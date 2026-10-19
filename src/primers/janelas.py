"""Janelas do alvo para o desenho: índice de k-meros das duas fitas e propriedades em lote (numpy).

GC, Tm de referência e ΔG da ponta saem iguais, bit a bit, às de primers.termo (mesmas tabelas, mesma ordem
de soma), para todas as janelas de um comprimento de uma vez; a unicidade fica O(N), sem str.find por janela.
"""

import math

import numpy as np
from Bio.SeqUtils.MeltingTemp import DNA_NN3

from primers.modelo import normalizar, rc
from primers.nn import R_CAL, dg37_ponta
from primers.termo import CondicoesReferencia

_COD_BASE = np.full(256, -1, dtype=np.int64)
for _i, _b in enumerate("ACGT"):
    _COD_BASE[ord(_b)] = _i
_COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _tabela_vizinhos() -> tuple[np.ndarray, np.ndarray]:
    """ΔH e ΔS (DNA_NN3) do vizinho x→y no índice 4·x + y, com a mesma busca de nn.dh_ds_vizinhos."""
    dh = np.zeros(16)
    ds = np.zeros(16)
    for i, x in enumerate("ACGT"):
        for j, y in enumerate("ACGT"):
            chave = f"{x}{y}/{_COMP[x]}{_COMP[y]}"
            if chave not in DNA_NN3:
                chave = chave[::-1]
            dh[4 * i + j], ds[4 * i + j] = DNA_NN3[chave]
    return dh, ds


_DH_VIZ, _DS_VIZ = _tabela_vizinhos()
# iniciação por base terminal, no índice A=0, C=1, G=2, T=3
_DH_INIT = np.array([DNA_NN3[f"init_{p}"][0] for p in ("A/T", "G/C", "G/C", "A/T")])
_DS_INIT = np.array([DNA_NN3[f"init_{p}"][1] for p in ("A/T", "G/C", "G/C", "A/T")])


def codigos_kmer(cod: np.ndarray, k: int, reverso: bool) -> np.ndarray:
    """Código base 4 do k-mero alvo[p:p+k] para cada p (com reverso, do complemento reverso dele).

    `cod` tem A=0, C=1, G=2, T=3 e −1 para base ambígua; k-mero com base ambígua recebe −1.
    """
    n = len(cod) - k + 1
    if n <= 0:
        return np.zeros(0, dtype=np.int64)
    valido = np.ones(n, dtype=bool)
    codigo = np.zeros(n, dtype=np.int64)
    for t in range(k):
        parte = cod[k - 1 - t : k - 1 - t + n] if reverso else cod[t : t + n]
        valido &= parte >= 0
        codigo = codigo * 4 + ((3 - parte) if reverso else parte)
    return np.where(valido, codigo, -1)


class IndiceAlvo:
    """Alvo codificado, k-meros das duas fitas e ΔG das pontas, para avaliar as janelas em lote."""

    def __init__(self, alvo: str, cfg: dict) -> None:
        self.alvo = normalizar(alvo)
        self.n = len(self.alvo)
        self.ref = CondicoesReferencia.de_config(cfg)
        self.bases_ponta = int(cfg["extremidade"]["bases"])
        self.k = max(1, min(int(cfg["comprimento"]["min"]), int(cfg["comprimento"]["loop_min"]), 31))

        cod = _COD_BASE[np.frombuffer(self.alvo.encode("ascii"), dtype=np.uint8)]
        valida = cod >= 0
        base = np.where(valida, cod, 0)
        self.amb_acum = np.concatenate(([0], np.cumsum(~valida)))
        self.gc_acum = np.concatenate(([0], np.cumsum((cod == 1) | (cod == 2))))
        self.dh_init = _DH_INIT[base]
        self.ds_init = _DS_INIT[base]
        self.dh_viz = _DH_VIZ[4 * base[:-1] + base[1:]]
        self.ds_viz = _DS_VIZ[4 * base[:-1] + base[1:]]

        # ΔG°37 das `bases_ponta` bases a partir de cada posição, no sentido da fita + e em rc
        cod_ponta = codigos_kmer(cod, self.bases_ponta, reverso=False)
        pos_ok = np.flatnonzero(cod_ponta >= 0)
        _, primeiro, inverso = np.unique(cod_ponta[pos_ok], return_index=True, return_inverse=True)
        trechos = [self.alvo[p : p + self.bases_ponta] for p in pos_ok[primeiro].tolist()]
        self.dg_ponta = np.full(len(cod_ponta), np.nan)
        self.dg_ponta_rc = np.full(len(cod_ponta), np.nan)
        if trechos:
            self.dg_ponta[pos_ok] = np.array([dg37_ponta(t) for t in trechos])[inverso]
            self.dg_ponta_rc[pos_ok] = np.array([dg37_ponta(rc(t)) for t in trechos])[inverso]

        # unicidade: k-meros das duas fitas; só os que aparecem 2 ou mais vezes guardam as posições
        self.cod_mais = codigos_kmer(cod, self.k, reverso=False)
        cod_menos = codigos_kmer(cod, self.k, reverso=True)
        todos = np.concatenate((self.cod_mais[self.cod_mais >= 0], cod_menos[cod_menos >= 0]))
        unicos, contagens = np.unique(todos, return_counts=True)
        repetidos = unicos[contagens >= 2]
        self.suspeito = np.isin(self.cod_mais, repetidos) & (self.cod_mais >= 0)
        self.ocorrencias: dict[int, list[tuple[int, int]]] = {}
        for p in np.flatnonzero(self.suspeito).tolist():
            self.ocorrencias.setdefault(int(self.cod_mais[p]), []).append((0, p))
        for q in np.flatnonzero(np.isin(cod_menos, repetidos) & (cod_menos >= 0)).tolist():
            self.ocorrencias.setdefault(int(cod_menos[q]), []).append((1, q))

    def unica(self, ini: int, comprimento: int) -> bool:
        """True se alvo[ini:ini+comprimento] não aparece em outro ponto do alvo nem no rc(alvo).

        Um palíndromo (janela igual ao próprio complemento reverso) conta como repetido.
        """
        if not self.suspeito[ini]:
            return True
        w = self.alvo[ini : ini + comprimento]
        w_rc = rc(w)
        for fita, p in self.ocorrencias[int(self.cod_mais[ini])]:
            if fita == 0:
                if p != ini and self.alvo[p : p + comprimento] == w:
                    return False
            else:
                j = p - comprimento + self.k
                if j >= 0 and self.alvo[j : j + comprimento] == w_rc:
                    return False
        return True

    def lote(self, comprimento: int, ini_min: int, ini_max: int) -> dict[str, np.ndarray]:
        """Propriedades das janelas [ini, ini+comprimento), ini de ini_min a ini_max (inclusive).

        "ambigua" (bool), "gc" (%), "tm"/"tm_rc" (Tm de referência do trecho e do rc dele, como
        termo.tm_referencia sem autocomplementaridade), "dg_ini"/"dg_fim" (ΔG das bases_ponta bases do
        início/fim da janela) e "dg_ini_rc"/"dg_fim_rc" (as mesmas bases em rc). Filtre por "ambigua" antes.
        """
        L = comprimento
        a, z = ini_min, ini_max + 1
        amb = (self.amb_acum[a + L : z + L] - self.amb_acum[a:z]) > 0
        gc = ((self.gc_acum[a + L : z + L] - self.gc_acum[a:z]) / L) * 100.0
        dh, ds, dh_rc, ds_rc = (np.zeros(z - a) for _ in range(4))
        for t in range(L - 1):
            dh += self.dh_viz[a + t : z + t]
            ds += self.ds_viz[a + t : z + t]
            dh_rc += self.dh_viz[a + L - 2 - t : z + L - 2 - t]
            ds_rc += self.ds_viz[a + L - 2 - t : z + L - 2 - t]
        dh += self.dh_init[a:z]
        ds += self.ds_init[a:z]
        dh += self.dh_init[a + L - 1 : z + L - 1]
        ds += self.ds_init[a + L - 1 : z + L - 1]
        dh_rc += self.dh_init[a + L - 1 : z + L - 1]
        ds_rc += self.ds_init[a + L - 1 : z + L - 1]
        dh_rc += self.dh_init[a:z]
        ds_rc += self.ds_init[a:z]
        log_c = R_CAL * math.log(self.ref.oligo_nM * 1e-9 / 4.0)
        sal = 16.6 * math.log10(self.ref.na_equivalente_M)
        b = self.bases_ponta
        return {
            "ambigua": amb,
            "gc": gc,
            "tm": 1000.0 * dh / (ds + log_c) - 273.15 + sal,
            "tm_rc": 1000.0 * dh_rc / (ds_rc + log_c) - 273.15 + sal,
            "dg_ini": self.dg_ponta[a:z],
            "dg_fim": self.dg_ponta[a + L - b : z + L - b],
            "dg_ini_rc": self.dg_ponta_rc[a:z],
            "dg_fim_rc": self.dg_ponta_rc[a + L - b : z + L - b],
        }
