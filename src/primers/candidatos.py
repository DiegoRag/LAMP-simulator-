"""Candidatos do desenho LAMP: regras por janela, estruturas na reação e teto por janela do alvo.

Janela com base ambígua (IUPAC ou N) não vira candidato: o desenho não propõe primer degenerado.
"""

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from primer3.thermoanalysis import ThermoAnalysis

from primers.janelas import IndiceAlvo
from primers.modelo import Regiao, rc
from primers.termo import Condicoes

TIPOS_BASE: tuple[str, ...] = ("F3", "F2", "F1", "B1c", "B2c", "B3c")

INFO_TIPO = {
    # tipo: (nome_primer, eh_rc, chave_tm, chave_conc, eh_loop)
    "F3": ("F3", False, "F3_B3", "externos", False),
    "F2": ("F2", False, "F2_B2", "internos", False),
    "F1": ("F1c", True, "F1c_B1c", "internos", False),
    "B1c": ("B1c", False, "F1c_B1c", "internos", False),
    "B2c": ("B2", True, "F2_B2", "internos", False),
    "B3c": ("B3", True, "F3_B3", "externos", False),
    "LF": ("LF", True, "loops", "loops", True),
    "LB": ("LB", False, "loops", "loops", True),
}

# Chave de IndiceAlvo.lote com o ΔG da ponta que a regra de extremidade avalia em cada tipo: a 3′ do
# primer (a 5′ em F1c e B1c), vista na janela da fita +; "_rc" quando o primer é o complemento reverso.
PONTA_TIPO = {
    "F3": "dg_fim", "F2": "dg_fim", "F1": "dg_fim_rc", "B1c": "dg_ini",
    "B2c": "dg_ini_rc", "B3c": "dg_ini_rc", "LF": "dg_ini_rc", "LB": "dg_fim",
}

# Regras na ordem em que são aplicadas; cada janela descartada conta só na primeira que a reprovou.
# "teto": passou em tudo, mas a sua janela do alvo já tinha candidatos_por_janela melhores.
REGRAS: tuple[str, ...] = (
    "ambigua", "gc_pct", "tm_c", "extremidade", "repetida", "grampo", "homodimero", "teto",
)


@dataclass(frozen=True, slots=True)
class Candidato:
    """Candidato individual a uma região de primer LAMP."""

    tipo: str
    regiao: Regiao
    primer: str
    tm_c: float
    desvio_tm: float


@lru_cache(maxsize=16)
def _analisador(cond: Condicoes, conc_nM: float) -> ThermoAnalysis:
    """primer3 nas condições da reação, os mesmos parâmetros de termo.grampo e termo.homodimero.

    Reaproveitado entre chamadas: o desenho testa grampo em centenas de milhares de janelas ACGT, e
    termo.grampo gastaria mais tempo normalizando a sequência do que o próprio primer3.
    """
    return ThermoAnalysis(
        mv_conc=cond.monovalentes_mM,
        dv_conc=cond.mg_mM,
        dntp_conc=cond.dntp_total_mM,
        dna_conc=conc_nM,
        temp_c=cond.temp_C,
    )


def reprova_grampo(primer: str, tipo: str, cfg: dict, cond: Condicoes) -> bool:
    """True se o grampo do primer (ACGT; primer3, na reação) tem Tm >= [reacao] temp_C − grampo_margem_C."""
    res = _analisador(cond, float(cfg["concentracao_nM"][INFO_TIPO[tipo][3]])).calc_hairpin(primer)
    limite = cond.temp_C - float(cfg["desenho"]["grampo_margem_C"])
    return bool(res.structure_found) and res.tm >= limite


def reprova_homodimero(primer: str, tipo: str, cfg: dict, cond: Condicoes) -> bool:
    """True se o homodímero (só ACGT; à concentração do primer no tubo) tem Tm >= temp_C − margem."""
    res = _analisador(cond, float(cfg["concentracao_nM"][INFO_TIPO[tipo][3]])).calc_homodimer(primer)
    limite = cond.temp_C - float(cfg["desenho"]["homodimero_margem_C"])
    return bool(res.structure_found) and res.tm >= limite


def gerar_candidatos(
    alvo: str,
    cfg: dict,
    tipos: tuple[str, ...] = TIPOS_BASE,
    dentro: Regiao | None = None,
    *,
    indice: IndiceAlvo | None = None,
    eliminados: dict[str, dict[str, int]] | None = None,
    testar_homodimero: bool = True,
) -> dict[str, list[Candidato]]:
    """Gera os candidatos de cada tipo, ordenados por (regiao.ini, regiao.fim).

    1. Janelas [ini, ini+L) para cada L da faixa do tipo; com `dentro`, só as contidas em `dentro`.
    2. Regras, nesta ordem: base ambígua, GC, Tm de referência, ponta e repetição no alvo ou no rc(alvo).
    3. Sem `dentro`: teto por janela. O alvo é dividido em janelas de [desenho] janela_nt nt pelo ini; em
       cada uma, na ordem (desvio_tm, ini, fim), ficam os primeiros candidatos_por_janela que passam no
       grampo (e no homodímero, se testar_homodimero). Com `dentro` (laços): sem teto.
    `eliminados`, se dado, recebe {tipo: {regra: janelas descartadas}} com todas as chaves de REGRAS.
    """
    if indice is None:
        indice = IndiceAlvo(alvo, cfg)
    cond = Condicoes.de_config(cfg)
    largura = int(cfg["desenho"]["janela_nt"])
    por_janela = int(cfg["desenho"]["candidatos_por_janela"]) if dentro is None else -1
    gc_min = float(cfg["gc_pct"]["min"])
    gc_max = float(cfg["gc_pct"]["max"])
    dg_max = float(cfg["extremidade"]["dg_max"])

    resultado: dict[str, list[Candidato]] = {}
    for tipo in tipos:
        _, eh_rc, chave_tm, _, eh_loop = INFO_TIPO[tipo]
        c_min = int(cfg["comprimento"]["loop_min" if eh_loop else "min"])
        c_max = int(cfg["comprimento"]["loop_max" if eh_loop else "max"])
        tm_min, tm_max = float(cfg["tm_c"][chave_tm][0]), float(cfg["tm_c"][chave_tm][1])
        centro = (tm_min + tm_max) / 2.0
        cont = dict.fromkeys(REGRAS, 0)

        partes: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
        for L in range(c_min, c_max + 1):
            ini_min, ini_max = (dentro.ini, dentro.fim - L) if dentro is not None else (0, indice.n - L)
            ini_max = min(ini_max, indice.n - L)
            if L < indice.bases_ponta or ini_max < ini_min:
                continue
            lt = indice.lote(L, ini_min, ini_max)
            tm = lt["tm_rc"] if eh_rc else lt["tm"]
            restante = ~lt["ambigua"]
            cont["ambigua"] += int(np.count_nonzero(lt["ambigua"]))
            for regra, ok in (
                ("gc_pct", (lt["gc"] >= gc_min) & (lt["gc"] <= gc_max)),
                ("tm_c", (tm >= tm_min) & (tm <= tm_max)),
                ("extremidade", lt[PONTA_TIPO[tipo]] <= dg_max),
            ):
                cont[regra] += int(np.count_nonzero(restante & ~ok))
                restante &= ok
            idx = np.nonzero(restante)[0]
            inis = idx + ini_min
            unica = np.ones(len(idx), dtype=bool)
            for j in np.nonzero(indice.suspeito[inis])[0]:
                unica[j] = indice.unica(int(inis[j]), L)
            cont["repetida"] += int(np.count_nonzero(~unica))
            partes.append((inis[unica], inis[unica] + L, tm[idx][unica]))

        inis = np.concatenate([p[0] for p in partes]) if partes else np.zeros(0, dtype=np.int64)
        fims = np.concatenate([p[1] for p in partes]) if partes else np.zeros(0, dtype=np.int64)
        tms = np.concatenate([p[2] for p in partes]) if partes else np.zeros(0)
        desvios = np.abs(tms - centro)
        janelas = inis // largura if dentro is None else np.zeros(len(inis), dtype=np.int64)

        ordem = np.lexsort((fims, inis, desvios, janelas))
        inis_o, fims_o = inis[ordem].tolist(), fims[ordem].tolist()
        tms_o, desvios_o = tms[ordem].tolist(), desvios[ordem].tolist()
        # início de cada janela do alvo na ordem (com `dentro`, uma janela só)
        janelas_o = janelas[ordem]
        cortes = np.flatnonzero(np.diff(janelas_o)) + 1 if len(ordem) else np.zeros(0, dtype=np.int64)
        limites = [0, *cortes.tolist(), len(ordem)]

        candidatos: list[Candidato] = []
        for a, z in zip(limites[:-1], limites[1:], strict=True):
            aceitos = 0
            for k in range(a, z):
                if aceitos == por_janela:
                    cont["teto"] += z - k
                    break
                ini, fim = inis_o[k], fims_o[k]
                primer = rc(indice.alvo[ini:fim]) if eh_rc else indice.alvo[ini:fim]
                if reprova_grampo(primer, tipo, cfg, cond):
                    cont["grampo"] += 1
                    continue
                if testar_homodimero and reprova_homodimero(primer, tipo, cfg, cond):
                    cont["homodimero"] += 1
                    continue
                aceitos += 1
                candidatos.append(Candidato(tipo, Regiao(ini, fim), primer, tms_o[k], desvios_o[k]))
        candidatos.sort(key=lambda c: (c.regiao.ini, c.regiao.fim))
        resultado[tipo] = candidatos
        if eliminados is not None:
            eliminados[tipo] = cont
    return resultado
