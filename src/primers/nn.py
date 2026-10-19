"""Modelo de vizinhos mais próximos (DNA/DNA) usado por primers.termo: tabelas e fórmulas da literatura.

Vizinhos: DNA_NN3 do Biopython = Allawi & SantaLucia (1997), mesmos valores da tabela unificada de
SantaLucia (1998). Todas as funções recebem sequências só com A, C, G e T (a expansão IUPAC fica em termo).
"""

import math

from Bio.SeqUtils.MeltingTemp import DNA_NN3

R_CAL = 1.987  # cal/(K·mol)
T37_K = 310.15

_COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

# SantaLucia (1998) PNAS 95:1460, tabela 2: ΔG°37 em NaCl 1 M (kcal/mol) por vizinho 5'->3'.
DG37_NN: dict[str, float] = {
    "AA": -1.00, "TT": -1.00, "AT": -0.88, "TA": -0.58,
    "CA": -1.45, "TG": -1.45, "GT": -1.44, "AC": -1.44,
    "CT": -1.28, "AG": -1.28, "GA": -1.30, "TC": -1.30,
    "CG": -2.17, "GC": -2.24, "GG": -1.84, "CC": -1.84,
}
DG37_INIT_GC: float = 0.98
DG37_INIT_AT: float = 1.03
DG37_SIMETRIA: float = 0.43

# SantaLucia & Hicks (2004) Annu Rev Biophys Biomol Struct 33:415, tabela 4: ΔG°37 do laço de
# grampo em NaCl 1 M (kcal/mol), com ΔH = 0.
DG37_LACO_GRAMPO: dict[int, float] = {
    3: 3.5, 4: 3.5, 5: 3.3, 6: 4.0, 7: 4.2, 8: 4.3, 9: 4.5, 10: 4.6,
    12: 5.0, 14: 5.1, 16: 5.3, 18: 5.5, 20: 5.7, 25: 6.1, 30: 6.3,
}


def autocomplementar(v: str) -> bool:
    """True se v é igual ao seu complemento reverso."""
    return v == "".join(_COMP[b] for b in reversed(v))


def fracao_gc(v: str) -> float:
    """Fração de G + C em v, de 0 a 1."""
    return sum(1 for b in v if b in "GC") / len(v)


def dh_ds_vizinhos(v: str) -> tuple[float, float]:
    """Soma de ΔH (kcal/mol) e ΔS (cal/(K·mol)) dos vizinhos de v na tabela DNA_NN3, sem iniciação."""
    dh = 0.0
    ds = 0.0
    for i in range(len(v) - 1):
        x = v[i]
        y = v[i + 1]
        chave = f"{x}{y}/{_COMP[x]}{_COMP[y]}"
        if chave not in DNA_NN3:
            chave = chave[::-1]
        h, s = DNA_NN3[chave]
        dh += h
        ds += s
    return dh, ds


def dh_ds_duplex(v: str) -> tuple[float, float]:
    """ΔH e ΔS do duplex perfeito v/complemento em NaCl 1 M: vizinhos, iniciação nas pontas e simetria."""
    dh, ds = dh_ds_vizinhos(v)
    for ponta in (v[0], v[-1]):
        h, s = DNA_NN3["init_G/C"] if ponta in "GC" else DNA_NN3["init_A/T"]
        dh += h
        ds += s
    if autocomplementar(v):
        h, s = DNA_NN3["sym"]
        dh += h
        ds += s
    return dh, ds


def dg37_ponta(v: str) -> float:
    """ΔG°37 (kcal/mol, NaCl 1 M) do duplex v/complemento: SantaLucia (1998) com iniciação e simetria."""
    dg = sum(DG37_NN[v[i : i + 2]] for i in range(len(v) - 1))
    dg += DG37_INIT_GC if v[0] in "GC" else DG37_INIT_AT
    dg += DG37_INIT_GC if v[-1] in "GC" else DG37_INIT_AT
    if autocomplementar(v):
        dg += DG37_SIMETRIA
    return dg


def correcao_owczarzy(f_gc: float, n_bp: int, mono_mM: float, mg_livre_mM: float) -> float:
    """Termo somado a 1/Tm (1/K): Owczarzy et al. (2008) eq. 16, ou (2004) eq. 22 sem Mg²⁺."""
    mono = mono_mM / 1000.0
    mg = mg_livre_mM / 1000.0
    if mono <= 0.0 and mg <= 0.0:
        raise ValueError("Sem cátions na reação: Na⁺ + K⁺ + NH₄⁺ ou o Mg²⁺ livre precisa ser > 0.")
    if mg <= 0.0 or (mono > 0.0 and math.sqrt(mg) / mono < 0.22):
        ln_mono = math.log(mono)
        return (4.29 * f_gc - 3.95) * 1e-5 * ln_mono + 9.40e-6 * ln_mono**2
    a, b, c, d = 3.92e-5, -9.11e-6, 6.26e-5, 1.42e-5
    e, f, g = -4.82e-4, 5.25e-4, 8.31e-5
    if mono > 0.0 and math.sqrt(mg) / mono < 6.0:
        ln_mono = math.log(mono)
        a = 3.92e-5 * (0.843 - 0.352 * math.sqrt(mono) * ln_mono)
        d = 1.42e-5 * (1.279 - 4.03e-3 * ln_mono - 8.03e-3 * ln_mono**2)
        g = 8.31e-5 * (0.486 - 0.258 * ln_mono + 5.25e-3 * ln_mono**3)
    ln_mg = math.log(mg)
    return (
        a
        + b * ln_mg
        + f_gc * (c + d * ln_mg)
        + (1.0 / (2.0 * (n_bp - 1))) * (e + f * ln_mg + g * ln_mg**2)
    )


def dg_laco_grampo(laco_nt: int) -> float:
    """ΔG°37 (kcal/mol, NaCl 1 M) do laço de grampo de laco_nt nt (SantaLucia & Hicks 2004)."""
    if laco_nt < 3:
        raise ValueError(f"Laço de grampo precisa de ao menos 3 nt (recebido: {laco_nt}).")
    if laco_nt in DG37_LACO_GRAMPO:
        return DG37_LACO_GRAMPO[laco_nt]
    if laco_nt > 30:
        return DG37_LACO_GRAMPO[30] + 2.44 * (R_CAL / 1000.0) * T37_K * math.log(laco_nt / 30.0)
    menor = max(t for t in DG37_LACO_GRAMPO if t < laco_nt)
    maior = min(t for t in DG37_LACO_GRAMPO if t > laco_nt)
    g0 = DG37_LACO_GRAMPO[menor]
    g1 = DG37_LACO_GRAMPO[maior]
    return g0 + (g1 - g0) * (laco_nt - menor) / (maior - menor)
