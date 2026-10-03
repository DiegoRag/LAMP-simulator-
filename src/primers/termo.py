"""Cálculos termodinâmicos para primers LAMP: Tm, GC, ΔG de extremidade, grampos e dímeros."""

from dataclasses import dataclass
import primer3
from Bio.SeqUtils.MeltingTemp import DNA_NN3
from primers.modelo import normalizar

_COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


@dataclass(frozen=True)
class Condicoes:
    """Condições físico-químicas de reação para cálculo de Tm e ΔG."""

    na_mM: float = 50.0
    mg_mM: float = 8.0
    dntp_mM: float = 1.4
    temp_C: float = 37.0

    @classmethod
    def de_config(cls, cfg: dict) -> "Condicoes":
        """Constrói Condicoes a partir de dicionário com seção 'reacao'."""
        reacao = cfg["reacao"]
        return cls(
            na_mM=float(reacao["na_mM"]),
            mg_mM=float(reacao["mg_mM"]),
            dntp_mM=float(reacao["dntp_mM"]),
            temp_C=float(reacao["temp_C"]),
        )


def gc_pct(seq: str) -> float:
    """Retorna o percentual de G+C na sequência, de 0 a 100."""
    s = normalizar(seq)
    if not s:
        return 0.0
    return (s.count("G") + s.count("C")) / len(s) * 100.0


def tm(seq: str, cond: Condicoes, conc_nM: float) -> float:
    """Calcula a temperatura de fusão (Tm em °C) via primer3.calc_tm."""
    s = normalizar(seq)
    return float(
        primer3.calc_tm(
            s,
            mv_conc=cond.na_mM,
            dv_conc=cond.mg_mM,
            dntp_conc=cond.dntp_mM,
            dna_conc=conc_nM,
        )
    )


@dataclass(frozen=True)
class Estrutura:
    """Resultado de predição de estrutura secundária (grampo ou dímero)."""

    tm_c: float
    dg_kcal: float
    encontrada: bool


def grampo(seq: str, cond: Condicoes, conc_nM: float = 50.0) -> Estrutura:
    """Prediz formação de grampo (hairpin) via primer3.calc_hairpin."""
    s = normalizar(seq)
    res = primer3.calc_hairpin(
        s,
        mv_conc=cond.na_mM,
        dv_conc=cond.mg_mM,
        dntp_conc=cond.dntp_mM,
        dna_conc=conc_nM,
        temp_c=cond.temp_C,
    )
    return Estrutura(
        tm_c=float(res.tm),
        dg_kcal=float(res.dg) / 1000.0,
        encontrada=bool(res.structure_found),
    )


def homodimero(seq: str, cond: Condicoes, conc_nM: float = 50.0) -> Estrutura:
    """Prediz formação de homodímero via primer3.calc_homodimer."""
    s = normalizar(seq)
    res = primer3.calc_homodimer(
        s,
        mv_conc=cond.na_mM,
        dv_conc=cond.mg_mM,
        dntp_conc=cond.dntp_mM,
        dna_conc=conc_nM,
        temp_c=cond.temp_C,
    )
    return Estrutura(
        tm_c=float(res.tm),
        dg_kcal=float(res.dg) / 1000.0,
        encontrada=bool(res.structure_found),
    )


def heterodimero(a: str, b: str, cond: Condicoes, conc_nM: float = 50.0) -> Estrutura:
    """Prediz formação de heterodímero entre duas sequências via primer3.calc_heterodimer."""
    s_a = normalizar(a)
    s_b = normalizar(b)
    res = primer3.calc_heterodimer(
        s_a,
        s_b,
        mv_conc=cond.na_mM,
        dv_conc=cond.mg_mM,
        dntp_conc=cond.dntp_mM,
        dna_conc=conc_nM,
        temp_c=cond.temp_C,
    )
    return Estrutura(
        tm_c=float(res.tm),
        dg_kcal=float(res.dg) / 1000.0,
        encontrada=bool(res.structure_found),
    )


def dg_extremidade(seq: str, lado: str, bases: int, cond: Condicoes) -> float:
    """Calcula o ΔG de empilhamento (kcal/mol) das bases terminais da extremidade.

    Sem termos de iniciação: o critério do guia de desenho trata da estabilidade
    de empilhamento da extremidade, e a iniciação é constante por duplex.

    Parâmetros:
        seq: sequência de DNA.
        lado: "3" para extremidade 3' (últimas bases), "5" para extremidade 5' (primeiras bases).
        bases: quantidade de bases da extremidade a considerar.
        cond: condições de reação (temperatura em Celsius).

    Retorna:
        ΔG em kcal/mol calculado por nearest-neighbor usando DNA_NN3 (SantaLucia & Hicks 2004).
    """
    if lado not in ("3", "5"):
        raise ValueError(f"lado deve ser '3' ou '5', recebido: '{lado}'")

    s = normalizar(seq)
    if len(s) < bases:
        raise ValueError(
            f"Comprimento da sequência ({len(s)}) é menor que o número de bases ({bases})."
        )
    if bases < 1:
        raise ValueError(f"Número de bases ({bases}) deve ser >= 1.")

    subseq = s[-bases:] if lado == "3" else s[:bases]

    dh_total = 0.0
    ds_total = 0.0

    for i in range(len(subseq) - 1):
        x = subseq[i]
        y = subseq[i + 1]
        cx = _COMP[x]
        cy = _COMP[y]
        chave = f"{x}{y}/{cx}{cy}"
        if chave not in DNA_NN3:
            chave = chave[::-1]
        dh, ds = DNA_NN3[chave]
        dh_total += dh
        ds_total += ds

    t_kelvin = cond.temp_C + 273.15
    return dh_total - t_kelvin * (ds_total / 1000.0)
