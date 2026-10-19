"""Termodinâmica de primers LAMP: Tm de referência (PrimerExplorer V5) e da reação, ponta e estruturas.

Tabelas e fórmulas do modelo de vizinhos ficam em primers.nn.
"""

import math
from dataclasses import dataclass

import primer3

from primers.modelo import BASES_IUPAC, expandir_iupac, normalizar
from primers.nn import (
    R_CAL,
    T37_K,
    autocomplementar,
    correcao_owczarzy,
    dg37_ponta,
    dg_laco_grampo,
    dh_ds_duplex,
    dh_ds_vizinhos,
    fracao_gc,
)

MAX_NT_ESTRUTURA: int = 60


@dataclass(frozen=True)
class Condicoes:
    """Condições da reação LAMP: Tm da reação (informativa), grampos e dímeros."""

    na_mM: float = 0.0
    k_mM: float = 50.0
    nh4_mM: float = 20.0
    mg_mM: float = 8.0
    dntp_cada_mM: float = 1.4
    temp_C: float = 65.0
    modo_iupac: str = "media"
    max_expansoes_iupac: int = 64

    @property
    def monovalentes_mM(self) -> float:
        """Na⁺ + K⁺ + NH₄⁺, em mM."""
        return round(self.na_mM + self.k_mM + self.nh4_mM, 9)

    @property
    def dntp_total_mM(self) -> float:
        """dNTP total (os quatro), em mM: 4 × dntp_cada_mM."""
        return round(4.0 * self.dntp_cada_mM, 9)

    @property
    def mg_livre_mM(self) -> float:
        """Mg²⁺ livre, em mM: Mg total menos o dNTP total (cada dNTP quela um Mg²⁺), nunca < 0."""
        return max(0.0, round(self.mg_mM - self.dntp_total_mM, 9))

    @classmethod
    def de_config(cls, cfg: dict) -> "Condicoes":
        """Constrói Condicoes a partir de cfg["reacao"] (e cfg["iupac"], opcional)."""
        reacao = cfg["reacao"]
        if "dntp_cada_mM" not in reacao:
            if "dntp_mM" in reacao:
                raise ValueError(
                    "A chave reacao.dntp_mM foi substituída por reacao.dntp_cada_mM "
                    "(mM de cada dNTP; o total é 4 × esse valor)."
                )
            raise ValueError("Falta a chave reacao.dntp_cada_mM (mM de cada dNTP).")
        iupac = cfg.get("iupac", {})
        return cls(
            na_mM=float(reacao.get("na_mM", 0.0)),
            k_mM=float(reacao.get("k_mM", 0.0)),
            nh4_mM=float(reacao.get("nh4_mM", 0.0)),
            mg_mM=float(reacao["mg_mM"]),
            dntp_cada_mM=float(reacao["dntp_cada_mM"]),
            temp_C=float(reacao["temp_C"]),
            modo_iupac=str(iupac.get("modo_tm", "media")),
            max_expansoes_iupac=int(iupac.get("max_expansoes", 64)),
        )


@dataclass(frozen=True)
class CondicoesReferencia:
    """Convenção do PrimerExplorer V5 para a Tm de referência (a das faixas do guia da Eiken)."""

    oligo_nM: float = 100.0
    na_mM: float = 50.0
    mg_mM: float = 4.0
    modo_iupac: str = "media"
    max_expansoes_iupac: int = 64

    @property
    def na_equivalente_M(self) -> float:
        """[Na⁺] + 4·√[Mg²⁺], em mol/L (conversão do PrimerExplorer)."""
        return self.na_mM / 1000.0 + 4.0 * math.sqrt(self.mg_mM / 1000.0)

    @classmethod
    def de_config(cls, cfg: dict) -> "CondicoesReferencia":
        """Constrói a convenção a partir de cfg["tm_referencia"] (padrões 100 nM, 50 mM, 4 mM)."""
        ref = cfg.get("tm_referencia", {})
        iupac = cfg.get("iupac", {})
        return cls(
            oligo_nM=float(ref.get("oligo_nM", 100.0)),
            na_mM=float(ref.get("na_mM", 50.0)),
            mg_mM=float(ref.get("mg_mM", 4.0)),
            modo_iupac=str(iupac.get("modo_tm", "media")),
            max_expansoes_iupac=int(iupac.get("max_expansoes", 64)),
        )


def gc_pct(seq: str) -> float:
    """Retorna o percentual ponderado de G+C na sequência, de 0 a 100, considerando IUPAC."""
    s = normalizar(seq)
    if not s:
        return 0.0
    soma_gc = sum(
        sum(1 for b in BASES_IUPAC[c] if b in ("G", "C")) / len(BASES_IUPAC[c])
        for c in s
    )
    return (soma_gc / len(s)) * 100.0


def _modo(modo_cond: str, modo_iupac: str | None) -> str:
    modo = modo_cond if modo_iupac is None else modo_iupac
    if modo not in ("media", "pior_caso"):
        raise ValueError(f"modo_iupac inválido: {modo!r}. Esperado 'media' ou 'pior_caso'.")
    return modo


def _agregar(valores: list[float], modo: str, pior: str) -> float:
    if len(valores) == 1:
        return valores[0]
    if modo == "pior_caso":
        return min(valores) if pior == "min" else max(valores)
    return sum(valores) / len(valores)


def _variantes_tm(seq: str, max_expansoes: int) -> list[str]:
    s = normalizar(seq)
    if not s:
        raise ValueError("Sequência vazia não possui Tm.")
    if len(s) < 2:
        raise ValueError("Sequência com menos de 2 nt não possui Tm.")
    return expandir_iupac(s, max_expansoes=max_expansoes)


def tm_referencia(seq: str, ref: CondicoesReferencia, *, modo_iupac: str | None = None) -> float:
    """Tm de referência (°C), PrimerExplorer V5: 1000·ΔH/(ΔS + R·ln(C/4)) − 273,15 + 16,6·log10(Na_eq)."""
    modo = _modo(ref.modo_iupac, modo_iupac)
    variantes = _variantes_tm(seq, ref.max_expansoes_iupac)
    na_eq = ref.na_equivalente_M
    if na_eq <= 0.0:
        raise ValueError("Na⁺ + 4·√Mg²⁺ da referência precisa ser > 0.")
    c_m = ref.oligo_nM * 1e-9
    tms: list[float] = []
    for v in variantes:
        dh, ds = dh_ds_duplex(v)
        fator = c_m if autocomplementar(v) else c_m / 4.0
        tms.append(1000.0 * dh / (ds + R_CAL * math.log(fator)) - 273.15 + 16.6 * math.log10(na_eq))
    return _agregar(tms, modo, "min")


def tm(seq: str, cond: Condicoes, conc_nM: float, *, modo_iupac: str | None = None) -> float:
    """Tm da reação (°C), informativa: duplex primer/alvo a conc_nM, sal de Owczarzy et al. (2008)."""
    modo = _modo(cond.modo_iupac, modo_iupac)
    variantes = _variantes_tm(seq, cond.max_expansoes_iupac)
    c_m = conc_nM * 1e-9
    tms: list[float] = []
    for v in variantes:
        dh, ds = dh_ds_duplex(v)
        # primer em excesso sobre o alvo: C_A − C_B/2 ≈ C do primer (SantaLucia 1998), sem o /4
        tm_1m = 1000.0 * dh / (ds + R_CAL * math.log(c_m))
        corr = correcao_owczarzy(fracao_gc(v), len(v), cond.monovalentes_mM, cond.mg_livre_mM)
        tms.append(1.0 / (1.0 / tm_1m + corr) - 273.15)
    return _agregar(tms, modo, "min")


def tm_dobramento(seq: str, laco_nt: int, cond: Condicoes, *, modo_iupac: str | None = None) -> float:
    """Tm da reação (°C) de F1c/B1c como dobramento: haste seq/complemento + laço, sem concentração."""
    modo = _modo(cond.modo_iupac, modo_iupac)
    variantes = _variantes_tm(seq, cond.max_expansoes_iupac)
    ds_laco = -1000.0 * dg_laco_grampo(laco_nt) / T37_K
    tms: list[float] = []
    for v in variantes:
        dh, ds = dh_ds_vizinhos(v)
        tm_1m = 1000.0 * dh / (ds + ds_laco)
        corr = correcao_owczarzy(fracao_gc(v), len(v), cond.monovalentes_mM, cond.mg_livre_mM)
        tms.append(1.0 / (1.0 / tm_1m + corr) - 273.15)
    return _agregar(tms, modo, "min")


@dataclass(frozen=True)
class Estrutura:
    """Resultado de predição de estrutura secundária (grampo ou dímero)."""

    tm_c: float
    dg_kcal: float
    encontrada: bool


def grampo(seq: str, cond: Condicoes, conc_nM: float = 50.0) -> Estrutura:
    """Prediz formação de grampo (hairpin) via primer3.calc_hairpin, na temperatura da reação."""
    s = normalizar(seq)
    if not s:
        raise ValueError("Sequência vazia para cálculo de estrutura secundária.")
    if len(s) > MAX_NT_ESTRUTURA:
        return Estrutura(tm_c=0.0, dg_kcal=0.0, encontrada=False)

    variantes = expandir_iupac(s, max_expansoes=cond.max_expansoes_iupac)
    estruturas: list[Estrutura] = []
    for v in variantes:
        try:
            res = primer3.calc_hairpin(
                v,
                mv_conc=cond.monovalentes_mM,
                dv_conc=cond.mg_mM,
                dntp_conc=cond.dntp_total_mM,
                dna_conc=conc_nM,
                temp_c=cond.temp_C,
            )
            estruturas.append(
                Estrutura(
                    tm_c=float(res.tm),
                    dg_kcal=float(res.dg) / 1000.0,
                    encontrada=bool(res.structure_found),
                )
            )
        except RuntimeError:
            estruturas.append(Estrutura(tm_c=0.0, dg_kcal=0.0, encontrada=False))

    return min(estruturas, key=lambda e: (not e.encontrada, e.dg_kcal))


def homodimero(seq: str, cond: Condicoes, conc_nM: float = 50.0) -> Estrutura:
    """Prediz formação de homodímero via primer3.calc_homodimer, na temperatura da reação."""
    s = normalizar(seq)
    if not s:
        raise ValueError("Sequência vazia para cálculo de estrutura secundária.")
    if len(s) > MAX_NT_ESTRUTURA:
        return Estrutura(tm_c=0.0, dg_kcal=0.0, encontrada=False)

    variantes = expandir_iupac(s, max_expansoes=cond.max_expansoes_iupac)
    estruturas: list[Estrutura] = []
    for v in variantes:
        try:
            res = primer3.calc_homodimer(
                v,
                mv_conc=cond.monovalentes_mM,
                dv_conc=cond.mg_mM,
                dntp_conc=cond.dntp_total_mM,
                dna_conc=conc_nM,
                temp_c=cond.temp_C,
            )
            estruturas.append(
                Estrutura(
                    tm_c=float(res.tm),
                    dg_kcal=float(res.dg) / 1000.0,
                    encontrada=bool(res.structure_found),
                )
            )
        except RuntimeError:
            estruturas.append(Estrutura(tm_c=0.0, dg_kcal=0.0, encontrada=False))

    return min(estruturas, key=lambda e: (not e.encontrada, e.dg_kcal))


def heterodimero(a: str, b: str, cond: Condicoes, conc_nM: float = 50.0) -> Estrutura:
    """Prediz heterodímero (melhor alinhamento em qualquer ponto) via primer3.calc_heterodimer."""
    s_a = normalizar(a)
    s_b = normalizar(b)
    if not s_a or not s_b:
        raise ValueError("Sequência vazia para cálculo de estrutura secundária.")
    if len(s_a) > MAX_NT_ESTRUTURA or len(s_b) > MAX_NT_ESTRUTURA:
        return Estrutura(tm_c=0.0, dg_kcal=0.0, encontrada=False)

    limite_cada = max(1, int(cond.max_expansoes_iupac ** 0.5))
    vars_a = expandir_iupac(s_a, max_expansoes=limite_cada)
    vars_b = expandir_iupac(s_b, max_expansoes=limite_cada)

    estruturas: list[Estrutura] = []
    for va in vars_a:
        for vb in vars_b:
            try:
                res = primer3.calc_heterodimer(
                    va,
                    vb,
                    mv_conc=cond.monovalentes_mM,
                    dv_conc=cond.mg_mM,
                    dntp_conc=cond.dntp_total_mM,
                    dna_conc=conc_nM,
                    temp_c=cond.temp_C,
                )
                estruturas.append(
                    Estrutura(
                        tm_c=float(res.tm),
                        dg_kcal=float(res.dg) / 1000.0,
                        encontrada=bool(res.structure_found),
                    )
                )
            except RuntimeError:
                estruturas.append(Estrutura(tm_c=0.0, dg_kcal=0.0, encontrada=False))

    return min(estruturas, key=lambda e: (not e.encontrada, e.dg_kcal))


def dimero_ponta3(a: str, b: str, cond: Condicoes, conc_nM: float = 50.0) -> Estrutura:
    """Dímero mais estável com a ponta 3′ de a ou de b pareada (calc_end_stability nos dois sentidos)."""
    s_a = normalizar(a)
    s_b = normalizar(b)
    if not s_a or not s_b:
        raise ValueError("Sequência vazia para cálculo de estrutura secundária.")
    if len(s_a) > MAX_NT_ESTRUTURA or len(s_b) > MAX_NT_ESTRUTURA:
        return Estrutura(tm_c=0.0, dg_kcal=0.0, encontrada=False)

    limite_cada = max(1, int(cond.max_expansoes_iupac ** 0.5))
    vars_a = expandir_iupac(s_a, max_expansoes=limite_cada)
    vars_b = expandir_iupac(s_b, max_expansoes=limite_cada)

    estruturas: list[Estrutura] = []
    for va in vars_a:
        for vb in vars_b:
            for x, y in ((va, vb), (vb, va)):
                try:
                    res = primer3.calc_end_stability(
                        x,
                        y,
                        mv_conc=cond.monovalentes_mM,
                        dv_conc=cond.mg_mM,
                        dntp_conc=cond.dntp_total_mM,
                        dna_conc=conc_nM,
                        temp_c=cond.temp_C,
                    )
                    estruturas.append(
                        Estrutura(
                            tm_c=float(res.tm),
                            dg_kcal=float(res.dg) / 1000.0,
                            encontrada=bool(res.structure_found),
                        )
                    )
                except RuntimeError:
                    estruturas.append(Estrutura(tm_c=0.0, dg_kcal=0.0, encontrada=False))

    return min(estruturas, key=lambda e: (not e.encontrada, e.dg_kcal))


def dg_extremidade(
    seq: str,
    lado: str,
    bases: int,
    cond: Condicoes | CondicoesReferencia,
    *,
    modo_iupac: str | None = None,
) -> float:
    """ΔG (kcal/mol) da ponta, PrimerExplorer V5: ΔG°37 em NaCl 1 M das `bases` terminais (nn.dg37_ponta)."""
    if lado not in ("3", "5"):
        raise ValueError(f"lado deve ser '3' ou '5', recebido: '{lado}'")

    s = normalizar(seq)
    if len(s) < bases:
        raise ValueError(
            f"Comprimento da sequência ({len(s)}) é menor que o número de bases ({bases})."
        )
    if bases < 1:
        raise ValueError(f"Número de bases ({bases}) deve ser >= 1.")

    modo = _modo(cond.modo_iupac, modo_iupac)
    subseq = s[-bases:] if lado == "3" else s[:bases]
    variantes = expandir_iupac(subseq, max_expansoes=cond.max_expansoes_iupac)
    return _agregar([dg37_ponta(v) for v in variantes], modo, "max")
