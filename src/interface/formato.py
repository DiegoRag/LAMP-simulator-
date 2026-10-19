"""Funções puras de formatação de texto para a interface web."""

from datetime import datetime

MENOS = "−"


def num(valor: float | int | None, casas: int = 2) -> str:
    """Formata valor numérico com vírgula decimal e sinal tipográfico, ou travessão se None."""
    if valor is None:
        return "—"
    s = f"{float(valor):.{casas}f}"
    if s.startswith("-") and float(s) == 0.0:
        s = s[1:]
    return s.replace("-", MENOS).replace(".", ",")


def num_curto(valor: float | int) -> str:
    """Formata número sem zeros à direita desnecessários, com vírgula e sinal tipográfico."""
    s = f"{float(valor):.10g}"
    if s.startswith("-") and float(s) == 0.0:
        s = s[1:]
    return s.replace("-", MENOS).replace(".", ",")


def milhar(n: int) -> str:
    """Formata inteiro usando ponto como separador de milhar."""
    s = f"{int(n):,}".replace(",", ".")
    if s.startswith("-"):
        return MENOS + s[1:]
    return s


def faixa(minimo: float | int, maximo: float | int) -> str:
    """Formata intervalo numérico separado por meia-risca."""
    return f"{num_curto(minimo)}–{num_curto(maximo)}"


def coord(ini: int, fim: int) -> str:
    """Formata intervalo semiaberto [ini, fim) base 0 como coordenadas base 1 inclusivas (ini+1–fim)."""
    return f"{ini + 1}–{fim}"


def pct(x: float | int, total: float | int) -> str:
    """Calcula porcentagem formatada com três casas e ponto decimal para estilo CSS."""
    if total == 0:
        return "0.000"
    return f"{float(x) / float(total) * 100:.3f}"


def data_hora(momento: datetime) -> str:
    """Formata data e hora no padrão dia/mês/ano hora:minuto."""
    return momento.strftime("%d/%m/%Y %H:%M")


def plural(n: int, singular: str, plural: str) -> str:
    """Retorna número seguido da forma singular se 1, ou plural caso contrário."""
    return f"{n} {singular}" if n == 1 else f"{n} {plural}"


def juntar(nomes: list[str]) -> str:
    """Junta lista de termos com vírgula e 'e' antes do último."""
    if not nomes:
        return ""
    if len(nomes) == 1:
        return nomes[0]
    if len(nomes) == 2:
        return f"{nomes[0]} e {nomes[1]}"
    return f"{', '.join(nomes[:-1])} e {nomes[-1]}"
