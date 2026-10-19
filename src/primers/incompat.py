"""Incompatibilidades de base entre um conjunto LAMP e o alvo, região por região.

É o que a análise entrega para a ponte com a simulação (esteira.ponte): quantas bases
de cada região não casam com o alvo e quantas dessas caem perto da ponta que estende.

Ponta que estende, por região (coordenadas da fita escolhida, base zero):
- F3, F2, LF, B3c, B2c, LB: a ponta 3' do próprio primer (casamento.distancia_3p).
- F1 e B1c: no haltere do LAMP, a fita se dobra e estende a partir da cópia de F1c
  (B1c) que o próprio primer trouxe; a ponta 3' dessa fita é complementar à ponta 5'
  de F1c (B1c). Por isso, para F1 e B1c conta a distância até a ponta 5' do primer,
  a mesma ponta cuja estabilidade a análise mede (extremidade 5').
A regra mora em casamento.distancia_ponta_extensao (fase 6), reexportada aqui.
"""

from dataclasses import dataclass

from primers.casamento import distancia_ponta_extensao, regioes_do_conjunto
from primers.modelo import ConjuntoLAMP
from primers.restricoes import Localizacao

__all__ = [
    "ESSENCIAIS",
    "LACOS",
    "Incompat",
    "IncompatRegiao",
    "distancia_ponta_extensao",
    "medir_incompat",
    "regioes_para_json",
]

ESSENCIAIS: tuple[str, ...] = ("F3", "F2", "F1", "B1c", "B2c", "B3c")
LACOS: tuple[str, ...] = ("LF", "LB")


@dataclass(frozen=True)
class IncompatRegiao:
    """Incompatibilidades de uma região do conjunto contra o alvo."""

    regiao: str  # "F3", "F2", "F1", "B1c", "B2c", "B3c", "LF" ou "LB"
    encontrada: bool  # False: nem exata nem com até [casamento] max_erros substituições
    total: int | None  # bases que não casam; None quando não encontrada
    ponta3: int | None  # dessas, quantas a menos de bases_ponta nt da ponta que estende


@dataclass(frozen=True)
class Incompat:
    """Incompatibilidades do conjunto inteiro contra o alvo."""

    por_regiao: tuple[IncompatRegiao, ...]  # na ordem de casamento.regioes_do_conjunto
    total: int | None  # soma das regiões; None se alguma região não foi encontrada
    ponta3: int | None  # idem, só as da ponta que estende
    nao_encontradas: tuple[str, ...]  # regiões sem posição no alvo


def medir_incompat(conjunto: ConjuntoLAMP, loc: Localizacao, bases_ponta: int) -> Incompat:
    """Conta as incompatibilidades de cada região a partir da localização no alvo."""
    if bases_ponta < 1:
        raise ValueError(f"bases_ponta deve ser >= 1, veio {bases_ponta}")
    itens: list[IncompatRegiao] = []
    for nome, _ in regioes_do_conjunto(conjunto):
        if nome in loc.regioes:
            itens.append(IncompatRegiao(nome, True, 0, 0))
        elif nome in loc.aproximadas:
            pos = loc.aproximadas[nome]
            perto = sum(
                1
                for d in pos.diferencas
                if distancia_ponta_extensao(nome, pos.regiao, d.posicao) < bases_ponta
            )
            itens.append(IncompatRegiao(nome, True, len(pos.diferencas), perto))
        else:
            itens.append(IncompatRegiao(nome, False, None, None))
    faltam = tuple(i.regiao for i in itens if not i.encontrada)
    if faltam:
        return Incompat(tuple(itens), None, None, faltam)
    total = sum(i.total or 0 for i in itens)
    ponta3 = sum(i.ponta3 or 0 for i in itens)
    return Incompat(tuple(itens), total, ponta3, ())


def regioes_para_json(incompat: Incompat) -> list[dict]:
    """Regiões em dicionários JSON: regiao, encontrada, total e ponta (a ponta que estende)."""
    return [
        {"regiao": i.regiao, "encontrada": i.encontrada, "total": i.total, "ponta": i.ponta3}
        for i in incompat.por_regiao
    ]
