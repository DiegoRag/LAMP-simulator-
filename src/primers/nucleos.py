"""Busca do desenho LAMP por núcleo (F2, B2c): todos os núcleos do alvo, monte limitado e diversidade.

Para cada F2 e cada B2c com fim(B2c) − ini(F2) na faixa F2_B2, o melhor par interno (F1, B1c) e o melhor
par externo (F3, B3c) saem independentes um do outro, porque nenhuma distância liga um par ao outro. Assim
cada núcleo custa duas buscas pequenas, e o alvo inteiro é percorrido sem teto de combinações.
"""

import heapq
from bisect import bisect_left, bisect_right
from collections.abc import Callable, Iterator
from dataclasses import dataclass

from primers.candidatos import TIPOS_BASE, Candidato
from primers.modelo import Regiao


@dataclass(frozen=True, slots=True)
class Nucleo:
    """Melhor conjunto completo de um núcleo (F2, B2c), com a pontuação da etapa 1."""

    pontuacao: float  # desvio_tm + pen_f2_b2
    chave: tuple[int, ...]  # ini de F3, F2, F1, B1c, B2c, B3c e depois os fim, nesta ordem
    desvio_tm: float  # soma dos seis desvio_tm, na ordem F3, F2, F1, B1c, B2c, B3c
    pen_f2_b2: float
    combo: tuple[Candidato, ...]  # (F3, F2, F1, B1c, B2c, B3c)


def sobreposicao(a: Regiao, b: Regiao) -> int:
    """Número de nucleotídeos em comum entre duas regiões (0 se disjuntas)."""
    return max(0, min(a.fim, b.fim) - max(a.ini, b.ini))


def mesmo_nucleo(a: Nucleo, b: Nucleo, limite_nt: int) -> bool:
    """True se F2 e B2c de a e b se sobrepõem, os dois ao mesmo tempo, em mais de limite_nt nt."""
    return (
        sobreposicao(a.combo[1].regiao, b.combo[1].regiao) > limite_nt
        and sobreposicao(a.combo[4].regiao, b.combo[4].regiao) > limite_nt
    )


def _no_intervalo(
    lista: list[Candidato], inis: list[int], lo: int, hi: int, pelo_fim: bool, l_min: int, l_max: int
) -> list[Candidato]:
    """Candidatos com ini (ou fim, se pelo_fim) em [lo, hi], na ordem (desvio_tm, ini, fim)."""
    if lo > hi:
        return []
    if pelo_fim:
        fatia = lista[bisect_left(inis, lo - l_max) : bisect_right(inis, hi - l_min)]
        sel = [c for c in fatia if lo <= c.regiao.fim <= hi]
    else:
        sel = lista[bisect_left(inis, lo) : bisect_right(inis, hi)]
    return sorted(sel, key=lambda c: (c.desvio_tm, c.regiao.ini, c.regiao.fim))


def _melhor_par(
    esq: list[Candidato], dir_: list[Candidato], ok: Callable[[Candidato, Candidato], bool]
) -> tuple[Candidato, Candidato] | None:
    """Par (a, b) de menor a.desvio_tm + b.desvio_tm com ok(a, b); listas na ordem (desvio_tm, ini, fim).

    Empate: fica o primeiro achado (a mais cedo na lista, depois b).
    """
    melhor: tuple[float, Candidato, Candidato] | None = None
    if not esq or not dir_:
        return None
    for a in esq:
        if melhor is not None and a.desvio_tm + dir_[0].desvio_tm >= melhor[0]:
            break
        for b in dir_:
            if melhor is not None and a.desvio_tm + b.desvio_tm >= melhor[0]:
                break
            if ok(a, b):
                melhor = (a.desvio_tm + b.desvio_tm, a, b)
                break
    return None if melhor is None else (melhor[1], melhor[2])


def percorrer_nucleos(cands: dict[str, list[Candidato]], cfg: dict) -> Iterator[Nucleo]:
    """Gera o melhor conjunto de cada núcleo (F2, B2c) válido, na ordem de F2 e depois de B2c.

    `cands[tipo]` vem de gerar_candidatos (ordenada por ini, fim). Distâncias de [distancias]:
    F3 com fim em [ini(F2) − F3_F2_gap[1], ini(F2) − F3_F2_gap[0]]; F1 com ini em ini(F2) + F2_F1;
    B2c com fim em ini(F2) + F2_B2; B1c com fim em fim(B2c) − B1c_B2c[1] … fim(B2c) − B1c_B2c[0];
    B3c com ini em fim(B2c) + B2c_B3c_gap; par interno com ini(B1c) − fim(F1) em F1_B1c_gap; par externo
    com fim(B3c) − ini(F3) em F3_B3_total (se a chave existir).
    """
    D = cfg["distancias"]
    l_min = int(cfg["comprimento"]["min"])
    l_max = int(cfg["comprimento"]["max"])
    pref = cfg["desenho"]["F2_B2_pref"]
    peso_f2_b2 = float(cfg["desenho"]["peso_f2_b2"])
    total = D.get("F3_B3_total")
    gap = D["F1_B1c_gap"]
    listas = {t: cands.get(t, []) for t in TIPOS_BASE}
    inis = {t: [c.regiao.ini for c in listas[t]] for t in TIPOS_BASE}

    def ok_interno(f1: Candidato, b1c: Candidato) -> bool:
        return gap[0] <= b1c.regiao.ini - f1.regiao.fim <= gap[1]

    def ok_externo(f3: Candidato, b3c: Candidato) -> bool:
        return total is None or total[0] <= b3c.regiao.fim - f3.regiao.ini <= total[1]

    lados_b2c: dict[Regiao, tuple[list[Candidato], list[Candidato]]] = {}
    for f2 in listas["F2"]:
        i2 = f2.regiao.ini
        f3s = _no_intervalo(
            listas["F3"], inis["F3"], i2 - D["F3_F2_gap"][1], i2 - D["F3_F2_gap"][0], True, l_min, l_max
        )
        f1s = _no_intervalo(
            listas["F1"], inis["F1"], i2 + D["F2_F1"][0], i2 + D["F2_F1"][1], False, l_min, l_max
        )
        if not f3s or not f1s:
            continue
        lo, hi = i2 + D["F2_B2"][0], i2 + D["F2_B2"][1]
        fatia = listas["B2c"][bisect_left(inis["B2c"], lo - l_max) : bisect_right(inis["B2c"], hi - l_min)]
        for b2c in fatia:
            f = b2c.regiao.fim
            if not lo <= f <= hi:
                continue
            if b2c.regiao not in lados_b2c:
                lados_b2c[b2c.regiao] = (
                    _no_intervalo(listas["B1c"], inis["B1c"], f - D["B1c_B2c"][1], f - D["B1c_B2c"][0],
                                  True, l_min, l_max),
                    _no_intervalo(listas["B3c"], inis["B3c"], f + D["B2c_B3c_gap"][0],
                                  f + D["B2c_B3c_gap"][1], False, l_min, l_max),
                )
            b1cs, b3cs = lados_b2c[b2c.regiao]
            interno = _melhor_par(f1s, b1cs, ok_interno)
            externo = _melhor_par(f3s, b3cs, ok_externo) if interno is not None else None
            if interno is None or externo is None:
                continue
            combo = (externo[0], f2, interno[0], interno[1], b2c, externo[1])
            desvio = 0.0
            for c in combo:
                desvio += c.desvio_tm
            dist = f - i2
            pen = peso_f2_b2 * max(0.0, float(pref[0]) - dist, dist - float(pref[1]))
            chave = tuple(c.regiao.ini for c in combo) + tuple(c.regiao.fim for c in combo)
            yield Nucleo(desvio + pen, chave, desvio, pen, combo)


def melhores_nucleos(cands: dict[str, list[Candidato]], cfg: dict) -> tuple[list[Nucleo], int]:
    """Os [desenho] nucleos_guardados melhores núcleos, na ordem (pontuacao, chave), e quantos havia."""
    contagem = 0

    def contar() -> Iterator[Nucleo]:
        nonlocal contagem
        for nucleo in percorrer_nucleos(cands, cfg):
            contagem += 1
            yield nucleo

    guardados = heapq.nsmallest(
        int(cfg["desenho"]["nucleos_guardados"]), contar(), key=lambda x: (x.pontuacao, x.chave)
    )
    return guardados, contagem


def escolher_distintos(
    nucleos: list[Nucleo], cfg: dict, aceita: Callable[[Nucleo], bool]
) -> list[Nucleo]:
    """Na ordem dada, guarda os núcleos que não são mesmo_nucleo de um já guardado e que aceita() aprova.

    Para em [desenho] pre_selecao núcleos. O teste de mesmo_nucleo vem antes de aceita (que pode custar
    caro, como o homodímero).
    """
    limite = int(cfg["desenho"]["sobreposicao_max_nt"])
    maximo = int(cfg["desenho"]["pre_selecao"])
    escolhidos: list[Nucleo] = []
    for nucleo in nucleos:
        if len(escolhidos) >= maximo:
            break
        if any(mesmo_nucleo(nucleo, outro, limite) for outro in escolhidos):
            continue
        if aceita(nucleo):
            escolhidos.append(nucleo)
    return escolhidos
