"""Cinética da amplificação: do número de cópias e da eficiência ao tempo da virada.

Modelo (PROVISÓRIO, constantes em config.CINETICA):
    t0 = atraso + max(0, log10_copias_ref - log10(cópias)) · min_por_log / eficiência
    k  = k_ref · eficiência
t0 é o meio da queda de pH (minutos) e k a velocidade da virada (1/min). Menos cópias
ou eficiência menor atrasam a virada; sem cópia nenhuma ou com eficiência zero, o tubo
não amplifica (A = 0). O número de cópias num tubo é sorteado por Poisson em volta da
média: com média de 1 cópia, cerca de 37 % dos tubos não recebem molde nenhum.
"""

import math

import numpy as np

from .config import CINETICA


def _numero(valor, nome):
    try:
        x = float(valor)
    except (TypeError, ValueError):
        raise ValueError(f"{nome} deve ser um número, veio {valor!r}") from None
    if x != x:  # NaN
        raise ValueError(f"{nome} não pode ser NaN")
    return x


def _validar(eficiencia, copias):
    if not 0.0 <= _numero(eficiencia, "eficiência") <= 1.0:
        raise ValueError(f"eficiência deve estar entre 0 e 1, veio {eficiencia!r}")
    if _numero(copias, "número de cópias") < 0.0:
        raise ValueError(f"número de cópias deve ser >= 0, veio {copias!r}")


def tempo_virada(eficiencia, copias, cinetica=CINETICA):
    """t0 determinístico (minutos), sem ruído; t0_max_min quando não há amplificação."""
    _validar(eficiencia, copias)
    teto = float(cinetica["t0_max_min"])
    if eficiencia <= 0.0 or copias < 1:
        return teto
    decadas = max(0.0, float(cinetica["log10_copias_ref"]) - math.log10(copias))
    t0 = float(cinetica["atraso_min"]) + decadas * float(cinetica["min_por_log"]) / float(eficiencia)
    return min(t0, teto)


def sortear_copias(media, rng):
    """Cópias de molde que caem num tubo: Poisson com a média dada (0 se a média é 0)."""
    if _numero(media, "média de cópias") < 0.0:
        raise ValueError(f"média de cópias deve ser >= 0, veio {media!r}")
    if media <= 0:
        return 0
    return int(rng.poisson(float(media)))


def parametros_cineticos(eficiencia, copias, rng, cinetica=CINETICA):
    """Sorteia A, t0 e k de um tubo. Sempre consome os mesmos 3 números do rng, nesta
    ordem: k_ref (uniforme), A (uniforme), ruído de t0 (normal)."""
    _validar(eficiencia, copias)
    k_ref = float(rng.uniform(*cinetica["k_por_min"]))
    a = float(rng.uniform(*cinetica["A"]))
    ruido = float(rng.normal(0.0, float(cinetica["sd_t0_min"])))
    teto = float(cinetica["t0_max_min"])
    k = k_ref * float(eficiencia)
    if eficiencia <= 0.0 or copias < 1:
        return {"A": 0.0, "t0": teto, "k": k}
    t_det = tempo_virada(eficiencia, copias, cinetica)
    if t_det >= teto:  # nem com o ruído a virada cai na janela: fica no teto, sem sorteio
        return {"A": a, "t0": teto, "k": k}
    return {"A": a, "t0": float(np.clip(t_det + ruido, 0.0, teto)), "k": k}
