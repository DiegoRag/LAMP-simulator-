"""Fase 5, sessão 2: logística normalizada (lampsim.quimica) e cinética por cópias e eficiência."""

import math

import numpy as np
import pytest

from lampsim.cinetica import parametros_cineticos, sortear_copias, tempo_virada
from lampsim.config import CINETICA
from lampsim.quimica import ph_no_tempo, progresso_logistico


def test_curvas_partem_do_mesmo_ph():
    # antes: um positivo com t0 cedo começava mais ácido que um negativo em t = 0
    cedo = ph_no_tempo(0.0, 8.7, 2.4, t0=3.0, k=1.0, A=1.0)
    negativo = ph_no_tempo(0.0, 8.7, 2.4, t0=3.0, k=1.0, A=0.0)
    assert cedo == negativo == 8.7


def test_progresso_sem_estouro_numerico():
    with np.errstate(over="raise"):
        p = progresso_logistico(np.array([0.0, 30.0]), t0=999.0, k=1.2, A=1.0)
    np.testing.assert_allclose(p, 0.0, atol=1e-12)
    # k = 0 (eficiência zero) não divide por zero
    np.testing.assert_allclose(progresso_logistico(np.array([0.0, 30.0]), t0=5.0, k=0.0, A=1.0), 0.0)


def test_tempo_virada_pelas_formulas():
    # t0 = 2 + (9 − log10 N) · 1,5 / eficiência
    assert tempo_virada(1.0, 1e5) == pytest.approx(8.0)
    assert tempo_virada(1.0, 10) == pytest.approx(14.0)
    assert tempo_virada(0.5, 1e5) == pytest.approx(14.0)
    assert tempo_virada(1.0, 1e10) == pytest.approx(2.0)  # acima da referência não acelera


def test_menos_copias_ou_menos_eficiencia_atrasam():
    copias = [1e6, 1e5, 1e4, 1e3, 1e2, 10, 1]
    t = [tempo_virada(1.0, n) for n in copias]
    assert all(b > a for a, b in zip(t, t[1:], strict=False))
    efs = [1.0, 0.8, 0.5, 0.3]
    t = [tempo_virada(e, 1e4) for e in efs]
    assert all(b > a for a, b in zip(t, t[1:], strict=False))


def test_conjunto_com_seis_erros_na_ponta_nao_vira():
    # eficiência 0,3^6 (seis incompatibilidades na ponta 3′): t0 bate no teto
    assert tempo_virada(0.3 ** 6, 1e5) == CINETICA["t0_max_min"]
    assert tempo_virada(0.0, 1e5) == CINETICA["t0_max_min"]
    assert tempo_virada(1.0, 0) == CINETICA["t0_max_min"]


def test_parametros_cineticos_sem_amplificacao():
    rng = np.random.default_rng(1)
    p = parametros_cineticos(1.0, 0, rng)
    assert p["A"] == 0.0 and p["t0"] == CINETICA["t0_max_min"]


def test_parametros_cineticos_reprodutiveis_e_nas_faixas():
    a = parametros_cineticos(0.8, 1000, np.random.default_rng(42))
    b = parametros_cineticos(0.8, 1000, np.random.default_rng(42))
    assert a == b
    lo, hi = CINETICA["k_por_min"]
    assert lo * 0.8 <= a["k"] <= hi * 0.8
    assert CINETICA["A"][0] <= a["A"] <= CINETICA["A"][1]
    assert abs(a["t0"] - tempo_virada(0.8, 1000)) < 5 * CINETICA["sd_t0_min"]


def test_sortear_copias_poisson():
    rng = np.random.default_rng(3)
    assert sortear_copias(0, rng) == 0
    n = [sortear_copias(1.0, np.random.default_rng(s)) for s in range(2000)]
    # P(0) = e^-1 ≈ 0,368 para média 1
    assert abs(n.count(0) / 2000 - math.exp(-1)) < 0.04
    grande = sortear_copias(1e5, np.random.default_rng(5))
    assert abs(grande - 1e5) < 5 * math.sqrt(1e5)


@pytest.mark.parametrize("ef, n", [(-0.1, 10), (1.5, 10), (float("nan"), 10), (0.5, -1)])
def test_entradas_invalidas(ef, n):
    with pytest.raises(ValueError):
        tempo_virada(ef, n)
