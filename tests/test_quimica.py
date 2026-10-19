import numpy as np
import pytest

from lampsim.quimica import fracao_acida, fracao_basica, ph_no_tempo, progresso_logistico

PKA = 8.0


def test_metade_acida_no_pka():
    assert fracao_acida(PKA, PKA) == pytest.approx(0.5)


def test_fracoes_somam_um():
    ph = np.linspace(5, 10, 51)
    np.testing.assert_allclose(fracao_acida(ph, PKA) + fracao_basica(ph, PKA), 1.0)


def test_mais_acido_quanto_menor_o_ph():
    ph = np.linspace(5, 10, 51)
    assert np.all(np.diff(fracao_acida(ph, PKA)) < 0)


def test_uma_unidade_de_ph_e_fator_dez():
    # Henderson-Hasselbalch: um pH abaixo do pKa -> razão ácida/básica = 10
    f = fracao_acida(PKA - 1, PKA)
    assert f / (1 - f) == pytest.approx(10.0)


def test_progresso_vai_de_zero_ate_A():
    # fase 5: a logística é normalizada para partir de zero exato em t = 0
    t = np.array([0.0, 8.0, 60.0])
    p = progresso_logistico(t, t0=8.0, k=1.0, A=0.9)
    assert p[0] == 0.0
    # 0,9 · (0,5 − L0) / (1 − L0), com L0 = 1 / (1 + e^8)
    assert p[1] == pytest.approx(0.449849, abs=1e-6)
    assert p[2] == pytest.approx(0.9)


def test_ph_cai_e_nunca_sobe():
    t = np.linspace(0, 30, 61)
    ph = ph_no_tempo(t, ph_inicial=8.7, queda_max=2.4, t0=7.0, k=0.9, A=1.0, deriva=0.001)
    assert np.all(np.diff(ph) <= 0)
    assert ph[-1] == pytest.approx(8.7 - 2.4 - 0.03, abs=1e-3)


def test_sem_reacao_e_sem_deriva_ph_constante():
    t = np.linspace(0, 30, 61)
    ph = ph_no_tempo(t, ph_inicial=8.7, queda_max=2.4, t0=7.0, k=0.9, A=0.0)
    np.testing.assert_allclose(ph, 8.7)
