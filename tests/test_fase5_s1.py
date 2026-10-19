"""Fase 5, sessão 1: espectro do vermelho de fenol, turbidez e matriz da câmera (lampsim.cor)."""

import colorsys

import numpy as np
import pytest

from lampsim.config import ESPECTRO
from lampsim.cor import (
    LAMBDA_NM,
    _bandas,
    absorbancia,
    absorbancia_pico_acida,
    cor_do_indicador,
    matriz_camera_para_xyz,
    rgb_para_hex,
    transmitancia_para_rgb,
)
from lampsim.quimica import fracao_acida

ROTAS = ["cie", "3bandas"]


def _banda(forma, lam):
    return float(_bandas(ESPECTRO[forma], np.array([float(lam)]))[0])


def test_razoes_de_yao_e_byrne_2001():
    # e1, e2 e e3 do vermelho de fenol em 433 e 558 nm (Environ. Sci. Technol. 35:1197)
    acida_433 = _banda("acida", 433)
    assert _banda("acida", 558) / acida_433 == pytest.approx(0.00244, rel=0.02)
    assert _banda("basica", 558) / acida_433 == pytest.approx(2.734, rel=0.01)
    assert _banda("basica", 433) / acida_433 == pytest.approx(0.1075, rel=0.02)


def test_ponto_isosbestico_perto_de_478_nm():
    # antes: picos de mesma absorbância e banda inventada em 420 nm
    lam = np.arange(440.0, 521.0, 1.0)
    d = _bandas(ESPECTRO["basica"], lam) - _bandas(ESPECTRO["acida"], lam)
    cruzamentos = lam[np.where(np.diff(np.sign(d)))[0]]
    assert len(cruzamentos) == 1 and 474 <= cruzamentos[0] <= 482


def test_absorbancia_no_pico_e_epsilon_c_l():
    assert absorbancia_pico_acida() == pytest.approx(23000.0 * 100e-6 * 0.5)
    a = absorbancia(1.0)[0]
    i433 = int(np.where(LAMBDA_NM == 435.0)[0][0])
    assert a[i433] == pytest.approx(1.15 * _banda("acida", 435) + ESPECTRO["fundo"])


def _matiz_graus(rgb):
    return colorsys.rgb_to_hsv(*rgb)[0] * 360.0


def test_ph_7_3_nao_sai_amarelo():
    # antes: pH 7,3 saía #F2DD75 (amarelo, matiz 51°)
    rgb = cor_do_indicador(fracao_acida(7.3, 8.0))[0]
    assert _matiz_graus(rgb) < 35.0
    assert rgb[0] - rgb[1] > 0.25
    # pH 6,3, fim de uma reação positiva: amarelo
    assert _matiz_graus(cor_do_indicador(fracao_acida(6.3, 8.0))[0]) > 45.0


@pytest.mark.parametrize("rota", ROTAS)
def test_forma_basica_e_viva_nao_acinzentada(rota):
    # antes: #C78AB9, malva-acinzentado (saturação 0,31)
    r, g, b = cor_do_indicador(0.0, rota=rota)[0]
    assert colorsys.rgb_to_hsv(r, g, b)[1] > 0.8
    assert r > g + 0.6 and b > g + 0.5


def test_cores_conhecidas_da_rota_cie():
    hexas = rgb_para_hex(cor_do_indicador(fracao_acida(np.array([6.3, 7.3, 8.7]), 8.0)))
    assert hexas == ["#FFEB4C", "#FFA866", "#FE00B1"]


def test_turbidez_escurece_por_igual():
    limpo = cor_do_indicador(0.5)[0]
    turvo = cor_do_indicador(0.5, turbidez=0.05)[0]
    assert np.all(turvo < limpo)
    np.testing.assert_allclose(cor_do_indicador(0.5, turbidez=0.0), cor_do_indicador(0.5))
    # uma turbidez por tubo
    dois = cor_do_indicador([0.5, 0.5], turbidez=[0.0, 0.05])
    np.testing.assert_allclose(dois[0], limpo)
    np.testing.assert_allclose(dois[1], turvo)


def test_matriz_da_camera():
    m = matriz_camera_para_xyz()
    assert m.shape == (3, 3)
    # o branco continua branco depois da matriz
    np.testing.assert_allclose(
        transmitancia_para_rgb(np.ones((1, LAMBDA_NM.size)), rota="3bandas"), 1.0, atol=1e-6
    )
    # com a matriz, a rota da câmera fica perto da colorimetria CIE
    f = np.linspace(0.0, 1.0, 5)
    assert np.abs(cor_do_indicador(f, rota="cie") - cor_do_indicador(f, rota="3bandas")).max() < 0.2


def test_rgb_para_hex():
    assert rgb_para_hex([[1.0, 0.0, 0.5]]) == ["#FF0080"]
    assert rgb_para_hex(np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])) == ["#000000", "#FFFFFF"]
