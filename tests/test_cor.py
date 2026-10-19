import numpy as np
import pytest

from lampsim.cor import LAMBDA_NM, cmf_cie1931, cor_do_indicador, transmitancia_para_rgb

ROTAS = ["cie", "3bandas"]


def test_cmf_y_tem_pico_perto_de_555nm():
    y = cmf_cie1931()[1]
    assert 545 <= LAMBDA_NM[np.argmax(y)] <= 565
    assert y.max() == pytest.approx(1.0, abs=0.02)


@pytest.mark.parametrize("rota", ROTAS)
def test_tubo_transparente_e_branco(rota):
    rgb = transmitancia_para_rgb(np.ones((1, LAMBDA_NM.size)), rota=rota)
    np.testing.assert_allclose(rgb, 1.0, atol=1e-6)


@pytest.mark.parametrize("rota", ROTAS)
def test_forma_acida_e_amarela(rota):
    r, g, b = cor_do_indicador(1.0, rota=rota)[0]
    assert r > b + 0.25 and g > b + 0.25


@pytest.mark.parametrize("rota", ROTAS)
def test_forma_basica_e_rosa(rota):
    r, g, b = cor_do_indicador(0.0, rota=rota)[0]
    assert r > g and b > g  # vermelho e azul passam, o verde é absorvido


@pytest.mark.parametrize("rota", ROTAS)
def test_cor_muda_num_so_sentido(rota):
    # da forma básica para a ácida: o azul só cai e o verde só sobe
    rgb = cor_do_indicador(np.linspace(0, 1, 41), rota=rota)
    assert np.all(np.diff(rgb[:, 2]) <= 1e-9)
    assert np.all(np.diff(rgb[:, 1]) >= -1e-9)


@pytest.mark.parametrize("rota", ROTAS)
def test_valores_entre_zero_e_um(rota):
    rgb = cor_do_indicador(np.linspace(0, 1, 11), rota=rota)
    assert rgb.min() >= 0.0 and rgb.max() <= 1.0


def test_rota_invalida():
    with pytest.raises(ValueError):
        cor_do_indicador(0.5, rota="xyz")
