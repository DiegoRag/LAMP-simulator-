"""Fase 5, sessão 3: diagnóstico pela curva, papel do tubo, eficiência e cópias, hash do lote."""

import numpy as np
import pytest

from lampsim import config
from lampsim.classes import (
    DIAGNOSTICOS,
    amostra_de_reacao,
    diagnosticar,
    simular,
    sortear_amostra,
)
from lampsim.manifesto import gerar_lote, hash_lote, parametros_do_simulador

T60 = np.arange(0.0, 60.01, 0.5)


@pytest.mark.parametrize(
    "f0, fl, papel, esperado",
    [
        (0.1, 0.9, "amostra", "positivo"),
        (0.1, 0.2, "amostra", "negativo"),
        (0.1, 0.5, "amostra", "positivo"),  # o limiar conta como positivo
        (0.1, 0.9, "controle_negativo", "invalido"),
        (0.1, 0.2, "controle_negativo", "controle_ok"),
        (0.1, 0.9, "controle_positivo", "controle_ok"),
        (0.1, 0.2, "controle_positivo", "invalido"),
        (0.8, 0.9, "amostra", "invalido"),  # já começa ácido: falha de reação
    ],
)
def test_tabela_do_diagnostico(f0, fl, papel, esperado):
    assert diagnosticar(f0, fl, papel) == esperado
    assert esperado in DIAGNOSTICOS


def test_negativo_e_controle_negativo_tem_a_mesma_optica():
    # antes: duas classes com as mesmas faixas e rótulos diferentes (erro mínimo de 50 %)
    a = sortear_amostra("negativo", 5)
    b = sortear_amostra("negativo", 5, papel="controle_negativo")
    np.testing.assert_array_equal(simular(a)["rgb"], simular(b)["rgb"])
    assert (a.diagnostico, b.diagnostico) == ("negativo", "controle_ok")


def test_controle_contaminado_que_ainda_nao_virou_e_controle_ok():
    # antes: aos 15 min metade dos controle_contaminado era igual a um negativo e saía "invalido"
    leitura = int(config.TEMPO["leitura_nominal_min"] / config.TEMPO["passo_min"])
    contaminados = [a for a in gerar_lote(100, 2026) if a.id.startswith("controle_contaminado")]
    vistos = set()
    for a in contaminados:
        f = simular(a)["fracao_acida"][leitura]
        assert a.diagnostico == ("invalido" if f >= 0.5 else "controle_ok")
        vistos.add(a.diagnostico)
    assert vistos == {"invalido", "controle_ok"}


def test_conjunto_perfeito_e_com_seis_erros_dao_curvas_diferentes():
    # antes: as duas curvas eram estatisticamente iguais
    perfeito = amostra_de_reacao(1.0, 1e5, 7)
    ruim = amostra_de_reacao(0.3 ** 6, 1e5, 7)
    assert perfeito.diagnostico == "positivo" and ruim.diagnostico == "negativo"
    assert perfeito.t0 == pytest.approx(8.404, abs=1e-3)  # 2 + (9 − log10 99 731) · 1,5 + ruído
    assert ruim.t0 == config.CINETICA["t0_max_min"]
    assert simular(ruim, t=T60)["fracao_acida"][-1] < 0.4
    assert simular(perfeito, t=T60)["fracao_acida"][-1] > 0.9


def test_sem_molde_nao_amplifica():
    a = amostra_de_reacao(1.0, 0, 7)
    assert a.copias == 0 and a.A == 0.0 and a.diagnostico == "negativo"
    s = simular(a, t=T60)
    assert np.all(np.abs(s["ph"] - (a.ph_inicial - a.deriva * T60)) < 1e-12)


def test_amostra_de_reacao_reprodutivel_e_com_papel():
    assert amostra_de_reacao(0.7, 100, 11) == amostra_de_reacao(0.7, 100, 11)
    c = amostra_de_reacao(1.0, 1e5, 3, papel="controle_positivo", id_amostra="cp")
    assert (c.id, c.classe, c.papel, c.diagnostico) == ("cp", "reacao", "controle_positivo", "controle_ok")
    with pytest.raises(ValueError):
        amostra_de_reacao(1.2, 100, 1)


def test_curva_parte_do_ph_inicial():
    a = sortear_amostra("positivo_forte", 1)
    assert simular(a)["ph"][0] == a.ph_inicial


def test_turbidez_acompanha_o_produto():
    a = sortear_amostra("positivo_forte", 2)
    s = simular(a)
    assert s["progresso"][0] == 0.0 and s["progresso"][-1] > 0.85


def test_hash_muda_com_parametros_de_quimica_e_cor(monkeypatch):
    # antes: o hash só via as amostras; dois lotes de cores diferentes tinham o mesmo hash
    lote = gerar_lote(2, 1)
    antes = hash_lote(lote)
    monkeypatch.setitem(config.ESPECTRO, "turbidez_max", 0.06)
    assert hash_lote(lote) != antes
    monkeypatch.setitem(config.ESPECTRO, "turbidez_max", 0.05)
    assert hash_lote(lote) == antes
    monkeypatch.setitem(config.QUIMICA, "pka", 7.9)
    assert hash_lote(lote) != antes


def test_parametros_do_simulador_tem_tudo():
    p = parametros_do_simulador()
    assert set(p) == {"versao_parametros", "tempo", "quimica", "espectro", "cinetica",
                      "mix", "leitura", "classes", "lote"}
