"""Fase 3, sessão 5: análise com as duas Tm, F1c/B1c como dobramento e dímeros pela ponta 3′."""

from pathlib import Path

import pytest

from primers.analise import analisar
from primers.entrada import carregar_config, ler_alvo, ler_conjunto
from primers.termo import (
    Condicoes,
    CondicoesReferencia,
    dimero_ponta3,
    heterodimero,
    tm,
    tm_dobramento,
    tm_referencia,
)

DADOS = Path(__file__).parent / "dados"
CONFIG = Path(__file__).parents[1] / "config" / "primers.toml"


def _linha(resultado, nome):
    return next(linha for linha in resultado.linhas if linha["nome"] == nome)


def test_linhas_com_tm_de_referencia_e_tm_da_reacao():
    cfg = carregar_config(CONFIG)
    cond = Condicoes.de_config(cfg)
    ref = CondicoesReferencia.de_config(cfg)
    leitura = ler_conjunto(DADOS / "conjunto_ok.csv")
    res = analisar(leitura, ler_alvo(DADOS / "alvo_ok.fasta"), cfg)

    f2 = _linha(res, "F2")
    assert f2["tm_c"] == pytest.approx(tm_referencia("CGCATATGGATTCTGTAC", ref))
    assert f2["tm_c"] == pytest.approx(51.2695, abs=1e-3)
    assert f2["tm_reacao_c"] == pytest.approx(tm("CGCATATGGATTCTGTAC", cond, 1600.0))
    assert f2["tm_reacao_c"] == pytest.approx(59.0914, abs=1e-3)
    assert _linha(res, "F3")["tm_reacao_c"] == pytest.approx(58.7531, abs=1e-3)
    fip = _linha(res, "FIP")
    assert fip["tm_c"] is None
    assert fip["tm_reacao_c"] is None


def test_f1c_e_b1c_na_reacao_como_dobramento_com_laco_da_geometria():
    cfg = carregar_config(CONFIG)
    cond = Condicoes.de_config(cfg)
    leitura = ler_conjunto(DADOS / "conjunto_ok.csv")
    alvo = ler_alvo(DADOS / "alvo_ok.fasta")
    f1c = leitura.conjunto.F1c
    b1c = leitura.conjunto.B1c

    res = analisar(leitura, alvo, cfg)  # F2_F1 = 50, B1c_B2c = 50, sem ligante
    assert _linha(res, "F1c")["tm_reacao_c"] == pytest.approx(tm_dobramento(f1c, 50, cond))
    assert _linha(res, "F1c")["tm_reacao_c"] == pytest.approx(72.1688, abs=1e-3)
    assert _linha(res, "B1c")["tm_reacao_c"] == pytest.approx(tm_dobramento(b1c, 50, cond))
    assert _linha(res, "B1c")["tm_reacao_c"] == pytest.approx(71.9396, abs=1e-3)

    alvo_55 = alvo[:60] + "AAAAA" + alvo[60:]  # F2_F1 passa a 55
    res_55 = analisar(leitura, alvo_55, cfg)
    assert _linha(res_55, "F1c")["tm_reacao_c"] == pytest.approx(71.8412, abs=1e-3)
    assert _linha(res_55, "B1c")["tm_reacao_c"] == pytest.approx(71.9396, abs=1e-3)

    sem_alvo = analisar(leitura, None, cfg)  # distância presumida de [dobramento]
    assert _linha(sem_alvo, "F1c")["tm_reacao_c"] == pytest.approx(72.1688, abs=1e-3)


def test_heterodimeros_priorizam_a_ponta_3():
    cfg = carregar_config(CONFIG)
    cond = Condicoes.de_config(cfg)
    leitura = ler_conjunto(DADOS / "conjunto_ok.csv")
    res = analisar(leitura, None, cfg)
    het = res.heterodimeros
    assert het["nomes"] == ["F3", "B3", "FIP", "BIP"]
    fip = leitura.primers["FIP"]
    bip = leitura.primers["BIP"]
    b3 = leitura.primers["B3"]
    assert het["dg_kcal"][2][3] == pytest.approx(dimero_ponta3(fip, bip, cond).dg_kcal)
    assert het["dg_kcal"][2][3] == pytest.approx(-1.0339, abs=1e-3)
    assert het["dg_qualquer_kcal"][1][3] == pytest.approx(heterodimero(b3, bip, cond).dg_kcal)
    assert het["dg_qualquer_kcal"][1][3] == pytest.approx(-3.8475, abs=1e-3)
    assert het["dg_kcal"][1][3] == pytest.approx(1.5274, abs=1e-3)


def test_laco_com_geometria_fora_de_ordem_usa_a_distancia_presumida():
    cfg = carregar_config(CONFIG)
    cond = Condicoes.de_config(cfg)
    leitura = ler_conjunto(DADOS / "conjunto_ok.csv")
    alvo = ler_alvo(DADOS / "alvo_ok.fasta")
    alvo_f1_antes = alvo[81:103] + alvo[:81] + alvo[103:]  # F1 antes do F2: F2_F1 = -53
    res = analisar(leitura, alvo_f1_antes, cfg)
    assert ("F2_F1", "distancia") in {(v.item, v.regra) for v in res.violacoes}
    f1c = leitura.conjunto.F1c
    assert _linha(res, "F1c")["tm_reacao_c"] == pytest.approx(tm_dobramento(f1c, 50, cond))
