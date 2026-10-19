"""Fase 3, sessão 3: condições da reação, Tm de referência e da reação, dobramento, ponta e estruturas."""

import itertools
import math
from dataclasses import replace

import primer3
import pytest

import primers.nn as nn
import primers.termo as termo
from primers.termo import (
    Condicoes,
    CondicoesReferencia,
    dg_extremidade,
    dg_laco_grampo,
    dimero_ponta3,
    grampo,
    heterodimero,
    tm,
    tm_dobramento,
    tm_referencia,
)

COND = Condicoes()
REF = CondicoesReferencia()
SEQ = "ACGTTGCAGCTAGCTTGCAG"
F1C = "ACCCCCAAGTGTTATTGGCAAGTT"


def test_condicoes_dntp_por_nucleotideo_e_mg_livre():
    reacao = {"na_mM": 0.0, "k_mM": 50.0, "nh4_mM": 20.0, "mg_mM": 8.0, "dntp_cada_mM": 1.4, "temp_C": 65.0}
    cfg = {"reacao": reacao}
    cond = Condicoes.de_config(cfg)
    assert cond == COND
    assert cond.monovalentes_mM == pytest.approx(70.0)
    assert cond.dntp_total_mM == pytest.approx(5.6)
    assert cond.mg_livre_mM == pytest.approx(2.4)
    assert replace(cond, mg_mM=5.0).mg_livre_mM == 0.0

    antigo = {"reacao": {"na_mM": 50.0, "mg_mM": 8.0, "dntp_mM": 1.4, "temp_C": 37.0}}
    with pytest.raises(ValueError, match="dntp_cada_mM"):
        Condicoes.de_config(antigo)


def test_thal_recebe_dntp_total_monovalentes_e_temperatura_da_reacao(monkeypatch):
    chamadas = []
    verdadeiro = primer3.calc_hairpin

    def falso(*args, **kwargs):
        chamadas.append(kwargs)
        return verdadeiro("ACGT")

    for nome in ("calc_hairpin", "calc_homodimer", "calc_heterodimer", "calc_end_stability"):
        monkeypatch.setattr(termo.primer3, nome, falso)
    grampo(SEQ, COND)
    termo.homodimero(SEQ, COND)
    heterodimero(SEQ, F1C, COND)
    dimero_ponta3(SEQ, F1C, COND)
    assert len(chamadas) == 5
    for kw in chamadas:
        assert kw["mv_conc"] == pytest.approx(70.0)
        assert kw["dv_conc"] == pytest.approx(8.0)
        assert kw["dntp_conc"] == pytest.approx(5.6)
        assert kw["temp_c"] == pytest.approx(65.0)


def test_tm_referencia_convencao_primerexplorer():
    assert tm_referencia(SEQ, REF) == pytest.approx(63.0682, abs=1e-3)
    assert tm_referencia("GCGAATTCGC", REF) == pytest.approx(36.4258, abs=1e-3)
    assert tm_referencia(SEQ, replace(REF, oligo_nM=1000.0)) == pytest.approx(66.4979, abs=1e-3)
    # o Mg entra como Na⁺ equivalente [Na⁺] + 4·√[Mg²⁺] no termo 16,6·log10
    delta = tm_referencia(SEQ, replace(REF, mg_mM=9.0)) - tm_referencia(SEQ, REF)
    esperado = 16.6 * math.log10((0.05 + 4 * math.sqrt(0.009)) / (0.05 + 4 * math.sqrt(0.004)))
    assert delta == pytest.approx(esperado, abs=1e-9)
    with pytest.raises(ValueError):
        tm_referencia("A", REF)


def test_tm_da_reacao_owczarzy_2008():
    assert tm(SEQ, COND, 1600.0) == pytest.approx(70.2726, abs=1e-3)
    assert tm(SEQ, COND, 200.0) == pytest.approx(67.2574, abs=1e-3)
    # 1,4 mM lido como total deixaria 6,6 mM de Mg²⁺ livre e a Tm alta
    errado = replace(COND, dntp_cada_mM=0.35)
    assert errado.mg_livre_mM == pytest.approx(6.6)
    assert tm(SEQ, errado, 1600.0) == pytest.approx(71.6633, abs=1e-3)
    # sem Mg²⁺ livre vale só a correção de monovalentes (Owczarzy 2004)
    assert tm(SEQ, Condicoes(mg_mM=0.0, dntp_cada_mM=0.0), 1600.0) == pytest.approx(65.0085, abs=1e-3)
    assert tm(SEQ, replace(COND, mg_mM=5.0), 1600.0) == pytest.approx(65.0085, abs=1e-3)
    # o primer3 2.3.1 perde o termo de comprimento da eq. 16 (divisão inteira)
    p3 = primer3.calc_tm(
        SEQ, mv_conc=70, dv_conc=8, dntp_conc=5.6, dna_conc=1600, salt_corrections_method="owczarzy"
    )
    assert p3 == pytest.approx(67.0221, abs=1e-3)
    # na convenção C/4 do primer3 (1600/4 = 400 nM), a eq. 16 completa dá 1,23 °C a mais
    assert tm(SEQ, COND, 400.0) - p3 == pytest.approx(1.2345, abs=1e-3)


def test_laco_e_dobramento_sem_termo_de_concentracao():
    assert dg_laco_grampo(3) == pytest.approx(3.5)
    assert dg_laco_grampo(11) == pytest.approx(4.8)
    assert dg_laco_grampo(30) == pytest.approx(6.3)
    assert dg_laco_grampo(50) == pytest.approx(7.0681, abs=1e-4)
    with pytest.raises(ValueError):
        dg_laco_grampo(2)
    assert tm_dobramento(F1C, 50, COND) == pytest.approx(76.3734, abs=1e-3)
    assert tm_dobramento(F1C, 40, COND) == pytest.approx(77.0820, abs=1e-3)
    assert tm_dobramento(F1C, 60, COND) == pytest.approx(75.7965, abs=1e-3)
    assert tm(F1C, COND, 1600.0) == pytest.approx(70.8526, abs=1e-3)


def test_ponta_na_convencao_primerexplorer():
    assert dg_extremidade("AAAAAA", "3", 6, REF) == pytest.approx(-2.94)
    assert dg_extremidade("CCCCGGGGGG", "3", 6, REF) == pytest.approx(-7.24)
    assert dg_extremidade("ACGCGTTTTT", "5", 6, REF) == pytest.approx(-6.97)
    seis = ["".join(p) for p in itertools.product("ACGT", repeat=6)]
    passam = sum(1 for s in seis if dg_extremidade(s, "3", 6, REF) <= -4.0)
    assert passam == 3193


def test_estruturas_avaliadas_na_temperatura_da_reacao():
    s = "TATAGGACTTTGTCCTTATTA"
    a37 = grampo(s, replace(COND, temp_C=37.0))
    a65 = grampo(s, COND)
    assert a37.dg_kcal == pytest.approx(-2.027, abs=0.01)
    assert a65.dg_kcal == pytest.approx(0.255, abs=0.01)


def test_dimero_com_ponta_3_extensivel():
    a = "ATATAGACCGTGCAGGCCATCAT"
    b_meio = "TCTCTGCCTGCACGGTCTCTCT"
    b_ponta = "TCTCTATGATGGCCTGCACGGTC"
    assert heterodimero(a, b_meio, COND).dg_kcal == pytest.approx(-9.51, abs=0.01)
    assert dimero_ponta3(a, b_meio, COND).dg_kcal == pytest.approx(-4.66, abs=0.01)
    assert dimero_ponta3(a, b_ponta, COND).dg_kcal == pytest.approx(-11.53, abs=0.01)
    assert heterodimero(a, b_ponta, COND).dg_kcal == pytest.approx(-11.53, abs=0.01)


def test_nn_tabelas_e_formulas():
    assert nn.DG37_NN["GC"] == -2.24
    assert nn.dg37_ponta("GGGGGG") == pytest.approx(-7.24)
    assert nn.dh_ds_vizinhos("GCGAATTCGC") == pytest.approx((-80.2, -212.4))
    # duas iniciações G·C (+0,1 / −2,8 cada) e a simetria (−1,4) do palíndromo
    assert nn.dh_ds_duplex("GCGAATTCGC") == pytest.approx((-80.0, -219.4))
    mono = (4.29 * 0.5 - 3.95) * 1e-5 * math.log(0.07) + 9.40e-6 * math.log(0.07) ** 2
    assert nn.correcao_owczarzy(0.5, 20, 70.0, 0.0) == pytest.approx(mono)
    assert nn.correcao_owczarzy(0.5, 20, 70.0, 2.4) == pytest.approx(6.5591e-05, rel=1e-4)
