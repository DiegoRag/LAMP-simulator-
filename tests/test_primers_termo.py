"""Testes unitários para cálculos termodinâmicos de primers LAMP."""

import tomllib
from pathlib import Path

import pytest

from primers.termo import (
    Condicoes,
    Estrutura,
    dg_extremidade,
    gc_pct,
    grampo,
    heterodimero,
    homodimero,
    tm,
)


def test_gc_pct():
    assert gc_pct("GGCC") == 100
    assert gc_pct("ATAT") == 0
    assert gc_pct("ACGT") == 50


def test_tm_monotonia():
    cond = Condicoes()
    tm_gc = tm("GCGCGGCCGCGCGGCCGCGC", cond, conc_nM=200)
    tm_at = tm("ATATATTAATATATTAATAT", cond, conc_nM=200)
    assert tm_gc > tm_at


def test_estruturas_secundarias():
    cond = Condicoes()
    s1 = "ACGTACGTGGCCAATTGCAT"
    s2 = "ATGCAATTGGCCACGTACGT"

    g = grampo(s1, cond)
    assert isinstance(g, Estrutura)
    assert -20 < g.dg_kcal < 5

    h = homodimero(s1, cond)
    assert isinstance(h, Estrutura)
    assert -20 < h.dg_kcal < 5

    # s2 é o complemento reverso de s1: duplex perfeito de 20 pares. A 65 °C (temperatura da
    # reação) o ΔG fica perto de -11,3 kcal/mol; em cal/mol seria da ordem de -11000.
    het = heterodimero(s1, s2, cond)
    assert isinstance(het, Estrutura)
    assert het.encontrada
    assert -20 < het.dg_kcal < -5


def test_dg_extremidade():
    cond = Condicoes()
    # convenção do PrimerExplorer V5: ΔG°37 em NaCl 1 M com iniciação (+0,98 G·C, +1,03 A·T)
    assert dg_extremidade("TTTTGGGGGG", "3", 6, cond) == pytest.approx(-7.24, abs=0.01)
    assert dg_extremidade("GGGGGGTTTT", "5", 6, cond) == pytest.approx(-7.24, abs=0.01)
    assert dg_extremidade("AAAAAA", "3", 6, cond) == pytest.approx(-2.94, abs=0.01)
    assert dg_extremidade("CCACGTAC", "3", 6, cond) == pytest.approx(-5.06, abs=0.01)

    with pytest.raises(ValueError):
        dg_extremidade("CCACGTAC", "x", 6, cond)


def test_condicoes_de_config():
    cfg = {"reacao": {"na_mM": 10, "k_mM": 40, "nh4_mM": 5, "mg_mM": 6, "dntp_cada_mM": 1.0, "temp_C": 60}}
    cond = Condicoes.de_config(cfg)
    assert cond.na_mM == 10
    assert cond.k_mM == 40
    assert cond.nh4_mM == 5
    assert cond.mg_mM == 6
    assert cond.dntp_cada_mM == 1.0
    assert cond.temp_C == 60
    assert cond.monovalentes_mM == 55
    assert cond.dntp_total_mM == 4.0
    assert cond.mg_livre_mM == 2.0


def test_config_toml():
    caminho = Path("config/primers.toml")
    if not caminho.exists():
        caminho = Path(__file__).resolve().parent.parent / "config" / "primers.toml"

    with open(caminho, "rb") as f:
        cfg = tomllib.load(f)

    assert "reacao" in cfg
    assert "tm_c" in cfg
    assert "distancias" in cfg
    assert "desenho" in cfg
    assert cfg["extremidade"]["dg_max"] == -4.0
    assert cfg["reacao"]["dntp_cada_mM"] == 1.4
    assert "dntp_mM" not in cfg["reacao"]
    assert cfg["reacao"]["temp_C"] == 65.0
    assert cfg["tm_referencia"] == {"oligo_nM": 100.0, "na_mM": 50.0, "mg_mM": 4.0}
    assert cfg["concentracao_nM"] == {"externos": 200.0, "internos": 1600.0, "loops": 400.0}
