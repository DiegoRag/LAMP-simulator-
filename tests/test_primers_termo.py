"""Testes unitários para cálculos termodinâmicos de primers LAMP."""

from pathlib import Path
import tomllib
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

    # s2 é o complemento reverso de s1: duplex perfeito de 20 pares. O ΔG tem de ser
    # fortemente negativo e estar em kcal/mol — em cal/mol seria da ordem de -24000.
    het = heterodimero(s1, s2, cond)
    assert isinstance(het, Estrutura)
    assert het.encontrada
    assert -40 < het.dg_kcal < -15


def test_dg_extremidade():
    cond = Condicoes()
    assert dg_extremidade("TTTTGGGGGG", "3", 6, cond) == pytest.approx(-9.140, abs=0.01)
    assert dg_extremidade("GGGGGGTTTT", "5", 6, cond) == pytest.approx(-9.140, abs=0.01)
    assert dg_extremidade("AAAAAA", "3", 6, cond) == pytest.approx(-5.073, abs=0.01)
    assert dg_extremidade("CCACGTAC", "3", 6, cond) == pytest.approx(-7.116, abs=0.01)

    with pytest.raises(ValueError):
        dg_extremidade("CCACGTAC", "x", 6, cond)


def test_condicoes_de_config():
    cfg = {"reacao": {"na_mM": 40, "mg_mM": 6, "dntp_mM": 1.0, "temp_C": 25}}
    cond = Condicoes.de_config(cfg)
    assert cond.na_mM == 40
    assert cond.mg_mM == 6
    assert cond.dntp_mM == 1.0
    assert cond.temp_C == 25


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
