"""Testes para o módulo de formatação de texto."""

from datetime import datetime

from interface.formato import coord, data_hora, faixa, juntar, milhar, num, num_curto, pct, plural


def test_formato_num():
    assert num(68.75) == "68,75"
    assert num(-6.96) == "−6,96"
    assert num(None) == "—"
    assert num(40.17, 1) == "40,2"
    assert num(-0.001) == "0,00"
    assert num(5) == "5,00"


def test_formato_num_curto():
    assert num_curto(59.0) == "59"
    assert num_curto(0.5) == "0,5"
    assert num_curto(-4.0) == "−4"
    assert num_curto(66.64) == "66,64"
    assert num_curto(1.4 * 4) == "5,6"


def test_formato_milhar():
    assert milhar(6644) == "6.644"
    assert milhar(200000) == "200.000"
    assert milhar(5) == "5"
    assert milhar(5635218) == "5.635.218"


def test_formato_faixa_e_coord():
    assert faixa(59.0, 61.0) == "59–61"
    assert faixa(0, 1000) == "0–1000"
    assert coord(60, 76) == "61–76"


def test_formato_pct():
    assert pct(60, 320) == "18.750"
    assert pct(0, 320) == "0.000"
    assert pct(320, 320) == "100.000"


def test_formato_data_hora():
    assert data_hora(datetime(2026, 10, 16, 14, 3)) == "16/10/2026 14:03"


def test_formato_plural_e_juntar():
    assert plural(1, "violação", "violações") == "1 violação"
    assert plural(0, "violação", "violações") == "0 violações"
    assert juntar(["F1"]) == "F1"
    assert juntar(["F1", "B1c"]) == "F1 e B1c"
    assert juntar(["F1", "B1c", "LB"]) == "F1, B1c e LB"
