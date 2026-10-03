"""Testes para o módulo de restrições de primers LAMP."""

from pathlib import Path
from primers.entrada import carregar_config, ler_alvo, ler_conjunto
from primers.modelo import Regiao
from primers.restricoes import localizar, verificar

DADOS = Path(__file__).parent / "dados"
CONFIG = Path(__file__).parents[1] / "config" / "primers.toml"


def test_conjunto_ok_sem_alvo():
    cfg = carregar_config(CONFIG)
    leitura = ler_conjunto(DADOS / "conjunto_ok.csv")
    violacoes = verificar(leitura.conjunto, None, cfg)
    assert {(v.item, v.regra) for v in violacoes} == set()


def test_conjunto_ok_com_alvo_ok():
    cfg = carregar_config(CONFIG)
    leitura = ler_conjunto(DADOS / "conjunto_ok.csv")
    alvo = ler_alvo(DADOS / "alvo_ok.fasta")

    violacoes = verificar(leitura.conjunto, alvo, cfg)
    assert {(v.item, v.regra) for v in violacoes} == set()

    regioes, _ = localizar(leitura.conjunto, alvo)
    assert regioes["F3"] == Regiao(0, 21)
    assert regioes["F2"] == Regiao(31, 49)
    assert regioes["F1"] == Regiao(81, 103)
    assert regioes["B1c"] == Regiao(108, 127)
    assert regioes["B2c"] == Regiao(157, 177)
    assert regioes["B3c"] == Regiao(187, 206)


def test_conjunto_violacoes_sem_alvo():
    cfg = carregar_config(CONFIG)
    leitura = ler_conjunto(DADOS / "conjunto_violacoes.csv")
    violacoes = verificar(leitura.conjunto, None, cfg)

    esperado = {
        ("B1c", "extremidade"),
        ("B1c", "gc_pct"),
        ("B1c", "tm_c"),
        ("F3", "gc_pct"),
        ("F3", "tm_c"),
    }
    assert {(v.item, v.regra) for v in violacoes} == esperado


def test_conjunto_ok_com_alvo_gap():
    cfg = carregar_config(CONFIG)
    leitura = ler_conjunto(DADOS / "conjunto_ok.csv")
    alvo = ler_alvo(DADOS / "alvo_gap.fasta")

    violacoes = verificar(leitura.conjunto, alvo, cfg)
    assert {(v.item, v.regra) for v in violacoes} == {("F3_F2_gap", "distancia")}
    assert violacoes[0].valor == 70


def test_conjunto_ok_com_alvo_sem_f3():
    cfg = carregar_config(CONFIG)
    leitura = ler_conjunto(DADOS / "conjunto_ok.csv")
    alvo_ok = ler_alvo(DADOS / "alvo_ok.fasta")
    alvo_sem_f3 = "T" * 21 + alvo_ok[21:]

    violacoes = verificar(leitura.conjunto, alvo_sem_f3, cfg)
    pares = {(v.item, v.regra) for v in violacoes}
    assert ("F3", "nao_encontrado") in pares
    assert not any(v.item == "F3_F2_gap" for v in violacoes)
