"""Testes para o módulo de casamento com tolerância por distância de Hamming."""

from pathlib import Path

from primers.casamento import Diferenca, Posicao, melhor_posicao, regioes_do_conjunto
from primers.entrada import ler_alvo
from primers.modelo import ConjuntoLAMP, Regiao

RAIZ = Path(__file__).resolve().parents[1]


def test_melhor_posicao_casamento_exato():
    assert melhor_posicao("CCGG", "AAAACCCCGGGGTTTT") == Posicao(Regiao(6, 10), ())


def test_melhor_posicao_empate_menor_inicio():
    assert melhor_posicao("CCAG", "AAAACCCCGGGGTTTT") == Posicao(
        Regiao(5, 9),
        (Diferenca(7, "A", "C"),),
    )


def test_melhor_posicao_vazia_ou_maior_que_alvo():
    assert melhor_posicao("ACGTACGT", "ACG") is None
    assert melhor_posicao("", "ACGT") is None


def test_melhor_posicao_alvo_ok():
    alvo = ler_alvo(RAIZ / "tests/dados/alvo_ok.fasta")
    assert melhor_posicao("TTACAGAACATACACGTCAG", alvo) == Posicao(
        Regiao(50, 70),
        (Diferenca(55, "G", "T"),),
    )


def test_regioes_do_conjunto():
    conjunto = ConjuntoLAMP(
        F3="AAAC",
        B3="CCCA",
        F2="GGGT",
        F1c="TTTG",
        B2="ACGG",
        B1c="CATG",
        LF="GACC",
    )
    assert regioes_do_conjunto(conjunto) == [
        ("F3", "AAAC"),
        ("F2", "GGGT"),
        ("F1", "CAAA"),
        ("B1c", "CATG"),
        ("B2c", "CCGT"),
        ("B3c", "TGGG"),
        ("LF", "GGTC"),
    ]
