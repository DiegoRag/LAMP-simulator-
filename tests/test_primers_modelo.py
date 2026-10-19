"""Testes unitários para o modelo de dados de primers LAMP."""

import pytest

from primers.modelo import ConjuntoLAMP, Regiao, normalizar, primers_de_regioes, rc


def test_rc():
    assert rc("ACGT") == "ACGT"
    assert rc("AAAC") == "GTTT"
    assert rc(rc("GATTACAGG")) == "GATTACAGG"


def test_normalizar():
    assert normalizar(" acgt \n") == "ACGT"
    assert normalizar("ACGU") == "ACGT"
    with pytest.raises(ValueError) as exc_info:
        normalizar("ACGX")
    assert "'X'" in str(exc_info.value)
    assert "posição 4" in str(exc_info.value)


def test_regiao_validacao():
    with pytest.raises(ValueError):
        Regiao(5, 5)
    with pytest.raises(ValueError):
        Regiao(-1, 3)


def test_fip_bip():
    conj_sem_ligante = ConjuntoLAMP(
        F3="AAAA",
        B3="CCCC",
        F2="CCCC",
        F1c="AAAA",
        B2="TTTT",
        B1c="GGGG",
        ligante="",
    )
    assert conj_sem_ligante.FIP == "AAAACCCC"
    assert conj_sem_ligante.BIP == "GGGGTTTT"

    conj_com_ligante = ConjuntoLAMP(
        F3="AAAA",
        B3="CCCC",
        F2="CCCC",
        F1c="AAAA",
        B2="TTTT",
        B1c="GGGG",
        ligante="TTTT",
    )
    assert conj_com_ligante.FIP == "AAAATTTTCCCC"


def test_convencao_orientacao():
    blocos = [
        ("F3", "ATGCGTACGTTAGCATGCAA"),
        ("F2", "GGATCCAAGCTTGAATTCGA"),
        ("LF", "TTGACCGGTACCAGTCGA"),
        ("F1", "CCTAGGACTAGTGCGGCCGC"),
        ("B1c", "GTCGACATCGATAGATCTCA"),
        ("LB", "AGCTAGCTAGGCATTCGT"),
        ("B2c", "CTGCAGTCTAGAGCTCGCAT"),
        ("B3c", "TCGAGCAAGGTACCATGGTT"),
    ]

    alvo_partes = []
    regioes = {}
    pos = 0
    for i, (nome, seq) in enumerate(blocos):
        if i > 0:
            alvo_partes.append("AAAAA")
            pos += len("AAAAA")
        ini = pos
        fim = pos + len(seq)
        regioes[nome] = Regiao(ini, fim)
        alvo_partes.append(seq)
        pos = fim

    alvo = "".join(alvo_partes)
    conj = primers_de_regioes(alvo, regioes)

    assert conj.F3 == "ATGCGTACGTTAGCATGCAA"
    assert conj.F2 == "GGATCCAAGCTTGAATTCGA"
    assert conj.B1c == "GTCGACATCGATAGATCTCA"
    assert conj.LB == "AGCTAGCTAGGCATTCGT"
    assert conj.F1c == "GCGGCCGCACTAGTCCTAGG"
    assert conj.B3 == "AACCATGGTACCTTGCTCGA"
    assert conj.B2 == "ATGCGAGCTCTAGACTGCAG"
    assert conj.LF == "TCGACTGGTACCGGTCAA"


def test_chave_obrigatoria_ausente():
    blocos = [
        ("F3", "ATGCGTACGTTAGCATGCAA"),
        ("F2", "GGATCCAAGCTTGAATTCGA"),
        ("LF", "TTGACCGGTACCAGTCGA"),
        ("F1", "CCTAGGACTAGTGCGGCCGC"),
        ("B1c", "GTCGACATCGATAGATCTCA"),
        ("LB", "AGCTAGCTAGGCATTCGT"),
        ("B2c", "CTGCAGTCTAGAGCTCGCAT"),
        ("B3c", "TCGAGCAAGGTACCATGGTT"),
    ]
    regioes = {}
    pos = 0
    for i, (nome, seq) in enumerate(blocos):
        if i > 0:
            pos += 5
        ini = pos
        fim = pos + len(seq)
        if nome != "B3c":
            regioes[nome] = Regiao(ini, fim)
        pos = fim

    with pytest.raises(KeyError) as exc_info:
        primers_de_regioes("A" * pos, regioes)
    assert "B3c" in str(exc_info.value)
