"""Fase 6, sessão 1: busca aproximada vetorizada e erros de base situados na ponta que estende."""

import random
import time
from pathlib import Path

import primers.incompat
from primers.casamento import (
    REGIAO_PARA_PRIMER,
    Diferenca,
    ErroDeBase,
    Posicao,
    distancia_ponta_extensao,
    erros_de_base,
    melhor_posicao,
    ocorrencias_compativeis,
    ponta_extensao,
)
from primers.entrada import carregar_config, ler_alvo, ler_conjunto
from primers.incompat import IncompatRegiao, medir_incompat, regioes_para_json
from primers.modelo import Regiao, bases_compativeis
from primers.restricoes import localizar_detalhado

RAIZ = Path(__file__).resolve().parents[1]
DADOS = RAIZ / "tests" / "dados"
CFG = carregar_config(RAIZ / "config" / "primers.toml")


def _erros(seq: str, texto: str, i: int) -> int:
    return sum(1 for j, b in enumerate(seq) if not bases_compativeis(b, texto[i + j]))


def _melhor_bruta(seq: str, texto: str, ini: int, fim: int, max_erros: int) -> Posicao | None:
    """A busca antiga, base a base em Python: a primeira posição com o menor número de erros."""
    melhor = None
    for i in range(ini, fim - len(seq) + 1):
        e = _erros(seq, texto, i)
        if e <= max_erros and (melhor is None or e < melhor[0]):
            melhor = (e, i)
    if melhor is None:
        return None
    i = melhor[1]
    difs = tuple(
        Diferenca(i + j, b, texto[i + j]) for j, b in enumerate(seq) if not bases_compativeis(b, texto[i + j])
    )
    return Posicao(Regiao(i, i + len(seq)), difs, "+")


def test_busca_vetorizada_da_o_mesmo_que_a_busca_base_a_base():
    rng = random.Random(6)
    for rodada in range(40):
        texto = "".join(rng.choices("ACGT", k=rng.randint(30, 400)))
        if rodada % 2:  # alvo com N e códigos IUPAC: cai na varredura vetorizada
            pos = rng.randrange(len(texto))
            texto = texto[:pos] + rng.choice("NRY") + texto[pos + 1 :]
        k = rng.randint(4, 22)
        if k > len(texto):
            continue
        i0 = rng.randrange(len(texto) - k + 1)
        seq = list(texto[i0 : i0 + k])
        for _ in range(rng.randint(0, 3)):
            j = rng.randrange(k)
            seq[j] = rng.choice("ACGTRYN")
        seq = "".join(seq)
        ini = rng.randrange(0, max(1, len(texto) - k))
        fim = rng.randint(ini + k, len(texto))
        esperado = _melhor_bruta(seq, texto, ini, fim, 3)
        assert melhor_posicao(seq, texto, max_erros=3, janela=Regiao(ini, fim)) == esperado
        exatas = [Regiao(i, i + k) for i in range(len(texto) - k + 1) if _erros(seq, texto, i) == 0]
        assert ocorrencias_compativeis(texto, seq) == exatas


def test_alvo_grande_com_n_nao_trava():
    # 2 Mb sintéticos com um N: antes, a busca base a base em Python levava segundos por região
    rng = random.Random(7)
    alvo = "".join(rng.choices("ACGT", k=2_000_000))
    alvo = alvo[:1000] + "N" + alvo[1001:]
    seq = alvo[1_500_000:1_500_021]
    seq_2erros = "T" + seq[1:10] + ("A" if seq[10] != "A" else "C") + seq[11:]
    inicio = time.perf_counter()
    pos = melhor_posicao(seq_2erros, alvo, max_erros=3)
    exatas = ocorrencias_compativeis(alvo, seq)
    assert time.perf_counter() - inicio < 2.0
    assert pos is not None and pos.regiao == Regiao(1_500_000, 1_500_021)
    assert len(pos.diferencas) == (seq[0] != "T") + 1
    assert Regiao(1_500_000, 1_500_021) in exatas


def test_ponta_que_estende_por_regiao():
    assert [ponta_extensao(r) for r in ("F3", "F2", "F1", "B1c", "B2c", "B3c", "LF", "LB")] == [
        "3'", "3'", "5'", "5'", "3'", "3'", "3'", "3'",
    ]
    assert REGIAO_PARA_PRIMER["F1"] == "F1c" and REGIAO_PARA_PRIMER["B2c"] == "B2"
    assert REGIAO_PARA_PRIMER["B3c"] == "B3" and REGIAO_PARA_PRIMER["LB"] == "LB"
    # a regra é uma só: a do incompat (fase 5) é a mesma função
    assert primers.incompat.distancia_ponta_extensao is distancia_ponta_extensao
    assert distancia_ponta_extensao("F1", Regiao(81, 103), 102) == 0
    assert distancia_ponta_extensao("B1c", Regiao(108, 127), 126) == 18


def test_erros_de_base_contam_no_primer():
    f2 = Posicao(Regiao(31, 49), (Diferenca(48, "C", "A"),))
    assert erros_de_base("F2", f2) == [
        ErroDeBase("F2", "F2", "+", 49, 18, 18, "C", "A", "3'", 0),
    ]
    # LF é lido pelo complemento reverso: base 61 da fita + = base 13 do primer, complementada
    lf = Posicao(Regiao(55, 73), (Diferenca(60, "G", "T"),))
    assert erros_de_base("LF", lf) == [ErroDeBase("LF", "LF", "+", 61, 13, 18, "C", "A", "3'", 5)]
    # F1: a última base da região na fita é a ponta 5' de F1c, a do laço
    f1 = Posicao(Regiao(81, 103), (Diferenca(102, "A", "G"), Diferenca(92, "R", "T")), fita="-")
    assert erros_de_base("F1", f1) == [
        ErroDeBase("F1", "F1c", "-", 103, 1, 22, "T", "C", "5'", 0),
        ErroDeBase("F1", "F1c", "-", 93, 11, 22, "Y", "A", "5'", 10),
    ]
    b1c = Posicao(Regiao(108, 127), (Diferenca(108, "C", "T"),))
    assert erros_de_base("B1c", b1c) == [ErroDeBase("B1c", "B1c", "+", 109, 1, 19, "C", "T", "5'", 0)]
    b2c = Posicao(Regiao(157, 177), (Diferenca(157, "G", "A"),))
    assert erros_de_base("B2c", b2c) == [ErroDeBase("B2c", "B2", "+", 158, 20, 20, "C", "T", "3'", 0)]
    assert erros_de_base("F3", Posicao(Regiao(0, 21), ())) == []


def _com_troca(alvo: str, pos: int, nova: str) -> str:
    return alvo[:pos] + nova + alvo[pos + 1 :]


def test_erro_no_meio_e_na_ponta_5_de_f1c_ficam_distintos():
    alvo = ler_alvo(DADOS / "alvo_ok.fasta")
    conjunto = ler_conjunto(DADOS / "conjunto_ok.csv").conjunto
    assert alvo[81:103] == "CCAATAACACTTGGGGGTTTAC"  # F1 = rc(F1c)

    loc_ponta = localizar_detalhado(conjunto, _com_troca(alvo, 102, "T"), CFG)
    loc_meio = localizar_detalhado(conjunto, _com_troca(alvo, 92, "G"), CFG)
    f1_ponta = next(i for i in medir_incompat(conjunto, loc_ponta, 5).por_regiao if i.regiao == "F1")
    f1_meio = next(i for i in medir_incompat(conjunto, loc_meio, 5).por_regiao if i.regiao == "F1")
    assert f1_ponta == IncompatRegiao("F1", True, 1, 1)
    assert f1_meio == IncompatRegiao("F1", True, 1, 0)

    # a mensagem da violação diz a ponta certa: 5' de F1c, não 3'
    assert [v.mensagem for v in loc_ponta.violacoes] == [
        "Região F1 não encontrada exata no alvo; melhor posição aproximada 82–103 "
        "(posição 103: C esperado, T no alvo (a 0 nt da ponta 5' de F1c))"
    ]
    assert [v.mensagem for v in loc_meio.violacoes] == [
        "Região F1 não encontrada exata no alvo; melhor posição aproximada 82–103 "
        "(posição 93: T esperado, G no alvo (a 10 nt da ponta 5' de F1c))"
    ]

    assert regioes_para_json(medir_incompat(conjunto, loc_ponta, 5))[2] == {
        "regiao": "F1",
        "encontrada": True,
        "total": 1,
        "ponta": 1,
    }
