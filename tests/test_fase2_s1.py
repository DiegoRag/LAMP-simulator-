"""Testes de unidade da Sessão 1 (Fase 2): modelo IUPAC/RNA, casamento, termo e contrato v2."""

import copy
import json
from pathlib import Path

import pytest

from esteira.contrato import validar_selecao
from primers.casamento import Diferenca, Posicao, melhor_posicao, ocorrencias_compativeis
from primers.modelo import ConjuntoLAMP, Regiao, bases_compativeis, normalizar, rc
from primers.termo import Condicoes, Estrutura, dg_extremidade, gc_pct, grampo, heterodimero, homodimero, tm

EXEMPLO = Path(__file__).parent / "dados" / "selecao_v1_exemplo.json"


@pytest.fixture
def exemplo():
    return json.loads(EXEMPLO.read_text(encoding="utf-8"))


def test_normalizar_rna_iupac_e_notacao_5p_3p():
    assert normalizar("5'-AUGC-3'") == "ATGC"
    assert normalizar("5′ - augc - 3′") == "ATGC"
    assert normalizar("5′–ACGT–3′") == "ACGT"
    assert normalizar("acgturykmswbdhvn") == "ACGTTRYKMSWBDHVN"
    assert normalizar("acgurykmswbdhvn") == "ACGTRYKMSWBDHVN"
    assert rc("ACGTRYKMSWBDHVN") == "NBDHVWSKMRYACGT"
    assert rc(rc("ACGTRYKMSWBDHVN")) == "ACGTRYKMSWBDHVN"


def test_erro_caractere_invalido_mostra_posicao_sem_despejar_sequencia():
    seq_longa = ("A" * 500) + "X" + ("C" * 500)
    with pytest.raises(ValueError) as exc:
        normalizar(seq_longa)
    assert "'X'" in str(exc.value)
    assert "posição 501" in str(exc.value)
    assert len(str(exc.value)) < 120

    with pytest.raises(ValueError) as exc_vazio:
        ConjuntoLAMP(
            F3="",
            B3="ACGT",
            F2="ACGT",
            F1c="ACGT",
            B2="ACGT",
            B1c="ACGT",
        )
    assert "Sequência vazia para o primer 'F3'" in str(exc_vazio.value)


def test_casamento_iupac_por_conjunto_de_bases():
    assert bases_compativeis("R", "A") is True
    assert bases_compativeis("R", "G") is True
    assert bases_compativeis("R", "C") is False
    assert bases_compativeis("N", "T") is True

    assert ocorrencias_compativeis("AAACCCTTTGGG", "AAAYCC") == [Regiao(0, 6)]

    assert melhor_posicao(
        "CGCATRTGGATTCTGTAC",
        "AAAAACGCATATGGATTCTGTACTTTTT",
        max_erros=0,
    ) == Posicao(regiao=Regiao(5, 23), diferencas=(), fita="+")

    assert melhor_posicao(
        "CGCATRTGGATTCTGTAC",
        "AAAAACGCATCTGGATTCTGTACTTTTT",
        max_erros=3,
    ) == Posicao(regiao=Regiao(5, 23), diferencas=(Diferenca(10, "R", "C"),), fita="+")


def test_termo_iupac_gc_tm_e_dg_extremidade():
    assert gc_pct("ACRS") == pytest.approx(62.5)
    assert gc_pct("NNNN") == pytest.approx(50.0)
    assert gc_pct("WWWW") == pytest.approx(0.0)
    assert gc_pct("SSSS") == pytest.approx(100.0)

    cond = Condicoes()
    seq_deg = "GCGCGGCCGCGCGGCCGCGY"
    v_c = "GCGCGGCCGCGCGGCCGCGC"
    v_t = "GCGCGGCCGCGCGGCCGCGT"
    tm_c = tm(v_c, cond, 200.0)
    tm_t = tm(v_t, cond, 200.0)

    assert tm(seq_deg, cond, 200.0, modo_iupac="media") == pytest.approx((tm_c + tm_t) / 2.0)
    assert tm(seq_deg, cond, 200.0, modo_iupac="pior_caso") == pytest.approx(min(tm_c, tm_t))

    dg_c = dg_extremidade("CCACGTAC", "3", 6, cond)
    dg_t = dg_extremidade("CCACGTAT", "3", 6, cond)

    assert dg_extremidade("CCACGTAY", "3", 6, cond, modo_iupac="media") == pytest.approx(
        (dg_c + dg_t) / 2.0
    )
    assert dg_extremidade("CCACGTAY", "3", 6, cond, modo_iupac="pior_caso") == pytest.approx(
        max(dg_c, dg_t)
    )


def test_termo_limite_60_nt_e_sequencia_vazia():
    cond = Condicoes()
    seq_64 = "ACGT" * 16

    assert grampo(seq_64, cond) == Estrutura(tm_c=0.0, dg_kcal=0.0, encontrada=False)
    assert homodimero(seq_64, cond) == Estrutura(tm_c=0.0, dg_kcal=0.0, encontrada=False)
    assert heterodimero(seq_64, "ACGTACGTACGTACGTACGT", cond) == Estrutura(
        tm_c=0.0, dg_kcal=0.0, encontrada=False
    )

    with pytest.raises(ValueError):
        tm("", cond, 200.0)
    with pytest.raises(ValueError):
        grampo("", cond)
    with pytest.raises(ValueError):
        homodimero("", cond)
    with pytest.raises(ValueError):
        heterodimero("", "ACGT", cond)


def test_contrato_v2_aceita_iupac_e_ids_flexiveis(exemplo):
    dados_v2 = copy.deepcopy(exemplo)
    dados_v2["esquema"] = "lamp-sim/selecao@2"
    dados_v2["conjuntos"][0]["conjunto_id"] = "fosA3"
    dados_v2["conjuntos"][0]["por_alelo"][0]["alelo_id"] = "AP018671.1"
    dados_v2["conjuntos"][0]["primers"]["FIP"] = "aagtttcgtgctgacgtgRGCTAAAGACACGCATATG"

    assert validar_selecao(dados_v2) == []

    dados_v1 = copy.deepcopy(dados_v2)
    dados_v1["esquema"] = "lamp-sim/selecao@1"
    problemas_v1 = validar_selecao(dados_v1)
    assert any("conjunto_id" in p for p in problemas_v1)
    assert any("alelo_id" in p for p in problemas_v1)
    assert any("FIP" in p for p in problemas_v1)
