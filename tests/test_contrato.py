"""Contrato selecao@1: o exemplo versionado passa, e cada desvio do formato é apontado."""

import json
from pathlib import Path

import pytest

from esteira.contrato import ESQUEMA, ler_selecao, validar_selecao

EXEMPLO = Path(__file__).parent / "dados" / "selecao_v1_exemplo.json"
MEDIDA = "conjuntos[0].por_alelo[0]"


@pytest.fixture
def exemplo():
    return json.loads(EXEMPLO.read_text(encoding="utf-8"))


def test_exemplo_versionado_e_valido(exemplo):
    assert exemplo["esquema"] == ESQUEMA == "lamp-sim/selecao@1"
    assert validar_selecao(exemplo) == []
    assert ler_selecao(EXEMPLO) == exemplo


def test_raiz_que_nao_e_objeto():
    assert validar_selecao([]) == ["a raiz deve ser um objeto JSON"]


def test_versao_desconhecida_e_recusada(exemplo):
    exemplo["esquema"] = "lamp-sim/selecao@99"
    assert validar_selecao(exemplo) == [
        "esquema 'lamp-sim/selecao@99' não reconhecido; esperado 'lamp-sim/selecao@1' ou 'lamp-sim/selecao@2'"
    ]


def test_nome_de_campo_errado_nao_passa_calado(exemplo):
    medida = exemplo["conjuntos"][0]["por_alelo"][0]
    medida["mismatches_3p"] = medida.pop("incompat_3p")
    assert validar_selecao(exemplo) == [
        f"{MEDIDA}: falta o campo 'incompat_3p'",
        f"{MEDIDA}: campo desconhecido 'mismatches_3p'",
    ]


def test_desenho_sem_aviso(exemplo):
    del exemplo["conjuntos"][0]["aviso"]
    assert validar_selecao(exemplo) == ["conjuntos[0]: conjunto de desenho sem aviso (ADR-0013, regra 2)"]


def test_conjunto_com_ponto_e_alelo_sem_ponto(exemplo):
    # No artigo, "fosA3" nomeia um conjunto e um alelo; o contrato não deixa os dois se confundirem.
    exemplo["conjuntos"][0]["conjunto_id"] = "fosA3"
    exemplo["conjuntos"][0]["por_alelo"][0]["alelo_id"] = "fosA.3"
    assert validar_selecao(exemplo) == [
        "conjuntos[0].conjunto_id 'fosA3' deve ser texto com ponto, como 'fosA.1'",
        f"{MEDIDA}.alelo_id 'fosA.3' deve ser texto sem ponto, como 'fosA3'",
    ]


def test_conjunto_repetido(exemplo):
    exemplo["conjuntos"][1]["conjunto_id"] = "sint.1"
    assert validar_selecao(exemplo) == ["conjuntos[1]: conjunto 'sint.1' repetido"]


@pytest.mark.parametrize("valor", [-1, True, 1.5, "2"])
def test_contagem_invalida(exemplo, valor):
    exemplo["conjuntos"][0]["por_alelo"][0]["incompat_total"] = valor
    assert validar_selecao(exemplo) == [f"{MEDIDA}.incompat_total deve ser inteiro >= 0 ou null"]


def test_ponta_3p_nao_passa_do_total(exemplo):
    medida = exemplo["conjuntos"][0]["por_alelo"][0]
    medida["incompat_total"], medida["incompat_3p"] = 1, 2
    assert validar_selecao(exemplo) == [f"{MEDIDA}: incompat_3p (2) maior que incompat_total (1)"]


@pytest.mark.parametrize("valor", [0, 0.5, 1, None])
def test_f_previsto_valido(exemplo, valor):
    exemplo["conjuntos"][0]["por_alelo"][0]["f_previsto"] = valor
    assert validar_selecao(exemplo) == []


@pytest.mark.parametrize("valor", [1.2, -0.1, float("nan"), True, "0.5"])
def test_f_previsto_invalido(exemplo, valor):
    exemplo["conjuntos"][0]["por_alelo"][0]["f_previsto"] = valor
    assert validar_selecao(exemplo) == [f"{MEDIDA}.f_previsto deve ser um número de 0 a 1, ou null"]


def test_primers_do_tubo(exemplo):
    primers = exemplo["conjuntos"][0]["primers"]
    primers["FIP"] = primers["FIP"].lower()
    primers["LoopF"] = primers.pop("BIP")
    assert validar_selecao(exemplo) == [
        "conjuntos[0].primers: falta o campo 'BIP'",
        "conjuntos[0].primers: campo desconhecido 'LoopF'",
        "conjuntos[0].primers.FIP deve ser uma sequência só de A, C, G e T",
    ]


def test_ler_selecao_lista_todos_os_problemas(exemplo, tmp_path):
    del exemplo["conjuntos"][0]["aviso"]
    exemplo["conjuntos"][1]["selecionado"] = "sim"
    caminho = tmp_path / "selecao.json"
    caminho.write_text(json.dumps(exemplo), encoding="utf-8")
    with pytest.raises(ValueError) as erro:
        ler_selecao(caminho)
    mensagem = str(erro.value)
    assert "conjuntos[0]: conjunto de desenho sem aviso (ADR-0013, regra 2)" in mensagem
    assert "conjuntos[1].selecionado deve ser true ou false" in mensagem
