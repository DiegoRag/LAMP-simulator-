"""Testes para o gerenciamento de sessão da interface web."""

from pathlib import Path

import pytest

from interface.sessao import Sessao
from primers.entrada import carregar_config

RAIZ = Path(__file__).resolve().parents[1]
cfg = carregar_config(RAIZ / "config" / "primers.toml")
cfg_web = {
    "entrada": {"tamanho_max_kb": 64, "alvo_desenho_max_pb": 2000},
    "analise": {"tm_como_calibracao": True},
    "heterodimeros": {"escala_kcal": [-10.0, -8.0, -6.0, -4.0]},
    "exemplos": [
        {
            "nome": "sintético",
            "conjunto": "tests/dados/conjunto_ok.csv",
            "alvo": "tests/dados/alvo_ok.fasta",
        },
        {
            "nome": "inexistente",
            "conjunto": "dados/nao_existe.csv",
            "alvo": "dados/nao_existe.fasta",
        },
    ],
}


def test_sessao_fluxo_completo():
    sessao = Sessao(RAIZ, cfg, cfg_web)
    assert sessao.exemplos_disponiveis() == ["sintético"]

    sessao.carregar_exemplo("sintético")
    assert sessao.conjunto is not None
    assert sessao.conjunto.nome == "conjunto_ok.csv"
    assert sessao.alvo is not None
    assert sessao.alvo.nome == "alvo_ok.fasta"
    assert sessao.alvo.descricao == "alvo_sintetico_ok"
    assert len(sessao.alvo.sequencia) == 216

    ex = sessao.analisar()
    assert ex.id == 1
    assert ex.tipo == "analise"
    assert ex.rotulo_conjunto == "conjunto_ok"
    assert ex.rotulo_alvo == "alvo_ok"
    assert len(ex.resultado.linhas) == 8
    # conjunto_ok reprova em Tm e ponta na convenção do PrimerExplorer V5 (fase 3)
    assert {(v.item, v.regra) for v in ex.resultado.violacoes} == {
        ("F3", "tm_c"),
        ("F3", "extremidade"),
        ("B3", "tm_c"),
        ("F2", "tm_c"),
        ("F1c", "tm_c"),
        ("F1c", "extremidade"),
        ("B2", "tm_c"),
        ("B1c", "tm_c"),
    }

    d = sessao.desenhar()
    assert d.id == 2
    assert d.tipo == "desenho"
    assert d.rotulo_conjunto == "—"
    assert len(d.resultado.conjuntos) == 4  # fase 4: só núcleos F2–B2c distintos

    p = sessao.analisar_proposta(2, 1)
    assert p.id == 3
    assert p.rotulo_conjunto == "proposta 1 do desenho 2"
    assert [linha["nome"] for linha in p.resultado.linhas] == [
        "F3",
        "B3",
        "F2",
        "F1c",
        "B2",
        "B1c",
        "FIP",
        "BIP",
    ]
    assert p.resultado.violacoes == []

    u_an = sessao.ultima("analise")
    assert u_an is not None
    assert u_an.id == 3

    u_des = sessao.ultima("desenho")
    assert u_des is not None
    assert u_des.id == 2

    assert sessao.execucao(1) is ex

    with pytest.raises(KeyError):
        sessao.execucao(99)

    with pytest.raises(KeyError):
        sessao.analisar_proposta(1, 1)

    with pytest.raises(KeyError):
        sessao.analisar_proposta(2, 9)


def test_sessao_validacoes_e_erros():
    sessao = Sessao(RAIZ, cfg, cfg_web)

    with pytest.raises(ValueError, match="^Carregue um conjunto de primers antes de analisar.$"):
        sessao.analisar()

    with pytest.raises(ValueError, match="^Carregue a sequência-alvo antes de desenhar.$"):
        sessao.desenhar()

    with pytest.raises(ValueError) as excinfo_hdr:
        sessao.carregar_conjunto("x.csv", b"nome,seq\nF3,ACGT\n")
    assert str(excinfo_hdr.value) == "x.csv: Cabeçalho inválido na linha 1: esperado 'nome,sequencia'"

    with pytest.raises(ValueError) as excinfo_tam:
        sessao.carregar_conjunto("grande.csv", b"x" * (64 * 1024 + 1))
    assert str(excinfo_tam.value) == "grande.csv: arquivo maior que 64 KB."

    with pytest.raises(ValueError) as excinfo_alvo:
        sessao.carregar_alvo("invalido.fasta", b">a\nACGX\n")
    msg_alvo = str(excinfo_alvo.value)
    assert msg_alvo.startswith("invalido.fasta: ")
    assert "entrada" not in msg_alvo

    with pytest.raises(ValueError) as excinfo_ex:
        sessao.carregar_exemplo("inexistente")
    assert str(excinfo_ex.value) == "Exemplo indisponível: inexistente."


def test_sessao_alvo_desenho_max_pb():
    cfg_web_curto = {
        "entrada": {"tamanho_max_kb": 64, "alvo_desenho_max_pb": 100},
        "analise": {"tm_como_calibracao": True},
        "heterodimeros": {"escala_kcal": [-10.0, -8.0, -6.0, -4.0]},
        "exemplos": [
            {
                "nome": "sintético",
                "conjunto": "tests/dados/conjunto_ok.csv",
                "alvo": "tests/dados/alvo_ok.fasta",
            }
        ],
    }
    sessao = Sessao(RAIZ, cfg, cfg_web_curto)
    sessao.carregar_exemplo("sintético")
    with pytest.raises(ValueError) as excinfo:
        sessao.desenhar()
    assert str(excinfo.value) == "O alvo tem 216 pb; o desenho pela interface aceita até 100 pb."
