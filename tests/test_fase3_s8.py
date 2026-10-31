"""Fase 3, sessão 8: a interface mostra a Tm de referência e a Tm da reação, e as novas condições."""

from pathlib import Path

from starlette.testclient import TestClient

from interface.app import criar_app
from interface.sessao import Sessao
from interface.telas import tela_analise, tela_parametros, texto_reacao
from primers.entrada import carregar_config
from primers.termo import Condicoes

RAIZ = Path(__file__).resolve().parents[1]
DADOS = RAIZ / "tests" / "dados"
ESCALA = [-10.0, -8.0, -6.0, -4.0]
cfg = carregar_config(RAIZ / "config" / "primers.toml")
cfg_web = {
    "entrada": {"tamanho_max_kb": 64, "alvo_desenho_max_pb": 2000},
    "analise": {"tm_como_calibracao": True},
    "heterodimeros": {"escala_kcal": ESCALA},
    "exemplos": [
        {"nome": "sintético", "conjunto": "tests/dados/conjunto_ok.csv", "alvo": "tests/dados/alvo_ok.fasta"},
    ],
}


def test_texto_reacao():
    assert texto_reacao(Condicoes()) == "K⁺ 50 mM, NH₄⁺ 20 mM, Mg²⁺ 8 mM, dNTP 1,4 mM cada, 65 °C"
    cond = Condicoes(na_mM=10.0, k_mM=0.0, nh4_mM=0.0, temp_C=60.0)
    assert texto_reacao(cond) == "Na⁺ 10 mM, Mg²⁺ 8 mM, dNTP 1,4 mM cada, 60 °C"


def test_tabela_de_primers_mostra_a_tm_da_reacao():
    sessao = Sessao(RAIZ, cfg, cfg_web)
    sessao.carregar_exemplo("sintético")
    dados = tela_analise(sessao.analisar(), cfg, tm_como_calibracao=True, escala_kcal=ESCALA)
    por_nome = {p["nome"]: p for p in dados["primers"]}
    assert (por_nome["F1c"]["tm"], por_nome["F1c"]["tm_reacao"]) == ("58,90", "72,17")
    assert (por_nome["B2"]["tm"], por_nome["B2"]["tm_reacao"]) == ("52,89", "59,66")
    assert por_nome["FIP"]["tm_reacao"] == "—"


def test_parametros_mostram_as_duas_convencoes():
    nova = tela_parametros(Sessao(RAIZ, cfg, cfg_web), cfg)
    assert [r["nome"] for r in nova["reacao"]] == [
        "Monovalentes",
        "Mg²⁺",
        "dNTP",
        "Temperatura da reação",
        "Tm de referência",
        "ΔG de ponta",
    ]
    assert nova["reacao"][1]["nota"] == "Livre: 2,4 mM, o total menos o dNTP total."
    assert nova["reacao"][4] == {
        "nome": "Tm de referência",
        "valor": "100 nM · Na⁺ 50 mM · Mg²⁺ 4 mM",
        "nota": "Convenção do PrimerExplorer V5, a das faixas da Eiken.",
        "destaque": True,
    }
    assert nova["regras"][2] == {
        "nome": "Ponta estável",
        "valor": "ΔG°37 das 6 bases ≤ −4 kcal/mol · 3′ de F3, F2, B2, B3, LF, LB; 5′ de F1c e B1c",
        "marca": "",
    }
    assert nova["regras"][3]["marca"] == "pendente"
    assert nova["regras"][4]["valor"] == "Só reportados, com a ponta 3′ pareada — limiar a definir com dado"

    sessao = Sessao(RAIZ, cfg, cfg_web)
    sessao.carregar_exemplo("sintético")
    sessao.analisar()
    cal = tela_parametros(sessao, cfg)["calibracao"]
    assert cal["linhas"][2]["pontos"][1]["cheio"] is False
    assert cal["descricao"] == (
        "Gráfico de Tm por tipo de primer: Externos, faixa 59–61: F3 54,94 (reação: 58,75), "
        "B3 56,61 (reação: 60,74); Internos, ponta 3′, faixa 59–61: F2 51,27 (reação: 59,09), "
        "B2 52,89 (reação: 59,66); Internos, ponta 5′, faixa 64–66: F1c 58,90 (reação: 72,17), "
        "B1c 57,17 (reação: 71,94)."
    )


def test_paginas_html_com_as_duas_tm():
    cliente = TestClient(criar_app(RAIZ))
    cliente.post(
        "/conjunto",
        files={"arquivo": ("conjunto_ok.csv", (DADOS / "conjunto_ok.csv").read_bytes(), "text/csv")},
        follow_redirects=False,
    )
    cliente.post(
        "/alvo",
        files={"arquivo": ("alvo_ok.fasta", (DADOS / "alvo_ok.fasta").read_bytes(), "text/plain")},
        follow_redirects=False,
    )
    cliente.post("/analisar", follow_redirects=False)
    aba_primers = cliente.get("/analise/1/primers").text
    assert '<span class="nc">Tm</span> reação' in aba_primers
    assert "72,17" in aba_primers
    aba_dimeros = cliente.get("/analise/1/dimeros").text
    assert "diagonal = homodímero pela ponta 3′ (a coluna Homodímero" in aba_dimeros
    parametros = cliente.get("/parametros").text
    assert "Tm de referência (PrimerExplorer V5)" in parametros
    assert "Tm da reação (informativa)" in parametros
    assert "Opção A adotada" in parametros
    assert "convenção do PrimerExplorer V5 (oligo 100 nM · Na⁺ 50 mM · Mg²⁺ 4 mM)" in parametros
    assert "Faixa do guia" in parametros
    assert "dNTP × 4" not in parametros
    assert "Decisão pendente" not in parametros
