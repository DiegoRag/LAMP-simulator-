"""Testes dos detalhes da análise em abas (Fase 2, s3)."""

import re
import tomllib
from pathlib import Path

import pytest
from starlette.testclient import TestClient

from interface.app import criar_app
from interface.sessao import Sessao
from interface.telas_fluxo import tela_detalhe, tela_detalhe_vazia
from primers.entrada import carregar_config

RAIZ = Path(__file__).resolve().parents[1]
DADOS = RAIZ / "tests" / "dados"

URLS_ABAS = [
    ("resumo", "/analise/1"),
    ("regioes", "/analise/1/regioes"),
    ("primers", "/analise/1/primers"),
    ("violacoes", "/analise/1/violacoes"),
    ("dimeros", "/analise/1/dimeros"),
]


@pytest.fixture
def cliente_com_analise():
    """Cliente Starlette com o fluxo de teste completo executado."""
    app = criar_app(RAIZ)
    cliente = TestClient(app)
    try:
        res1 = cliente.post(
            "/testar/primers",
            files={"arquivo": ("conjunto_ok.csv", (DADOS / "conjunto_ok.csv").read_bytes(), "text/csv")},
            follow_redirects=False,
        )
        assert res1.status_code == 303
        assert res1.headers["location"] == "/testar/alvo"

        res2 = cliente.post(
            "/testar/alvo",
            files={"arquivo": ("alvo_ok.fasta", (DADOS / "alvo_ok.fasta").read_bytes(), "text/plain")},
            follow_redirects=False,
        )
        assert res2.status_code == 303
        assert res2.headers["location"] == "/testar/resultado/1"

        yield cliente
    finally:
        app.state.analisador.encerrar()


def test_cinco_urls_respondem_200_com_aria_current_e_marca(cliente_com_analise):
    """Verifica que as 5 abas respondem 200, com aria-current só na própria aba e body class marca."""
    for chave_atual, url in URLS_ABAS:
        res = cliente_com_analise.get(url)
        assert res.status_code == 200
        assert 'body class="marca' in res.text

        # A barra superior deve ter o link Análise com aria-current="page"
        assert '<a class="marca-nav-link" href="/analise" aria-current="page">Análise</a>' in res.text

        # Na navegação de abas, apenas a aba atual tem aria-current="page"
        for chave_outra, url_outra in URLS_ABAS:
            padrao = rf'href="{re.escape(url_outra)}"[^>]*aria-current="page"'
            if chave_outra == chave_atual:
                assert re.search(padrao, res.text) is not None, f"Aba ativa {url_outra} sem aria-current"
            else:
                assert re.search(padrao, res.text) is None, f"Aba inativa {url_outra} com aria-current"


def test_rotas_validacao_e_exportar(cliente_com_analise):
    """Verifica 404 para aba desconhecida, 404 para id inexistente e 200 para exportar."""
    res_qualquer = cliente_com_analise.get("/analise/1/qualquer")
    assert res_qualquer.status_code == 404

    res_99 = cliente_com_analise.get("/analise/99")
    assert res_99.status_code == 404

    res_exp = cliente_com_analise.get("/analise/1/exportar")
    assert res_exp.status_code == 200
    assert res_exp.headers["content-type"] == "application/zip"


def test_rota_analise_redirecionamento_e_vazia(cliente_com_analise):
    """Verifica /analise com e sem execução na sessão."""
    res_com = cliente_com_analise.get("/analise", follow_redirects=False)
    assert res_com.status_code == 303
    assert res_com.headers["location"] == "/analise/1"

    app_vazio = criar_app(RAIZ)
    cliente_vazio = TestClient(app_vazio)
    try:
        res_sem = cliente_vazio.get("/analise", follow_redirects=False)
        assert res_sem.status_code == 200
        assert "Nenhuma análise ainda" in res_sem.text
        assert 'body class="marca' in res_sem.text
    finally:
        app_vazio.state.analisador.encerrar()


def test_tm_calibracao_propagacao_e_alternancia(cliente_com_analise):
    """Verifica que tm_calibracao propaga para os links de aba e alterna no link de calibração."""
    res = cliente_com_analise.get("/analise/1/primers?tm_calibracao=0")
    assert res.status_code == 200

    for _, sufixo in [
        ("resumo", ""),
        ("regioes", "regioes"),
        ("primers", "primers"),
        ("violacoes", "violacoes"),
        ("dimeros", "dimeros"),
    ]:
        url_com_sufixo = f"/analise/1/{sufixo}" if sufixo else "/analise/1"
        assert f'href="{url_com_sufixo}?tm_calibracao=0"' in res.text

    assert 'href="/analise/1/primers?tm_calibracao=1"' in res.text


def test_tela_detalhe_direto():
    """Verifica tela_detalhe e tela_detalhe_vazia diretamente."""
    cfg_primers = carregar_config(RAIZ / "config" / "primers.toml")
    with open(RAIZ / "config" / "web.toml", "rb") as f:
        cfg_web = tomllib.load(f)

    sessao = Sessao(RAIZ, cfg_primers, cfg_web)
    sessao.carregar_conjunto("conjunto_ok.csv", (DADOS / "conjunto_ok.csv").read_bytes())
    sessao.carregar_alvo("alvo_ok.fasta", (DADOS / "alvo_ok.fasta").read_bytes())
    execucao = sessao.analisar()

    dados = tela_detalhe(execucao, "resumo", cfg_primers, cfg_web)
    assert dados["aba"] == "resumo"
    assert dados["par"] == "conjunto_ok × alvo_ok"
    assert dados["url_resultado"] == f"/testar/resultado/{execucao.id}"
    assert dados["url_exportar"] == f"/analise/{execucao.id}/exportar"
    assert dados["url_outro"] == "/testar"

    abas = dados["abas"]
    assert len(abas) == 5
    assert [a["chave"] for a in abas] == ["resumo", "regioes", "primers", "violacoes", "dimeros"]
    assert [a["rotulo"] for a in abas] == ["Resumo", "Regiões no alvo", "Primers", "Violações", "Dímeros"]
    assert abas[0]["href"] == f"/analise/{execucao.id}"
    assert abas[1]["href"] == f"/analise/{execucao.id}/regioes"
    assert abas[2]["href"] == f"/analise/{execucao.id}/primers"
    assert abas[3]["href"] == f"/analise/{execucao.id}/violacoes"
    assert abas[4]["href"] == f"/analise/{execucao.id}/dimeros"

    assert abas[0]["atual"] is True
    assert all(not a["atual"] for a in abas[1:])

    # Contagem só na aba violacoes
    assert isinstance(abas[3]["contagem"], int)
    assert abas[3]["contagem"] > 0
    assert all(a["contagem"] is None for i, a in enumerate(abas) if i != 3)

    # tm_calibracao=0 na aba primers
    dados_tm = tela_detalhe(execucao, "primers", cfg_primers, cfg_web, tm_calibracao=0)
    assert dados_tm["href_alternar_calibracao"] == f"/analise/{execucao.id}/primers?tm_calibracao=1"
    for a in dados_tm["abas"]:
        assert a["href"].endswith("?tm_calibracao=0")
    assert dados_tm["abas"][2]["atual"] is True

    # tela_detalhe_vazia
    vazia = tela_detalhe_vazia()
    assert vazia == {
        "titulo": "Nenhuma análise ainda",
        "texto": "Envie um conjunto de primers e um gene-alvo para ver a análise detalhada aqui.",
        "url_acao": "/testar",
        "rotulo_acao": "Testar primers",
    }


def test_selo_aba_violacoes_coerente_com_kpi():
    """Com conjunto_violacoes.csv × alvo_ok.fasta, o selo da aba Violações bate com o KPI."""
    cfg_primers = carregar_config(RAIZ / "config" / "primers.toml")
    with open(RAIZ / "config" / "web.toml", "rb") as f:
        cfg_web = tomllib.load(f)

    sessao = Sessao(RAIZ, cfg_primers, cfg_web)
    sessao.carregar_conjunto("conjunto_violacoes.csv", (DADOS / "conjunto_violacoes.csv").read_bytes())
    sessao.carregar_alvo("alvo_ok.fasta", (DADOS / "alvo_ok.fasta").read_bytes())
    execucao = sessao.analisar()

    dados = tela_detalhe(execucao, "resumo", cfg_primers, cfg_web)
    kpi_viol = next(k for k in dados["kpis"] if k["rotulo"] == "Violações")
    aba_viol = next(a for a in dados["abas"] if a["chave"] == "violacoes")
    assert aba_viol["contagem"] == int(kpi_viol["valor"])
    assert aba_viol["contagem"] == 7
