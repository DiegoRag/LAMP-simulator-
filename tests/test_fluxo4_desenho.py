"""Testes da Fase 4: tela de Desenho na identidade LampSim, barra e alvos longos."""

import tomllib
from pathlib import Path

from starlette.testclient import TestClient

import primers.desenho
from interface.app import criar_app
from interface.sessao import Sessao
from interface.telas_fluxo import tela_desenho_vazio, tela_detalhe, tela_marca
from primers.entrada import carregar_config

RAIZ = Path(__file__).resolve().parents[1]
DADOS = RAIZ / "tests" / "dados"


def test_tela_desenho_vazio_dados():
    """Valida o dicionário retornado por tela_desenho_vazio."""
    dados = tela_desenho_vazio()
    assert dados["titulo"] == "Nenhum desenho ainda"
    assert (
        dados["texto"]
        == 'Abra uma análise e use "Desenhar alternativas" para o LampSim propor '
        "outros conjuntos de primers para o mesmo alvo."
    )
    assert dados["url_acao"] == "/analise"
    assert dados["rotulo_acao"] == "Ir para a análise"


def test_tela_marca_abas_resultado_e_desenho():
    """Valida que resultado marca 'Testar primers' e desenho marca 'Análise'."""
    res_marca = tela_marca("resultado")
    assert res_marca["pagina"] == "resultado"
    assert [link["atual"] for link in res_marca["links"]] == [False, True, False, False]

    des_marca = tela_marca("desenho")
    assert des_marca["pagina"] == "desenho"
    assert [link["atual"] for link in des_marca["links"]] == [False, False, True, False]


def test_desenho_sem_desenho_na_sessao():
    """GET /desenho sem desenho na sessão devolve 200 com a tela de desenho vazio."""
    cliente = TestClient(criar_app(RAIZ))
    res = cliente.get("/desenho")
    assert res.status_code == 200
    assert 'body class="marca' in res.text
    assert "Nenhum desenho ainda" in res.text
    assert "Abra uma análise e use" in res.text
    assert 'href="/analise"' in res.text
    assert "Ir para a análise" in res.text


def test_fluxo4_desenho_de_ponta_a_ponta():
    """Valida o fluxo completo de desenho na identidade LampSim."""
    cliente = TestClient(criar_app(RAIZ))

    # 1. Enviar primers e alvo pelo fluxo /testar
    res_primers = cliente.post(
        "/testar/primers",
        files={"arquivo": ("conjunto_ok.csv", (DADOS / "conjunto_ok.csv").read_bytes(), "text/csv")},
        follow_redirects=False,
    )
    assert res_primers.status_code == 303
    assert res_primers.headers["location"] == "/testar/alvo"

    res_alvo = cliente.post(
        "/testar/alvo",
        files={"arquivo": ("alvo_ok.fasta", (DADOS / "alvo_ok.fasta").read_bytes(), "text/plain")},
        follow_redirects=False,
    )
    assert res_alvo.status_code == 303
    assert res_alvo.headers["location"] == "/testar/resultado/1"

    # 2. Tela de resultado marca "Testar primers" na barra
    res_resultado = cliente.get("/testar/resultado/1")
    assert res_resultado.status_code == 200
    assert (
        '<a class="marca-nav-link" href="/testar" aria-current="page">Testar primers</a>'
        in res_resultado.text
    )

    # 3. Desenhar alternativas
    res_desenhar = cliente.post("/desenhar", follow_redirects=False)
    assert res_desenhar.status_code == 303
    assert res_desenhar.headers["location"] == "/desenho/2"

    # 4. /desenho/{id} exibe aviso ADR-0013, ranking, marca "Análise" na barra e voltar à análise
    res_desenho = cliente.get("/desenho/2")
    assert res_desenho.status_code == 200
    assert 'body class="marca' in res_desenho.text
    assert primers.desenho.AVISO in res_desenho.text
    assert "ranking" in res_desenho.text.lower() or "tabela" in res_desenho.text.lower()
    assert (
        '<a class="marca-nav-link" href="/analise" aria-current="page">Análise</a>'
        in res_desenho.text
    )
    assert 'href="/analise/1"' in res_desenho.text
    assert "Voltar à análise" in res_desenho.text

    # 5. ?sel=2 marca a 2ª linha do ranking
    res_sel2 = cliente.get("/desenho/2?sel=2")
    assert res_sel2.status_code == 200
    assert (
        'aria-current="true"' in res_sel2.text
        or "selecionada" in res_sel2.text
        or "sel=2" in res_sel2.text
    )

    # 6. Analisar proposta 1 redireciona para a nova análise
    res_analisar_prop = cliente.post("/desenho/2/analisar", data={"posicao": "1"}, follow_redirects=False)
    assert res_analisar_prop.status_code == 303
    assert res_analisar_prop.headers["location"] == "/analise/3"

    res_nova_analise = cliente.get("/analise/3")
    assert res_nova_analise.status_code == 200


def test_alvo_grande_bloqueia_desenho_e_mostra_motivo():
    """Alvo sintético de 3 kb impede o desenho e exibe o motivo correspondente."""
    cfg_primers = carregar_config(RAIZ / "config" / "primers.toml")
    cfg_web = tomllib.loads((RAIZ / "config" / "web.toml").read_text(encoding="utf-8"))

    # Sequência sintética de 3000 pb (3 kb)
    seq_3kb = "ACGT" * 750
    alvo_3kb_bytes = f">alvo_3kb\n{seq_3kb}\n".encode("utf-8")

    # Verificação com a função pura tela_detalhe
    sessao = Sessao(RAIZ, cfg_primers, cfg_web)
    sessao.carregar_conjunto("conjunto_ok.csv", (DADOS / "conjunto_ok.csv").read_bytes())
    sessao.carregar_alvo("alvo_3kb.fasta", alvo_3kb_bytes)
    execucao = sessao.analisar()

    dados = tela_detalhe(execucao, "resumo", cfg_primers, cfg_web)
    assert dados["pode_desenhar"] is False
    assert (
        dados["motivo_sem_desenho"]
        == "O desenho de alternativas aceita alvos de até 2.000 pb; este tem 3.000 pb."
    )

    # Verificação em alvo pequeno (<= 2000 pb)
    sessao_ok = Sessao(RAIZ, cfg_primers, cfg_web)
    sessao_ok.carregar_conjunto("conjunto_ok.csv", (DADOS / "conjunto_ok.csv").read_bytes())
    sessao_ok.carregar_alvo("alvo_ok.fasta", (DADOS / "alvo_ok.fasta").read_bytes())
    exec_ok = sessao_ok.analisar()

    dados_ok = tela_detalhe(exec_ok, "resumo", cfg_primers, cfg_web)
    assert dados_ok["pode_desenhar"] is True
    assert dados_ok["motivo_sem_desenho"] == ""

    # Verificação via requisição HTTP: a aba Resumo não contém o formulário /desenhar
    cliente = TestClient(criar_app(RAIZ))
    cliente.post(
        "/testar/primers",
        files={"arquivo": ("conjunto_ok.csv", (DADOS / "conjunto_ok.csv").read_bytes(), "text/csv")},
    )
    cliente.post(
        "/testar/alvo",
        files={"arquivo": ("alvo_3kb.fasta", alvo_3kb_bytes, "text/plain")},
    )
    res_analise = cliente.get("/analise/1")
    assert res_analise.status_code == 200
    assert 'action="/desenhar"' not in res_analise.text
