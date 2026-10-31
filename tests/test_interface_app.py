"""Testes de ponta a ponta da interface web local do lamp-sim."""

import io
import json
import re
import tomllib
import zipfile
from pathlib import Path

from starlette.testclient import TestClient

import primers.desenho
from interface.app import criar_app

RAIZ = Path(__file__).resolve().parents[1]
DADOS = RAIZ / "tests" / "dados"


def test_inicio() -> None:
    cliente = TestClient(criar_app(RAIZ))
    res = cliente.get("/")
    assert res.status_code == 200
    assert "Exames genéticos em qualquer lugar, a qualquer hora." in res.text

    res_dados = cliente.get("/dados")
    assert res_dados.status_code == 200
    assert "<title>LAMP Primers — Início</title>" in res_dados.text


def test_fluxo_analise() -> None:
    cliente = TestClient(criar_app(RAIZ))
    res = cliente.post(
        "/conjunto",
        files={"arquivo": ("conjunto_ok.csv", (DADOS / "conjunto_ok.csv").read_bytes(), "text/csv")},
        follow_redirects=False,
    )
    assert res.status_code == 303
    assert res.headers["location"] == "/dados"

    res = cliente.get("/dados")
    assert "conjunto_ok.csv" in res.text

    res = cliente.post(
        "/alvo",
        files={"arquivo": ("alvo_ok.fasta", (DADOS / "alvo_ok.fasta").read_bytes(), "text/plain")},
        follow_redirects=False,
    )
    assert res.status_code == 303
    assert res.headers["location"] == "/dados"

    res = cliente.get("/dados")
    assert "alvo_ok.fasta" in res.text

    res = cliente.post("/analisar", follow_redirects=False)
    assert res.status_code == 303
    assert res.headers["location"] == "/analise/1"

    res = cliente.get("/analise/1")
    assert res.status_code == 200
    assert "<title>LampSim · Análise: conjunto_ok × alvo_ok</title>" in res.text

    res_primers = cliente.get("/analise/1/primers")
    assert res_primers.status_code == 200
    assert "GTGACAGACTCTCCTTAGATA" in res_primers.text
    assert "GTAAACCCCCAAGTGTTATTGG" in res_primers.text
    assert "CGCATATGGATTCTGTAC" in res_primers.text

    res = cliente.get("/analise", follow_redirects=False)
    assert res.status_code == 303
    assert res.headers["location"] == "/analise/1"

    res = cliente.get("/analise/1?tm_calibracao=0")
    assert res.status_code == 200


def test_exportar_analise() -> None:
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

    res = cliente.get("/analise/1/exportar")
    assert res.status_code == 200
    assert res.headers["content-type"] == "application/zip"
    assert "analise_1.zip" in res.headers["content-disposition"]
    with zipfile.ZipFile(io.BytesIO(res.content)) as zf:
        # fase 6: com conjunto completo e alvo, o zip leva o selecao.json (contrato v2)
        assert zf.namelist() == [
            "primers.csv",
            "heterodimeros.csv",
            "violacoes.csv",
            "resultado.json",
            "selecao.json",
        ]


def test_fluxo_desenho() -> None:
    cliente = TestClient(criar_app(RAIZ))
    res = cliente.post("/exemplo", data={"nome": "conjunto_ok — sintético"}, follow_redirects=False)
    assert res.status_code == 303
    assert res.headers["location"] == "/dados"

    res = cliente.post("/desenhar", follow_redirects=False)
    assert res.status_code == 303
    assert res.headers["location"] == "/desenho/1"

    res = cliente.get("/desenho/1")
    assert res.status_code == 200
    assert "<title>LampSim · Desenho: alvo_ok</title>" in res.text
    assert primers.desenho.AVISO in res.text
    assert "TGACAGACTCTCCTTAGATAGCTAA" in res.text
    assert "ACCCCCAAGTGTTATTGGCAAGTT" in res.text

    res = cliente.get("/desenho/1?sel=2")
    assert res.status_code == 200
    assert "CTCCGTCCACCAGTGTGAATCG" in res.text  # fase 4: a proposta 2 vem de outro núcleo F2–B2c

    res = cliente.get("/desenho/1?agrupar=1")
    assert res.status_code == 200

    res = cliente.get("/desenho/1/exportar")
    assert res.status_code == 200
    with zipfile.ZipFile(io.BytesIO(res.content)) as zf:
        assert zf.namelist() == ["conjuntos.csv", "resultado.json"]

    res = cliente.post("/desenho/1/analisar", data={"posicao": "1"}, follow_redirects=False)
    assert res.status_code == 303
    assert res.headers["location"] == "/analise/2"

    res = cliente.get("/analise/2")
    assert res.status_code == 200


def test_mensagens_uma_vez() -> None:
    cliente = TestClient(criar_app(RAIZ))
    res = cliente.post("/analisar", follow_redirects=False)
    assert res.status_code == 303
    assert res.headers["location"] == "/dados"

    res = cliente.get("/dados")
    assert "Carregue um conjunto de primers antes de analisar." in res.text

    res = cliente.get("/dados")
    assert "Carregue um conjunto de primers antes de analisar." not in res.text

    res = cliente.post("/desenhar", follow_redirects=False)
    assert res.status_code == 303
    assert res.headers["location"] == "/dados"

    res = cliente.get("/dados")
    assert "Carregue a sequência-alvo antes de desenhar." in res.text

    res = cliente.post(
        "/conjunto",
        files={"arquivo": ("x.csv", b"nome,seq\nF3,ACGT\n", "text/csv")},
        follow_redirects=False,
    )
    assert res.status_code == 303
    assert res.headers["location"] == "/dados"

    res = cliente.get("/dados")
    assert "x.csv: Cabeçalho inválido na linha 1" in res.text


def test_nao_encontrado() -> None:
    cliente = TestClient(criar_app(RAIZ))
    res = cliente.get("/analise/99")
    assert res.status_code == 404
    assert "<title>LAMP Primers — Não encontrado</title>" in res.text
    assert "Execução não encontrada" in res.text
    assert cliente.get("/desenho/99").status_code == 404

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

    assert cliente.get("/desenho/1").status_code == 404
    assert cliente.post("/desenho/99/analisar", data={"posicao": "1"}).status_code == 404


def test_nao_encontrado_estatico() -> None:
    cliente = TestClient(criar_app(RAIZ))
    res = cliente.get("/static/nao-existe.css")
    assert res.status_code == 404
    assert "O endereço pedido não existe" in res.text


def test_parametros() -> None:
    cliente = TestClient(criar_app(RAIZ))
    res = cliente.get("/parametros")
    assert res.status_code == 200
    assert (
        "<title>LampSim · Parâmetros</title>" in res.text
        or "<title>LAMP Primers — Parâmetros</title>" in res.text
    )
    assert "config/primers.toml" in res.text


def test_estaticos() -> None:
    cliente = TestClient(criar_app(RAIZ))
    assert cliente.get("/static/estilo.css").status_code == 200
    assert cliente.get("/static/app.js").status_code == 200


def test_docs_desativados() -> None:
    cliente = TestClient(criar_app(RAIZ))
    assert cliente.get("/docs").status_code == 404


def test_config_web() -> None:
    with open(RAIZ / "config" / "web.toml", "rb") as f:
        cfg = tomllib.load(f)
    assert cfg["servidor"]["host"] == "127.0.0.1"


CABECALHO = {"X-Requested-With": "lamp-carregamento"}


def test_destino_em_json() -> None:
    """Verifica que respostas POST com cabeçalho de carregamento retornam destino em JSON."""
    cliente = TestClient(criar_app(RAIZ))
    res = cliente.post("/exemplo", data={"nome": "conjunto_ok — sintético"}, headers=CABECALHO)
    assert res.status_code == 200
    assert res.json() == {"destino": "/dados"}

    res = cliente.post("/analisar", headers=CABECALHO)
    assert res.status_code == 200
    assert res.json() == {"destino": "/analise/1"}

    res = cliente.post("/desenhar", headers=CABECALHO)
    assert res.status_code == 200
    assert res.json() == {"destino": "/desenho/2"}

    res = cliente.post("/desenho/2/analisar", data={"posicao": "1"}, headers=CABECALHO)
    assert res.status_code == 200
    assert res.json() == {"destino": "/analise/3"}


def test_destino_em_json_com_erro() -> None:
    """Verifica mensagem de erro ao requisitar destino em JSON sem dados carregados."""
    cliente = TestClient(criar_app(RAIZ))
    res = cliente.post("/analisar", headers=CABECALHO)
    assert res.status_code == 200
    assert res.json() == {"destino": "/dados"}

    res = cliente.get("/dados")
    assert "Carregue um conjunto de primers antes de analisar." in res.text

    res = cliente.get("/dados")
    assert "Carregue um conjunto de primers antes de analisar." not in res.text


def test_envio_de_arquivo_em_json() -> None:
    """Verifica upload de conjunto com cabeçalho de carregamento retornando destino em JSON."""
    cliente = TestClient(criar_app(RAIZ))
    res = cliente.post(
        "/conjunto",
        files={"arquivo": ("conjunto_ok.csv", (DADOS / "conjunto_ok.csv").read_bytes(), "text/csv")},
        headers=CABECALHO,
    )
    assert res.status_code == 200
    assert res.json() == {"destino": "/dados"}

    res = cliente.get("/dados")
    assert "conjunto_ok.csv" in res.text


def test_estaticos_novos() -> None:
    """Verifica disponibilidade dos novos arquivos estáticos."""
    cliente = TestClient(criar_app(RAIZ))
    assert cliente.get("/static/helice.js").status_code == 200
    assert cliente.get("/static/fontes/open-sans-latin.woff2").status_code == 200


def test_estaticos_revalidados() -> None:
    """Verifica que respostas de estáticos trazem Cache-Control: no-cache."""
    cliente = TestClient(criar_app(RAIZ))
    res = cliente.get("/static/estilo.css")
    assert res.status_code == 200
    assert res.headers["cache-control"] == "no-cache"


def test_pagina_tubos() -> None:
    """Verifica resposta 200 e elementos obrigatórios da página /tubos."""
    cliente = TestClient(criar_app(RAIZ))
    res = cliente.get("/tubos")
    assert res.status_code == 200
    assert "data-tubos-bandeja" in res.text
    assert res.text.count('data-tubo="') == 9
    assert 'id="tubos-dados"' in res.text
    assert "/static/tubos.js" in res.text

    m = re.search(r'<script type="application/json" id="tubos-dados">(.*?)</script>', res.text, re.DOTALL)
    assert m is not None
    dados = json.loads(m.group(1))
    assert len(dados["tubos"]) == 9
    # fase 5: sem análise, a página mostra as curvas manuais, em minutos e com aviso
    assert dados["modo"] == "manual"
    assert "data-tubos-aviso" in res.text
    assert "0,0 min" in res.text


def test_links_navegacao() -> None:
    """Verifica links de navegação da barra superior sem link para /tubos."""
    cliente = TestClient(criar_app(RAIZ))
    for rota in ("/dados", "/tubos"):
        res = cliente.get(rota)
        assert res.status_code == 200
        assert '<a class="link-topo" href="/testar">Testar primers</a>' in res.text
        assert '<a class="link-topo" href="/analise">Análise</a>' in res.text
        assert '<a class="link-topo" href="/">Início</a>' in res.text
        assert '<a class="link-topo" href="/parametros">Parâmetros</a>' in res.text
        assert 'class="link-topo" href="/tubos"' not in res.text

    res_param = cliente.get("/parametros")
    assert res_param.status_code == 200
    assert (
        '<a class="marca-nav-link" href="/parametros" aria-current="page">Parâmetros</a>'
        in res_param.text
    )
