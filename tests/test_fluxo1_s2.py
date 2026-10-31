"""Testes da spec s2: rotas novas, sessão por cookie e isolamento."""

from pathlib import Path

from starlette.testclient import TestClient

from interface.app import criar_app

RAIZ = Path(__file__).resolve().parents[1]
DADOS = RAIZ / "tests" / "dados"


def test_cookie_sessao_e_estaticos() -> None:
    app = criar_app(RAIZ)
    cliente = TestClient(app)

    res = cliente.get("/")
    assert res.status_code == 200
    set_cookie = res.headers.get("set-cookie", "")
    assert "lampsim_sessao=" in set_cookie
    set_cookie_lower = set_cookie.lower()
    assert "httponly" in set_cookie_lower
    assert "samesite=lax" in set_cookie_lower
    assert "path=/" in set_cookie_lower
    assert "max-age=7200" in set_cookie_lower

    res2 = cliente.get("/")
    assert res2.status_code == 200
    assert "set-cookie" not in res2.headers

    res_estatico = cliente.get("/static/estilo.css")
    assert res_estatico.status_code == 200
    assert "set-cookie" not in res_estatico.headers


def test_isolamento_entre_clientes() -> None:
    app = criar_app(RAIZ)
    cliente_a = TestClient(app)
    cliente_b = TestClient(app)

    # Cliente A carrega primers no novo fluxo
    res_a = cliente_a.post(
        "/testar/primers",
        files={"arquivo": ("conjunto_ok.csv", (DADOS / "conjunto_ok.csv").read_bytes(), "text/csv")},
        follow_redirects=False,
    )
    assert res_a.status_code == 303
    assert res_a.headers["location"] == "/testar/alvo"

    # Cliente B não tem conjunto, então /testar/alvo redireciona para /testar
    res_b = cliente_b.get("/testar/alvo", follow_redirects=False)
    assert res_b.status_code == 303
    assert res_b.headers["location"] == "/testar"

    # Erro de A em rota antiga não vaza para B
    res_err_a = cliente_a.post(
        "/conjunto",
        files={"arquivo": ("vazio.csv", b"", "text/csv")},
        follow_redirects=False,
    )
    assert res_err_a.status_code == 303
    assert res_err_a.headers["location"] == "/dados"

    res_dados_b = cliente_b.get("/dados")
    assert res_dados_b.status_code == 200
    assert "vazio.csv" not in res_dados_b.text

    res_dados_a = cliente_a.get("/dados")
    assert res_dados_a.status_code == 200
    assert "vazio.csv" in res_dados_a.text


def test_fluxo_passo_primers_validacoes() -> None:
    app = criar_app(RAIZ)
    cliente = TestClient(app)

    # Sem arquivo -> 400
    res_sem = cliente.post("/testar/primers")
    assert res_sem.status_code == 400
    assert "location" not in res_sem.headers
    assert "Escolha um arquivo antes de enviar." in res_sem.text
    assert 'role="alert"' in res_sem.text

    # Arquivo com nome vazio -> 400
    res_vazio = cliente.post("/testar/primers", files={"arquivo": ("", b"dados", "text/csv")})
    assert res_vazio.status_code == 400
    assert "Escolha um arquivo antes de enviar." in res_vazio.text

    # Arquivo inválido -> 400
    res_invalido = cliente.post(
        "/testar/primers",
        files={"arquivo": ("ruim.csv", b"isto nao e um csv", "text/csv")},
    )
    assert res_invalido.status_code == 400
    assert "location" not in res_invalido.headers
    assert "ruim.csv" in res_invalido.text

    # Arquivo válido -> 303 /testar/alvo
    res_ok = cliente.post(
        "/testar/primers",
        files={"arquivo": ("conjunto_ok.csv", (DADOS / "conjunto_ok.csv").read_bytes(), "text/csv")},
        follow_redirects=False,
    )
    assert res_ok.status_code == 303
    assert res_ok.headers["location"] == "/testar/alvo"

    res_alvo_page = cliente.get("/testar/alvo")
    assert res_alvo_page.status_code == 200
    assert "Passo 2 de 2" in res_alvo_page.text
    assert "conjunto_ok.csv" in res_alvo_page.text


def test_fluxo_passo_alvo_e_resultado() -> None:
    app = criar_app(RAIZ)
    cliente = TestClient(app)

    # Sem conjunto na sessão -> 303 /testar
    res_sem_conjunto = cliente.post(
        "/testar/alvo",
        files={"arquivo": ("alvo_ok.fasta", (DADOS / "alvo_ok.fasta").read_bytes(), "text/plain")},
        follow_redirects=False,
    )
    assert res_sem_conjunto.status_code == 303
    assert res_sem_conjunto.headers["location"] == "/testar"

    # Carrega conjunto
    cliente.post(
        "/testar/primers",
        files={"arquivo": ("conjunto_ok.csv", (DADOS / "conjunto_ok.csv").read_bytes(), "text/csv")},
        follow_redirects=False,
    )

    # POST alvo sem arquivo -> 400
    res_sem_arq = cliente.post("/testar/alvo")
    assert res_sem_arq.status_code == 400
    assert "Escolha um arquivo antes de enviar." in res_sem_arq.text

    # POST alvo com arquivo correto -> 303 /testar/resultado/1
    res_alvo = cliente.post(
        "/testar/alvo",
        files={"arquivo": ("alvo_ok.fasta", (DADOS / "alvo_ok.fasta").read_bytes(), "text/plain")},
        follow_redirects=False,
    )
    assert res_alvo.status_code == 303
    assert res_alvo.headers["location"] == "/testar/resultado/1"

    # GET /testar/resultado/1
    res_resultado = cliente.get("/testar/resultado/1")
    assert res_resultado.status_code == 200
    assert "conjunto_ok × alvo_ok" in res_resultado.text
    assert (
        "Estes primers devem amplificar este alvo, com ressalvas: 2 regras do guia não cumpridas."
        in res_resultado.text
    )
    assert 'id="tubos-dados"' in res_resultado.text
    assert "/analise/1/exportar" in res_resultado.text

    # GET /testar/resultado/99 -> 404
    assert cliente.get("/testar/resultado/99").status_code == 404

    # ID de desenho -> 404
    res_desenho = cliente.post("/desenhar", follow_redirects=False)
    assert res_desenho.status_code == 303
    assert res_desenho.headers["location"] == "/desenho/2"
    assert cliente.get("/testar/resultado/2").status_code == 404

    # Isolamento de resultado entre clientes: cliente_b não vê o resultado 1
    cliente_b = TestClient(app)
    assert cliente_b.get("/testar/resultado/1").status_code == 404

    # GET /dados -> 200 com título
    res_dados = cliente.get("/dados")
    assert res_dados.status_code == 200
    assert "<title>LAMP Primers — Início</title>" in res_dados.text
