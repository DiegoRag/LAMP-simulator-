"""Testes dos nomes e detalhes de tubos em linguagem simples (conserto s1)."""

import tomllib
from pathlib import Path

from starlette.testclient import TestClient

import esteira.ponte
import primers.entrada
from interface.app import criar_app
from interface.sessao import Sessao
from interface.telas import _serie_simulada, tela_tubos

RAIZ = Path(__file__).resolve().parents[1]
DADOS = RAIZ / "tests" / "dados"


def test_serie_simulada_detalhes() -> None:
    """Valida detalhes na série simulada padrão e com item customizado."""
    tubos_padrao = _serie_simulada({})
    assert len(tubos_padrao) == 9
    assert all(bool(t.get("detalhe")) for t in tubos_padrao)

    tubos_custom = _serie_simulada({"serie": [{"nome": "X", "copias": 5}]})
    assert len(tubos_custom) == 1
    assert tubos_custom[0]["detalhe"] == ""


def test_tela_tubos_detalhe_copiado() -> None:
    """Valida cópia de detalhes para os tubos simulados em tela_tubos."""
    cfg_primers = primers.entrada.carregar_config(RAIZ / "config" / "primers.toml")
    cfg_esteira = esteira.ponte.carregar_config(RAIZ / "config" / "esteira.toml")
    cfg_web = tomllib.loads((RAIZ / "config" / "web.toml").read_text(encoding="utf-8"))

    sessao = Sessao(RAIZ, cfg_primers, cfg_web)
    sessao.carregar_conjunto("conjunto_ok.csv", (DADOS / "conjunto_ok.csv").read_bytes())
    sessao.carregar_alvo("alvo_ok.fasta", (DADOS / "alvo_ok.fasta").read_bytes())
    execucao = sessao.analisar()

    dados_tubos = tela_tubos(cfg_web, cfg_esteira, execucao, "simulacao")
    assert dados_tubos["modo"] == "simulacao"
    assert len(dados_tubos["tubos"]) == 9
    assert all("detalhe" in t and bool(t["detalhe"]) for t in dados_tubos["tubos"])
    assert dados_tubos["tubos"][0]["nome"] == "Carga muito alta"
    assert dados_tubos["tubos"][0]["detalhe"] == "1 milhão de cópias do gene"


def test_resultado_rota_tubos_nomes() -> None:
    """Valida a rota /testar/resultado/1 com nomes para leigos e ausência de 10⁶."""
    app = criar_app(RAIZ)
    cliente = TestClient(app)
    try:
        conjunto_bytes = (DADOS / "conjunto_ok.csv").read_bytes()
        alvo_bytes = (DADOS / "alvo_ok.fasta").read_bytes()

        res_p1 = cliente.post(
            "/testar/primers",
            files={"arquivo": ("conjunto_ok.csv", conjunto_bytes, "text/csv")},
            follow_redirects=False,
        )
        assert res_p1.status_code == 303
        assert res_p1.headers["location"] == "/testar/alvo"

        res_p2 = cliente.post(
            "/testar/alvo",
            files={"arquivo": ("alvo_ok.fasta", alvo_bytes, "text/plain")},
            follow_redirects=False,
        )
        assert res_p2.status_code == 303
        assert res_p2.headers["location"] == "/testar/resultado/1"

        res = cliente.get("/testar/resultado/1")
        assert res.status_code == 200
        html = res.text

        assert "Carga muito alta" in html
        assert "1 milhão de cópias do gene" in html
        assert "Controle sem DNA" in html
        assert "deve ficar rosa" in html
        assert "10⁶" not in html
    finally:
        app.state.analisador.encerrar()
