"""Testes do tempo limite de análise no fluxo guiado (spec s2 do conserto)."""

import tomllib
from pathlib import Path

import pytest
from starlette.testclient import TestClient

import primers.entrada
from interface.analise_isolada import AnalisadorIsolado, TempoEsgotado, tempo_max_s
from interface.app import criar_app
from interface.sessao import Execucao

RAIZ = Path(__file__).resolve().parents[1]
DADOS = RAIZ / "tests" / "dados"
TEMPO_CURTO_S = 0.001
TEMPO_LONGO_S = 30.0


def test_tempo_max_s_config() -> None:
    """Valida a leitura de tempo_max_s com fallbacks e o arquivo web.toml real."""
    assert tempo_max_s({}) == 45.0
    assert tempo_max_s({"analise": {"tempo_max_s": 1.5}}) == 1.5
    assert tempo_max_s({"analise": {"tempo_max_s": True}}) == 45.0
    assert tempo_max_s({"analise": {"tempo_max_s": 0}}) == 45.0
    assert tempo_max_s({"analise": {"tempo_max_s": -1}}) == 45.0
    assert tempo_max_s({"analise": {"tempo_max_s": "3"}}) == 45.0

    cfg_real = tomllib.loads((RAIZ / "config" / "web.toml").read_text(encoding="utf-8"))
    assert tempo_max_s(cfg_real) == 45.0


def test_analisador_isolado_fluxo() -> None:
    """Testa o AnalisadorIsolado diretamente com análise normal, erro e timeout."""
    cfg_primers = primers.entrada.carregar_config(RAIZ / "config" / "primers.toml")
    cfg_web = tomllib.loads((RAIZ / "config" / "web.toml").read_text(encoding="utf-8"))

    conjunto_bytes = (DADOS / "conjunto_ok.csv").read_bytes()
    alvo_bytes = (DADOS / "alvo_ok.fasta").read_bytes()

    analisador = AnalisadorIsolado(RAIZ, cfg_primers, cfg_web, tempo_max=TEMPO_LONGO_S)
    try:
        # Sucesso inicial dentro do limite de 30 s
        execucao = analisador.analisar("conjunto_ok.csv", conjunto_bytes, "alvo_ok.fasta", alvo_bytes)
        assert isinstance(execucao, Execucao)
        assert execucao.tipo == "analise"

        # FASTA inválido levanta ValueError comum, não TempoEsgotado
        with pytest.raises(ValueError) as excinfo:
            analisador.analisar(
                "conjunto_ok.csv",
                conjunto_bytes,
                "invalido.fasta",
                b"isto nao e fasta",
            )
        assert not isinstance(excinfo.value, TempoEsgotado)

        # Limite curto estoura e levanta TempoEsgotado
        analisador.tempo_max = TEMPO_CURTO_S
        with pytest.raises(TempoEsgotado):
            analisador.analisar("conjunto_ok.csv", conjunto_bytes, "alvo_ok.fasta", alvo_bytes)

        # Voltando o limite para 30 s, o processo recriado atende normalmente
        analisador.tempo_max = TEMPO_LONGO_S
        execucao_recriada = analisador.analisar(
            "conjunto_ok.csv",
            conjunto_bytes,
            "alvo_ok.fasta",
            alvo_bytes,
        )
        assert isinstance(execucao_recriada, Execucao)
        assert execucao_recriada.tipo == "analise"
    finally:
        analisador.encerrar()


def test_analise_tempo_limite_app() -> None:
    """Testa a interrupção por limite de tempo integrada à rota POST /testar/alvo."""
    app = criar_app(RAIZ)
    cliente = TestClient(app)
    try:
        conjunto_bytes = (DADOS / "conjunto_ok.csv").read_bytes()
        alvo_bytes = (DADOS / "alvo_ok.fasta").read_bytes()

        # Passo 1: carrega primers com sucesso
        res1 = cliente.post(
            "/testar/primers",
            files={"arquivo": ("conjunto_ok.csv", conjunto_bytes, "text/csv")},
            follow_redirects=False,
        )
        assert res1.status_code == 303
        assert res1.headers["location"] == "/testar/alvo"

        # Configura limite mínimo para forçar timeout
        app.state.analisador.tempo_max = TEMPO_CURTO_S

        # Passo 2: envio do alvo deve falhar com 400 e mensagem de interrupção
        res2 = cliente.post(
            "/testar/alvo",
            files={"arquivo": ("alvo_ok.fasta", alvo_bytes, "text/plain")},
            follow_redirects=False,
        )
        assert res2.status_code == 400
        assert "location" not in res2.headers
        assert 'role="alert"' in res2.text
        assert "Passo 2 de 2" in res2.text
        assert "foi interrompida" in res2.text

        # Volta o limite para 30 s e reenvia o alvo
        app.state.analisador.tempo_max = TEMPO_LONGO_S
        res3 = cliente.post(
            "/testar/alvo",
            files={"arquivo": ("alvo_ok.fasta", alvo_bytes, "text/plain")},
            follow_redirects=False,
        )
        # O id começa em 1 pois a execução interrompida não foi registrada na sessão
        assert res3.status_code == 303
        assert res3.headers["location"] == "/testar/resultado/1"
    finally:
        app.state.analisador.encerrar()
