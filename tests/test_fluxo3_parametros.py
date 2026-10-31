"""Testes da tela de Parâmetros na identidade LampSim, em abas com URL própria (fase 3)."""

from pathlib import Path

from starlette.testclient import TestClient

from interface.app import criar_app
from interface.sessao import Sessao
from interface.telas_fluxo import tela_parametros_marca
from primers.entrada import carregar_config

RAIZ = Path(__file__).resolve().parents[1]
DADOS = RAIZ / "tests" / "dados"
CFG_PRIMERS = carregar_config(RAIZ / "config" / "primers.toml")
CFG_WEB = {
    "entrada": {"tamanho_max_kb": 64, "alvo_desenho_max_pb": 2000},
    "analise": {"tm_como_calibracao": True},
    "heterodimeros": {"escala_kcal": [-10.0, -8.0, -6.0, -4.0]},
    "exemplos": [],
}

URLS_PARAMETROS = (
    ("calibracao", "/parametros", "Calibração de Tm"),
    ("reacao", "/parametros/reacao", "Condições de reação"),
    ("distancias", "/parametros/distancias", "Distâncias"),
    ("regras", "/parametros/regras", "Regras por primer"),
    ("desenho", "/parametros/desenho", "Desenho"),
)


def test_tela_parametros_marca_estrutura():
    """Verifica contrato da função pura tela_parametros_marca."""
    sessao = Sessao(RAIZ, CFG_PRIMERS, CFG_WEB)
    for chave_alvo, _, _ in URLS_PARAMETROS:
        dados = tela_parametros_marca(sessao, CFG_PRIMERS, chave_alvo)
        assert dados["aba"] == chave_alvo
        assert dados["url_testar"] == "/testar"
        assert len(dados["abas"]) == 5
        for item in dados["abas"]:
            assert set(item.keys()) == {"chave", "rotulo", "href", "atual"}
            if item["chave"] == chave_alvo:
                assert item["atual"] is True
            else:
                assert item["atual"] is False
            if item["chave"] == "calibracao":
                assert item["href"] == "/parametros"
            else:
                assert item["href"] == f"/parametros/{item['chave']}"


def test_cinco_urls_status_200_body_marca_e_aria_current():
    """As 5 URLs respondem 200 com body marca e aria-current na barra e na aba ativa."""
    cliente = TestClient(criar_app(RAIZ))
    for chave_atual, url, _ in URLS_PARAMETROS:
        res = cliente.get(url)
        assert res.status_code == 200
        assert 'body class="marca' in res.text
        assert (
            '<a class="marca-nav-link" href="/parametros" aria-current="page">Parâmetros</a>'
            in res.text
        )

        for chave_aba, url_aba, _ in URLS_PARAMETROS:
            trecho = f'href="{url_aba}" aria-current="page"'
            if chave_aba == chave_atual:
                if chave_aba == "calibracao":
                    # Na aba calibracao, a URL /parametros aparece na barra e na aba
                    assert res.text.count(trecho) == 2
                else:
                    assert res.text.count(trecho) == 1
            else:
                if chave_aba != "calibracao":
                    assert trecho not in res.text

        # Exatamente duas marcações de aria-current="page" no documento inteiro:
        # uma na barra principal ("Parâmetros") e uma na aba ativa.
        assert res.text.count('aria-current="page"') == 2


def test_chave_desconhecida_retorna_404():
    """Acesso a aba inexistente em /parametros/{aba} devolve 404."""
    cliente = TestClient(criar_app(RAIZ))
    res = cliente.get("/parametros/xyz")
    assert res.status_code == 404
    assert cliente.get("/parametros/nao_existe").status_code == 404


def test_calibracao_sem_analise_mostra_link_testar():
    """Sem análise na sessão, a aba Calibração aponta para /testar."""
    cliente = TestClient(criar_app(RAIZ))
    res = cliente.get("/parametros")
    assert res.status_code == 200
    assert 'href="/testar"' in res.text
    assert "Testar primers" in res.text


def test_calibracao_e_distancias_com_analise_do_fluxo_testar():
    """Após analisar pelo fluxo /testar, mostra o gráfico em Calibração e referência em Distâncias."""
    cliente = TestClient(criar_app(RAIZ))

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
    assert "/testar/resultado/" in res_alvo.headers["location"]

    # Aba Calibração com análise: gráfico cal-grafico presente
    res_cal = cliente.get("/parametros")
    assert res_cal.status_code == 200
    assert "cal-grafico" in res_cal.text

    # Aba Distâncias com análise: coluna da referência com valores medidos
    res_dist = cliente.get("/parametros/distancias")
    assert res_dist.status_code == 200
    assert "conjunto_ok × alvo_ok" in res_dist.text
    assert "146" in res_dist.text
