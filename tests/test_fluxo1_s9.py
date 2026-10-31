"""Testes de ponta a ponta do fluxo guiado da fase 1 (spec s9)."""

import json
import re
from pathlib import Path

from starlette.testclient import TestClient

from interface.app import criar_app

RAIZ = Path(__file__).resolve().parents[1]
DADOS = RAIZ / "tests" / "dados"


def test_fluxo_guiado_ponta_a_ponta() -> None:
    cliente = TestClient(criar_app(RAIZ))

    # 1. GET / -> landing page
    res = cliente.get("/")
    assert res.status_code == 200
    assert 'href="/testar"' in res.text
    assert "Exames genéticos em qualquer lugar, a qualquer hora." in res.text
    assert "/static/marca/lampsim_horizontal_colorido.svg" in res.text
    assert "lampsim_completo_negativo.svg" not in res.text
    assert "/static/marca.css" in res.text
    assert "/static/fluxo.js" in res.text
    assert '<a class="marca-nav-link" href="/" aria-current="page">Início</a>' in res.text

    # 2. GET /testar -> passo 1
    res = cliente.get("/testar")
    assert res.status_code == 200
    assert "Passo 1 de 2" in res.text
    assert "data-envio-automatico" in res.text
    assert 'accept=".csv,.fasta,.fa,.fna"' in res.text
    assert 'name="arquivo"' in res.text
    assert '<a class="marca-nav-link" href="/testar" aria-current="page">Testar primers</a>' in res.text
    assert 'role="alert"' not in res.text

    # 3. POST /testar/primers com conjunto_ok.csv
    conjunto_bytes = (DADOS / "conjunto_ok.csv").read_bytes()
    res = cliente.post(
        "/testar/primers",
        files={"arquivo": ("conjunto_ok.csv", conjunto_bytes, "text/csv")},
        follow_redirects=False,
    )
    assert res.status_code == 303
    assert res.headers["location"] == "/testar/alvo"

    # 4. GET /testar/alvo -> passo 2
    res = cliente.get("/testar/alvo")
    assert res.status_code == 200
    assert "Passo 2 de 2" in res.text
    assert "conjunto_ok.csv" in res.text
    assert 'data-espera="fluxo"' in res.text
    assert "data-espera-fluxo" in res.text
    assert 'accept=".fasta,.fa,.fna"' in res.text

    # 5. POST /testar/alvo com alvo_ok.fasta
    alvo_bytes = (DADOS / "alvo_ok.fasta").read_bytes()
    res = cliente.post(
        "/testar/alvo",
        files={"arquivo": ("alvo_ok.fasta", alvo_bytes, "text/plain")},
        follow_redirects=False,
    )
    assert res.status_code == 303
    location = res.headers["location"]
    assert location.startswith("/testar/resultado/")
    id_execucao = location.rsplit("/", 1)[-1]

    # 6. GET /testar/resultado/{id}
    res = cliente.get(location)
    assert res.status_code == 200
    assert "conjunto_ok × alvo_ok" in res.text
    frase_esperada = (
        "Estes primers devem amplificar este alvo, com ressalvas: "
        "2 regras do guia não cumpridas."
    )
    assert frase_esperada in res.text
    assert "Na simulação, detectam a partir de 1 cópia." in res.text
    assert "6 de 6" in res.text
    assert "100 %" in res.text
    assert "veredito-alerta" in res.text
    assert "/static/marca-resultado.css" in res.text
    assert "/static/tubos.js" in res.text
    assert res.text.count('data-tubo="') == 9

    m = re.search(r'<script type="application/json" id="tubos-dados">(.*?)</script>', res.text, re.DOTALL)
    assert m is not None
    dados_tubos = json.loads(m.group(1))
    assert dados_tubos["modo"] == "simulacao"
    assert len(dados_tubos["tubos"]) == 9

    # 7. Links do resultado
    res_analise = cliente.get(f"/analise/{id_execucao}")
    assert res_analise.status_code == 200

    res_exportar = cliente.get(f"/analise/{id_execucao}/exportar")
    assert res_exportar.status_code == 200
    assert "application/zip" in res_exportar.headers.get("content-type", "")

    # 8. Fluxo de erro do passo 2
    res_erro = cliente.post(
        "/testar/alvo",
        files={"arquivo": ("alvo_ruim.fasta", b"", "text/plain")},
        follow_redirects=False,
    )
    assert res_erro.status_code == 400
    assert "location" not in res_erro.headers
    assert 'role="alert"' in res_erro.text
    assert "Passo 2 de 2" in res_erro.text

    # 9. Navbar antiga em /dados
    res_dados = cliente.get("/dados")
    assert res_dados.status_code == 200
    assert 'class="marca-nav-antiga"' in res_dados.text
    assert 'href="/testar"' in res_dados.text
    assert 'class="trilho"' not in res_dados.text
