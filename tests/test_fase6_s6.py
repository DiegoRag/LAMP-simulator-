"""Fase 6, sessão 6: as páginas mostram o KPI "a calibrar", a régua com zero, a fita e a base 1."""

from pathlib import Path

from fastapi.testclient import TestClient

from interface.app import criar_app

RAIZ = Path(__file__).resolve().parents[1]
DADOS = RAIZ / "tests" / "dados"
TEMPLATES = RAIZ / "src" / "interface" / "templates"
ALVO = (DADOS / "alvo_ok.fasta").read_text(encoding="utf-8").splitlines()[1]


def _analisar(conjunto: bytes, alvo: bytes | None = None) -> TestClient:
    cliente = TestClient(criar_app(RAIZ))
    cliente.post("/conjunto", files={"arquivo": ("c.csv", conjunto, "text/csv")}, follow_redirects=False)
    if alvo is not None:
        cliente.post("/alvo", files={"arquivo": ("a.fasta", alvo, "text/plain")}, follow_redirects=False)
    cliente.post("/analisar", follow_redirects=False)
    return cliente


def test_analise_mostra_fita_base_1_e_link_dos_tubos():
    cliente = _analisar((DADOS / "conjunto_ok.csv").read_bytes(), (DADOS / "alvo_ok.fasta").read_bytes())
    html = cliente.get("/analise/1").text
    assert "fita + · 216 pb · coordenadas a partir de 1 · F1, B2c, B3c e LF são lidas" in html
    assert '<a class="botao botao-secundario" href="/tubos?execucao=1">Ver nos tubos</a>' in html
    assert "e 6 de Tm a calibrar" in html
    assert "medida-tem-viol" in html  # 2 violações de ponta, fora as de Tm
    assert "regua-mini-zero" not in html  # nenhuma distância negativa


def test_regua_com_distancia_negativa_mostra_o_zero():
    conjunto = (DADOS / "conjunto_ok.csv").read_text(encoding="utf-8").replace(
        "GTGACAGACTCTCCTTAGATA", ALVO[20:41]
    )
    cliente = _analisar(conjunto.encode("utf-8"), (DADOS / "alvo_ok.fasta").read_bytes())
    html = cliente.get("/analise/1").text
    assert '<span class="regua-mini-zero" style="left: 15.241%"></span>' in html
    assert "<strong>−10 nt</strong>" in html


def test_sem_alvo_nao_ha_link_dos_tubos():
    cliente = _analisar((DADOS / "conjunto_ok.csv").read_bytes())
    assert "Ver nos tubos" not in cliente.get("/analise/1").text


def test_kpi_a_calibrar_tem_estilo_proprio():
    analise = (TEMPLATES / "analise.html").read_text(encoding="utf-8")
    assert "medida-tem-calibrar" in analise and "medida-calibrar" in analise
    css = (RAIZ / "src" / "interface" / "static" / "estilo.css").read_text(encoding="utf-8")
    for classe in (".medida-tem-calibrar::before", ".medida-calibrar", ".regua-mini-zero"):
        assert classe in css


def test_parametros_marcam_as_faixas_de_tm():
    cliente = _analisar((DADOS / "conjunto_ok.csv").read_bytes(), (DADOS / "alvo_ok.fasta").read_bytes())
    html = cliente.get("/parametros").text
    assert '<span class="pilula pilula-calibrar">faixas a calibrar</span>' in html
    assert 'faixa 59–61 °C · 200 nM na reação <span class="chip-marca">a calibrar</span>' in html
    assert "Em nt, na fita em que o conjunto foi localizado." in html
    assert "Em nt, na fita +." not in html


def test_desenho_mostra_base_1_e_nota_de_tm():
    cliente = TestClient(criar_app(RAIZ))
    cliente.post("/exemplo", data={"nome": "conjunto_ok — sintético"}, follow_redirects=False)
    cliente.post("/desenhar", follow_redirects=False)
    html = cliente.get("/desenho/1").text
    assert "fita + · 216 pb · coordenadas a partir de 1 · autoverificada" in html
    assert 'ainda estão a calibrar. <a href="/parametros">Ver em Parâmetros</a></p>' in html
