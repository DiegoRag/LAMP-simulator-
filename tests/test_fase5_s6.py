"""Fase 5, sessão 6: página Tubos alimentada pela simulação do conjunto analisado."""

import json
import re
from datetime import datetime
from pathlib import Path

from starlette.testclient import TestClient

from esteira.ponte import carregar_config, eixo_minutos, simular_tubos
from interface.app import criar_app
from interface.sessao import Execucao
from interface.telas import PADRAO_SERIE, tela_tubos, texto_eficiencia
from primers.analise import Resultado

RAIZ = Path(__file__).resolve().parents[1]
DADOS = RAIZ / "tests" / "dados"
CFG_ESTEIRA = carregar_config(RAIZ / "config" / "esteira.toml")
TROCA = {"A": "C", "C": "A", "G": "T", "T": "G"}


def _dados(texto: str) -> dict:
    m = re.search(r'<script type="application/json" id="tubos-dados">(.*?)</script>', texto, re.DOTALL)
    assert m is not None
    return json.loads(m.group(1))


def _cliente_com_analise(alvo_fasta: bytes) -> TestClient:
    cliente = TestClient(criar_app(RAIZ))
    conjunto = (DADOS / "conjunto_ok.csv").read_bytes()
    cliente.post("/conjunto", files={"arquivo": ("conjunto_ok.csv", conjunto, "text/csv")})
    cliente.post("/alvo", files={"arquivo": ("alvo.fasta", alvo_fasta, "text/plain")})
    res = cliente.post("/analisar", follow_redirects=False)
    assert res.status_code == 303 and res.headers["location"] == "/analise/1"
    return cliente


def _alvo_mutado(posicoes: list[int]) -> bytes:
    linhas = (DADOS / "alvo_ok.fasta").read_text(encoding="utf-8").split()
    bases = list("".join(linhas[1:]))
    for p in posicoes:
        bases[p] = TROCA[bases[p]]
    return (linhas[0] + "\n" + "".join(bases) + "\n").encode("utf-8")


def test_texto_eficiencia():
    assert texto_eficiencia(None) == "—"
    assert texto_eficiencia(1.0) == "100 %"
    assert texto_eficiencia(0.614125) == "61 %"
    assert texto_eficiencia(0.015) == "1,5 %"
    assert texto_eficiencia(0.3 ** 6) == "0,073 %"
    assert texto_eficiencia(0.0) == "0 %"


def test_tubos_vem_da_simulacao_do_conjunto_analisado():
    cliente = _cliente_com_analise((DADOS / "alvo_ok.fasta").read_bytes())
    res = cliente.get("/tubos")
    assert res.status_code == 200
    dados = _dados(res.text)
    assert dados["modo"] == "simulacao" and dados["aviso"] is None
    assert dados["origem"]["id"] == 1 and dados["origem"]["eficiencia"] == 1.0
    assert dados["origem"]["incompat_texto"] == "0 no total, 0 na ponta que estende"  # fase 6
    assert len(dados["t_min"]) == 121 and dados["t_min"][-1] == 60.0
    # as curvas embutidas são as do lampsim para a série de config/web.toml
    esperado = simular_tubos(1.0, PADRAO_SERIE, eixo_minutos(60.0, 0.5), 2026)
    assert [t["ph"] for t in dados["tubos"]] == [t["ph"] for t in esperado]
    assert [t["cores"] for t in dados["tubos"]] == [t["cores"] for t in esperado]
    assert [t["diagnostico"] for t in dados["tubos"]] == ["positivo"] * 7 + ["controle_ok"] * 2
    assert dados["tubos"][7]["t_positivo_texto"] == "não vira"
    assert res.text.count('data-tubo="') == 9
    assert "Leitura aos 15 min" in res.text
    assert 'href="/tubos?execucao=1&amp;modo=manual"' in res.text
    assert '<option value="60" selected>60×</option>' in res.text


def test_conjunto_com_seis_erros_na_ponta_nao_vira():
    # o exemplo da revisão, agora na página: antes as curvas eram as mesmas do conjunto perfeito
    cliente = _cliente_com_analise(_alvo_mutado([20, 48, 102, 108, 157, 187]))
    dados = _dados(cliente.get("/tubos").text)
    assert dados["modo"] == "simulacao"
    assert dados["origem"]["incompat_texto"] == "6 no total, 6 na ponta que estende"  # fase 6
    assert dados["origem"]["eficiencia_texto"] == "0,073 %"
    diags = [t["diagnostico"] for t in dados["tubos"]]
    assert diags == ["negativo"] * 7 + ["controle_ok"] * 2
    assert all(t["t_positivo_min"] is None for t in dados["tubos"][:8])


def test_regiao_essencial_ausente_avisa():
    cliente = _cliente_com_analise(_alvo_mutado([31, 33, 35, 37, 39]))
    dados = _dados(cliente.get("/tubos").text)
    assert dados["modo"] == "simulacao" and dados["origem"]["eficiencia"] == 0.0
    assert dados["aviso"].startswith("Região essencial sem posição no alvo (F2)")


def test_modo_manual_opcional_com_aviso():
    cliente = _cliente_com_analise((DADOS / "alvo_ok.fasta").read_bytes())
    res = cliente.get("/tubos?modo=manual")
    dados = _dados(res.text)
    assert dados["modo"] == "manual" and dados["aviso"].startswith("Curvas manuais")
    assert 'href="/tubos?execucao=1"' in res.text  # volta aos tubos simulados
    assert "data-tubos-editor" in res.text
    assert dados["tubos"][1]["inicio"] == 24.0  # minutos, não segundos


def test_execucao_inexistente_ou_de_desenho_da_404():
    cliente = TestClient(criar_app(RAIZ))
    assert cliente.get("/tubos?execucao=99").status_code == 404


def test_analise_sem_alvo_cai_no_manual():
    execucao = Execucao(
        id=5,
        tipo="analise",
        criada_em=datetime(2026, 10, 18, 12, 0),
        rotulo_conjunto="x",
        rotulo_alvo="—",
        alvo=None,
        leitura=None,
        resultado=Resultado(linhas=[], heterodimeros={}, violacoes=[], avisos=[], regioes=None),
    )
    dados = tela_tubos({}, CFG_ESTEIRA, execucao)
    assert dados["modo"] == "manual" and dados["pode_simular"] is False
    assert dados["aviso"].startswith("A análise 5 não tem alvo")
    assert dados["url_manual"] == "/tubos?execucao=5&modo=manual"
