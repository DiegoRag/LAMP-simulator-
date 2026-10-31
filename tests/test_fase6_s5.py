"""Fase 6, sessão 5: a exportação da Análise leva o resultado.json completo e o selecao.json (v2)."""

import io
import json
import zipfile
from pathlib import Path

from fastapi.testclient import TestClient

from esteira.contrato import validar_selecao
from esteira.ponte import eixo_minutos, simular_tubos, tubos_da_selecao
from interface.app import criar_app, selecao_da_execucao
from primers.desenho import AVISO

RAIZ = Path(__file__).resolve().parents[1]
DADOS = RAIZ / "tests" / "dados"


def _enviar(cliente: TestClient, alvo: bool = True) -> None:
    cliente.post(
        "/conjunto",
        files={"arquivo": ("conjunto_ok.csv", (DADOS / "conjunto_ok.csv").read_bytes(), "text/csv")},
        follow_redirects=False,
    )
    if alvo:
        cliente.post(
            "/alvo",
            files={"arquivo": ("alvo_ok.fasta", (DADOS / "alvo_ok.fasta").read_bytes(), "text/plain")},
            follow_redirects=False,
        )


def _zip(cliente: TestClient, url: str) -> dict:
    res = cliente.get(url)
    assert res.status_code == 200
    with zipfile.ZipFile(io.BytesIO(res.content)) as zf:
        return {nome: zf.read(nome).decode("utf-8") for nome in zf.namelist()}


def test_exportacao_traz_selecao_v2_e_resultado_completo():
    cliente = TestClient(criar_app(RAIZ))
    _enviar(cliente)
    cliente.post("/analisar", follow_redirects=False)
    arquivos = _zip(cliente, "/analise/1/exportar")
    assert list(arquivos) == [
        "primers.csv",
        "heterodimeros.csv",
        "violacoes.csv",
        "resultado.json",
        "selecao.json",
    ]

    selecao = json.loads(arquivos["selecao.json"])
    assert validar_selecao(selecao) == []
    assert selecao["esquema"] == "lamp-sim/selecao@2"
    assert selecao["gerado_por"]["experimento"] == "interface/analise/1"
    conjunto = selecao["conjuntos"][0]
    assert (conjunto["conjunto_id"], conjunto["origem"]) == ("conjunto_ok", "usuario")
    assert conjunto["selecionado"] is True
    assert "aviso" not in conjunto
    assert conjunto["primers"]["FIP"] == "GTAAACCCCCAAGTGTTATTGGCGCATATGGATTCTGTAC"
    assert conjunto["primers"]["ligante"] == ""
    medida = conjunto["por_alelo"][0]
    assert (medida["alelo_id"], medida["acesso"], medida["f_previsto"]) == ("alvo_sintetico_ok", None, 1.0)
    assert (medida["incompat_total"], medida["incompat_3p"]) == (0, 0)

    resultado = json.loads(arquivos["resultado.json"])
    assert resultado["alvo"]["id"] == "alvo_sintetico_ok"
    assert resultado["alvo"]["arquivo"] == "alvo_ok.fasta"
    assert resultado["fita"] == "+"
    assert resultado["conjunto"]["FIP"]["origem"] == "partes informadas"
    assert resultado["eficiencia_prevista"] == {
        "f_previsto": 1.0,
        "modelo": "esteira.ponte.eficiencia_prevista com config/esteira.toml (provisório)",
        "parametros": {"fator_ponta3": 0.3, "fator_outra": 0.85, "piso_laco": 0.5},
    }

    # o lado leitor do contrato simula os mesmos tubos que a eficiência daria direto
    serie = [{"nome": "10³ cópias", "copias": 1000.0, "papel": "amostra"}]
    t = eixo_minutos(30.0, 0.5)
    assert tubos_da_selecao(selecao, serie, t, 7) == simular_tubos(1.0, serie, t, 7)


def test_sem_alvo_nao_ha_selecao():
    cliente = TestClient(criar_app(RAIZ))
    _enviar(cliente, alvo=False)
    cliente.post("/analisar", follow_redirects=False)
    arquivos = _zip(cliente, "/analise/1/exportar")
    assert list(arquivos) == ["primers.csv", "heterodimeros.csv", "violacoes.csv", "resultado.json"]
    resultado = json.loads(arquivos["resultado.json"])
    assert resultado["alvo"] is None and resultado["eficiencia_prevista"] is None
    app = criar_app(RAIZ)
    sessao = app.state.sessoes.obter(None)[0].sessao
    sessao.carregar_conjunto("conjunto_ok.csv", (DADOS / "conjunto_ok.csv").read_bytes())
    execucao = sessao.analisar()
    assert execucao.origem == "usuario"
    assert selecao_da_execucao(execucao, sessao.cfg_primers, {"eficiencia": {}}) is None


def test_proposta_do_desenho_vai_como_desenho_com_aviso():
    cliente = TestClient(criar_app(RAIZ))
    cliente.post("/exemplo", data={"nome": "conjunto_ok — sintético"}, follow_redirects=False)
    cliente.post("/desenhar", follow_redirects=False)
    res = cliente.post("/desenho/1/analisar", data={"posicao": "1"}, follow_redirects=False)
    assert res.headers["location"] == "/analise/2"
    selecao = json.loads(_zip(cliente, "/analise/2/exportar")["selecao.json"])
    conjunto = selecao["conjuntos"][0]
    assert (conjunto["origem"], conjunto["aviso"]) == ("desenho", AVISO)
    assert conjunto["conjunto_id"] == "proposta 1 do desenho 1"
    assert validar_selecao(selecao) == []


def test_desenho_mostra_que_as_faixas_de_tm_estao_a_calibrar():
    cliente = TestClient(criar_app(RAIZ))
    cliente.post("/exemplo", data={"nome": "conjunto_ok — sintético"}, follow_redirects=False)
    cliente.post("/desenhar", follow_redirects=False)
    res = cliente.get("/desenho/1")
    assert res.status_code == 200
    assert "filtram os candidatos do desenho e ainda estão a calibrar." in res.text
