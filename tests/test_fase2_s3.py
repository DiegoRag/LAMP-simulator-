"""Testes da Fase 2 — Sessão 3: Interface Web para Entrada Aberta."""

import tomllib
from pathlib import Path

from starlette.testclient import TestClient

import primers.entrada
from interface.app import criar_app
from interface.sessao import Sessao
from interface.telas import tela_analise, tela_inicio

RAIZ = Path(__file__).resolve().parents[1]
CFG_PRIMERS = primers.entrada.carregar_config(RAIZ / "config" / "primers.toml")
CFG_WEB = tomllib.loads((RAIZ / "config" / "web.toml").read_text(encoding="utf-8"))
ESCALA_TESTE = CFG_WEB["heterodimeros"]["escala_kcal"]


def test_sessao_e_tela_inicio_csv_ponto_e_virgula_e_nomes_com_prefixo_sufixo():
    sessao = Sessao(RAIZ, CFG_PRIMERS, CFG_WEB)
    conteudo_csv = (
        b'"name";"sequence"\n'
        b'"fosA-F3";"5\'-GTGACAGACTCTCCTTAGATA-3\'"\n'
        b'"B3_v2";"TGTTGCTGGATTGTCTGTT"\n'
        b'"alvo_F2";"CGCATATGGATTCTGTAC"\n'
        b'"fosA.F1c";"GTAAACCCCCAAGTGTTATTGG"\n'
        b'"B2";"GCAATAGTATAAAGTCTCGC"\n'
        b'"B1c_1";"AAGGGTAATTCTCCGTCCA"\n'
        b'"Loop-F";"CGTATCCAGCGACAGGCA"\n'
    )
    sessao.carregar_conjunto("primers.csv", conteudo_csv)
    dados = tela_inicio(sessao)
    assert dados["conjunto"]["reconhecidos"] == [
        {"texto": "F3 ← fosA-F3", "lado": "F"},
        {"texto": "B3 ← B3_v2", "lado": "B"},
        {"texto": "F2 ← alvo_F2", "lado": "F"},
        {"texto": "F1c ← fosA.F1c", "lado": "F"},
        {"texto": "B2", "lado": "B"},
        {"texto": "B1c ← B1c_1", "lado": "B"},
        {"texto": "LF ← Loop-F", "lado": "L"},
    ]


def test_sessao_e_tela_inicio_fasta_multi_registro():
    sessao = Sessao(RAIZ, CFG_PRIMERS, CFG_WEB)
    p_conj = RAIZ / "tests" / "dados" / "conjunto_ok.csv"
    sessao.carregar_conjunto(p_conj.name, p_conj.read_bytes())

    seq_ok = (RAIZ / "tests" / "dados" / "alvo_ok.fasta").read_bytes()
    linhas_alvo = seq_ok.split(b"\n", 1)
    seq_ok_dna = linhas_alvo[1].strip()

    conteudo_multi = (
        b">alelo_vazio\n"
        + b"A" * 200
        + b"\n>alelo_correto\n"
        + seq_ok_dna
        + b"\n"
    )
    sessao.carregar_alvo("alvos.fasta", conteudo_multi)
    dados_ini = tela_inicio(sessao)
    assert dados_ini["alvo"]["resumo"] == "FASTA · 2 registros · 216 pb"

    ex = sessao.analisar()
    # Tm e ponta de conjunto_ok reprovam na convenção da fase 3; aqui só a geometria importa
    assert [v for v in ex.resultado.violacoes if v.regra not in ("tm_c", "extremidade")] == []
    assert ex.alvo.descricao == "alelo_correto"
    assert any("2 registros" in a for a in ex.resultado.avisos)


def test_tela_analise_fip_bip_sem_alvo_com_violacoes():
    sessao = Sessao(RAIZ, CFG_PRIMERS, CFG_WEB)
    conteudo_csv = (
        b"nome,sequencia\n"
        b"F3,GTGACAGACTCTCCTTAGATA\n"
        b"B3,TGTTGCTGGATTGTCTGTT\n"
        b"FIP,GTAAACCCCCAAGTGTTATTGGCGCATATGGATTCTGTAC\n"
        b"BIP,TATATAAATTCTCCGTCCAGCAATAGTATAAAGTCTCGC\n"
    )
    sessao.carregar_conjunto("primers_quimera.csv", conteudo_csv)
    ex = sessao.analisar()
    dados = tela_analise(ex, CFG_PRIMERS, tm_como_calibracao=True, escala_kcal=ESCALA_TESTE)

    assert dados["total_violacoes"] == 9
    primers_map = {p["nome"]: p for p in dados["primers"]}
    assert primers_map["BIP"]["situacao"] == "2 + Tm"
    assert primers_map["BIP"]["situacao_tipo"] == "viol"
    assert primers_map["FIP"]["situacao"] == "1 + Tm"


def test_app_nenhum_runtime_error_chega_ao_usuario(monkeypatch):
    app = criar_app(RAIZ)
    cliente = TestClient(app)

    resp_ex = cliente.post(
        "/exemplo", data={"nome": "conjunto_ok — sintético"}, follow_redirects=False
    )
    assert resp_ex.status_code == 303

    def falha_analise(*args, **kwargs):
        raise RuntimeError("falha simulada do primer3")

    monkeypatch.setattr("interface.sessao.analisar", falha_analise)

    resp_analisar = cliente.post("/analisar", follow_redirects=False)
    assert resp_analisar.status_code == 303
    assert resp_analisar.headers["location"] == "/dados"

    resp_inicio = cliente.get("/dados")
    assert resp_inicio.status_code == 200
    assert "falha simulada do primer3" in resp_inicio.text
