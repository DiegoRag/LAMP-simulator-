"""Fase 6, sessão 4: telas coerentes com "a calibrar", réguas com negativos, base 1 e erros na ponta."""

import copy
from pathlib import Path

import pytest

import primers.casamento
from interface.sessao import Sessao
from interface.telas import (
    regua_distancia,
    resumo_violacoes,
    tela_analise,
    tela_desenho,
    tela_inicio,
    tela_parametros,
    texto_erro_de_base,
    texto_faixas_tm,
    texto_fita,
)
from primers.casamento import ErroDeBase
from primers.entrada import carregar_config
from primers.modelo import rc
from primers.restricoes import Violacao

RAIZ = Path(__file__).resolve().parents[1]
DADOS = RAIZ / "tests" / "dados"
CFG = carregar_config(RAIZ / "config" / "primers.toml")
CFG_WEB = {
    "entrada": {"tamanho_max_kb": 64, "alvo_desenho_max_pb": 2000},
    "analise": {"tm_como_calibracao": True},
    "exemplos": [
        {"nome": "sintético", "conjunto": "tests/dados/conjunto_ok.csv", "alvo": "tests/dados/alvo_ok.fasta"}
    ],
}
ESCALA = [-10.0, -8.0, -6.0, -4.0]
ALVO = (DADOS / "alvo_ok.fasta").read_text(encoding="utf-8").splitlines()[1]


def _v(regra: str) -> Violacao:
    return Violacao(item="F3", regra=regra, valor=1.0, faixa="", mensagem="")


def _sessao(cfg_web: dict = CFG_WEB) -> Sessao:
    return Sessao(RAIZ, CFG, cfg_web)


def test_resumo_de_violacoes_com_e_sem_a_calibrar():
    tm, ponta = _v("tm_c"), _v("extremidade")
    assert resumo_violacoes([tm, tm, ponta], True) == {
        "n_violacoes": 1,
        "n_tm": 2,
        "detalhe": "e 2 de Tm a calibrar",
        "tipo": "viol",
        "texto": "1 violação e 2 de Tm a calibrar",
    }
    # só Tm fora da faixa: o KPI não fica vermelho, fica "a calibrar"
    assert resumo_violacoes([tm], True) == {
        "n_violacoes": 0,
        "n_tm": 1,
        "detalhe": "e 1 de Tm a calibrar",
        "tipo": "calibrar",
        "texto": "0 violações e 1 de Tm a calibrar",
    }
    assert resumo_violacoes([], True)["tipo"] == "ok"
    assert resumo_violacoes([tm, ponta], False) == {
        "n_violacoes": 2,
        "n_tm": 1,
        "detalhe": "1 de Tm",
        "tipo": "viol",
        "texto": "2 violações",
    }


@pytest.mark.parametrize(
    ("d", "d_min", "d_max", "esquerda", "largura", "marca", "zero"),
    [
        (-12, 0, 60, "17.085", "62.696", "4.545", "17.085"),  # régua de −16,35 a 79,35
        (30, 0, 1000, "4.545", "72.727", "6.727", ""),  # faixa larga: de −62,5 a 1.312,5
        (400, 120, 180, "31.818", "13.636", "95.455", ""),  # bem acima da faixa: de −20 a 420
        (5, 5, 5, "77.273", "0.000", "77.273", ""),
    ],
)
def test_regua_cabe_a_marca(d, d_min, d_max, esquerda, largura, marca, zero):
    assert regua_distancia(d, d_min, d_max) == {
        "banda_esquerda": esquerda,
        "banda_largura": largura,
        "marca": marca,
        "zero": zero,
    }


def test_regua_nunca_sai_dos_limites():
    for d in range(-300, 1500, 7):
        for d_min, d_max in ((0, 60), (40, 60), (120, 180), (0, 1000), (-20, 10)):
            r = regua_distancia(d, d_min, d_max)
            assert 4.5 <= float(r["marca"]) <= 95.5
            assert float(r["banda_esquerda"]) >= 0.0
            assert float(r["banda_esquerda"]) + float(r["banda_largura"]) <= 100.001


def test_distancia_negativa_na_analise():
    # F3 deslocado 20 nt para dentro do alvo: sobrepõe F2 em 10 nt (F3 → F2 = −10)
    conjunto = (DADOS / "conjunto_ok.csv").read_text(encoding="utf-8").replace(
        "GTGACAGACTCTCCTTAGATA", ALVO[20:41]
    )
    sessao = _sessao()
    sessao.carregar_conjunto("sobreposto.csv", conjunto.encode("utf-8"))
    sessao.carregar_alvo("alvo_ok.fasta", (DADOS / "alvo_ok.fasta").read_bytes())
    dados = tela_analise(sessao.analisar(), CFG, True, ESCALA)
    assert dados["distancias"][0] == {
        "nome": "F3 → F2",
        "medida": "−10",
        "faixa": "0–60",
        "dentro": False,
        "texto": "10 abaixo do mínimo",
        "regua": {"banda_esquerda": "15.241", "banda_largura": "64.171", "marca": "4.545", "zero": "15.241"},
    }
    assert dados["mapa"]["faixas_fora"][0]["texto"] == "F3 → F2: −10 nt · mínimo 0"
    assert "Distância −10 nt" in [i["texto"] for i in dados["grupos_violacoes"][0]["itens"]]


def test_texto_do_erro_de_base():
    na_ponta_5 = ErroDeBase("F1", "F1c", "+", 103, 1, 22, "G", "A", "5'", 0)
    assert texto_erro_de_base(na_ponta_5, 5) == (
        "base 1 de 22 do F1c: G no primer, o alvo pede A, a 0 nt da ponta 5′ (a que fecha o laço), "
        "dentro das 5 bases da ponta; posição 103"
    )
    longe = ErroDeBase("B2c", "B2", "-", 1234, 3, 20, "C", "T", "3'", 17)
    assert texto_erro_de_base(longe, 5) == (
        "base 3 de 20 do B2: C no primer, o alvo pede T, a 17 nt da ponta 3′; posição 1.234"
    )
    assert (texto_fita("+"), texto_fita("-"), texto_fita(None)) == ("+", "−", "+")


def test_erro_na_ponta_5_de_f1c_na_tela():
    sessao = _sessao()
    sessao.carregar_conjunto("conjunto_ok.csv", (DADOS / "conjunto_ok.csv").read_bytes())
    alvo_mut = ALVO[:102] + "T" + ALVO[103:]  # última base de F1 = ponta 5′ de F1c
    sessao.carregar_alvo("alvo_mut.fasta", f">mut\n{alvo_mut}\n".encode())
    dados = tela_analise(sessao.analisar(), CFG, True, ESCALA)
    assert dados["kpis"][2]["detalhe"] == "F1c com 1 dif., 1 na ponta"
    itens = dados["grupos_violacoes"][0]["itens"]
    item = next(i for i in itens if i["item"] == "F1c" and i["texto"].startswith("Não"))
    assert item["texto"] == (
        "Não localizado: melhor posição 82–103 (fita +), 1 diferença (base 1 de 22 do F1c: G no primer, "
        "o alvo pede A, a 0 nt da ponta 5′ (a que fecha o laço), dentro das 5 bases da ponta; posição 103)"
    )


def test_a_pagina_nao_repete_a_busca(monkeypatch):
    sessao = _sessao()
    sessao.carregar_conjunto("conjunto_ok.csv", (DADOS / "conjunto_ok.csv").read_bytes())
    sessao.carregar_alvo("alvo_mut.fasta", f">mut\n{ALVO[:48]}A{ALVO[49:]}\n".encode())
    execucao = sessao.analisar()

    def proibida(*args, **kwargs):
        raise AssertionError("a tela não deve buscar de novo no alvo")

    monkeypatch.setattr(primers.casamento, "melhor_posicao", proibida)
    monkeypatch.setattr(primers.casamento, "_buscar_em_fita", proibida)
    dados = tela_analise(execucao, CFG, True, ESCALA)
    assert dados["kpis"][2]["detalhe"] == "F2 com 1 dif., 1 na ponta"


def test_mapa_em_base_1_e_com_a_fita():
    sessao = _sessao()
    sessao.carregar_conjunto("conjunto_ok.csv", (DADOS / "conjunto_ok.csv").read_bytes())
    sessao.carregar_alvo("alvo_rc.fasta", f">rc\n{rc(ALVO)}\n".encode())
    execucao = sessao.analisar()
    dados = tela_analise(execucao, CFG, True, ESCALA)
    assert dados["mapa"]["fita"] == "−"
    assert dados["mapa"]["legenda"] == (
        "fita − (complemento reverso do alvo) · 216 pb · coordenadas a partir de 1 · "
        "F1, B2c, B3c e LF são lidas pelo complemento reverso do primer"
    )
    assert dados["mapa"]["ticks"][0]["valor"] == "1"
    assert dados["mapa"]["descricao"].startswith(
        "Mapa linear do alvo de 216 pb, fita −, coordenadas a partir de 1, com as regiões F3 em 1–21, "
    )
    assert dados["url_tubos"] == f"/tubos?execucao={execucao.id}"


def test_historico_segue_a_chave_do_config():
    for chave, texto, tipo in (
        (True, "2 violações e 6 de Tm a calibrar · 6 de 6 regiões", "viol"),
        (False, "8 violações · 6 de 6 regiões", "viol"),
    ):
        cfg_web = copy.deepcopy(CFG_WEB)
        cfg_web["analise"]["tm_como_calibracao"] = chave
        sessao = _sessao(cfg_web)
        sessao.carregar_exemplo("sintético")
        sessao.analisar()
        hist = tela_inicio(sessao)["historico"][0]
        assert (hist["resultado"], hist["resultado_tipo"]) == (texto, tipo)


def test_sem_alvo_nao_ha_link_dos_tubos():
    sessao = _sessao()
    sessao.carregar_conjunto("conjunto_ok.csv", (DADOS / "conjunto_ok.csv").read_bytes())
    assert tela_analise(sessao.analisar(), CFG, True, ESCALA)["url_tubos"] == ""


def test_parametros_rotulam_tm_e_dntp():
    assert texto_faixas_tm(CFG) == "F3/B3 59–61, F2/B2 59–61, F1c/B1c 64–66 e loops 64–66 °C"
    sessao = _sessao()
    sessao.carregar_exemplo("sintético")
    sessao.analisar()
    dados = tela_parametros(sessao, CFG)
    assert dados["tm_a_calibrar"] is True
    assert dados["regras"][-1] == {
        "nome": "Tm de referência",
        "valor": "F3/B3 59–61, F2/B2 59–61, F1c/B1c 64–66 e loops 64–66 °C",
        "marca": "a calibrar",
    }
    assert {linha["marca"] for linha in dados["calibracao"]["linhas"]} == {"a calibrar"}

    cfg_web = copy.deepcopy(CFG_WEB)
    cfg_web["analise"]["tm_como_calibracao"] = False
    sem = tela_parametros(_sessao(cfg_web), CFG)
    assert sem["tm_a_calibrar"] is False and sem["regras"][-1]["marca"] == "conferir"

    # dNTP total (4 × 2,5 = 10 mM) acima do Mg²⁺ total (8 mM): a tela avisa em vez de mostrar 0 calado
    cfg_total = copy.deepcopy(CFG)
    cfg_total["reacao"]["dntp_cada_mM"] = 2.5
    mg = tela_parametros(_sessao(), cfg_total)["reacao"][1]
    assert mg == {
        "nome": "Mg²⁺",
        "valor": "8 mM",
        "nota": (
            "Livre: 0 mM: o dNTP total (10 mM) não é menor que o Mg²⁺ total (8 mM). "
            "Confira se reacao.dntp_cada_mM está em mM de cada dNTP, e não do total."
        ),
        "destaque": True,
    }


def test_desenho_avisa_que_as_faixas_de_tm_estao_a_calibrar():
    sessao = _sessao()
    sessao.carregar_exemplo("sintético")
    execucao = sessao.desenhar()
    assert tela_desenho(execucao, CFG, 1, False)["nota_tm"] == ""
    assert tela_desenho(execucao, CFG, 1, False, tm_a_calibrar=True)["nota_tm"] == (
        "As faixas de Tm de referência (F3/B3 59–61, F2/B2 59–61, F1c/B1c 64–66 e loops 64–66 °C) "
        "filtram os candidatos do desenho e ainda estão a calibrar."
    )
