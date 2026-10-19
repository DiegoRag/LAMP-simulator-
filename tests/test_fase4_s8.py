"""Fase 4, sessão 8: tela do desenho com diagnóstico, faixa relaxada e o resumo da busca."""

import random
import tomllib
from pathlib import Path

from starlette.testclient import TestClient

from interface.app import criar_app
from interface.sessao import Sessao
from interface.telas import tela_desenho, tela_parametros
from primers.entrada import carregar_config

RAIZ = Path(__file__).resolve().parents[1]
cfg = carregar_config(RAIZ / "config" / "primers.toml")
cfg_web = tomllib.loads((RAIZ / "config" / "web.toml").read_text(encoding="utf-8"))


def _fasta_rico_em_at() -> bytes:
    rng = random.Random(7)
    seq = "".join(rng.choices("ACGT", weights=(32, 18, 18, 32), k=1000))  # GC de 33,9 %
    return f">alvo_at_sintetico\n{seq}\n".encode()


def test_tela_do_desenho_vazio_oferece_a_faixa_relaxada():
    sessao = Sessao(RAIZ, cfg, cfg_web)
    sessao.carregar_alvo("at.fasta", _fasta_rico_em_at())
    vazio = tela_desenho(sessao.desenhar(), cfg, sel=1, agrupar=False)
    assert vazio["vazio"] is True
    assert vazio["relaxado"] is False
    assert vazio["pode_relaxar"] is True
    assert vazio["mensagem_vazio"].startswith("Há candidatos para as seis regiões (o mais escasso é F1")
    assert "truncada" not in vazio
    rel = tela_desenho(sessao.desenhar(relaxado=True), cfg, sel=1, agrupar=False)
    assert rel["vazio"] is False
    assert rel["relaxado"] is True
    assert rel["pode_relaxar"] is False
    assert rel["mensagem_vazio"] == ""
    assert len(rel["lista"]) == 3


def test_parametros_do_desenho():
    desenho = tela_parametros(Sessao(RAIZ, cfg, cfg_web), cfg)["desenho"]
    assert [d["nome"] for d in desenho] == [
        "Teto por janela",
        "Núcleos guardados",
        "Pré-seleção",
        "Peso do dímero",
        "Resultados",
        "Ligante FIP/BIP",
        "Propostas distintas",
    ]
    assert desenho[0]["valor"] == "2 por tipo a cada 10 nt"
    assert desenho[1]["valor"] == "até 5.000, do alvo inteiro"
    assert desenho[6]["valor"] == "F2 e B2c não sobrepõem mais de 10 nt ao mesmo tempo"


def test_pagina_do_desenho_vazio_e_relaxado():
    cliente = TestClient(criar_app(RAIZ))
    cliente.post("/alvo", files={"arquivo": ("at.fasta", _fasta_rico_em_at(), "text/plain")})
    cliente.post("/desenhar", follow_redirects=False)
    res = cliente.get("/desenho/1")
    assert res.status_code == 200
    assert "o mais escasso é F1, com 24" in res.text
    assert 'name="relaxado" value="1"' in res.text
    assert "Desenhar com as faixas relaxadas" in res.text
    assert "núcleos F2–B2c examinados" in res.text
    assert "Alvo inteiro percorrido" in res.text
    cliente.post("/desenhar", data={"relaxado": "1"}, follow_redirects=False)
    res = cliente.get("/desenho/2")
    assert res.status_code == 200
    assert "Faixas relaxadas, fora do guia" in res.text
    assert 'name="relaxado"' not in res.text
    assert "Cada proposta vem de um núcleo F2–B2c diferente." in res.text
