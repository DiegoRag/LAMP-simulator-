"""Fase 4, sessão 7: sessão e rota /desenhar com a faixa relaxada."""

import io
import json
import random
import tomllib
import zipfile
from pathlib import Path

from starlette.testclient import TestClient

from interface.app import criar_app
from interface.sessao import Sessao
from primers.entrada import carregar_config

RAIZ = Path(__file__).resolve().parents[1]


def _fasta_rico_em_at() -> bytes:
    rng = random.Random(7)
    seq = "".join(rng.choices("ACGT", weights=(32, 18, 18, 32), k=1000))  # GC de 33,9 %
    return f">alvo_at_sintetico\n{seq}\n".encode()


def test_sessao_desenha_com_faixa_relaxada():
    cfg = carregar_config(RAIZ / "config" / "primers.toml")
    cfg_web = tomllib.loads((RAIZ / "config" / "web.toml").read_text(encoding="utf-8"))
    sessao = Sessao(RAIZ, cfg, cfg_web)
    sessao.carregar_alvo("at.fasta", _fasta_rico_em_at())
    normal = sessao.desenhar()
    assert normal.resultado.conjuntos == []
    assert normal.resultado.estatisticas["relaxado"] is False
    rel = sessao.desenhar(relaxado=True)
    assert rel.id == normal.id + 1
    assert len(rel.resultado.conjuntos) == 3
    assert rel.resultado.estatisticas["relaxado"] is True


def test_rota_desenhar_com_relaxado():
    cliente = TestClient(criar_app(RAIZ))
    res = cliente.post(
        "/alvo", files={"arquivo": ("at.fasta", _fasta_rico_em_at(), "text/plain")}, follow_redirects=False
    )
    assert res.status_code == 303
    res = cliente.post("/desenhar", follow_redirects=False)
    assert res.headers["location"] == "/desenho/1"
    res = cliente.post("/desenhar", data={"relaxado": "1"}, follow_redirects=False)
    assert res.status_code == 303
    assert res.headers["location"] == "/desenho/2"
    for id_exec, relaxado, n in ((1, False, 0), (2, True, 3)):
        res = cliente.get(f"/desenho/{id_exec}/exportar")
        assert res.status_code == 200
        with zipfile.ZipFile(io.BytesIO(res.content)) as zf:
            dados = json.loads(zf.read("resultado.json"))
        assert dados["estatisticas"]["relaxado"] is relaxado
        assert len(dados["conjuntos"]) == n
