"""Fase 4, sessão 6: CLI desenhar com --relaxado e diagnóstico do desenho vazio, em ASCII."""

import json
import random
from pathlib import Path

from primers.cli import main

CONFIG = (Path(__file__).parents[1] / "config" / "primers.toml").resolve()


def _alvo_rico_em_at() -> str:
    rng = random.Random(7)
    return "".join(rng.choices("ACGT", weights=(32, 18, 18, 32), k=1000))  # GC de 33,9 %


def test_cli_desenhar_vazio_e_relaxado(monkeypatch, tmp_path, capsys):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "at.fasta").write_text(f">alvo_at_sintetico\n{_alvo_rico_em_at()}\n", encoding="utf-8")
    base = ["desenhar", "--alvo", "at.fasta", "--config", str(CONFIG)]

    assert main([*base, "--saida", "saida/normal"]) == 0
    saida = capsys.readouterr().out
    assert saida.isascii()
    assert "Conjuntos: 0," in saida
    assert "nucleos distintos: 0, faixas relaxadas: False." in saida
    assert (
        "Diagnostico: Ha candidatos para as seis regioes (o mais escasso e F1, com 24), mas nenhuma "
        "combinacao respeita as distancias de [distancias]. Tente o desenho com as faixas relaxadas de "
        "[desenho_relaxado]."
    ) in saida

    assert main([*base, "--saida", "saida/relaxado", "--relaxado"]) == 0
    saida = capsys.readouterr().out
    assert saida.isascii()
    assert "Conjuntos: 3," in saida
    assert "faixas relaxadas: True." in saida
    assert "Diagnostico" not in saida
    dados = json.loads((tmp_path / "saida" / "relaxado" / "resultado.json").read_text(encoding="utf-8"))
    assert dados["estatisticas"]["relaxado"] is True
    assert len(dados["conjuntos"]) == 3
    assert "busca_truncada" not in dados["estatisticas"]
