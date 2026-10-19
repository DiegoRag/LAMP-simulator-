"""Fase 6, sessão 7: o notebook 01_simulador roda com o modelo de classes da fase 5 (papel e diagnóstico)."""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

NOTEBOOK = Path(__file__).resolve().parents[1] / "notebooks" / "01_simulador.ipynb"
CELULA_AMBIENTE = "fb5f1b9a"  # monta o Drive no Colab e acha a pasta do código: fica fora do teste


def _celulas_de_codigo() -> list[dict]:
    caderno = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    assert caderno["nbformat"] == 4
    return [c for c in caderno["cells"] if c["cell_type"] == "code"]


def test_notebook_nao_usa_mais_o_diagnostico_da_classe():
    for celula in _celulas_de_codigo():
        fonte = "".join(celula["source"])
        assert "['diagnostico']" not in fonte and '["diagnostico"]' not in fonte
        compile(fonte, f"01_simulador#{celula['id']}", "exec")
    fontes = "".join("".join(c["source"]) for c in _celulas_de_codigo())
    assert "from lampsim.config import LOTE, QUIMICA, TEMPO, VERSAO_PARAMETROS" in fontes
    assert "a.diagnostico" in fontes and "papel=np.array([a.papel for a in lote])" in fontes


def test_notebook_roda_de_ponta_a_ponta(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(plt, "show", lambda *a, **k: None)
    ambiente = {"__name__": "__main__", "DADOS": tmp_path}
    for celula in _celulas_de_codigo():
        if celula["id"] == CELULA_AMBIENTE:
            continue
        exec(compile("".join(celula["source"]), f"01_simulador#{celula['id']}", "exec"), ambiente)
    plt.close("all")
    saida = capsys.readouterr().out
    assert "Reprodutível: mesma semente, mesmo hash." in saida
    linhas = {linha.split()[0]: linha for linha in saida.splitlines() if linha.strip()}
    assert linhas["positivo_forte"].endswith("positivo 20")
    assert linhas["negativo"].endswith("negativo 20")
    assert linhas["controle_negativo"].endswith("controle_ok 20")
    assert linhas["falha_reacao"].endswith("invalido 20")
    assert "controle_negativo" in linhas["controle_contaminado"]
    assert len(list((tmp_path / "lotes").glob("manifesto_*.json"))) == 1
