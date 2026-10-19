import json

import numpy as np
import pytest

from lampsim.classes import (
    NOMES_CLASSES,
    diagnostico_da_curva,
    eixo_tempo,
    simular,
    sortear_amostra,
)
from lampsim.config import CLASSES, LOTE, MIX, TEMPO
from lampsim.manifesto import gerar_lote, hash_lote, montar_manifesto, salvar_manifesto


def test_cinco_classes_opticas_e_seis_grupos_no_lote():
    # fase 5: negativo e controle negativo são a mesma classe óptica; o papel é outro campo
    assert NOMES_CLASSES == (
        "positivo_forte", "positivo_perto_lod", "positivo_tardio", "negativo", "falha_reacao",
    )
    assert [g for g, _, _ in LOTE] == [
        "positivo_forte", "positivo_perto_lod", "negativo",
        "controle_negativo", "controle_contaminado", "falha_reacao",
    ]


def test_eixo_tempo():
    t = eixo_tempo()
    assert t[0] == 0 and t[-1] == TEMPO["duracao_min"]
    assert np.diff(t) == pytest.approx(TEMPO["passo_min"])


@pytest.mark.parametrize("classe", NOMES_CLASSES)
def test_parametros_dentro_das_faixas(classe):
    for semente in range(50):
        a = sortear_amostra(classe, semente)
        for p in ("ph_inicial", "queda_max", "deriva"):
            lo, hi = CLASSES[classe].get(p, MIX[p])
            assert lo <= getattr(a, p) <= hi
        lo, hi = CLASSES[classe]["eficiencia"]
        assert lo <= a.eficiencia <= hi
        assert a.copias >= 0


def test_mesma_semente_mesma_amostra():
    assert sortear_amostra("negativo", 7) == sortear_amostra("negativo", 7)


def test_classe_invalida():
    with pytest.raises(ValueError):
        sortear_amostra("talvez", 1)


def test_papel_invalido():
    with pytest.raises(ValueError):
        sortear_amostra("negativo", 1, papel="talvez")


def _f(classe, semente, papel="amostra"):
    return simular(sortear_amostra(classe, semente, papel=papel))["fracao_acida"]


@pytest.mark.parametrize("semente", range(20))
def test_comportamento_de_cada_classe(semente):
    leitura = int(TEMPO["leitura_nominal_min"] / TEMPO["passo_min"])
    # positivo forte: já amarelo na leitura nominal
    assert _f("positivo_forte", semente)[leitura] > 0.9
    # negativo, como amostra ou como controle: continua rosado até o fim
    assert _f("negativo", semente)[-1] < 0.4
    assert _f("negativo", semente, "controle_negativo")[-1] < 0.4
    # falha de reação: já começa ácido
    assert _f("falha_reacao", semente)[0] > 0.75


def test_diagnostico_sai_da_curva():
    for a in gerar_lote(20, 2026):
        s = simular(a)
        assert a.diagnostico == diagnostico_da_curva(s["t"], s["fracao_acida"], a.papel)


def test_simulacao_tem_formatos_coerentes():
    s = simular(sortear_amostra("positivo_perto_lod", 3), rota="3bandas")
    n = eixo_tempo().size
    assert s["ph"].shape == (n,) and s["rgb"].shape == (n, 3) and s["progresso"].shape == (n,)


def test_lote_reprodutivel():
    a, b = gerar_lote(5, 2026), gerar_lote(5, 2026)
    assert hash_lote(a) == hash_lote(b)
    assert hash_lote(a) != hash_lote(gerar_lote(5, 2027))
    assert len(a) == 30 and len({x.id for x in a}) == 30


def test_manifesto_salvo(tmp_path):
    lote = gerar_lote(2, 1)
    caminho = salvar_manifesto(montar_manifesto(lote, 1, 2), tmp_path)
    dados = json.loads(caminho.read_text(encoding="utf-8"))
    assert dados["hash_lote"] == hash_lote(lote)
    assert dados["n_total"] == 12 and len(dados["amostras"]) == 12
    assert dados["parametros"]["espectro"]["concentracao_uM"] == 100.0
