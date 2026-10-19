"""Fase 4, sessão 9: desenho de ponta a ponta em alvos sintéticos longos, ricos em AT e sem solução."""

import random
from pathlib import Path

from primers.desenho import AVISO_RELAXADO, config_relaxada, desenhar
from primers.entrada import carregar_config
from primers.restricoes import verificar

CONFIG = Path(__file__).parents[1] / "config" / "primers.toml"


def _aleatorio(n: int, semente: int, pesos: tuple[int, int, int, int] = (1, 1, 1, 1)) -> str:
    rng = random.Random(semente)
    return "".join(rng.choices("ACGT", weights=pesos, k=n))


def test_alvo_de_30_kb_tem_5_propostas():
    # antes da fase 4: 0 propostas em 162 s (o teto global de 200 caía em lugares sem relação)
    res = desenhar(_aleatorio(30_000, 42), carregar_config(CONFIG))
    assert [c.regioes["F2"].ini for c in res.conjuntos] == [19462, 9577, 8379, 21894, 16801]
    assert res.estatisticas["conjuntos_examinados"] == 30043


def test_faixa_relaxada_para_alvo_rico_em_at():
    cfg = carregar_config(CONFIG)
    alvo = _aleatorio(1000, 7, (32, 18, 18, 32))  # GC de 33,9 %
    res = desenhar(alvo, cfg)
    assert res.conjuntos == []
    assert res.estatisticas["candidatos_por_tipo"]["F1"] == 24
    assert res.estatisticas["diagnostico"] == (
        "Há candidatos para as seis regiões (o mais escasso é F1, com 24), mas nenhuma combinação "
        "respeita as distâncias de [distancias]. Tente o desenho com as faixas relaxadas de "
        "[desenho_relaxado]."
    )
    rel = desenhar(alvo, cfg, relaxado=True)
    assert len(rel.conjuntos) == 3
    assert (rel.estatisticas["relaxado"], rel.estatisticas["diagnostico"]) == (True, "")
    for c in rel.conjuntos:
        assert c.avisos[-1] == AVISO_RELAXADO
        assert verificar(c.conjunto, alvo, config_relaxada(cfg)) == []
    assert any(verificar(c.conjunto, alvo, cfg) for c in rel.conjuntos)  # fora do guia


def test_diagnostico_de_ponta_a_ponta():
    cfg = carregar_config(CONFIG)
    res = desenhar("A" * 300, cfg)
    assert res.conjuntos == []
    assert res.estatisticas["diagnostico"] == (
        "Nenhum candidato a F3, F2, F1, B1c, B2c e B3c: a regra que mais eliminou janelas foi GC fora "
        "da faixa de [gc_pct]. Tente o desenho com as faixas relaxadas de [desenho_relaxado]."
    )
    assert desenhar("ACGT", cfg).estatisticas["diagnostico"].startswith(
        "Nenhuma janela para F3, F2, F1, B1c, B2c e B3c: o alvo é curto demais."
    )
    assert desenhar("N" * 300, cfg).estatisticas["eliminados_por_tipo"]["F2"]["ambigua"] == 2236
