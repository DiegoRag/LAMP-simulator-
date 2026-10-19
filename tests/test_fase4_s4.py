"""Fase 4, sessão 4: desenho no alvo inteiro, propostas distintas, laços nos dímeros e FIP longo."""

import random
from pathlib import Path

import pytest

from primers.candidatos import reprova_homodimero
from primers.desenho import desenhar, pior_dimero_tubo
from primers.entrada import carregar_config
from primers.termo import Condicoes, dimero_ponta3

CONFIG = Path(__file__).parents[1] / "config" / "primers.toml"


def test_alvo_de_10_kb():
    cfg = carregar_config(CONFIG)
    cond = Condicoes.de_config(cfg)
    rng = random.Random(42)
    res = desenhar("".join(rng.choices("ACGT", k=10_000)), cfg)
    # antes da fase 4: 5 propostas, 4 delas com o mesmo F2 (8059), e só 15 combinações examinadas
    assert [c.regioes["F2"].ini for c in res.conjuntos] == [9577, 1940, 8379, 4523, 7222]
    assert [round(c.pontuacao, 4) for c in res.conjuntos] == [0.1115, 0.1845, 0.2004, 0.2102, 0.3073]
    e = res.estatisticas
    assert (e["conjuntos_examinados"], e["nucleos_guardados"], e["nucleos_distintos"]) == (10025, 5000, 50)
    assert (e["relaxado"], e["diagnostico"]) == (False, "")
    assert "busca_truncada" not in e
    assert sum(e["eliminados_por_tipo"]["F2"].values()) + e["candidatos_por_tipo"]["F2"] == 79836

    def sobre(r, s) -> int:
        return max(0, min(r.fim, s.fim) - max(r.ini, s.ini))

    for i, a in enumerate(res.conjuntos):  # núcleos distintos
        for b in res.conjuntos[i + 1 :]:
            f2, b2 = sobre(a.regioes["F2"], b.regioes["F2"]), sobre(a.regioes["B2c"], b.regioes["B2c"])
            assert not (f2 > 10 and b2 > 10)
    for c in res.conjuntos:  # homodímero testado em todo primer proposto
        for tipo, primer in (("F3", c.conjunto.F3), ("F2", c.conjunto.F2), ("F1", c.conjunto.F1c),
                             ("B1c", c.conjunto.B1c), ("B2c", c.conjunto.B2), ("B3c", c.conjunto.B3)):
            assert reprova_homodimero(primer, tipo, cfg, cond) is False
        dist = c.regioes["B2c"].fim - c.regioes["F2"].ini
        pen = 0.5 * max(0, 120 - dist, dist - 160)
        assert c.pontuacao == pytest.approx(c.desvio_tm + pen + 0.5 * max(0.0, -c.pior_heterodimero_kcal))

    c1 = res.conjuntos[0]  # dímeros com os laços
    assert c1.avisos == []
    assert set(c1.regioes) == {"F3", "F2", "F1", "B1c", "B2c", "B3c", "LF", "LB"}
    tubo = [c1.conjunto.F3, c1.conjunto.B3, c1.conjunto.FIP, c1.conjunto.BIP, c1.conjunto.LF, c1.conjunto.LB]
    pares = [dimero_ponta3(tubo[i], tubo[j], cond).dg_kcal for i in range(6) for j in range(i + 1, 6)]
    assert c1.pior_heterodimero_kcal == pytest.approx(min(pares))


def test_fip_acima_do_limite_do_primer3():
    cond = Condicoes.de_config(carregar_config(CONFIG))
    fip = "ACGTTGCAAGTCCTGATGCAGTCAGT" * 2 + "GATTCAGCACTGCATG"  # 68 nt
    b3 = "AAGGCGTTGTTGCTGGATT"
    pior, avisos = pior_dimero_tubo({"FIP": fip, "B3": b3}, cond, 60)
    assert avisos == ["FIP com 68 nt: dímeros avaliados nos 60 nt da ponta 3′ (limite do primer3)"]
    assert pior == pytest.approx(dimero_ponta3(fip[-60:], b3, cond).dg_kcal)
    assert pior_dimero_tubo({"B3": b3}, cond, 60) == (0.0, [])
