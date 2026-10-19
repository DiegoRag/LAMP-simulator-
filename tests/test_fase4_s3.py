"""Fase 4, sessão 3: busca por núcleo F2–B2c no alvo inteiro, monte limitado e diversidade."""

import copy
import random
from pathlib import Path

import pytest

from primers.candidatos import Candidato, gerar_candidatos
from primers.entrada import carregar_config, ler_alvo
from primers.modelo import Regiao
from primers.nucleos import (
    Nucleo,
    escolher_distintos,
    melhores_nucleos,
    mesmo_nucleo,
    percorrer_nucleos,
    sobreposicao,
)
from primers.restricoes import calcular_distancias

DADOS = Path(__file__).parent / "dados"
CONFIG = Path(__file__).parents[1] / "config" / "primers.toml"
TIPOS = ("F3", "F2", "F1", "B1c", "B2c", "B3c")


def _nucleo(f2: Regiao, b2c: Regiao) -> Nucleo:
    regs = (Regiao(0, 18), f2, Regiao(200, 220), Regiao(230, 250), b2c, Regiao(400, 420))
    combo = tuple(Candidato(t, r, "A" * (r.fim - r.ini), 60.0, 0.0) for t, r in zip(TIPOS, regs, strict=True))
    return Nucleo(0.0, tuple(r.ini for r in regs) + tuple(r.fim for r in regs), 0.0, 0.0, combo)


def _cands_1kb(cfg):
    rng = random.Random(42)
    return gerar_candidatos("".join(rng.choices("ACGT", k=1000)), cfg)


def test_sobreposicao_e_mesmo_nucleo():
    assert sobreposicao(Regiao(10, 30), Regiao(25, 40)) == 5
    assert sobreposicao(Regiao(10, 30), Regiao(30, 40)) == 0
    a = _nucleo(Regiao(100, 120), Regiao(260, 280))
    assert mesmo_nucleo(a, _nucleo(Regiao(105, 125), Regiao(262, 282)), 10) is True  # 15 e 18 nt
    assert mesmo_nucleo(a, _nucleo(Regiao(110, 130), Regiao(262, 282)), 10) is False  # F2 só 10 nt
    assert mesmo_nucleo(a, _nucleo(Regiao(100, 120), Regiao(300, 320)), 10) is False  # B2c disjunto
    assert mesmo_nucleo(a, _nucleo(Regiao(110, 130), Regiao(262, 282)), 0) is True


def test_nucleos_respeitam_as_distancias():
    cfg = carregar_config(CONFIG)
    D = cfg["distancias"]
    nucleos = list(percorrer_nucleos(_cands_1kb(cfg), cfg))
    assert len(nucleos) == 666
    for nu in nucleos:
        regioes = {t: c.regiao for t, c in zip(TIPOS, nu.combo, strict=True)}
        for nome, valor in calcular_distancias(regioes).items():
            assert D[nome][0] <= valor <= D[nome][1], (nome, valor)
        assert nu.desvio_tm == pytest.approx(sum(c.desvio_tm for c in nu.combo))
        dist = regioes["B2c"].fim - regioes["F2"].ini
        assert nu.pen_f2_b2 == pytest.approx(0.5 * max(0, 120 - dist, dist - 160))
        assert nu.pontuacao == pytest.approx(nu.desvio_tm + nu.pen_f2_b2)
        assert nu.chave == tuple(r.ini for r in regioes.values()) + tuple(r.fim for r in regioes.values())
    pares = [(nu.combo[1].regiao, nu.combo[4].regiao) for nu in nucleos]
    assert len(set(pares)) == len(pares)  # um núcleo por par (F2, B2c), na ordem de F2
    assert [p[0].ini for p in pares] == sorted(p[0].ini for p in pares)


def test_cada_nucleo_tem_o_melhor_conjunto():
    cfg = carregar_config(CONFIG)
    D = cfg["distancias"]
    cands = gerar_candidatos(ler_alvo(DADOS / "alvo_ok.fasta"), cfg)
    nucleos = list(percorrer_nucleos(cands, cfg))
    assert len(nucleos) == 6
    for nu in nucleos:
        f2, b2c = nu.combo[1], nu.combo[4]
        melhor = None
        for f3 in cands["F3"]:
            for f1 in cands["F1"]:
                for b1c in cands["B1c"]:
                    for b3c in cands["B3c"]:
                        reg = {"F3": f3.regiao, "F2": f2.regiao, "F1": f1.regiao, "B1c": b1c.regiao,
                               "B2c": b2c.regiao, "B3c": b3c.regiao}
                        if all(D[k][0] <= v <= D[k][1] for k, v in calcular_distancias(reg).items()):
                            soma = f3.desvio_tm + f1.desvio_tm + b1c.desvio_tm + b3c.desvio_tm
                            melhor = soma if melhor is None else min(melhor, soma)
        assert nu.desvio_tm == pytest.approx(melhor + f2.desvio_tm + b2c.desvio_tm)


def test_monte_limitado_e_distintos():
    cfg = carregar_config(CONFIG)
    cands = _cands_1kb(cfg)
    todos = sorted(percorrer_nucleos(cands, cfg), key=lambda x: (x.pontuacao, x.chave))
    guardados, total = melhores_nucleos(cands, cfg)
    assert total == 666
    assert guardados == todos  # 666 < nucleos_guardados = 5000
    cfg_3 = copy.deepcopy(cfg)
    cfg_3["desenho"]["nucleos_guardados"] = 3
    assert melhores_nucleos(cands, cfg_3) == (todos[:3], 666)
    assert round(todos[0].pontuacao, 4) == 0.067

    distintos = escolher_distintos(guardados, cfg, lambda nu: True)
    assert len(distintos) == 50  # pre_selecao
    assert distintos[0] == guardados[0]
    assert all(not mesmo_nucleo(a, b, 10) for i, a in enumerate(distintos) for b in distintos[i + 1 :])
    sem_o_primeiro = escolher_distintos(guardados, cfg, lambda nu: nu is not guardados[0])
    assert guardados[0] not in sem_o_primeiro  # aceita() é respeitado
    assert len(sem_o_primeiro) == 50
