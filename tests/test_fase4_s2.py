"""Fase 4, sessão 2: candidatos com teto por janela, unicidade nas duas fitas e filtros de estrutura."""

import copy
import random
from collections import Counter
from pathlib import Path

from primers.candidatos import (
    INFO_TIPO,
    REGRAS,
    TIPOS_BASE,
    gerar_candidatos,
    reprova_grampo,
    reprova_homodimero,
)
from primers.entrada import carregar_config, ler_alvo
from primers.modelo import Regiao, rc
from primers.restricoes import verificar_primers
from primers.termo import Condicoes, grampo, homodimero

DADOS = Path(__file__).parent / "dados"
CONFIG = Path(__file__).parents[1] / "config" / "primers.toml"


def test_teto_por_janela_no_alvo_inteiro():
    cfg = carregar_config(CONFIG)
    cond = Condicoes.de_config(cfg)
    rng = random.Random(42)
    eliminados: dict = {}
    cands = gerar_candidatos("".join(rng.choices("ACGT", k=1000)), cfg, eliminados=eliminados)
    assert {t: len(v) for t, v in cands.items()} == {
        "F3": 157, "F2": 157, "F1": 136, "B1c": 128, "B2c": 157, "B3c": 157,
    }
    assert eliminados["F3"] == {
        "ambigua": 0, "gc_pct": 2427, "tm_c": 4612, "extremidade": 129,
        "repetida": 0, "grampo": 46, "homodimero": 0, "teto": 465,
    }
    for tipo in TIPOS_BASE:
        por_janela = Counter(c.regiao.ini // 10 for c in cands[tipo])
        assert max(por_janela.values()) <= 2
        assert min(por_janela) < 10 and max(por_janela) > 90  # começo e fim do alvo, sem truncar
        assert sum(eliminados[tipo].values()) + len(cands[tipo]) == 7836  # cada janela de 18–25 nt, 1 vez
        assert list(eliminados[tipo]) == list(REGRAS)
        nome, _, _, chave_conc, _ = INFO_TIPO[tipo]
        conc = cfg["concentracao_nM"][chave_conc]
        for c in cands[tipo]:
            assert verificar_primers({nome: c.primer}, cfg) == []
            g = grampo(c.primer, cond, conc)
            h = homodimero(c.primer, cond, conc)
            assert not (g.encontrada and g.tm_c >= 65.0)
            assert not (h.encontrada and h.tm_c >= 65.0)


def test_eliminados_no_alvo_ok():
    cfg = carregar_config(CONFIG)
    eliminados: dict = {}
    cands = gerar_candidatos(ler_alvo(DADOS / "alvo_ok.fasta"), cfg, eliminados=eliminados)
    assert len(cands["F2"]) == 35
    assert eliminados["F2"] == {
        "ambigua": 0, "gc_pct": 540, "tm_c": 824, "extremidade": 43,
        "repetida": 0, "grampo": 0, "homodimero": 0, "teto": 122,
    }
    assert eliminados["F1"]["teto"] == 64


def test_reprova_grampo_e_homodimero():
    cfg = carregar_config(CONFIG)
    cond = Condicoes.de_config(cfg)
    assert reprova_grampo("GACCGCGGTTTTCCGCGGTCATAGT", "F2", cfg, cond) is True  # grampo de Tm 93 °C
    assert reprova_homodimero("GACCGCGGTTTTCCGCGGTCATAGT", "F2", cfg, cond) is False  # 60,9 °C
    assert reprova_homodimero("CTGACCGCGGAAACCGCGGTCATGA", "F2", cfg, cond) is True  # 65,3 °C a 1,6 µM
    assert reprova_grampo("AGACACGCATATGGATTCTGTAC", "F2", cfg, cond) is False
    assert reprova_homodimero("AGACACGCATATGGATTCTGTAC", "F2", cfg, cond) is False
    cfg_margem = copy.deepcopy(cfg)
    cfg_margem["desenho"]["homodimero_margem_C"] = -1.0  # limite passa a 66 °C
    assert reprova_homodimero("CTGACCGCGGAAACCGCGGTCATGA", "F2", cfg_margem, cond) is False


def test_repeticao_invertida_base_ambigua_e_lacos():
    cfg = carregar_config(CONFIG)
    alvo = ler_alvo(DADOS / "alvo_ok.fasta")
    assert Regiao(26, 49) in [c.regiao for c in gerar_candidatos(alvo, cfg, tipos=("F2",))["F2"]]
    elim: dict = {}
    cands = gerar_candidatos(alvo + rc(alvo[26:49]), cfg, tipos=("F2",), eliminados=elim)
    assert Regiao(26, 49) not in [c.regiao for c in cands["F2"]]  # a cópia só existe na fita oposta
    assert elim["F2"]["repetida"] == 2

    elim = {}
    cands = gerar_candidatos(alvo[:30] + "N" + alvo[31:], cfg, eliminados=elim)
    assert all(not c.regiao.ini <= 30 < c.regiao.fim for lista in cands.values() for c in lista)
    assert {t: elim[t]["ambigua"] for t in TIPOS_BASE} == dict.fromkeys(TIPOS_BASE, 172)

    elim = {}
    lf = gerar_candidatos(alvo, cfg, tipos=("LF",), dentro=Regiao(40, 110), eliminados=elim)["LF"]
    assert len(lf) == 29
    assert elim["LF"]["teto"] == 0  # laços (com `dentro`) não têm teto
    assert all(40 <= c.regiao.ini and c.regiao.fim <= 110 for c in lf)
    assert lf[0].regiao == Regiao(60, 83)
    assert lf[0].primer == rc(alvo[60:83])
