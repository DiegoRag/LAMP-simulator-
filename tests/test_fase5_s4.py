"""Fase 5, sessão 4: incompatibilidades por região medidas pela análise (primers.incompat)."""

from pathlib import Path

import pytest

from primers.analise import analisar
from primers.entrada import carregar_config, ler_alvos, ler_conjunto
from primers.incompat import distancia_ponta_extensao, medir_incompat
from primers.modelo import Regiao

RAIZ = Path(__file__).resolve().parents[1]
DADOS = RAIZ / "tests" / "dados"
CFG = carregar_config(RAIZ / "config" / "primers.toml")
# alvo sintético: F3 0–21, F2 31–49, F1 81–103, B1c 108–127, B2c 157–177, B3c 187–206 (base 0)
ALVO = ler_alvos(DADOS / "alvo_ok.fasta")[0].sequencia
TROCA = {"A": "C", "C": "A", "G": "T", "T": "G"}


def _mutar(posicoes: list[int]) -> str:
    bases = list(ALVO)
    for p in posicoes:
        bases[p] = TROCA[bases[p]]
    return "".join(bases)


def _analisar(alvo: str | None):
    leitura = ler_conjunto(DADOS / "conjunto_ok.csv", alvo=alvo, ligante="", max_erros=3, circular=False)
    return analisar(leitura, alvo, CFG)


def test_config_tem_bases_da_ponta():
    assert CFG["casamento"]["bases_ponta_3p"] == 5


def test_distancia_ate_a_ponta_que_estende():
    # F3, F2, LB: ponta 3' do primer no fim da região; B2c, B3c, LF: no início
    assert distancia_ponta_extensao("F3", Regiao(0, 21), 20) == 0
    assert distancia_ponta_extensao("F3", Regiao(0, 21), 15) == 5
    assert distancia_ponta_extensao("B2c", Regiao(157, 177), 157) == 0
    assert distancia_ponta_extensao("B3c", Regiao(187, 206), 190) == 3
    # F1 e B1c: a fita do haltere estende a partir da cópia da ponta 5' de F1c (B1c)
    assert distancia_ponta_extensao("F1", Regiao(81, 103), 102) == 0
    assert distancia_ponta_extensao("F1", Regiao(81, 103), 81) == 21
    assert distancia_ponta_extensao("B1c", Regiao(108, 127), 108) == 0
    assert distancia_ponta_extensao("B1c", Regiao(108, 127), 126) == 18


def test_conjunto_perfeito_tem_zero_incompatibilidades():
    inc = _analisar(ALVO).incompat
    assert (inc.total, inc.ponta3, inc.nao_encontradas) == (0, 0, ())
    assert [i.regiao for i in inc.por_regiao] == ["F3", "F2", "F1", "B1c", "B2c", "B3c"]


def test_seis_erros_na_ponta_3():
    # o exemplo da revisão: uma troca na base terminal de cada uma das seis regiões
    inc = _analisar(_mutar([20, 48, 102, 108, 157, 187])).incompat
    assert (inc.total, inc.ponta3) == (6, 6)
    assert all((i.total, i.ponta3) == (1, 1) for i in inc.por_regiao)


def test_erros_longe_da_ponta():
    # F3 a 15 nt da ponta; F1 e B1c na ponta 3' do primer, mas longe da ponta que estende
    inc = _analisar(_mutar([5, 81, 126])).incompat
    assert (inc.total, inc.ponta3) == (3, 0)
    assert [(i.regiao, i.total) for i in inc.por_regiao if i.total] == [("F3", 1), ("F1", 1), ("B1c", 1)]


def test_limite_de_cinco_bases():
    # F3: posição 16 está a 4 nt da ponta (conta), posição 15 a 5 nt (não conta)
    inc = _analisar(_mutar([15, 16])).incompat
    f3 = inc.por_regiao[0]
    assert (f3.regiao, f3.total, f3.ponta3) == ("F3", 2, 1)


def test_regiao_nao_encontrada():
    inc = _analisar(_mutar([31, 33, 35, 37, 39])).incompat  # 5 trocas em F2 > max_erros
    assert inc.nao_encontradas == ("F2",)
    assert inc.total is None and inc.ponta3 is None
    f2 = inc.por_regiao[1]
    assert (f2.regiao, f2.encontrada, f2.total, f2.ponta3) == ("F2", False, None, None)


def test_sem_alvo_nao_ha_medida():
    assert _analisar(None).incompat is None


def test_bases_ponta_invalida():
    leitura = ler_conjunto(DADOS / "conjunto_ok.csv", alvo=ALVO, ligante="", max_erros=3, circular=False)
    from primers.restricoes import localizar_detalhado

    loc = localizar_detalhado(leitura.conjunto, ALVO, CFG)
    with pytest.raises(ValueError):
        medir_incompat(leitura.conjunto, loc, 0)
