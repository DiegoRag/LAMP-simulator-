"""Fase 3, sessão 4: faixas de Tm e regra da ponta na convenção do PrimerExplorer V5."""

from pathlib import Path

import pytest

from primers.entrada import carregar_config, ler_conjunto
from primers.restricoes import avaliar_primer_quimera, verificar, verificar_primers
from primers.termo import Condicoes, CondicoesReferencia, tm, tm_referencia

DADOS = Path(__file__).parent / "dados"
CONFIG = Path(__file__).parents[1] / "config" / "primers.toml"

# Conjunto sintético desenhado sobre alvo_ok.fasta já na convenção nova (proposta 1 do desenho).
BOM = {
    "F3": "TGACAGACTCTCCTTAGATAGCTAA",
    "B3": "AAGGCGTTGTTGCTGGATT",
    "F2": "AGACACGCATATGGATTCTGTAC",
    "F1c": "ACCCCCAAGTGTTATTGGCAAGTT",
    "B2": "ATAAAGTCTCGCCTTACTTAACCC",
    "B1c": "GGCAAGGGTAATTCTCCGTCCAC",
}


def test_faixas_usam_a_tm_de_referencia_e_nao_a_da_reacao():
    cfg = carregar_config(CONFIG)
    assert verificar_primers(BOM, cfg) == []
    ref = CondicoesReferencia.de_config(cfg)
    assert tm_referencia(BOM["F2"], ref) == pytest.approx(59.9846, abs=1e-3)
    assert tm_referencia(BOM["F1c"], ref) == pytest.approx(65.0956, abs=1e-3)
    # nas condições da reação as mesmas Tm ficariam acima das faixas
    cond = Condicoes.de_config(cfg)
    assert tm(BOM["F2"], cond, 1600.0) == pytest.approx(65.8889, abs=1e-3)
    assert tm(BOM["F1c"], cond, 1600.0) == pytest.approx(70.8526, abs=1e-3)


def test_valor_da_violacao_de_tm_e_a_tm_de_referencia():
    cfg = carregar_config(CONFIG)
    leitura = ler_conjunto(DADOS / "conjunto_ok.csv")
    v_f2 = next(v for v in verificar(leitura.conjunto, None, cfg) if v.item == "F2" and v.regra == "tm_c")
    assert v_f2.valor == pytest.approx(51.2695, abs=1e-3)
    assert v_f2.faixa == "[59.0, 61.0]"


def test_ponta_fraca_reprova_na_regra_nova():
    cfg = carregar_config(CONFIG)
    viols = verificar_primers({"F3": "TGACAGACTCTCCTTAGATAAAAAA"}, cfg)
    v_ext = [v for v in viols if v.regra == "extremidade"]
    assert len(v_ext) == 1
    assert v_ext[0].valor == pytest.approx(-2.94)
    assert v_ext[0].faixa == "<= -4.0"
    assert verificar_primers({"F3": BOM["F3"]}, cfg) == []


def test_quimera_escolhe_corte_pela_tm_de_referencia():
    cfg = carregar_config(CONFIG)
    partes, viols = avaliar_primer_quimera("FIP", BOM["F1c"] + BOM["F2"], cfg)
    assert partes == {"F1c": BOM["F1c"], "F2": BOM["F2"], "ligante": ""}
    assert viols == []
    partes_b, viols_b = avaliar_primer_quimera("BIP", BOM["B1c"] + BOM["B2"], cfg)
    assert partes_b == {"B1c": BOM["B1c"], "B2": BOM["B2"], "ligante": ""}
    assert viols_b == []


def test_restricoes_nao_dependem_das_condicoes_da_reacao():
    cfg = carregar_config(CONFIG)
    cfg_sem_reacao = {k: v for k, v in cfg.items() if k != "reacao"}
    assert verificar_primers(BOM, cfg_sem_reacao) == []
