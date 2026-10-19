"""Fase 3, sessão 6: desenho com a Tm de referência e ponta 3′; saída com as duas Tm."""

import csv
import json
from pathlib import Path

import pytest

from primers.analise import analisar
from primers.desenho import desenhar, gerar_candidatos
from primers.entrada import carregar_config, ler_alvo, ler_conjunto
from primers.modelo import Regiao
from primers.restricoes import verificar
from primers.saida import exportar
from primers.termo import Condicoes, CondicoesReferencia, dimero_ponta3, tm_referencia

DADOS = Path(__file__).parent / "dados"
CONFIG = Path(__file__).parents[1] / "config" / "primers.toml"


def test_candidatos_usam_a_tm_de_referencia():
    cfg = carregar_config(CONFIG)
    ref = CondicoesReferencia.de_config(cfg)
    cands = gerar_candidatos(ler_alvo(DADOS / "alvo_ok.fasta"), cfg, tipos=("F2", "F1"))
    assert len(cands["F2"]) == 35  # fase 4: teto por janela
    for c in cands["F2"]:
        assert c.tm_c == pytest.approx(tm_referencia(c.primer, ref))
        assert 59.0 <= c.tm_c <= 61.0
        assert c.desvio_tm == pytest.approx(abs(c.tm_c - 60.0))
    for c in cands["F1"]:
        assert 64.0 <= c.tm_c <= 66.0


def test_ranking_usa_o_dimero_com_ponta_3():
    cfg = carregar_config(CONFIG)
    cond = Condicoes.de_config(cfg)
    res = desenhar(ler_alvo(DADOS / "alvo_ok.fasta"), cfg)
    c1 = res.conjuntos[0]
    tubo = [c1.conjunto.F3, c1.conjunto.B3, c1.conjunto.FIP, c1.conjunto.BIP]
    pares = [dimero_ponta3(tubo[i], tubo[j], cond).dg_kcal for i in range(4) for j in range(i + 1, 4)]
    assert c1.pior_heterodimero_kcal == pytest.approx(min(pares))
    assert c1.pior_heterodimero_kcal == pytest.approx(-0.0282, abs=1e-3)
    assert c1.desvio_tm == pytest.approx(0.2762, abs=1e-3)
    assert c1.pontuacao == pytest.approx(0.2762 + 0.5 * 0.0282, abs=1e-3)


def test_proposta_1_no_alvo_ok():
    cfg = carregar_config(CONFIG)
    alvo = ler_alvo(DADOS / "alvo_ok.fasta")
    res = desenhar(alvo, cfg)
    assert res.estatisticas["candidatos_por_tipo"] == {  # fase 4: teto por janela
        "F3": 35, "F2": 35, "F1": 22, "B1c": 20, "B2c": 32, "B3c": 32,
    }
    c1 = res.conjuntos[0]
    assert c1.regioes == {
        "F3": Regiao(1, 26),
        "F2": Regiao(26, 49),
        "F1": Regiao(75, 99),
        "B1c": Regiao(105, 128),
        "B2c": Regiao(145, 169),
        "B3c": Regiao(194, 213),
    }
    assert c1.conjunto.F2 == "AGACACGCATATGGATTCTGTAC"
    assert c1.conjunto.F1c == "ACCCCCAAGTGTTATTGGCAAGTT"
    assert c1.avisos == ["sem LF candidato", "sem LB candidato"]
    assert verificar(c1.conjunto, alvo, cfg) == []


def test_saida_exporta_tm_da_reacao_e_as_duas_convencoes(tmp_path):
    cfg = carregar_config(CONFIG)
    leitura = ler_conjunto(DADOS / "conjunto_ok.csv")
    res = analisar(leitura, ler_alvo(DADOS / "alvo_ok.fasta"), cfg)
    exportar(res, tmp_path, cfg)

    with open(tmp_path / "primers.csv", encoding="utf-8") as f:
        linhas = list(csv.DictReader(f))
    assert list(linhas[0].keys())[4:6] == ["tm_c", "tm_reacao_c"]
    f2 = next(linha for linha in linhas if linha["nome"] == "F2")
    assert (f2["tm_c"], f2["tm_reacao_c"]) == ("51.27", "59.09")

    dados = json.loads((tmp_path / "resultado.json").read_text(encoding="utf-8"))
    assert dados["condicoes"] == {
        "reacao": {
            "na_mM": 0.0,
            "k_mM": 50.0,
            "nh4_mM": 20.0,
            "monovalentes_mM": 70.0,
            "mg_mM": 8.0,
            "dntp_cada_mM": 1.4,
            "dntp_total_mM": 5.6,
            "mg_livre_mM": 2.4,
            "temp_C": 65.0,
        },
        "referencia": {"convencao": "PrimerExplorer V5", "oligo_nM": 100.0, "na_mM": 50.0, "mg_mM": 4.0},
    }
    assert set(dados["heterodimeros"].keys()) == {"nomes", "dg_kcal", "dg_qualquer_kcal"}
