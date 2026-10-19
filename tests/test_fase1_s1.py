"""Testes da Sessão 1 (Fase 1): Casamento aproximado, localização nas duas fitas,
ocorrências múltiplas, alvos circulares e distâncias.
"""

from pathlib import Path

from primers.casamento import Posicao, distancia_3p, melhor_posicao
from primers.entrada import carregar_config, ler_alvo, ler_conjunto
from primers.modelo import Regiao, rc
from primers.restricoes import Violacao, calcular_distancias, localizar_detalhado, verificar

DADOS = Path(__file__).parent / "dados"
CONFIG = Path(__file__).parents[1] / "config" / "primers.toml"


def _geometricas(violacoes: list[Violacao]) -> list[Violacao]:
    """Só as violações de geometria: as de Tm e de ponta de conjunto_ok mudaram na fase 3."""
    return [v for v in violacoes if v.regra not in ("tm_c", "extremidade")]


def test_localizar_fita_menos_alvo_invertido():
    cfg = carregar_config(CONFIG)
    leitura = ler_conjunto(DADOS / "conjunto_ok.csv")
    conjunto = leitura.conjunto
    alvo_ok = ler_alvo(DADOS / "alvo_ok.fasta")
    alvo_rc = rc(alvo_ok)

    loc = localizar_detalhado(conjunto, alvo_rc, cfg)
    violacoes = verificar(conjunto, alvo_rc, cfg)

    assert loc.fita == "-"
    assert _geometricas(violacoes) == []
    assert loc.regioes["F3"] == Regiao(0, 21)
    assert loc.regioes["F2"] == Regiao(31, 49)
    assert loc.regioes["F1"] == Regiao(81, 103)
    assert loc.regioes["B1c"] == Regiao(108, 127)
    assert loc.regioes["B2c"] == Regiao(157, 177)
    assert loc.regioes["B3c"] == Regiao(187, 206)
    assert any("fita complementar (−)" in a for a in loc.avisos)


def test_multiplas_ocorrencias_escolhe_combinacao_coerente():
    cfg = carregar_config(CONFIG)
    leitura = ler_conjunto(DADOS / "conjunto_ok.csv")
    conjunto = leitura.conjunto
    alvo_ok = ler_alvo(DADOS / "alvo_ok.fasta")

    f2 = conjunto.F2
    prefixo = f2 + ("A" * 115)
    alvo_rep = prefixo + alvo_ok

    assert _geometricas(verificar(conjunto, alvo_rep, cfg)) == []
    loc = localizar_detalhado(conjunto, alvo_rep, cfg)
    assert loc.regioes["F2"] == Regiao(164, 182)
    assert calcular_distancias(loc.regioes)["F2_B2"] == 146
    assert any("Região F2 possui 2 ocorrências" in a for a in loc.avisos)


def test_casamento_limiar_janela_e_duas_fitas():
    assert melhor_posicao("ACGTACGTACGTACGTAC", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", max_erros=3) is None
    assert melhor_posicao("AAACCCGGGTTTAAACCC", rc("AAACCCGGGTTTAAACCC"), ambas_fitas=True) == Posicao(
        Regiao(0, 18), (), fita="-"
    )

    # F2: primer 5'->3' na fita ativa, extremidade 3' em regiao.fim - 1
    reg_f2 = Regiao(31, 49)
    assert distancia_3p("F2", reg_f2, 48) == 0
    assert distancia_3p("F2", reg_f2, 45) == 3

    # LF: rc(primer) na fita ativa, extremidade 3' em regiao.ini
    reg_lf = Regiao(55, 73)
    assert distancia_3p("LF", reg_lf, 55) == 0
    assert distancia_3p("LF", reg_lf, 60) == 5


def test_regiao_com_um_erro_mantem_distancias_e_informa_3p():
    cfg = carregar_config(CONFIG)
    leitura = ler_conjunto(DADOS / "conjunto_ok.csv")
    conjunto = leitura.conjunto
    alvo_ok = ler_alvo(DADOS / "alvo_ok.fasta")

    # Troca a última base de F2 (posição 48 de alvo_ok, base 'C') por 'A'
    assert alvo_ok[48] == "C"
    alvo_1erro = alvo_ok[:48] + "A" + alvo_ok[49:]

    loc = localizar_detalhado(conjunto, alvo_1erro, cfg)
    violacoes = verificar(conjunto, alvo_1erro, cfg)

    assert "F2" not in loc.regioes
    assert "F2" in loc.aproximadas
    assert loc.aproximadas["F2"].regiao == Regiao(31, 49)
    assert calcular_distancias(loc.regioes_completas)["F2_B2"] == 146

    v_f2 = next(v for v in violacoes if v.item == "F2" and v.regra == "nao_encontrado")
    assert "posição 49" in v_f2.mensagem
    assert "a 0 nt da ponta 3'" in v_f2.mensagem


def test_f2_b2_ate_180_e_f3_b3_total():
    cfg = carregar_config(CONFIG)
    leitura = ler_conjunto(DADOS / "conjunto_ok.csv")
    conjunto = leitura.conjunto
    alvo_ok = ler_alvo(DADOS / "alvo_ok.fasta")

    # Alvo com F2_B2 = 173 (inserção de 27 A's entre F1 e B1c na posição 105)
    alvo_173 = alvo_ok[:105] + ("A" * 27) + alvo_ok[105:]
    assert _geometricas(verificar(conjunto, alvo_173, cfg)) == []

    # Alvo com F3_F2_gap = 60, F2_B2 = 170, B2c_B3c_gap = 60 -> F3_B3_total = 330
    alvo_330 = (
        alvo_ok[:25]
        + ("A" * 50)
        + alvo_ok[25:105]
        + ("A" * 24)
        + alvo_ok[105:180]
        + ("A" * 50)
        + alvo_ok[180:]
    )
    violacoes = verificar(conjunto, alvo_330, cfg)
    assert _geometricas(violacoes) == [
        Violacao(
            item="F3_B3_total",
            regra="distancia",
            valor=330.0,
            faixa="[120, 280]",
            mensagem="Distância F3_B3_total (330) fora da faixa [120, 280]",
        )
    ]


def test_alvo_circular_cruzando_origem():
    cfg = carregar_config(CONFIG)
    leitura = ler_conjunto(DADOS / "conjunto_ok.csv")
    conjunto = leitura.conjunto
    alvo_ok = ler_alvo(DADOS / "alvo_ok.fasta")

    alvo_circ = alvo_ok[100:] + alvo_ok[:100]

    # Com circular=False: amplicon quebrado na origem
    assert len(_geometricas(verificar(conjunto, alvo_circ, cfg))) > 0

    # Com circular=True: localizado perfeitamente
    loc = localizar_detalhado(conjunto, alvo_circ, cfg, circular=True)
    assert loc.cruza_origem is True
    assert any("cruza a origem do alvo circular" in a for a in loc.avisos)
    assert _geometricas(verificar(conjunto, alvo_circ, cfg, circular=True)) == []
    assert calcular_distancias(loc.regioes) == {
        "F2_B2": 146,
        "F2_F1": 50,
        "B1c_B2c": 50,
        "F3_F2_gap": 10,
        "B2c_B3c_gap": 10,
        "F1_B1c_gap": 5,
        "F3_B3_total": 206,
    }
