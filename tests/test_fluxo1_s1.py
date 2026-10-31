"""Testes da Sessão 1 (Fase 1): Casamento aproximado, localização nas duas fitas,
ocorrências múltiplas, alvos circulares e distâncias.
"""

import threading
import tomllib
from pathlib import Path

import pytest

from interface.sessao import Sessao
from interface.sessoes import PADRAO_SESSOES, Sessoes, config_sessoes
from primers.casamento import Posicao, distancia_3p, melhor_posicao
from primers.entrada import carregar_config, ler_alvo, ler_conjunto
from primers.modelo import Regiao, rc
from primers.restricoes import Violacao, calcular_distancias, localizar_detalhado, verificar

RAIZ = Path(__file__).resolve().parents[1]
PADRAO = {"cookie": "lampsim_sessao", "validade_min": 120.0, "max_sessoes": 500,
          "cookie_seguro": False}
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


# --- Registro de sessões por cookie (fase do fluxo guiado, s1) ---


def _criar(validade_min: float = 10, max_sessoes: int = 100) -> tuple[Sessoes, list[float]]:
    tempo = [1000.0]
    return Sessoes(lambda: Sessao(RAIZ, {}, {}), validade_min, max_sessoes,
                   lambda: tempo[0]), tempo


@pytest.mark.parametrize("token", [None, "", "desconhecido"])
def test_sessoes_token_novo(token: str | None) -> None:
    sessoes, _ = _criar()
    entrada, nova = sessoes.obter(token)
    assert nova is True
    assert entrada.token != token
    assert len(entrada.token) >= 40


def test_sessoes_mesmo_token_mesma_sessao() -> None:
    sessoes, _ = _criar()
    primeira, _ = sessoes.obter(None)
    segunda, nova = sessoes.obter(primeira.token)
    assert nova is False
    assert segunda.sessao is primeira.sessao


def test_sessoes_tokens_diferentes_sessoes_diferentes() -> None:
    sessoes, _ = _criar()
    a, _ = sessoes.obter(None)
    b, _ = sessoes.obter(None)
    assert a.token != b.token
    assert a.sessao is not b.sessao


def test_sessoes_expiracao() -> None:
    sessoes, tempo = _criar(validade_min=10)
    entrada, _ = sessoes.obter(None)
    tempo[0] += 599
    mesma, nova = sessoes.obter(entrada.token)
    assert nova is False and mesma.sessao is entrada.sessao
    tempo[0] += 599  # o uso renovou o prazo
    mesma, nova = sessoes.obter(entrada.token)
    assert nova is False and mesma.sessao is entrada.sessao
    tempo[0] += 601
    outra, nova = sessoes.obter(entrada.token)
    assert nova is True
    assert outra.sessao is not entrada.sessao
    assert outra.token != entrada.token
    assert len(sessoes) == 1


def test_sessoes_limite_descarta_menos_usada() -> None:
    sessoes, tempo = _criar(max_sessoes=2)
    a, _ = sessoes.obter(None)
    tempo[0] += 1
    b, _ = sessoes.obter(None)
    tempo[0] += 1
    sessoes.obter(a.token)  # a passa a ser a mais recente
    tempo[0] += 1
    c, _ = sessoes.obter(None)
    assert len(sessoes) == 2
    assert sessoes.obter(a.token)[1] is False
    assert sessoes.obter(c.token)[1] is False
    assert sessoes.obter(b.token)[1] is True


def test_sessoes_mensagem_sai_uma_vez() -> None:
    sessoes, _ = _criar()
    entrada, _ = sessoes.obter(None)
    assert sessoes.tomar_mensagem(entrada.token) is None
    sessoes.definir_mensagem(entrada.token, {"tipo": "erro", "texto": "falhou"})
    assert sessoes.tomar_mensagem(entrada.token) == {"tipo": "erro", "texto": "falhou"}
    assert sessoes.tomar_mensagem(entrada.token) is None


def test_sessoes_mensagem_token_inexistente() -> None:
    sessoes, _ = _criar()
    sessoes.definir_mensagem("nao-existe", {"tipo": "info", "texto": "x"})
    assert sessoes.tomar_mensagem("nao-existe") is None
    assert len(sessoes) == 0


def test_config_sessoes_vazia_e_padrao() -> None:
    assert config_sessoes({}) == PADRAO
    assert PADRAO_SESSOES["cookie"] == PADRAO["cookie"]


def test_config_sessoes_valida() -> None:
    cfg = {"sessoes": {"cookie": "meu", "validade_min": 30, "max_sessoes": 7,
                       "cookie_seguro": True}}
    assert config_sessoes(cfg) == {"cookie": "meu", "validade_min": 30.0, "max_sessoes": 7,
                                   "cookie_seguro": True}


def test_config_sessoes_invalida_cai_no_padrao() -> None:
    cfg = {"sessoes": {"cookie": "", "validade_min": "10", "max_sessoes": True,
                       "cookie_seguro": "sim"}}
    assert config_sessoes(cfg) == PADRAO
    assert config_sessoes({"sessoes": {"validade_min": -5, "max_sessoes": -1}}) == PADRAO
    assert config_sessoes({"sessoes": {"validade_min": True, "max_sessoes": 2.5}}) == PADRAO
    assert config_sessoes({"sessoes": "texto"}) == PADRAO


def test_config_sessoes_do_web_toml_real() -> None:
    with (RAIZ / "config" / "web.toml").open("rb") as f:
        cfg = tomllib.load(f)
    assert config_sessoes(cfg)["cookie"] == "lampsim_sessao"


def test_sessoes_concorrencia_simples() -> None:
    sessoes = Sessoes(lambda: Sessao(RAIZ, {}, {}), 120, 100)
    erros: list[Exception] = []

    def trabalhar() -> None:
        try:
            for _ in range(50):
                sessoes.obter(None)
        except Exception as e:  # noqa: BLE001
            erros.append(e)

    threads = [threading.Thread(target=trabalhar) for _ in range(8)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert not erros
    assert len(sessoes) <= 100
