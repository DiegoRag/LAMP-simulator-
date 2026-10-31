"""Testes de alvo contendo bases N e janela de cadeia para evitar explosão combinatória."""

import random
import time
from pathlib import Path

from primers.analise import analisar
from primers.casamento import (
    base_casa_no_alvo,
    ocorrencias_aproximadas,
    ocorrencias_compativeis,
)
from primers.entrada import carregar_config, ler_alvo, ler_conjunto
from primers.modelo import Regiao

RAIZ = Path(__file__).resolve().parents[1]
DADOS = RAIZ / "tests" / "dados"
CONFIG = RAIZ / "config" / "primers.toml"


def test_alvo_com_n_e_flancos_aleatorios():
    """Alvo com 3 kb aleatórios + N*90 + alvo_ok + N*90 + 3 kb aleatórios."""
    cfg = carregar_config(CONFIG)
    leitura = ler_conjunto(DADOS / "conjunto_ok.csv")
    alvo_ok = ler_alvo(DADOS / "alvo_ok.fasta")

    rng = random.Random(42)
    flanco_esq = "".join(rng.choices("ACGT", k=3000))
    flanco_dir = "".join(rng.choices("ACGT", k=3000))
    alvo = flanco_esq + ("N" * 90) + alvo_ok + ("N" * 90) + flanco_dir

    deslocamento = len(flanco_esq) + 90

    t0 = time.perf_counter()
    resultado = analisar(leitura, alvo, cfg)
    duracao = time.perf_counter() - t0

    assert duracao < 5.0, f"Análise demorou {duracao:.2f} s (esperado < 5 s)"

    # Localiza as 6 regiões nas mesmas posições relativas do alvo_ok
    regioes_esperadas = {
        "F3": Regiao(0 + deslocamento, 21 + deslocamento),
        "F2": Regiao(31 + deslocamento, 49 + deslocamento),
        "F1": Regiao(81 + deslocamento, 103 + deslocamento),
        "B1c": Regiao(108 + deslocamento, 127 + deslocamento),
        "B2c": Regiao(157 + deslocamento, 177 + deslocamento),
        "B3c": Regiao(187 + deslocamento, 206 + deslocamento),
    }
    assert resultado.regioes is not None
    for nome, reg in regioes_esperadas.items():
        assert resultado.regioes[nome] == reg

    # Tem o aviso de bases N
    aviso_esperado = (
        "O alvo tem 180 bases N (desconhecidas) em 2 trecho(s); "
        "nenhum primer é localizado sobre elas."
    )
    assert aviso_esperado in resultado.avisos


def test_ocorrencias_compativeis_com_n():
    """N no alvo não casa; primer com N contra alvo ACGT continua casando."""
    # Nenhuma ocorrência de sub contra texto quando alvo é N
    assert ocorrencias_compativeis("ACGT" * 5, "N" * 40) == []
    assert ocorrencias_compativeis("N" * 40, "ACGT" * 5) == []

    # Primer com N contra alvo ACGT continua casando
    alvo_acgt = "ACGT" * 10
    ocorr = ocorrencias_compativeis(alvo_acgt, "ACGN")
    assert len(ocorr) == 10

    # Teste escalar base_casa_no_alvo
    assert not base_casa_no_alvo("A", "N")
    assert not base_casa_no_alvo("N", "N")
    assert base_casa_no_alvo("N", "A")
    assert base_casa_no_alvo("A", "A")
    assert not base_casa_no_alvo("A", "C")


def test_ocorrencias_aproximadas_com_n():
    """No texto com bloco de N, cada N conta como erro/incompatibilidade."""
    # seq = "ACGT", texto = "ACNT" -> 1 erro (N conta como erro)
    assert ocorrencias_aproximadas("ACGT", "ACNT", max_erros=1, max_ini=1) == [(0, 1)]
    assert ocorrencias_aproximadas("ACGT", "ACNT", max_erros=0, max_ini=1) == []

    # Bloco com 3 Ns no meio
    seq = "ACGTACGTACGT"  # 12 nt
    texto = "ACGTNNNTACGT"  # 3 Ns
    assert ocorrencias_aproximadas(seq, texto, max_erros=2, max_ini=1) == []
    assert ocorrencias_aproximadas(seq, texto, max_erros=3, max_ini=1) == [(0, 3)]


def test_alvo_repetido_duas_vezes():
    """Alvo_ok repetido 2 vezes separado por 5 kb aleatórios: < 5 s e escolhe a 1a cópia."""
    cfg = carregar_config(CONFIG)
    leitura = ler_conjunto(DADOS / "conjunto_ok.csv")
    alvo_ok = ler_alvo(DADOS / "alvo_ok.fasta")

    rng = random.Random(123)
    separador = "".join(rng.choices("ACGT", k=5000))
    alvo_duplo = alvo_ok + separador + alvo_ok

    t0 = time.perf_counter()
    resultado = analisar(leitura, alvo_duplo, cfg)
    duracao = time.perf_counter() - t0

    assert duracao < 5.0, f"Análise demorou {duracao:.2f} s (esperado < 5 s)"

    # Escolhe a primeira cópia (mesmas posições de alvo_ok)
    assert resultado.regioes is not None
    assert resultado.regioes["F3"] == Regiao(0, 21)
    assert resultado.regioes["F2"] == Regiao(31, 49)
    assert resultado.regioes["F1"] == Regiao(81, 103)
    assert resultado.regioes["B1c"] == Regiao(108, 127)
    assert resultado.regioes["B2c"] == Regiao(157, 177)
    assert resultado.regioes["B3c"] == Regiao(187, 206)
