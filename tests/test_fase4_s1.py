"""Fase 4, sessão 1: janelas do alvo em lote, índice de k-meros e unicidade nas duas fitas."""

import random
from pathlib import Path

import numpy as np

from primers.entrada import carregar_config, ler_alvo
from primers.janelas import IndiceAlvo, codigos_kmer
from primers.modelo import rc
from primers.termo import CondicoesReferencia, dg_extremidade, gc_pct, tm_referencia

DADOS = Path(__file__).parent / "dados"
CONFIG = Path(__file__).parents[1] / "config" / "primers.toml"


def test_lote_igual_a_termo_bit_a_bit():
    cfg = carregar_config(CONFIG)
    ref = CondicoesReferencia.de_config(cfg)
    alvo = ler_alvo(DADOS / "alvo_ok.fasta")
    indice = IndiceAlvo(alvo, cfg)
    for L in (15, 18, 25):
        lt = indice.lote(L, 0, len(alvo) - L)
        assert len(lt["tm"]) == len(alvo) - L + 1
        for i in range(len(alvo) - L + 1):
            w = alvo[i : i + L]
            assert not lt["ambigua"][i]
            assert lt["gc"][i] == gc_pct(w)
            assert lt["tm"][i] == tm_referencia(w, ref)
            assert lt["tm_rc"][i] == tm_referencia(rc(w), ref)
            assert lt["dg_fim"][i] == dg_extremidade(w, "3", 6, ref)
            assert lt["dg_ini"][i] == dg_extremidade(w, "5", 6, ref)
            assert lt["dg_ini_rc"][i] == dg_extremidade(rc(w), "3", 6, ref)
            assert lt["dg_fim_rc"][i] == dg_extremidade(rc(w), "5", 6, ref)
    trecho = indice.lote(20, 100, 104)
    assert len(trecho["tm"]) == 5
    assert trecho["tm"][4] == tm_referencia(alvo[104:124], ref)


def test_codigos_kmer():
    cod = np.array([0, 1, 2, 3])  # ACGT
    assert codigos_kmer(cod, 2, reverso=False).tolist() == [1, 6, 11]  # AC, CG, GT
    assert codigos_kmer(cod, 2, reverso=True).tolist() == [11, 6, 1]  # rc: GT, CG, AC
    assert codigos_kmer(np.array([0, -1, 2, 3]), 2, reverso=False).tolist() == [-1, -1, 11]
    assert codigos_kmer(cod, 5, reverso=False).tolist() == []


def test_unica_nas_duas_fitas_contra_forca_bruta():
    cfg = carregar_config(CONFIG)
    rng = random.Random(3)
    base = "".join(rng.choices("ACGT", k=400))
    alvo = (
        base[:200] + base[50:70]  # repetição direta de [50, 70)
        + base[200:260] + rc(base[100:122])  # repetição invertida de [100, 122)
        + base[260:300] + "GAATTCGGATCCGGATCCGAATTC"  # palíndromo de 24 nt
        + "ACACACACACACACACACACACAC" + base[300:]  # repetição em tandem
    )
    indice = IndiceAlvo(alvo, cfg)
    assert indice.k == 15
    repetidas = 0
    for L in (15, 22):
        janelas = [alvo[j : j + L] for j in range(len(alvo) - L + 1)]
        for i, w in enumerate(janelas):
            outra = any(x == w for j, x in enumerate(janelas) if j != i)
            no_rc = rc(w) in janelas
            assert indice.unica(i, L) == (not outra and not no_rc), (i, L)
            repetidas += outra or no_rc
    assert repetidas == 64
    assert indice.unica(50, 20) is False  # cópia direta
    assert indice.unica(100, 22) is False  # cópia invertida, só visível no rc(alvo)
    assert indice.unica(alvo.index("GAATTCGGATCC"), 24) is False  # palíndromo
    assert indice.unica(0, 20) is True


def test_janela_ambigua_e_alvo_longo():
    cfg = carregar_config(CONFIG)
    alvo = ler_alvo(DADOS / "alvo_ok.fasta")
    lt = IndiceAlvo(alvo[:30] + "N" + alvo[31:], cfg).lote(20, 0, 40)
    assert lt["ambigua"].tolist() == [i <= 30 < i + 20 for i in range(41)]
    rng = random.Random(11)
    base = "".join(rng.choices("ACGT", k=100_000))
    indice = IndiceAlvo(base + rc(base[10:35]), cfg)  # cópia invertida a 100 kb da original
    assert indice.unica(10, 25) is False
    assert indice.unica(100_000, 25) is False
    assert indice.unica(40_000, 25) is True
    assert int(np.count_nonzero(indice.suspeito)) < 100
