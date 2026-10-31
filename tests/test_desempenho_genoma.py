"""Testes de desempenho e exatidão para genoma inteiro (s1)."""

import random
import time
from pathlib import Path

import pytest
from Bio.Seq import Seq

from primers.casamento import (
    Diferenca,
    Posicao,
    _buscar_em_fita,
    _indice_kmers,
    base_casa_no_alvo,
    ocorrencias_aproximadas,
)
from primers.entrada import carregar_config, ler_alvo, ler_conjunto
from primers.modelo import BASES_IUPAC, Regiao, bases_compativeis, normalizar, rc
from primers.restricoes import localizar_detalhado, verificar

DADOS = Path(__file__).parent / "dados"
CONFIG = Path(__file__).parents[1] / "config" / "primers.toml"


def _ocorrencias_aproximadas_referencia(
    seq: str, texto: str, max_erros: int, max_ini: int
) -> list[tuple[int, int]]:
    """Referência ingênua para teste de ocorrências aproximadas."""
    k = len(seq)
    if k == 0 or k > len(texto) or max_ini <= 0 or max_erros < 0:
        return []
    lim = min(max_ini, len(texto) - k + 1)
    res: list[tuple[int, int]] = []
    for idx in range(lim):
        sub = texto[idx : idx + k]
        erros = sum(1 for a, b in zip(seq, sub, strict=True) if not base_casa_no_alvo(a, b))
        if erros <= max_erros:
            res.append((idx, erros))
    res.sort(key=lambda x: (x[1], x[0]))
    return res


def _buscar_em_fita_referencia(seq: str, texto: str, max_erros: int) -> Posicao | None:
    """Referência ingênua para busca de melhor posição em fita."""
    k = len(seq)
    menor_dif: int | None = None
    melhor_pos: Posicao | None = None
    for idx in range(len(texto) - k + 1):
        sub = texto[idx : idx + k]
        erros = sum(1 for a, b in zip(seq, sub, strict=True) if not base_casa_no_alvo(a, b))
        if erros <= max_erros:
            if menor_dif is None or erros < menor_dif:
                menor_dif = erros
                difs = tuple(
                    Diferenca(posicao=idx + j, esperado=seq[j], no_alvo=texto[idx + j])
                    for j in range(k)
                    if not base_casa_no_alvo(seq[j], texto[idx + j])
                )
                melhor_pos = Posicao(regiao=Regiao(idx, idx + k), diferencas=difs, fita="+")
                if menor_dif == 0:
                    break
    return melhor_pos


def test_ocorrencias_aproximadas_30_casos_aleatorios():
    """Valida ocorrencias_aproximadas contra referência ingênua em 30 casos sintéticos."""
    rng = random.Random(12345)
    bases_iupac_lista = list("ACGTRYSWKMBDHVN")

    for caso in range(30):
        tam_texto = rng.randint(300, 3000)
        tam_seq = rng.randint(10, 25)
        max_erros = rng.randint(0, 3)

        if caso % 2 == 0:
            texto_chars = [rng.choice("ACGT") for _ in range(tam_texto)]
        else:
            texto_chars = [rng.choice(bases_iupac_lista) for _ in range(tam_texto)]

        seq_chars = [rng.choice("ACGT") for _ in range(tam_seq)]

        max_ini = rng.randint(tam_seq, tam_texto)
        lim = min(max_ini, tam_texto - tam_seq + 1)
        pos_plantio = rng.randint(0, lim - 1)

        sub_plantada = list(seq_chars)
        erros_plantar = rng.randint(0, 3)
        pos_mutar = rng.sample(range(tam_seq), min(erros_plantar, tam_seq))
        for p in pos_mutar:
            candidatos_incomp = [b for b in "ACGT" if not bases_compativeis(b, seq_chars[p])]
            if candidatos_incomp:
                sub_plantada[p] = rng.choice(candidatos_incomp)

        texto_chars[pos_plantio : pos_plantio + tam_seq] = sub_plantada

        texto = "".join(texto_chars)
        seq = "".join(seq_chars)

        esperado = _ocorrencias_aproximadas_referencia(seq, texto, max_erros, max_ini)
        obtido = ocorrencias_aproximadas(seq, texto, max_erros, max_ini)

        assert obtido == esperado, f"Divergência no caso {caso}: obtido={obtido}, esperado={esperado}"


def test_normalizar_genoma_sintetico_e_invalido():
    """Valida normalizar em genoma longo sintético e detecção de caractere inválido."""
    rng = random.Random(42)
    tam_genoma = 5000
    alfabeto = "ACGTMRYSWKMBDHVN"
    seq_valida = "".join(rng.choice(alfabeto) for _ in range(tam_genoma))

    norm = normalizar(seq_valida)
    assert len(norm) == tam_genoma

    pos_invalida = rng.randint(1, tam_genoma)
    seq_invalida = seq_valida[: pos_invalida - 1] + "X" + seq_valida[pos_invalida:]

    with pytest.raises(ValueError) as excinfo:
        normalizar(seq_invalida)

    msg_esperada = (
        f"Caractere inválido 'X' na posição {pos_invalida} "
        "(esperado base IUPAC A, C, G, T, U ou N)."
    )
    assert msg_esperada in str(excinfo.value)


def test_rc_devolve_o_mesmo_que_biopython():
    """Valida que rc retorna o mesmo que Bio.Seq.reverse_complement."""
    rng = random.Random(999)
    for _ in range(25):
        tam = rng.randint(10, 500)
        seq = "".join(rng.choice("ACGT") for _ in range(tam))
        assert rc(seq) == str(Seq(seq).reverse_complement())

    seq_iupac = "ACGTUNRYSWKMBDHV"
    assert rc(seq_iupac) == str(Seq(seq_iupac.replace("U", "T")).reverse_complement())


def test_buscar_em_fita_ramo_pombos_20_casos():
    """Valida ramo dos pombos de _buscar_em_fita em 20 casos aleatórios com texto de 5000 nt."""
    rng = random.Random(54321)

    for caso in range(20):
        tam_texto = 5000
        texto_chars = [rng.choice("ACGT") for _ in range(tam_texto)]

        max_erros = rng.randint(1, 3)
        k = rng.randint((max_erros + 1) * 3, 30)

        seq_chars = [rng.choice("ACGT") for _ in range(k)]
        pos_plantio = rng.randint(0, tam_texto - k)

        erros_plantar = rng.randint(0, max_erros)
        sub_plantada = list(seq_chars)
        pos_mutar = rng.sample(range(k), erros_plantar)
        for p in pos_mutar:
            opcoes = [b for b in "ACGT" if b != seq_chars[p]]
            sub_plantada[p] = rng.choice(opcoes)

        texto_chars[pos_plantio : pos_plantio + k] = sub_plantada

        texto = "".join(texto_chars)
        seq = "".join(seq_chars)

        esperado = _buscar_em_fita_referencia(seq, texto, max_erros)
        obtido = _buscar_em_fita(seq, texto, "+", None, max_erros)

        assert obtido == esperado, f"Divergência no caso {caso}: obtido={obtido}, esperado={esperado}"


def test_verificar_com_loc_reaproveitado():
    """Valida que verificar(..., loc=loc) devolve o mesmo que sem loc para dados padrão."""
    cfg = carregar_config(CONFIG)
    leitura = ler_conjunto(DADOS / "conjunto_ok.csv")
    assert leitura.conjunto is not None
    alvo = ler_alvo(DADOS / "alvo_ok.fasta")

    loc = localizar_detalhado(leitura.conjunto, alvo, cfg)
    v_sem_loc = verificar(leitura.conjunto, alvo, cfg)
    v_com_loc = verificar(leitura.conjunto, alvo, cfg, loc=loc)

    assert v_com_loc == v_sem_loc


def test_ler_conjunto_fip_bip_inteiros_em_alvo_200kb(tmp_path):
    """Valida separação de FIP/BIP inteiros em alvo sintético de 200 kb com amplicon plantado."""
    rng = random.Random(42)
    seq_alvo_ok = ler_alvo(DADOS / "alvo_ok.fasta")
    leitura_ref = ler_conjunto(DADOS / "conjunto_ok.csv")
    assert leitura_ref.conjunto is not None
    conj = leitura_ref.conjunto

    tam_total = 200_000
    pos_plantio = rng.randint(10_000, tam_total - len(seq_alvo_ok) - 10_000)
    flanco_esq = "".join(rng.choice("ACGT") for _ in range(pos_plantio))
    flanco_dir = "".join(rng.choice("ACGT") for _ in range(tam_total - pos_plantio - len(seq_alvo_ok)))
    alvo_200k = flanco_esq + seq_alvo_ok + flanco_dir
    assert len(alvo_200k) == tam_total

    fip_int = conj.F1c + conj.F2
    bip_int = conj.B1c + conj.B2
    csv_conteudo = (
        "nome,sequencia\n"
        f"F3,{conj.F3}\n"
        f"B3,{conj.B3}\n"
        f"FIP,{fip_int}\n"
        f"BIP,{bip_int}\n"
    )
    arq_csv = tmp_path / "conjunto_fip_bip_inteiros.csv"
    arq_csv.write_text(csv_conteudo, encoding="utf-8")

    leitura = ler_conjunto(arq_csv, alvo=alvo_200k)
    assert leitura.conjunto is not None
    assert leitura.conjunto.F1c == conj.F1c
    assert leitura.conjunto.F2 == conj.F2
    assert leitura.conjunto.B1c == conj.B1c
    assert leitura.conjunto.B2 == conj.B2
    assert "FIP" in leitura.separados
    assert "BIP" in leitura.separados


def test_ler_conjunto_fip_bip_caminho_aproximado_em_alvo_200kb(tmp_path):
    """Valida separação de FIP/BIP pelo caminho aproximado com 1 erro plantado no F2."""
    rng = random.Random(99)
    seq_alvo_ok = ler_alvo(DADOS / "alvo_ok.fasta")
    leitura_ref = ler_conjunto(DADOS / "conjunto_ok.csv")
    assert leitura_ref.conjunto is not None
    conj = leitura_ref.conjunto

    tam_total = 200_000
    pos_plantio = rng.randint(10_000, tam_total - len(seq_alvo_ok) - 10_000)
    flanco_esq = "".join(rng.choice("ACGT") for _ in range(pos_plantio))
    flanco_dir = "".join(rng.choice("ACGT") for _ in range(tam_total - pos_plantio - len(seq_alvo_ok)))
    alvo_200k = flanco_esq + seq_alvo_ok + flanco_dir

    f2_chars = list(conj.F2)
    # Planta 1 erro na base 8 do F2 (troca G por A, que não casa com o alvo naquele ponto)
    f2_chars[8] = "A" if f2_chars[8] != "A" else "C"
    f2_mutado = "".join(f2_chars)

    fip_int_mutado = conj.F1c + f2_mutado
    bip_int = conj.B1c + conj.B2
    csv_conteudo = (
        "nome,sequencia\n"
        f"F3,{conj.F3}\n"
        f"B3,{conj.B3}\n"
        f"FIP,{fip_int_mutado}\n"
        f"BIP,{bip_int}\n"
    )
    arq_csv = tmp_path / "conjunto_fip_bip_aproximado.csv"
    arq_csv.write_text(csv_conteudo, encoding="utf-8")

    leitura = ler_conjunto(arq_csv, alvo=alvo_200k)
    assert leitura.conjunto is not None
    assert leitura.conjunto.F1c == conj.F1c
    assert leitura.conjunto.F2 == f2_mutado
    assert leitura.conjunto.B1c == conj.B1c
    assert leitura.conjunto.B2 == conj.B2
    assert "FIP" in leitura.separados
    assert "BIP" in leitura.separados
    assert any("FIP separado por casamento aproximado" in aviso for aviso in leitura.avisos)


def test_ocorrencias_aproximadas_indice_kmers_15_casos():
    """Valida ocorrencias_aproximadas com blocos de 4–5 nt (max_erros=3) contra referência."""
    rng = random.Random(1337)
    bases_iupac_lista = list("ACGTRYSWKMBDHV")

    for caso in range(15):
        tam_texto = rng.randint(25000, 60000)
        tam_seq = rng.randint(16, 22)
        max_erros = 3

        texto_chars = [rng.choice("ACGT") for _ in range(tam_texto)]

        # Casos com e sem blocos de 'N' no texto
        if caso % 2 == 1:
            num_blocos_n = rng.randint(1, 3)
            for _ in range(num_blocos_n):
                tam_n = rng.randint(10, 80)
                pos_n = rng.randint(0, tam_texto - tam_n)
                texto_chars[pos_n : pos_n + tam_n] = ["N"] * tam_n

        # Sequências com e sem IUPAC
        if caso % 3 == 0:
            seq_chars = [rng.choice(bases_iupac_lista) for _ in range(tam_seq)]
        else:
            seq_chars = [rng.choice("ACGT") for _ in range(tam_seq)]

        max_ini = rng.randint(tam_seq, tam_texto)
        lim = min(max_ini, tam_texto - tam_seq + 1)
        pos_plantio = rng.randint(0, lim - 1)

        # Monta sequência a ser plantada no texto (garante que texto só tem ACGT e N)
        sub_plantada: list[str] = []
        for ch in seq_chars:
            if ch in "ACGT":
                sub_plantada.append(ch)
            else:
                sub_plantada.append(rng.choice(list(BASES_IUPAC[ch])))

        erros_plantar = rng.randint(0, max_erros)
        pos_mutar = rng.sample(range(tam_seq), min(erros_plantar, tam_seq))
        for p in pos_mutar:
            candidatos_incomp = [b for b in "ACGT" if not bases_compativeis(b, seq_chars[p])]
            if candidatos_incomp:
                sub_plantada[p] = rng.choice(candidatos_incomp)

        texto_chars[pos_plantio : pos_plantio + tam_seq] = sub_plantada

        texto = "".join(texto_chars)
        seq = "".join(seq_chars)

        esperado = _ocorrencias_aproximadas_referencia(seq, texto, max_erros, max_ini)
        obtido = ocorrencias_aproximadas(seq, texto, max_erros, max_ini)

        assert obtido == esperado, f"Divergência no caso {caso}: obtido={obtido}, esperado={esperado}"


def test_indice_kmers_texto_com_r_devolve_none_e_mantem_referencia():
    """Valida que texto com 'R' faz _indice_kmers devolver None e mantém resultado correto."""
    rng = random.Random(888)
    tam_texto = 30000
    texto_chars = [rng.choice("ACGT") for _ in range(tam_texto)]
    # Coloca 'R' no meio do texto
    texto_chars[tam_texto // 2] = "R"
    texto = "".join(texto_chars)

    assert _indice_kmers(texto) is None
    assert _indice_kmers(texto, 4) is None
    assert _indice_kmers(texto, 6) is None

    seq = "".join(rng.choice("ACGT") for _ in range(20))
    max_erros = 2
    max_ini = tam_texto

    esperado = _ocorrencias_aproximadas_referencia(seq, texto, max_erros, max_ini)
    obtido = ocorrencias_aproximadas(seq, texto, max_erros, max_ini)
    assert obtido == esperado


def test_indice_kmers_curto_e_longo_estrutura():
    """Valida tamanhos e tipos das estruturas devolvidas por _indice_kmers para k=4 e k=6."""
    texto = "ACGT" * 6000  # 24.000 nt
    idx_4 = _indice_kmers(texto, 4)
    assert idx_4 is not None
    pos_4, fat_4 = idx_4
    assert len(fat_4) == 4**4 + 1
    assert len(pos_4) == len(texto) - 4 + 1

    idx_6 = _indice_kmers(texto, 6)
    assert idx_6 is not None
    pos_6, fat_6 = idx_6
    assert len(fat_6) == 4**6 + 1
    assert len(pos_6) == len(texto) - 6 + 1


def test_desempenho_ocorrencias_aproximadas_2mb_50_chamadas():
    """Valida que 50 chamadas em texto aleatório de 2 Mb terminam em menos de 5 s."""
    rng = random.Random(9999)
    tam_texto = 2_000_000
    texto = "".join(rng.choices("ACGT", k=tam_texto))
    # Sequência de 18 nt com max_erros = 3 (gera blocos de 4 e 5 nt)
    seq = "AAAAATTTTTCCCCCGGG"
    max_erros = 3
    max_ini = tam_texto

    t0 = time.monotonic()
    for _ in range(50):
        res = ocorrencias_aproximadas(seq, texto, max_erros, max_ini)
    tempo_total = time.monotonic() - t0

    # o texto é aleatório: pode haver ocorrências ao acaso; cada uma tem de ter os erros certos
    for ini, erros in res:
        sub = texto[ini : ini + len(seq)]
        assert erros == sum(1 for a, b in zip(seq, sub, strict=True) if not base_casa_no_alvo(a, b))
        assert erros <= max_erros
    assert tempo_total < 5.0, f"Tempo excessivo: {tempo_total:.2f} s (esperado < 5 s)"
