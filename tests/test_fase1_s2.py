"""Testes para a Sessão 2 (Fase 1): separação de FIP/BIP, geometria e detecção de ligantes."""

from pathlib import Path

import pytest

from primers.analise import analisar
from primers.entrada import carregar_config, ler_alvo, ler_conjunto
from primers.modelo import rc

DADOS = Path(__file__).parent / "dados"
CONFIG = Path(__file__).parents[1] / "config" / "primers.toml"


def test_separar_fip_com_snp_em_f1(tmp_path):
    alvo_ok = ler_alvo(DADOS / "alvo_ok.fasta")
    assert alvo_ok[90] == "C"
    alvo_snp = alvo_ok[:90] + "A" + alvo_ok[91:]

    f1c = "GTAAACCCCCAAGTGTTATTGG"
    f2 = "CGCATATGGATTCTGTAC"
    b1c = "AAGGGTAATTCTCCGTCCA"
    b2 = "GCAATAGTATAAAGTCTCGC"

    f_csv = tmp_path / "fip_snp.csv"
    f_csv.write_text(
        f"nome,sequencia\n"
        f"F3,GTGACAGACTCTCCTTAGATA\n"
        f"B3,TGTTGCTGGATTGTCTGTT\n"
        f"FIP,{f1c + f2}\n"
        f"BIP,{b1c + b2}\n",
        encoding="utf-8",
    )

    leitura = ler_conjunto(f_csv, alvo=alvo_snp)

    assert leitura.conjunto is not None
    assert leitura.conjunto.F1c == "GTAAACCCCCAAGTGTTATTGG"
    assert leitura.conjunto.F2 == "CGCATATGGATTCTGTAC"
    assert leitura.conjunto.B1c == "AAGGGTAATTCTCCGTCCA"
    assert leitura.conjunto.B2 == "GCAATAGTATAAAGTCTCGC"
    assert any("casamento aproximado" in a and "FIP" in a for a in leitura.avisos)


def test_separar_fip_bip_com_ligante_tttt_notomi(tmp_path):
    alvo_ok = ler_alvo(DADOS / "alvo_ok.fasta")
    f1c = "GTAAACCCCCAAGTGTTATTGG"
    f2 = "CGCATATGGATTCTGTAC"
    b1c = "AAGGGTAATTCTCCGTCCA"
    b2 = "GCAATAGTATAAAGTCTCGC"

    f_csv = tmp_path / "fip_bip_tttt.csv"
    f_csv.write_text(
        f"nome,sequencia\n"
        f"F3,GTGACAGACTCTCCTTAGATA\n"
        f"B3,TGTTGCTGGATTGTCTGTT\n"
        f"FIP,{f1c + 'TTTT' + f2}\n"
        f"BIP,{b1c + 'TTTT' + b2}\n",
        encoding="utf-8",
    )

    leitura = ler_conjunto(f_csv, alvo=alvo_ok, ligante="")

    assert leitura.conjunto is not None
    assert leitura.conjunto.ligante == "TTTT"
    assert leitura.conjunto.F1c == "GTAAACCCCCAAGTGTTATTGG"
    assert leitura.conjunto.F2 == "CGCATATGGATTCTGTAC"
    assert leitura.conjunto.B1c == "AAGGGTAATTCTCCGTCCA"
    assert leitura.conjunto.B2 == "GCAATAGTATAAAGTCTCGC"


def test_detectar_ligante_quando_partes_e_fip_informados(tmp_path):
    f1c = "GTAAACCCCCAAGTGTTATTGG"
    f2 = "CGCATATGGATTCTGTAC"
    b1c = "AAGGGTAATTCTCCGTCCA"
    b2 = "GCAATAGTATAAAGTCTCGC"

    f_csv = tmp_path / "partes_e_fip_tttt.csv"
    f_csv.write_text(
        f"nome,sequencia\n"
        f"F3,GTGACAGACTCTCCTTAGATA\n"
        f"B3,TGTTGCTGGATTGTCTGTT\n"
        f"F1c,{f1c}\n"
        f"F2,{f2}\n"
        f"B1c,{b1c}\n"
        f"B2,{b2}\n"
        f"FIP,{f1c + 'TTTT' + f2}\n"
        f"BIP,{b1c + 'TTTT' + b2}\n",
        encoding="utf-8",
    )

    leitura = ler_conjunto(f_csv, ligante="")

    assert leitura.conjunto is not None
    assert leitura.conjunto.ligante == "TTTT"


def test_separar_fip_bip_em_alvo_na_fita_menos(tmp_path):
    alvo_ok = ler_alvo(DADOS / "alvo_ok.fasta")
    alvo_menos = rc(alvo_ok)
    cfg = carregar_config(CONFIG)

    f1c = "GTAAACCCCCAAGTGTTATTGG"
    f2 = "CGCATATGGATTCTGTAC"
    b1c = "AAGGGTAATTCTCCGTCCA"
    b2 = "GCAATAGTATAAAGTCTCGC"

    f_csv = tmp_path / "fita_menos.csv"
    f_csv.write_text(
        f"nome,sequencia\n"
        f"F3,GTGACAGACTCTCCTTAGATA\n"
        f"B3,TGTTGCTGGATTGTCTGTT\n"
        f"FIP,{f1c + f2}\n"
        f"BIP,{b1c + b2}\n",
        encoding="utf-8",
    )

    leitura = ler_conjunto(f_csv, alvo=alvo_menos)
    assert leitura.conjunto is not None
    assert leitura.conjunto.F1c == f1c
    assert leitura.conjunto.F2 == f2
    assert leitura.conjunto.B1c == b1c
    assert leitura.conjunto.B2 == b2

    res = analisar(leitura, alvo_menos, cfg)
    assert res.fita == "-"
    # Tm e ponta de conjunto_ok reprovam na convenção da fase 3; aqui só a geometria importa
    assert [v for v in res.violacoes if v.regra not in ("tm_c", "extremidade")] == []
    assert any("fita complementar" in a for a in res.avisos)


def test_ligantes_diferentes_fip_bip_erro(tmp_path):
    f1c = "GTAAACCCCCAAGTGTTATTGG"
    f2 = "CGCATATGGATTCTGTAC"
    b1c = "AAGGGTAATTCTCCGTCCA"
    b2 = "GCAATAGTATAAAGTCTCGC"

    f_csv = tmp_path / "ligantes_diferentes.csv"
    f_csv.write_text(
        f"nome,sequencia\n"
        f"F3,GTGACAGACTCTCCTTAGATA\n"
        f"B3,TGTTGCTGGATTGTCTGTT\n"
        f"F1c,{f1c}\n"
        f"F2,{f2}\n"
        f"B1c,{b1c}\n"
        f"B2,{b2}\n"
        f"FIP,{f1c + 'TTTT' + f2}\n"
        f"BIP,{b1c + 'TTT' + b2}\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError) as exc:
        ler_conjunto(f_csv, ligante="")
    assert "Ligante do FIP" in str(exc.value)
