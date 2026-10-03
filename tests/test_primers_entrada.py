"""Testes para o módulo de entrada de primers LAMP."""

from pathlib import Path
import pytest
from primers.entrada import carregar_config, ler_alvo, ler_conjunto

DADOS = Path(__file__).parent / "dados"


def test_conjunto_ok_csv_fasta():
    leitura_csv = ler_conjunto(DADOS / "conjunto_ok.csv")
    leitura_fasta = ler_conjunto(DADOS / "conjunto_ok.fasta")

    assert leitura_csv.primers == leitura_fasta.primers
    assert leitura_csv.primers["FIP"] == leitura_csv.primers["F1c"] + leitura_csv.primers["F2"]
    assert leitura_csv.conjunto is not None


def test_validacoes_csv(tmp_path):
    # CSV com nome flp vira chave LF
    f_flp = tmp_path / "flp.csv"
    f_flp.write_text("nome,sequencia\nflp,ACGTACGTACGTACGT\n", encoding="utf-8")
    l_flp = ler_conjunto(f_flp)
    assert "LF" in l_flp.primers
    assert l_flp.primers["LF"] == "ACGTACGTACGTACGT"

    # Nome XYZ -> ValueError
    f_xyz = tmp_path / "xyz.csv"
    f_xyz.write_text("nome,sequencia\nXYZ,ACGTACGTACGTACGT\n", encoding="utf-8")
    with pytest.raises(ValueError) as exc_xyz:
        ler_conjunto(f_xyz)
    assert "XYZ" in str(exc_xyz.value)

    # Nome repetido -> ValueError
    f_rep = tmp_path / "repetido.csv"
    f_rep.write_text("nome,sequencia\nF3,ACGTACGTACGTACGT\nF3,ACGTACGTACGTACGT\n", encoding="utf-8")
    with pytest.raises(ValueError):
        ler_conjunto(f_rep)

    # Sequência ACGU na terceira linha -> ValueError com "linha 3"
    f_acgu = tmp_path / "acgu.csv"
    f_acgu.write_text("nome,sequencia\nF3,ACGTACGTACGTACGT\nF2,ACGUACGTACGTACGT\n", encoding="utf-8")
    with pytest.raises(ValueError) as exc_acgu:
        ler_conjunto(f_acgu)
    assert "linha 3" in str(exc_acgu.value)

    # F1c sem F2 -> ValueError
    f_f1c = tmp_path / "f1c_sem_f2.csv"
    f_f1c.write_text("nome,sequencia\nF1c,ACGTACGTACGTACGT\n", encoding="utf-8")
    with pytest.raises(ValueError):
        ler_conjunto(f_f1c)

    # Partes + FIP diferente -> ValueError com "FIP"
    f_fip_dif = tmp_path / "fip_dif.csv"
    f_fip_dif.write_text(
        "nome,sequencia\nF1c,ACGTACGTACGTACGT\nF2,ACGTACGTACGTACGT\nFIP,TTTTTTTTTTTTTTTT\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError) as exc_fip:
        ler_conjunto(f_fip_dif)
    assert "FIP" in str(exc_fip.value)


def test_fip_bip_sem_partes_sem_alvo(tmp_path):
    f1c = "GTAAACCCCCAAGTGTTATTGG"
    f2 = "CGCATATGGATTCTGTAC"
    b1c = "AAGGGTAATTCTCCGTCCA"
    b2 = "GCAATAGTATAAAGTCTCGC"

    f_csv = tmp_path / "fip_bip_sem_alvo.csv"
    f_csv.write_text(
        f"nome,sequencia\n"
        f"F3,GTGACAGACTCTCCTTAGATA\n"
        f"B3,TGTTGCTGGATTGTCTGTT\n"
        f"FIP,{f1c + f2}\n"
        f"BIP,{b1c + b2}\n",
        encoding="utf-8",
    )
    leitura = ler_conjunto(f_csv)
    assert leitura.conjunto is None
    assert "FIP" in leitura.primers
    assert any("FIP sem F1c/F2" in a for a in leitura.avisos)


def test_fip_bip_sem_partes_com_alvo(tmp_path):
    f1c = "GTAAACCCCCAAGTGTTATTGG"
    f2 = "CGCATATGGATTCTGTAC"
    b1c = "AAGGGTAATTCTCCGTCCA"
    b2 = "GCAATAGTATAAAGTCTCGC"

    f_csv = tmp_path / "fip_bip_com_alvo.csv"
    f_csv.write_text(
        f"nome,sequencia\n"
        f"F3,GTGACAGACTCTCCTTAGATA\n"
        f"B3,TGTTGCTGGATTGTCTGTT\n"
        f"FIP,{f1c + f2}\n"
        f"BIP,{b1c + b2}\n",
        encoding="utf-8",
    )
    alvo = ler_alvo(DADOS / "alvo_ok.fasta")
    leitura = ler_conjunto(f_csv, alvo=alvo)

    assert leitura.conjunto is not None
    assert leitura.conjunto.F1c == "GTAAACCCCCAAGTGTTATTGG"
    assert leitura.conjunto.F2 == "CGCATATGGATTCTGTAC"
    assert leitura.conjunto.B1c == "AAGGGTAATTCTCCGTC"
    assert leitura.conjunto.B2 == "CAGCAATAGTATAAAGTCTCGC"
    assert any("ambíguo" in a and "FIP" in a for a in leitura.avisos)
    assert any("ambíguo" in a and "BIP" in a for a in leitura.avisos)


def test_ler_alvo_dois_registros(tmp_path):
    f_fasta = tmp_path / "multiplos.fasta"
    f_fasta.write_text(">reg1\nACGT\n>reg2\nTGCA\n", encoding="utf-8")
    with pytest.raises(ValueError):
        ler_alvo(f_fasta)
