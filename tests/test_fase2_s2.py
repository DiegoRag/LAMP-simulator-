"""Testes da Fase 2 (Sessão 2): entrada aberta, multi-FASTA, CSV flexível,
avaliação de quimeras e análise segura."""

from pathlib import Path

import pytest

from primers.analise import analisar, analisar_registros
from primers.entrada import carregar_config, ler_alvo, ler_alvos, ler_conjunto

DADOS = Path(__file__).parent / "dados"
CONFIG = Path(__file__).resolve().parents[1] / "config" / "primers.toml"


def test_csv_excel_brasileiro_ponto_e_virgula_aspas_name_sequence_e_5p_3p(tmp_path):
    f_csv = tmp_path / "excel_br.csv"
    conteudo = (
        '\ufeff"name";"sequence"\n'
        '"fosA-F3";"5\'-GTGACAGACTCTCCTTAGATA-3\'"\n'
        '"B3_v2";"5′ – TGTTGCTGGATTGTCTGTT – 3′"\n'
        '"alvo_F2";"CGCATATGGATTCTGTAC"\n'
        '"fosA.F1c";"GTAAACCCCCAAGTGTTATTGG"\n'
        '"B2";"GCAATAGTATAAAGTCTCGC"\n'
        '"B1c_1";"AAGGGTAATTCTCCGTCCA"\n'
        '"Loop-F";"CGTATCCAGCGACAGGCA"\n'
    )
    f_csv.write_text(conteudo, encoding="utf-8")

    leitura = ler_conjunto(f_csv)
    assert leitura.conjunto is not None
    assert leitura.conjunto.F3 == "GTGACAGACTCTCCTTAGATA"
    assert leitura.conjunto.B3 == "TGTTGCTGGATTGTCTGTT"
    assert leitura.conjunto.LF == "CGTATCCAGCGACAGGCA"


def test_fasta_com_bom_sem_cabecalho_e_vazio(tmp_path):
    # FASTA de alvo com BOM e U convertido para T
    f_bom = tmp_path / "alvo_bom.fasta"
    f_bom.write_text("\ufeff>alvo1\nAUGCACGT\n", encoding="utf-8")
    assert ler_alvo(f_bom) == "ATGCACGT"

    # FASTA de alvo sem '>'
    f_sem_maior = tmp_path / "sem_cabecalho.fasta"
    f_sem_maior.write_text("ACGTACGT\n", encoding="utf-8")
    with pytest.raises(ValueError) as exc_sem:
        ler_alvo(f_sem_maior)
    assert "falta cabeçalho '>'" in str(exc_sem.value)

    # FASTA de primers vazio
    f_vazio = tmp_path / "vazio.fasta"
    f_vazio.write_text("   \n\n", encoding="utf-8")
    with pytest.raises(ValueError) as exc_vazio:
        ler_conjunto(f_vazio)
    assert "vazio" in str(exc_vazio.value).lower()

    # FASTA de primers com registro de sequência vazia
    f_reg_vazio = tmp_path / "reg_vazio.fasta"
    f_reg_vazio.write_text(">F3\n\n", encoding="utf-8")
    with pytest.raises(ValueError) as exc_reg:
        ler_conjunto(f_reg_vazio)
    assert "Sequência vazia" in str(exc_reg.value)


def test_fasta_multi_registro_escolha_e_analise_painel(tmp_path):
    cfg = carregar_config(CONFIG)
    seq_ok = ler_alvo(DADOS / "alvo_ok.fasta")

    # Garante que a posição 48 seja alterada ('C' para 'A')
    assert seq_ok[48] == "C"
    seq_snp = seq_ok[:48] + "A" + seq_ok[49:]
    seq_incompativel = "A" * 220

    f_multi = tmp_path / "multi.fasta"
    conteudo = (
        f">var_incompativel\n{seq_incompativel}\n"
        f">var_correta\n{seq_ok}\n"
        f">var_snp\n{seq_snp}\n"
    )
    f_multi.write_text(conteudo, encoding="utf-8")

    regs = ler_alvos(f_multi)
    assert len(regs) == 3

    leitura = ler_conjunto(DADOS / "conjunto_ok.csv")
    # Escolha automática baseada no conjunto de primers
    assert ler_alvo(f_multi, conjunto=leitura.conjunto) == regs[1].sequencia
    # Escolha explícita por ID
    assert ler_alvo(f_multi, registro="var_snp") == regs[2].sequencia

    painel = analisar_registros(DADOS / "conjunto_ok.csv", regs, cfg)
    # Tm e ponta de conjunto_ok reprovam na convenção da fase 3; aqui só a geometria importa
    def geometricas(res):
        return [v for v in res.violacoes if v.regra not in ("tm_c", "extremidade")]

    assert geometricas(painel["var_correta"]) == []
    assert len(geometricas(painel["var_snp"])) == 1
    assert geometricas(painel["var_snp"])[0].item == "F2"
    assert geometricas(painel["var_snp"])[0].regra == "nao_encontrado"
    assert len(geometricas(painel["var_incompativel"])) == 6


def test_fip_bip_inteiros_sem_alvo_avaliados_por_cortes_plausiveis(tmp_path):
    cfg = carregar_config(CONFIG)

    # Caso 1: 4 primers sem alvo e sem violações
    f_ok4 = tmp_path / "ok4.csv"
    f_ok4.write_text(
        "nome,sequencia\n"
        "F3,GTGACAGACTCTCCTTAGATA\n"
        "B3,TGTTGCTGGATTGTCTGTT\n"
        "FIP,GTAAACCCCCAAGTGTTATTGGCGCATATGGATTCTGTAC\n"
        "BIP,AAGGGTAATTCTCCGTCCAGCAATAGTATAAAGTCTCGC\n",
        encoding="utf-8",
    )
    leitura_ok = ler_conjunto(f_ok4)
    res_ok = analisar(leitura_ok, None, cfg)
    # O corte de FIP com menos violações agora é F1c = GTAAACCCCCAAGTGTTAT + F2 = TGGCGCATATGGATTCTGTAC
    assert {
        (v.item, v.regra) for v in res_ok.violacoes
    } == {
        ("F3", "tm_c"),
        ("F3", "extremidade"),
        ("B3", "tm_c"),
        ("F1c", "tm_c"),
        ("F1c", "extremidade"),
        ("B2", "tm_c"),
        ("B1c", "tm_c"),
    }

    # Caso 2: 4 primers sem alvo com violações em F3 e B1c (dentro do BIP)
    f_viol4 = tmp_path / "viol4.csv"
    f_viol4.write_text(
        "nome,sequencia\n"
        "F3,ATTATAAATTTAGATTAT\n"
        "B3,TGTTGCTGGATTGTCTGTT\n"
        "FIP,GTAAACCCCCAAGTGTTATTGGCGCATATGGATTCTGTAC\n"
        "BIP,TATATAAATTCTCCGTCCAGCAATAGTATAAAGTCTCGC\n",
        encoding="utf-8",
    )
    leitura_viol = ler_conjunto(f_viol4)
    res_viol = analisar(leitura_viol, None, cfg)
    assert {
        (v.item, v.regra) for v in res_viol.violacoes
    } == {
        ("F3", "gc_pct"),
        ("F3", "tm_c"),
        ("F3", "extremidade"),
        ("B3", "tm_c"),
        ("F1c", "tm_c"),
        ("F1c", "extremidade"),
        ("B2", "tm_c"),
        ("B1c", "gc_pct"),
        ("B1c", "tm_c"),
        ("B1c", "extremidade"),
    }


def test_fip_acima_de_60_nt_gera_aviso_sem_runtime_error(tmp_path):
    cfg = carregar_config(CONFIG)

    f_longo = tmp_path / "longo.csv"
    fip_longo = (
        "GTAAACCCCCAAGTGTTATTGG" + "TTTT" + "CGCATATGGATTCTGTAC" + "ACGTACGTACGTACGTACGT"
    )
    assert len(fip_longo) == 64

    f_longo.write_text(
        f"nome,sequencia\n"
        f"F3,GTGACAGACTCTCCTTAGATA\n"
        f"B3,TGTTGCTGGATTGTCTGTT\n"
        f"FIP,{fip_longo}\n"
        f"BIP,AAGGGTAATTCTCCGTCCAGCAATAGTATAAAGTCTCGC\n",
        encoding="utf-8",
    )
    leitura = ler_conjunto(f_longo)
    res = analisar(leitura, None, cfg)

    # Verifica linha do FIP
    linha_fip = next(linha for linha in res.linhas if linha["nome"] == "FIP")
    assert linha_fip["grampo_tm_c"] is None
    assert linha_fip["grampo_dg_kcal"] is None
    assert linha_fip["homodimero_dg_kcal"] is None
    assert "acima de 60 nt" in linha_fip["avisos"]

    # Verifica avisos globais e violações
    assert any("acima do limite de 60 nt do Primer3" in a for a in res.avisos)
    assert any(v.item == "FIP" and v.regra == "comprimento" for v in res.violacoes)
