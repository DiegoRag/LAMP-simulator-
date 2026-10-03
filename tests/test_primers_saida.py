"""Testes para o módulo de saída e CLI de primers LAMP."""

import csv
import json
from pathlib import Path
from primers.analise import analisar
from primers.cli import main
from primers.entrada import carregar_config, ler_alvo, ler_conjunto
from primers.saida import exportar

DADOS = Path(__file__).parent / "dados"
CONFIG = Path(__file__).parents[1] / "config" / "primers.toml"


def test_exportar_analise_ok(tmp_path):
    cfg = carregar_config(CONFIG)
    leitura_ok = ler_conjunto(DADOS / "conjunto_ok.csv")
    alvo_ok = ler_alvo(DADOS / "alvo_ok.fasta")
    res = analisar(leitura_ok, alvo_ok, cfg)

    pasta = tmp_path / "a"
    exportar(res, pasta, cfg)

    assert (pasta / "primers.csv").exists()
    assert (pasta / "heterodimeros.csv").exists()
    assert (pasta / "violacoes.csv").exists()
    assert (pasta / "resultado.json").exists()

    with open(pasta / "primers.csv", "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        chaves_esperadas = [
            "nome",
            "sequencia",
            "comprimento",
            "gc_pct",
            "tm_c",
            "dg_extremidade_kcal",
            "extremidade",
            "grampo_tm_c",
            "grampo_dg_kcal",
            "homodimero_dg_kcal",
            "avisos",
        ]
        assert reader.fieldnames == chaves_esperadas
        linhas = list(reader)
        nomes = [l["nome"] for l in linhas]
        assert nomes == ["F3", "B3", "F2", "F1c", "B2", "B1c", "FIP", "BIP"]
        linha_fip = next(l for l in linhas if l["nome"] == "FIP")
        assert linha_fip["tm_c"] == ""

    with open(pasta / "violacoes.csv", "r", encoding="utf-8") as f:
        reader_viol = csv.reader(f)
        linhas_viol = list(reader_viol)
        assert len(linhas_viol) == 1
        assert linhas_viol[0] == ["item", "regra", "valor", "faixa", "mensagem"]

    with open(pasta / "resultado.json", "r", encoding="utf-8") as f:
        dados_json = json.load(f)

    assert set(dados_json.keys()) == {
        "versao_codigo",
        "condicoes",
        "config",
        "avisos",
        "primers",
        "heterodimeros",
        "violacoes",
        "regioes",
    }
    assert dados_json["regioes"]["F2"] == [31, 49]


def test_heterodimeros_csv(tmp_path):
    cfg = carregar_config(CONFIG)
    leitura_ok = ler_conjunto(DADOS / "conjunto_ok.csv")
    alvo_ok = ler_alvo(DADOS / "alvo_ok.fasta")
    res = analisar(leitura_ok, alvo_ok, cfg)

    pasta = tmp_path / "het"
    exportar(res, pasta, cfg)

    with open(pasta / "heterodimeros.csv", "r", encoding="utf-8") as f:
        linhas = [l.strip() for l in f.readlines() if l.strip()]

    assert linhas[0] == "nome,F3,B3,FIP,BIP"
    assert len(linhas) == 5


def test_exportar_deterministico(tmp_path):
    cfg = carregar_config(CONFIG)
    leitura_ok = ler_conjunto(DADOS / "conjunto_ok.csv")
    alvo_ok = ler_alvo(DADOS / "alvo_ok.fasta")
    res = analisar(leitura_ok, alvo_ok, cfg)

    pasta1 = tmp_path / "d1"
    pasta2 = tmp_path / "d2"
    exportar(res, pasta1, cfg)
    exportar(res, pasta2, cfg)

    for nome in ("primers.csv", "heterodimeros.csv", "violacoes.csv", "resultado.json"):
        assert (pasta1 / nome).read_bytes() == (pasta2 / nome).read_bytes()


def test_cli_analisar(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    entrada = (DADOS / "conjunto_ok.csv").resolve()
    cfg_path = CONFIG.resolve()

    codigo_sucesso = main(
        [
            "analisar",
            "--entrada",
            str(entrada),
            "--config",
            str(cfg_path),
            "--saida",
            "saida/x",
        ]
    )
    assert codigo_sucesso == 0
    assert (tmp_path / "saida" / "x" / "primers.csv").exists()

    codigo_erro = main(
        [
            "analisar",
            "--entrada",
            str(entrada),
            "--config",
            str(cfg_path),
            "--saida",
            "fora",
        ]
    )
    assert codigo_erro == 2
