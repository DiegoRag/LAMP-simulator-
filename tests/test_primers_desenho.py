"""Testes para o módulo de desenho de primers LAMP."""

import copy
import json
from pathlib import Path

import pytest

from primers.cli import main
from primers.desenho import (
    AVISO,
    AVISO_ASCII,
    desenhar,
    gerar_candidatos,
)
from primers.entrada import carregar_config, ler_alvo
from primers.modelo import Regiao
from primers.restricoes import Violacao, calcular_distancias, verificar
from primers.saida import exportar_desenho

DADOS = Path(__file__).parent / "dados"
CONFIG = Path(__file__).parents[1] / "config" / "primers.toml"


def test_gerar_candidatos_alvo_ok():
    cfg = carregar_config(CONFIG)
    alvo = ler_alvo(DADOS / "alvo_ok.fasta")
    cands = gerar_candidatos(alvo, cfg)

    # fase 4: teto por janela (2 a cada 10 nt) no lugar do teto global; antes 157/157/86/87/153/153
    assert len(cands["F3"]) == 35
    assert len(cands["F2"]) == 35
    assert len(cands["F1"]) == 22
    assert len(cands["B1c"]) == 20
    assert len(cands["B2c"]) == 32
    assert len(cands["B3c"]) == 32

    regioes_esperadas = {
        "F3": Regiao(1, 26),
        "F2": Regiao(26, 49),
        "F1": Regiao(75, 99),
        "B1c": Regiao(105, 128),
        "B2c": Regiao(145, 169),
        "B3c": Regiao(194, 213),
    }

    for tipo, regiao in regioes_esperadas.items():
        regioes_encontradas = [c.regiao for c in cands[tipo]]
        assert regiao in regioes_encontradas, f"Região {regiao} não encontrada em cands[{tipo}]"


def test_desenhar_alvo_ok():
    cfg = carregar_config(CONFIG)
    alvo = ler_alvo(DADOS / "alvo_ok.fasta")
    res = desenhar(alvo, cfg)

    assert len(res.conjuntos) >= 1
    # fase 4: a busca percorre o alvo inteiro (não há mais busca_truncada)
    assert res.estatisticas["diagnostico"] == ""
    assert res.estatisticas["nucleos_distintos"] == 4

    posicoes = [c.posicao for c in res.conjuntos]
    assert posicoes == list(range(1, len(res.conjuntos) + 1))

    pontuacoes = [c.pontuacao for c in res.conjuntos]
    for i in range(len(pontuacoes) - 1):
        assert pontuacoes[i] <= pontuacoes[i + 1]

    cfg_copia = copy.deepcopy(cfg)
    cfg_copia["desenho"]["peso_dimero"] = 0
    res_copia = desenhar(alvo, cfg_copia)

    assert len(res_copia.conjuntos) >= 1
    assert res_copia.conjuntos[0].desvio_tm == pytest.approx(0.2762, abs=1e-3)


def test_autoverificacao_conjuntos():
    cfg = carregar_config(CONFIG)
    alvo = ler_alvo(DADOS / "alvo_ok.fasta")
    res = desenhar(alvo, cfg)

    assert len(res.conjuntos) >= 1
    for c in res.conjuntos:
        assert verificar(c.conjunto, alvo, cfg) == []


def test_determinismo(tmp_path):
    cfg = carregar_config(CONFIG)
    alvo = ler_alvo(DADOS / "alvo_ok.fasta")

    res1 = desenhar(alvo, cfg)
    res2 = desenhar(alvo, cfg)

    pasta1 = tmp_path / "d1"
    pasta2 = tmp_path / "d2"

    exportar_desenho(res1, alvo, pasta1, cfg)
    exportar_desenho(res2, alvo, pasta2, cfg)

    assert (pasta1 / "conjuntos.csv").read_bytes() == (pasta2 / "conjuntos.csv").read_bytes()
    assert (pasta1 / "resultado.json").read_bytes() == (pasta2 / "resultado.json").read_bytes()


def test_aviso_e_privacidade(tmp_path):
    assert "ADR-0013" in AVISO
    assert "ADR-0013" in AVISO_ASCII
    assert AVISO_ASCII.isascii()

    cfg = carregar_config(CONFIG)
    alvo = ler_alvo(DADOS / "alvo_ok.fasta")
    res = desenhar(alvo, cfg)

    pasta = tmp_path / "saida_aviso"
    exportar_desenho(res, alvo, pasta, cfg)

    primeira_linha_csv = (pasta / "conjuntos.csv").read_text(encoding="utf-8").splitlines()[0]
    assert primeira_linha_csv.startswith("# Candidatos computacionais")

    conteudo_json_str = (pasta / "resultado.json").read_text(encoding="utf-8")
    dados_json = json.loads(conteudo_json_str)

    assert "ADR-0013" in dados_json["aviso"]
    assert alvo not in conteudo_json_str


def test_alvo_sem_solucao():
    cfg = carregar_config(CONFIG)
    alvo = "A" * 300
    res = desenhar(alvo, cfg)

    assert res.conjuntos == []
    assert all(qtd == 0 for qtd in res.estatisticas["candidatos_por_tipo"].values())


def test_guarda_autoverificacao(monkeypatch):
    cfg = carregar_config(CONFIG)
    alvo = ler_alvo(DADOS / "alvo_ok.fasta")

    def mock_verificar(conjunto, alvo, cfg):
        return [
            Violacao(
                item="F3",
                regra="teste",
                valor=0.0,
                faixa="[0, 1]",
                mensagem="Violação simulada de teste",
            )
        ]

    monkeypatch.setattr("primers.desenho.verificar", mock_verificar)

    with pytest.raises(RuntimeError) as exc_info:
        desenhar(alvo, cfg)

    assert "bug no desenho" in str(exc_info.value)


def test_calcular_distancias():
    regioes = {
        "F3": Regiao(0, 21),
        "F2": Regiao(31, 49),
        "F1": Regiao(81, 103),
        "B1c": Regiao(108, 127),
        "B2c": Regiao(157, 177),
        "B3c": Regiao(187, 206),
    }
    dist = calcular_distancias(regioes)
    esperado = {
        "F2_B2": 146,
        "F2_F1": 50,
        "B1c_B2c": 50,
        "F3_F2_gap": 10,
        "B2c_B3c_gap": 10,
        "F1_B1c_gap": 5,
        "F3_B3_total": 206,
    }
    assert dist == esperado


def test_cli_desenhar(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    alvo_path = (DADOS / "alvo_ok.fasta").resolve()
    cfg_path = CONFIG.resolve()

    codigo_sucesso = main(
        [
            "desenhar",
            "--alvo",
            str(alvo_path),
            "--config",
            str(cfg_path),
            "--saida",
            "saida/d",
        ]
    )
    assert codigo_sucesso == 0
    assert (tmp_path / "saida" / "d" / "conjuntos.csv").exists()
    assert (tmp_path / "saida" / "d" / "resultado.json").exists()

    codigo_erro = main(
        [
            "desenhar",
            "--alvo",
            str(alvo_path),
            "--config",
            str(cfg_path),
            "--saida",
            "fora",
        ]
    )
    assert codigo_erro == 2


def test_escolha_f3_b3c_otimos_por_nucleo():
    cfg = carregar_config(CONFIG)
    alvo = ler_alvo(DADOS / "alvo_ok.fasta")
    res = desenhar(alvo, cfg)
    cands = gerar_candidatos(alvo, cfg)
    D = cfg["distancias"]

    assert len(res.conjuntos) >= 1
    for c in res.conjuntos:
        reg_f2 = c.regioes["F2"]
        reg_f3 = c.regioes["F3"]
        reg_b2c = c.regioes["B2c"]
        reg_b3c = c.regioes["B3c"]

        validos_f3 = [
            cand
            for cand in cands["F3"]
            if D["F3_F2_gap"][0] <= reg_f2.ini - cand.regiao.fim <= D["F3_F2_gap"][1]
        ]
        assert len(validos_f3) > 0
        melhor_f3 = min(validos_f3, key=lambda x: (x.desvio_tm, x.regiao.ini, x.regiao.fim))
        assert reg_f3 == melhor_f3.regiao
        assert c.conjunto.F3 == melhor_f3.primer

        validos_b3c = [
            cand
            for cand in cands["B3c"]
            if D["B2c_B3c_gap"][0] <= cand.regiao.ini - reg_b2c.fim <= D["B2c_B3c_gap"][1]
        ]
        assert len(validos_b3c) > 0
        melhor_b3c = min(validos_b3c, key=lambda x: (x.desvio_tm, x.regiao.ini, x.regiao.fim))
        assert reg_b3c == melhor_b3c.regiao
        assert c.conjunto.B3 == melhor_b3c.primer

