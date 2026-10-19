"""Fase 6, sessão 2: resultado.json completo (alvo, fita, corte do FIP/BIP, erros, incompat)."""

import hashlib
import json
from pathlib import Path

from primers.analise import analisar
from primers.cli import main
from primers.desenho import desenhar
from primers.entrada import carregar_config, ler_alvo, ler_conjunto
from primers.saida import (
    COORDENADAS,
    COORDENADAS_DESENHO,
    CSV_COORDENADAS,
    descrever_alvo,
    descrever_conjunto,
    exportar,
    exportar_desenho,
    versao_do_codigo,
)

RAIZ = Path(__file__).resolve().parents[1]
DADOS = RAIZ / "tests" / "dados"
CONFIG = RAIZ / "config" / "primers.toml"
CFG = carregar_config(CONFIG)
ALVO = ler_alvo(DADOS / "alvo_ok.fasta")
ZEROS = [
    {"regiao": r, "encontrada": True, "total": 0, "ponta": 0} for r in ("F3", "F2", "F1", "B1c", "B2c", "B3c")
]


def _json(pasta: Path) -> dict:
    return json.loads((pasta / "resultado.json").read_text(encoding="utf-8"))


def test_resultado_json_reproduz_a_analise(tmp_path):
    leitura = ler_conjunto(DADOS / "conjunto_ok.csv", alvo=ALVO)
    res = analisar(leitura, ALVO, CFG)
    alvo_json = descrever_alvo(ALVO, id_registro="alvo_sintetico_ok", arquivo="alvo_ok.fasta")
    eficiencia = {"f_previsto": 1.0, "modelo": "teste", "parametros": {"fator_ponta3": 0.3}}
    exportar(res, tmp_path, CFG, leitura=leitura, alvo=alvo_json, eficiencia=eficiencia)
    dados = _json(tmp_path)

    assert dados["coordenadas"] == COORDENADAS
    assert dados["alvo"] == {
        "id": "alvo_sintetico_ok",
        "arquivo": "alvo_ok.fasta",
        "comprimento": 216,
        "sha256": hashlib.sha256(ALVO.encode("utf-8")).hexdigest(),
        "circular": False,
    }
    assert ALVO not in (tmp_path / "resultado.json").read_text(encoding="utf-8")  # nunca a sequência
    assert dados["fita"] == "+"
    assert dados["conjunto"]["FIP"] == {
        "sequencia": "GTAAACCCCCAAGTGTTATTGGCGCATATGGATTCTGTAC",
        "origem": "partes informadas",
        "F1c": [1, 22],
        "ligante": None,
        "F2": [23, 40],
    }
    assert dados["conjunto"]["BIP"]["B1c"] == [1, 19] and dados["conjunto"]["BIP"]["B2"] == [20, 39]
    assert dados["conjunto"]["ligante"] == "" and dados["conjunto"]["LF"] is None
    assert dados["regioes_aproximadas"] == {}
    assert dados["incompat"] == {
        "total": 0,
        "ponta": 0,
        "bases_ponta": 5,
        "nao_encontradas": [],
        "por_regiao": ZEROS,
    }
    assert dados["eficiencia_prevista"] == eficiencia


def test_fip_inteiro_registra_o_corte(tmp_path):
    c = ler_conjunto(DADOS / "conjunto_ok.csv").conjunto
    csv = f"nome,sequencia\nF3,{c.F3}\nB3,{c.B3}\nFIP,{c.F1c}TTTT{c.F2}\nBIP,{c.B1c}TTTT{c.B2}\n"
    (tmp_path / "inteiros.csv").write_text(csv, encoding="utf-8")
    leitura = ler_conjunto(tmp_path / "inteiros.csv", alvo=ALVO)
    assert leitura.separados == ("FIP", "BIP")
    conjunto = descrever_conjunto(leitura)
    assert conjunto is not None
    assert conjunto["ligante"] == "TTTT"
    assert conjunto["FIP"] == {
        "sequencia": "GTAAACCCCCAAGTGTTATTGGTTTTCGCATATGGATTCTGTAC",
        "origem": "FIP inteiro separado nas partes",
        "F1c": [1, 22],
        "ligante": [23, 26],
        "F2": [27, 44],
    }
    assert conjunto["BIP"] == {
        "sequencia": "AAGGGTAATTCTCCGTCCATTTTGCAATAGTATAAAGTCTCGC",
        "origem": "BIP inteiro separado nas partes",
        "B1c": [1, 19],
        "ligante": [20, 23],
        "B2": [24, 43],
    }
    assert ler_conjunto(DADOS / "conjunto_ok.csv", alvo=ALVO).separados == ()


def test_erros_de_base_e_incompat_no_json(tmp_path):
    assert ALVO[48] == "C"  # última base de F2
    alvo_mut = ALVO[:48] + "A" + ALVO[49:]
    leitura = ler_conjunto(DADOS / "conjunto_ok.csv", alvo=alvo_mut)
    exportar(analisar(leitura, alvo_mut, CFG), tmp_path, CFG, leitura=leitura)
    dados = _json(tmp_path)
    assert dados["regioes_aproximadas"] == {
        "F2": {
            "regiao": [31, 49],
            "erros": [
                {
                    "regiao": "F2",
                    "primer": "F2",
                    "fita": "+",
                    "posicao": 49,
                    "base_no_primer": 18,
                    "comprimento": 18,
                    "primer_tem": "C",
                    "alvo_pede": "A",
                    "ponta": "3'",
                    "distancia_ponta": 0,
                }
            ],
        }
    }
    assert dados["incompat"]["total"] == 1 and dados["incompat"]["ponta"] == 1
    assert dados["incompat"]["por_regiao"][1] == {"regiao": "F2", "encontrada": True, "total": 1, "ponta": 1}
    assert dados["alvo"] is None and dados["eficiencia_prevista"] is None


def test_sem_alvo_os_campos_ficam_nulos(tmp_path):
    leitura = ler_conjunto(DADOS / "conjunto_ok.csv")
    exportar(analisar(leitura, None, CFG), tmp_path, CFG)
    dados = _json(tmp_path)
    nulos = ("alvo", "fita", "conjunto", "regioes", "regioes_aproximadas", "incompat", "eficiencia_prevista")
    for chave in nulos:
        assert dados[chave] is None, chave


def test_cli_grava_o_alvo_e_o_corte(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    codigo = main(
        [
            "analisar",
            "--entrada",
            str(DADOS / "conjunto_ok.csv"),
            "--alvo",
            str(DADOS / "alvo_ok.fasta"),
            "--config",
            str(CONFIG),
            "--saida",
            "saida/c",
        ]
    )
    assert codigo == 0
    dados = _json(tmp_path / "saida" / "c")
    assert dados["alvo"]["id"] == "alvo_sintetico_ok"
    assert dados["alvo"]["arquivo"] == "alvo_ok.fasta"
    assert dados["conjunto"]["FIP"]["origem"] == "partes informadas"
    assert dados["incompat"]["total"] == 0
    assert dados["eficiencia_prevista"] is None  # a CLI do primers não conhece a esteira


def test_desenho_rotula_as_coordenadas(tmp_path):
    exportar_desenho(desenhar(ALVO, CFG), ALVO, tmp_path, CFG)
    linhas = (tmp_path / "conjuntos.csv").read_text(encoding="utf-8").splitlines()
    assert linhas[0].startswith("# Candidatos computacionais")
    assert linhas[1] == CSV_COORDENADAS
    assert linhas[2].startswith("posicao,pontuacao")
    assert _json(tmp_path)["coordenadas"] == COORDENADAS_DESENHO
    assert isinstance(versao_do_codigo(), str) and versao_do_codigo()
