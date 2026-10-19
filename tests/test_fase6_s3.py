"""Fase 6, sessão 3: contrato selecao@2 com incompat_regioes, e a esteira fechada (escreve e lê)."""

import copy
import hashlib
import json
from pathlib import Path

import pytest

from esteira.contrato import ESQUEMA_V2, REGIOES_INCOMPAT, ler_selecao, validar_selecao
from esteira.ponte import (
    acesso_do_id,
    assinatura,
    carregar_config,
    eixo_minutos,
    medida_da_selecao,
    selecao_da_analise,
    simular_tubos,
    tubos_da_selecao,
)
from primers.analise import analisar
from primers.entrada import carregar_config as carregar_config_primers
from primers.entrada import ler_alvo, ler_conjunto
from primers.incompat import Incompat, IncompatRegiao

RAIZ = Path(__file__).resolve().parents[1]
DADOS = RAIZ / "tests" / "dados"
CFG = carregar_config(RAIZ / "config" / "esteira.toml")
CFG_PRIMERS = carregar_config_primers(RAIZ / "config" / "primers.toml")
GERADO = {"versao_codigo": "teste", "experimento": "fase6", "config_sha256": "0"}
PRIMERS = {"F3": "ACGT", "B3": "ACGT", "FIP": "ACGT", "BIP": "ACGT"}
ESSENCIAIS = ("F3", "F2", "F1", "B1c", "B2c", "B3c")
MEDIDA = "conjuntos[0].por_alelo[0]"
SERIE = [
    {"nome": "10⁵ cópias", "copias": 1e5, "papel": "amostra"},
    {"nome": "Controle negativo", "copias": 0, "papel": "controle_negativo"},
]


def _inc(f1_total: int, f1_ponta: int) -> Incompat:
    itens = tuple(
        IncompatRegiao(r, True, f1_total, f1_ponta) if r == "F1" else IncompatRegiao(r, True, 0, 0)
        for r in ESSENCIAIS
    )
    return Incompat(itens, f1_total, f1_ponta, ())


def _selecao(incompat, **extra) -> dict:
    argumentos = {
        "conjunto_id": "sintetico",
        "origem": "usuario",
        "primers": PRIMERS,
        "incompat": incompat,
        "alelo_id": "alvo_sintetico",
        "acesso": None,
        "cfg": CFG,
        "gerado_por": GERADO,
    } | extra
    return selecao_da_analise(**argumentos)


def test_selecao_da_analise_de_conjunto_ok():
    leitura = ler_conjunto(DADOS / "conjunto_ok.csv")
    alvo = ler_alvo(DADOS / "alvo_ok.fasta")
    res = analisar(leitura, alvo, CFG_PRIMERS)
    dados = _selecao(res.incompat, alelo_id="alvo_sintetico_ok")
    assert dados["esquema"] == ESQUEMA_V2 == "lamp-sim/selecao@2"
    assert validar_selecao(dados) == []
    conjunto = dados["conjuntos"][0]
    assert conjunto["selecionado"] is True
    assert conjunto["motivo"] == (
        "eficiência prevista 1 (modelo provisório de config/esteira.toml); "
        "0 incompatibilidades, 0 na ponta que estende"
    )
    assert conjunto["por_alelo"] == [
        {
            "alelo_id": "alvo_sintetico_ok",
            "acesso": None,
            "incompat_total": 0,
            "incompat_3p": 0,
            "f_previsto": 1.0,
            "incompat_regioes": [
                {"regiao": r, "encontrada": True, "total": 0, "ponta": 0} for r in ESSENCIAIS
            ],
        }
    ]


def test_erro_no_meio_e_na_ponta_5_de_f1c_ficam_distintos_no_contrato():
    meio = _selecao(_inc(1, 0))["conjuntos"][0]["por_alelo"][0]
    ponta = _selecao(_inc(1, 1))["conjuntos"][0]["por_alelo"][0]
    assert meio["incompat_regioes"][2] == {"regiao": "F1", "encontrada": True, "total": 1, "ponta": 0}
    assert ponta["incompat_regioes"][2] == {"regiao": "F1", "encontrada": True, "total": 1, "ponta": 1}
    assert meio["f_previsto"] == pytest.approx(0.85)
    assert ponta["f_previsto"] == pytest.approx(0.3)


def test_a_esteira_fecha_selecao_escrita_e_lida(tmp_path):
    dados = _selecao(_inc(1, 1))
    caminho = tmp_path / "selecao.json"
    caminho.write_text(json.dumps(dados, ensure_ascii=False, indent=2), encoding="utf-8")
    lidos = ler_selecao(caminho)
    t = eixo_minutos(30.0, 0.5)
    tubos = tubos_da_selecao(lidos, SERIE, t, 2026)
    assert tubos == simular_tubos(0.3, SERIE, t, 2026)
    assert medida_da_selecao(lidos, "sintetico", "alvo_sintetico")["f_previsto"] == pytest.approx(0.3)


def test_leitura_sem_medida_ou_com_nome_errado():
    sem_medida = _selecao(None)
    assert sem_medida["conjuntos"][0]["selecionado"] is False
    assert sem_medida["conjuntos"][0]["motivo"] == "sem medida contra o alvo: nada a prever"
    assert "incompat_regioes" not in sem_medida["conjuntos"][0]["por_alelo"][0]
    with pytest.raises(ValueError, match="sem f_previsto"):
        tubos_da_selecao(sem_medida, SERIE, eixo_minutos(10.0, 0.5), 1)
    with pytest.raises(KeyError):
        medida_da_selecao(sem_medida, conjunto_id="outro")
    with pytest.raises(KeyError):
        medida_da_selecao(sem_medida, alelo_id="outro")


def test_regiao_essencial_ausente():
    itens = [
        IncompatRegiao(r, False, None, None) if r == "F2" else IncompatRegiao(r, True, 0, 0)
        for r in ESSENCIAIS
    ]
    sem_f2 = _selecao(Incompat(tuple(itens), None, None, ("F2",)))
    medida = sem_f2["conjuntos"][0]["por_alelo"][0]
    assert (medida["incompat_total"], medida["incompat_3p"], medida["f_previsto"]) == (None, None, 0.0)
    assert sem_f2["conjuntos"][0]["selecionado"] is False
    assert sem_f2["conjuntos"][0]["motivo"].endswith("; sem posição no alvo: F2")


def test_validacao_de_incompat_regioes():
    base = _selecao(_inc(1, 1))
    assert REGIOES_INCOMPAT == ("F3", "F2", "F1", "B1c", "B2c", "B3c", "LF", "LB")

    def com(mudar) -> list[str]:
        dados = copy.deepcopy(base)
        mudar(dados["conjuntos"][0]["por_alelo"][0])
        return validar_selecao(dados)

    def soma_errada(m):
        m["incompat_total"] = 2

    def regiao_estranha(m):
        m["incompat_regioes"][0]["regiao"] = "F4"

    def repetida(m):
        m["incompat_regioes"][1]["regiao"] = "F3"

    def ponta_demais(m):
        m["incompat_regioes"][0]["ponta"] = 1

    def ausente_com_numero(m):
        m["incompat_regioes"][0]["encontrada"] = False

    def ausente_com_total(m):
        m["incompat_regioes"][0].update(encontrada=False, total=None, ponta=None)

    onde = f"{MEDIDA}.incompat_regioes"
    assert com(soma_errada) == [f"{MEDIDA}: incompat_total (2) difere da soma das regiões (1)"]
    assert com(regiao_estranha) == [f"{onde}[0].regiao 'F4' deve ser uma de {list(REGIOES_INCOMPAT)}"]
    assert com(repetida) == [f"{onde}[1]: região 'F3' repetida"]
    assert com(ponta_demais) == [f"{onde}[0]: ponta (1) maior que total (0)"]
    assert com(ausente_com_numero) == [f"{onde}[0]: região não encontrada leva total e ponta null"]
    assert com(ausente_com_total) == [
        f"{MEDIDA}: com região não encontrada, incompat_total e incompat_3p devem ser null"
    ]


def test_v1_continua_congelada():
    dados = _selecao(_inc(0, 0))
    v1 = copy.deepcopy(dados)
    v1["esquema"] = "lamp-sim/selecao@1"
    v1["conjuntos"][0]["conjunto_id"] = "sint.1"
    problemas = validar_selecao(v1)
    assert f"{MEDIDA}: campo desconhecido 'incompat_regioes'" in problemas
    assert "conjuntos[0].origem 'usuario' deve ser 'artigo' ou 'desenho'" in problemas
    dados["conjuntos"][0]["origem"] = "bancada"
    assert validar_selecao(dados) == [
        "conjuntos[0].origem 'bancada' deve ser 'artigo', 'desenho' ou 'usuario'"
    ]


def test_selecao_invalida_nao_e_gravada(tmp_path):
    with pytest.raises(ValueError, match="falta o campo 'BIP'"):
        _selecao(_inc(0, 0), primers={"F3": "ACGT", "B3": "ACGT", "FIP": "ACGT"})
    with pytest.raises(ValueError, match="sem aviso"):
        _selecao(_inc(0, 0), origem="desenho")
    dados = _selecao(_inc(0, 0))
    dados["conjuntos"][0]["selecionado"] = "sim"
    caminho = tmp_path / "s.json"
    caminho.write_text(json.dumps(dados), encoding="utf-8")
    with pytest.raises(ValueError, match="não segue lamp-sim/selecao@2"):
        ler_selecao(caminho)


def test_acesso_e_assinatura():
    assert acesso_do_id("AP018671.1") == "AP018671.1"
    assert acesso_do_id("NZ_CP012345.1") == "NZ_CP012345.1"
    assert acesso_do_id("MN908947.3") == "MN908947.3"
    assert acesso_do_id("AP018671") is None
    assert acesso_do_id("alvo_sintetico_ok") is None
    assert acesso_do_id(None) is None
    assert assinatura("abc", "teste", {"b": 1, "a": "é"}) == {
        "versao_codigo": "abc",
        "experimento": "teste",
        "config_sha256": hashlib.sha256('{"a": "é", "b": 1}'.encode()).hexdigest(),
    }
