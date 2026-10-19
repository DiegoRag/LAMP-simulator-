"""Fase 5, sessão 5: ponte análise -> eficiência -> tubos simulados (esteira.ponte)."""

from pathlib import Path

import pytest

from esteira.contrato import ESQUEMA_V2, validar_selecao
from esteira.ponte import (
    carregar_config,
    eficiencia_prevista,
    eixo_minutos,
    medida_por_alelo,
    paleta_progresso,
    sementes_dos_tubos,
    simular_tubos,
)
from primers.incompat import Incompat, IncompatRegiao

RAIZ = Path(__file__).resolve().parents[1]
CFG = carregar_config(RAIZ / "config" / "esteira.toml")
ESSENCIAIS = ["F3", "F2", "F1", "B1c", "B2c", "B3c"]


def _inc(itens: list[IncompatRegiao]) -> Incompat:
    faltam = tuple(i.regiao for i in itens if not i.encontrada)
    if faltam:
        return Incompat(tuple(itens), None, None, faltam)
    return Incompat(tuple(itens), sum(i.total for i in itens), sum(i.ponta3 for i in itens), ())


PERFEITO = _inc([IncompatRegiao(n, True, 0, 0) for n in ESSENCIAIS])
SEIS_NA_PONTA = _inc([IncompatRegiao(n, True, 1, 1) for n in ESSENCIAIS])


def test_config_da_ponte():
    assert CFG["eficiencia"] == {"fator_ponta3": 0.3, "fator_outra": 0.85, "piso_laco": 0.5}
    assert carregar_config() == CFG  # caminho padrão: config/esteira.toml da raiz


def test_eficiencia_prevista():
    assert eficiencia_prevista(None, CFG) is None
    assert eficiencia_prevista(PERFEITO, CFG) == 1.0
    assert eficiencia_prevista(SEIS_NA_PONTA, CFG) == pytest.approx(0.3 ** 6)
    tres_fora = _inc([IncompatRegiao(n, True, int(n in ("F3", "F1", "B1c")), 0) for n in ESSENCIAIS])
    assert eficiencia_prevista(tres_fora, CFG) == pytest.approx(0.85 ** 3)
    mista = _inc([IncompatRegiao("F3", True, 3, 1)] + [IncompatRegiao(n, True, 0, 0) for n in ESSENCIAIS[1:]])
    assert eficiencia_prevista(mista, CFG) == pytest.approx(0.3 * 0.85 ** 2)


def test_regiao_essencial_ausente_zera_e_laco_ausente_tem_piso():
    sem_f2 = _inc([IncompatRegiao(n, n != "F2", 0 if n != "F2" else None, 0 if n != "F2" else None)
                   for n in ESSENCIAIS])
    assert eficiencia_prevista(sem_f2, CFG) == 0.0
    limpas = [IncompatRegiao(n, True, 0, 0) for n in ESSENCIAIS]
    sem_lf = _inc(limpas + [IncompatRegiao("LF", False, None, None)])
    assert eficiencia_prevista(sem_lf, CFG) == 0.5
    lb_ruim = _inc(limpas + [IncompatRegiao("LB", True, 1, 1)])
    assert eficiencia_prevista(lb_ruim, CFG) == 0.5  # 0,3 fica no piso de 0,5


def test_medida_segue_o_contrato_v2():
    medida = medida_por_alelo(SEIS_NA_PONTA, "alvo_sintetico", None, CFG)
    assert medida == {
        "alelo_id": "alvo_sintetico",
        "acesso": None,
        "incompat_total": 6,
        "incompat_3p": 6,
        "f_previsto": pytest.approx(0.000729),
    }
    selecao = {
        "esquema": ESQUEMA_V2,
        "gerado_por": {"versao_codigo": "teste", "experimento": "fase5", "config_sha256": "0"},
        "conjuntos": [{
            "conjunto_id": "sintetico",
            "origem": "desenho",
            "aviso": "candidato computacional",
            "primers": {"F3": "ACGT", "B3": "ACGT", "FIP": "ACGT", "BIP": "ACGT"},
            "por_alelo": [medida, medida_por_alelo(None, "sem_alvo", None, CFG)],
            "selecionado": True,
            "motivo": "teste",
        }],
    }
    assert validar_selecao(selecao) == []


def test_eixo_e_sementes():
    t = eixo_minutos(60.0, 0.5)
    assert t.size == 121 and t[0] == 0.0 and t[-1] == 60.0
    assert sementes_dos_tubos(2026, 3) == [479243620, 454024514, 1818751028]


SERIE = [
    {"nome": "10⁵ cópias", "copias": 1e5, "papel": "amostra"},
    {"nome": "Controle negativo", "copias": 0, "papel": "controle_negativo"},
    {"nome": "Controle positivo", "copias": 1e5, "papel": "controle_positivo", "eficiencia": 1.0},
]


def test_tubos_do_conjunto_perfeito_e_do_conjunto_com_seis_erros():
    t = eixo_minutos(60.0, 0.5)
    bom = simular_tubos(1.0, SERIE, t, 2026)
    ruim = simular_tubos(0.3 ** 6, SERIE, t, 2026)
    assert [x["diagnostico"] for x in bom] == ["positivo", "controle_ok", "controle_ok"]
    assert [x["diagnostico"] for x in ruim] == ["negativo", "controle_ok", "controle_ok"]
    assert bom[0]["ph"][-1] < 7.0 and ruim[0]["ph"][-1] > 8.4
    assert ruim[0]["t0_min"] is None and bom[1]["t0_min"] is None
    # o controle positivo tem molde próprio: não depende do conjunto
    assert bom[2]["ph"] == ruim[2]["ph"] and ruim[2]["eficiencia"] == 1.0
    assert all(len(x["ph"]) == 121 and len(x["cores"]) == 121 for x in bom)
    assert bom[0]["copias_media"] == 1e5 and bom[1]["copias"] == 0
    assert simular_tubos(1.0, SERIE, t, 2026) == bom  # mesma semente, mesmos tubos


def test_papel_invalido():
    with pytest.raises(ValueError):
        simular_tubos(1.0, [{"nome": "x", "copias": 1, "papel": "talvez"}], eixo_minutos(10, 0.5), 1)


def test_paleta_das_curvas_manuais():
    paleta = paleta_progresso(8.7, 2.4)
    assert len(paleta) == 101
    assert (paleta[0], paleta[50], paleta[100]) == ("#FE00B1", "#FB8B6B", "#F3DF48")
