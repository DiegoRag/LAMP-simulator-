"""Fase 4, sessão 5: chaves novas do config, faixa relaxada e diagnóstico do desenho vazio."""

from pathlib import Path

from primers.candidatos import REGRAS, TIPOS_BASE
from primers.diagnostico import config_relaxada, diagnosticar
from primers.entrada import carregar_config

CONFIG = Path(__file__).parents[1] / "config" / "primers.toml"
SUGESTAO = " Tente o desenho com as faixas relaxadas de [desenho_relaxado]."


def _stats(cands: dict[str, int], dominante: dict[str, str]) -> dict:
    elim = {t: dict.fromkeys(REGRAS, 0) for t in TIPOS_BASE}
    for t, regra in dominante.items():
        elim[t][regra] = 100
        elim[t]["extremidade"] = 7
    return {"candidatos_por_tipo": cands, "eliminados_por_tipo": elim}


def test_chaves_do_desenho():
    cfg = carregar_config(CONFIG)
    d = cfg["desenho"]
    assert (d["janela_nt"], d["candidatos_por_janela"], d["nucleos_guardados"]) == (10, 2, 5000)
    assert (d["sobreposicao_max_nt"], d["pre_selecao"], d["n_resultados"]) == (10, 50, 5)
    assert (d["grampo_margem_C"], d["homodimero_margem_C"], d["peso_dimero"]) == (0.0, 0.0, 0.5)
    assert "candidatos_por_regiao" not in d  # teto global trocado pelo teto por janela
    assert "conjuntos_examinados" not in d  # a busca não para mais num número de combinações
    assert cfg["desenho_relaxado"] == {"gc_pct": [30.0, 70.0], "tm_folga_C": 1.0}
    assert cfg["dimeros"]["max_nt_estrutura"] == 60


def test_config_relaxada():
    cfg = carregar_config(CONFIG)
    rel = config_relaxada(cfg)
    assert rel["gc_pct"] == {"min": 30.0, "max": 70.0}
    assert rel["tm_c"] == {
        "F1c_B1c": [63.0, 67.0], "F2_B2": [58.0, 62.0], "F3_B3": [58.0, 62.0], "loops": [63.0, 67.0],
    }
    assert cfg["gc_pct"] == {"min": 40.0, "max": 65.0}  # o cfg original não muda
    assert rel["desenho"] == cfg["desenho"]


def test_diagnostico_sem_candidatos():
    todos = dict.fromkeys(TIPOS_BASE, 0)
    st = _stats(todos, dict.fromkeys(TIPOS_BASE, "gc_pct"))
    assert diagnosticar(st, False) == (
        "Nenhum candidato a F3, F2, F1, B1c, B2c e B3c: a regra que mais eliminou janelas foi "
        "GC fora da faixa de [gc_pct]." + SUGESTAO
    )
    alguns = {"F3": 5, "F2": 5, "F1": 0, "B1c": 0, "B2c": 4, "B3c": 0}
    st = _stats(alguns, {"F1": "tm_c", "B1c": "tm_c", "B3c": "repetida"})
    assert diagnosticar(st, True) == (
        "Nenhum candidato a F1 e B1c: a regra que mais eliminou janelas foi Tm de referência fora da "
        "faixa de [tm_c]. Nenhum candidato a B3c: a regra que mais eliminou janelas foi sequência "
        "repetida no alvo ou no complemento reverso."
    )


def test_diagnostico_empate_alvo_curto_e_geometria():
    st = _stats(dict.fromkeys(TIPOS_BASE, 0), {})
    for t in TIPOS_BASE:
        st["eliminados_por_tipo"][t]["tm_c"] = 9
        st["eliminados_por_tipo"][t]["gc_pct"] = 9  # empate: vale a regra que vem antes em REGRAS
    assert diagnosticar(st, True).endswith("foi GC fora da faixa de [gc_pct].")
    vazio = _stats(dict.fromkeys(TIPOS_BASE, 0), {})
    assert diagnosticar(vazio, True) == (
        "Nenhuma janela para F3, F2, F1, B1c, B2c e B3c: o alvo é curto demais."
    )
    cheio = _stats({"F3": 9, "F2": 9, "F1": 4, "B1c": 3, "B2c": 3, "B3c": 9}, {})
    assert diagnosticar(cheio, False) == (
        "Há candidatos para as seis regiões (o mais escasso é B1c, com 3), mas nenhuma combinação "
        "respeita as distâncias de [distancias]." + SUGESTAO
    )
