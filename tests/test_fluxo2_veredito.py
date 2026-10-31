"""Testes dos vereditos em linguagem simples (alvo não encontrado e amplificação lenta)."""

from types import SimpleNamespace

from interface.telas_fluxo import (
    COMPLEMENTO_ALEM_DISSO,
    COMPLEMENTO_NAO_ENCONTRADOS,
    FRASE_AMPLIFICA,
    FRASE_AMPLIFICAM_DEVAGAR,
    FRASE_NAO_ENCONTRADOS,
    FRASE_NAO_SIMULADO,
    TOM_ALERTA,
    TOM_NEGATIVO,
    TOM_POSITIVO,
    veredito,
)

EFICIENCIA_BAIXA = 0.15
COPIAS_ZERO = 0.0
COPIAS_UMA = 1.0
COPIAS_POSITIVA = 100.0
COPIAS_AMOSTRA = 1000.0
COPIAS_CONTROLE_POS = 10000.0
LEITURA_TEXTO_15 = "15 min"


def test_veredito_modo_manual_alvo_presente_sem_regioes():
    """Valida tom negativo e frase de primers não encontrados quando regiões está vazio ou None."""
    alvo = SimpleNamespace(nome="alvo_teste.fasta")
    tubos_manual = {
        "modo": "manual",
        "origem": None,
        "tubos": [],
        "leitura_texto": LEITURA_TEXTO_15,
        "aviso": "Curvas manuais.",
    }

    # Caso com regioes None
    exec_sem_regioes_none = SimpleNamespace(
        alvo=alvo,
        resultado=SimpleNamespace(regioes=None, violacoes=[]),
    )
    v1 = veredito(exec_sem_regioes_none, tubos_manual, tm_como_calibracao=True)
    assert v1["tom"] == TOM_NEGATIVO
    assert v1["frase"] == FRASE_NAO_ENCONTRADOS
    assert v1["complemento"] == COMPLEMENTO_NAO_ENCONTRADOS
    assert v1["sensibilidade"] is None

    # Caso com regioes dicionário vazio
    exec_sem_regioes_vazio = SimpleNamespace(
        alvo=alvo,
        resultado=SimpleNamespace(regioes={}, violacoes=[]),
    )
    v2 = veredito(exec_sem_regioes_vazio, tubos_manual, tm_como_calibracao=True)
    assert v2["tom"] == TOM_NEGATIVO
    assert v2["frase"] == FRASE_NAO_ENCONTRADOS
    assert v2["complemento"] == COMPLEMENTO_NAO_ENCONTRADOS
    assert v2["sensibilidade"] is None


def test_veredito_modo_manual_alvo_presente_com_regioes():
    """Valida preservação de FRASE_NAO_SIMULADO quando regiões foram localizadas no alvo."""
    alvo = SimpleNamespace(nome="alvo_teste.fasta")
    exec_com_regioes = SimpleNamespace(
        alvo=alvo,
        resultado=SimpleNamespace(regioes={"F3": SimpleNamespace()}, violacoes=[]),
    )
    tubos_manual = {
        "modo": "manual",
        "origem": None,
        "tubos": [],
        "leitura_texto": LEITURA_TEXTO_15,
        "aviso": "Curvas manuais selecionadas.",
    }
    v = veredito(exec_com_regioes, tubos_manual, tm_como_calibracao=True)
    assert v["tom"] == TOM_ALERTA
    assert v["frase"] == FRASE_NAO_SIMULADO
    assert v["complemento"] == "Curvas manuais selecionadas."
    assert v["sensibilidade"] is None


def test_veredito_simulacao_eficiencia_baixa_todas_amostras_negativas():
    """Valida tom alerta e frase de amplificação lenta quando nenhuma amostra positiva a tempo."""
    tubos_lentos = {
        "modo": "simulacao",
        "origem": {"eficiencia": EFICIENCIA_BAIXA, "nao_encontradas": []},
        "tubos": [
            {"papel": "amostra", "diagnostico": "negativo", "copias_media": COPIAS_AMOSTRA},
            {"papel": "controle_negativo", "diagnostico": "controle_ok", "copias_media": COPIAS_ZERO},
            {"papel": "controle_positivo", "diagnostico": "controle_ok", "copias_media": COPIAS_CONTROLE_POS},
        ],
        "leitura_texto": LEITURA_TEXTO_15,
        "aviso": "",
    }

    # Sem violações: complemento deve ser None
    exec_sem_viol = SimpleNamespace(
        alvo=SimpleNamespace(nome="alvo.fasta"),
        resultado=SimpleNamespace(violacoes=[]),
    )
    v_sem_viol = veredito(exec_sem_viol, tubos_lentos, tm_como_calibracao=True)
    assert v_sem_viol["tom"] == TOM_ALERTA
    assert v_sem_viol["frase"] == FRASE_AMPLIFICAM_DEVAGAR.format(leitura=LEITURA_TEXTO_15)
    assert "amplificam devagar" in v_sem_viol["frase"]
    assert "aos 15 min" in v_sem_viol["frase"]
    assert v_sem_viol["complemento"] is None
    assert v_sem_viol["sensibilidade"] is None

    # Com violações: complemento deve detalhar as regras não cumpridas
    exec_com_viol = SimpleNamespace(
        alvo=SimpleNamespace(nome="alvo.fasta"),
        resultado=SimpleNamespace(
            violacoes=[SimpleNamespace(regra="gc"), SimpleNamespace(regra="dg")]
        ),
    )
    v_com_viol = veredito(exec_com_viol, tubos_lentos, tm_como_calibracao=True)
    assert v_com_viol["tom"] == TOM_ALERTA
    assert v_com_viol["frase"] == FRASE_AMPLIFICAM_DEVAGAR.format(leitura=LEITURA_TEXTO_15)
    assert v_com_viol["complemento"] == COMPLEMENTO_ALEM_DISSO.format(
        regras="2 regras do guia não cumpridas"
    )
    assert v_com_viol["sensibilidade"] is None


def test_veredito_simulacao_com_amostra_positiva_mantem_comportamento():
    """Valida que simulações com amostras positivas continuam gerando veredito padrão."""
    tubos_com_positiva = {
        "modo": "simulacao",
        "origem": {"eficiencia": EFICIENCIA_BAIXA, "nao_encontradas": []},
        "tubos": [
            {"papel": "amostra", "diagnostico": "positivo", "copias_media": COPIAS_POSITIVA},
            {"papel": "amostra", "diagnostico": "negativo", "copias_media": COPIAS_UMA},
            {"papel": "controle_negativo", "diagnostico": "controle_ok", "copias_media": COPIAS_ZERO},
        ],
        "leitura_texto": LEITURA_TEXTO_15,
        "aviso": "",
    }
    exec_sem_viol = SimpleNamespace(
        alvo=SimpleNamespace(nome="alvo.fasta"),
        resultado=SimpleNamespace(violacoes=[]),
    )
    v = veredito(exec_sem_viol, tubos_com_positiva, tm_como_calibracao=True)
    assert v["tom"] == TOM_POSITIVO
    assert v["frase"] == FRASE_AMPLIFICA
    assert v["sensibilidade"] == "Na simulação, detectam a partir de 100 cópias."
    assert v["complemento"] is None
