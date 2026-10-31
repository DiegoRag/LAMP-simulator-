"""Testes da Fase 1 - Sessão 3: coordenadas, tela de análise e telas do fluxo guiado."""

import tomllib
from pathlib import Path
from types import SimpleNamespace

from Bio import SeqIO

from esteira.ponte import carregar_config as carregar_config_esteira
from esteira.ponte import paleta_progresso
from interface.formato import coord
from interface.sessao import Sessao
from interface.telas import tela_analise, tela_tubos
from interface.telas_fluxo import (
    ACEITA_ALVO,
    ACEITA_PRIMERS,
    AVISO_SIMULADO,
    ESPERA_FLUXO_MS_PADRAO,
    RODAPE,
    espera_fluxo_ms,
    resumo_primers,
    tela_landing,
    tela_marca,
    tela_passo_alvo,
    tela_passo_primers,
    tela_resultado,
    veredito,
)
from primers.entrada import carregar_config as carregar_config_primers
from primers.modelo import rc

RAIZ = Path(__file__).resolve().parents[1]
DADOS = RAIZ / "tests" / "dados"

CFG_PRIMERS = carregar_config_primers(RAIZ / "config" / "primers.toml")
CFG_WEB_TESTE = {
    "entrada": {"tamanho_max_kb": 64, "alvo_desenho_max_pb": 2000},
    "analise": {"tm_como_calibracao": True},
    "heterodimeros": {"escala_kcal": [-10.0, -8.0, -6.0, -4.0]},
    "exemplos": [
        {
            "nome": "sintético",
            "conjunto": "tests/dados/conjunto_ok.csv",
            "alvo": "tests/dados/alvo_ok.fasta",
        },
    ],
}
ESCALA_TESTE = [-10.0, -8.0, -6.0, -4.0]

PH_REF = 8.7
QUEDA_REF = 2.4
TOTAL_PASSOS = 2
PASSO_UM = 1
PASSO_DOIS = 2
COPIAS_MIL = 1000
QTD_TUBOS_LANDING = 5
ESPERA_ZERO = 0
ESPERA_INVALIDA_NEGATIVA = -5


def test_coord_base_1_genbank():
    """Valida formatação de coordenadas no padrão GenBank (base 1)."""
    assert coord(0, 21) == "1–21"
    assert coord(31, 49) == "32–49"
    assert coord(187, 206) == "188–206"


def test_tela_analise_regiao_com_um_erro_mantem_7_distancias_e_mostra_ponta_3p():
    """Valida preservação de distâncias e exibição de ponta 3' em alvo com substituição."""
    sessao = Sessao(RAIZ, CFG_PRIMERS, CFG_WEB_TESTE)
    conteudo_conj = (DADOS / "conjunto_ok.csv").read_bytes()
    sessao.carregar_conjunto("conjunto_ok.csv", conteudo_conj)

    registros = list(SeqIO.parse(str(DADOS / "alvo_ok.fasta"), "fasta"))
    seq_orig = str(registros[0].seq)
    assert seq_orig[48] == "C"
    seq_mut = seq_orig[:48] + "A" + seq_orig[49:]
    alvo_bytes = f">{registros[0].id}\n{seq_mut}\n".encode("utf-8")
    sessao.carregar_alvo("alvo_mut.fasta", alvo_bytes)

    ex = sessao.analisar()
    dados = tela_analise(ex, CFG_PRIMERS, tm_como_calibracao=True, escala_kcal=ESCALA_TESTE)

    assert dados["kpis"][2] == {
        "rotulo": "Regiões localizadas",
        "valor": "5 de 6",
        "detalhe": "F2 com 1 dif., 1 na ponta",
        "tipo": "neutro",
    }
    assert dados["kpis"][3]["valor"] == "7 de 7"
    assert len(dados["distancias"]) == 7

    nomes_dist = [d["nome"] for d in dados["distancias"]]
    for esperado in ("F2 → B2", "F2 → F1", "F3 → F2", "F3 → B3"):
        assert esperado in nomes_dist

    itens_viol = dados["grupos_violacoes"][0]["itens"]
    viol_f2 = next(item for item in itens_viol if item["item"] == "F2")
    assert viol_f2["texto"] == (
        "Não localizado: melhor posição 32–49 (fita +), 1 diferença (base 18 de 18 do F2: C no primer, "
        "o alvo pede A, a 0 nt da ponta 3′, dentro das 5 bases da ponta; posição 49)"
    )


def test_tela_analise_alvo_fita_menos():
    """Valida análise quando o alvo está na fita complementar reversa."""
    sessao = Sessao(RAIZ, CFG_PRIMERS, CFG_WEB_TESTE)
    conteudo_conj = (DADOS / "conjunto_ok.csv").read_bytes()
    sessao.carregar_conjunto("conjunto_ok.csv", conteudo_conj)

    registros = list(SeqIO.parse(str(DADOS / "alvo_ok.fasta"), "fasta"))
    seq_rc = rc(str(registros[0].seq))
    alvo_bytes = f">{registros[0].id}_rc\n{seq_rc}\n".encode("utf-8")
    sessao.carregar_alvo("alvo_rc.fasta", alvo_bytes)

    ex = sessao.analisar()
    dados = tela_analise(ex, CFG_PRIMERS, tm_como_calibracao=True, escala_kcal=ESCALA_TESTE)

    assert "fita −" in dados["subtitulo"]
    assert dados["kpis"][2]["valor"] == "6 de 6"
    assert dados["kpis"][3]["valor"] == "7 de 7"
    assert any("fita complementar (−)" in a for a in dados["avisos"])


def _montar_sessao_real():
    with open(RAIZ / "config" / "web.toml", "rb") as f:
        cfg_web = tomllib.load(f)
    cfg_primers = carregar_config_primers(RAIZ / "config" / "primers.toml")
    cfg_esteira = carregar_config_esteira(RAIZ / "config" / "esteira.toml")
    sessao = Sessao(RAIZ, cfg_primers, cfg_web)
    return sessao, cfg_primers, cfg_web, cfg_esteira


def test_tela_marca():
    """Valida estrutura de links, página atual e rodapé da marca."""
    marca = tela_marca("testar")
    assert [(link["rotulo"], link["href"]) for link in marca["links"]] == [
        ("Início", "/"),
        ("Testar primers", "/testar"),
        ("Análise", "/analise"),
        ("Parâmetros", "/parametros"),
    ]
    assert [link["atual"] for link in marca["links"]] == [False, True, False, False]
    assert marca["pagina"] == "testar"
    assert marca["aviso"] == AVISO_SIMULADO
    assert marca["rodape"] == RODAPE

    marca_res = tela_marca("resultado")
    assert [link["atual"] for link in marca_res["links"]] == [False, True, False, False]

    marca_des = tela_marca("desenho")
    assert [link["atual"] for link in marca_des["links"]] == [False, False, True, False]


def test_tela_landing():
    """Valida paleta e tubos ilustrativos da landing com config e fallback."""
    _, _, cfg_web, _ = _montar_sessao_real()
    paleta_esperada = paleta_progresso(PH_REF, QUEDA_REF)

    for cfg in (cfg_web, {}):
        tela = tela_landing(cfg)
        assert tela["cor_negativo"] == paleta_esperada[0]
        assert tela["cor_positivo"] == paleta_esperada[-1]

    tubos_landing = tela_landing(cfg_web)["tubos"]
    assert len(tubos_landing) == QTD_TUBOS_LANDING
    assert tubos_landing[-1]["positivo"] is False
    assert tubos_landing[-1]["rotulo"] == "Controle sem DNA"
    assert all(t["positivo"] for t in tubos_landing[:-1])


def test_resumo_primers():
    """Valida resumo de primers quando ausente e após carregamento."""
    sessao, _, _, _ = _montar_sessao_real()
    assert resumo_primers(sessao) is None

    sessao.carregar_conjunto("conjunto_ok.csv", (DADOS / "conjunto_ok.csv").read_bytes())
    resumo = resumo_primers(sessao)
    assert resumo is not None
    assert resumo["nome_arquivo"] == "conjunto_ok.csv"
    assert resumo["resumo"].startswith("CSV · ")


def test_tela_passo_primers_e_alvo():
    """Valida os contextos dos passos 1 e 2 do fluxo guiado."""
    sessao, _, _, _ = _montar_sessao_real()

    p1 = tela_passo_primers(sessao)
    assert p1["passo"] == PASSO_UM
    assert p1["total"] == TOTAL_PASSOS
    assert p1["acao"] == "/testar/primers"
    assert p1["aceita"] == ACEITA_PRIMERS
    assert p1["id_campo"] == "arquivo-primers"
    assert p1["rotulo_botao"] == "Escolher arquivo de primers"
    assert p1["primers"] is None

    sessao.carregar_conjunto("conjunto_ok.csv", (DADOS / "conjunto_ok.csv").read_bytes())
    p1_carregado = tela_passo_primers(sessao)
    assert p1_carregado["primers"] is not None
    assert p1_carregado["primers"]["nome_arquivo"] == "conjunto_ok.csv"

    p2 = tela_passo_alvo(sessao)
    assert p2["passo"] == PASSO_DOIS
    assert p2["total"] == TOTAL_PASSOS
    assert p2["acao"] == "/testar/alvo"
    assert p2["aceita"] == ACEITA_ALVO
    assert p2["id_campo"] == "arquivo-alvo"
    assert p2["rotulo_botao"] == "Escolher arquivo do gene-alvo"
    assert p2["primers"]["nome_arquivo"] == "conjunto_ok.csv"
    assert p2["espera_ms"] == ESPERA_FLUXO_MS_PADRAO
    p2_custom = tela_passo_alvo(sessao, {"carregamento": {"fluxo_ms": ESPERA_ZERO}})
    assert p2_custom["espera_ms"] == ESPERA_ZERO


def test_espera_fluxo_ms():
    """Valida leitura do tempo de espera com padrão, valor zero e entradas inválidas."""
    sessao, _, _, _ = _montar_sessao_real()
    assert espera_fluxo_ms({}) == ESPERA_FLUXO_MS_PADRAO
    assert espera_fluxo_ms({"carregamento": {"fluxo_ms": ESPERA_ZERO}}) == ESPERA_ZERO
    for invalido in (True, ESPERA_INVALIDA_NEGATIVA, "3000"):
        assert espera_fluxo_ms({"carregamento": {"fluxo_ms": invalido}}) == ESPERA_FLUXO_MS_PADRAO
    assert tela_passo_alvo(sessao)["espera_ms"] == ESPERA_FLUXO_MS_PADRAO


def test_veredito_par_real():
    """Valida veredito do par sintético conjunto_ok x alvo_ok."""
    sessao, _, cfg_web, cfg_esteira = _montar_sessao_real()
    sessao.carregar_conjunto("conjunto_ok.csv", (DADOS / "conjunto_ok.csv").read_bytes())
    sessao.carregar_alvo("alvo_ok.fasta", (DADOS / "alvo_ok.fasta").read_bytes())
    execucao = sessao.analisar()
    tubos = tela_tubos(cfg_web, cfg_esteira, execucao, "simulacao")

    v = veredito(execucao, tubos, tm_como_calibracao=True)
    assert v["tom"] == "alerta"
    assert v["frase"] == (
        "Estes primers devem amplificar este alvo, com ressalvas: 2 regras do guia não cumpridas."
    )
    assert v["sensibilidade"] == "Na simulação, detectam a partir de 1 cópia."
    assert v["complemento"] is None


def test_veredito_ramos_manuais():
    """Valida os diferentes ramos de veredito usando contextos simulados."""
    exec_sem_viol = SimpleNamespace(resultado=SimpleNamespace(violacoes=[]))

    tubos_manual = {
        "modo": "manual",
        "origem": {"eficiencia": 1.0, "nao_encontradas": []},
        "tubos": [],
        "leitura_texto": "15 min",
        "aviso": "Curvas manuais de teste.",
    }
    v_manual = veredito(exec_sem_viol, tubos_manual, tm_como_calibracao=True)
    assert v_manual["tom"] == "alerta"
    assert v_manual["frase"] == "Não foi possível simular esta reação."
    assert v_manual["complemento"] == "Curvas manuais de teste."
    assert v_manual["sensibilidade"] is None

    tubos_ef0 = {
        "modo": "simulacao",
        "origem": {"eficiencia": 0.0, "nao_encontradas": []},
        "tubos": [{"papel": "amostra", "diagnostico": "negativo", "copias_media": 10.0}],
        "leitura_texto": "15 min",
        "aviso": "",
    }
    v_ef0 = veredito(exec_sem_viol, tubos_ef0, tm_como_calibracao=True)
    assert v_ef0["tom"] == "negativo"
    assert v_ef0["frase"] == "Estes primers não devem amplificar este alvo."
    assert v_ef0["complemento"] is None

    tubos_ef0_regioes = {
        "modo": "simulacao",
        "origem": {"eficiencia": 0.0, "nao_encontradas": ["F3", "B3"]},
        "tubos": [],
        "leitura_texto": "15 min",
        "aviso": "",
    }
    v_ef0_reg = veredito(exec_sem_viol, tubos_ef0_regioes, tm_como_calibracao=True)
    assert v_ef0_reg["tom"] == "negativo"
    assert v_ef0_reg["complemento"] == "Regiões não encontradas no alvo: F3 e B3."

    tubos_controle_invalido = {
        "modo": "simulacao",
        "origem": {"eficiencia": 1.0, "nao_encontradas": []},
        "tubos": [
            {"papel": "amostra", "diagnostico": "positivo", "copias_media": 100.0},
            {"papel": "controle_negativo", "diagnostico": "invalido", "copias_media": 0.0},
        ],
        "leitura_texto": "15 min",
        "aviso": "",
    }
    v_ctrl = veredito(exec_sem_viol, tubos_controle_invalido, tm_como_calibracao=True)
    assert v_ctrl["tom"] == "alerta"
    assert v_ctrl["frase"] == "Simulação inválida: um controle não se comportou como esperado."

    tubos_ok = {
        "modo": "simulacao",
        "origem": {"eficiencia": 1.0, "nao_encontradas": []},
        "tubos": [{"papel": "amostra", "diagnostico": "positivo", "copias_media": 10.0}],
        "leitura_texto": "15 min",
        "aviso": "",
    }
    v_ok = veredito(exec_sem_viol, tubos_ok, tm_como_calibracao=True)
    assert v_ok["tom"] == "positivo"
    assert v_ok["frase"] == "Estes primers devem amplificar este alvo."
    assert v_ok["sensibilidade"] == "Na simulação, detectam a partir de 10 cópias."

    tubos_nenhuma = {
        "modo": "simulacao",
        "origem": {"eficiencia": 1.0, "nao_encontradas": []},
        "tubos": [{"papel": "amostra", "diagnostico": "negativo", "copias_media": float(COPIAS_MIL)}],
        "leitura_texto": "15 min",
        "aviso": "",
    }
    v_nenhuma = veredito(exec_sem_viol, tubos_nenhuma, tm_como_calibracao=True)
    assert v_nenhuma["tom"] == "alerta"
    assert v_nenhuma["frase"] == (
        "Estes primers encontram o alvo, mas amplificam devagar: na simulação nenhuma amostra "
        "ficou amarela aos 15 min."
    )
    assert v_nenhuma["sensibilidade"] is None

    tubos_mil = {
        "modo": "simulacao",
        "origem": {"eficiencia": 1.0, "nao_encontradas": []},
        "tubos": [
            {"papel": "amostra", "diagnostico": "positivo", "copias_media": float(COPIAS_MIL * 10)},
            {"papel": "amostra", "diagnostico": "positivo", "copias_media": float(COPIAS_MIL)},
            {"papel": "amostra", "diagnostico": "negativo", "copias_media": 1.0},
        ],
        "leitura_texto": "15 min",
        "aviso": "",
    }
    v_mil = veredito(exec_sem_viol, tubos_mil, tm_como_calibracao=True)
    assert v_mil["sensibilidade"] == "Na simulação, detectam a partir de 1.000 cópias."


def test_tela_resultado():
    """Valida a montagem do dicionário de contexto da tela de resultado."""
    sessao, cfg_primers, cfg_web, cfg_esteira = _montar_sessao_real()
    sessao.carregar_conjunto("conjunto_ok.csv", (DADOS / "conjunto_ok.csv").read_bytes())
    sessao.carregar_alvo("alvo_ok.fasta", (DADOS / "alvo_ok.fasta").read_bytes())
    execucao = sessao.analisar()
    tubos = tela_tubos(cfg_web, cfg_esteira, execucao, "simulacao")

    tela = tela_resultado(execucao, tubos, cfg_primers, cfg_web)
    assert tela["id"] == 1
    assert tela["par"] == "conjunto_ok × alvo_ok"
    assert tela["conjunto"] == "conjunto_ok"
    assert tela["alvo"] == "alvo_ok"
    assert [n["rotulo"] for n in tela["numeros"]] == [
        "Regiões localizadas",
        "Violações",
        "Eficiência prevista",
    ]
    assert [n["valor"] for n in tela["numeros"]] == ["6 de 6", "2", "100 %"]
    assert tela["numeros"][1]["detalhe"] == "e 6 de Tm a calibrar"
    assert tela["url_analise"] == "/analise/1"
    assert tela["url_exportar"] == "/analise/1/exportar"
    assert tela["url_outro"] == "/testar"
    assert tela["simulacao"] is True
    assert tela["tubos"] is tubos
