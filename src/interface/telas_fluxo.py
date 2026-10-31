"""Dados das telas do fluxo guiado (landing, passos e resultado), sem FastAPI."""

from esteira.ponte import paleta_progresso
from interface.formato import juntar, milhar, plural
from interface.sessao import Execucao, Sessao
from interface.telas import resumo_violacoes, tela_analise, tela_inicio, tela_parametros

AVISO_SIMULADO = "Resultados simulados e provisórios até a validação em bancada."
RODAPE = "Projeto LAMP · Senac + FMABC"
ACEITA_PRIMERS = ".csv,.fasta,.fa,.fna"
ACEITA_ALVO = ".fasta,.fa,.fna"
ERRO_SEM_ARQUIVO = "Escolha um arquivo antes de enviar."
LINKS_MARCA = (
    ("Início", "/", "landing"),
    ("Testar primers", "/testar", "testar"),
    ("Análise", "/analise", "analise"),
    ("Parâmetros", "/parametros", "parametros"),
)

ABAS_DETALHE = (
    ("resumo", "Resumo", ""),
    ("regioes", "Regiões no alvo", "regioes"),
    ("primers", "Primers", "primers"),
    ("violacoes", "Violações", "violacoes"),
    ("dimeros", "Dímeros", "dimeros"),
)

ABAS_PARAMETROS = (
    ("calibracao", "Calibração de Tm", ""),
    ("reacao", "Condições de reação", "reacao"),
    ("distancias", "Distâncias", "distancias"),
    ("regras", "Regras por primer", "regras"),
    ("desenho", "Desenho", "desenho"),
)

PH_INICIAL_PADRAO = 8.7
QUEDA_PH_PADRAO = 2.4
ESPERA_FLUXO_MS_PADRAO = 3000

TUBOS_LANDING = (
    {"rotulo": "Carga alta", "positivo": True, "atraso_ms": 0},
    {"rotulo": "Carga média", "positivo": True, "atraso_ms": 450},
    {"rotulo": "Carga baixa", "positivo": True, "atraso_ms": 900},
    {"rotulo": "Uma única cópia", "positivo": True, "atraso_ms": 1350},
    {"rotulo": "Controle sem DNA", "positivo": False, "atraso_ms": 0},
)

PASSO_TOTAL = 2
ACAO_PRIMERS = "/testar/primers"
ACAO_ALVO = "/testar/alvo"
ID_CAMPO_PRIMERS = "arquivo-primers"
ID_CAMPO_ALVO = "arquivo-alvo"
ROTULO_BOTAO_PRIMERS = "Escolher arquivo de primers"
ROTULO_BOTAO_ALVO = "Escolher arquivo do gene-alvo"

TOM_POSITIVO = "positivo"
TOM_ALERTA = "alerta"
TOM_NEGATIVO = "negativo"
FRASE_NAO_SIMULADO = "Não foi possível simular esta reação."
FRASE_NAO_ENCONTRADOS = "Estes primers não foram encontrados neste gene-alvo."
COMPLEMENTO_NAO_ENCONTRADOS = (
    "Eles não devem amplificar este alvo. Confira se o arquivo do gene-alvo é o certo."
)
FRASE_NAO_AMPLIFICA = "Estes primers não devem amplificar este alvo."
FRASE_CONTROLE_INVALIDO = "Simulação inválida: um controle não se comportou como esperado."
FRASE_AMPLIFICAM_DEVAGAR = (
    "Estes primers encontram o alvo, mas amplificam devagar: na simulação nenhuma amostra "
    "ficou amarela aos {leitura}."
)
COMPLEMENTO_ALEM_DISSO = "Além disso: {regras}."
FRASE_AMPLIFICA = "Estes primers devem amplificar este alvo."
FRASE_RESSALVAS = "Estes primers devem amplificar este alvo, com ressalvas: {regras}."
ROTULO_REGIOES = "Regiões localizadas"
ROTULO_VIOLACOES = "Violações"
ROTULO_EFICIENCIA = "Eficiência prevista"
DETALHE_EFICIENCIA = "modelo provisório"
TRAVESSAO = "—"
URL_OUTRO = "/testar"


PAGINA_DA_ABA = {
    "resultado": "testar",
    "desenho": "analise",
}


def tela_marca(pagina: str) -> dict:
    """Contexto da marca (barra de links, aviso e rodapé) para a página informada."""
    aba = PAGINA_DA_ABA.get(pagina, pagina)
    return {
        "pagina": pagina,
        "links": [
            {"rotulo": rotulo, "href": href, "atual": aba == chave}
            for rotulo, href, chave in LINKS_MARCA
        ],
        "aviso": AVISO_SIMULADO,
        "rodape": RODAPE,
    }


def tela_landing(cfg_web: dict) -> dict:
    """Dados da landing: cores científicas do vermelho de fenol e os tubos da ilustração."""
    manual = cfg_web.get("tubos", {}).get("manual", {})
    ph_inicial = manual.get("ph_inicial", PH_INICIAL_PADRAO)
    queda_ph = manual.get("queda_ph", QUEDA_PH_PADRAO)
    paleta = paleta_progresso(ph_inicial, queda_ph)
    return {
        "cor_negativo": paleta[0],
        "cor_positivo": paleta[-1],
        "tubos": [dict(t) for t in TUBOS_LANDING],
    }


def resumo_primers(sessao: Sessao) -> dict | None:
    """Nome do arquivo e resumo do conjunto carregado, ou None se não há conjunto."""
    if sessao.conjunto is None:
        return None
    conjunto = tela_inicio(sessao)["conjunto"]
    return {"nome_arquivo": conjunto["nome_arquivo"], "resumo": conjunto["resumo"]}


def tela_passo_primers(sessao: Sessao) -> dict:
    """Dados do passo 1 (envio do conjunto de primers)."""
    return {
        "passo": 1,
        "total": PASSO_TOTAL,
        "acao": ACAO_PRIMERS,
        "aceita": ACEITA_PRIMERS,
        "id_campo": ID_CAMPO_PRIMERS,
        "rotulo_botao": ROTULO_BOTAO_PRIMERS,
        "primers": resumo_primers(sessao),
    }


def espera_fluxo_ms(cfg_web: dict) -> int:
    """Tempo mínimo em ms da camada de espera do passo 2 antes do envio."""
    carregamento = cfg_web.get("carregamento") if isinstance(cfg_web, dict) else None
    if not isinstance(carregamento, dict):
        return ESPERA_FLUXO_MS_PADRAO
    valor = carregamento.get("fluxo_ms")
    if isinstance(valor, bool):
        return ESPERA_FLUXO_MS_PADRAO
    if isinstance(valor, (int, float)) and valor >= 0:
        return int(valor)
    return ESPERA_FLUXO_MS_PADRAO


def tela_passo_alvo(sessao: Sessao, cfg_web: dict | None = None) -> dict:
    """Dados do passo 2 (envio do gene-alvo); a rota só chama com conjunto carregado."""
    cfg = cfg_web if cfg_web is not None else sessao.cfg_web
    return {
        "passo": 2,
        "total": PASSO_TOTAL,
        "acao": ACAO_ALVO,
        "aceita": ACEITA_ALVO,
        "id_campo": ID_CAMPO_ALVO,
        "rotulo_botao": ROTULO_BOTAO_ALVO,
        "primers": resumo_primers(sessao),
        "espera_ms": espera_fluxo_ms(cfg or {}),
    }


def _sensibilidade(tubos: dict) -> str:
    """Frase com a menor quantidade de cópias que positivou, ou aviso de que nenhuma positivou."""
    positivas = [
        round(t["copias_media"])
        for t in tubos.get("tubos", [])
        if t.get("papel") == "amostra" and t.get("diagnostico") == "positivo"
    ]
    if not positivas:
        return f"Na simulação, nenhuma amostra positivou aos {tubos['leitura_texto']}."
    n = min(positivas)
    unidade = "cópia" if n == 1 else "cópias"
    return f"Na simulação, detectam a partir de {milhar(n)} {unidade}."


def veredito(execucao: Execucao, tubos: dict, tm_como_calibracao: bool) -> dict:
    """Veredito em linguagem simples a partir do contexto dos tubos simulados e das violações."""
    if tubos["modo"] != "simulacao":
        if getattr(execucao, "alvo", None) is not None and not getattr(
            getattr(execucao, "resultado", None), "regioes", None
        ):
            return {
                "tom": TOM_NEGATIVO,
                "frase": FRASE_NAO_ENCONTRADOS,
                "complemento": COMPLEMENTO_NAO_ENCONTRADOS,
                "sensibilidade": None,
            }
        return {
            "tom": TOM_ALERTA,
            "frase": FRASE_NAO_SIMULADO,
            "complemento": tubos.get("aviso"),
            "sensibilidade": None,
        }
    sensibilidade = _sensibilidade(tubos)
    origem = tubos.get("origem") or {}
    if origem.get("eficiencia") == 0:
        nao_encontradas = origem.get("nao_encontradas") or []
        complemento = None
        if nao_encontradas:
            complemento = f"Regiões não encontradas no alvo: {juntar(list(nao_encontradas))}."
        return {
            "tom": TOM_NEGATIVO,
            "frase": FRASE_NAO_AMPLIFICA,
            "complemento": complemento,
            "sensibilidade": sensibilidade,
        }
    if any(
        t.get("papel") != "amostra" and t.get("diagnostico") == "invalido"
        for t in tubos.get("tubos", [])
    ):
        return {
            "tom": TOM_ALERTA,
            "frase": FRASE_CONTROLE_INVALIDO,
            "complemento": None,
            "sensibilidade": sensibilidade,
        }
    amostras = [t for t in tubos.get("tubos", []) if t.get("papel") == "amostra"]
    tem_amostra_positiva = any(t.get("diagnostico") == "positivo" for t in amostras)
    eficiencia = origem.get("eficiencia")
    if (
        amostras
        and not tem_amostra_positiva
        and eficiencia is not None
        and eficiencia > 0
    ):
        n = resumo_violacoes(execucao.resultado.violacoes, tm_como_calibracao)["n_violacoes"]
        complemento = (
            COMPLEMENTO_ALEM_DISSO.format(
                regras=plural(n, "regra do guia não cumprida", "regras do guia não cumpridas")
            )
            if n > 0
            else None
        )
        return {
            "tom": TOM_ALERTA,
            "frase": FRASE_AMPLIFICAM_DEVAGAR.format(leitura=tubos["leitura_texto"]),
            "complemento": complemento,
            "sensibilidade": None,
        }
    n = resumo_violacoes(execucao.resultado.violacoes, tm_como_calibracao)["n_violacoes"]
    if n > 0:
        regras = plural(n, "regra do guia não cumprida", "regras do guia não cumpridas")
        return {
            "tom": TOM_ALERTA,
            "frase": FRASE_RESSALVAS.format(regras=regras),
            "complemento": None,
            "sensibilidade": sensibilidade,
        }
    return {
        "tom": TOM_POSITIVO,
        "frase": FRASE_AMPLIFICA,
        "complemento": None,
        "sensibilidade": sensibilidade,
    }


def _kpi(kpis: list[dict], rotulo: str) -> dict:
    """Pega o KPI pelo rótulo."""
    return next(k for k in kpis if k["rotulo"] == rotulo)


def tela_resultado(execucao: Execucao, tubos: dict, cfg_primers: dict, cfg_web: dict) -> dict:
    """Dados da tela de resultado: veredito, três números, links e o contexto dos tubos."""
    tm = bool(cfg_web.get("analise", {}).get("tm_como_calibracao", True))
    kpis = tela_analise(execucao, cfg_primers, tm, cfg_web["heterodimeros"]["escala_kcal"])["kpis"]
    kpi_reg = _kpi(kpis, ROTULO_REGIOES)
    kpi_viol = _kpi(kpis, ROTULO_VIOLACOES)
    origem = tubos.get("origem")
    texto_ef = origem["eficiencia_texto"] if origem else TRAVESSAO
    return {
        "id": execucao.id,
        "conjunto": execucao.rotulo_conjunto,
        "alvo": execucao.rotulo_alvo,
        "par": f"{execucao.rotulo_conjunto} × {execucao.rotulo_alvo}",
        "veredito": veredito(execucao, tubos, tm),
        "numeros": [
            {"rotulo": ROTULO_REGIOES, "valor": kpi_reg["valor"], "detalhe": kpi_reg["detalhe"]},
            {"rotulo": ROTULO_VIOLACOES, "valor": kpi_viol["valor"], "detalhe": kpi_viol["detalhe"]},
            {"rotulo": ROTULO_EFICIENCIA, "valor": texto_ef, "detalhe": DETALHE_EFICIENCIA},
        ],
        "url_analise": f"/analise/{execucao.id}",
        "url_exportar": f"/analise/{execucao.id}/exportar",
        "url_outro": URL_OUTRO,
        "simulacao": tubos["modo"] == "simulacao",
        "tubos": tubos,
    }


def tela_detalhe_vazia() -> dict:
    """Aba Análise sem nenhuma análise na sessão."""
    return {
        "titulo": "Nenhuma análise ainda",
        "texto": "Envie um conjunto de primers e um gene-alvo para ver a análise detalhada aqui.",
        "url_acao": "/testar",
        "rotulo_acao": "Testar primers",
    }


def tela_desenho_vazio() -> dict:
    """Aba Desenho sem nenhum desenho na sessão."""
    return {
        "titulo": "Nenhum desenho ainda",
        "texto": (
            'Abra uma análise e use "Desenhar alternativas" para o LampSim propor '
            "outros conjuntos de primers para o mesmo alvo."
        ),
        "url_acao": "/analise",
        "rotulo_acao": "Ir para a análise",
    }


def tela_detalhe(
    execucao: Execucao,
    aba: str,
    cfg_primers: dict,
    cfg_web: dict,
    tm_calibracao: int | None = None,
) -> dict:
    """Contexto da análise detalhada (visual LampSim): o de tela_analise mais as abas."""
    if tm_calibracao is None:
        tm = bool(cfg_web.get("analise", {}).get("tm_como_calibracao", True))
    else:
        tm = bool(tm_calibracao)

    ctx = tela_analise(
        execucao,
        cfg_primers,
        tm,
        cfg_web["heterodimeros"]["escala_kcal"],
    )

    url_base_aba = f"/analise/{execucao.id}" if aba == "resumo" else f"/analise/{execucao.id}/{aba}"
    ctx["href_alternar_calibracao"] = f"{url_base_aba}?tm_calibracao={0 if tm else 1}"

    if execucao.resultado is not None:
        total_viol = resumo_violacoes(execucao.resultado.violacoes, tm)["n_violacoes"]
    else:
        total_viol = 0

    abas = []
    for chave_aba, rotulo_aba, sufixo in ABAS_DETALHE:
        base_href = f"/analise/{execucao.id}" if not sufixo else f"/analise/{execucao.id}/{sufixo}"
        href = f"{base_href}?tm_calibracao={tm_calibracao}" if tm_calibracao is not None else base_href
        abas.append({
            "chave": chave_aba,
            "rotulo": rotulo_aba,
            "href": href,
            "atual": chave_aba == aba,
            "contagem": total_viol if chave_aba == "violacoes" else None,
        })

    par = (
        f"{execucao.rotulo_conjunto} × {execucao.rotulo_alvo}"
        if execucao.alvo is not None
        else execucao.rotulo_conjunto
    )

    limite = int(cfg_web.get("entrada", {}).get("alvo_desenho_max_pb", 2000))
    if execucao.alvo is not None:
        n = len(execucao.alvo.sequencia)
        pode_desenhar = n <= limite
        motivo_sem_desenho = (
            ""
            if pode_desenhar
            else (
                f"O desenho de alternativas aceita alvos de até {milhar(limite)} pb; "
                f"este tem {milhar(n)} pb."
            )
        )
    else:
        pode_desenhar = False
        motivo_sem_desenho = ""

    ctx.update({
        "pode_desenhar": pode_desenhar,
        "motivo_sem_desenho": motivo_sem_desenho,
        "aba": aba,
        "abas": abas,
        "url_resultado": f"/testar/resultado/{execucao.id}",
        "url_exportar": f"/analise/{execucao.id}/exportar",
        "url_outro": URL_OUTRO,
        "par": par,
    })
    return ctx


def tela_parametros_marca(sessao: Sessao, cfg_primers: dict, aba: str) -> dict:
    """Contexto de Parâmetros no visual LampSim: o de telas.tela_parametros mais as abas."""
    ctx = tela_parametros(sessao, cfg_primers)
    abas = [
        {
            "chave": chave,
            "rotulo": rotulo,
            "href": f"/parametros/{sufixo}" if sufixo else "/parametros",
            "atual": chave == aba,
        }
        for chave, rotulo, sufixo in ABAS_PARAMETROS
    ]
    ctx.update({
        "aba": aba,
        "abas": abas,
        "url_testar": "/testar",
    })
    return ctx
