"""Geração dos dicionários de contexto para as telas da interface web."""

import csv
import io
import math
from pathlib import Path

from esteira.ponte import eficiencia_prevista, eixo_minutos, paleta_progresso, simular_tubos
from interface.formato import coord, data_hora, faixa, juntar, milhar, num, num_curto, pct, plural
from interface.sessao import Execucao, Sessao
from lampsim.config import LEITURA, PAPEIS, QUIMICA, TEMPO
from primers.casamento import ErroDeBase, Posicao, erros_de_base, regioes_do_conjunto
from primers.desenho import AVISO, ConjuntoDesenhado, ResultadoDesenho
from primers.entrada import reconhecer_nome_primer
from primers.modelo import Regiao
from primers.restricoes import calcular_distancias
from primers.termo import Condicoes, CondicoesReferencia

AVISO_CURTO = "Candidatos computacionais, não validados experimentalmente (ADR-0013)."

NOMES_DISTANCIA = {
    "F3_F2_gap": "F3 → F2",
    "F2_F1": "F2 → F1",
    "F1_B1c_gap": "F1 → B1c",
    "B1c_B2c": "B1c → B2c",
    "B2c_B3c_gap": "B2c → B3c",
    "F2_B2": "F2 → B2",
    "F3_B3_total": "F3 → B3",
}

PRIMER_PARA_REGIAO = {
    "F3": "F3",
    "F2": "F2",
    "F1c": "F1",
    "B1c": "B1c",
    "B2": "B2c",
    "B3": "B3c",
    "LF": "LF",
    "LB": "LB",
}

REGIAO_PARA_PRIMER = {
    "F3": "F3",
    "F2": "F2",
    "F1": "F1c",
    "B1c": "B1c",
    "B2c": "B2",
    "B3c": "B3",
    "LF": "LF",
    "LB": "LB",
}

FAIXA_TM_PRIMER = {
    "F3": "F3_B3",
    "B3": "F3_B3",
    "F2": "F2_B2",
    "B2": "F2_B2",
    "F1c": "F1c_B1c",
    "B1c": "F1c_B1c",
    "LF": "loops",
    "LB": "loops",
}

PASSOS_EIXO = [
    10,
    20,
    25,
    40,
    50,
    100,
    200,
    250,
    500,
    1000,
    2000,
    2500,
    5000,
    10000,
    20000,
    25000,
    50000,
    100000,
    200000,
    250000,
    500000,
    1000000,
]

TUBOS_PADRAO_TEMPO = {
    "duracao_min": 30.0,
    "duracao_max_min": 60.0,
    "passo_min": 0.5,
    "velocidade": 60.0,
}
TUBOS_PADRAO_VELOCIDADES = [30.0, 60.0, 120.0, 300.0]
TUBOS_PADRAO_SEMENTE = 2026
TUBOS_PADRAO_MANUAL = {"ph_inicial": 8.7, "queda_ph": 2.4}

PADRAO_SERIE = [
    {"nome": "10⁶ cópias", "copias": 1000000.0, "papel": "amostra"},
    {"nome": "10⁵ cópias", "copias": 100000.0, "papel": "amostra"},
    {"nome": "10⁴ cópias", "copias": 10000.0, "papel": "amostra"},
    {"nome": "10³ cópias", "copias": 1000.0, "papel": "amostra"},
    {"nome": "100 cópias", "copias": 100.0, "papel": "amostra"},
    {"nome": "10 cópias", "copias": 10.0, "papel": "amostra"},
    {"nome": "1 cópia", "copias": 1.0, "papel": "amostra"},
    {"nome": "Controle negativo", "copias": 0.0, "papel": "controle_negativo"},
    {"nome": "Controle positivo", "copias": 100000.0, "papel": "controle_positivo", "eficiencia": 1.0},
]

PADRAO_TUBOS = [
    {"nome": "Amostra 1", "curva": "linear", "inicio": 0.0, "duracao": 30.0, "alcance": 100.0},
    {"nome": "Amostra 2", "curva": "sigmoide", "inicio": 24.0, "duracao": 6.0, "alcance": 100.0},
    {"nome": "Amostra 3", "curva": "sigmoide", "inicio": 8.0, "duracao": 12.0, "alcance": 100.0},
    {"nome": "Amostra 4", "curva": "sigmoide", "inicio": 16.0, "duracao": 8.0, "alcance": 100.0},
    {"nome": "Amostra 5", "curva": "linear", "inicio": 10.0, "duracao": 20.0, "alcance": 100.0},
    {"nome": "Amostra 6", "curva": "sigmoide", "inicio": 4.0, "duracao": 6.0, "alcance": 100.0},
    {"nome": "Amostra 7", "curva": "sigmoide", "inicio": 18.0, "duracao": 12.0, "alcance": 60.0},
    {"nome": "Amostra 8", "curva": "linear", "inicio": 0.0, "duracao": 30.0, "alcance": 0.0},
    {"nome": "Referência", "curva": "linear", "inicio": 0.0, "duracao": 30.0, "alcance": 0.0},
]

DIAGNOSTICO_TEXTO = {
    "positivo": "positivo",
    "negativo": "negativo",
    "controle_ok": "controle ok",
    "invalido": "inválido",
}
PAPEL_TEXTO = {
    "amostra": "amostra",
    "controle_negativo": "controle negativo",
    "controle_positivo": "controle positivo",
}
AVISO_TUBOS_MANUAL = (
    "Curvas manuais: digitadas à mão em config/web.toml, não vêm do conjunto analisado nem do "
    "simulador. Servem só para ilustrar."
)
AVISO_TUBOS_SEM_ANALISE = (
    "Nenhuma análise nesta sessão, então não há tubos simulados. Estas curvas são digitadas à "
    "mão em config/web.toml: analise um conjunto contra um alvo para ver a simulação."
)
AVISO_TUBOS_SEM_MEDIDA = (
    "A análise {id} não tem alvo ou não tem o conjunto completo, então não há eficiência "
    "prevista. Estas curvas são digitadas à mão em config/web.toml e não vêm de simulação."
)


FAIXAS_TM_TEXTO = (("F3/B3", "F3_B3"), ("F2/B2", "F2_B2"), ("F1c/B1c", "F1c_B1c"), ("loops", "loops"))


def resumo_violacoes(violacoes: list, tm_como_calibracao: bool) -> dict:
    """Contagem das violações coerente com "a calibrar".

    Com a chave ligada, Tm fora da faixa não conta como violação: aparece à parte, como
    "a calibrar", no KPI da Análise e no histórico do Início, igual à tabela de primers.
    Devolve n_violacoes, n_tm, detalhe (texto do KPI), tipo ("viol", "calibrar" ou "ok")
    e texto (resumo do histórico).
    """
    n_tm = sum(1 for v in violacoes if v.regra == "tm_c")
    if tm_como_calibracao:
        n_viol = len(violacoes) - n_tm
        detalhe = f"e {n_tm} de Tm a calibrar" if n_tm else ""
        tipo = "viol" if n_viol else ("calibrar" if n_tm else "ok")
    else:
        n_viol = len(violacoes)
        detalhe = f"{n_tm} de Tm" if n_tm else ""
        tipo = "viol" if n_viol else "ok"
    texto = plural(n_viol, "violação", "violações")
    if tm_como_calibracao and n_tm:
        texto += f" {detalhe}"
    return {"n_violacoes": n_viol, "n_tm": n_tm, "detalhe": detalhe, "tipo": tipo, "texto": texto}


def regua_distancia(d: float, d_min: float, d_max: float) -> dict:
    """Posições (porcentagem) da faixa e da marca numa régua que cabe a medida e a faixa.

    A régua vai de min(0, d, d_min) a max(1,25 · d_max, d, 1), com 5 % de folga de cada lado:
    a marca nunca sai da régua, nem com distância negativa (regiões sobrepostas) nem com
    faixa larga. "zero" marca a posição do 0 quando a régua tem parte negativa ("" se não).
    """
    ini = min(0.0, float(d), float(d_min))
    fim = max(float(d_max) * 1.25, float(d), 1.0)
    folga = (fim - ini) * 0.05
    ini -= folga
    fim += folga
    span = fim - ini
    return {
        "banda_esquerda": pct(d_min - ini, span),
        "banda_largura": pct(d_max - d_min, span),
        "marca": pct(d - ini, span),
        "zero": pct(-ini, span) if min(d, d_min) < 0 else "",
    }


def texto_fita(fita: str | None) -> str:
    """"+" ou "−" (sinal tipográfico) para a fita em que o conjunto foi localizado."""
    return "−" if fita == "-" else "+"


def texto_erro_de_base(e: ErroDeBase, bases_ponta: int) -> str:
    """Uma diferença de base contada no primer, com a distância até a ponta que estende."""
    if e.ponta == "5'":
        ponta = "da ponta 5′ (a que fecha o laço)"
    else:
        ponta = "da ponta 3′"
    texto = (
        f"base {e.base_no_primer} de {e.comprimento} do {e.primer}: {e.primer_tem} no primer, "
        f"o alvo pede {e.alvo_pede}, a {e.distancia_ponta} nt {ponta}"
    )
    if e.distancia_ponta < bases_ponta:
        texto += f", dentro das {bases_ponta} bases da ponta"
    return f"{texto}; posição {milhar(e.posicao)}"


def texto_faixas_tm(cfg_primers: dict) -> str:
    """Faixas de Tm de referência do config, como "F3/B3 59–61, ..., loops 64–66 °C"."""
    partes = [
        f"{rotulo} {faixa(cfg_primers['tm_c'][chave][0], cfg_primers['tm_c'][chave][1])}"
        for rotulo, chave in FAIXAS_TM_TEXTO
        if chave in cfg_primers.get("tm_c", {})
    ]
    return f"{juntar(partes)} °C" if partes else ""


def lado_de(item: str) -> str:
    """Determina o lado ('F', 'B', 'L', 'X') a partir do nome de primer ou região."""
    item_upper = item.upper()
    if item in ("F3", "F2", "F1c", "F1") or item_upper in ("F3", "F2", "F1C", "F1"):
        return "F"
    if item in ("B3", "B2", "B1c", "B2c", "B3c") or item_upper in ("B3", "B2", "B1C", "B2C", "B3C"):
        return "B"
    if item in ("LF", "LB") or item_upper in ("LF", "LB", "FLP", "BLP", "LOOPF", "LOOPB"):
        return "L"
    return "X"


def marcas_de_eixo(comprimento: int) -> list[dict]:
    """Gera marcas de eixo linear para um determinado comprimento."""
    if comprimento <= 0:
        return []
    passo = PASSOS_EIXO[-1]
    for p in PASSOS_EIXO:
        if comprimento / p <= 8:
            passo = p
            break

    valores: list[int] = []
    v = 0
    while v < comprimento:
        valores.append(v)
        v += passo
    valores.append(comprimento)

    total_pts = len(valores)
    ticks: list[dict] = []
    for i, val in enumerate(valores):
        if i == 0:
            alinhamento = "inicio"
        elif i == total_pts - 1:
            alinhamento = "fim"
        else:
            alinhamento = "meio"
        ticks.append({
            "valor": milhar(val) if val > 0 else "1",  # coordenadas a partir de 1
            "esquerda": pct(val, comprimento),
            "alinhamento": alinhamento,
        })
    return ticks


def tela_base(pagina: str, sessao: Sessao, versao: tuple[str, str]) -> dict:
    """Monta o contexto básico comum a todas as páginas, incluindo etapas e tempos de espera."""
    rot_conj = Path(sessao.conjunto.nome).stem if sessao.conjunto else None
    rot_alvo = Path(sessao.alvo.nome).stem if sessao.alvo else None

    if rot_conj and rot_alvo:
        projeto = f"{rot_conj} × {rot_alvo}"
    elif rot_conj:
        projeto = rot_conj
    elif rot_alvo:
        projeto = rot_alvo
    else:
        projeto = "Nenhum arquivo carregado"

    dados_feitos = sessao.conjunto is not None or sessao.alvo is not None
    analise_feita = any(ex.tipo == "analise" for ex in sessao.execucoes)
    desenho_feito = any(ex.tipo == "desenho" for ex in sessao.execucoes)

    if pagina == "inicio":
        estado_dados = "atual"
    elif dados_feitos:
        estado_dados = "feita"
    else:
        estado_dados = "pendente"

    if pagina == "analise":
        estado_analise = "atual"
    elif analise_feita:
        estado_analise = "feita"
    else:
        estado_analise = "pendente"

    if pagina == "desenho":
        estado_desenho = "atual"
    elif desenho_feito:
        estado_desenho = "feita"
    else:
        estado_desenho = "pendente"

    etapas = [
        {"numero": 1, "nome": "Dados", "href": "/", "estado": estado_dados},
        {"numero": 2, "nome": "Análise", "href": "/analise", "estado": estado_analise},
        {"numero": 3, "nome": "Desenho", "href": "/desenho", "estado": estado_desenho},
    ]

    cfg_carregamento = sessao.cfg_web.get("carregamento", {})
    espera = {
        "analise_ms": int(cfg_carregamento.get("analise_ms", 2600)),
        "desenho_ms": int(cfg_carregamento.get("desenho_ms", 3200)),
        "envio_ms": int(cfg_carregamento.get("envio_ms", 900)),
    }

    return {
        "pagina": pagina,
        "projeto": projeto,
        "ramo": versao[0],
        "commit": versao[1],
        "aviso": AVISO_CURTO,
        "etapas": etapas,
        "espera": espera,
    }


def tela_inicio(sessao: Sessao) -> dict:
    """Monta os dados da tela inicial de carregamento e histórico."""
    conjunto_dict = None
    if sessao.conjunto is not None:
        nome_arq = sessao.conjunto.nome
        ext = Path(nome_arq).suffix.lower()
        texto = sessao.conjunto.conteudo.decode("utf-8-sig", errors="replace")

        nomes_originais: list[str] = []
        if ext == ".csv":
            linhas_nao_vazias = [linha for linha in texto.splitlines() if linha.strip()]
            if linhas_nao_vazias:
                primeira_linha = linhas_nao_vazias[0]
                if primeira_linha.count(";") > primeira_linha.count(","):
                    delim = ";"
                elif "\t" in primeira_linha and "," not in primeira_linha:
                    delim = "\t"
                else:
                    delim = ","

                leitor = csv.reader(io.StringIO(texto), delimiter=delim)
                col_nome = 0
                cabecalho_passou = False
                for partes in leitor:
                    partes = [col.strip() for col in partes]
                    if not partes or not any(partes):
                        continue
                    if not cabecalho_passou:
                        cabecalho_passou = True
                        for idx_c, col in enumerate(partes[:2]):
                            col_lower = col.lower()
                            if col_lower in (
                                "nome", "name", "primer", "nome_primer", "primer_name", "id"
                            ):
                                col_nome = idx_c
                                break
                            elif col_lower in ("sequencia", "sequência", "sequence", "seq"):
                                col_nome = 1 if idx_c == 0 else 0
                        continue
                    if len(partes) > col_nome and partes[col_nome]:
                        nomes_originais.append(partes[col_nome])
        else:
            for linha in texto.splitlines():
                linha_limpa = linha.strip()
                if linha_limpa.startswith(">"):
                    resto = linha_limpa[1:].strip()
                    partes = resto.split()
                    if partes:
                        nomes_originais.append(partes[0])

        reconhecidos: list[dict] = []
        canons_presentes: set[str] = set()
        for nome_orig in nomes_originais:
            reconhecido = reconhecer_nome_primer(nome_orig)
            canon = reconhecido if reconhecido is not None else nome_orig
            canons_presentes.add(canon)
            if canon != nome_orig:
                texto_item = f"{canon} ← {nome_orig}"
            else:
                texto_item = canon
            reconhecidos.append({
                "texto": texto_item,
                "lado": lado_de(canon),
            })

        n_prim = len(nomes_originais)
        tipo_arq = "CSV" if ext == ".csv" else "FASTA"
        resumo = f"{tipo_arq} · {plural(n_prim, 'primer', 'primers')}"

        partes_fip_bip = {"F1c", "F2", "B1c", "B2"}.issubset(canons_presentes)
        tem_fip = "FIP" in canons_presentes
        tem_bip = "BIP" in canons_presentes

        if partes_fip_bip and not tem_fip and not tem_bip:
            resumo += " · FIP e BIP montados a partir das partes"
        elif tem_fip and tem_bip and "F1c" not in canons_presentes and "B1c" not in canons_presentes:
            resumo += " · FIP e BIP inteiros"

        conjunto_dict = {
            "nome_arquivo": nome_arq,
            "resumo": resumo,
            "reconhecidos": reconhecidos,
            "avisos": list(sessao.conjunto.avisos),
        }

    alvo_dict = None
    if sessao.alvo is not None:
        n_reg = max(1, len(getattr(sessao.alvo, "registros", ()) or ()))
        rotulo_regs = plural(n_reg, "registro", "registros")
        alvo_dict = {
            "nome_arquivo": sessao.alvo.nome,
            "resumo": f"FASTA · {rotulo_regs} · {milhar(len(sessao.alvo.sequencia))} pb",
            "cabecalho": sessao.alvo.descricao,
        }

    limite_pb = int(sessao.cfg_web.get("entrada", {}).get("alvo_desenho_max_pb", 2000))
    pode_analisar = sessao.conjunto is not None

    if sessao.alvo is None:
        pode_desenhar = False
        motivo_sem_desenho = "Carregue a sequência-alvo para desenhar."
    elif len(sessao.alvo.sequencia) > limite_pb:
        pode_desenhar = False
        motivo_sem_desenho = (
            f"O alvo tem {milhar(len(sessao.alvo.sequencia))} pb; "
            f"o desenho pela interface aceita até {milhar(limite_pb)} pb."
        )
    else:
        pode_desenhar = True
        motivo_sem_desenho = ""

    tm_como_calibracao = bool(sessao.cfg_web.get("analise", {}).get("tm_como_calibracao", True))
    historico: list[dict] = []
    for ex in reversed(sessao.execucoes):
        if ex.tipo == "analise":
            resumo_v = resumo_violacoes(ex.resultado.violacoes, tm_como_calibracao)
            res_str = resumo_v["texto"]
            if (
                ex.alvo is not None
                and ex.resultado.regioes is not None
                and ex.leitura is not None
                and ex.leitura.conjunto is not None
            ):
                localizadas = len(ex.resultado.regioes)
                total_reg = len(regioes_do_conjunto(ex.leitura.conjunto))
                res_str += f" · {localizadas} de {total_reg} regiões"
            res_tipo = resumo_v["tipo"]
            href = f"/analise/{ex.id}"
            historico.append({
                "tipo": "Análise",
                "conjunto": ex.rotulo_conjunto,
                "alvo": ex.rotulo_alvo,
                "resultado": res_str,
                "resultado_tipo": res_tipo,
                "data": data_hora(ex.criada_em),
                "href": href,
            })
        elif ex.tipo == "desenho":
            k_prop = len(ex.resultado.conjuntos)
            res_str = plural(k_prop, "proposta", "propostas")
            res_tipo = "info" if k_prop > 0 else "calibrar"
            href = f"/desenho/{ex.id}"
            historico.append({
                "tipo": "Desenho",
                "conjunto": "—",
                "alvo": ex.rotulo_alvo,
                "resultado": res_str,
                "resultado_tipo": res_tipo,
                "data": data_hora(ex.criada_em),
                "href": href,
            })

    return {
        "exemplos": sessao.exemplos_disponiveis(),
        "conjunto": conjunto_dict,
        "alvo": alvo_dict,
        "pode_analisar": pode_analisar,
        "pode_desenhar": pode_desenhar,
        "motivo_sem_desenho": motivo_sem_desenho,
        "historico": historico,
    }


def texto_reacao(cond: Condicoes) -> str:
    """Resumo das condições da reação: íons monovalentes presentes, Mg²⁺, dNTP de cada e temperatura."""
    partes = [
        f"{rotulo} {num_curto(valor)} mM"
        for rotulo, valor in (("Na⁺", cond.na_mM), ("K⁺", cond.k_mM), ("NH₄⁺", cond.nh4_mM))
        if valor > 0
    ]
    partes.append(f"Mg²⁺ {num_curto(cond.mg_mM)} mM")
    partes.append(f"dNTP {num_curto(cond.dntp_cada_mM)} mM cada")
    partes.append(f"{num_curto(cond.temp_C)} °C")
    return ", ".join(partes)


def tela_analise(
    execucao: Execucao,
    cfg_primers: dict,
    tm_como_calibracao: bool,
    escala_kcal: list[float],
) -> dict:
    """Monta os dados da tela de análise de propriedades físico-químicas."""
    r = execucao.resultado
    titulo = (
        f"{execucao.rotulo_conjunto} × {execucao.rotulo_alvo}"
        if execucao.alvo is not None
        else execucao.rotulo_conjunto
    )

    reacao_str = texto_reacao(Condicoes.de_config(cfg_primers))

    fita_res = getattr(r, "fita", None)
    if execucao.alvo is not None:
        tam_alvo = len(execucao.alvo.sequencia)
        if fita_res == "-":
            subtitulo = (
                f"{plural(len(r.linhas), 'primer', 'primers')} · "
                f"alvo {execucao.rotulo_alvo}, {milhar(tam_alvo)} pb · fita − · {reacao_str}"
            )
        else:
            subtitulo = (
                f"{plural(len(r.linhas), 'primer', 'primers')} · "
                f"alvo {execucao.rotulo_alvo}, {milhar(tam_alvo)} pb · {reacao_str}"
            )
    else:
        subtitulo = (
            f"{plural(len(r.linhas), 'primer', 'primers')} · "
            f"sem alvo: regiões e distâncias não calculadas · {reacao_str}"
        )

    href_alternar_calibracao = f"/analise/{execucao.id}?tm_calibracao={0 if tm_como_calibracao else 1}"
    rotulo_alternar_calibracao = (
        "Mostrar Tm fora da faixa como violação"
        if tm_como_calibracao
        else "Mostrar Tm fora da faixa como calibração"
    )

    # As posições aproximadas vêm prontas da análise (localizar_detalhado): a página não
    # repete a busca, que em genoma inteiro custava segundos a cada abertura
    melhores_posicoes: dict[str, Posicao] = (
        dict(r.regioes_aproximadas)
        if isinstance(getattr(r, "regioes_aproximadas", None), dict)
        else {}
    )
    bases_ponta = int(cfg_primers.get("casamento", {}).get("bases_ponta_3p", 5))

    regioes_para_dist: dict[str, Regiao] | None
    if r.regioes is None:
        regioes_para_dist = None
    else:
        regioes_para_dist = dict(r.regioes)
        for reg_k, pos_v in melhores_posicoes.items():
            if reg_k not in regioes_para_dist:
                regioes_para_dist[reg_k] = pos_v.regiao

    # KPIs
    total_viol = len(r.violacoes)
    resumo_v = resumo_violacoes(r.violacoes, tm_como_calibracao)

    if execucao.alvo is None:
        val_reg = "—"
        det_reg = "sem alvo"
    elif r.regioes is None:
        val_reg = "—"
        det_reg = "conjunto incompleto"
    else:
        total_reg = (
            len(regioes_do_conjunto(execucao.leitura.conjunto))
            if (execucao.leitura and execucao.leitura.conjunto)
            else 0
        )
        val_reg = f"{len(r.regioes)} de {total_reg}"
        itens_det_reg: list[str] = []
        if execucao.leitura and execucao.leitura.conjunto:
            for reg_chave, _ in regioes_do_conjunto(execucao.leitura.conjunto):
                if reg_chave not in r.regioes and reg_chave in melhores_posicoes:
                    primer_n = REGIAO_PARA_PRIMER.get(reg_chave, reg_chave)
                    erros_reg = erros_de_base(reg_chave, melhores_posicoes[reg_chave])
                    n_ponta = sum(1 for e in erros_reg if e.distancia_ponta < bases_ponta)
                    item_reg = f"{primer_n} com {len(erros_reg)} dif."
                    if n_ponta:
                        item_reg += f", {n_ponta} na ponta"
                    itens_det_reg.append(item_reg)
        det_reg = ", ".join(itens_det_reg)

    if execucao.alvo is None:
        val_dist = "—"
        det_dist = "sem alvo"
    elif regioes_para_dist is None:
        val_dist = "—"
        det_dist = "conjunto incompleto"
    else:
        dist_calc = calcular_distancias(regioes_para_dist)
        cfg_dist = cfg_primers.get("distancias", {})
        dentro_cnt = 0
        fora_nomes: list[str] = []
        for k_dist in NOMES_DISTANCIA:
            if k_dist in dist_calc and k_dist in cfg_dist:
                v_d = dist_calc[k_dist]
                d_min, d_max = cfg_dist[k_dist]
                if d_min <= v_d <= d_max:
                    dentro_cnt += 1
                else:
                    fora_nomes.append(f"{NOMES_DISTANCIA[k_dist].replace(' ', '')} fora")
        val_dist = f"{dentro_cnt} de {len(dist_calc)}"
        det_dist = ", ".join(fora_nomes)

    kpis = [
        {"rotulo": "Primers lidos", "valor": str(len(r.linhas)), "detalhe": "", "tipo": "neutro"},
        {
            "rotulo": "Violações",
            "valor": str(resumo_v["n_violacoes"]),
            "detalhe": resumo_v["detalhe"],
            "tipo": resumo_v["tipo"],
        },
        {"rotulo": "Regiões localizadas", "valor": val_reg, "detalhe": det_reg, "tipo": "neutro"},
        {"rotulo": "Distâncias", "valor": val_dist, "detalhe": det_dist, "tipo": "neutro"},
    ]

    # Aviso de calibração
    aviso_calibracao = ""
    viols_tm = [v for v in r.violacoes if v.regra == "tm_c"]
    if tm_como_calibracao and viols_tm:
        primers_tm = {v.item for v in viols_tm}
        n_p = len(primers_tm)
        x_str = "1 primer fica" if n_p == 1 else f"{n_p} primers ficam"
        maior_dist = 0.0
        for v in viols_tm:
            chave_tm = FAIXA_TM_PRIMER.get(v.item, "F3_B3")
            faixa_tm = cfg_primers["tm_c"][chave_tm]
            tm_val = float(v.valor) if v.valor is not None else 0.0
            if tm_val < faixa_tm[0]:
                dist_borda = faixa_tm[0] - tm_val
            elif tm_val > faixa_tm[1]:
                dist_borda = tm_val - faixa_tm[1]
            else:
                dist_borda = 0.0
            if dist_borda > maior_dist:
                maior_dist = dist_borda
        aviso_calibracao = (
            f"As faixas de Tm ainda não estão calibradas. {x_str} fora da faixa de Tm, "
            f"a até {num(maior_dist, 1)} °C dela. "
            f"Leia as violações de Tm como comparação, não como reprovação."
        )

    # Mapa
    mapa = None
    if execucao.alvo is not None and r.regioes is not None:
        tam_alvo = len(execucao.alvo.sequencia)
        segmentos: list[dict] = []
        for reg_k, reg_v in r.regioes.items():
            segmentos.append({
                "nome": reg_k,
                "regiao": reg_k,
                "ini": reg_v.ini,
                "fim": reg_v.fim,
                "esquerda": pct(reg_v.ini, tam_alvo),
                "largura": pct(reg_v.fim - reg_v.ini, tam_alvo),
                "fantasma": False,
            })
        for reg_k, pos_v in melhores_posicoes.items():
            segmentos.append({
                "nome": reg_k,
                "regiao": reg_k,
                "ini": pos_v.regiao.ini,
                "fim": pos_v.regiao.fim,
                "esquerda": pct(pos_v.regiao.ini, tam_alvo),
                "largura": pct(pos_v.regiao.fim - pos_v.regiao.ini, tam_alvo),
                "fantasma": True,
            })
        segmentos.sort(key=lambda s: (s["ini"], s["fim"]))
        ticks = marcas_de_eixo(tam_alvo)

        faixas_fora: list[dict] = []
        dist_calc = calcular_distancias(regioes_para_dist) if regioes_para_dist else {}
        cfg_dist = cfg_primers.get("distancias", {})
        for k_d, nome_fmt in NOMES_DISTANCIA.items():
            if k_d in dist_calc and k_d in cfg_dist:
                medida = dist_calc[k_d]
                d_min, d_max = cfg_dist[k_d]
                if medida < d_min or medida > d_max:
                    a: int | None = None
                    b: int | None = None
                    if regioes_para_dist is not None:
                        if k_d == "F2_B2" and "F2" in regioes_para_dist and "B2c" in regioes_para_dist:
                            a, b = regioes_para_dist["F2"].ini, regioes_para_dist["B2c"].fim
                        elif k_d == "F2_F1" and "F2" in regioes_para_dist and "F1" in regioes_para_dist:
                            a, b = regioes_para_dist["F2"].ini, regioes_para_dist["F1"].ini
                        elif k_d == "B1c_B2c" and "B1c" in regioes_para_dist and "B2c" in regioes_para_dist:
                            a, b = regioes_para_dist["B1c"].fim, regioes_para_dist["B2c"].fim
                        elif k_d == "F3_F2_gap" and "F3" in regioes_para_dist and "F2" in regioes_para_dist:
                            a, b = regioes_para_dist["F3"].fim, regioes_para_dist["F2"].ini
                        elif (
                            k_d == "B2c_B3c_gap" and "B2c" in regioes_para_dist and "B3c" in regioes_para_dist
                        ):
                            a, b = regioes_para_dist["B2c"].fim, regioes_para_dist["B3c"].ini
                        elif k_d == "F1_B1c_gap" and "F1" in regioes_para_dist and "B1c" in regioes_para_dist:
                            a, b = regioes_para_dist["F1"].fim, regioes_para_dist["B1c"].ini
                        elif (
                            k_d == "F3_B3_total" and "F3" in regioes_para_dist and "B3c" in regioes_para_dist
                        ):
                            a, b = regioes_para_dist["F3"].ini, regioes_para_dist["B3c"].fim

                    if a is not None and b is not None:
                        lim_txt = (
                            f"máximo {num_curto(d_max)}" if medida > d_max else f"mínimo {num_curto(d_min)}"
                        )
                        texto_fora = f"{nome_fmt}: {milhar(medida)} nt · {lim_txt}"
                        faixas_fora.append({
                            "texto": texto_fora,
                            "esquerda": pct(min(a, b), tam_alvo),
                            "largura": pct(abs(b - a), tam_alvo),
                        })

        sinal_fita = texto_fita(fita_res)
        partes_desc: list[str] = []
        for s in segmentos:
            item_desc = f"{s['nome']} em {coord(s['ini'], s['fim'])}"
            if s["fantasma"]:
                item_desc += " (melhor posição, não exata)"
            partes_desc.append(item_desc)
        desc_mapa = (
            f"Mapa linear do alvo de {milhar(tam_alvo)} pb, fita {sinal_fita}, coordenadas a partir de 1, "
            "com as regiões " + ", ".join(partes_desc) + "."
        )
        if fita_res == "-":
            legenda_fita = "fita − (complemento reverso do alvo)"
        else:
            legenda_fita = "fita +"
        legenda_mapa = (
            f"{legenda_fita} · {milhar(tam_alvo)} pb · coordenadas a partir de 1 · "
            "F1, B2c, B3c e LF são lidas pelo complemento reverso do primer"
        )

        mapa = {
            "comprimento": tam_alvo,
            "fita": sinal_fita,
            "legenda": legenda_mapa,
            "segmentos": segmentos,
            "ticks": ticks,
            "faixas_fora": faixas_fora,
            "descricao": desc_mapa,
        }

    # Distâncias
    distancias_lista: list[dict] = []
    if regioes_para_dist is not None:
        dist_calc = calcular_distancias(regioes_para_dist)
        cfg_dist = cfg_primers.get("distancias", {})
        for k_dist, nome_dist in NOMES_DISTANCIA.items():
            if k_dist in dist_calc and k_dist in cfg_dist:
                d = dist_calc[k_dist]
                d_min, d_max = cfg_dist[k_dist]
                dentro = d_min <= d <= d_max
                if dentro:
                    txt_dist = "dentro"
                elif d > d_max:
                    txt_dist = f"{d - d_max} acima do máximo"
                else:
                    txt_dist = f"{d_min - d} abaixo do mínimo"
                distancias_lista.append({
                    "nome": nome_dist,
                    "medida": milhar(d),
                    "faixa": faixa(d_min, d_max),
                    "dentro": dentro,
                    "texto": txt_dist,
                    "regua": regua_distancia(d, d_min, d_max),
                })

    # Primers
    ligante = (
        execucao.leitura.conjunto.ligante
        if (execucao.leitura and execucao.leitura.conjunto)
        else ""
    )
    primers_linhas: list[dict] = []
    for linha in r.linhas:
        nome = linha["nome"]
        if (
            nome == "FIP"
            and execucao.leitura
            and "F1c" in execucao.leitura.primers
            and "F2" in execucao.leitura.primers
        ):
            seq = execucao.leitura.primers["F1c"]
            seq2 = ligante + execucao.leitura.primers["F2"]
        elif (
            nome == "BIP"
            and execucao.leitura
            and "B1c" in execucao.leitura.primers
            and "B2" in execucao.leitura.primers
        ):
            seq = execucao.leitura.primers["B1c"]
            seq2 = ligante + execucao.leitura.primers["B2"]
        else:
            seq = linha["sequencia"]
            seq2 = ""

        if nome in ("FIP", "BIP") or execucao.alvo is None or r.regioes is None:
            regiao_str = "—"
        else:
            reg_k = PRIMER_PARA_REGIAO.get(nome)
            if reg_k and reg_k in r.regioes:
                regiao_str = coord(r.regioes[reg_k].ini, r.regioes[reg_k].fim)
            else:
                regiao_str = "não localizado"

        regiao_viol = regiao_str == "não localizado"
        nt_viol = any(v.item == nome and v.regra == "comprimento" for v in r.violacoes)
        gc_viol = any(v.item == nome and v.regra == "gc_pct" for v in r.violacoes)
        tm_viol = any(v.item == nome and v.regra == "tm_c" for v in r.violacoes)

        if tm_viol:
            tm_estado = "calibrar" if tm_como_calibracao else "viol"
        else:
            tm_estado = "ok"

        ext_val = linha.get("extremidade")
        if ext_val == "3'":
            ponta = "3′"
        elif ext_val == "5'":
            ponta = "5′"
        else:
            ponta = ""

        if linha.get("grampo_tm_c") is None:
            grampo_str = "—"
        else:
            grampo_str = f"{num(linha['grampo_tm_c'], 1)} / {num(linha['grampo_dg_kcal'])}"

        partes_presentes = (
            ("F1c" in execucao.leitura.primers and "F2" in execucao.leitura.primers)
            if (nome == "FIP" and execucao.leitura)
            else (
                ("B1c" in execucao.leitura.primers and "B2" in execucao.leitura.primers)
                if (nome == "BIP" and execucao.leitura)
                else False
            )
        )

        if nome in ("FIP", "BIP") and partes_presentes:
            situacao = "só métricas"
            situacao_tipo = "neutro"
        else:
            if nome in ("FIP", "BIP"):
                itens_associados = ("FIP", "F1c", "F2") if nome == "FIP" else ("BIP", "B1c", "B2")
                viols_do_primer = [v for v in r.violacoes if v.item in itens_associados]
                padrao_sem_viols = ("só métricas", "neutro")
            else:
                reg_k = PRIMER_PARA_REGIAO.get(nome)
                viols_do_primer = [
                    v
                    for v in r.violacoes
                    if v.item == nome or (v.regra == "nao_encontrado" and v.item == reg_k)
                ]
                padrao_sem_viols = ("ok", "ok")

            if not viols_do_primer:
                situacao, situacao_tipo = padrao_sem_viols
            else:
                tem_tm = any(v.regra == "tm_c" for v in viols_do_primer)
                outras = sum(1 for v in viols_do_primer if v.regra != "tm_c")
                if tm_como_calibracao:
                    if outras == 0 and not tem_tm:
                        situacao = "ok"
                        situacao_tipo = "ok"
                    elif outras == 0 and tem_tm:
                        situacao = "Tm a calibrar"
                        situacao_tipo = "calibrar"
                    elif outras > 0 and tem_tm:
                        situacao = f"{outras} + Tm"
                        situacao_tipo = "viol"
                    else:
                        situacao = plural(outras, "violação", "violações")
                        situacao_tipo = "viol"
                else:
                    total_p = outras + (1 if tem_tm else 0)
                    if total_p == 0:
                        situacao = "ok"
                        situacao_tipo = "ok"
                    else:
                        situacao = plural(total_p, "violação", "violações")
                        situacao_tipo = "viol"

        primers_linhas.append({
            "nome": nome,
            "lado": lado_de(nome),
            "seq": seq,
            "seq2": seq2,
            "regiao": regiao_str,
            "regiao_viol": regiao_viol,
            "nt": str(linha["comprimento"]),
            "nt_viol": nt_viol,
            "gc": num(linha["gc_pct"]),
            "gc_viol": gc_viol,
            "tm": num(linha["tm_c"]),
            "tm_estado": tm_estado,
            "tm_reacao": num(linha.get("tm_reacao_c")),
            "dg": num(linha["dg_extremidade_kcal"]),
            "ponta": ponta,
            "grampo": grampo_str,
            "homodimero": num(linha["homodimero_dg_kcal"]),
            "situacao": situacao,
            "situacao_tipo": situacao_tipo,
        })

    tem_fip_bip = any(linha["nome"] in ("FIP", "BIP") for linha in r.linhas)
    if tem_fip_bip:
        lig_desc = f"com o ligante {ligante}" if ligante else "sem ligante"
        nota_fip_bip = (
            f"FIP e BIP são montados das partes, {lig_desc}. "
            f"Tm e ΔG de ponta não se aplicam a um primer quimérico."
        )
    else:
        nota_fip_bip = ""

    # Grupos de violações
    itens_outras: list[dict] = []
    itens_tm: list[dict] = []

    for v in r.violacoes:
        if v.regra == "tm_c":
            chave_tm = FAIXA_TM_PRIMER.get(v.item, "F3_B3")
            faixa_tm = cfg_primers["tm_c"][chave_tm]
            itens_tm.append({
                "item": v.item,
                "lado": lado_de(v.item),
                "texto": f"{num(v.valor)} °C",
                "faixa": f"faixa {faixa(faixa_tm[0], faixa_tm[1])}",
            })
        elif v.regra == "comprimento":
            if v.item in ("FIP", "BIP") and getattr(v, "faixa", None):
                faixa_txt = f"faixa {v.faixa.strip('[]').replace(', ', '–')}"
            elif v.item in ("LF", "LB"):
                c_min = cfg_primers["comprimento"]["loop_min"]
                c_max = cfg_primers["comprimento"]["loop_max"]
                faixa_txt = f"faixa {faixa(c_min, c_max)}"
            else:
                c_min = cfg_primers["comprimento"]["min"]
                c_max = cfg_primers["comprimento"]["max"]
                faixa_txt = f"faixa {faixa(c_min, c_max)}"
            itens_outras.append({
                "item": v.item,
                "lado": lado_de(v.item),
                "texto": f"Comprimento {int(v.valor)} nt",
                "faixa": faixa_txt,
            })
        elif v.regra == "gc_pct":
            gc_min = cfg_primers["gc_pct"]["min"]
            gc_max = cfg_primers["gc_pct"]["max"]
            itens_outras.append({
                "item": v.item,
                "lado": lado_de(v.item),
                "texto": f"GC {num(v.valor)} %",
                "faixa": f"faixa {faixa(gc_min, gc_max)}",
            })
        elif v.regra == "extremidade":
            ponta_txt = "5′" if v.item in ("F1c", "B1c") else "3′"
            dg_max = cfg_primers["extremidade"]["dg_max"]
            itens_outras.append({
                "item": v.item,
                "lado": lado_de(v.item),
                "texto": f"ΔG da ponta {ponta_txt} {num(v.valor)} kcal/mol",
                "faixa": f"máximo {num_curto(dg_max)}",
            })
        elif v.regra == "nao_encontrado":
            primer_item = REGIAO_PARA_PRIMER.get(v.item, v.item)
            txt_item = "Não localizado"
            if v.item in melhores_posicoes:
                mp = melhores_posicoes[v.item]
                erros_v = erros_de_base(v.item, mp)
                det_difs = " | ".join(texto_erro_de_base(e, bases_ponta) for e in erros_v)
                txt_item += (
                    f": melhor posição {coord(mp.regiao.ini, mp.regiao.fim)} (fita {texto_fita(mp.fita)}), "
                    f"{plural(len(erros_v), 'diferença', 'diferenças')} ({det_difs})"
                )
            itens_outras.append({
                "item": primer_item,
                "lado": lado_de(primer_item),
                "texto": txt_item,
                "faixa": "",
            })
        elif v.regra == "distancia":
            nome_sem_espacos = NOMES_DISTANCIA.get(v.item, v.item).replace(" ", "")
            d_min, d_max = cfg_primers["distancias"][v.item]
            itens_outras.append({
                "item": nome_sem_espacos,
                "lado": "X",
                "texto": f"Distância {milhar(int(v.valor))} nt",
                "faixa": f"faixa {faixa(d_min, d_max)}",
            })
        elif v.regra == "loop_fora_da_regiao":
            itens_outras.append({
                "item": v.item,
                "lado": "L",
                "texto": "Fora da região do loop",
                "faixa": f"região {v.faixa}",
            })

    grupos_violacoes: list[dict] = []
    if itens_outras:
        grupos_violacoes.append({
            "titulo": f"Composição e alvo · {len(itens_outras)}",
            "tipo": "viol",
            "itens": itens_outras,
        })
    if itens_tm:
        if tm_como_calibracao:
            titulo_tm = f"Tm — faixa a calibrar · {len(itens_tm)}"
            tipo_tm = "calibrar"
        else:
            titulo_tm = f"Tm fora da faixa · {len(itens_tm)}"
            tipo_tm = "viol"
        grupos_violacoes.append({
            "titulo": titulo_tm,
            "tipo": tipo_tm,
            "itens": itens_tm,
        })

    # Heterodímeros
    nomes_het = r.heterodimeros.get("nomes", [])
    matriz_het = r.heterodimeros.get("dg_kcal", [])
    e0, e1, e2, e3 = escala_kcal[0], escala_kcal[1], escala_kcal[2], escala_kcal[3]

    def faixa_het(v_kcal: float) -> int:
        if v_kcal <= e0:
            return 0
        if v_kcal <= e1:
            return 1
        if v_kcal <= e2:
            return 2
        if v_kcal <= e3:
            return 3
        return 4

    linhas_het: list[dict] = []
    for i, nome_a in enumerate(nomes_het):
        celulas = []
        for j, _ in enumerate(nomes_het):
            v_val = matriz_het[i][j]
            celulas.append({
                "valor": num(v_val),
                "faixa": faixa_het(v_val),
                "diagonal": i == j,
            })
        linhas_het.append({"nome": nome_a, "celulas": celulas})

    legenda_het = [
        {"texto": f"> {num_curto(e3)}", "faixa": 4},
        {"texto": f"{num_curto(e2)} a {num_curto(e3)}", "faixa": 3},
        {"texto": f"{num_curto(e1)} a {num_curto(e2)}", "faixa": 2},
        {"texto": f"{num_curto(e0)} a {num_curto(e1)}", "faixa": 1},
        {"texto": f"≤ {num_curto(e0)}", "faixa": 0},
    ]

    resumo_het = ""
    if nomes_het and matriz_het:
        min_geral = None
        pos_geral = None
        for i in range(len(nomes_het)):
            for j in range(len(nomes_het)):
                val = matriz_het[i][j]
                if min_geral is None or val < min_geral:
                    min_geral = val
                    pos_geral = (i, j)

        if pos_geral is not None and min_geral is not None:
            i_min, j_min = pos_geral
            if i_min == j_min:
                a_nome = nomes_het[i_min]
                resumo_het = f"Par mais estável: {a_nome} × {a_nome}, {num(min_geral)} (homodímero)."
                min_dist = None
                pos_dist = None
                for i in range(len(nomes_het)):
                    for j in range(i + 1, len(nomes_het)):
                        val_d = matriz_het[i][j]
                        if min_dist is None or val_d < min_dist:
                            min_dist = val_d
                            pos_dist = (i, j)
                if pos_dist is not None and min_dist is not None:
                    a_d = nomes_het[pos_dist[0]]
                    b_d = nomes_het[pos_dist[1]]
                    resumo_het += f" Entre primers distintos: {a_d} × {b_d}, {num(min_dist)}."
            else:
                i_ord = min(i_min, j_min)
                j_ord = max(i_min, j_min)
                a_nome = nomes_het[i_ord]
                b_nome = nomes_het[j_ord]
                resumo_het = f"Par mais estável: {a_nome} × {b_nome}, {num(min_geral)}."

    return {
        "id": execucao.id,
        "titulo": titulo,
        "subtitulo": subtitulo,
        "tm_como_calibracao": tm_como_calibracao,
        "href_alternar_calibracao": href_alternar_calibracao,
        "rotulo_alternar_calibracao": rotulo_alternar_calibracao,
        "pode_desenhar": execucao.alvo is not None,
        "url_tubos": f"/tubos?execucao={execucao.id}" if getattr(r, "incompat", None) is not None else "",
        "kpis": kpis,
        "aviso_calibracao": aviso_calibracao,
        "mapa": mapa,
        "distancias": distancias_lista,
        "primers": primers_linhas,
        "nota_fip_bip": nota_fip_bip,
        "total_violacoes": total_viol,
        "grupos_violacoes": grupos_violacoes,
        "heterodimeros": {
            "nomes": nomes_het,
            "linhas": linhas_het,
            "legenda": legenda_het,
            "resumo": resumo_het,
        },
        "avisos": list(r.avisos),
    }


def tela_desenho(
    execucao: Execucao,
    cfg_primers: dict,
    sel: int,
    agrupar: bool,
    tm_a_calibrar: bool = False,
) -> dict:
    """Monta os dados da tela de resultados do desenho de primers.

    tm_a_calibrar ([analise] tm_como_calibracao do web.toml): avisa que as faixas de Tm,
    que no desenho filtram os candidatos, ainda estão a calibrar.
    """
    rotulo_alvo = execucao.rotulo_alvo
    tam_alvo = len(execucao.alvo.sequencia) if execucao.alvo else 0
    cabecalho = execucao.alvo.descricao if execucao.alvo else ""

    res_desenho: ResultadoDesenho = execucao.resultado  # type: ignore[assignment]
    estatisticas = res_desenho.estatisticas
    cands_tipo = estatisticas.get("candidatos_por_tipo", {})

    candidatos: list[dict] = []
    for tipo in ("F3", "F2", "F1", "B1c", "B2c", "B3c"):
        n_c = cands_tipo.get(tipo, 0)
        lado = "F" if tipo in ("F3", "F2", "F1") else "B"
        candidatos.append({"texto": f"{tipo} {n_c}", "lado": lado})

    examinados = milhar(estatisticas.get("conjuntos_examinados", 0))
    relaxado = bool(estatisticas.get("relaxado", False))

    conjuntos = res_desenho.conjuntos
    vazio = len(conjuntos) == 0
    mensagem_vazio = ""
    if vazio:
        mensagem_vazio = str(estatisticas.get("diagnostico", "")) or (
            "Nenhuma proposta passou nas regras do config para este alvo. "
            "Veja em Parâmetros as faixas usadas."
        )
    pode_relaxar = vazio and not relaxado

    grupos: list[list[ConjuntoDesenhado]] = []
    mapa_grupos: dict[tuple, list[ConjuntoDesenhado]] = {}
    for c in conjuntos:
        chave = (
            num(c.pontuacao),
            c.regioes.get("F3"),
            c.regioes.get("F2"),
            c.regioes.get("B2c"),
            c.regioes.get("B3c"),
        )
        if chave not in mapa_grupos:
            lista_grupo: list[ConjuntoDesenhado] = []
            mapa_grupos[chave] = lista_grupo
            grupos.append(lista_grupo)
        mapa_grupos[chave].append(c)

    posicoes_existentes = {c.posicao for c in conjuntos}
    if sel not in posicoes_existentes:
        sel = conjuntos[0].posicao if conjuntos else 1

    if agrupar:
        for grp in grupos:
            if any(c.posicao == sel for c in grp):
                sel = grp[0].posicao
                break

    a_param = 1 if agrupar else 0
    href_alternar_agrupar = f"/desenho/{execucao.id}?sel=1&agrupar={1 - a_param}"

    prop1 = next((c for c in conjuntos if c.posicao == 1), conjuntos[0] if conjuntos else None)

    def comp_regiao(reg_nome: str, prop: ConjuntoDesenhado) -> tuple[str, bool]:
        r = prop.regioes.get(reg_nome)
        r1 = prop1.regioes.get(reg_nome) if prop1 else None
        if r is None:
            return "—", (r1 is not None)
        if r1 is None:
            return coord(r.ini, r.fim), True
        if r == r1:
            return coord(r.ini, r.fim), False
        d = r.ini - r1.ini if r.ini != r1.ini else r.fim - r1.fim
        sinal = "+" if d >= 0 else "−"
        return f"{coord(r.ini, r.fim)} ({sinal}{abs(d)})", True

    def formatar_avisos(avisos_orig: list[str]) -> str:
        avs = []
        for av in avisos_orig:
            if av == "sem LF candidato":
                avs.append("sem LF")
            elif av == "sem LB candidato":
                avs.append("sem LB")
            else:
                avs.append(av)
        return ", ".join(avs)

    lista_itens: list[dict] = []
    itens_exibidos = [grp[0] for grp in grupos] if agrupar else conjuntos

    for c in itens_exibidos:
        if agrupar:
            chave_c = (
                num(c.pontuacao),
                c.regioes.get("F3"),
                c.regioes.get("F2"),
                c.regioes.get("B2c"),
                c.regioes.get("B3c"),
            )
            tam_grp = len(mapa_grupos.get(chave_c, []))
            extra = f" + {plural(tam_grp - 1, 'variante', 'variantes')}" if tam_grp > 1 else ""
        else:
            extra = ""

        f1_txt, f1_dif = comp_regiao("F1", c)
        b1c_txt, b1c_dif = comp_regiao("B1c", c)
        lb_txt, lb_dif = comp_regiao("LB", c)

        lista_itens.append({
            "num": str(c.posicao),
            "extra": extra,
            "selecionada": c.posicao == sel,
            "href": f"/desenho/{execucao.id}?sel={c.posicao}&agrupar={a_param}",
            "pontuacao": num(c.pontuacao),
            "desvio_tm": f"{num(c.desvio_tm)} °C",
            "pior_dimero": num(c.pior_heterodimero_kcal),
            "f1": f1_txt,
            "f1_dif": f1_dif,
            "b1c": b1c_txt,
            "b1c_dif": b1c_dif,
            "lb": lb_txt,
            "lb_dif": lb_dif,
            "avisos": formatar_avisos(c.avisos),
        })

    nota_empate = ""
    if len(conjuntos) >= 2:
        p0 = conjuntos[0]
        pont0 = num(p0.pontuacao)
        iguais_base = all(
            num(c.pontuacao) == pont0
            and c.regioes.get("F3") == p0.regioes.get("F3")
            and c.regioes.get("F2") == p0.regioes.get("F2")
            and c.regioes.get("B2c") == p0.regioes.get("B2c")
            and c.regioes.get("B3c") == p0.regioes.get("B3c")
            for c in conjuntos
        )
        if iguais_base:
            diferentes: list[str] = []
            for reg_k in ("F1", "B1c", "LF", "LB"):
                primeira = p0.regioes.get(reg_k)
                if any(c.regioes.get(reg_k) != primeira for c in conjuntos):
                    diferentes.append(reg_k)
            nota_empate = (
                f"As {len(conjuntos)} propostas empatam em {pont0} e compartilham F3, F2, B2c e B3c — "
                f"diferem só por deslocamentos nas bordas de {juntar(diferentes)}. "
                f"Na prática, é um desenho só."
            )

    proposta_dict = None
    if not vazio:
        prop_sel = next(c for c in conjuntos if c.posicao == sel)
        c_sel = prop_sel.conjunto
        ligante = c_sel.ligante
        regioes = prop_sel.regioes

        segs_prop: list[dict] = []
        for reg_k in ("F3", "F2", "F1", "B1c", "B2c", "B3c"):
            if reg_k in regioes:
                r_obj = regioes[reg_k]
                segs_prop.append({
                    "nome": reg_k,
                    "regiao": reg_k,
                    "ini": r_obj.ini,
                    "fim": r_obj.fim,
                    "esquerda": pct(r_obj.ini, tam_alvo),
                    "largura": pct(r_obj.fim - r_obj.ini, tam_alvo),
                    "fantasma": False,
                })
        if "LF" in regioes:
            r_lf = regioes["LF"]
            segs_prop.append({
                "nome": "LF",
                "regiao": "LF",
                "ini": r_lf.ini,
                "fim": r_lf.fim,
                "esquerda": pct(r_lf.ini, tam_alvo),
                "largura": pct(r_lf.fim - r_lf.ini, tam_alvo),
                "fantasma": False,
            })
        elif "F2" in regioes and "F1" in regioes and regioes["F2"].fim < regioes["F1"].ini:
            ini_f = regioes["F2"].fim
            fim_f = regioes["F1"].ini
            segs_prop.append({
                "nome": "sem LF",
                "regiao": "LF",
                "ini": ini_f,
                "fim": fim_f,
                "esquerda": pct(ini_f, tam_alvo),
                "largura": pct(fim_f - ini_f, tam_alvo),
                "fantasma": True,
            })

        if "LB" in regioes:
            r_lb = regioes["LB"]
            segs_prop.append({
                "nome": "LB",
                "regiao": "LB",
                "ini": r_lb.ini,
                "fim": r_lb.fim,
                "esquerda": pct(r_lb.ini, tam_alvo),
                "largura": pct(r_lb.fim - r_lb.ini, tam_alvo),
                "fantasma": False,
            })
        elif "B1c" in regioes and "B2c" in regioes and regioes["B1c"].fim < regioes["B2c"].ini:
            ini_b = regioes["B1c"].fim
            fim_b = regioes["B2c"].ini
            segs_prop.append({
                "nome": "sem LB",
                "regiao": "LB",
                "ini": ini_b,
                "fim": fim_b,
                "esquerda": pct(ini_b, tam_alvo),
                "largura": pct(fim_b - ini_b, tam_alvo),
                "fantasma": True,
            })

        segs_prop.sort(key=lambda s: (s["ini"], s["fim"]))
        ticks_prop = marcas_de_eixo(tam_alvo)

        desc_segs = ", ".join(f"{s['nome']} {coord(s['ini'], s['fim'])}" for s in segs_prop)
        desc_prop = f"Mapa do alvo de {tam_alvo} pb da proposta {sel}: {desc_segs}."

        primers_prop = [
            {
                "nome": "F3",
                "lado": "F",
                "seq": c_sel.F3,
                "seq2": "",
                "regiao": coord(regioes["F3"].ini, regioes["F3"].fim),
                "nt": str(len(c_sel.F3)),
                "copiar": c_sel.F3,
            },
            {
                "nome": "B3",
                "lado": "B",
                "seq": c_sel.B3,
                "seq2": "",
                "regiao": f"B3c {coord(regioes['B3c'].ini, regioes['B3c'].fim)}",
                "nt": str(len(c_sel.B3)),
                "copiar": c_sel.B3,
            },
            {
                "nome": "FIP",
                "lado": "X",
                "seq": c_sel.F1c,
                "seq2": ligante + c_sel.F2,
                "regiao": (
                    f"F1 {coord(regioes['F1'].ini, regioes['F1'].fim)} · "
                    f"F2 {coord(regioes['F2'].ini, regioes['F2'].fim)}"
                ),
                "nt": str(len(c_sel.FIP)),
                "copiar": c_sel.FIP,
            },
            {
                "nome": "BIP",
                "lado": "X",
                "seq": c_sel.B1c,
                "seq2": ligante + c_sel.B2,
                "regiao": (
                    f"B1c {coord(regioes['B1c'].ini, regioes['B1c'].fim)} · "
                    f"B2c {coord(regioes['B2c'].ini, regioes['B2c'].fim)}"
                ),
                "nt": str(len(c_sel.BIP)),
                "copiar": c_sel.BIP,
            },
        ]

        if c_sel.LF is not None and "LF" in regioes:
            primers_prop.append({
                "nome": "LF",
                "lado": "L",
                "seq": c_sel.LF,
                "seq2": "",
                "regiao": coord(regioes["LF"].ini, regioes["LF"].fim),
                "nt": str(len(c_sel.LF)),
                "copiar": c_sel.LF,
            })
        else:
            primers_prop.append({
                "nome": "LF",
                "lado": "L",
                "seq": "—",
                "seq2": "",
                "regiao": "sem candidato",
                "nt": "—",
                "copiar": "",
            })

        if c_sel.LB is not None and "LB" in regioes:
            primers_prop.append({
                "nome": "LB",
                "lado": "L",
                "seq": c_sel.LB,
                "seq2": "",
                "regiao": coord(regioes["LB"].ini, regioes["LB"].fim),
                "nt": str(len(c_sel.LB)),
                "copiar": c_sel.LB,
            })
        else:
            primers_prop.append({
                "nome": "LB",
                "lado": "L",
                "seq": "—",
                "seq2": "",
                "regiao": "sem candidato",
                "nt": "—",
                "copiar": "",
            })

        partes_nota: list[str] = []
        loop_min = cfg_primers["comprimento"]["loop_min"]
        loop_max = cfg_primers["comprimento"]["loop_max"]

        if c_sel.LF is None and "F2" in regioes and "F1" in regioes:
            partes_nota.append(
                f"Sem LF: nenhuma janela de {faixa(loop_min, loop_max)} nt entre F2 e F1 "
                f"(posições {coord(regioes['F2'].fim, regioes['F1'].ini)}) passou nas regras de loop."
            )
        if c_sel.LB is None and "B1c" in regioes and "B2c" in regioes:
            partes_nota.append(
                f"Sem LB: nenhuma janela de {faixa(loop_min, loop_max)} nt entre B1c e B2c "
                f"(posições {coord(regioes['B1c'].fim, regioes['B2c'].ini)}) passou nas regras de loop."
            )
        lig_desc = f"com o ligante {ligante}" if ligante else "sem ligante"
        partes_nota.append(f"FIP = F1c + F2 e BIP = B1c + B2, {lig_desc}.")

        proposta_dict = {
            "num": str(sel),
            "mapa": {
                "comprimento": tam_alvo,
                "segmentos": segs_prop,
                "ticks": ticks_prop,
                "descricao": desc_prop,
            },
            "primers": primers_prop,
            "nota": " ".join(partes_nota),
        }

    return {
        "id": execucao.id,
        "alvo": rotulo_alvo,
        "comprimento": tam_alvo,
        "cabecalho": cabecalho,
        "aviso": AVISO,
        "candidatos": candidatos,
        "examinados": examinados,
        "relaxado": relaxado,
        "pode_relaxar": pode_relaxar,
        "agrupar": agrupar,
        "sel": sel,
        "href_alternar_agrupar": href_alternar_agrupar,
        "nota_empate": nota_empate,
        "nota_tm": (
            f"As faixas de Tm de referência ({texto_faixas_tm(cfg_primers)}) filtram os candidatos "
            "do desenho e ainda estão a calibrar."
            if tm_a_calibrar
            else ""
        ),
        "vazio": vazio,
        "mensagem_vazio": mensagem_vazio,
        "lista": lista_itens,
        "proposta": proposta_dict,
    }


def tela_parametros(sessao: Sessao, cfg_primers: dict) -> dict:
    """Monta os dados da tela de parâmetros: condições, Tm de referência e Tm da reação."""
    ref_exec: Execucao | None = None
    referencia = None
    for ex in reversed(sessao.execucoes):
        if ex.tipo == "analise" and any(linha.get("tm_c") is not None for linha in ex.resultado.linhas):
            ref_exec = ex
            rotulo = (
                f"{ex.rotulo_conjunto} × {ex.rotulo_alvo}"
                if ex.alvo is not None
                else ex.rotulo_conjunto
            )
            referencia = {"rotulo": rotulo, "id": ex.id}
            break

    cond = Condicoes.de_config(cfg_primers)
    conv = CondicoesReferencia.de_config(cfg_primers)
    tm_a_calibrar = bool(sessao.cfg_web.get("analise", {}).get("tm_como_calibracao", True))
    marca_tm = "a calibrar" if tm_a_calibrar else ""

    calibracao = None

    if ref_exec is not None:
        cats_def = [
            ("Externos", ("F3", "B3"), "F3_B3", "externos"),
            ("Internos, ponta 3′", ("F2", "B2"), "F2_B2", "internos"),
            ("Internos, ponta 5′", ("F1c", "B1c"), "F1c_B1c", "internos"),
            ("Loops", ("LF", "LB"), "loops", "loops"),
        ]

        linhas_por_nome = {
            linha["nome"]: linha for linha in ref_exec.resultado.linhas if linha.get("tm_c") is not None
        }
        cats_ativas: list[dict] = []
        todos_valores: list[float] = []
        primers_medidos: list[dict] = []

        for cat_nome, prims_nomes, tm_key, conc_key in cats_def:
            prims_cat: list[dict] = []
            conc_nM = float(cfg_primers["concentracao_nM"][conc_key])
            faixa_tm = cfg_primers["tm_c"][tm_key]
            a, b = float(faixa_tm[0]), float(faixa_tm[1])

            for p_nome in prims_nomes:
                if p_nome in linhas_por_nome:
                    l_dados = linhas_por_nome[p_nome]
                    t_ref = float(l_dados["tm_c"])
                    t_reac = float(l_dados["tm_reacao_c"])
                    todos_valores.extend([t_ref, t_reac])
                    info_p = {
                        "nome": p_nome,
                        "lado": lado_de(p_nome),
                        "tm_ref": t_ref,
                        "tm_reacao": t_reac,
                        "borda_inf": a,
                        "borda_sup": b,
                    }
                    prims_cat.append(info_p)
                    primers_medidos.append(info_p)

            if prims_cat:
                todos_valores.extend([a, b])
                cats_ativas.append({
                    "nome": cat_nome,
                    "faixa_a": a,
                    "faixa_b": b,
                    "conc_nM": conc_nM,
                    "primers": prims_cat,
                })

        if cats_ativas and todos_valores:
            val_min = min(todos_valores)
            val_max = max(todos_valores)
            lo = int(2 * math.floor((val_min - 2) / 2))
            hi = int(2 * math.ceil((val_max + 2) / 2))
            span = hi - lo

            def pos_val(v: float) -> str:
                return pct(v - lo, span)

            linhas_cal: list[dict] = []
            for c_info in cats_ativas:
                pontos: list[dict] = []
                for p_info in c_info["primers"]:
                    pontos.append({
                        "nome": p_info["nome"],
                        "lado": p_info["lado"],
                        "cheio": True,
                        "esquerda": pos_val(p_info["tm_ref"]),
                        "rotulo": True,
                    })
                    pontos.append({
                        "nome": p_info["nome"],
                        "lado": p_info["lado"],
                        "cheio": False,
                        "esquerda": pos_val(p_info["tm_reacao"]),
                        "rotulo": False,
                    })

                linhas_cal.append({
                    "nome": c_info["nome"],
                    "faixa": faixa(c_info["faixa_a"], c_info["faixa_b"]),
                    "marca": marca_tm,
                    "conc": f"{num_curto(c_info['conc_nM'])} nM",
                    "banda_esquerda": pos_val(c_info["faixa_a"]),
                    "banda_largura": pct(c_info["faixa_b"] - c_info["faixa_a"], span),
                    "pontos": pontos,
                })

            valores_ticks = list(range(lo, hi + 1, 2))
            ticks_cal: list[dict] = []
            for i, val in enumerate(valores_ticks):
                if i == 0:
                    alinhamento = "inicio"
                elif i == len(valores_ticks) - 1:
                    alinhamento = "fim"
                else:
                    alinhamento = "meio"
                v_str = f"{val} °C" if i == len(valores_ticks) - 1 else str(val)
                ticks_cal.append({
                    "valor": v_str,
                    "esquerda": pos_val(val),
                    "alinhamento": alinhamento,
                })

            total_med = len(primers_medidos)
            dentro_count = sum(
                1 for p in primers_medidos if p["borda_inf"] <= p["tm_ref"] <= p["borda_sup"]
            )
            if dentro_count == total_med:
                inicio_txt = f"Todos os {total_med} primers medidos ficam dentro da faixa."
            elif dentro_count == 0:
                inicio_txt = f"Nenhum dos {total_med} primers medidos fica dentro da faixa."
            else:
                inicio_txt = f"{dentro_count} dos {total_med} primers medidos ficam dentro da faixa."

            duplex = [p for p in primers_medidos if p["nome"] not in ("F1c", "B1c")]
            texto_cal = inicio_txt
            if duplex:
                dif_media = sum(p["tm_reacao"] - p["tm_ref"] for p in duplex) / len(duplex)
                sentido = "acima" if dif_media >= 0 else "abaixo"
                texto_cal += (
                    f" Na reação, a Tm dos primers em duplex fica em média {num(abs(dif_media), 1)} °C "
                    f"{sentido} da Tm de referência."
                )
            if len(duplex) < total_med:
                texto_cal += " F1c e B1c entram na reação como dobramento, fora da média."

            partes_desc_cal: list[str] = []
            for c_info in cats_ativas:
                p_strs = [
                    f"{p['nome']} {num(p['tm_ref'])} (reação: {num(p['tm_reacao'])})"
                    for p in c_info["primers"]
                ]
                partes_desc_cal.append(
                    f"{c_info['nome']}, faixa {faixa(c_info['faixa_a'], c_info['faixa_b'])}: "
                    f"{', '.join(p_strs)}"
                )
            desc_cal = f"Gráfico de Tm por tipo de primer: {'; '.join(partes_desc_cal)}."

            calibracao = {
                "linhas": linhas_cal,
                "ticks": ticks_cal,
                "texto": texto_cal,
                "descricao": desc_cal,
            }

    reacao = [
        {
            "nome": "Monovalentes",
            "valor": f"{num_curto(cond.monovalentes_mM)} mM",
            "nota": (
                f"Na⁺ {num_curto(cond.na_mM)} + K⁺ {num_curto(cond.k_mM)} + "
                f"NH₄⁺ {num_curto(cond.nh4_mM)} mM. Conferir: a mistura WarmStart é proprietária; "
                "o tampão isotérmico 1X da NEB tem 50 mM de KCl e 10 mM de (NH₄)₂SO₄."
            ),
            "destaque": False,
        },
        {
            "nome": "Mg²⁺",
            "valor": f"{num_curto(cond.mg_mM)} mM",
            "nota": (
                f"Livre: {num_curto(cond.mg_livre_mM)} mM, o total menos o dNTP total."
                if cond.mg_livre_mM > 0
                else (
                    f"Livre: 0 mM: o dNTP total ({num_curto(cond.dntp_total_mM)} mM) não é menor que o "
                    f"Mg²⁺ total ({num_curto(cond.mg_mM)} mM). Confira se reacao.dntp_cada_mM está em mM "
                    "de cada dNTP, e não do total."
                )
            ),
            "destaque": cond.mg_livre_mM <= 0,
        },
        {
            "nome": "dNTP",
            "valor": f"{num_curto(cond.dntp_cada_mM)} mM cada",
            "nota": f"Total de {num_curto(cond.dntp_total_mM)} mM, o valor enviado ao primer3.",
            "destaque": False,
        },
        {
            "nome": "Temperatura da reação",
            "valor": f"{num_curto(cond.temp_C)} °C",
            "nota": "Grampos e dímeros são avaliados nesta temperatura.",
            "destaque": False,
        },
        {
            "nome": "Tm de referência",
            "valor": (
                f"{num_curto(conv.oligo_nM)} nM · Na⁺ {num_curto(conv.na_mM)} mM · "
                f"Mg²⁺ {num_curto(conv.mg_mM)} mM"
            ),
            "nota": "Convenção do PrimerExplorer V5, a das faixas da Eiken.",
            "destaque": True,
        },
        {
            "nome": "ΔG de ponta",
            "valor": "37 °C · NaCl 1 M",
            "nota": (
                "Bases da ponta com iniciação, tabela unificada de SantaLucia (1998), "
                "como no PrimerExplorer V5."
            ),
            "destaque": False,
        },
    ]

    distancias_param: list[dict] = []
    cfg_dist = cfg_primers.get("distancias", {})
    dist_calc = (
        calcular_distancias(ref_exec.resultado.regioes)
        if (ref_exec and ref_exec.resultado.regioes)
        else {}
    )
    for k_dist, faixa_val in cfg_dist.items():
        nome_dist = NOMES_DISTANCIA.get(k_dist, k_dist)
        faixa_str = faixa(faixa_val[0], faixa_val[1])
        if k_dist in dist_calc:
            val_d = dist_calc[k_dist]
            dentro = faixa_val[0] <= val_d <= faixa_val[1]
            valor_str = milhar(val_d)
        else:
            valor_str = "—"
            dentro = None
        distancias_param.append({
            "nome": nome_dist,
            "faixa": faixa_str,
            "valor": valor_str,
            "dentro": dentro,
        })

    c_min = cfg_primers["comprimento"]["min"]
    c_max = cfg_primers["comprimento"]["max"]
    l_min = cfg_primers["comprimento"]["loop_min"]
    l_max = cfg_primers["comprimento"]["loop_max"]
    gc_min = cfg_primers["gc_pct"]["min"]
    gc_max = cfg_primers["gc_pct"]["max"]
    bases = cfg_primers["extremidade"]["bases"]
    dg_max = cfg_primers["extremidade"]["dg_max"]
    ext_c = num_curto(cfg_primers["concentracao_nM"]["externos"])
    int_c = num_curto(cfg_primers["concentracao_nM"]["internos"])
    loop_c = num_curto(cfg_primers["concentracao_nM"]["loops"])

    regras = [
        {
            "nome": "Comprimento",
            "valor": f"{faixa(c_min, c_max)} nt · loops {faixa(l_min, l_max)} nt",
            "marca": "",
        },
        {
            "nome": "GC",
            "valor": f"{faixa(gc_min, gc_max)} %",
            "marca": "conferir",
        },
        {
            "nome": "Ponta estável",
            "valor": (
                f"ΔG°37 das {bases} bases ≤ {num_curto(dg_max)} kcal/mol · "
                f"3′ de F3, F2, B2, B3, LF, LB; 5′ de F1c e B1c"
            ),
            "marca": "",
        },
        {
            "nome": "Concentração",
            "valor": f"externos {ext_c} nM · internos {int_c} nM · loops {loop_c} nM",
            "marca": "pendente",
        },
        {
            "nome": "Dímeros",
            "valor": "Só reportados, com a ponta 3′ pareada — limiar a definir com dado",
            "marca": "",
        },
        {
            "nome": "Tm de referência",
            "valor": texto_faixas_tm(cfg_primers),
            "marca": marca_tm or "conferir",
        },
    ]

    d_cfg = cfg_primers.get("desenho", {})
    janela = d_cfg.get("janela_nt", 10)
    por_janela = d_cfg.get("candidatos_por_janela", 2)
    guardados = d_cfg.get("nucleos_guardados", 5000)
    sobre = d_cfg.get("sobreposicao_max_nt", 10)
    pre_sel = d_cfg.get("pre_selecao", 50)
    peso = d_cfg.get("peso_dimero", 0.5)
    n_res = d_cfg.get("n_resultados", 5)
    ligante = d_cfg.get("ligante", "")
    lig_val = ligante if ligante else "nenhum"

    desenho_param = [
        {
            "nome": "Teto por janela",
            "valor": f"{por_janela} por tipo a cada {janela} nt",
            "marca": "",
        },
        {
            "nome": "Núcleos guardados",
            "valor": f"até {milhar(guardados)}, do alvo inteiro",
            "marca": "",
        },
        {
            "nome": "Pré-seleção",
            "valor": f"{pre_sel} conjuntos para a etapa de dímeros",
            "marca": "",
        },
        {
            "nome": "Peso do dímero",
            "valor": f"{num_curto(peso)} por kcal/mol",
            "marca": "provisório",
        },
        {
            "nome": "Resultados",
            "valor": str(n_res),
            "marca": "",
        },
        {
            "nome": "Ligante FIP/BIP",
            "valor": lig_val,
            "marca": "",
        },
        {
            "nome": "Propostas distintas",
            "valor": f"F2 e B2c não sobrepõem mais de {sobre} nt ao mesmo tempo",
            "marca": "",
        },
    ]

    return {
        "config": "config/primers.toml",
        "tm_a_calibrar": tm_a_calibrar,
        "referencia": referencia,
        "calibracao": calibracao,
        "reacao": reacao,
        "distancias": distancias_param,
        "regras": regras,
        "desenho": desenho_param,
    }


def _numero_cfg(cfg: object, chave: str, padrao: float) -> float:
    """Número positivo e finito do config, ou o padrão quando falta ou é inválido."""
    if not isinstance(cfg, dict):
        return float(padrao)
    try:
        valor = float(cfg.get(chave, padrao))
    except (TypeError, ValueError):
        return float(padrao)
    if not math.isfinite(valor) or valor <= 0.0:
        return float(padrao)
    return valor


def _exemplares_manuais(cfg_tubos: dict, duracao_padrao: float) -> list[dict]:
    """Curvas manuais de [[tubos.exemplares]], com tempos em minutos e valores saneados."""
    exemplares_cfg = cfg_tubos.get("exemplares")
    if not exemplares_cfg or not isinstance(exemplares_cfg, list):
        return [dict(t) for t in PADRAO_TUBOS]
    tubos = []
    for ex in exemplares_cfg:
        if not isinstance(ex, dict):
            continue
        curva = ex.get("curva")
        if curva not in ("linear", "sigmoide"):
            curva = "linear"
        try:
            inicio = max(0.0, float(ex.get("inicio", 0.0)))
        except (ValueError, TypeError):
            inicio = 0.0
        try:
            duracao = float(ex.get("duracao", duracao_padrao))
            if duracao <= 0.0:
                duracao = 1.0
        except (ValueError, TypeError):
            duracao = 1.0
        try:
            alcance = min(100.0, max(0.0, float(ex.get("alcance", 100.0))))
        except (ValueError, TypeError):
            alcance = 100.0
        tubos.append({
            "nome": str(ex.get("nome", "Amostra")),
            "curva": curva,
            "inicio": inicio,
            "duracao": duracao,
            "alcance": alcance,
        })
    return tubos or [dict(t) for t in PADRAO_TUBOS]


def _serie_simulada(cfg_tubos: dict) -> list[dict]:
    """Tubos de [[tubos.serie]]: nome, média de cópias, papel e eficiência própria opcional."""
    serie_cfg = cfg_tubos.get("serie")
    if not serie_cfg or not isinstance(serie_cfg, list):
        return [dict(t) for t in PADRAO_SERIE]
    serie = []
    for item in serie_cfg:
        if not isinstance(item, dict):
            continue
        try:
            copias = float(item.get("copias", -1.0))
        except (TypeError, ValueError):
            continue
        if not math.isfinite(copias) or copias < 0.0:
            continue
        papel = item.get("papel", "amostra")
        tubo = {
            "nome": str(item.get("nome", "Tubo")),
            "copias": copias,
            "papel": papel if papel in PAPEIS else "amostra",
        }
        if "eficiencia" in item:
            try:
                ef = float(item["eficiencia"])
            except (TypeError, ValueError):
                ef = -1.0
            if 0.0 <= ef <= 1.0:
                tubo["eficiencia"] = ef
        serie.append(tubo)
    return serie or [dict(t) for t in PADRAO_SERIE]


def texto_eficiencia(ef: float | None) -> str:
    """Eficiência como porcentagem: '100 %', '61 %', '1,5 %', '0,073 %', '0 %' ou '—'."""
    if ef is None:
        return "—"
    v = ef * 100.0
    if v >= 10.0:
        return f"{num(v, 0)} %"
    if v >= 1.0:
        return f"{num(v, 1)} %"
    if v > 0.0:
        return f"{num(v, 3)} %"
    return "0 %"


def _primeiro_positivo(t_min: list[float], ph: list[float]) -> float | None:
    """Primeiro minuto em que a fração ácida chega ao limiar de positivo do lampsim."""
    pka = float(QUIMICA["pka"])
    limiar = float(LEITURA["limiar_positivo"])
    for t, v in zip(t_min, ph, strict=True):
        if 1.0 / (1.0 + 10.0 ** (v - pka)) >= limiar:
            return float(t)
    return None


def tela_tubos(
    cfg_web: dict,
    cfg_esteira: dict,
    execucao: Execucao | None = None,
    modo: str = "simulacao",
) -> dict:
    """Contexto da página /tubos: tubos simulados com o conjunto analisado, ou curvas manuais.

    modo "simulacao" usa a eficiência prevista pela análise (esteira.ponte) e o lampsim; cai
    para as curvas manuais, com aviso, quando a execução não tem incompatibilidades medidas.
    modo "manual" mostra as curvas de [[tubos.exemplares]], sempre com aviso.
    """
    cfg_tubos = cfg_web.get("tubos", {}) if isinstance(cfg_web, dict) else {}
    if not isinstance(cfg_tubos, dict):
        cfg_tubos = {}

    duracao_max = _numero_cfg(cfg_tubos, "duracao_max_min", TUBOS_PADRAO_TEMPO["duracao_max_min"])
    duracao = min(_numero_cfg(cfg_tubos, "duracao_min", TUBOS_PADRAO_TEMPO["duracao_min"]), duracao_max)
    passo = _numero_cfg(cfg_tubos, "passo_min", TUBOS_PADRAO_TEMPO["passo_min"])
    velocidade = _numero_cfg(cfg_tubos, "velocidade", TUBOS_PADRAO_TEMPO["velocidade"])
    velocidades_cfg = cfg_tubos.get("velocidades", TUBOS_PADRAO_VELOCIDADES)
    if not isinstance(velocidades_cfg, list):
        velocidades_cfg = TUBOS_PADRAO_VELOCIDADES
    velocidades = {velocidade}
    for v in velocidades_cfg:
        if isinstance(v, int | float) and not isinstance(v, bool) and math.isfinite(v) and v > 0:
            velocidades.add(float(v))
    try:
        semente = int(cfg_tubos.get("semente", TUBOS_PADRAO_SEMENTE))
    except (TypeError, ValueError):
        semente = TUBOS_PADRAO_SEMENTE
    cfg_manual = cfg_tubos.get("manual", {})
    ph_inicial = _numero_cfg(cfg_manual, "ph_inicial", TUBOS_PADRAO_MANUAL["ph_inicial"])
    queda_ph = _numero_cfg(cfg_manual, "queda_ph", TUBOS_PADRAO_MANUAL["queda_ph"])
    paleta = paleta_progresso(ph_inicial, queda_ph)

    origem = None
    eficiencia = None
    if execucao is not None and execucao.tipo == "analise":
        incompat = getattr(execucao.resultado, "incompat", None)
        eficiencia = eficiencia_prevista(incompat, cfg_esteira)
        if incompat is None:
            incompat_texto = "não medidas"
        elif incompat.total is None:
            incompat_texto = f"região sem posição no alvo: {juntar(list(incompat.nao_encontradas))}"
        else:
            incompat_texto = f"{incompat.total} no total, {incompat.ponta3} na ponta que estende"
        origem = {
            "id": execucao.id,
            "conjunto": execucao.rotulo_conjunto,
            "alvo": execucao.rotulo_alvo,
            "url_analise": f"/analise/{execucao.id}",
            "eficiencia": eficiencia,
            "eficiencia_texto": texto_eficiencia(eficiencia),
            "incompat_total": None if incompat is None else incompat.total,
            "incompat_3p": None if incompat is None else incompat.ponta3,
            "incompat_texto": incompat_texto,
            "nao_encontradas": [] if incompat is None else list(incompat.nao_encontradas),
        }

    pode_simular = eficiencia is not None
    url_simulacao = f"/tubos?execucao={origem['id']}" if origem else "/tubos"
    url_manual = f"{url_simulacao}{'&' if origem else '?'}modo=manual"
    leitura_min = float(TEMPO["leitura_nominal_min"])

    base = {
        "duracao_min": duracao,
        "duracao_max_min": duracao_max,
        "passo_min": passo,
        "leitura_min": leitura_min,
        "leitura_texto": f"{num_curto(leitura_min)} min",
        "pka": float(QUIMICA["pka"]),
        "velocidade": velocidade,
        "velocidades": sorted(velocidades),
        "paleta": paleta,
        "cor_inicial": paleta[0],
        "cor_final": paleta[-1],
        "origem": origem,
        "pode_simular": pode_simular,
        "url_simulacao": url_simulacao,
        "url_manual": url_manual,
    }

    if modo == "manual" or not pode_simular:
        if modo == "manual":
            aviso = AVISO_TUBOS_MANUAL
        elif origem is None:
            aviso = AVISO_TUBOS_SEM_ANALISE
        else:
            aviso = AVISO_TUBOS_SEM_MEDIDA.format(id=origem["id"])
        return base | {
            "modo": "manual",
            "aviso": aviso,
            "t_min": [],
            "tubos": _exemplares_manuais(cfg_tubos, duracao),
        }

    t = eixo_minutos(duracao_max, passo)
    t_min = [round(float(x), 3) for x in t]
    tubos = simular_tubos(eficiencia, _serie_simulada(cfg_tubos), t, semente)
    for tubo in tubos:
        t_pos = _primeiro_positivo(t_min, tubo["ph"])
        tubo["t_positivo_min"] = t_pos
        tubo["t_positivo_texto"] = "não vira" if t_pos is None else f"{num(t_pos, 1)} min"
        tubo["diagnostico_texto"] = DIAGNOSTICO_TEXTO.get(tubo["diagnostico"], tubo["diagnostico"])
        tubo["papel_texto"] = PAPEL_TEXTO.get(tubo["papel"], tubo["papel"])
        tubo["copias_media_texto"] = milhar(round(tubo["copias_media"]))
        tubo["copias_texto"] = milhar(tubo["copias"])
        tubo["eficiencia_texto"] = texto_eficiencia(tubo["eficiencia"])

    faltam = origem["nao_encontradas"] if origem else []
    if eficiencia == 0.0:
        aviso = (
            f"Região essencial sem posição no alvo ({juntar(faltam)}): eficiência prevista 0, "
            "nenhum tubo de amostra amplifica."
            if faltam
            else "Eficiência prevista 0: nenhum tubo de amostra amplifica."
        )
    elif faltam:
        aviso = f"Laço sem posição no alvo ({juntar(faltam)}): conta como metade da eficiência."
    else:
        aviso = None
    return base | {"modo": "simulacao", "aviso": aviso, "t_min": t_min, "tubos": tubos}
