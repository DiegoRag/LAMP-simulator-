"""Módulo de análise termodinâmica e estrutural de primers LAMP."""

import re
from dataclasses import dataclass
from pathlib import Path

from primers.casamento import Posicao
from primers.entrada import ORDEM_CANONICA, Leitura, ler_conjunto
from primers.incompat import Incompat, medir_incompat
from primers.modelo import Regiao, normalizar, rc
from primers.restricoes import (
    Violacao,
    calcular_distancias,
    localizar_detalhado,
    verificar,
    verificar_primers,
)
from primers.termo import (
    Condicoes,
    CondicoesReferencia,
    dg_extremidade,
    dimero_ponta3,
    gc_pct,
    grampo,
    heterodimero,
    homodimero,
    tm,
    tm_dobramento,
    tm_referencia,
)

PRIMERS_TUBO = ["F3", "B3", "FIP", "BIP", "LF", "LB"]


@dataclass
class Resultado:
    """Resultado da análise termodinâmica e de restrições de um conjunto de primers."""

    linhas: list[dict]
    heterodimeros: dict
    violacoes: list[Violacao]
    avisos: list[str]
    regioes: dict[str, Regiao] | None
    fita: str | None = None
    regioes_aproximadas: dict[str, Posicao] | None = None
    incompat: Incompat | None = None  # só com conjunto completo e alvo (fase 5)


def analisar(leitura: Leitura, alvo: str | None, cfg: dict) -> Resultado:
    """Analisa propriedades físico-químicas, termodinâmicas e restrições dos primers.

    Parâmetros:
        leitura: Leitura contendo primers e conjunto LAMP.
        alvo: sequência-alvo normalizada (opcional).
        cfg: dicionário com a configuração do primers.toml.

    Retorna:
        Resultado consolidado com métricas por primer, heterodímeros, violações e regiões.
    """
    if not leitura.primers:
        raise ValueError("Nenhum primer informado para análise.")
    for nome, seq_bruta in leitura.primers.items():
        seq_norm = normalizar(seq_bruta)
        if not seq_norm:
            raise ValueError(f"Sequência vazia para o primer '{nome}'.")
    if alvo is not None:
        alvo = normalizar(alvo)
        if not alvo:
            raise ValueError("Sequência-alvo vazia.")

    cond = Condicoes.de_config(cfg)
    ref = CondicoesReferencia.de_config(cfg)
    conc_cfg = cfg["concentracao_nM"]
    bases_ext = int(cfg["extremidade"]["bases"])
    max_nt_estrutura = int(cfg.get("dimeros", {}).get("max_nt_estrutura", 60))

    # Regiões e localização (antes das linhas: o laço de F1c/B1c depende delas)
    if leitura.conjunto is not None and alvo is not None:
        loc = localizar_detalhado(leitura.conjunto, alvo, cfg)
        regioes = loc.regioes
        fita = loc.fita
        regioes_aproximadas = loc.aproximadas
        dists = calcular_distancias(loc.regioes_completas)
        bases_ponta = int(cfg.get("casamento", {}).get("bases_ponta_3p", 5))
        incompat = medir_incompat(leitura.conjunto, loc, bases_ponta)
    else:
        loc = None
        regioes = None
        fita = None
        regioes_aproximadas = None
        dists = {}
        incompat = None

    if leitura.conjunto is not None:
        ligante = leitura.conjunto.ligante
    else:
        ligante = str(cfg.get("desenho", {}).get("ligante", ""))
    dist_padrao = int(cfg.get("dobramento", {}).get("distancia_nt_padrao", 50))
    d_f2_f1 = dists.get("F2_F1", dist_padrao)
    d_b1c_b2c = dists.get("B1c_B2c", dist_padrao)
    # geometria fora de ordem (laço < 3 nt) não serve para o laço: vale a distância presumida
    laco_f1c = len(ligante) + (d_f2_f1 if len(ligante) + d_f2_f1 >= 3 else dist_padrao)
    laco_b1c = len(ligante) + (d_b1c_b2c if len(ligante) + d_b1c_b2c >= 3 else dist_padrao)

    linhas: list[dict] = []
    avisos_estruturas: list[str] = []

    for nome in ORDEM_CANONICA:
        if nome not in leitura.primers:
            continue

        seq = leitura.primers[nome]
        comprimento = len(seq)
        gc = gc_pct(seq)

        if nome in ("F3", "B3"):
            conc = float(conc_cfg["externos"])
        elif nome in ("LF", "LB"):
            conc = float(conc_cfg["loops"])
        else:  # F2, F1c, B2, B1c, FIP, BIP
            conc = float(conc_cfg["internos"])

        aviso_linha = ""
        if nome == "FIP" and "F1c" not in leitura.primers:
            aviso_linha = "FIP sem F1c/F2"
        elif nome == "BIP" and "B1c" not in leitura.primers:
            aviso_linha = "BIP sem B1c/B2"

        if comprimento > max_nt_estrutura:
            g_tm = None
            g_dg = None
            h_dg = None
            aviso_nt = f"acima de {max_nt_estrutura} nt (estruturas secundárias omitidas)"
            aviso_linha = f"{aviso_linha}; {aviso_nt}" if aviso_linha else aviso_nt
            avisos_estruturas.append(
                f"Primer {nome} tem {comprimento} nt (acima do limite de "
                f"{max_nt_estrutura} nt do Primer3); grampo e dímeros não calculados."
            )
        else:
            g = grampo(seq, cond, conc_nM=conc)
            h = homodimero(seq, cond, conc_nM=conc)
            g_tm = g.tm_c if g.encontrada else None
            g_dg = g.dg_kcal if g.encontrada else None
            h_dg = h.dg_kcal if h.encontrada else None

        if nome in ("FIP", "BIP"):
            tm_val = None
            tm_reacao_val = None
            dg_ext = None
            extremidade = None
        else:
            tm_val = tm_referencia(seq, ref)
            if nome == "F1c":
                tm_reacao_val = tm_dobramento(seq, laco_f1c, cond)
            elif nome == "B1c":
                tm_reacao_val = tm_dobramento(seq, laco_b1c, cond)
            else:
                tm_reacao_val = tm(seq, cond, conc_nM=conc)
            if nome in ("F1c", "B1c"):
                extremidade = "5'"
                dg_ext = dg_extremidade(seq, "5", bases_ext, ref)
            else:
                extremidade = "3'"
                dg_ext = dg_extremidade(seq, "3", bases_ext, ref)

        linhas.append(
            {
                "nome": nome,
                "sequencia": seq,
                "comprimento": comprimento,
                "gc_pct": gc,
                "tm_c": tm_val,
                "tm_reacao_c": tm_reacao_val,
                "dg_extremidade_kcal": dg_ext,
                "extremidade": extremidade,
                "grampo_tm_c": g_tm,
                "grampo_dg_kcal": g_dg,
                "homodimero_dg_kcal": h_dg,
                "avisos": aviso_linha,
            }
        )

    # Heterodímeros entre primers de tubo presentes
    nomes_tubo = [n for n in PRIMERS_TUBO if n in leitura.primers]
    matriz_ponta3: list[list[float]] = []
    matriz_qualquer: list[list[float]] = []
    for a in nomes_tubo:
        linha_p3: list[float] = []
        linha_q: list[float] = []
        seq_a = leitura.primers[a]
        for b in nomes_tubo:
            seq_b = leitura.primers[b]
            if len(seq_a) > max_nt_estrutura or len(seq_b) > max_nt_estrutura:
                linha_p3.append(0.0)
                linha_q.append(0.0)
            else:
                linha_p3.append(dimero_ponta3(seq_a, seq_b, cond).dg_kcal)
                linha_q.append(heterodimero(seq_a, seq_b, cond).dg_kcal)
        matriz_ponta3.append(linha_p3)
        matriz_qualquer.append(linha_q)

    heterodimeros = {"nomes": nomes_tubo, "dg_kcal": matriz_ponta3, "dg_qualquer_kcal": matriz_qualquer}

    # Avisos
    avisos = list(leitura.avisos)
    for a_est in avisos_estruturas:
        if a_est not in avisos:
            avisos.append(a_est)
    if alvo is not None and leitura.conjunto is None:
        avisos.append("alvo informado, mas o conjunto está incompleto")

    # Regiões e localização
    if loc is not None:
        for av in loc.avisos:
            if av not in avisos:
                avisos.append(av)

    if alvo is not None and "N" in alvo:
        n_n = alvo.count("N")
        trechos_n = len(re.findall(r"N+", alvo))
        txt_milhar = f"{n_n:,}".replace(",", ".")
        aviso_n = (
            f"O alvo tem {txt_milhar} bases N (desconhecidas) em {trechos_n} trecho(s); "
            "nenhum primer é localizado sobre elas."
        )
        if aviso_n not in avisos:
            avisos.append(aviso_n)

    # Violações
    if leitura.conjunto is not None:
        if alvo is not None and fita == "-":
            v_test = verificar(leitura.conjunto, alvo, cfg, loc=loc)
            if any(v.regra == "nao_encontrado" for v in v_test):
                violacoes = verificar(leitura.conjunto, rc(alvo), cfg)
            else:
                violacoes = v_test
        else:
            violacoes = verificar(leitura.conjunto, alvo, cfg, loc=loc)
    else:
        violacoes = verificar_primers(leitura.primers, cfg)

    return Resultado(
        linhas=linhas,
        heterodimeros=heterodimeros,
        violacoes=violacoes,
        avisos=avisos,
        regioes=regioes,
        fita=fita,
        regioes_aproximadas=regioes_aproximadas,
        incompat=incompat,
    )


def analisar_registros(
    leitura_ou_caminho: Leitura | str | Path,
    registros: list,
    cfg: dict,
) -> dict[str, Resultado]:
    """Analisa um arquivo de conjunto de primers contra cada registro de um FASTA multi-registro."""
    ligante_cfg = str(cfg.get("desenho", {}).get("ligante", ""))
    max_erros_cfg = int(cfg.get("casamento", {}).get("max_erros", 3))
    circular_cfg = bool(cfg.get("alvo", {}).get("circular", False))
    resultados: dict[str, Resultado] = {}
    for reg in registros:
        if isinstance(leitura_ou_caminho, Leitura):
            leitura_reg = leitura_ou_caminho
        else:
            leitura_reg = ler_conjunto(
                leitura_ou_caminho,
                alvo=reg.sequencia,
                ligante=ligante_cfg,
                max_erros=max_erros_cfg,
                circular=circular_cfg,
            )
        resultados[reg.id] = analisar(leitura_reg, reg.sequencia, cfg)
    return resultados
