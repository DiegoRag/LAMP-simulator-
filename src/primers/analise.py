"""Módulo de análise termodinâmica e estrutural de primers LAMP."""

from dataclasses import dataclass
from primers.entrada import Leitura, ORDEM_CANONICA
from primers.modelo import Regiao
from primers.restricoes import Violacao, localizar, verificar, verificar_primers
from primers.termo import Condicoes, dg_extremidade, gc_pct, grampo, heterodimero, homodimero, tm

PRIMERS_TUBO = ["F3", "B3", "FIP", "BIP", "LF", "LB"]


@dataclass
class Resultado:
    """Resultado da análise termodinâmica e de restrições de um conjunto de primers."""

    linhas: list[dict]
    heterodimeros: dict
    violacoes: list[Violacao]
    avisos: list[str]
    regioes: dict[str, Regiao] | None


def analisar(leitura: Leitura, alvo: str | None, cfg: dict) -> Resultado:
    """Analisa propriedades físico-químicas, termodinâmicas e restrições dos primers.

    Parâmetros:
        leitura: Leitura contendo primers e conjunto LAMP.
        alvo: sequência-alvo normalizada (opcional).
        cfg: dicionário com a configuração do primers.toml.

    Retorna:
        Resultado consolidado com métricas por primer, heterodímeros, violações e regiões.
    """
    cond = Condicoes.de_config(cfg)
    conc_cfg = cfg["concentracao_nM"]
    bases_ext = int(cfg["extremidade"]["bases"])

    linhas: list[dict] = []
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

        g = grampo(seq, cond, conc_nM=conc)
        h = homodimero(seq, cond, conc_nM=conc)

        if nome in ("FIP", "BIP"):
            tm_val = None
            dg_ext = None
            extremidade = None
        else:
            tm_val = tm(seq, cond, conc_nM=conc)
            if nome in ("F1c", "B1c"):
                extremidade = "5'"
                dg_ext = dg_extremidade(seq, "5", bases_ext, cond)
            else:
                extremidade = "3'"
                dg_ext = dg_extremidade(seq, "3", bases_ext, cond)

        aviso_linha = ""
        if nome == "FIP" and "F1c" not in leitura.primers:
            aviso_linha = "FIP sem F1c/F2"
        elif nome == "BIP" and "B1c" not in leitura.primers:
            aviso_linha = "BIP sem B1c/B2"

        linhas.append(
            {
                "nome": nome,
                "sequencia": seq,
                "comprimento": comprimento,
                "gc_pct": gc,
                "tm_c": tm_val,
                "dg_extremidade_kcal": dg_ext,
                "extremidade": extremidade,
                # Sem estrutura, o primer3 devolve Tm = 0 e dG = 0; exportar 0.00 sugeriria um
                # grampo que derrete a 0 grau. Celula vazia = nenhuma estrutura encontrada.
                "grampo_tm_c": g.tm_c if g.encontrada else None,
                "grampo_dg_kcal": g.dg_kcal if g.encontrada else None,
                "homodimero_dg_kcal": h.dg_kcal if h.encontrada else None,
                "avisos": aviso_linha,
            }
        )

    # Heterodímeros entre primers de tubo presentes
    nomes_tubo = [n for n in PRIMERS_TUBO if n in leitura.primers]
    matriz_het: list[list[float]] = []
    for a in nomes_tubo:
        linha_het: list[float] = []
        seq_a = leitura.primers[a]
        for b in nomes_tubo:
            seq_b = leitura.primers[b]
            het = heterodimero(seq_a, seq_b, cond)
            linha_het.append(het.dg_kcal)
        matriz_het.append(linha_het)

    heterodimeros = {"nomes": nomes_tubo, "dg_kcal": matriz_het}

    # Violações
    if leitura.conjunto is not None:
        violacoes = verificar(leitura.conjunto, alvo, cfg)
    else:
        violacoes = verificar_primers(leitura.primers, cfg)

    # Regiões
    if leitura.conjunto is not None and alvo is not None:
        regioes, _ = localizar(leitura.conjunto, alvo)
    else:
        regioes = None

    # Avisos
    avisos = list(leitura.avisos)
    if alvo is not None and leitura.conjunto is None:
        avisos.append("alvo informado, mas o conjunto está incompleto")

    return Resultado(
        linhas=linhas,
        heterodimeros=heterodimeros,
        violacoes=violacoes,
        avisos=avisos,
        regioes=regioes,
    )
