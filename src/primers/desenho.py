"""Módulo de desenho de conjuntos de primers LAMP."""

from bisect import bisect_left, bisect_right
from dataclasses import dataclass
from primers.modelo import ConjuntoLAMP, Regiao, normalizar, primers_de_regioes, rc
from primers.restricoes import calcular_distancias, verificar, verificar_primers
from primers.termo import Condicoes, heterodimero, tm

AVISO = (
    "Candidatos computacionais, não validados experimentalmente. A definição do conjunto "
    "que vai ao tubo é da equipe biomédica (ADR-0013)."
)

AVISO_ASCII = (
    "Candidatos computacionais, nao validados experimentalmente. A definicao do conjunto "
    "que vai ao tubo e da equipe biomedica (ADR-0013)."
)

INFO_TIPO = {
    # tipo: (nome_primer, eh_rc, chave_tm, chave_conc, eh_loop)
    "F3": ("F3", False, "F3_B3", "externos", False),
    "F2": ("F2", False, "F2_B2", "internos", False),
    "F1": ("F1c", True, "F1c_B1c", "internos", False),
    "B1c": ("B1c", False, "F1c_B1c", "internos", False),
    "B2c": ("B2", True, "F2_B2", "internos", False),
    "B3c": ("B3", True, "F3_B3", "externos", False),
    "LF": ("LF", True, "loops", "loops", True),
    "LB": ("LB", False, "loops", "loops", True),
}


@dataclass(frozen=True)
class Candidato:
    """Candidato individual a uma região de primer LAMP."""

    tipo: str
    regiao: Regiao
    primer: str
    tm_c: float
    desvio_tm: float


@dataclass
class ConjuntoDesenhado:
    """Conjunto LAMP desenhado e ranqueado."""

    posicao: int
    pontuacao: float
    desvio_tm: float
    pior_heterodimero_kcal: float
    conjunto: ConjuntoLAMP
    regioes: dict[str, Regiao]
    avisos: list[str]


@dataclass
class ResultadoDesenho:
    """Resultado completo da etapa de desenho de primers LAMP."""

    conjuntos: list[ConjuntoDesenhado]
    estatisticas: dict


def gerar_candidatos(
    alvo: str,
    cfg: dict,
    tipos: tuple[str, ...] = ("F3", "F2", "F1", "B1c", "B2c", "B3c"),
    dentro: Regiao | None = None,
) -> dict[str, list[Candidato]]:
    """Gera candidatos de primers para os tipos de região especificados.

    1. Para cada tipo, cada comprimento da faixa e cada ini: janela [ini, ini+L).
       Com dentro, só janelas inteiramente contidas em dentro.
    2. Descarte a janela se a sequência da região ocorrer mais de uma vez no alvo.
    3. Mantenha só se verificar_primers({nome_do_primer: primer}, cfg) == [].
    4. Teto: se sobrarem mais que cfg["desenho"]["candidatos_por_regiao"],
       fique com os de menor desvio_tm, desempate por (regiao.ini, regiao.fim).
    5. Devolva cada lista ordenada por (regiao.ini, regiao.fim).
    """
    alvo = normalizar(alvo)
    cond = Condicoes.de_config(cfg)
    teto = int(cfg["desenho"]["candidatos_por_regiao"])
    len_alvo = len(alvo)

    resultado: dict[str, list[Candidato]] = {}

    for tipo in tipos:
        nome_primer, eh_rc, chave_tm, chave_conc, eh_loop = INFO_TIPO[tipo]

        if eh_loop:
            c_min = int(cfg["comprimento"]["loop_min"])
            c_max = int(cfg["comprimento"]["loop_max"])
        else:
            c_min = int(cfg["comprimento"]["min"])
            c_max = int(cfg["comprimento"]["max"])

        conc = float(cfg["concentracao_nM"][chave_conc])
        tm_cfg = cfg["tm_c"][chave_tm]
        centro_tm = (float(tm_cfg[0]) + float(tm_cfg[1])) / 2.0

        candidatos_tipo: list[Candidato] = []

        for L in range(c_min, c_max + 1):
            if dentro is not None:
                ini_min = dentro.ini
                ini_max = dentro.fim - L
            else:
                ini_min = 0
                ini_max = len_alvo - L

            if ini_max < ini_min:
                continue

            for ini in range(ini_min, ini_max + 1):
                fim = ini + L
                seq = alvo[ini:fim]

                # 2. Descarte a janela se a sequência ocorrer mais de uma vez no alvo
                if alvo.find(seq) != ini or alvo.find(seq, ini + 1) != -1:
                    continue

                primer_seq = rc(seq) if eh_rc else seq

                # 3. Mantenha só se verificar_primers == []
                if verificar_primers({nome_primer: primer_seq}, cfg):
                    continue

                tm_val = tm(primer_seq, cond, conc)
                desvio = abs(tm_val - centro_tm)

                candidatos_tipo.append(
                    Candidato(
                        tipo=tipo,
                        regiao=Regiao(ini, fim),
                        primer=primer_seq,
                        tm_c=tm_val,
                        desvio_tm=desvio,
                    )
                )

        # 4. Teto
        if len(candidatos_tipo) > teto:
            candidatos_tipo.sort(key=lambda c: (c.desvio_tm, c.regiao.ini, c.regiao.fim))
            candidatos_tipo = candidatos_tipo[:teto]

        # 5. Ordenação por (regiao.ini, regiao.fim)
        candidatos_tipo.sort(key=lambda c: (c.regiao.ini, c.regiao.fim))
        resultado[tipo] = candidatos_tipo

    return resultado


def desenhar(alvo: str, cfg: dict, n: int | None = None) -> ResultadoDesenho:
    """Desenha conjuntos de primers LAMP para o alvo conforme restrições e ranking."""
    alvo = normalizar(alvo)
    cands = gerar_candidatos(alvo, cfg)

    D = cfg["distancias"]
    L_min = int(cfg["comprimento"]["min"])
    L_max = int(cfg["comprimento"]["max"])
    max_examinados = int(cfg["desenho"]["conjuntos_examinados"])

    f2_list = cands.get("F2", [])
    f1_list = cands.get("F1", [])
    b1c_list = cands.get("B1c", [])
    b2c_list = cands.get("B2c", [])
    f3_list = cands.get("F3", [])
    b3c_list = cands.get("B3c", [])

    f1_inis = [c.regiao.ini for c in f1_list]
    b1c_inis = [c.regiao.ini for c in b1c_list]
    b2c_inis = [c.regiao.ini for c in b2c_list]
    f3_inis = [c.regiao.ini for c in f3_list]
    b3c_inis = [c.regiao.ini for c in b3c_list]

    conjuntos_examinados = 0
    busca_truncada = False
    combinacoes: list[tuple[Candidato, Candidato, Candidato, Candidato, Candidato, Candidato]] = []

    # Cache do melhor F3 por ini(F2)
    cache_f3: dict[int, Candidato | None] = {}
    for c_f2 in f2_list:
        ini_f2 = c_f2.regiao.ini
        if ini_f2 not in cache_f3:
            min_fim_f3 = ini_f2 - D["F3_F2_gap"][1]
            max_fim_f3 = ini_f2 - D["F3_F2_gap"][0]
            if min_fim_f3 > max_fim_f3:
                cache_f3[ini_f2] = None
            else:
                idx_start = bisect_left(f3_inis, min_fim_f3 - L_max)
                idx_end = bisect_right(f3_inis, max_fim_f3 - L_min)
                validos_f3 = [
                    c for c in f3_list[idx_start:idx_end]
                    if min_fim_f3 <= c.regiao.fim <= max_fim_f3
                ]
                cache_f3[ini_f2] = (
                    min(validos_f3, key=lambda c: (c.desvio_tm, c.regiao.ini, c.regiao.fim))
                    if validos_f3
                    else None
                )

    # Cache do melhor B3c por fim(B2c)
    cache_b3c: dict[int, Candidato | None] = {}
    for c_b2c in b2c_list:
        fim_b2c = c_b2c.regiao.fim
        if fim_b2c not in cache_b3c:
            min_ini_b3c = fim_b2c + D["B2c_B3c_gap"][0]
            max_ini_b3c = fim_b2c + D["B2c_B3c_gap"][1]
            if min_ini_b3c > max_ini_b3c:
                cache_b3c[fim_b2c] = None
            else:
                idx_start = bisect_left(b3c_inis, min_ini_b3c)
                idx_end = bisect_right(b3c_inis, max_ini_b3c)
                validos_b3c = b3c_list[idx_start:idx_end]
                cache_b3c[fim_b2c] = (
                    min(validos_b3c, key=lambda c: (c.desvio_tm, c.regiao.ini, c.regiao.fim))
                    if validos_b3c
                    else None
                )

    try:
        for c_f2 in f2_list:
            ini_f2 = c_f2.regiao.ini
            c_f3 = cache_f3[ini_f2]
            if c_f3 is None:
                continue

            # F1 com ini(F1) - ini(F2) em D["F2_F1"]
            min_ini_f1 = ini_f2 + D["F2_F1"][0]
            max_ini_f1 = ini_f2 + D["F2_F1"][1]
            idx_f1_start = bisect_left(f1_inis, min_ini_f1)
            idx_f1_end = bisect_right(f1_inis, max_ini_f1)

            for c_f1 in f1_list[idx_f1_start:idx_f1_end]:
                fim_f1 = c_f1.regiao.fim

                # B1c com ini(B1c) - fim(F1) em D["F1_B1c_gap"] e
                # fim(B1c) <= ini(F2) + D["F2_B2"][1] - D["B1c_B2c"][0] (poda)
                min_ini_b1c = fim_f1 + D["F1_B1c_gap"][0]
                max_ini_b1c = fim_f1 + D["F1_B1c_gap"][1]
                idx_b1c_start = bisect_left(b1c_inis, min_ini_b1c)
                idx_b1c_end = bisect_right(b1c_inis, max_ini_b1c)
                max_fim_b1c = ini_f2 + D["F2_B2"][1] - D["B1c_B2c"][0]

                for c_b1c in b1c_list[idx_b1c_start:idx_b1c_end]:
                    if c_b1c.regiao.fim > max_fim_b1c:
                        continue
                    fim_b1c = c_b1c.regiao.fim

                    # B2c com fim(B2c) - fim(B1c) em D["B1c_B2c"] e fim(B2c) - ini(F2) em D["F2_B2"]
                    min_fim_b2c = max(fim_b1c + D["B1c_B2c"][0], ini_f2 + D["F2_B2"][0])
                    max_fim_b2c = min(fim_b1c + D["B1c_B2c"][1], ini_f2 + D["F2_B2"][1])
                    if min_fim_b2c > max_fim_b2c:
                        continue

                    idx_b2c_start = bisect_left(b2c_inis, min_fim_b2c - L_max)
                    idx_b2c_end = bisect_right(b2c_inis, max_fim_b2c - L_min)

                    for c_b2c in b2c_list[idx_b2c_start:idx_b2c_end]:
                        if not (min_fim_b2c <= c_b2c.regiao.fim <= max_fim_b2c):
                            continue
                        fim_b2c = c_b2c.regiao.fim

                        c_b3c = cache_b3c[fim_b2c]
                        if c_b3c is None:
                            continue

                        combinacoes.append((c_f3, c_f2, c_f1, c_b1c, c_b2c, c_b3c))
                        conjuntos_examinados += 1
                        if conjuntos_examinados >= max_examinados:
                            busca_truncada = True
                            raise StopIteration
    except StopIteration:
        pass

    stats = {
        "candidatos_por_tipo": {t: len(cands.get(t, [])) for t in ("F3", "F2", "F1", "B1c", "B2c", "B3c")},
        "conjuntos_examinados": conjuntos_examinados,
        "busca_truncada": busca_truncada,
    }

    if not combinacoes:
        return ResultadoDesenho(conjuntos=[], estatisticas=stats)

    # 3. Ranking em duas etapas
    # Etapa 1: desvio_tm total e chave_coord
    candidatos_etapa1 = []
    for combo in combinacoes:
        c_f3, c_f2, c_f1, c_b1c, c_b2c, c_b3c = combo
        desvio_tm_total = (
            c_f3.desvio_tm
            + c_f2.desvio_tm
            + c_f1.desvio_tm
            + c_b1c.desvio_tm
            + c_b2c.desvio_tm
            + c_b3c.desvio_tm
        )
        chave_coord = (
            c_f3.regiao.ini,
            c_f2.regiao.ini,
            c_f1.regiao.ini,
            c_b1c.regiao.ini,
            c_b2c.regiao.ini,
            c_b3c.regiao.ini,
            c_f3.regiao.fim,
            c_f2.regiao.fim,
            c_f1.regiao.fim,
            c_b1c.regiao.fim,
            c_b2c.regiao.fim,
            c_b3c.regiao.fim,
        )
        candidatos_etapa1.append((desvio_tm_total, chave_coord, combo))

    candidatos_etapa1.sort(key=lambda x: (x[0], x[1]))
    pre_selecao_k = int(cfg["desenho"]["pre_selecao"])
    pre_selecionados = candidatos_etapa1[:pre_selecao_k]

    # Etapa 2: heterodímeros e pontuação final
    cond = Condicoes.de_config(cfg)
    ligante = cfg["desenho"]["ligante"]
    peso_dimero = float(cfg["desenho"]["peso_dimero"])
    candidatos_etapa2 = []

    for desvio_tm_total, chave_coord, combo in pre_selecionados:
        c_f3, c_f2, c_f1, c_b1c, c_b2c, c_b3c = combo
        regioes_base = {
            "F3": c_f3.regiao,
            "F2": c_f2.regiao,
            "F1": c_f1.regiao,
            "B1c": c_b1c.regiao,
            "B2c": c_b2c.regiao,
            "B3c": c_b3c.regiao,
        }
        dists = calcular_distancias(regioes_base)
        if any(d_nome in dists and not (D[d_nome][0] <= dists[d_nome] <= D[d_nome][1]) for d_nome in D):
            continue
        conjunto_base = primers_de_regioes(alvo, regioes_base, ligante=ligante)
        primers_tubo_4 = [conjunto_base.F3, conjunto_base.B3, conjunto_base.FIP, conjunto_base.BIP]
        dgs = []
        for i in range(len(primers_tubo_4)):
            for j in range(i + 1, len(primers_tubo_4)):
                het = heterodimero(primers_tubo_4[i], primers_tubo_4[j], cond)
                dgs.append(het.dg_kcal)
        pior_het = min(dgs)
        pontuacao = desvio_tm_total + peso_dimero * max(0.0, -pior_het)

        candidatos_etapa2.append(
            {
                "pontuacao": pontuacao,
                "desvio_tm": desvio_tm_total,
                "pior_heterodimero_kcal": pior_het,
                "chave_coord": chave_coord,
                "combo": combo,
                "regioes_base": regioes_base,
            }
        )

    candidatos_etapa2.sort(key=lambda x: (x["pontuacao"], x["chave_coord"]))
    limite_n = int(cfg["desenho"]["n_resultados"]) if n is None else n
    finais = candidatos_etapa2[:limite_n]

    # 4. Loops e 5. Autoverificação
    conjuntos_desenhados: list[ConjuntoDesenhado] = []
    for idx, item in enumerate(finais, start=1):
        c_f3, c_f2, c_f1, c_b1c, c_b2c, c_b3c = item["combo"]
        regioes = dict(item["regioes_base"])
        avisos: list[str] = []

        # LF: dentro de [fim(F2), ini(F1))
        fim_f2 = c_f2.regiao.fim
        ini_f1 = c_f1.regiao.ini
        if fim_f2 < ini_f1:
            dentro_lf = Regiao(fim_f2, ini_f1)
            cands_lf = gerar_candidatos(alvo, cfg, tipos=("LF",), dentro=dentro_lf).get("LF", [])
            if cands_lf:
                best_lf = min(cands_lf, key=lambda c: (c.desvio_tm, c.regiao.ini, c.regiao.fim))
                regioes["LF"] = best_lf.regiao
            else:
                avisos.append("sem LF candidato")
        else:
            avisos.append("sem LF candidato")

        # LB: dentro de [fim(B1c), ini(B2c))
        fim_b1c = c_b1c.regiao.fim
        ini_b2c = c_b2c.regiao.ini
        if fim_b1c < ini_b2c:
            dentro_lb = Regiao(fim_b1c, ini_b2c)
            cands_lb = gerar_candidatos(alvo, cfg, tipos=("LB",), dentro=dentro_lb).get("LB", [])
            if cands_lb:
                best_lb = min(cands_lb, key=lambda c: (c.desvio_tm, c.regiao.ini, c.regiao.fim))
                regioes["LB"] = best_lb.regiao
            else:
                avisos.append("sem LB candidato")
        else:
            avisos.append("sem LB candidato")

        conjunto_final = primers_de_regioes(alvo, regioes, ligante=ligante)

        # 5. Autoverificação
        violacoes = verificar(conjunto_final, alvo, cfg)
        if violacoes:
            raise RuntimeError(f"bug no desenho: {violacoes}")

        conjuntos_desenhados.append(
            ConjuntoDesenhado(
                posicao=idx,
                pontuacao=item["pontuacao"],
                desvio_tm=item["desvio_tm"],
                pior_heterodimero_kcal=item["pior_heterodimero_kcal"],
                conjunto=conjunto_final,
                regioes=regioes,
                avisos=avisos,
            )
        )

    return ResultadoDesenho(conjuntos=conjuntos_desenhados, estatisticas=stats)
