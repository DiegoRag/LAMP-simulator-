"""Módulo de verificação de restrições de desenho para primers LAMP."""

from dataclasses import dataclass

from primers.casamento import (
    REGIAO_PARA_PRIMER,
    Posicao,
    distancia_ponta_extensao,
    janela_esperada,
    melhor_posicao,
    ocorrencias_compativeis,
    ponta_extensao,
    regioes_do_conjunto,
)
from primers.modelo import ConjuntoLAMP, Regiao, normalizar, rc
from primers.termo import CondicoesReferencia, dg_extremidade, gc_pct, tm_referencia

ORDEM_VERIFICACAO = ["F3", "B3", "F2", "F1c", "B2", "B1c", "LF", "LB"]
JANELA_CADEIA_NT_PADRAO: int = 1000


@dataclass(frozen=True)
class Violacao:
    """Representa uma violação de restrição de desenho."""

    item: str
    regra: str
    valor: float | None
    faixa: str
    mensagem: str


@dataclass(frozen=True)
class Localizacao:
    """Resultado detalhado da localização das regiões de um conjunto LAMP no alvo."""

    fita: str  # "+" (alvo original) ou "-" (rc(alvo))
    regioes: dict[str, Regiao]  # regiões com casamento exato na fita escolhida
    aproximadas: dict[str, Posicao]  # regiões não exatas localizadas com <= max_erros na janela
    regioes_completas: dict[str, Regiao]  # união de regioes exatas + regiao de cada aproximada
    violacoes: list[Violacao]  # violações nao_encontrado
    avisos: list[str]  # avisos de fita "-", múltiplas ocorrências, posição aproximada e origem circular
    cruza_origem: bool = False  # True quando circular=True e r.fim > len(alvo)


def _verificar_partes_individuais(primers: dict[str, str], cfg: dict) -> list[Violacao]:
    """Aplica comprimento, GC, Tm de referência e extremidade aos primers individuais presentes."""
    ref = CondicoesReferencia.de_config(cfg)
    violacoes: list[Violacao] = []

    for nome in ORDEM_VERIFICACAO:
        if nome not in primers:
            continue

        seq = primers[nome]

        # 1. Comprimento
        if nome in ("LF", "LB"):
            c_min = int(cfg["comprimento"]["loop_min"])
            c_max = int(cfg["comprimento"]["loop_max"])
        else:
            c_min = int(cfg["comprimento"]["min"])
            c_max = int(cfg["comprimento"]["max"])

        c = len(seq)
        if c < c_min or c > c_max:
            faixa_c = f"[{c_min}, {c_max}]"
            violacoes.append(
                Violacao(
                    item=nome,
                    regra="comprimento",
                    valor=float(c),
                    faixa=faixa_c,
                    mensagem=f"Comprimento {c} fora da faixa {faixa_c}",
                )
            )

        # 2. GC%
        gc = gc_pct(seq)
        gc_min = float(cfg["gc_pct"]["min"])
        gc_max = float(cfg["gc_pct"]["max"])
        if gc < gc_min or gc > gc_max:
            faixa_gc = f"[{gc_min:.1f}, {gc_max:.1f}]"
            violacoes.append(
                Violacao(
                    item=nome,
                    regra="gc_pct",
                    valor=float(gc),
                    faixa=faixa_gc,
                    mensagem=f"GC {gc:.2f}% fora da faixa {faixa_gc}",
                )
            )

        # 3. Tm de referência (convenção do PrimerExplorer V5, a das faixas de [tm_c])
        if nome in ("F3", "B3"):
            tm_cfg = cfg["tm_c"]["F3_B3"]
        elif nome in ("F2", "B2"):
            tm_cfg = cfg["tm_c"]["F2_B2"]
        elif nome in ("F1c", "B1c"):
            tm_cfg = cfg["tm_c"]["F1c_B1c"]
        else:
            tm_cfg = cfg["tm_c"]["loops"]

        tm_val = tm_referencia(seq, ref)
        tm_min, tm_max = float(tm_cfg[0]), float(tm_cfg[1])
        if tm_val < tm_min or tm_val > tm_max:
            faixa_tm = f"[{tm_min:.1f}, {tm_max:.1f}]"
            violacoes.append(
                Violacao(
                    item=nome,
                    regra="tm_c",
                    valor=float(tm_val),
                    faixa=faixa_tm,
                    mensagem=f"Tm {tm_val:.2f}°C fora da faixa {faixa_tm}",
                )
            )

        # 4. Extremidade
        lado = "5" if nome in ("F1c", "B1c") else "3"
        bases = int(cfg["extremidade"]["bases"])
        dg_max = float(cfg["extremidade"]["dg_max"])
        dg = dg_extremidade(seq, lado, bases, ref)
        if dg > dg_max:
            faixa_dg = f"<= {dg_max:.1f}"
            violacoes.append(
                Violacao(
                    item=nome,
                    regra="extremidade",
                    valor=float(dg),
                    faixa=faixa_dg,
                    mensagem=(
                        f"ΔG da extremidade {lado}' ({dg:.2f} kcal/mol) "
                        f"maior que o máximo ({dg_max:.1f} kcal/mol)"
                    ),
                )
            )

    return violacoes


def avaliar_primer_quimera(
    nome_int: str,
    seq_int: str,
    cfg: dict,
) -> tuple[dict[str, str] | None, list[Violacao]]:
    """Avalia FIP ou BIP não separado testando todos os cortes plausíveis nas partes."""
    seq_norm = normalizar(seq_int)
    nome_5, nome_3 = ("F1c", "F2") if nome_int == "FIP" else ("B1c", "B2")
    c_min = int(cfg["comprimento"]["min"])
    c_max = int(cfg["comprimento"]["max"])
    bases_ext = int(cfg["extremidade"]["bases"])
    ref = CondicoesReferencia.de_config(cfg)
    tm_5_cfg = cfg["tm_c"]["F1c_B1c"]
    tm_3_cfg = cfg["tm_c"]["F2_B2"]
    centro_5 = (float(tm_5_cfg[0]) + float(tm_5_cfg[1])) / 2.0
    centro_3 = (float(tm_3_cfg[0]) + float(tm_3_cfg[1])) / 2.0
    lig_padrao = str(cfg.get("desenho", {}).get("ligante", ""))
    ligs_cfg = cfg.get("desenho", {}).get("ligantes_candidatos", ["", "TTTT", "TTT", "TTTTT"])
    cands_lig = (lig_padrao, *[str(x) for x in ligs_cfg if str(x) != lig_padrao])
    max_lig = max((len(x) for x in cands_lig), default=0)

    cortes: list[tuple[tuple, dict[str, str], list[Violacao]]] = []
    for lig in cands_lig:
        len_lig = len(lig)
        for i in range(c_min, len(seq_norm) - c_min - len_lig + 1):
            if len_lig > 0 and seq_norm[i : i + len_lig] != lig:
                continue
            p5 = seq_norm[:i]
            p3 = seq_norm[i + len_lig :]
            if c_min <= len(p5) <= c_max and c_min <= len(p3) <= c_max:
                v_corte = _verificar_partes_individuais({nome_5: p5, nome_3: p3}, cfg)
                desvio_tm = abs(tm_referencia(p5, ref) - centro_5) + abs(
                    tm_referencia(p3, ref) - centro_3
                )
                pref_lig = 0 if lig == lig_padrao else (1 if lig == "" else 2)
                desvio_comp = abs(len(p5) - len(p3))
                chave = (len(v_corte), desvio_tm, pref_lig, desvio_comp, i)
                cortes.append((chave, {nome_5: p5, nome_3: p3, "ligante": lig}, v_corte))

    if cortes:
        _, melhores_partes, melhores_viols = min(cortes, key=lambda x: x[0])
        return melhores_partes, melhores_viols

    faixa_int = f"[{2 * c_min}, {2 * c_max + max_lig}]"
    viols: list[Violacao] = [
        Violacao(
            item=nome_int,
            regra="comprimento",
            valor=float(len(seq_norm)),
            faixa=faixa_int,
            mensagem=(
                f"Comprimento {len(seq_norm)} fora da faixa {faixa_int} "
                f"para separar {nome_int} em {nome_5} e {nome_3}"
            ),
        )
    ]
    if len(seq_norm) >= 2 * bases_ext:
        meio = len(seq_norm) // 2
        p5 = seq_norm[:meio]
        p3 = seq_norm[meio:]
        v_meio = _verificar_partes_individuais({nome_5: p5, nome_3: p3}, cfg)
        viols.extend([v for v in v_meio if v.regra != "comprimento"])

    return None, viols


def verificar_primers(primers: dict[str, str], cfg: dict) -> list[Violacao]:
    """Aplica comprimento, GC, Tm e extremidade aos primers presentes."""
    violacoes = _verificar_partes_individuais(primers, cfg)

    if "FIP" in primers and ("F1c" not in primers or "F2" not in primers):
        _, viols_fip = avaliar_primer_quimera("FIP", primers["FIP"], cfg)
        violacoes.extend(viols_fip)

    if "BIP" in primers and ("B1c" not in primers or "B2" not in primers):
        _, viols_bip = avaliar_primer_quimera("BIP", primers["BIP"], cfg)
        violacoes.extend(viols_bip)

    return violacoes


def calcular_distancias(regioes: dict[str, Regiao]) -> dict[str, int]:
    """Calcula as 7 distâncias: F2_B2, F2_F1, B1c_B2c, F3_F2_gap, B2c_B3c_gap, F1_B1c_gap, F3_B3_total.

    Devolve apenas as distâncias cujas regiões necessárias estão todas presentes.
    """
    dist: dict[str, int] = {}
    if "B2c" in regioes and "F2" in regioes:
        dist["F2_B2"] = regioes["B2c"].fim - regioes["F2"].ini
    if "F1" in regioes and "F2" in regioes:
        dist["F2_F1"] = regioes["F1"].ini - regioes["F2"].ini
    if "B2c" in regioes and "B1c" in regioes:
        dist["B1c_B2c"] = regioes["B2c"].fim - regioes["B1c"].fim
    if "F2" in regioes and "F3" in regioes:
        dist["F3_F2_gap"] = regioes["F2"].ini - regioes["F3"].fim
    if "B3c" in regioes and "B2c" in regioes:
        dist["B2c_B3c_gap"] = regioes["B3c"].ini - regioes["B2c"].fim
    if "B1c" in regioes and "F1" in regioes:
        dist["F1_B1c_gap"] = regioes["B1c"].ini - regioes["F1"].fim
    if "B3c" in regioes and "F3" in regioes:
        dist["F3_B3_total"] = regioes["B3c"].fim - regioes["F3"].ini
    return dist


def _todas_ocorrencias(texto: str, sub: str, limite_ini: int | None = None) -> list[Regiao]:
    """Retorna todas as ocorrências compatíveis de sub em texto com pos < limite_ini."""
    return ocorrencias_compativeis(texto, sub, limite_ini=limite_ini)


def _tam_lis(seq: list[int]) -> int:
    """Comprimento da maior subsequência estritamente crescente."""
    if not seq:
        return 0
    dp = [1] * len(seq)
    for i in range(len(seq)):
        for j in range(i):
            if seq[j] < seq[i]:
                dp[i] = max(dp[i], dp[j] + 1)
    return max(dp)


def localizar_detalhado(
    conjunto: ConjuntoLAMP,
    alvo: str,
    cfg: dict | None = None,
    *,
    circular: bool | None = None,
) -> Localizacao:
    """Resultado detalhado da localização das regiões de um conjunto LAMP no alvo."""
    eh_circular = (
        circular if circular is not None else bool((cfg or {}).get("alvo", {}).get("circular", False))
    )
    casamento_cfg = (cfg or {}).get("casamento") or {}
    max_erros = int(casamento_cfg.get("max_erros", 3))
    janela = int(casamento_cfg.get("janela_cadeia_nt", JANELA_CADEIA_NT_PADRAO))
    distancias_cfg = (cfg or {}).get(
        "distancias",
        {
            "F2_B2": [120, 180],
            "F2_F1": [40, 60],
            "B1c_B2c": [40, 60],
            "F3_F2_gap": [0, 60],
            "B2c_B3c_gap": [0, 60],
            "F1_B1c_gap": [0, 1000],
            "F3_B3_total": [120, 280],
        },
    )

    mapa_regioes = regioes_do_conjunto(conjunto)
    mapa_nomes = [nome for nome, _ in mapa_regioes]
    essenciais = ("F3", "F2", "F1", "B1c", "B2c", "B3c")

    N = len(alvo)
    alvo_rc = rc(alvo)  # uma vez só: em genoma inteiro, cada rc custa uma passada no alvo
    mult_ocorrencias: dict[str, tuple[int, int]] = {}
    for reg_nome, seq_reg in mapa_regioes:
        if eh_circular and N > 0:
            n_m = len(_todas_ocorrencias(alvo + alvo, seq_reg, limite_ini=N))
            n_rc = len(_todas_ocorrencias(alvo_rc + alvo_rc, seq_reg, limite_ini=N))
        else:
            n_m = len(_todas_ocorrencias(alvo, seq_reg, limite_ini=N))
            n_rc = len(_todas_ocorrencias(alvo_rc, seq_reg, limite_ini=N))
        mult_ocorrencias[reg_nome] = (n_m, n_rc)

    resultados_fita: list[
        tuple[tuple, str, dict[str, Regiao], dict[str, Posicao], dict[str, Regiao]]
    ] = []

    for fita_id, seq_fita in (("+", alvo), ("-", alvo_rc)):
        texto_busca = (seq_fita + seq_fita) if (eh_circular and N > 0) else seq_fita
        ocorr: dict[str, list[Regiao]] = {}
        for reg_nome, seq_reg in mapa_regioes:
            ocorr[reg_nome] = _todas_ocorrencias(texto_busca, seq_reg, limite_ini=N)

        presentes = [nome for nome in essenciais if ocorr.get(nome)]
        if not presentes:
            regioes_exatas: dict[str, Regiao] = {}
        else:
            todas_cadeias: list[dict[str, Regiao]] = []
            for start_idx in range(len(presentes)):
                nome_anc = presentes[start_idx]
                for r_anc in ocorr[nome_anc]:
                    cands_por_m: list[tuple[str, list[Regiao]]] = []
                    for m in range(start_idx, len(presentes)):
                        nm = presentes[m]
                        if m == start_idx:
                            cands_m = [r_anc]
                        else:
                            cands_m = [
                                r for r in ocorr[nm]
                                if r.ini >= r_anc.ini and r.ini - r_anc.ini <= janela
                            ]
                            if eh_circular:
                                cands_m.extend(
                                    Regiao(r.ini + N, r.fim + N)
                                    for r in ocorr[nm]
                                    if r.ini < r_anc.ini
                                    and (r.ini + N) - r_anc.ini <= janela
                                    and (r.fim + N) - r_anc.ini <= N
                                )
                        cands_por_m.append((nm, cands_m))

                    def backtrack(
                        idx_m: int,
                        last_fim: int,
                        atual: list[tuple[str, Regiao]],
                        cands_por_m=cands_por_m,
                        todas_cadeias=todas_cadeias,
                    ):
                        if idx_m == len(cands_por_m):
                            if atual:
                                todas_cadeias.append(dict(atual))
                            return
                        reg_m_nome, reg_m_cands = cands_por_m[idx_m]
                        for r in reg_m_cands:
                            if r.ini >= last_fim:
                                atual.append((reg_m_nome, r))
                                backtrack(idx_m + 1, r.fim, atual)
                                atual.pop()
                        backtrack(idx_m + 1, last_fim, atual)

                    backtrack(0, 0, [])

            if todas_cadeias:

                def chave_combo(combo: dict[str, Regiao]):
                    dists = calcular_distancias(combo)
                    n_viol_dist = sum(
                        1
                        for d_k, (d_min, d_max) in distancias_cfg.items()
                        if d_k in dists and not (d_min <= dists[d_k] <= d_max)
                    )
                    pen_dist = sum(
                        max(0, d_min - dists[d_k], dists[d_k] - d_max)
                        for d_k, (d_min, d_max) in distancias_cfg.items()
                        if d_k in dists
                    )
                    span = (
                        (max(r.fim for r in combo.values()) - min(r.ini for r in combo.values()))
                        if combo
                        else 0
                    )
                    coords_key = tuple(
                        (combo[k].ini, combo[k].fim) for k in essenciais if k in combo
                    )
                    return (-len(combo), n_viol_dist, pen_dist, span, coords_key)

                melhor_combo = min(todas_cadeias, key=chave_combo)
                regioes_exatas = dict(melhor_combo)
            else:
                regioes_exatas = {}

            for nome in presentes:
                if nome not in regioes_exatas and ocorr[nome]:
                    regioes_exatas[nome] = ocorr[nome][0]

        # Seleção de LF e LB
        if "LF" in mapa_nomes and ocorr.get("LF"):
            cands_lf = list(ocorr["LF"])
            if eh_circular and any(r.fim > N for r in regioes_exatas.values()):
                cands_lf.extend(
                    Regiao(r.ini + N, r.fim + N)
                    for r in ocorr["LF"]
                    if r.fim + N <= 2 * N
                )
            if "F2" in regioes_exatas and "F1" in regioes_exatas:
                ini_gap = regioes_exatas["F2"].fim
                fim_gap = regioes_exatas["F1"].ini

                def score_lf(r: Regiao, ini_gap=ini_gap, fim_gap=fim_gap):
                    contida = 0 if (r.ini >= ini_gap and r.fim <= fim_gap) else 1
                    dist = max(0, ini_gap - r.ini) + max(0, r.fim - fim_gap)
                    return (contida, dist, r.ini, r.fim)

                regioes_exatas["LF"] = min(cands_lf, key=score_lf)
            else:
                regioes_exatas["LF"] = min(cands_lf, key=lambda r: (r.ini, r.fim))

        if "LB" in mapa_nomes and ocorr.get("LB"):
            cands_lb = list(ocorr["LB"])
            if eh_circular and any(r.fim > N for r in regioes_exatas.values()):
                cands_lb.extend(
                    Regiao(r.ini + N, r.fim + N)
                    for r in ocorr["LB"]
                    if r.fim + N <= 2 * N
                )
            if "B1c" in regioes_exatas and "B2c" in regioes_exatas:
                ini_gap = regioes_exatas["B1c"].fim
                fim_gap = regioes_exatas["B2c"].ini

                def score_lb(r: Regiao, ini_gap=ini_gap, fim_gap=fim_gap):
                    contida = 0 if (r.ini >= ini_gap and r.fim <= fim_gap) else 1
                    dist = max(0, ini_gap - r.ini) + max(0, r.fim - fim_gap)
                    return (contida, dist, r.ini, r.fim)

                regioes_exatas["LB"] = min(cands_lb, key=score_lb)
            else:
                regioes_exatas["LB"] = min(cands_lb, key=lambda r: (r.ini, r.fim))

        # Busca aproximada para regiões ausentes
        tam_busca = len(texto_busca) if eh_circular else N
        aproximadas_fita: dict[str, Posicao] = {}
        for reg_nome, seq_reg in mapa_regioes:
            if reg_nome not in regioes_exatas:
                jan = janela_esperada(reg_nome, regioes_exatas, tam_busca, cfg)
                mp = melhor_posicao(
                    seq_reg,
                    texto_busca[:tam_busca],
                    max_erros=max_erros,
                    janela=jan,
                    ambas_fitas=False,
                )
                if mp is not None:
                    aproximadas_fita[reg_nome] = Posicao(
                        regiao=mp.regiao, diferencas=mp.diferencas, fita=fita_id
                    )

        regioes_completas = dict(regioes_exatas) | {k: v.regiao for k, v in aproximadas_fita.items()}

        # Pontuação da fita
        inis_essenciais = [regioes_completas[k].ini for k in essenciais if k in regioes_completas]
        n_essenciais_ord = _tam_lis(inis_essenciais)
        dists_comp = calcular_distancias(regioes_completas)
        pen_dist_completas = sum(
            max(0, d_min - dists_comp[d_k], dists_comp[d_k] - d_max)
            for d_k, (d_min, d_max) in distancias_cfg.items()
            if d_k in dists_comp
        )
        erros_aprox = sum(len(p.diferencas) for p in aproximadas_fita.values())
        desempate = 1 if fita_id == "+" else 0

        score_fita = (
            n_essenciais_ord,
            len(regioes_exatas),
            len(regioes_completas),
            -erros_aprox,
            -pen_dist_completas,
            desempate,
        )
        resultados_fita.append(
            (score_fita, fita_id, regioes_exatas, aproximadas_fita, regioes_completas)
        )

    # Escolha da melhor fita
    resultados_fita.sort(key=lambda x: x[0], reverse=True)
    _, fita_esc, regioes_esc, aproximadas_esc, completas_esc = resultados_fita[0]

    cruza_origem = bool(eh_circular and any(r.fim > N for r in completas_esc.values()))

    avisos: list[str] = []
    if fita_esc == "-":
        avisos.append(
            "Conjunto localizado na fita complementar (−); "
            "coordenadas referem-se a rc(alvo) no sentido 5'→3'."
        )
    if cruza_origem:
        avisos.append("Amplicon cruza a origem do alvo circular.")
    for reg_nome, _ in mapa_regioes:
        n_m, n_rc = mult_ocorrencias.get(reg_nome, (0, 0))
        if n_m + n_rc > 1:
            avisos.append(
                f"Região {reg_nome} possui {n_m + n_rc} ocorrências no alvo "
                f"(fita +: {n_m}, fita −: {n_rc}); foi escolhida a ocorrência coerente com o conjunto."
            )
    if aproximadas_esc:
        nomes_aprox = ", ".join(aproximadas_esc.keys())
        avisos.append(f"Posição aproximada usada para calcular distâncias de: {nomes_aprox}.")

    violacoes: list[Violacao] = []
    for reg_nome, _ in mapa_regioes:
        if reg_nome not in regioes_esc:
            if reg_nome in aproximadas_esc:
                mp = aproximadas_esc[reg_nome]
                ponta = f"ponta {ponta_extensao(reg_nome)} de {REGIAO_PARA_PRIMER[reg_nome]}"
                det = "; ".join(
                    f"posição {d.posicao + 1}: {d.esperado} esperado, {d.no_alvo} no alvo "
                    f"(a {distancia_ponta_extensao(reg_nome, mp.regiao, d.posicao)} nt da {ponta})"
                    for d in mp.diferencas
                )
                mensagem = (
                    f"Região {reg_nome} não encontrada exata no alvo; "
                    f"melhor posição aproximada {mp.regiao.ini + 1}–{mp.regiao.fim} ({det})"
                )
                violacoes.append(
                    Violacao(
                        item=reg_nome,
                        regra="nao_encontrado",
                        valor=None,
                        faixa="",
                        mensagem=mensagem,
                    )
                )
            else:
                mensagem = f"Região {reg_nome} não encontrada no alvo"
                violacoes.append(
                    Violacao(
                        item=reg_nome,
                        regra="nao_encontrado",
                        valor=None,
                        faixa="",
                        mensagem=mensagem,
                    )
                )

    return Localizacao(
        fita=fita_esc,
        regioes=regioes_esc,
        aproximadas=aproximadas_esc,
        regioes_completas=completas_esc,
        violacoes=violacoes,
        avisos=avisos,
        cruza_origem=cruza_origem,
    )


def localizar(
    conjunto: ConjuntoLAMP,
    alvo: str,
    cfg: dict | None = None,
    *,
    circular: bool | None = None,
) -> tuple[dict[str, Regiao], list[Violacao]]:
    """Mantém retrocompatibilidade retornando `(loc.regioes, loc.violacoes)`."""
    loc = localizar_detalhado(conjunto, alvo, cfg, circular=circular)
    return loc.regioes, loc.violacoes


def verificar(
    conjunto: ConjuntoLAMP,
    alvo: str | None,
    cfg: dict,
    *,
    circular: bool | None = None,
    loc: Localizacao | None = None,
) -> list[Violacao]:
    """Verifica restrições de desenho para o conjunto e, se informado o alvo, distâncias e loops."""
    primers_dict = {
        "F3": conjunto.F3,
        "B3": conjunto.B3,
        "F2": conjunto.F2,
        "F1c": conjunto.F1c,
        "B2": conjunto.B2,
        "B1c": conjunto.B1c,
    }
    if conjunto.LF is not None:
        primers_dict["LF"] = conjunto.LF
    if conjunto.LB is not None:
        primers_dict["LB"] = conjunto.LB

    violacoes = verificar_primers(primers_dict, cfg)

    if alvo is None:
        return violacoes

    loc_efetivo = loc if loc is not None else localizar_detalhado(conjunto, alvo, cfg, circular=circular)
    violacoes.extend(loc_efetivo.violacoes)

    regioes_verif = loc_efetivo.regioes_completas

    # Distâncias na ordem do arquivo de configuração
    distancias_cfg = cfg.get("distancias", {})
    dist_calculadas = calcular_distancias(regioes_verif)
    for dist_nome, faixa in distancias_cfg.items():
        if dist_nome in dist_calculadas:
            val = dist_calculadas[dist_nome]
            d_min, d_max = float(faixa[0]), float(faixa[1])
            if val < d_min or val > d_max:
                faixa_str = f"[{int(d_min)}, {int(d_max)}]"
                violacoes.append(
                    Violacao(
                        item=dist_nome,
                        regra="distancia",
                        valor=float(val),
                        faixa=faixa_str,
                        mensagem=f"Distância {dist_nome} ({val}) fora da faixa {faixa_str}",
                    )
                )

    # Verificação de loops
    if (
        conjunto.LF is not None
        and "LF" in regioes_verif
        and "F2" in regioes_verif
        and "F1" in regioes_verif
    ):
        fim_f2 = regioes_verif["F2"].fim
        ini_f1 = regioes_verif["F1"].ini
        r_lf = regioes_verif["LF"]
        if not (r_lf.ini >= fim_f2 and r_lf.fim <= ini_f1):
            faixa_str = f"[{fim_f2}, {ini_f1})"
            violacoes.append(
                Violacao(
                    item="LF",
                    regra="loop_fora_da_regiao",
                    valor=float(r_lf.ini),
                    faixa=faixa_str,
                    mensagem=f"LF [{r_lf.ini}, {r_lf.fim}) fora da região [{fim_f2}, {ini_f1})",
                )
            )

    if (
        conjunto.LB is not None
        and "LB" in regioes_verif
        and "B1c" in regioes_verif
        and "B2c" in regioes_verif
    ):
        fim_b1c = regioes_verif["B1c"].fim
        ini_b2c = regioes_verif["B2c"].ini
        r_lb = regioes_verif["LB"]
        if not (r_lb.ini >= fim_b1c and r_lb.fim <= ini_b2c):
            faixa_str = f"[{fim_b1c}, {ini_b2c})"
            violacoes.append(
                Violacao(
                    item="LB",
                    regra="loop_fora_da_regiao",
                    valor=float(r_lb.ini),
                    faixa=faixa_str,
                    mensagem=f"LB [{r_lb.ini}, {r_lb.fim}) fora da região [{fim_b1c}, {ini_b2c})",
                )
            )

    return violacoes
