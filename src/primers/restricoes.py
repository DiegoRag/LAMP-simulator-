"""Módulo de verificação de restrições de desenho para primers LAMP."""

from dataclasses import dataclass
from primers.modelo import ConjuntoLAMP, Regiao, rc
from primers.termo import Condicoes, dg_extremidade, gc_pct, tm

ORDEM_VERIFICACAO = ["F3", "B3", "F2", "F1c", "B2", "B1c", "LF", "LB"]


@dataclass(frozen=True)
class Violacao:
    """Representa uma violação de restrição de desenho."""

    item: str
    regra: str
    valor: float | None
    faixa: str
    mensagem: str


def verificar_primers(primers: dict[str, str], cfg: dict) -> list[Violacao]:
    """Aplica comprimento, GC, Tm e extremidade aos primers presentes.

    Ignora FIP e BIP.
    """
    cond = Condicoes.de_config(cfg)
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

        # 3. Tm
        if nome in ("F3", "B3"):
            conc = float(cfg["concentracao_nM"]["externos"])
            tm_cfg = cfg["tm_c"]["F3_B3"]
        elif nome in ("F2", "B2"):
            conc = float(cfg["concentracao_nM"]["internos"])
            tm_cfg = cfg["tm_c"]["F2_B2"]
        elif nome in ("F1c", "B1c"):
            conc = float(cfg["concentracao_nM"]["internos"])
            tm_cfg = cfg["tm_c"]["F1c_B1c"]
        else:
            conc = float(cfg["concentracao_nM"]["loops"])
            tm_cfg = cfg["tm_c"]["loops"]

        tm_val = tm(seq, cond, conc)
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
        dg = dg_extremidade(seq, lado, bases, cond)
        if dg > dg_max:
            faixa_dg = f"<= {dg_max:.1f}"
            violacoes.append(
                Violacao(
                    item=nome,
                    regra="extremidade",
                    valor=float(dg),
                    faixa=faixa_dg,
                    mensagem=f"ΔG da extremidade {lado}' ({dg:.2f} kcal/mol) maior que o máximo ({dg_max:.1f} kcal/mol)",
                )
            )

    return violacoes


def localizar(conjunto: ConjuntoLAMP, alvo: str) -> tuple[dict[str, Regiao], list[Violacao]]:
    """Localiza regiões do conjunto no alvo e identifica as não encontradas."""
    mapa_regioes: list[tuple[str, str]] = [
        ("F3", conjunto.F3),
        ("F2", conjunto.F2),
        ("F1", rc(conjunto.F1c)),
        ("B1c", conjunto.B1c),
        ("B2c", rc(conjunto.B2)),
        ("B3c", rc(conjunto.B3)),
    ]
    if conjunto.LF is not None:
        mapa_regioes.append(("LF", rc(conjunto.LF)))
    if conjunto.LB is not None:
        mapa_regioes.append(("LB", conjunto.LB))

    regioes: dict[str, Regiao] = {}
    violacoes: list[Violacao] = []

    for reg_nome, seq_fita_mais in mapa_regioes:
        pos = alvo.find(seq_fita_mais)
        if pos >= 0:
            regioes[reg_nome] = Regiao(pos, pos + len(seq_fita_mais))
        else:
            violacoes.append(
                Violacao(
                    item=reg_nome,
                    regra="nao_encontrado",
                    valor=None,
                    faixa="",
                    mensagem=f"Região {reg_nome} não encontrada no alvo",
                )
            )

    return regioes, violacoes


def verificar(conjunto: ConjuntoLAMP, alvo: str | None, cfg: dict) -> list[Violacao]:
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

    regioes, v_loc = localizar(conjunto, alvo)
    violacoes.extend(v_loc)

    # Distâncias na ordem do arquivo de configuração
    distancias_cfg = cfg.get("distancias", {})
    for dist_nome, faixa in distancias_cfg.items():
        d_min, d_max = float(faixa[0]), float(faixa[1])
        val = None

        if dist_nome == "F2_B2" and "B2c" in regioes and "F2" in regioes:
            val = regioes["B2c"].fim - regioes["F2"].ini
        elif dist_nome == "F2_F1" and "F1" in regioes and "F2" in regioes:
            val = regioes["F1"].ini - regioes["F2"].ini
        elif dist_nome == "B1c_B2c" and "B2c" in regioes and "B1c" in regioes:
            val = regioes["B2c"].fim - regioes["B1c"].fim
        elif dist_nome == "F3_F2_gap" and "F2" in regioes and "F3" in regioes:
            val = regioes["F2"].ini - regioes["F3"].fim
        elif dist_nome == "B2c_B3c_gap" and "B3c" in regioes and "B2c" in regioes:
            val = regioes["B3c"].ini - regioes["B2c"].fim
        elif dist_nome == "F1_B1c_gap" and "B1c" in regioes and "F1" in regioes:
            val = regioes["B1c"].ini - regioes["F1"].fim

        if val is not None:
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
    if conjunto.LF is not None and "LF" in regioes and "F2" in regioes and "F1" in regioes:
        fim_f2 = regioes["F2"].fim
        ini_f1 = regioes["F1"].ini
        r_lf = regioes["LF"]
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

    if conjunto.LB is not None and "LB" in regioes and "B1c" in regioes and "B2c" in regioes:
        fim_b1c = regioes["B1c"].fim
        ini_b2c = regioes["B2c"].ini
        r_lb = regioes["LB"]
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
