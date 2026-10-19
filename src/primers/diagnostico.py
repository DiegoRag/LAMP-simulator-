"""Desenho vazio: que regra eliminou os candidatos, e a faixa relaxada para tentar de novo."""

import copy

from primers.candidatos import REGRAS, TIPOS_BASE

AVISO_RELAXADO = "faixas relaxadas de [desenho_relaxado]: GC e Tm fora do guia da Eiken"

NOMES_REGRA = {
    "ambigua": "base ambígua (IUPAC ou N) na janela",
    "gc_pct": "GC fora da faixa de [gc_pct]",
    "tm_c": "Tm de referência fora da faixa de [tm_c]",
    "extremidade": "ponta menos estável que [extremidade] dg_max",
    "repetida": "sequência repetida no alvo ou no complemento reverso",
    "grampo": "grampo estável na temperatura da reação",
    "homodimero": "homodímero estável na temperatura da reação",
    "teto": "teto de candidatos por janela",
}


def config_relaxada(cfg: dict) -> dict:
    """Cópia do cfg com [gc_pct] e [tm_c] trocados pelas faixas de [desenho_relaxado]."""
    rel = cfg["desenho_relaxado"]
    novo = copy.deepcopy(cfg)
    novo["gc_pct"] = {"min": float(rel["gc_pct"][0]), "max": float(rel["gc_pct"][1])}
    folga = float(rel["tm_folga_C"])
    novo["tm_c"] = {k: [float(v[0]) - folga, float(v[1]) + folga] for k, v in cfg["tm_c"].items()}
    return novo


def _juntar(itens: list[str]) -> str:
    return itens[0] if len(itens) == 1 else ", ".join(itens[:-1]) + " e " + itens[-1]


def diagnosticar(estatisticas: dict, relaxado: bool) -> str:
    """Explica um desenho vazio: a regra que mais eliminou candidatos, ou a geometria que não fechou.

    Lê estatisticas["candidatos_por_tipo"] e estatisticas["eliminados_por_tipo"] (de desenhar).
    """
    cands = estatisticas["candidatos_por_tipo"]
    elim = estatisticas["eliminados_por_tipo"]
    vazios = [t for t in TIPOS_BASE if cands[t] == 0]
    if vazios:
        grupos: dict[str, list[str]] = {}
        for t in vazios:
            regra = max(REGRAS, key=lambda r: (elim[t][r], -REGRAS.index(r)))
            grupos.setdefault(regra if elim[t][regra] > 0 else "", []).append(t)
        partes = []
        for regra, tipos in grupos.items():
            if regra:
                partes.append(
                    f"Nenhum candidato a {_juntar(tipos)}: a regra que mais eliminou janelas foi "
                    f"{NOMES_REGRA[regra]}."
                )
            else:
                partes.append(f"Nenhuma janela para {_juntar(tipos)}: o alvo é curto demais.")
        texto = " ".join(partes)
    else:
        escasso = min(TIPOS_BASE, key=lambda t: (cands[t], TIPOS_BASE.index(t)))
        texto = (
            f"Há candidatos para as seis regiões (o mais escasso é {escasso}, com {cands[escasso]}), "
            "mas nenhuma combinação respeita as distâncias de [distancias]."
        )
    if not relaxado:
        texto += " Tente o desenho com as faixas relaxadas de [desenho_relaxado]."
    return texto
