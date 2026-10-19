"""Módulo de exportação de resultados em CSV e JSON."""

import csv
import hashlib
import json
import subprocess
from dataclasses import asdict
from pathlib import Path

from primers.analise import Resultado
from primers.casamento import erros_de_base
from primers.desenho import AVISO, ResultadoDesenho
from primers.entrada import Leitura
from primers.incompat import regioes_para_json
from primers.modelo import ConjuntoLAMP
from primers.termo import Condicoes, CondicoesReferencia

COORDENADAS = (
    "regioes: base 0, intervalo semiaberto [ini, fim), na fita indicada em 'fita' "
    "(no GenBank, de ini + 1 até fim); posicao e base_no_primer dos erros de base: base 1"
)
COORDENADAS_DESENHO = (
    "regioes: base 0, intervalo semiaberto [ini, fim), na fita +; no GenBank, de ini + 1 até fim"
)
CSV_COORDENADAS = (
    "# colunas _ini e _fim: base 0, intervalo semiaberto [ini, fim), na fita +; "
    "no GenBank, de ini + 1 até fim"
)

COLUNAS_PRIMERS = [
    "nome",
    "sequencia",
    "comprimento",
    "gc_pct",
    "tm_c",
    "tm_reacao_c",
    "dg_extremidade_kcal",
    "extremidade",
    "grampo_tm_c",
    "grampo_dg_kcal",
    "homodimero_dg_kcal",
    "avisos",
]

COLUNAS_VIOLACOES = [
    "item",
    "regra",
    "valor",
    "faixa",
    "mensagem",
]


def _formatar_campo_primers(v: object) -> str:
    if v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.2f}"
    return str(v)


def _formatar_valor_violacao(valor: float | None, regra: str) -> str:
    if valor is None:
        return ""
    if valor.is_integer() and regra in ("comprimento", "distancia"):
        return f"{int(valor)}"
    return f"{valor:.2f}"


def versao_do_codigo() -> str:
    """Commit do Git do repositório (git rev-parse HEAD), ou "desconhecida"."""
    try:
        repo_dir = Path(__file__).resolve().parents[2]
        proc = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=repo_dir,
            capture_output=True,
            text=True,
            check=False,
        )
        if proc.returncode == 0 and proc.stdout.strip():
            return proc.stdout.strip()
    except Exception:
        pass
    return "desconhecida"


def descrever_alvo(
    sequencia: str,
    *,
    id_registro: str | None = None,
    arquivo: str | None = None,
    circular: bool = False,
) -> dict:
    """Identificação do alvo para o resultado.json, sem a sequência (só o tamanho e o hash)."""
    return {
        "id": id_registro,
        "arquivo": arquivo,
        "comprimento": len(sequencia),
        "sha256": hashlib.sha256(sequencia.encode("utf-8")).hexdigest(),
        "circular": bool(circular),
    }


def _corte(nome: str, parte5: str, ligante: str, parte3: str, nomes: tuple[str, str], partes: bool) -> dict:
    """Como o primer interno foi cortado em partes, em posições base 1 do primer inteiro."""
    n5, nl = len(parte5), len(ligante)
    return {
        "sequencia": parte5 + ligante + parte3,
        "origem": "partes informadas" if partes else f"{nome} inteiro separado nas partes",
        nomes[0]: [1, n5],
        "ligante": [n5 + 1, n5 + nl] if nl else None,
        nomes[1]: [n5 + nl + 1, n5 + nl + len(parte3)],
    }


def descrever_conjunto(leitura: Leitura | None) -> dict | None:
    """Primers do conjunto e o corte do FIP e do BIP usados na análise (None sem conjunto completo)."""
    if leitura is None or leitura.conjunto is None:
        return None
    c: ConjuntoLAMP = leitura.conjunto
    return {
        "F3": c.F3,
        "B3": c.B3,
        "F2": c.F2,
        "F1c": c.F1c,
        "B2": c.B2,
        "B1c": c.B1c,
        "LF": c.LF,
        "LB": c.LB,
        "ligante": c.ligante,
        "FIP": _corte("FIP", c.F1c, c.ligante, c.F2, ("F1c", "F2"), "FIP" not in leitura.separados),
        "BIP": _corte("BIP", c.B1c, c.ligante, c.B2, ("B1c", "B2"), "BIP" not in leitura.separados),
    }


def exportar(
    resultado: Resultado,
    pasta: str | Path,
    cfg: dict,
    *,
    leitura: Leitura | None = None,
    alvo: dict | None = None,
    eficiencia: dict | None = None,
) -> None:
    """Exporta resultados da análise para pasta informada gerando 4 arquivos determinísticos.

    Arquivos gerados:
        - primers.csv: métricas físico-químicas e termodinâmicas por primer.
        - heterodimeros.csv: matriz quadrada de ΔG entre primers de tubo.
        - violacoes.csv: lista de violações encontradas.
        - resultado.json: estrutura consolidada completa. Para reproduzir a análise só pelo
          arquivo, ele traz também o alvo (descrever_alvo: id, tamanho e hash, nunca a
          sequência), a fita, os primers com o corte do FIP e do BIP (descrever_conjunto),
          as posições aproximadas com os erros de base, as incompatibilidades e, quando quem
          chama a passa (a esteira), a eficiência prevista.
    """
    pasta_p = Path(pasta)
    pasta_p.mkdir(parents=True, exist_ok=True)

    # 1. primers.csv
    caminho_primers = pasta_p / "primers.csv"
    with open(caminho_primers, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(COLUNAS_PRIMERS)
        for linha in resultado.linhas:
            row = [_formatar_campo_primers(linha[k]) for k in COLUNAS_PRIMERS]
            writer.writerow(row)

    # 2. heterodimeros.csv
    caminho_het = pasta_p / "heterodimeros.csv"
    nomes_tubo = resultado.heterodimeros["nomes"]
    matriz_het = resultado.heterodimeros["dg_kcal"]
    with open(caminho_het, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(["nome"] + nomes_tubo)
        for nome, linha_dg in zip(nomes_tubo, matriz_het, strict=True):
            row = [nome] + [f"{x:.2f}" for x in linha_dg]
            writer.writerow(row)

    # 3. violacoes.csv
    caminho_viol = pasta_p / "violacoes.csv"
    with open(caminho_viol, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(COLUNAS_VIOLACOES)
        for v in resultado.violacoes:
            val_str = _formatar_valor_violacao(v.valor, v.regra)
            writer.writerow([v.item, v.regra, val_str, v.faixa, v.mensagem])

    # 4. resultado.json
    versao_codigo = versao_do_codigo()

    cond = Condicoes.de_config(cfg)
    ref = CondicoesReferencia.de_config(cfg)
    condicoes_dict = {
        "reacao": {
            "na_mM": cond.na_mM,
            "k_mM": cond.k_mM,
            "nh4_mM": cond.nh4_mM,
            "monovalentes_mM": cond.monovalentes_mM,
            "mg_mM": cond.mg_mM,
            "dntp_cada_mM": cond.dntp_cada_mM,
            "dntp_total_mM": cond.dntp_total_mM,
            "mg_livre_mM": cond.mg_livre_mM,
            "temp_C": cond.temp_C,
        },
        "referencia": {
            "convencao": "PrimerExplorer V5",
            "oligo_nM": ref.oligo_nM,
            "na_mM": ref.na_mM,
            "mg_mM": ref.mg_mM,
        },
    }

    regioes_json = None
    if resultado.regioes is not None:
        regioes_json = {k: [r.ini, r.fim] for k, r in resultado.regioes.items()}

    aproximadas_json = None
    if resultado.regioes_aproximadas is not None:
        aproximadas_json = {
            k: {"regiao": [p.regiao.ini, p.regiao.fim], "erros": [asdict(e) for e in erros_de_base(k, p)]}
            for k, p in resultado.regioes_aproximadas.items()
        }

    incompat_json = None
    if resultado.incompat is not None:
        incompat_json = {
            "total": resultado.incompat.total,
            "ponta": resultado.incompat.ponta3,
            "bases_ponta": int(cfg.get("casamento", {}).get("bases_ponta_3p", 5)),
            "nao_encontradas": list(resultado.incompat.nao_encontradas),
            "por_regiao": regioes_para_json(resultado.incompat),
        }

    violacoes_json = [
        {
            "item": v.item,
            "regra": v.regra,
            "valor": v.valor,
            "faixa": v.faixa,
            "mensagem": v.mensagem,
        }
        for v in resultado.violacoes
    ]

    dados = {
        "versao_codigo": versao_codigo,
        "coordenadas": COORDENADAS,
        "alvo": alvo,
        "fita": resultado.fita,
        "conjunto": descrever_conjunto(leitura),
        "condicoes": condicoes_dict,
        "config": cfg,
        "avisos": resultado.avisos,
        "primers": resultado.linhas,
        "heterodimeros": resultado.heterodimeros,
        "violacoes": violacoes_json,
        "regioes": regioes_json,
        "regioes_aproximadas": aproximadas_json,
        "incompat": incompat_json,
        "eficiencia_prevista": eficiencia,
    }

    caminho_json = pasta_p / "resultado.json"
    with open(caminho_json, "w", encoding="utf-8") as f:
        json.dump(dados, f, ensure_ascii=False, indent=2)


COLUNAS_CONJUNTOS_DESENHO = [
    "posicao",
    "pontuacao",
    "desvio_tm",
    "pior_heterodimero_kcal",
    "F3",
    "B3",
    "FIP",
    "BIP",
    "LF",
    "LB",
    "F1c",
    "F2",
    "B1c",
    "B2",
    "F3_ini",
    "F3_fim",
    "F2_ini",
    "F2_fim",
    "F1_ini",
    "F1_fim",
    "B1c_ini",
    "B1c_fim",
    "B2c_ini",
    "B2c_fim",
    "B3c_ini",
    "B3c_fim",
    "LF_ini",
    "LF_fim",
    "LB_ini",
    "LB_fim",
    "avisos",
]


def exportar_desenho(resultado: ResultadoDesenho, alvo: str, pasta: str | Path, cfg: dict) -> None:
    """Exporta resultados do desenho de primers LAMP gerando conjuntos.csv e resultado.json."""
    pasta_p = Path(pasta)
    pasta_p.mkdir(parents=True, exist_ok=True)

    # 1. conjuntos.csv
    caminho_csv = pasta_p / "conjuntos.csv"
    with open(caminho_csv, "w", encoding="utf-8", newline="") as f:
        f.write(f"# {AVISO}\n")
        f.write(f"{CSV_COORDENADAS}\n")
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(COLUNAS_CONJUNTOS_DESENHO)
        for c in resultado.conjuntos:
            row = [
                str(c.posicao),
                f"{c.pontuacao:.2f}",
                f"{c.desvio_tm:.2f}",
                f"{c.pior_heterodimero_kcal:.2f}",
                c.conjunto.F3,
                c.conjunto.B3,
                c.conjunto.FIP,
                c.conjunto.BIP,
                c.conjunto.LF or "",
                c.conjunto.LB or "",
                c.conjunto.F1c,
                c.conjunto.F2,
                c.conjunto.B1c,
                c.conjunto.B2,
            ]
            for r_nome in ("F3", "F2", "F1", "B1c", "B2c", "B3c", "LF", "LB"):
                if r_nome in c.regioes and c.regioes[r_nome] is not None:
                    row.append(str(c.regioes[r_nome].ini))
                    row.append(str(c.regioes[r_nome].fim))
                else:
                    row.append("")
                    row.append("")
            row.append("; ".join(c.avisos))
            writer.writerow(row)

    # 2. resultado.json
    versao_codigo = versao_do_codigo()

    conjuntos_json = []
    for c in resultado.conjuntos:
        primers_dict = {
            "F3": c.conjunto.F3,
            "B3": c.conjunto.B3,
            "F2": c.conjunto.F2,
            "F1c": c.conjunto.F1c,
            "B2": c.conjunto.B2,
            "B1c": c.conjunto.B1c,
            "FIP": c.conjunto.FIP,
            "BIP": c.conjunto.BIP,
        }
        if c.conjunto.LF is not None:
            primers_dict["LF"] = c.conjunto.LF
        if c.conjunto.LB is not None:
            primers_dict["LB"] = c.conjunto.LB

        regioes_dict = {}
        for r_nome in ("F3", "F2", "F1", "B1c", "B2c", "B3c", "LF", "LB"):
            if r_nome in c.regioes and c.regioes[r_nome] is not None:
                regioes_dict[r_nome] = [c.regioes[r_nome].ini, c.regioes[r_nome].fim]

        conjuntos_json.append(
            {
                "posicao": c.posicao,
                "pontuacao": c.pontuacao,
                "desvio_tm": c.desvio_tm,
                "pior_heterodimero_kcal": c.pior_heterodimero_kcal,
                "primers": primers_dict,
                "regioes": regioes_dict,
                "avisos": c.avisos,
            }
        )

    dados = {
        "versao_codigo": versao_codigo,
        "aviso": AVISO,
        "coordenadas": COORDENADAS_DESENHO,
        "config": cfg,
        "alvo": {
            "comprimento": len(alvo),
            "sha256": hashlib.sha256(alvo.encode("utf-8")).hexdigest(),
        },
        "estatisticas": resultado.estatisticas,
        "conjuntos": conjuntos_json,
    }

    caminho_json = pasta_p / "resultado.json"
    with open(caminho_json, "w", encoding="utf-8") as f:
        json.dump(dados, f, ensure_ascii=False, indent=2)
