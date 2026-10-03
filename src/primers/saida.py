"""Módulo de exportação de resultados em CSV e JSON."""

import csv
import json
from pathlib import Path
import subprocess
from primers.analise import Resultado
from primers.termo import Condicoes

COLUNAS_PRIMERS = [
    "nome",
    "sequencia",
    "comprimento",
    "gc_pct",
    "tm_c",
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


def exportar(resultado: Resultado, pasta: str | Path, cfg: dict) -> None:
    """Exporta resultados da análise para pasta informada gerando 4 arquivos determinísticos.

    Arquivos gerados:
        - primers.csv: métricas físico-químicas e termodinâmicas por primer.
        - heterodimeros.csv: matriz quadrada de ΔG entre primers de tubo.
        - violacoes.csv: lista de violações encontradas.
        - resultado.json: estrutura consolidada completa.
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
        for nome, linha_dg in zip(nomes_tubo, matriz_het):
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
            versao_codigo = proc.stdout.strip()
        else:
            versao_codigo = "desconhecida"
    except Exception:
        versao_codigo = "desconhecida"

    cond = Condicoes.de_config(cfg)
    condicoes_dict = {
        "na_mM": cond.na_mM,
        "mg_mM": cond.mg_mM,
        "dntp_mM": cond.dntp_mM,
        "temp_C": cond.temp_C,
    }

    regioes_json = None
    if resultado.regioes is not None:
        regioes_json = {k: [r.ini, r.fim] for k, r in resultado.regioes.items()}

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
        "condicoes": condicoes_dict,
        "config": cfg,
        "avisos": resultado.avisos,
        "primers": resultado.linhas,
        "heterodimeros": resultado.heterodimeros,
        "violacoes": violacoes_json,
        "regioes": regioes_json,
    }

    caminho_json = pasta_p / "resultado.json"
    with open(caminho_json, "w", encoding="utf-8") as f:
        json.dump(dados, f, ensure_ascii=False, indent=2)
