"""Roda uma chamada do agy e troca para o modelo reserva quando a cota acaba.

Uso:
    python3 agy_retry.py SAIDA.json --reserva MODELO -- agy -p "..." --model PRINCIPAL ...

Sai com 0 em sucesso real (status SUCCESS, resposta não vazia, nenhuma ação negada),
2 em falha da tarefa e 3 quando a cota acabou e não havia reserva.
"""

import json
import subprocess
import sys
from pathlib import Path

MAX_TENTATIVAS = 3


def rodar(cmd: list[str], saida: Path) -> dict:
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    saida.write_text(proc.stdout + proc.stderr, encoding="utf-8")
    primeira = proc.stdout.strip().splitlines()[0] if proc.stdout.strip() else ""
    try:
        dados = json.loads(primeira)
    except json.JSONDecodeError:
        dados = {"status": "ERRO"}
    dados["bruto"] = (proc.stdout + proc.stderr)[-2000:]
    return dados


def cota_esgotada(dados: dict) -> bool:
    texto = json.dumps(dados, ensure_ascii=False)
    return "RESOURCE_EXHAUSTED" in texto or "quota reached" in texto


def trocar_modelo(cmd: list[str], modelo: str) -> list[str]:
    novo = list(cmd)
    if "--model" in novo:
        novo[novo.index("--model") + 1] = modelo
    else:
        novo += ["--model", modelo]
    return novo


def main() -> int:
    args = sys.argv[1:]
    if "--" not in args:
        print(__doc__)
        return 2
    corte = args.index("--")
    opcoes, cmd = args[:corte], args[corte + 1 :]
    saida = Path(opcoes[0])
    reserva = opcoes[opcoes.index("--reserva") + 1] if "--reserva" in opcoes else None

    dados: dict = {}
    for _tentativa in range(MAX_TENTATIVAS):
        dados = rodar(cmd, saida)
        if cota_esgotada(dados):
            if reserva is None:
                print(f"{saida}: cota esgotada, sem reserva")
                return 3
            print(f"{saida}: cota esgotada, recomeçando com {reserva}")
            cmd = trocar_modelo(cmd, reserva)
            reserva = None
            continue
        break

    negadas = dados.get("denied_actions") or []
    resposta = str(dados.get("response") or "").strip()
    ok = dados.get("status") == "SUCCESS" and resposta and not negadas
    uso = (dados.get("usage") or {}).get("total_tokens")
    print(f"{saida}: status={dados.get('status')} negadas={len(negadas)} tokens={uso} "
          f"conversa={dados.get('conversation_id')}")
    return 0 if ok else 2


if __name__ == "__main__":
    sys.exit(main())
