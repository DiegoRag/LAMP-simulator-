#!/bin/bash
# Roda um executor por spec, em paralelo, cada um na sua pasta lamp-sim-wK.
# Uso: executar-fase.sh N [K ...]   (sem K, roda todas as specs sK.md da fase)
set -uo pipefail
LAMP="$HOME/Documents/Lamp"
ORQ="$LAMP/orquestracao"
N="$1"; shift
FASE="$ORQ/fase-$N"
PRINCIPAL="${PRINCIPAL:-claude-sonnet-5-5-medium}"
RESERVA="${RESERVA:-gemini-3.8-flash-medium}"

python3 -c "import json,os;json.load(open(os.path.expanduser('~/.gemini/antigravity-cli/settings.json')))" \
  || { echo "settings.json do agy inválido: parei"; exit 1; }

if [ "$#" -gt 0 ]; then KS="$*"; else
  KS=$(ls "$FASE"/s*.md | sed -E 's#.*/s([0-9]+)\.md#\1#' | sort -n | tr '\n' ' ')
fi

for K in $KS; do
  W="$LAMP/lamp-sim-w$K"
  mkdir -p "$W/.delegar"
  cp "$FASE/contrato.md" "$W/.delegar/contrato.md"
  cp "$FASE/s$K.md" "$W/.delegar/s$K.md"
  PROMPT="Você é o executor da spec s$K. Leia .delegar/contrato.md e .delegar/s$K.md e implemente \
exatamente o que a spec pede, seguindo o contrato. As outras sessões implementam o resto do contrato \
ao mesmo tempo: importe direto o que o contrato promete, sem try/except ImportError nem cópias de \
reserva. Você NÃO pode executar comandos (serão negados): use só as ferramentas de leitura e escrita \
de arquivos. Não pare para esperar nada. Apenas escreva/edite os arquivos listados na spec. Não \
execute comandos, testes nem instalações. Não altere arquivos fora da lista. Ao terminar, liste os \
arquivos que alterou."
  ( cd "$W" && python3 "$ORQ/agy_retry.py" "$W/.delegar/saida.json" --reserva "$RESERVA" -- \
      agy -p "$PROMPT" --model "$PRINCIPAL" --output-format json --print-timeout 30m \
      > "$FASE/log-s$K.txt" 2>&1; echo "s$K terminou (código $?)" ) &
  echo "s$K iniciado em $W"
done
wait
echo "--- resumo ---"
cat "$FASE"/log-s*.txt
