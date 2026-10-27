#!/bin/bash
# Copia o repositório (sem .git, venv, dados e caches) para cada pasta de destino.
# Uso: sincronizar.sh DEST [DEST...]   (caminhos absolutos ou relativos a ~/Documents/Lamp)
set -euo pipefail
LAMP="$HOME/Documents/Lamp"
ORIGEM="$LAMP/lamp-sim/"
for d in "$@"; do
  case "$d" in /*) dest="$d" ;; *) dest="$LAMP/$d" ;; esac
  mkdir -p "$dest"
  rsync -a --delete \
    --exclude .git --exclude .venv --exclude dados/ --exclude saida/ --exclude _dados_local/ \
    --exclude .pytest_cache --exclude .ruff_cache --exclude __pycache__ --exclude graphify-out/ \
    --exclude .delegar/ "$ORIGEM" "$dest/"
  echo "sincronizado: $dest"
done
