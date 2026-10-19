"""Lotes de amostras reprodutíveis e o manifesto que os descreve.

Um lote é definido pela semente-mestre, pelo número de amostras por grupo e pelos
parâmetros do simulador (química, cor, cinética, leitura e classes). O manifesto
guarda tudo isso, os parâmetros sorteados de cada amostra e um hash que identifica o
lote: mesma entrada, mesmo hash; parâmetros diferentes, hash diferente.
Os dados gerados (imagens, tabelas) não vão para o git; o manifesto vai.
"""

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from . import __version__
from . import config as _config
from .classes import sortear_amostra


def parametros_do_simulador():
    """Todos os parâmetros que definem o resultado, lidos do config no momento da chamada."""
    return {
        "versao_parametros": _config.VERSAO_PARAMETROS,
        "tempo": _config.TEMPO,
        "quimica": _config.QUIMICA,
        "espectro": _config.ESPECTRO,
        "cinetica": _config.CINETICA,
        "mix": _config.MIX,
        "leitura": _config.LEITURA,
        "classes": _config.CLASSES,
        "lote": _config.LOTE,
    }


def gerar_lote(n_por_classe, semente_mestre, grupos=None):
    """Gera n amostras por grupo de config.LOTE (grupo, classe óptica, papel). Cada amostra
    recebe sua própria semente, derivada da semente-mestre (numpy SeedSequence), então o
    lote é reprodutível."""
    grupos = _config.LOTE if grupos is None else grupos
    filhas = np.random.SeedSequence(semente_mestre).spawn(n_por_classe * len(grupos))
    lote, i = [], 0
    for grupo, classe, papel in grupos:
        for j in range(n_por_classe):
            semente = int(filhas[i].generate_state(1)[0])
            lote.append(sortear_amostra(classe, semente, id_amostra=f"{grupo}-{j:04d}", papel=papel))
            i += 1
    return lote


def hash_lote(lote, parametros=None):
    """Hash SHA-256 (16 primeiros dígitos) dos parâmetros do simulador e de todas as amostras."""
    parametros = parametros_do_simulador() if parametros is None else parametros
    conteudo = json.dumps(
        {"parametros": parametros, "amostras": [a.como_dict() for a in lote]}, sort_keys=True
    )
    return hashlib.sha256(conteudo.encode("utf-8")).hexdigest()[:16]


def montar_manifesto(lote, semente_mestre, n_por_classe):
    parametros = parametros_do_simulador()
    return {
        "lampsim_versao": __version__,
        "versao_parametros": _config.VERSAO_PARAMETROS,
        "semente_mestre": semente_mestre,
        "n_por_classe": n_por_classe,
        "n_total": len(lote),
        "hash_lote": hash_lote(lote, parametros),
        "gerado_em_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "parametros": parametros,
        "amostras": [a.como_dict() for a in lote],
    }


def salvar_manifesto(manifesto, pasta):
    """Grava o manifesto como JSON; o nome do arquivo leva o hash do lote."""
    pasta = Path(pasta)
    pasta.mkdir(parents=True, exist_ok=True)
    caminho = pasta / f"manifesto_{manifesto['hash_lote']}.json"
    caminho.write_text(json.dumps(manifesto, ensure_ascii=False, indent=2), encoding="utf-8")
    return caminho
