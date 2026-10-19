"""Fronteiras entre os pacotes (roadmap da esteira): primers e lampsim nunca se importam.

A esteira pode importar os dois; nenhum dos dois importa a esteira. Assim cada pacote roda
sozinho, e a única ligação entre eles é o contrato selecao.json (esteira.contrato). A interface
web (`interface`) pode importar todos, e nenhum pacote importa a interface.
"""

import ast
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[1] / "src"

PROIBIDOS = {
    "primers": {"lampsim", "esteira", "interface"},
    "lampsim": {"primers", "esteira", "interface"},
    "esteira": {"interface"},
    "interface": set(),
}


def _importados(arquivo: Path) -> set[str]:
    """Pacotes de primeiro nível importados de forma absoluta pelo arquivo."""
    arvore = ast.parse(arquivo.read_text(encoding="utf-8"))
    nomes: set[str] = set()
    for no in ast.walk(arvore):
        if isinstance(no, ast.Import):
            nomes.update(alias.name.split(".")[0] for alias in no.names)
        elif isinstance(no, ast.ImportFrom) and no.level == 0 and no.module:
            nomes.add(no.module.split(".")[0])
    return nomes


@pytest.mark.parametrize("pacote", sorted(PROIBIDOS))
def test_pacote_nao_importa_o_vizinho(pacote):
    arquivos = sorted((SRC / pacote).rglob("*.py"))
    assert arquivos, f"src/{pacote} sem nenhum módulo"
    achados = [
        f"{arquivo.relative_to(SRC).as_posix()} importa {nome}"
        for arquivo in arquivos
        for nome in sorted(_importados(arquivo) & PROIBIDOS[pacote])
    ]
    assert achados == []
