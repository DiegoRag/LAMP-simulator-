"""Junta numa pasta de integração os arquivos que cada executor mudou em relação à base.

Uso: python3 integrar.py N BASE [K ...]
Cria orquestracao/int-N como cópia de BASE e copia por cima os arquivos alterados ou criados em
cada lamp-sim-wK. Lista quem mexeu em quê e acusa arquivos fora da lista da spec e conflitos.
"""

import filecmp
import re
import shutil
import sys
from pathlib import Path

LAMP = Path.home() / "Documents" / "Lamp"
ORQ = LAMP / "orquestracao"
IGNORAR = {".delegar", "__pycache__", ".pytest_cache", ".ruff_cache"}


def arquivos(raiz: Path) -> set[str]:
    return {
        str(p.relative_to(raiz))
        for p in raiz.rglob("*")
        if p.is_file() and not (IGNORAR & set(p.relative_to(raiz).parts))
    }


def permitidos(spec: Path) -> set[str]:
    """Caminhos da seção 'PODE editar' da spec (entre crases)."""
    texto = spec.read_text(encoding="utf-8")
    m = re.search(r"##\s*PODE editar(.*?)(\n##\s)", texto, re.S)
    return set(re.findall(r"`([^`]+\.[a-z]+)`", m.group(1))) if m else set()


def main() -> int:
    n, base = sys.argv[1], Path(sys.argv[2])
    fase = ORQ / f"fase-{n}"
    ks = sys.argv[3:] or sorted(
        (re.sub(r"\D", "", p.stem) for p in fase.glob("s*.md")), key=int
    )
    destino = ORQ / f"int-{n}"
    if destino.exists():
        shutil.rmtree(destino)
    shutil.copytree(base, destino)

    dono: dict[str, str] = {}
    problemas = []
    base_arqs = arquivos(base)
    for k in ks:
        w = LAMP / f"lamp-sim-w{k}"
        mudados = sorted(
            a for a in arquivos(w)
            if a not in base_arqs or not filecmp.cmp(base / a, w / a, shallow=False)
        )
        apagados = sorted(a for a in base_arqs if not (w / a).exists())
        lista = permitidos(fase / f"s{k}.md")
        print(f"s{k}: {len(mudados)} alterado(s) {mudados}")
        for a in mudados:
            if lista and a not in lista:
                problemas.append(f"s{k} mexeu fora da lista: {a}")
            if a in dono:
                problemas.append(f"conflito em {a}: {dono[a]} e s{k}")
            dono[a] = f"s{k}"
            (destino / a).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(w / a, destino / a)
        for a in apagados:
            problemas.append(f"s{k} apagou {a}")
    print("--- problemas ---" if problemas else "--- sem problemas ---")
    for p in problemas:
        print(p)
    return 1 if problemas else 0


if __name__ == "__main__":
    sys.exit(main())
