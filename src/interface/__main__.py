"""Ponto de entrada para execução da interface web local."""

import tomllib
from pathlib import Path

import uvicorn

from interface.app import criar_app


def main() -> None:
    """Inicia o servidor web local da interface de primers."""
    raiz = Path(__file__).resolve().parents[2]
    cfg_web = tomllib.loads((raiz / "config" / "web.toml").read_text(encoding="utf-8"))
    host = cfg_web["servidor"]["host"]
    porta = cfg_web["servidor"]["porta"]
    print(f"LAMP Primers em http://{host}:{porta}  (Ctrl+C para parar)")
    uvicorn.run(criar_app(raiz), host=host, port=int(porta))


if __name__ == "__main__":
    main()
