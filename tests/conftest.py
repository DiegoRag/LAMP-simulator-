"""Configuração comum dos testes do lamp-sim.

Os testes do primers precisam de primer3-py e biopython; os da interface web precisam deles e de
fastapi, jinja2, python-multipart e httpx2. Num ambiente sem eles, como o Colab do simulador, que
instala só o necessário para o lampsim, esses arquivos ficam fora da coleta em vez de derrubarem a
suíte inteira: cada pacote continua rodando sozinho. Na CI e no computador as dependências de
requirements.txt estão instaladas e tudo roda.
"""

import importlib.util

_DEPENDENCIAS = {
    "test_primers_*.py": ("primer3", "Bio"),
    "test_interface_*.py": ("primer3", "Bio", "fastapi", "jinja2", "python_multipart", "httpx2"),
}

collect_ignore_glob = [
    padrao
    for padrao, modulos in _DEPENDENCIAS.items()
    if not all(importlib.util.find_spec(m) for m in modulos)
]


def pytest_report_header(config):
    if collect_ignore_glob:
        return "fora da coleta por falta de dependência: " + ", ".join(collect_ignore_glob)
    return None
