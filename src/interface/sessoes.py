"""Registro de sessões por navegador, identificadas por um cookie aleatório."""

import secrets
import threading
import time
from collections import OrderedDict
from collections.abc import Callable
from dataclasses import dataclass

from interface.sessao import Sessao

PADRAO_SESSOES = {"cookie": "lampsim_sessao", "validade_min": 120, "max_sessoes": 500,
                  "cookie_seguro": False}
BYTES_TOKEN = 32
SEGUNDOS_POR_MINUTO = 60


def _numero_positivo(valor: object) -> bool:
    return isinstance(valor, (int, float)) and not isinstance(valor, bool) and valor > 0


def config_sessoes(cfg_web: dict) -> dict:
    """Seção [sessoes] do web.toml com fallback em PADRAO_SESSOES, tipos normalizados.

    Devolve {"cookie": str, "validade_min": float, "max_sessoes": int, "cookie_seguro": bool}.
    Valor ausente, de tipo errado ou não positivo cai no padrão (bool não conta como número).
    """
    bruta = cfg_web.get("sessoes", {}) if isinstance(cfg_web, dict) else {}
    if not isinstance(bruta, dict):
        bruta = {}
    cookie = bruta.get("cookie")
    if not isinstance(cookie, str) or not cookie:
        cookie = PADRAO_SESSOES["cookie"]
    validade = bruta.get("validade_min")
    if not _numero_positivo(validade):
        validade = PADRAO_SESSOES["validade_min"]
    maximo = bruta.get("max_sessoes")
    if not isinstance(maximo, int) or isinstance(maximo, bool) or maximo <= 0:
        maximo = PADRAO_SESSOES["max_sessoes"]
    seguro = bruta.get("cookie_seguro")
    if not isinstance(seguro, bool):
        seguro = PADRAO_SESSOES["cookie_seguro"]
    return {"cookie": cookie, "validade_min": float(validade), "max_sessoes": int(maximo),
            "cookie_seguro": seguro}


@dataclass
class EntradaSessao:
    """Sessão de um navegador, com o instante do último uso e a mensagem pendente."""

    token: str
    sessao: Sessao
    ultimo_uso: float
    mensagem: dict[str, str] | None = None   # {"tipo": "erro"|"info", "texto": str}


class Sessoes:
    """Registro thread-safe de uma Sessao por token de cookie."""

    def __init__(self, fabrica: Callable[[], Sessao], validade_min: float, max_sessoes: int,
                 relogio: Callable[[], float] = time.monotonic) -> None:
        self._fabrica = fabrica
        self._validade_s = validade_min * SEGUNDOS_POR_MINUTO
        self._max_sessoes = max_sessoes
        self._relogio = relogio
        self._entradas: OrderedDict[str, EntradaSessao] = OrderedDict()
        self._trava = threading.Lock()

    def obter(self, token: str | None) -> tuple[EntradaSessao, bool]:
        """Entrada do token e se ela é nova.

        Antes de tudo descarta as expiradas (sem uso há mais de validade_min * 60 s pelo relógio).
        Token conhecido → atualiza ultimo_uso, devolve (entrada, False).
        Token None, vazio, desconhecido ou expirado → cria entrada com token novo
        secrets.token_urlsafe(32) (nunca reaproveita o token do cliente), Sessao de fabrica(),
        devolve (entrada, True). Se passar de max_sessoes, descarta a menos usada recentemente.
        """
        with self._trava:
            agora = self._relogio()
            expirados = [t for t, e in self._entradas.items()
                         if agora - e.ultimo_uso > self._validade_s]
            for t in expirados:
                del self._entradas[t]
            if token and token in self._entradas:
                entrada = self._entradas[token]
                entrada.ultimo_uso = agora
                self._entradas.move_to_end(token)
                return entrada, False
            novo = secrets.token_urlsafe(BYTES_TOKEN)
            while novo in self._entradas:
                novo = secrets.token_urlsafe(BYTES_TOKEN)
            entrada = EntradaSessao(novo, self._fabrica(), agora)
            self._entradas[novo] = entrada
            while len(self._entradas) > self._max_sessoes:
                self._entradas.popitem(last=False)
            return entrada, True

    def tomar_mensagem(self, token: str) -> dict[str, str] | None:
        """Devolve e apaga a mensagem pendente da sessão (None se não houver ou token sumiu)."""
        with self._trava:
            entrada = self._entradas.get(token)
            if entrada is None:
                return None
            mensagem, entrada.mensagem = entrada.mensagem, None
            return mensagem

    def definir_mensagem(self, token: str, mensagem: dict[str, str]) -> None:
        """Guarda a mensagem pendente da sessão do token (ignora token inexistente)."""
        with self._trava:
            entrada = self._entradas.get(token)
            if entrada is not None:
                entrada.mensagem = mensagem

    def __len__(self) -> int:
        with self._trava:
            return len(self._entradas)
