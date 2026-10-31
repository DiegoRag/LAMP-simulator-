"""Execução de análise de primers em processo isolado com limite de tempo.

Permite interromper análises que excedem o tempo máximo configurado (por exemplo,
em genes-alvo excessivamente longos) sem travar a interface nem consumir CPU indefinidamente.
"""

import atexit
import multiprocessing
import multiprocessing.connection
import threading
from pathlib import Path
from typing import Any

from interface.sessao import Execucao

TEMPO_MAX_PADRAO_S = 45.0
PARTIDA_MAX_S = 60.0  # tempo para o processo novo importar tudo e avisar que está pronto
ESPERA_JOIN_S = 5.0
ESPERA_ENCERRAR_S = 2.0
ESPERA_TERMINATE_S = 1.0
MAX_TENTATIVAS = 2

MENSAGEM_TEMPO = (
    "A análise passou de {limite} s e foi interrompida. O gene-alvo pode ser grande "
    "demais: envie só o gene (não o genoma inteiro) ou um FASTA com menos registros."
)


class TempoEsgotado(ValueError):
    """Tempo limite de análise esgotado no processo isolado."""


def tempo_max_s(cfg_web: dict) -> float:
    """[analise] tempo_max_s com fallback; número (não bool) > 0, senão o padrão."""
    valor = cfg_web.get("analise", {}).get("tempo_max_s", TEMPO_MAX_PADRAO_S)
    if isinstance(valor, bool):
        return TEMPO_MAX_PADRAO_S
    if isinstance(valor, (int, float)) and valor > 0:
        return float(valor)
    return TEMPO_MAX_PADRAO_S


def _trabalhador(conexao: Any, raiz: str, cfg_primers: dict, cfg_web: dict) -> None:
    """Laço do processo filho: importa, avisa ("pronto",), e para cada pedido
    (nome_conjunto, bytes_conjunto, nome_alvo, bytes_alvo) cria uma Sessao nova, carrega os dois,
    roda analisar() e devolve ("ok", execucao) ou ("erro", str(e)) para ValueError/RuntimeError.
    Pedido None encerra o laço."""
    import sys
    from pathlib import Path

    raiz_path = Path(raiz)
    src_dir = str(raiz_path / "src")
    if src_dir not in sys.path:
        sys.path.insert(0, src_dir)

    from interface.sessao import Sessao

    try:
        conexao.send(("pronto",))
    except (EOFError, BrokenPipeError, OSError):
        return

    while True:
        try:
            pedido = conexao.recv()
        except (EOFError, KeyboardInterrupt):
            break
        if pedido is None:
            break
        nome_conjunto, bytes_conjunto, nome_alvo, bytes_alvo = pedido
        try:
            sessao = Sessao(raiz_path, cfg_primers, cfg_web)
            sessao.carregar_conjunto(nome_conjunto, bytes_conjunto)
            sessao.carregar_alvo(nome_alvo, bytes_alvo)
            execucao = sessao.analisar()
            conexao.send(("ok", execucao))
        except (ValueError, RuntimeError) as e:
            try:
                conexao.send(("erro", str(e)))
            except (EOFError, BrokenPipeError, OSError):
                break
        except Exception as e:
            try:
                conexao.send(("erro", str(e)))
            except (EOFError, BrokenPipeError, OSError):
                break


class AnalisadorIsolado:
    """Gerencia processo isolado para análise de primers com limite de tempo."""

    def __init__(self, raiz: Path, cfg_primers: dict, cfg_web: dict, tempo_max: float) -> None:
        self.raiz = Path(raiz)
        self.cfg_primers = cfg_primers
        self.cfg_web = cfg_web
        self.tempo_max: float = float(tempo_max)
        self._ctx = multiprocessing.get_context("spawn")
        self._processo: Any = None
        self._conexao: Any = None
        self._lock = threading.Lock()
        self._encerrado = False
        atexit.register(self.encerrar)

    def _limpar_processo(self) -> None:
        proc = self._processo
        conn = self._conexao
        self._processo = None
        self._conexao = None
        if conn is not None:
            try:
                conn.close()
            except Exception:
                pass
        if proc is not None:
            try:
                if proc.is_alive():
                    proc.terminate()
                proc.join(timeout=ESPERA_ENCERRAR_S)
            except Exception:
                pass

    def _criar_processo(self) -> None:
        self._limpar_processo()
        pai, filho = self._ctx.Pipe()
        proc = self._ctx.Process(
            target=_trabalhador,
            args=(filho, str(self.raiz), self.cfg_primers, self.cfg_web),
            daemon=True,
        )
        proc.start()
        filho.close()
        if pai.poll(PARTIDA_MAX_S):
            try:
                msg = pai.recv()
            except (EOFError, OSError):
                msg = None
            if msg == ("pronto",):
                self._processo = proc
                self._conexao = pai
                return
        try:
            if proc.is_alive():
                proc.terminate()
            proc.join(timeout=ESPERA_ENCERRAR_S)
        except Exception:
            pass
        pai.close()
        raise RuntimeError("Falha ao inicializar o processo de análise.")

    def _iniciar(self) -> None:
        with self._lock:
            if self._encerrado:
                return
            if self._processo is not None and self._processo.is_alive() and self._conexao is not None:
                return
            try:
                self._criar_processo()
            except Exception:
                pass

    def analisar(
        self,
        nome_conjunto: str,
        bytes_conjunto: bytes,
        nome_alvo: str,
        bytes_alvo: bytes,
    ) -> Execucao:
        """Roda a análise no processo filho e devolve a Execucao (id a ser trocado por
        Sessao.registrar). ValueError com a mensagem do filho em erro de entrada;
        TempoEsgotado(MENSAGEM_TEMPO.format(limite=...)) se passar de tempo_max."""
        with self._lock:
            if self._encerrado:
                raise RuntimeError("Analisador já encerrado.")

            tentativa = 0
            while tentativa < MAX_TENTATIVAS:
                tentativa += 1
                if self._processo is None or not self._processo.is_alive() or self._conexao is None:
                    self._criar_processo()

                pedido = (nome_conjunto, bytes_conjunto, nome_alvo, bytes_alvo)
                try:
                    self._conexao.send(pedido)
                except (EOFError, BrokenPipeError, OSError):
                    self._limpar_processo()
                    if tentativa < MAX_TENTATIVAS:
                        continue
                    raise RuntimeError("Falha no processo de análise.") from None

                limite_tempo = self.tempo_max
                chegou = self._conexao.poll(limite_tempo)
                if not chegou:
                    proc = self._processo
                    conn = self._conexao
                    self._processo = None
                    self._conexao = None
                    if proc is not None:
                        try:
                            if proc.is_alive():
                                proc.terminate()
                            proc.join(timeout=ESPERA_JOIN_S)
                        except Exception:
                            pass
                    if conn is not None:
                        try:
                            conn.close()
                        except Exception:
                            pass
                    threading.Thread(target=self._iniciar, daemon=True).start()
                    limite_fmt = (
                        int(limite_tempo)
                        if isinstance(limite_tempo, (int, float)) and limite_tempo == int(limite_tempo)
                        else limite_tempo
                    )
                    raise TempoEsgotado(MENSAGEM_TEMPO.format(limite=limite_fmt))

                try:
                    status, resultado = self._conexao.recv()
                except (EOFError, BrokenPipeError, OSError):
                    self._limpar_processo()
                    if tentativa < MAX_TENTATIVAS:
                        continue
                    raise RuntimeError("Falha no processo de análise.") from None

                if status == "ok":
                    return resultado
                if status == "erro":
                    raise ValueError(resultado)
                raise RuntimeError(f"Resposta inesperada do processo trabalhador: {status}")

            raise RuntimeError("Falha no processo de análise.")

    def encerrar(self) -> None:
        """Encerra o processo trabalhador e libera recursos."""
        with self._lock:
            self._encerrado = True
            proc = self._processo
            conn = self._conexao
            self._processo = None
            self._conexao = None
            if conn is not None:
                try:
                    conn.send(None)
                except Exception:
                    pass
            if proc is not None:
                try:
                    proc.join(timeout=ESPERA_ENCERRAR_S)
                except Exception:
                    pass
                if proc.is_alive():
                    try:
                        proc.terminate()
                        proc.join(timeout=ESPERA_TERMINATE_S)
                    except Exception:
                        pass
            if conn is not None:
                try:
                    conn.close()
                except Exception:
                    pass
