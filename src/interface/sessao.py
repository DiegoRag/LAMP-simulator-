"""Gerenciamento de sessão em memória para a interface web de primers."""

import tempfile
import threading
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from interface.formato import milhar
from primers.analise import Resultado, analisar
from primers.desenho import ConjuntoDesenhado, ResultadoDesenho, desenhar
from primers.entrada import (
    ORDEM_CANONICA,
    Leitura,
    RegistroAlvo,
    escolher_registro_alvo,
    ler_alvos,
    ler_conjunto,
)
from primers.restricoes import avaliar_primer_quimera


@dataclass(frozen=True)
class ArquivoCarregado:
    """Metadados e conteúdo bruto do conjunto de primers carregado."""

    nome: str
    conteudo: bytes
    avisos: tuple[str, ...]


@dataclass(frozen=True)
class AlvoCarregado:
    """Metadados e sequência normalizada do alvo genômico carregado."""

    nome: str
    descricao: str
    sequencia: str
    registros: tuple[RegistroAlvo, ...] = ()


@dataclass
class Execucao:
    """Registro histórico de uma análise ou de um desenho de primers."""

    id: int
    tipo: str  # "analise" ou "desenho"
    criada_em: datetime
    rotulo_conjunto: str
    rotulo_alvo: str
    alvo: AlvoCarregado | None
    leitura: Leitura | None  # apenas na análise
    resultado: Resultado | ResultadoDesenho
    origem: str = "usuario"  # do conjunto analisado: "usuario" (carregado) ou "desenho" (proposta)


def _salvar_temporario(dir_temp: str, nome_arquivo: str, conteudo: bytes) -> Path:
    ext = Path(nome_arquivo).suffix.lower()
    caminho = Path(dir_temp) / f"entrada{ext}"
    caminho.write_bytes(conteudo)
    return caminho


def _tratar_erro_leitor(erro: Exception, caminho_temp: Path, nome_arquivo: str) -> ValueError:
    msg = str(erro)
    msg = msg.replace(str(caminho_temp), nome_arquivo)
    msg = msg.replace(caminho_temp.as_posix(), nome_arquivo)
    msg = msg.replace(caminho_temp.name, nome_arquivo)
    return ValueError(f"{nome_arquivo}: {msg}")


class Sessao:
    """Estado em memória da interface web por sessão/servidor."""

    def __init__(self, raiz: Path, cfg_primers: dict, cfg_web: dict) -> None:
        """Inicializa a sessão com configurações e listas de execuções vazias."""
        self.raiz = Path(raiz)
        self.cfg_primers = cfg_primers
        self.cfg_web = cfg_web
        self.conjunto: ArquivoCarregado | None = None
        self.alvo: AlvoCarregado | None = None
        self.execucoes: list[Execucao] = []
        self._proximo_id: int = 1
        self._lock = threading.RLock()

    def _validar_tamanho(self, nome_arquivo: str, conteudo: bytes) -> None:
        tamanho_max_kb = self.cfg_web.get("entrada", {}).get("tamanho_max_kb", 8192)
        if len(conteudo) > tamanho_max_kb * 1024:
            raise ValueError(f"{nome_arquivo}: arquivo maior que {tamanho_max_kb} KB.")

    def carregar_conjunto(self, nome_arquivo: str, conteudo: bytes) -> None:
        """Valida e carrega o conjunto de primers na sessão."""
        with self._lock:
            self._validar_tamanho(nome_arquivo, conteudo)
            with tempfile.TemporaryDirectory() as tmp_dir:
                caminho_temp = _salvar_temporario(tmp_dir, nome_arquivo, conteudo)
                ligante = self.cfg_primers.get("desenho", {}).get("ligante", "")
                max_erros = int(self.cfg_primers.get("casamento", {}).get("max_erros", 3))
                circular = bool(self.cfg_primers.get("alvo", {}).get("circular", False))
                try:
                    leitura = ler_conjunto(
                        caminho_temp,
                        alvo=None,
                        ligante=ligante,
                        max_erros=max_erros,
                        circular=circular,
                    )
                except Exception as e:
                    raise _tratar_erro_leitor(e, caminho_temp, nome_arquivo) from e

            self.conjunto = ArquivoCarregado(
                nome=nome_arquivo,
                conteudo=conteudo,
                avisos=tuple(leitura.avisos),
            )

    def carregar_alvo(self, nome_arquivo: str, conteudo: bytes) -> None:
        """Valida e carrega o alvo FASTA na sessão."""
        with self._lock:
            self._validar_tamanho(nome_arquivo, conteudo)
            with tempfile.TemporaryDirectory() as tmp_dir:
                caminho_temp = _salvar_temporario(tmp_dir, nome_arquivo, conteudo)
                try:
                    registros = ler_alvos(caminho_temp)
                except Exception as e:
                    raise _tratar_erro_leitor(e, caminho_temp, nome_arquivo) from e

            escolhido = registros[0]
            if len(registros) > 1 and self.conjunto is not None:
                try:
                    with tempfile.TemporaryDirectory() as tmp_c_dir:
                        c_temp = _salvar_temporario(tmp_c_dir, self.conjunto.nome, self.conjunto.conteudo)
                        ligante = self.cfg_primers.get("desenho", {}).get("ligante", "")
                        max_erros = int(self.cfg_primers.get("casamento", {}).get("max_erros", 3))
                        circular = bool(self.cfg_primers.get("alvo", {}).get("circular", False))
                        leitura_tmp = ler_conjunto(
                            c_temp,
                            alvo=None,
                            ligante=ligante,
                            max_erros=max_erros,
                            circular=circular,
                        )
                        escolhido = escolher_registro_alvo(
                            registros,
                            leitura_tmp.conjunto or leitura_tmp.primers,
                        )
                except Exception:
                    escolhido = registros[0]

            self.alvo = AlvoCarregado(
                nome=nome_arquivo,
                descricao=escolhido.descricao,
                sequencia=escolhido.sequencia,
                registros=tuple(registros),
            )

    def exemplos_disponiveis(self) -> list[str]:
        """Lista os nomes dos exemplos configurados cujos arquivos existem em disco."""
        with self._lock:
            disponiveis: list[str] = []
            for ex in self.cfg_web.get("exemplos", []):
                p_conj = self.raiz / ex.get("conjunto", "")
                p_alvo = self.raiz / ex.get("alvo", "")
                if p_conj.is_file() and p_alvo.is_file():
                    disponiveis.append(ex["nome"])
            return disponiveis

    def carregar_exemplo(self, nome: str) -> None:
        """Carrega conjunto e alvo de um exemplo pré-configurado."""
        with self._lock:
            disponiveis = self.exemplos_disponiveis()
            if nome not in disponiveis:
                raise ValueError(f"Exemplo indisponível: {nome}.")

            exemplo = next(ex for ex in self.cfg_web.get("exemplos", []) if ex["nome"] == nome)
            p_conj = self.raiz / exemplo["conjunto"]
            p_alvo = self.raiz / exemplo["alvo"]

            self.carregar_conjunto(p_conj.name, p_conj.read_bytes())
            self.carregar_alvo(p_alvo.name, p_alvo.read_bytes())

    def analisar(self) -> Execucao:
        """Executa a análise do conjunto carregado contra o alvo carregado."""
        with self._lock:
            if self.conjunto is None:
                raise ValueError("Carregue um conjunto de primers antes de analisar.")

            ligante = self.cfg_primers.get("desenho", {}).get("ligante", "")
            max_erros = int(self.cfg_primers.get("casamento", {}).get("max_erros", 3))
            circular = bool(self.cfg_primers.get("alvo", {}).get("circular", False))

            with tempfile.TemporaryDirectory() as tmp_dir:
                caminho_temp = _salvar_temporario(tmp_dir, self.conjunto.nome, self.conjunto.conteudo)

                if self.alvo is not None and len(self.alvo.registros) > 1:
                    try:
                        leitura_prelim = ler_conjunto(
                            caminho_temp,
                            alvo=None,
                            ligante=ligante,
                            max_erros=max_erros,
                            circular=circular,
                        )
                    except Exception as e:
                        raise _tratar_erro_leitor(e, caminho_temp, self.conjunto.nome) from e

                    primers_ref = leitura_prelim.conjunto or leitura_prelim.primers
                    escolhido = escolher_registro_alvo(list(self.alvo.registros), primers_ref)
                    alvo_efetivo = AlvoCarregado(
                        nome=self.alvo.nome,
                        descricao=escolhido.descricao,
                        sequencia=escolhido.sequencia,
                        registros=self.alvo.registros,
                    )
                    self.alvo = alvo_efetivo

                alvo_seq = self.alvo.sequencia if self.alvo else None

                try:
                    leitura = ler_conjunto(
                        caminho_temp,
                        alvo=alvo_seq,
                        ligante=ligante,
                        max_erros=max_erros,
                        circular=circular,
                    )
                except Exception as e:
                    raise _tratar_erro_leitor(e, caminho_temp, self.conjunto.nome) from e

            try:
                resultado = analisar(leitura, alvo_seq, self.cfg_primers)
            except RuntimeError as e:
                raise ValueError(str(e)) from e

            if self.alvo is not None and len(self.alvo.registros) > 1:
                desc_id = (
                    self.alvo.descricao.split()[0]
                    if self.alvo.descricao.split()
                    else self.alvo.descricao
                )
                aviso = (
                    f"FASTA com {len(self.alvo.registros)} registros; "
                    f"selecionado '{desc_id}' para a análise."
                )
                if aviso not in resultado.avisos:
                    resultado.avisos.insert(0, aviso)

            novas_viols = list(resultado.violacoes)
            if "FIP" in leitura.primers and ("F1c" not in leitura.primers or "F2" not in leitura.primers):
                if not any(v.item in ("FIP", "F1c", "F2") for v in novas_viols):
                    _, viols_fip = avaliar_primer_quimera(
                        "FIP", leitura.primers["FIP"], self.cfg_primers
                    )
                    novas_viols.extend(viols_fip)
            if "BIP" in leitura.primers and ("B1c" not in leitura.primers or "B2" not in leitura.primers):
                if not any(v.item in ("BIP", "B1c", "B2") for v in novas_viols):
                    _, viols_bip = avaliar_primer_quimera(
                        "BIP", leitura.primers["BIP"], self.cfg_primers
                    )
                    novas_viols.extend(viols_bip)
            if len(novas_viols) != len(resultado.violacoes):
                resultado = Resultado(
                    linhas=resultado.linhas,
                    heterodimeros=resultado.heterodimeros,
                    violacoes=novas_viols,
                    avisos=resultado.avisos,
                    regioes=resultado.regioes,
                    fita=resultado.fita,
                    regioes_aproximadas=resultado.regioes_aproximadas,
                    incompat=resultado.incompat,
                )

            rotulo_conjunto = Path(self.conjunto.nome).stem
            rotulo_alvo = Path(self.alvo.nome).stem if self.alvo else "—"

            execucao = Execucao(
                id=self._proximo_id,
                tipo="analise",
                criada_em=datetime.now(),
                rotulo_conjunto=rotulo_conjunto,
                rotulo_alvo=rotulo_alvo,
                alvo=self.alvo,
                leitura=leitura,
                resultado=resultado,
            )
            self._proximo_id += 1
            self.execucoes.append(execucao)
            return execucao

    def desenhar(self, relaxado: bool = False) -> Execucao:
        """Desenha conjuntos de primers para a sequência-alvo carregada (relaxado: [desenho_relaxado])."""
        with self._lock:
            if self.alvo is None:
                raise ValueError("Carregue a sequência-alvo antes de desenhar.")

            n_pb = len(self.alvo.sequencia)
            limite_pb = int(self.cfg_web.get("entrada", {}).get("alvo_desenho_max_pb", 2000))
            if n_pb > limite_pb:
                raise ValueError(
                    f"O alvo tem {milhar(n_pb)} pb; "
                    f"o desenho pela interface aceita até {milhar(limite_pb)} pb."
                )

            resultado = desenhar(self.alvo.sequencia, self.cfg_primers, relaxado=relaxado)
            rotulo_alvo = Path(self.alvo.nome).stem

            execucao = Execucao(
                id=self._proximo_id,
                tipo="desenho",
                criada_em=datetime.now(),
                rotulo_conjunto="—",
                rotulo_alvo=rotulo_alvo,
                alvo=self.alvo,
                leitura=None,
                resultado=resultado,
            )
            self._proximo_id += 1
            self.execucoes.append(execucao)
            return execucao

    def analisar_proposta(self, id_desenho: int, posicao: int) -> Execucao:
        """Analisa uma proposta gerada pelo desenho como se fosse um conjunto lido."""
        with self._lock:
            exec_desenho = self.execucao(id_desenho)
            if exec_desenho.tipo != "desenho" or not isinstance(exec_desenho.resultado, ResultadoDesenho):
                raise KeyError(f"Execução {id_desenho} não é um desenho.")

            prop: ConjuntoDesenhado | None = None
            for c in exec_desenho.resultado.conjuntos:
                if c.posicao == posicao:
                    prop = c
                    break
            if prop is None:
                raise KeyError(f"Proposta de posição {posicao} não encontrada no desenho {id_desenho}.")

            c_lamp = prop.conjunto
            primers_dict: dict[str, str] = {}
            for nome_primer in ORDEM_CANONICA:
                if nome_primer == "F3":
                    primers_dict["F3"] = c_lamp.F3
                elif nome_primer == "B3":
                    primers_dict["B3"] = c_lamp.B3
                elif nome_primer == "F2":
                    primers_dict["F2"] = c_lamp.F2
                elif nome_primer == "F1c":
                    primers_dict["F1c"] = c_lamp.F1c
                elif nome_primer == "B2":
                    primers_dict["B2"] = c_lamp.B2
                elif nome_primer == "B1c":
                    primers_dict["B1c"] = c_lamp.B1c
                elif nome_primer == "FIP":
                    primers_dict["FIP"] = c_lamp.FIP
                elif nome_primer == "BIP":
                    primers_dict["BIP"] = c_lamp.BIP
                elif nome_primer == "LF" and c_lamp.LF is not None:
                    primers_dict["LF"] = c_lamp.LF
                elif nome_primer == "LB" and c_lamp.LB is not None:
                    primers_dict["LB"] = c_lamp.LB

            leitura = Leitura(primers=primers_dict, conjunto=c_lamp, avisos=[])
            alvo_seq = exec_desenho.alvo.sequencia if exec_desenho.alvo else None
            try:
                resultado = analisar(leitura, alvo_seq, self.cfg_primers)
            except RuntimeError as e:
                raise ValueError(str(e)) from e

            rotulo_conjunto = f"proposta {posicao} do desenho {id_desenho}"
            execucao = Execucao(
                id=self._proximo_id,
                tipo="analise",
                criada_em=datetime.now(),
                rotulo_conjunto=rotulo_conjunto,
                rotulo_alvo=exec_desenho.rotulo_alvo,
                alvo=exec_desenho.alvo,
                leitura=leitura,
                resultado=resultado,
                origem="desenho",
            )
            self._proximo_id += 1
            self.execucoes.append(execucao)
            return execucao

    def execucao(self, id: int) -> Execucao:
        """Busca execução pelo identificador numérico ou levanta KeyError."""
        with self._lock:
            for ex in self.execucoes:
                if ex.id == id:
                    return ex
            raise KeyError(f"Execução com id {id} não encontrada.")

    def ultima(self, tipo: str) -> Execucao | None:
        """Retorna a execução mais recente do tipo ('analise' ou 'desenho'), ou None."""
        with self._lock:
            for ex in reversed(self.execucoes):
                if ex.tipo == tipo:
                    return ex
            return None
