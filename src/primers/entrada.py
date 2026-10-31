"""Módulo de entrada para primers LAMP: leitura de CSV, FASTA e arquivos de configuração."""

import csv
import io
import re
import tomllib
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO

from primers.casamento import _so_acgt, ocorrencias_aproximadas, ocorrencias_compativeis
from primers.modelo import ConjuntoLAMP, Regiao, normalizar, rc

LIGANTES_PADRAO: tuple[str, ...] = ("", "TTTT", "TTT", "TTTTT")

MAPA_NOMES = {
    "f3": "F3",
    "b3": "B3",
    "f2": "F2",
    "f1c": "F1c",
    "b2": "B2",
    "b1c": "B1c",
    "fip": "FIP",
    "bip": "BIP",
    "lf": "LF",
    "lb": "LB",
    "flp": "LF",
    "loopf": "LF",
    "loop-f": "LF",
    "loop_f": "LF",
    "blp": "LB",
    "loopb": "LB",
    "loop-b": "LB",
    "loop_b": "LB",
}

ORDEM_CANONICA = ["F3", "B3", "F2", "F1c", "B2", "B1c", "FIP", "BIP", "LF", "LB"]

CABECALHOS_NOME = {"nome", "name", "primer", "id", "oligo", "identificador"}
CABECALHOS_SEQ = {
    "sequencia",
    "sequência",
    "sequence",
    "oligoseq",
    "sequencia 5'-3'",
    "sequencia 5'->3'",
    "sequence 5'-3'",
}


@dataclass(frozen=True)
class RegistroAlvo:
    """Registro individual lido de um arquivo FASTA de alvo."""

    id: str
    descricao: str
    sequencia: str


def reconhecer_nome_primer(nome_raw: str) -> str | None:
    """Reconhece o papel canônico de um primer a partir de um nome flexível."""
    s = nome_raw.strip().strip("\"'").strip()
    if not s:
        return None
    s_low = s.lower()
    if s_low in MAPA_NOMES:
        return MAPA_NOMES[s_low]

    raw_tokens = [t for t in re.split(r"[-_./:\s()\[\]]+", s_low) if t]
    tokens: list[str] = []
    i = 0
    while i < len(raw_tokens):
        if i + 1 < len(raw_tokens) and raw_tokens[i] == "loop" and raw_tokens[i + 1] in ("f", "b"):
            tokens.append(f"loop{raw_tokens[i + 1]}")
            i += 2
        else:
            tokens.append(raw_tokens[i])
            i += 1

    papeis = {MAPA_NOMES[t] for t in tokens if t in MAPA_NOMES}
    if not papeis:
        padrao = r"(f1c|b1c|fip|bip|loopf|loopb|flp|blp|f3|b3|f2|b2|lf|lb)v?\d+"
        for t in tokens:
            m = re.fullmatch(padrao, t)
            if m:
                papeis.add(MAPA_NOMES[m.group(1)])

    if len(papeis) == 1:
        return next(iter(papeis))
    return None


def carregar_config(caminho: str | Path) -> dict:
    """Carrega configuração TOML via tomllib."""
    with open(caminho, "rb") as f:
        return tomllib.load(f)


def _ler_registros_fasta(caminho: str | Path) -> list:
    """Lê registros de arquivo FASTA validando cabeçalho, remoção de BOM e não-vazio."""
    texto = Path(caminho).read_text(encoding="utf-8-sig").lstrip("\ufeff")
    if not texto.strip():
        raise ValueError(f"Arquivo FASTA vazio: {caminho}")
    primeira = next(ln.strip() for ln in texto.splitlines() if ln.strip())
    if not primeira.startswith(">"):
        raise ValueError(f"Arquivo FASTA inválido (falta cabeçalho '>'): {caminho}")
    registros = list(SeqIO.parse(io.StringIO(texto), "fasta"))
    if not registros:
        raise ValueError(f"Nenhum registro FASTA encontrado em {caminho}")
    for reg in registros:
        if not reg.id.strip():
            raise ValueError(f"Cabeçalho FASTA sem identificador em {caminho}")
    return registros


def ler_alvos(caminho: str | Path) -> list[RegistroAlvo]:
    """Lê todos os registros de um arquivo FASTA de alvo normalizando as sequências."""
    regs = _ler_registros_fasta(caminho)
    resultado: list[RegistroAlvo] = []
    for reg in regs:
        try:
            seq_norm = normalizar(str(reg.seq))
        except ValueError as e:
            raise ValueError(f"Sequência inválida no registro '{reg.id}': {e}") from e
        if not seq_norm:
            raise ValueError(f"Sequência vazia no registro '{reg.id}' do arquivo {caminho}")
        resultado.append(RegistroAlvo(id=reg.id, descricao=reg.description, sequencia=seq_norm))
    return resultado


def _contem_compativel(fita: str, s: str) -> bool:
    """Verifica se s ocorre em fita considerando compatibilidade IUPAC."""
    if not s or not fita or len(s) > len(fita):
        return False
    if _so_acgt(s) and _so_acgt(fita):
        return fita.find(s) != -1
    return bool(ocorrencias_compativeis(fita, s))


def escolher_registro_alvo(
    registros: list[RegistroAlvo],
    primers_ou_conjunto: ConjuntoLAMP | dict[str, str] | None = None,
    *,
    registro: str | int | None = None,
) -> RegistroAlvo:
    """Seleciona um registro alvo por índice, id ou maior compatibilidade com primers."""
    if not registros:
        raise ValueError("Lista de registros de alvo vazia.")
    if registro is not None:
        if isinstance(registro, int):
            if 0 <= registro < len(registros):
                return registros[registro]
            raise ValueError(
                f"Índice de registro {registro} fora do intervalo (0..{len(registros) - 1})."
            )
        if isinstance(registro, str):
            for r in registros:
                if r.id == registro or r.id.split(".")[0] == registro or registro in r.descricao:
                    return r
            raise ValueError(f"Registro '{registro}' não encontrado no FASTA.")

    if len(registros) == 1 or primers_ou_conjunto is None:
        return registros[0]

    seqs_busca: list[str] = []
    if isinstance(primers_ou_conjunto, ConjuntoLAMP):
        seqs_busca = [
            primers_ou_conjunto.F3,
            primers_ou_conjunto.F2,
            rc(primers_ou_conjunto.F1c),
            primers_ou_conjunto.B1c,
            rc(primers_ou_conjunto.B2),
            rc(primers_ou_conjunto.B3),
        ]
        if primers_ou_conjunto.LF:
            seqs_busca.append(rc(primers_ou_conjunto.LF))
        if primers_ou_conjunto.LB:
            seqs_busca.append(primers_ou_conjunto.LB)
    elif isinstance(primers_ou_conjunto, dict):
        for k in ("F3", "B3", "F2", "F1c", "B2", "B1c", "LF", "LB"):
            if k in primers_ou_conjunto:
                seqs_busca.append(normalizar(primers_ou_conjunto[k]))
        for k in ("FIP", "BIP"):
            if k in primers_ou_conjunto:
                s_int = normalizar(primers_ou_conjunto[k])
                seqs_busca.append(s_int)
                if len(s_int) >= 36:
                    seqs_busca.extend([s_int[:18], s_int[-18:]])

    melhor_idx = 0
    melhor_score = (-1, 0)
    for idx, reg in enumerate(registros):
        fita_mais = reg.sequencia
        fita_menos = rc(reg.sequencia)
        acertos_mais = sum(1 for s in seqs_busca if _contem_compativel(fita_mais, s))
        acertos_menos = sum(1 for s in seqs_busca if _contem_compativel(fita_menos, s))
        acertos = max(acertos_mais, acertos_menos)
        score = (acertos, -idx)
        if score > melhor_score:
            melhor_score = score
            melhor_idx = idx

    return registros[melhor_idx]


def ler_alvo(
    caminho: str | Path,
    *,
    registro: str | int | None = None,
    conjunto: ConjuntoLAMP | dict[str, str] | None = None,
) -> str:
    """Lê sequência-alvo de arquivo FASTA (1 ou vários registros)."""
    registros = ler_alvos(caminho)
    escolhido = escolher_registro_alvo(registros, conjunto, registro=registro)
    return escolhido.sequencia


@dataclass
class Leitura:
    """Resultado da leitura de um conjunto de primers."""

    primers: dict[str, str]
    conjunto: ConjuntoLAMP | None
    avisos: list[str]
    separados: tuple[str, ...] = ()  # "FIP" e/ou "BIP" lidos inteiros e separados nas partes pelo alvo


@dataclass(frozen=True)
class _CandidatoCorte:
    i: int
    lig_cand: str
    fita_id: str
    p5: str
    p3: str
    r_esq: Regiao
    r_dir: Regiao
    e_esq: int
    e_dir: int
    dist_loop: int
    score_geom: tuple
    chave: tuple


def _buscar_ancora(
    seq: str,
    s_busca: str,
    max_e: int,
    max_ini: int,
) -> list[tuple[int, int]]:
    """Encontra ocorrências de seq em s_busca até max_e erros."""
    k = len(seq)
    if k == 0 or k > len(s_busca):
        return []
    lim = min(max_ini, len(s_busca) - k + 1)
    exatas: list[tuple[int, int]] = []
    p = 0
    while p < lim:
        idx = s_busca.find(seq, p, len(s_busca))
        if idx == -1 or idx >= lim:
            break
        exatas.append((idx, 0))
        p = idx + 1
    if exatas:
        return exatas

    tem_iupac = not (_so_acgt(seq) and _so_acgt(s_busca))
    if tem_iupac:
        exatas_iupac = [(r.ini, 0) for r in ocorrencias_compativeis(s_busca, seq, limite_ini=lim)]
        if exatas_iupac:
            return exatas_iupac

    return ocorrencias_aproximadas(seq, s_busca, max_e, lim)


def ler_conjunto(
    caminho: str | Path,
    alvo: str | None = None,
    ligante: str = "",
    *,
    max_erros: int = 3,
    circular: bool = False,
) -> Leitura:
    """Lê conjunto de primers de arquivo CSV ou FASTA.

    Parâmetros:
        caminho: caminho para arquivo .csv, .fasta, .fa ou .fna.
        alvo: sequência-alvo normalizada (opcional, para separar FIP/BIP).
        ligante: sequência do ligante (opcional).

    Retorna:
        Instância de Leitura com primers na ordem canônica, conjunto e avisos.
    """
    p = Path(caminho)
    ext = p.suffix.lower()
    if ext not in (".csv", ".fasta", ".fa", ".fna"):
        raise ValueError(f"Extensão não suportada: '{p.suffix}'. Esperado .csv, .fasta, .fa ou .fna")

    primers_lidos: dict[str, str] = {}
    avisos: list[str] = []
    separados: list[str] = []

    if ext == ".csv":
        texto = p.read_text(encoding="utf-8-sig").lstrip("\ufeff")
        if not texto.strip():
            raise ValueError("Arquivo CSV vazio ou sem cabeçalho 'nome,sequencia'")
        linhas = texto.splitlines()
        primeira_ln = next((ln.strip() for ln in linhas if ln.strip()), "")
        if ";" in primeira_ln and primeira_ln.count(";") > primeira_ln.count(","):
            delim = ";"
        elif "\t" in primeira_ln and "," not in primeira_ln:
            delim = "\t"
        else:
            delim = ","

        header_encontrado = False
        col_nome, col_seq = 0, 1
        for num_linha, linha in enumerate(linhas, start=1):
            linha_limpa = linha.strip()
            if not linha_limpa:
                continue

            partes = [
                col.strip().strip("\"'").strip()
                for col in next(
                    csv.reader(
                        [linha_limpa],
                        delimiter=delim,
                        quotechar='"',
                        skipinitialspace=True,
                    )
                )
            ]
            while len(partes) > 2 and partes[-1] == "":
                partes.pop()

            if not header_encontrado:
                partes_lower = [col.lower() for col in partes]
                if (
                    len(partes_lower) == 2
                    and partes_lower[0] in CABECALHOS_NOME
                    and partes_lower[1] in CABECALHOS_SEQ
                ):
                    col_nome, col_seq = 0, 1
                elif (
                    len(partes_lower) == 2
                    and partes_lower[0] in CABECALHOS_SEQ
                    and partes_lower[1] in CABECALHOS_NOME
                ):
                    col_nome, col_seq = 1, 0
                else:
                    raise ValueError(
                        f"Cabeçalho inválido na linha {num_linha}: esperado 'nome,sequencia'"
                    )
                header_encontrado = True
                continue

            if len(partes) != 2:
                raise ValueError(
                    f"Linha {num_linha} inválida: esperados 2 campos (nome, sequencia), "
                    f"recebidos {len(partes)}"
                )

            nome_raw, seq_raw = partes[col_nome], partes[col_seq]
            nome_canonico = reconhecer_nome_primer(nome_raw)
            if nome_canonico is None:
                raise ValueError(f"Nome de primer desconhecido '{nome_raw}' na linha {num_linha}")
            if nome_canonico in primers_lidos:
                raise ValueError(f"Nome de primer repetido '{nome_canonico}' na linha {num_linha}")

            try:
                seq_norm = normalizar(seq_raw)
            except ValueError as e:
                raise ValueError(f"Sequência inválida na linha {num_linha}: {e}") from e
            if not seq_norm:
                raise ValueError(f"Sequência vazia na linha {num_linha} para o primer '{nome_raw}'")

            primers_lidos[nome_canonico] = seq_norm

        if not header_encontrado:
            raise ValueError("Arquivo CSV vazio ou sem cabeçalho 'nome,sequencia'")
        if not primers_lidos:
            raise ValueError("Nenhum primer encontrado no arquivo CSV.")

    else:
        registros = _ler_registros_fasta(caminho)
        for reg in registros:
            nome_raw = reg.id
            nome_canonico = reconhecer_nome_primer(nome_raw)
            if nome_canonico is None:
                raise ValueError(f"Nome de primer desconhecido no registro '{nome_raw}'")
            if nome_canonico in primers_lidos:
                raise ValueError(f"Nome de primer repetido '{nome_canonico}' no registro '{nome_raw}'")

            try:
                seq_norm = normalizar(str(reg.seq))
            except ValueError as e:
                raise ValueError(f"Sequência inválida no registro '{nome_raw}': {e}") from e
            if not seq_norm:
                raise ValueError(f"Sequência vazia no registro '{nome_raw}'")

            primers_lidos[nome_canonico] = seq_norm

        if not primers_lidos:
            raise ValueError(f"Nenhum primer encontrado no arquivo FASTA: {caminho}")

    tem_f1c = "F1c" in primers_lidos
    tem_f2 = "F2" in primers_lidos
    if tem_f1c != tem_f2:
        if tem_f1c:
            raise ValueError("F1c informado sem F2")
        raise ValueError("F2 informado sem F1c")

    tem_b1c = "B1c" in primers_lidos
    tem_b2 = "B2" in primers_lidos
    if tem_b1c != tem_b2:
        if tem_b1c:
            raise ValueError("B1c informado sem B2")
        raise ValueError("B2 informado sem B1c")

    # Detecção de ligante quando partes e primers inteiros são informados juntos
    lig_fip: str | None = None
    if tem_f1c and tem_f2 and "FIP" in primers_lidos:
        fip_esperado = primers_lidos["F1c"] + ligante + primers_lidos["F2"]
        if primers_lidos["FIP"] == fip_esperado:
            lig_fip = ligante
        elif (
            primers_lidos["FIP"].startswith(primers_lidos["F1c"])
            and primers_lidos["FIP"].endswith(primers_lidos["F2"])
            and len(primers_lidos["FIP"]) >= len(primers_lidos["F1c"]) + len(primers_lidos["F2"])
        ):
            meio_fip = primers_lidos["FIP"][
                len(primers_lidos["F1c"]) : len(primers_lidos["FIP"]) - len(primers_lidos["F2"])
            ]
            if meio_fip in LIGANTES_PADRAO:
                lig_fip = meio_fip
            else:
                raise ValueError(
                    f"FIP informado ({primers_lidos['FIP']}) difere de F1c + ligante + F2 ({fip_esperado})"
                )
        else:
            raise ValueError(
                f"FIP informado ({primers_lidos['FIP']}) difere de F1c + ligante + F2 ({fip_esperado})"
            )

    lig_bip: str | None = None
    if tem_b1c and tem_b2 and "BIP" in primers_lidos:
        bip_esperado = primers_lidos["B1c"] + ligante + primers_lidos["B2"]
        if primers_lidos["BIP"] == bip_esperado:
            lig_bip = ligante
        elif (
            primers_lidos["BIP"].startswith(primers_lidos["B1c"])
            and primers_lidos["BIP"].endswith(primers_lidos["B2"])
            and len(primers_lidos["BIP"]) >= len(primers_lidos["B1c"]) + len(primers_lidos["B2"])
        ):
            meio_bip = primers_lidos["BIP"][
                len(primers_lidos["B1c"]) : len(primers_lidos["BIP"]) - len(primers_lidos["B2"])
            ]
            if meio_bip in LIGANTES_PADRAO:
                lig_bip = meio_bip
            else:
                raise ValueError(
                    f"BIP informado ({primers_lidos['BIP']}) difere de B1c + ligante + B2 ({bip_esperado})"
                )
        else:
            raise ValueError(
                f"BIP informado ({primers_lidos['BIP']}) difere de B1c + ligante + B2 ({bip_esperado})"
            )

    if lig_fip is not None and lig_bip is not None and lig_fip != lig_bip:
        raise ValueError(f"Ligante do FIP ({lig_fip!r}) difere do ligante do BIP ({lig_bip!r})")

    ligante_efetivo = ligante
    if lig_fip is not None:
        ligante_efetivo = lig_fip
    elif lig_bip is not None:
        ligante_efetivo = lig_bip

    if tem_f1c and tem_f2:
        primers_lidos["FIP"] = primers_lidos["F1c"] + ligante_efetivo + primers_lidos["F2"]

    if tem_b1c and tem_b2:
        primers_lidos["BIP"] = primers_lidos["B1c"] + ligante_efetivo + primers_lidos["B2"]

    precisa_separar_fip = ("FIP" in primers_lidos) and (not tem_f1c)
    precisa_separar_bip = ("BIP" in primers_lidos) and (not tem_b1c)

    if precisa_separar_fip or precisa_separar_bip:
        if not alvo:
            if precisa_separar_fip:
                avisos.append("FIP sem F1c/F2 informado sem alvo para separação.")
            if precisa_separar_bip:
                avisos.append("BIP sem B1c/B2 informado sem alvo para separação.")
        else:
            fitas_candidatas = [("+", alvo), ("-", rc(alvo))]
            fita_preferida: str | None = None
            fitas_dados: list[dict] = []
            ancorou_f3_b3: dict[str, bool] = {}

            for fita_id, s_fita in fitas_candidatas:
                s_busca = (s_fita + s_fita) if (circular and len(s_fita) > 0) else s_fita
                lim_ini = len(s_fita) if circular else len(s_busca)

                tem_f3 = "F3" in primers_lidos
                tem_b3 = "B3" in primers_lidos

                r_f3: Regiao | None = None
                r_b3c: Regiao | None = None
                ancorado = False
                pos_f3_escolhida: int | None = None
                pos_b3c_escolhida: int | None = None

                if tem_f3 and tem_b3:
                    f3_seq = primers_lidos["F3"]
                    b3c_seq = rc(primers_lidos["B3"])
                    ocorr_f3 = _buscar_ancora(f3_seq, s_busca, max_erros, lim_ini)
                    ocorr_b3c = _buscar_ancora(b3c_seq, s_busca, max_erros, len(s_busca) - len(b3c_seq) + 1)

                    pares_ancora = []
                    for p_f, e_f in ocorr_f3:
                        for p_b, e_b in ocorr_b3c:
                            if p_f < p_b:
                                pares_ancora.append((e_f + e_b, p_b - p_f, p_f, p_b))

                    if pares_ancora:
                        pares_ancora.sort()
                        _, _, p_f, p_b = pares_ancora[0]
                        ancorado = True
                        pos_f3_escolhida = p_f
                        pos_b3c_escolhida = p_b
                        r_f3 = Regiao(p_f, p_f + len(f3_seq))
                        r_b3c = Regiao(p_b, p_b + len(b3c_seq))
                elif tem_f3:
                    f3_seq = primers_lidos["F3"]
                    ocorr_f3 = _buscar_ancora(f3_seq, s_busca, max_erros, lim_ini)
                    if ocorr_f3:
                        p_f = ocorr_f3[0][0]
                        r_f3 = Regiao(p_f, p_f + len(f3_seq))
                elif tem_b3:
                    b3c_seq = rc(primers_lidos["B3"])
                    ocorr_b3c = _buscar_ancora(b3c_seq, s_busca, max_erros, len(s_busca) - len(b3c_seq) + 1)
                    if ocorr_b3c:
                        p_b = ocorr_b3c[0][0]
                        r_b3c = Regiao(p_b, p_b + len(b3c_seq))

                ancorou_f3_b3[fita_id] = ancorado

                if (
                    ancorado
                    and pos_f3_escolhida is not None
                    and pos_b3c_escolhida is not None
                    and (0 <= pos_b3c_escolhida - pos_f3_escolhida <= 400)
                ):
                    jan_ini = max(0, pos_f3_escolhida)
                    jan_fim = min(len(s_busca), pos_b3c_escolhida + len(primers_lidos["B3"]))
                else:
                    jan_ini = 0
                    jan_fim = len(s_busca)

                fitas_dados.append(
                    {
                        "fita_id": fita_id,
                        "s_busca": s_busca,
                        "jan_ini": jan_ini,
                        "jan_fim": jan_fim,
                        "r_f3": r_f3,
                        "r_b3c": r_b3c,
                    }
                )

            if ancorou_f3_b3.get("+") and not ancorou_f3_b3.get("-"):
                fita_preferida = "+"
            elif ancorou_f3_b3.get("-") and not ancorou_f3_b3.get("+"):
                fita_preferida = "-"

            for tipo, nome_int in [("FIP", "FIP"), ("BIP", "BIP")]:
                if tipo == "FIP" and not precisa_separar_fip:
                    continue
                if tipo == "BIP" and not precisa_separar_bip:
                    continue

                seq_int = primers_lidos[tipo]
                cands_lig = (ligante_efetivo, *[x for x in LIGANTES_PADRAO if x != ligante_efetivo])

                candidatos: list[_CandidatoCorte] = []

                for lig_cand in cands_lig:
                    len_lig = len(lig_cand)
                    for i in range(10, len(seq_int) - 10 - len_lig + 1):
                        if len_lig > 0 and seq_int[i : i + len_lig] != lig_cand:
                            continue
                        p5 = seq_int[:i]
                        p3 = seq_int[i + len_lig :]

                        if tipo == "FIP":
                            seq_esq = p3  # F2
                            seq_dir = rc(p5)  # F1
                        else:
                            seq_esq = p5  # B1c
                            seq_dir = rc(p3)  # B2c

                        k_esq = len(seq_esq)
                        k_dir = len(seq_dir)

                        for fd in fitas_dados:
                            fita_id = fd["fita_id"]
                            s_busca = fd["s_busca"]
                            jan_ini = fd["jan_ini"]
                            jan_fim = fd["jan_fim"]
                            r_f3 = fd["r_f3"]
                            r_b3c = fd["r_b3c"]

                            tem_iupac_esq = not (_so_acgt(seq_esq) and _so_acgt(s_busca))
                            pos_esq_exatas: list[int] = []
                            if not tem_iupac_esq:
                                curr = jan_ini
                                while curr + k_esq <= jan_fim:
                                    idx = s_busca.find(seq_esq, curr, jan_fim)
                                    if idx == -1:
                                        break
                                    if not circular or idx < len(alvo):
                                        pos_esq_exatas.append(idx)
                                    curr = idx + 1
                            else:
                                for r in ocorrencias_compativeis(s_busca[jan_ini:jan_fim], seq_esq):
                                    idx = jan_ini + r.ini
                                    if not circular or idx < len(alvo):
                                        pos_esq_exatas.append(idx)

                            tem_iupac_dir = not (_so_acgt(seq_dir) and _so_acgt(s_busca))
                            pos_dir_exatas: list[int] = []
                            if not tem_iupac_dir:
                                curr = jan_ini
                                while curr + k_dir <= jan_fim:
                                    idx = s_busca.find(seq_dir, curr, jan_fim)
                                    if idx == -1:
                                        break
                                    pos_dir_exatas.append(idx)
                                    curr = idx + 1
                            else:
                                for r in ocorrencias_compativeis(s_busca[jan_ini:jan_fim], seq_dir):
                                    pos_dir_exatas.append(jan_ini + r.ini)

                            pares_encontrados: list[tuple[Regiao, Regiao, int, int, int]] = []
                            for p_e in pos_esq_exatas:
                                for p_d in pos_dir_exatas:
                                    if p_e < p_d:
                                        r_e = Regiao(p_e, p_e + k_esq)
                                        r_d = Regiao(p_d, p_d + k_dir)
                                        d_loop = (
                                            (r_d.ini - r_e.ini)
                                            if tipo == "FIP"
                                            else (r_d.fim - r_e.fim)
                                        )
                                        pares_encontrados.append((r_e, r_d, 0, 0, d_loop))

                            if not pares_encontrados and max_erros > 0:
                                lim_esq = min(
                                    jan_fim - k_esq + 1,
                                    len(alvo) if circular else jan_fim - k_esq + 1,
                                )
                                cands_esq = sorted(
                                    (
                                        (jan_ini + p, e)
                                        for p, e in ocorrencias_aproximadas(
                                            seq_esq,
                                            s_busca[jan_ini:jan_fim],
                                            max_erros,
                                            lim_esq - jan_ini,
                                        )
                                    ),
                                    key=lambda t: t[0],
                                )
                                erros_dir = {
                                    i + jan_ini: e
                                    for i, e in ocorrencias_aproximadas(
                                        seq_dir,
                                        s_busca[jan_ini:jan_fim],
                                        max_erros,
                                        jan_fim - jan_ini - k_dir + 1,
                                    )
                                }
                                for p_e, e_e in cands_esq:
                                    for d_loop in range(25, 86):
                                        if tipo == "FIP":
                                            p_d = p_e + d_loop
                                        else:
                                            p_d = p_e + k_esq - k_dir + d_loop
                                        if p_d <= p_e or p_d + k_dir > jan_fim:
                                            continue
                                        e_d = erros_dir.get(p_d)
                                        if e_d is None:
                                            continue
                                        if e_e + e_d <= max_erros:
                                            r_e = Regiao(p_e, p_e + k_esq)
                                            r_d = Regiao(p_d, p_d + k_dir)
                                            pares_encontrados.append((r_e, r_d, e_e, e_d, d_loop))

                            for r_e, r_d, e_e, e_d, d_loop in pares_encontrados:
                                fita_pen = (
                                    0 if (fita_preferida is None or fita_id == fita_preferida) else 1
                                )
                                e_total = e_e + e_d
                                sobreposicao = max(0, r_e.fim - r_d.ini)
                                pen_dist = max(0, 40 - d_loop, d_loop - 60)
                                pen_comp = (
                                    max(0, 18 - len(p5), len(p5) - 25)
                                    + max(0, 18 - len(p3), len(p3) - 25)
                                )
                                if tipo == "FIP":
                                    if r_f3 is not None:
                                        gap = r_e.ini - r_f3.fim
                                        pen_gap_ext = max(0, 0 - gap, gap - 60)
                                    else:
                                        pen_gap_ext = 0
                                else:
                                    if r_b3c is not None:
                                        gap = r_b3c.ini - r_d.fim
                                        pen_gap_ext = max(0, 0 - gap, gap - 60)
                                    else:
                                        pen_gap_ext = 0
                                pref_lig = (
                                    0
                                    if lig_cand == ligante_efetivo
                                    else (1 if lig_cand == "" else 2)
                                )
                                desvio_geom = abs(d_loop - 50) + abs(len(p5) - len(p3))
                                pref_fita_mais = 0 if fita_id == "+" else 1

                                score_geom = (
                                    fita_pen,
                                    e_total,
                                    sobreposicao,
                                    pen_dist,
                                    pen_comp,
                                    pen_gap_ext,
                                    pref_lig,
                                    desvio_geom,
                                    pref_fita_mais,
                                )
                                chave = (*score_geom, i)
                                candidatos.append(
                                    _CandidatoCorte(
                                        i=i,
                                        lig_cand=lig_cand,
                                        fita_id=fita_id,
                                        p5=p5,
                                        p3=p3,
                                        r_esq=r_e,
                                        r_dir=r_d,
                                        e_esq=e_e,
                                        e_dir=e_d,
                                        dist_loop=d_loop,
                                        score_geom=score_geom,
                                        chave=chave,
                                    )
                                )

                if not candidatos:
                    avisos.append(f"{nome_int} não pôde ser separado a partir do alvo.")
                else:
                    melhor = min(candidatos, key=lambda c: c.chave)
                    separados.append(nome_int)
                    if tipo == "FIP":
                        primers_lidos["F1c"] = melhor.p5
                        primers_lidos["F2"] = melhor.p3
                    else:
                        primers_lidos["B1c"] = melhor.p5
                        primers_lidos["B2"] = melhor.p3

                    ligante_efetivo = melhor.lig_cand

                    if melhor.score_geom[1] > 0:
                        e_tot = melhor.score_geom[1]
                        avisos.append(
                            f"{nome_int} separado por casamento aproximado "
                            f"({e_tot} diferença(s) no alvo) em i={melhor.i}."
                        )

                    cortes_exatos = sorted(
                        {
                            c.i
                            for c in candidatos
                            if c.fita_id == melhor.fita_id
                            and c.lig_cand == melhor.lig_cand
                            and c.score_geom[1] == 0
                            and c.score_geom[2] == 0
                        }
                    )
                    if len(cortes_exatos) > 1:
                        melhores_score = [c for c in candidatos if c.score_geom == melhor.score_geom]
                        if len(melhores_score) > 1:
                            avisos.append(
                                f"Corte do {nome_int} ambíguo (empate geométrico em {cortes_exatos}); "
                                f"usando i={melhor.i}. Informe as partes separadas para desambiguar."
                            )
                        else:
                            avisos.append(
                                f"Corte do {nome_int} ambíguo: cortes válidos em {cortes_exatos}; "
                                f"escolhido i={melhor.i} pela geometria. "
                                "Recomenda-se informar as partes separadas."
                            )

    primers_ordenados = {k: primers_lidos[k] for k in ORDEM_CANONICA if k in primers_lidos}

    essenciais = ("F3", "B3", "F2", "F1c", "B2", "B1c")
    if all(k in primers_ordenados for k in essenciais):
        conjunto = ConjuntoLAMP(
            F3=primers_ordenados["F3"],
            B3=primers_ordenados["B3"],
            F2=primers_ordenados["F2"],
            F1c=primers_ordenados["F1c"],
            B2=primers_ordenados["B2"],
            B1c=primers_ordenados["B1c"],
            LF=primers_ordenados.get("LF"),
            LB=primers_ordenados.get("LB"),
            ligante=ligante_efetivo,
        )
    else:
        conjunto = None
        avisos.append("Conjunto incompleto: faltam primers essenciais (F3, B3, F2, F1c, B2, B1c).")

    return Leitura(primers=primers_ordenados, conjunto=conjunto, avisos=avisos, separados=tuple(separados))
