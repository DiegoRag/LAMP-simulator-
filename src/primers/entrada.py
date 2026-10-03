"""Módulo de entrada para primers LAMP: leitura de CSV, FASTA e arquivos de configuração."""

from dataclasses import dataclass
from pathlib import Path
import tomllib
from Bio import SeqIO
from primers.modelo import ConjuntoLAMP, normalizar, rc

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
    "blp": "LB",
    "loopb": "LB",
}

ORDEM_CANONICA = ["F3", "B3", "F2", "F1c", "B2", "B1c", "FIP", "BIP", "LF", "LB"]


def carregar_config(caminho: str | Path) -> dict:
    """Carrega configuração TOML via tomllib."""
    with open(caminho, "rb") as f:
        return tomllib.load(f)


def ler_alvo(caminho: str | Path) -> str:
    """Lê sequência-alvo de arquivo FASTA contendo exatamente um registro."""
    registros = list(SeqIO.parse(caminho, "fasta"))
    if len(registros) != 1:
        raise ValueError(
            f"O arquivo {caminho} deve conter exatamente um registro FASTA, encontrado(s): {len(registros)}"
        )
    return normalizar(str(registros[0].seq))


@dataclass
class Leitura:
    """Resultado da leitura de um conjunto de primers."""

    primers: dict[str, str]
    conjunto: ConjuntoLAMP | None
    avisos: list[str]


def ler_conjunto(
    caminho: str | Path,
    alvo: str | None = None,
    ligante: str = "",
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

    if ext == ".csv":
        linhas = p.read_text(encoding="utf-8-sig").splitlines()
        header_encontrado = False
        for num_linha, linha in enumerate(linhas, start=1):
            linha_limpa = linha.strip()
            if not linha_limpa:
                continue

            partes = [col.strip() for col in linha_limpa.split(",")]
            if not header_encontrado:
                partes_lower = [col.lower() for col in partes]
                if len(partes_lower) != 2 or partes_lower[0] != "nome" or partes_lower[1] != "sequencia":
                    raise ValueError(f"Cabeçalho inválido na linha {num_linha}: esperado 'nome,sequencia'")
                header_encontrado = True
                continue

            if len(partes) != 2:
                raise ValueError(
                    f"Linha {num_linha} inválida: esperados 2 campos (nome, sequencia), recebidos {len(partes)}"
                )

            nome_raw, seq_raw = partes[0], partes[1]
            nome_key = nome_raw.lower()
            if nome_key not in MAPA_NOMES:
                raise ValueError(f"Nome de primer desconhecido '{nome_raw}' na linha {num_linha}")

            nome_canonico = MAPA_NOMES[nome_key]
            if nome_canonico in primers_lidos:
                raise ValueError(f"Nome de primer repetido '{nome_canonico}' na linha {num_linha}")

            try:
                seq_norm = normalizar(seq_raw)
            except ValueError as e:
                raise ValueError(f"Sequência inválida na linha {num_linha}: {e}") from e

            primers_lidos[nome_canonico] = seq_norm

        if not header_encontrado:
            raise ValueError("Arquivo CSV vazio ou sem cabeçalho 'nome,sequencia'")

    else:
        registros = list(SeqIO.parse(caminho, "fasta"))
        for reg in registros:
            nome_raw = reg.id
            nome_key = nome_raw.lower()
            if nome_key not in MAPA_NOMES:
                raise ValueError(f"Nome de primer desconhecido no registro '{nome_raw}'")

            nome_canonico = MAPA_NOMES[nome_key]
            if nome_canonico in primers_lidos:
                raise ValueError(f"Nome de primer repetido '{nome_canonico}' no registro '{nome_raw}'")

            try:
                seq_norm = normalizar(str(reg.seq))
            except ValueError as e:
                raise ValueError(f"Sequência inválida no registro '{nome_raw}': {e}") from e

            primers_lidos[nome_canonico] = seq_norm

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

    if tem_f1c and tem_f2:
        fip_esperado = primers_lidos["F1c"] + ligante + primers_lidos["F2"]
        if "FIP" in primers_lidos:
            if primers_lidos["FIP"] != fip_esperado:
                raise ValueError(
                    f"FIP informado ({primers_lidos['FIP']}) difere de F1c + ligante + F2 ({fip_esperado})"
                )
        primers_lidos["FIP"] = fip_esperado

    if tem_b1c and tem_b2:
        bip_esperado = primers_lidos["B1c"] + ligante + primers_lidos["B2"]
        if "BIP" in primers_lidos:
            if primers_lidos["BIP"] != bip_esperado:
                raise ValueError(
                    f"BIP informado ({primers_lidos['BIP']}) difere de B1c + ligante + B2 ({bip_esperado})"
                )
        primers_lidos["BIP"] = bip_esperado

    if "FIP" in primers_lidos and not tem_f1c:
        fip_seq = primers_lidos["FIP"]
        if not alvo:
            avisos.append("FIP sem F1c/F2 informado sem alvo para separação.")
        else:
            cortes_validos = []
            max_i = len(fip_seq) - 10 - len(ligante)
            for i in range(10, max_i + 1):
                p5 = fip_seq[:i]
                lig = fip_seq[i : i + len(ligante)]
                p3 = fip_seq[i + len(ligante) :]
                if lig != ligante:
                    continue
                pos_rc_p5 = alvo.find(rc(p5))
                pos_p3 = alvo.find(p3)
                if pos_rc_p5 >= 0 and pos_p3 >= 0 and pos_rc_p5 > pos_p3:
                    cortes_validos.append(i)

            if not cortes_validos:
                avisos.append("FIP não pôde ser separado a partir do alvo.")
            else:
                menor_i = min(cortes_validos)
                primers_lidos["F1c"] = fip_seq[:menor_i]
                primers_lidos["F2"] = fip_seq[menor_i + len(ligante) :]
                if len(cortes_validos) > 1:
                    avisos.append(
                        f"Corte do FIP ambíguo: cortes válidos em {cortes_validos}; "
                        f"usando menor i={menor_i}. Recomenda-se informar as partes separadas."
                    )

    if "BIP" in primers_lidos and not tem_b1c:
        bip_seq = primers_lidos["BIP"]
        if not alvo:
            avisos.append("BIP sem B1c/B2 informado sem alvo para separação.")
        else:
            cortes_validos = []
            max_i = len(bip_seq) - 10 - len(ligante)
            for i in range(10, max_i + 1):
                p5 = bip_seq[:i]
                lig = bip_seq[i : i + len(ligante)]
                p3 = bip_seq[i + len(ligante) :]
                if lig != ligante:
                    continue
                pos_p5 = alvo.find(p5)
                pos_rc_p3 = alvo.find(rc(p3))
                if pos_p5 >= 0 and pos_rc_p3 >= 0 and pos_rc_p3 > pos_p5:
                    cortes_validos.append(i)

            if not cortes_validos:
                avisos.append("BIP não pôde ser separado a partir do alvo.")
            else:
                menor_i = min(cortes_validos)
                primers_lidos["B1c"] = bip_seq[:menor_i]
                primers_lidos["B2"] = bip_seq[menor_i + len(ligante) :]
                if len(cortes_validos) > 1:
                    avisos.append(
                        f"Corte do BIP ambíguo: cortes válidos em {cortes_validos}; "
                        f"usando menor i={menor_i}. Recomenda-se informar as partes separadas."
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
            ligante=ligante,
        )
    else:
        conjunto = None
        avisos.append("Conjunto incompleto: faltam primers essenciais (F3, B3, F2, F1c, B2, B1c).")

    return Leitura(primers=primers_ordenados, conjunto=conjunto, avisos=avisos)
