"""Modelo de dados e convenção de orientação para primers LAMP."""

import itertools
import math
import re
from dataclasses import dataclass

from Bio.Seq import Seq

BASES_IUPAC: dict[str, frozenset[str]] = {
    "A": frozenset({"A"}),
    "C": frozenset({"C"}),
    "G": frozenset({"G"}),
    "T": frozenset({"T"}),
    "R": frozenset({"A", "G"}),
    "Y": frozenset({"C", "T"}),
    "S": frozenset({"C", "G"}),
    "W": frozenset({"A", "T"}),
    "K": frozenset({"G", "T"}),
    "M": frozenset({"A", "C"}),
    "B": frozenset({"C", "G", "T"}),
    "D": frozenset({"A", "G", "T"}),
    "H": frozenset({"A", "C", "T"}),
    "V": frozenset({"A", "C", "G"}),
    "N": frozenset({"A", "C", "G", "T"}),
}

ALFABETO_IUPAC: frozenset[str] = frozenset(BASES_IUPAC.keys())

_RE_5P = re.compile(r"^\s*5(?:['′’`´]\s*[-–—:]?|[-–—])\s*")
_RE_3P = re.compile(r"\s*(?:[-–—:]?\s*3['′’`´]|[-–—]\s*3)\s*$")


def limpar_notacao_extremidade(seq: str) -> str:
    """Remove marcadores opcionais de ponta 5' no início e 3' no fim da sequência."""
    s = seq.strip()
    s = _RE_5P.sub("", s)
    s = _RE_3P.sub("", s)
    return s.strip()


def normalizar(seq: str) -> str:
    """Limpa notação 5'/3' e espaços, converte U em T e valida bases IUPAC.

    Levanta ValueError indicando o caractere inválido e sua posição (base 1),
    sem despejar a sequência inteira na mensagem.
    """
    s = "".join(limpar_notacao_extremidade(seq).split()).upper()
    for pos, char in enumerate(s, start=1):
        if char not in ALFABETO_IUPAC and char != "U":
            raise ValueError(
                f"Caractere inválido '{char}' na posição {pos} "
                "(esperado base IUPAC A, C, G, T, U ou N)."
            )
    return s.replace("U", "T")


def rc(seq: str) -> str:
    """Retorna o complemento reverso da sequência normalizada usando Bio.Seq."""
    return str(Seq(normalizar(seq)).reverse_complement())


def bases_compativeis(a: str, b: str) -> bool:
    """Retorna True se duas bases normalizadas IUPAC compartilham ao menos uma base canônica."""
    return bool(BASES_IUPAC[a] & BASES_IUPAC[b])


def expandir_iupac(seq: str, max_expansoes: int = 64) -> list[str]:
    """Expande uma sequência com códigos IUPAC em variantes concretas ACGT.

    Se o produto das degenerescências for <= max_expansoes, retorna todas as combinações
    em ordem lexicográfica. Caso contrário, retorna as duas sequências extremas
    determinísticas (mais rica em A/T e mais rica em G/C).
    """
    s = normalizar(seq)
    if not s:
        return [""]
    if all(c in "ACGT" for c in s):
        return [s]

    opcoes = [sorted(BASES_IUPAC[c]) for c in s]
    total = math.prod(len(op) for op in opcoes)
    if total <= max_expansoes:
        return ["".join(p) for p in itertools.product(*opcoes)]

    seq_at = "".join(min(BASES_IUPAC[c], key=lambda b: (b in "GC", b)) for c in s)
    seq_gc = "".join(min(BASES_IUPAC[c], key=lambda b: (b in "AT", b)) for c in s)
    return [seq_at] if seq_at == seq_gc else [seq_at, seq_gc]


@dataclass(frozen=True)
class Regiao:
    """Região genômica na fita +, coordenadas base zero no intervalo semiaberto [ini, fim)."""

    ini: int
    fim: int

    def __post_init__(self) -> None:
        if self.ini < 0:
            raise ValueError(f"Coordenada ini ({self.ini}) deve ser >= 0.")
        if self.fim <= self.ini:
            raise ValueError(
                f"Coordenada fim ({self.fim}) deve ser estritamente maior que ini ({self.ini})."
            )


@dataclass(frozen=True)
class ConjuntoLAMP:
    """Conjunto de primers LAMP com orientação 5'->3' padronizada."""

    F3: str
    B3: str
    F2: str
    F1c: str
    B2: str
    B1c: str
    LF: str | None = None
    LB: str | None = None
    ligante: str = ""

    def __post_init__(self) -> None:
        for campo in ("F3", "B3", "F2", "F1c", "B2", "B1c"):
            val = getattr(self, campo)
            norm = normalizar(val)
            if not norm:
                raise ValueError(f"Sequência vazia para o primer '{campo}'.")
            object.__setattr__(self, campo, norm)
        for campo in ("LF", "LB"):
            val = getattr(self, campo)
            if val is not None:
                norm = normalizar(val)
                if not norm:
                    raise ValueError(f"Sequência vazia para o primer '{campo}'.")
                object.__setattr__(self, campo, norm)
        if self.ligante is not None:
            object.__setattr__(self, "ligante", normalizar(self.ligante))

    @property
    def FIP(self) -> str:
        """Primer interno forward: F1c + ligante + F2."""
        return self.F1c + self.ligante + self.F2

    @property
    def BIP(self) -> str:
        """Primer interno backward: B1c + ligante + B2."""
        return self.B1c + self.ligante + self.B2


def primers_de_regioes(
    alvo: str,
    regioes: dict[str, Regiao],
    ligante: str = "",
) -> ConjuntoLAMP:
    """Extrai os primers a partir das regiões na fita + do alvo conforme convenção LAMP.

    Aplica a tabela de orientação:
      F3  = seq(F3)
      B3  = rc(seq(B3c))
      F2  = seq(F2)
      F1c = rc(seq(F1))
      B1c = seq(B1c)
      B2  = rc(seq(B2c))
      LF  = rc(seq(LF))   (opcional)
      LB  = seq(LB)       (opcional)

    Levanta KeyError com o nome da chave obrigatória ausente.
    """
    obrigatorias = ("F3", "F2", "F1", "B1c", "B2c", "B3c")
    for chave in obrigatorias:
        if chave not in regioes:
            raise KeyError(chave)

    def seq(nome: str) -> str:
        r = regioes[nome]
        return alvo[r.ini : r.fim]

    f3 = seq("F3")
    b3 = rc(seq("B3c"))
    f2 = seq("F2")
    f1c = rc(seq("F1"))
    b1c = seq("B1c")
    b2 = rc(seq("B2c"))

    lf = rc(seq("LF")) if "LF" in regioes and regioes["LF"] is not None else None
    lb = seq("LB") if "LB" in regioes and regioes["LB"] is not None else None

    return ConjuntoLAMP(
        F3=f3,
        B3=b3,
        F2=f2,
        F1c=f1c,
        B2=b2,
        B1c=b1c,
        LF=lf,
        LB=lb,
        ligante=ligante,
    )
