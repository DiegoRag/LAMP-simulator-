"""Modelo de dados e convenção de orientação para primers LAMP."""

from dataclasses import dataclass
from Bio.Seq import Seq


def normalizar(seq: str) -> str:
    """Remove espaços em branco, converte para maiúsculas e valida bases ACGT.

    Levanta ValueError se encontrar qualquer caractere fora de ACGT, citando o caractere.
    """
    s = "".join(seq.split()).upper()
    for char in s:
        if char not in "ACGT":
            raise ValueError(f"Caractere inválido '{char}' na sequência: {seq}")
    return s


def rc(seq: str) -> str:
    """Retorna o complemento reverso da sequência usando Bio.Seq."""
    return str(Seq(seq).reverse_complement())


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
        for campo in ("F3", "B3", "F2", "F1c", "B2", "B1c", "ligante"):
            val = getattr(self, campo)
            if val is not None:
                object.__setattr__(self, campo, normalizar(val))
        for campo in ("LF", "LB"):
            val = getattr(self, campo)
            if val is not None:
                object.__setattr__(self, campo, normalizar(val))

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
