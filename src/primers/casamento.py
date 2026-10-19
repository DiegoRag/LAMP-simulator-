"""Casamento com tolerância entre as regiões de um conjunto LAMP e o alvo, por distância de Hamming.

Diz onde uma região que não aparece exata mais se aproxima do alvo e em quais bases difere.
Limitação declarada: só substituições, sem inserção nem deleção (indels).

Fora do caso só ACGT (que usa str.find e a casa dos pombos), a comparação base a base é
vetorizada em numpy: cada base vira uma máscara de 4 bits (A=1, C=2, G=4, T=8; IUPAC é a
união) e duas bases casam quando as máscaras se cruzam, como em modelo.bases_compativeis.
"""

import re
from dataclasses import dataclass

import numpy as np

from primers.modelo import BASES_IUPAC, ConjuntoLAMP, Regiao, bases_compativeis, rc

_NAO_ACGT = re.compile(r"[^ACGT]")
REGIOES_RC: frozenset[str] = frozenset({"F1", "B2c", "B3c", "LF"})
ORDEM_REGIOES: tuple[str, ...] = ("F3", "F2", "LF", "F1", "B1c", "LB", "B2c", "B3c")
REGIAO_PARA_PRIMER: dict[str, str] = {
    "F3": "F3",
    "F2": "F2",
    "F1": "F1c",
    "B1c": "B1c",
    "B2c": "B2",
    "B3c": "B3",
    "LF": "LF",
    "LB": "LB",
}
# F1 e B1c: no haltere, a fita estende pela ponta que pareia com a ponta 5' de F1c (B1c)
REGIOES_PONTA_5: frozenset[str] = frozenset({"F1", "B1c"})

_BITS = {"A": 1, "C": 2, "G": 4, "T": 8}
_MASCARA = np.zeros(256, dtype=np.uint8)
for _letra, _bases in BASES_IUPAC.items():
    _MASCARA[ord(_letra)] = sum(_BITS[b] for b in _bases)


@dataclass(frozen=True)
class Diferenca:
    """Diferença de base entre uma sequência e o alvo em determinada posição."""

    posicao: int  # coordenada no alvo (na fita escolhida), base zero
    esperado: str  # base da região derivada do primer, na fita escolhida
    no_alvo: str  # base do alvo nessa posição


@dataclass(frozen=True)
class Posicao:
    """Posição de melhor casamento de uma região no alvo com suas diferenças."""

    regiao: Regiao  # primers.modelo.Regiao
    diferencas: tuple[Diferenca, ...]  # vazia = casamento exato
    fita: str = "+"  # "+" para alvo, "-" para rc(alvo)


@dataclass(frozen=True)
class ErroDeBase:
    """Uma diferença de base contada do lado do primer, com a distância até a ponta que estende."""

    regiao: str  # "F3", "F2", "F1", "B1c", "B2c", "B3c", "LF" ou "LB"
    primer: str  # primer que traz a região: "F1c" para F1, "B2" para B2c, "B3" para B3c
    fita: str  # "+" ou "-": fita do alvo em que a região foi localizada
    posicao: int  # base 1, na fita `fita` do alvo
    base_no_primer: int  # base 1, contada a partir da ponta 5' do primer
    comprimento: int  # nt do primer (igual ao da região)
    primer_tem: str  # base do primer nessa posição, lida 5'->3' no primer
    alvo_pede: str  # base que o primer precisaria ter ali para casar com o alvo
    ponta: str  # "3'" ou "5'": a ponta de onde sai a extensão
    distancia_ponta: int  # nt entre a base e essa ponta; 0 = a própria base terminal


def _so_acgt(texto: str) -> bool:
    """True se o texto só tem A, C, G e T (para no primeiro caractere diferente)."""
    return _NAO_ACGT.search(texto) is None


def _mascaras(texto: str) -> np.ndarray:
    """Máscara de 4 bits de cada base do texto (0 para caractere fora do IUPAC)."""
    return _MASCARA[np.frombuffer(texto.encode("ascii", "replace"), dtype=np.uint8)]


def _inicios_sem_erro(seq_m: np.ndarray, alvo_m: np.ndarray, n_inicios: int) -> np.ndarray:
    """Inícios i em [0, n_inicios) em que seq_m casa base a base com alvo_m[i:], sem erro.

    Filtra base a base: depois das primeiras bases sobra pouca posição para conferir.
    """
    idx = np.flatnonzero(alvo_m[:n_inicios] & seq_m[0])
    for j in range(1, len(seq_m)):
        if idx.size == 0:
            break
        idx = idx[(alvo_m[idx + j] & seq_m[j]) != 0]
    return idx


def _erros_nos_inicios(seq_m: np.ndarray, alvo_m: np.ndarray, inicios: np.ndarray) -> np.ndarray:
    """Número de incompatibilidades de seq_m contra alvo_m[i:i+k] em cada início i dado."""
    erros = np.zeros(inicios.size, dtype=np.int32)
    for j, m in enumerate(seq_m):
        erros += (alvo_m[inicios + j] & m) == 0
    return erros


def ocorrencias_compativeis(
    texto: str,
    sub: str,
    limite_ini: int | None = None,
) -> list[Regiao]:
    """Retorna todas as ocorrências em que sub casa com texto sem incompatibilidades IUPAC."""
    if not sub or not texto or len(sub) > len(texto):
        return []
    k = len(sub)
    lim = len(texto) if limite_ini is None else min(len(texto), limite_ini)
    max_start = min(lim, len(texto) - k + 1)
    if max_start <= 0:
        return []

    res: list[Regiao] = []
    if _so_acgt(sub) and _so_acgt(texto):
        pos = 0
        while pos < max_start:
            idx = texto.find(sub, pos)
            if idx == -1 or idx >= max_start:
                break
            res.append(Regiao(idx, idx + k))
            pos = idx + 1
        return res

    inicios = _inicios_sem_erro(_mascaras(sub), _mascaras(texto[: max_start + k - 1]), max_start)
    return [Regiao(int(i), int(i) + k) for i in inicios]


def regioes_do_conjunto(conjunto: ConjuntoLAMP) -> list[tuple[str, str]]:
    """Devolve as regiões na fita + na convenção de restricoes.localizar."""
    mapa_regioes: list[tuple[str, str]] = [
        ("F3", conjunto.F3),
        ("F2", conjunto.F2),
        ("F1", rc(conjunto.F1c)),
        ("B1c", conjunto.B1c),
        ("B2c", rc(conjunto.B2)),
        ("B3c", rc(conjunto.B3)),
    ]
    if conjunto.LF is not None:
        mapa_regioes.append(("LF", rc(conjunto.LF)))
    if conjunto.LB is not None:
        mapa_regioes.append(("LB", conjunto.LB))
    return mapa_regioes


def ponta_extensao(reg_nome: str) -> str:
    """Ponta do primer de onde sai a extensão: "5'" para F1 e B1c (haltere), "3'" para as demais."""
    return "5'" if reg_nome in REGIOES_PONTA_5 else "3'"


def distancia_3p(reg_nome: str, regiao: Regiao, posicao: int) -> int:
    """Distância em nt (base 0) entre `posicao` e a base 3' terminal do primer correspondente.

    - Para reg_nome em REGIOES_RC ("F1", "B2c", "B3c", "LF"): a sequência na fita ativa é rc(primer),
      logo a extremidade 3' do primer está em regiao.ini -> retorna `posicao - regiao.ini`.
    - Para reg_nome em ("F3", "F2", "B1c", "LB"): a sequência na fita ativa é o próprio primer 5'->3',
      logo a extremidade 3' do primer está em regiao.fim - 1 -> retorna `(regiao.fim - 1) - posicao`.
    """
    if reg_nome in REGIOES_RC:
        return posicao - regiao.ini
    return (regiao.fim - 1) - posicao


def distancia_ponta_extensao(reg_nome: str, regiao: Regiao, posicao: int) -> int:
    """Distância em nt (base 0) entre `posicao` e a ponta que estende a partir da região.

    F1 = rc(F1c) na fita escolhida: a ponta 5' de F1c fica em fim - 1. B1c está como o primer:
    a ponta 5' fica em ini. As demais regiões estendem pela ponta 3' do primer (distancia_3p).
    """
    if reg_nome == "F1":
        return (regiao.fim - 1) - posicao
    if reg_nome == "B1c":
        return posicao - regiao.ini
    return distancia_3p(reg_nome, regiao, posicao)


def erros_de_base(reg_nome: str, pos: Posicao) -> list[ErroDeBase]:
    """As diferenças de uma posição aproximada, contadas no primer e situadas em relação à ponta."""
    k = pos.regiao.fim - pos.regiao.ini
    erros: list[ErroDeBase] = []
    for d in pos.diferencas:
        deslocamento = d.posicao - pos.regiao.ini  # base 0 dentro da região, na fita escolhida
        if reg_nome in REGIOES_RC:  # o primer é rc(região): conta do outro lado e complementa
            base_no_primer = k - deslocamento
            primer_tem, alvo_pede = rc(d.esperado), rc(d.no_alvo)
        else:
            base_no_primer = deslocamento + 1
            primer_tem, alvo_pede = d.esperado, d.no_alvo
        erros.append(
            ErroDeBase(
                regiao=reg_nome,
                primer=REGIAO_PARA_PRIMER.get(reg_nome, reg_nome),
                fita=pos.fita,
                posicao=d.posicao + 1,
                base_no_primer=base_no_primer,
                comprimento=k,
                primer_tem=primer_tem,
                alvo_pede=alvo_pede,
                ponta=ponta_extensao(reg_nome),
                distancia_ponta=distancia_ponta_extensao(reg_nome, pos.regiao, d.posicao),
            )
        )
    return erros


def janela_esperada(
    reg_nome: str,
    regioes_ancoradas: dict[str, Regiao],
    tam_alvo: int,
    cfg: dict | None = None,
) -> Regiao | None:
    """Calcula a janela geométrica [ini, fim) coerente para procurar `reg_nome` dado `regioes_ancoradas`."""
    if tam_alvo <= 0:
        return None
    if not regioes_ancoradas:
        return Regiao(0, tam_alvo)

    if reg_nome == "LF" and "F2" in regioes_ancoradas and "F1" in regioes_ancoradas:
        ini_w = max(0, regioes_ancoradas["F2"].fim)
        fim_w = min(tam_alvo, regioes_ancoradas["F1"].ini)
        return Regiao(ini_w, fim_w) if fim_w > ini_w else None

    if reg_nome == "LB" and "B1c" in regioes_ancoradas and "B2c" in regioes_ancoradas:
        ini_w = max(0, regioes_ancoradas["B1c"].fim)
        fim_w = min(tam_alvo, regioes_ancoradas["B2c"].ini)
        return Regiao(ini_w, fim_w) if fim_w > ini_w else None

    if reg_nome not in ORDEM_REGIOES:
        return Regiao(0, tam_alvo)

    idx = ORDEM_REGIOES.index(reg_nome)
    esq = [r for k, r in regioes_ancoradas.items() if k in ORDEM_REGIOES and ORDEM_REGIOES.index(k) < idx]
    dir_ = [r for k, r in regioes_ancoradas.items() if k in ORDEM_REGIOES and ORDEM_REGIOES.index(k) > idx]

    f3_b3_max = int((cfg or {}).get("distancias", {}).get("F3_B3_total", [120, 280])[1])
    margem = f3_b3_max + 80
    ini_w = max(r.ini for r in esq) if esq else max(0, min(r.ini for r in dir_) - margem)
    fim_w = min(r.fim for r in dir_) if dir_ else min(tam_alvo, max(r.fim for r in esq) + margem)
    ini_w = max(0, min(ini_w, tam_alvo))
    fim_w = max(0, min(fim_w, tam_alvo))
    if fim_w <= ini_w:
        return None
    return Regiao(ini_w, fim_w)


def _buscar_em_fita(
    seq: str,
    texto: str,
    fita_id: str,
    janela_fita: Regiao | None,
    max_erros: int | None,
) -> Posicao | None:
    k = len(seq)
    if k == 0 or k > len(texto):
        return None

    ini_b = 0 if janela_fita is None else max(0, janela_fita.ini)
    fim_b = len(texto) if janela_fita is None else min(len(texto), janela_fita.fim)

    if fim_b - ini_b < k:
        return None

    # Busca exata rápida em C
    pos_exato = texto.find(seq, ini_b, fim_b)
    if pos_exato != -1:
        return Posicao(regiao=Regiao(pos_exato, pos_exato + k), diferencas=(), fita=fita_id)

    # Se puramente ACGT e max_erros <= 0, como a busca exata falhou, não há casamento com 0 erros.
    if max_erros is not None and max_erros <= 0 and _so_acgt(seq) and _so_acgt(texto):
        return None

    # Filtro de candidatos por casa dos pombos (para sequências longas puramente ACGT)
    if (
        max_erros is not None
        and (fim_b - ini_b) > 2000
        and k >= (max_erros + 1) * 3
        and _so_acgt(seq)
        and _so_acgt(texto)
    ):
        m = max_erros
        set_i: set[int] = set()
        for r in range(m + 1):
            o_r = (r * k) // (m + 1)
            o_r1 = ((r + 1) * k) // (m + 1)
            bloco = seq[o_r:o_r1]
            lim_j = fim_b - k + o_r1
            curr = ini_b + o_r
            while True:
                j = texto.find(bloco, curr, lim_j)
                if j == -1:
                    break
                i = j - o_r
                if ini_b <= i <= fim_b - k:
                    set_i.add(i)
                curr = j + 1
        candidatos_i = sorted(set_i)
    else:
        # Vetorizado (numpy), com o mesmo resultado do laço abaixo: a primeira posição com o
        # menor número de erros, se não passar do teto. Com teto m, uma janela com até m erros
        # tem um dos m + 1 blocos da sequência sem erro (casa dos pombos): só esses inícios
        # são conferidos.
        teto = max(0, max_erros if max_erros is not None else k)
        alvo_m = _mascaras(texto[ini_b:fim_b])
        seq_m = _mascaras(seq)
        n = fim_b - ini_b - k + 1
        if teto < k:
            blocos = [((r * k) // (teto + 1), ((r + 1) * k) // (teto + 1)) for r in range(teto + 1)]
            inicios = np.unique(
                np.concatenate([_inicios_sem_erro(seq_m[a:b], alvo_m[a:], n) for a, b in blocos])
            )
        else:
            inicios = np.arange(n)
        if inicios.size == 0:
            return None
        erros = _erros_nos_inicios(seq_m, alvo_m, inicios)
        melhor = int(np.argmin(erros))
        if int(erros[melhor]) > teto:
            return None
        i = ini_b + int(inicios[melhor])
        difs = tuple(
            Diferenca(posicao=i + idx, esperado=seq[idx], no_alvo=texto[i + idx])
            for idx in range(k)
            if not bases_compativeis(seq[idx], texto[i + idx])
        )
        return Posicao(regiao=Regiao(i, i + k), diferencas=difs, fita=fita_id)

    menor_dif: int | None = None
    melhor_pos: Posicao | None = None

    for i in candidatos_i:
        teto = max_erros if max_erros is not None else k
        if menor_dif is not None:
            teto = min(teto, menor_dif - 1)

        erros = 0
        ultrapassou = False
        for idx in range(k):
            if not bases_compativeis(seq[idx], texto[i + idx]):
                erros += 1
                if erros > teto:
                    ultrapassou = True
                    break

        if not ultrapassou:
            difs = tuple(
                Diferenca(posicao=i + idx, esperado=seq[idx], no_alvo=texto[i + idx])
                for idx in range(k)
                if not bases_compativeis(seq[idx], texto[i + idx])
            )
            menor_dif = len(difs)
            melhor_pos = Posicao(regiao=Regiao(i, i + k), diferencas=difs, fita=fita_id)
            if menor_dif == 0:
                break

    return melhor_pos


def melhor_posicao(
    seq_fita_mais: str,
    alvo: str,
    *,
    max_erros: int | None = 3,
    janela: Regiao | None = None,
    ambas_fitas: bool = False,
) -> Posicao | None:
    """Encontra a melhor janela de casamento por distância de Hamming até `max_erros`."""
    if not ambas_fitas:
        return _buscar_em_fita(seq_fita_mais, alvo, "+", janela, max_erros)

    p_mais = _buscar_em_fita(seq_fita_mais, alvo, "+", janela, max_erros)
    p_menos = _buscar_em_fita(seq_fita_mais, rc(alvo), "-", janela, max_erros)

    if p_mais is None and p_menos is None:
        return None
    if p_mais is None:
        return p_menos
    if p_menos is None:
        return p_mais
    if len(p_menos.diferencas) < len(p_mais.diferencas):
        return p_menos
    return p_mais
