"""Amostras simuladas: sorteio de parâmetros, simulação e diagnóstico pela curva.

Uma amostra nasce de duas entradas físicas, a eficiência do conjunto de primers
contra o alvo (0 a 1) e a média de cópias de molde no tubo, mais o mix sorteado.
As classes ópticas de config.CLASSES são só faixas dessas entradas. O papel do tubo
(amostra, controle negativo, controle positivo) é um campo à parte, e o diagnóstico
sai da curva simulada no tempo de leitura, nunca do nome da classe.
"""

from dataclasses import asdict, dataclass

import numpy as np

from .cinetica import parametros_cineticos, sortear_copias
from .config import CINETICA, CLASSES, ESPECTRO, LEITURA, MIX, PAPEIS, QUIMICA, TEMPO
from .cor import cor_do_indicador
from .quimica import fracao_acida, ph_no_tempo, progresso_logistico

NOMES_CLASSES = tuple(CLASSES)
DIAGNOSTICOS = ("positivo", "negativo", "controle_ok", "invalido")


@dataclass(frozen=True)
class Amostra:
    id: str
    classe: str  # classe óptica de CLASSES, ou "reacao" se veio de eficiência e cópias
    papel: str  # amostra, controle_negativo ou controle_positivo
    diagnostico: str  # derivado da curva no tempo de leitura: positivo, negativo, controle_ok, invalido
    semente: int
    copias: int  # cópias de molde no tubo (Poisson em volta da média)
    eficiencia: float
    ph_inicial: float
    queda_max: float
    A: float
    t0: float
    k: float
    deriva: float

    def como_dict(self):
        return asdict(self)


def eixo_tempo(tempo=TEMPO):
    """Minutos amostrados, de 0 até a duração, inclusive."""
    n = int(round(tempo["duracao_min"] / tempo["passo_min"])) + 1
    return np.linspace(0.0, tempo["duracao_min"], n)


def diagnosticar(fracao_inicial, fracao_leitura, papel, leitura=LEITURA):
    """Diagnóstico de um tubo a partir da fração ácida em t = 0 e no tempo de leitura."""
    if papel not in PAPEIS:
        raise ValueError(f"papel desconhecido: {papel!r}; opções: {list(PAPEIS)}")
    if fracao_inicial >= leitura["limiar_falha"]:
        return "invalido"
    amplificou = fracao_leitura >= leitura["limiar_positivo"]
    if papel == "amostra":
        return "positivo" if amplificou else "negativo"
    if papel == "controle_negativo":
        return "invalido" if amplificou else "controle_ok"
    return "controle_ok" if amplificou else "invalido"


def diagnostico_da_curva(t, fracao, papel, leitura_min=None, leitura=LEITURA):
    """Diagnóstico a partir de uma curva (t em minutos, fração ácida), lida em leitura_min."""
    t = np.asarray(t, dtype=float)
    fracao = np.asarray(fracao, dtype=float)
    t_l = TEMPO["leitura_nominal_min"] if leitura_min is None else float(leitura_min)
    return diagnosticar(float(fracao[0]), float(np.interp(t_l, t, fracao)), papel, leitura)


def _montar(rng, id_amostra, classe, papel, semente, eficiencia, media_copias, faixas_mix, cinetica):
    if papel not in PAPEIS:
        raise ValueError(f"papel desconhecido: {papel!r}; opções: {list(PAPEIS)}")
    ph_inicial = float(rng.uniform(*faixas_mix["ph_inicial"]))
    queda_max = float(rng.uniform(*faixas_mix["queda_max"]))
    deriva = float(rng.uniform(*faixas_mix["deriva"]))
    copias = sortear_copias(media_copias, rng)
    cin = parametros_cineticos(eficiencia, copias, rng, cinetica)
    t_l = np.array([0.0, TEMPO["leitura_nominal_min"]])
    ph = ph_no_tempo(t_l, ph_inicial, queda_max, cin["t0"], cin["k"], cin["A"], deriva)
    f = fracao_acida(ph, QUIMICA["pka"])
    return Amostra(
        id=id_amostra,
        classe=classe,
        papel=papel,
        diagnostico=diagnosticar(float(f[0]), float(f[1]), papel),
        semente=int(semente),
        copias=int(copias),
        eficiencia=float(eficiencia),
        ph_inicial=ph_inicial,
        queda_max=queda_max,
        A=cin["A"],
        t0=cin["t0"],
        k=cin["k"],
        deriva=deriva,
    )


def amostra_de_reacao(eficiencia, copias_media, semente, papel="amostra", id_amostra=None,
                      mix=MIX, cinetica=CINETICA):
    """Amostra de um tubo com eficiência (0 a 1) e média de cópias dadas, reprodutível pela semente."""
    rng = np.random.default_rng(semente)
    return _montar(rng, id_amostra or f"reacao-{semente}", "reacao", papel, semente,
                   float(eficiencia), copias_media, mix, cinetica)


def sortear_amostra(classe, semente, id_amostra=None, classes=CLASSES, papel="amostra"):
    """Sorteia uma amostra da classe óptica, de forma reprodutível."""
    if classe not in classes:
        raise ValueError(f"classe desconhecida: {classe!r}; opções: {list(classes)}")
    faixas = classes[classe]
    rng = np.random.default_rng(semente)
    eficiencia = float(rng.uniform(*faixas["eficiencia"]))
    if faixas["log10_copias"] is None:
        media = 0.0
    else:
        media = float(10.0 ** rng.uniform(*faixas["log10_copias"]))
    faixas_mix = {chave: faixas.get(chave, MIX[chave]) for chave in MIX}
    return _montar(rng, id_amostra or f"{classe}-{semente}", classe, papel, semente,
                   eficiencia, media, faixas_mix, CINETICA)


def simular(amostra, t=None, rota="cie", pka=QUIMICA["pka"]):
    """Trajetória da amostra: tempo, pH, fração ácida, progresso e cor RGB em cada instante."""
    t = eixo_tempo() if t is None else np.asarray(t, dtype=float)
    progresso = progresso_logistico(t, amostra.t0, amostra.k, amostra.A)
    ph = ph_no_tempo(t, amostra.ph_inicial, amostra.queda_max, amostra.t0,
                     amostra.k, amostra.A, amostra.deriva)
    f = fracao_acida(ph, pka)
    turbidez = ESPECTRO["turbidez_max"] * progresso
    return {"t": t, "ph": ph, "fracao_acida": f, "progresso": progresso,
            "rgb": cor_do_indicador(f, rota=rota, turbidez=turbidez)}
