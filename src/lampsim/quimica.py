"""Química da reação: pH ao longo do tempo e forma do indicador.

Modelo (LAMP-S02-ESP):
- A amplificação libera prótons e o pH cai seguindo uma curva logística,
  normalizada para partir de zero em t = 0:
      L(t) = 1 / (1 + exp(-k (t - t0)))
      progresso(t) = A · (L(t) - L(0)) / (1 - L(0))
      pH(t) = pH_inicial - queda_max · progresso(t) - deriva · t
  Sem a normalização, um tubo com t0 cedo já começaria mais ácido que um
  negativo no instante zero.
- O indicador existe em duas formas; pela equação de Henderson-Hasselbalch,
  a fração na forma ácida (amarela) é
      f_acida(pH) = 1 / (1 + 10^(pH - pKa))
"""

import numpy as np


def fracao_acida(ph, pka):
    """Fração do indicador na forma ácida (amarela), entre 0 e 1."""
    return 1.0 / (1.0 + 10.0 ** (np.asarray(ph, dtype=float) - pka))


def fracao_basica(ph, pka):
    """Fração do indicador na forma básica (rosa), entre 0 e 1."""
    return 1.0 - fracao_acida(ph, pka)


def _logistica(x):
    """1 / (1 + exp(-x)) sem estouro numérico para |x| grande."""
    return 0.5 * (1.0 + np.tanh(0.5 * np.asarray(x, dtype=float)))


def progresso_logistico(t, t0, k, A):
    """Fração da queda de pH já ocorrida no tempo t: 0 em t = 0, tende a A."""
    t = np.asarray(t, dtype=float)
    l0 = _logistica(-k * t0)
    return A * (_logistica(k * (t - t0)) - l0) / (1.0 - l0)


def ph_no_tempo(t, ph_inicial, queda_max, t0, k, A, deriva=0.0):
    """pH do mix ao longo do tempo (minutos)."""
    t = np.asarray(t, dtype=float)
    return ph_inicial - queda_max * progresso_logistico(t, t0, k, A) - deriva * t
