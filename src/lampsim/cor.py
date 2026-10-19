"""Da fração do indicador até a cor RGB.

Caminho físico: espectro de absorção da mistura -> transmitância (Beer-Lambert)
-> luz que chega à câmera -> cor.

Duas rotas de conversão para RGB, para a decisão de faseamento registrada na
LAMP-S02-ESP (Seção 2.4):
- "cie": colorimetria padrão. Integra o espectro com as funções de
  correspondência de cor CIE 1931 (aproximação analítica de Wyman, Sloan e
  Shirley, 2013), converte XYZ para sRGB.
- "3bandas": um sensor de câmera com três filtros gaussianos (vermelho, verde,
  azul). O RGB do sensor passa para XYZ pela matriz da câmera (mínimos
  quadrados entre as sensibilidades do sensor e as funções CIE) e daí para
  sRGB, como no processamento de uma câmera; a gama do sRGB só entra no fim.

Nas duas rotas a cor é normalizada pelo branco (tubo transparente = RGB 1,1,1),
ou seja, um balanço de branco ideal. As imperfeições de câmera entram só na
fase 2.
"""

import numpy as np

from .config import ESPECTRO

LAMBDA_NM = np.arange(380.0, 781.0, 5.0)
_DL = 5.0

_XYZ_PARA_SRGB = np.array([
    [3.2406, -1.5372, -0.4986],
    [-0.9689, 1.8758, 0.0415],
    [0.0557, -0.2040, 1.0570],
])


def _gauss_assim(lam, mu, s1, s2):
    s = np.where(lam < mu, s1, s2)
    return np.exp(-0.5 * ((lam - mu) / s) ** 2)


def cmf_cie1931(lam=LAMBDA_NM):
    """Funções de correspondência de cor CIE 1931 (x̄, ȳ, z̄), forma analítica."""
    x = (1.056 * _gauss_assim(lam, 599.8, 37.9, 31.0)
         + 0.362 * _gauss_assim(lam, 442.0, 16.0, 26.7)
         - 0.065 * _gauss_assim(lam, 501.1, 20.4, 26.2))
    y = (0.821 * _gauss_assim(lam, 568.8, 46.9, 40.5)
         + 0.286 * _gauss_assim(lam, 530.9, 16.3, 31.1))
    z = (1.217 * _gauss_assim(lam, 437.0, 11.8, 36.0)
         + 0.681 * _gauss_assim(lam, 459.0, 26.0, 13.8))
    return np.stack([x, y, z])


def sensibilidade_3bandas(lam=LAMBDA_NM):
    """Sensibilidade simplificada de um sensor RGB (linhas: R, G, B)."""
    r = np.exp(-0.5 * ((lam - 600.0) / 30.0) ** 2)
    g = np.exp(-0.5 * ((lam - 540.0) / 30.0) ** 2)
    b = np.exp(-0.5 * ((lam - 460.0) / 25.0) ** 2)
    return np.stack([r, g, b])


def matriz_camera_para_xyz(lam=LAMBDA_NM, temp_k=None):
    """Matriz 3x3 M tal que XYZ ≈ M @ RGB_do_sensor, por mínimos quadrados.

    O ajuste é feito sobre a luz do iluminante em cada comprimento de onda: procura M
    com M @ (sensibilidade · luz) o mais perto possível de (CIE · luz).
    """
    temp = ESPECTRO["temperatura_cor_iluminante_k"] if temp_k is None else temp_k
    luz = iluminante_planck(temp, lam)
    sens = sensibilidade_3bandas(lam) * luz  # (3, n)
    alvo = cmf_cie1931(lam) * luz            # (3, n)
    solucao, *_ = np.linalg.lstsq(sens.T, alvo.T, rcond=None)  # (3, 3): sens.T @ X ≈ alvo.T
    return solucao.T


def iluminante_planck(temp_k, lam=LAMBDA_NM):
    """Espectro de corpo negro (aproximação de luz branca), normalizado."""
    h, c, kb = 6.62607015e-34, 2.99792458e8, 1.380649e-23
    lm = lam * 1e-9
    e = 1.0 / (lm ** 5 * (np.exp(h * c / (lm * kb * temp_k)) - 1.0))
    return e / e.max()


def _bandas(bandas, lam):
    """Soma das bandas (pico_nm, sigma_azul_nm, sigma_vermelho_nm, amplitude)."""
    total = np.zeros_like(np.asarray(lam, dtype=float))
    for mu, s_azul, s_verm, a in bandas:
        total = total + a * _gauss_assim(lam, mu, s_azul, s_verm)
    return total


def absorbancia_pico_acida(espectro=ESPECTRO):
    """Absorbância da forma ácida pura no pico de 433 nm: ε · c · l (adimensional)."""
    return (espectro["epsilon_433_acida_M_cm"] * espectro["concentracao_uM"] * 1e-6
            * espectro["caminho_cm"])


def absorbancia(f_acida, espectro=ESPECTRO, lam=LAMBDA_NM, turbidez=None):
    """Absorbância da mistura para cada fração ácida. Saída: (n, n_lambda).

    turbidez: absorbância neutra extra de cada tubo (n,), por exemplo o precipitado
    de pirofosfato de magnésio; None = nenhuma.
    """
    f = np.atleast_1d(np.asarray(f_acida, dtype=float))[:, None]
    escala = absorbancia_pico_acida(espectro)
    a = escala * (f * _bandas(espectro["acida"], lam)
                  + (1.0 - f) * _bandas(espectro["basica"], lam)) + espectro["fundo"]
    if turbidez is not None:
        a = a + np.atleast_1d(np.asarray(turbidez, dtype=float))[:, None]
    return a


def _gama_srgb(lin):
    lin = np.clip(lin, 0.0, 1.0)
    return np.where(lin <= 0.0031308, 12.92 * lin, 1.055 * lin ** (1 / 2.4) - 0.055)


def transmitancia_para_rgb(transmitancia, rota="cie", espectro=ESPECTRO, lam=LAMBDA_NM):
    """Converte transmitância (n, n_lambda) em RGB com gama sRGB, valores 0-1."""
    T = np.atleast_2d(transmitancia)
    luz = iluminante_planck(espectro["temperatura_cor_iluminante_k"], lam)
    if rota == "cie":
        sens = cmf_cie1931(lam)
        xyz = (T * luz) @ sens.T * _DL
        xyz_branco = (luz @ sens.T) * _DL
    elif rota == "3bandas":
        sens = sensibilidade_3bandas(lam)
        m = matriz_camera_para_xyz(lam, espectro["temperatura_cor_iluminante_k"])
        xyz = ((T * luz) @ sens.T * _DL) @ m.T
        xyz_branco = ((luz @ sens.T) * _DL) @ m.T
    else:
        raise ValueError(f"rota desconhecida: {rota!r} (use 'cie' ou '3bandas')")
    lin = xyz @ _XYZ_PARA_SRGB.T
    lin_branco = xyz_branco @ _XYZ_PARA_SRGB.T
    return _gama_srgb(lin / lin_branco)


def cor_do_indicador(f_acida, rota="cie", espectro=ESPECTRO, turbidez=None):
    """Cor RGB (0-1) do tubo para cada fração ácida do indicador. Saída: (n, 3)."""
    T = 10.0 ** (-absorbancia(f_acida, espectro, turbidez=turbidez))
    return transmitancia_para_rgb(T, rota=rota, espectro=espectro)


def rgb_para_hex(rgb):
    """Lista de cores '#RRGGBB' (maiúsculas) a partir de RGB 0-1 com forma (n, 3)."""
    bytes_ = np.clip(np.round(np.atleast_2d(rgb) * 255.0), 0, 255).astype(int)
    return [f"#{r:02X}{g:02X}{b:02X}" for r, g, b in bytes_]
