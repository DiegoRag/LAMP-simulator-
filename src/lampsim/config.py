"""Parâmetros do simulador — versão v1 (fase 5 da correção).

Todos os números do simulador ficam neste arquivo, e só nele. Assim, conferir
com a especificação (LAMP-S02-ESP) ou recalibrar com os tempos do artigo é
editar um lugar só.

Situação dos valores:
- QUIMICA["pka"]: vermelho de fenol, pKa ≈ 8,0. Yao e Byrne (2001) dão
  pK° = 5,798 + 666,7/T, ou seja, 8,03 a 25 °C em força iônica zero. O
  indicador em si é PREMISSA: a composição do mix ainda não foi confirmada
  pela FMABC.
- ESPECTRO: as razões de absortividade das duas formas são de literatura
  (Yao e Byrne, 2001); larguras de banda são ajustadas a essas razões e aos
  pontos isosbésticos; concentração e caminho óptico são A CALIBRAR.
- CINETICA e faixas de CLASSES: PROVISÓRIAS (a calibrar com os tempos do
  artigo). Escolhidas para que um positivo forte mude de cor bem antes da
  leitura e um positivo perto do limite de detecção fique perto dela.
"""

VERSAO_PARAMETROS = "v1-fase5"

# Linha do tempo da simulação (minutos)
TEMPO = {
    "duracao_min": 30.0,
    "passo_min": 0.5,
    "leitura_nominal_min": 15.0,  # tempo de referência para o diagnóstico
}

QUIMICA = {
    "pka": 8.0,  # vermelho de fenol
}

# Espectro de absorção do vermelho de fenol.
# Cada banda: (pico_nm, sigma_azul_nm, sigma_vermelho_nm, amplitude), uma gaussiana
# com largura própria de cada lado do pico. A amplitude é relativa a ε433 da forma
# ácida (HI⁻ = 1). Vínculos de literatura (Yao e Byrne, 2001, Environ. Sci. Technol.
# 35:1197, λ1 = 433 nm e λ2 = 558 nm):
#   e1 = ε558(HI⁻)/ε433(HI⁻) = 0,00244 -> sigma 36,0 nm da banda ácida
#   e2 = ε558(I²⁻)/ε433(HI⁻) = 2,734   -> amplitude da banda básica
#   e3 = ε433(I²⁻)/ε433(HI⁻) = 0,1075  -> absorção da forma básica em 433 nm
# Pontos isosbésticos usados no ajuste: 478 nm (477 nm em Rovati et al., 2012;
# ~480 nm em meio de cultura) e 365 nm (Rovati et al., 2012). A banda de 340 nm da
# forma básica é o que o ajuste pede para cumprir e3 e o isosbéstico de 365 nm.
# O sigma vermelho de 14,5 nm da banda básica vem da meia-largura de ~17 nm do lado
# vermelho no espectro em pH 10 de Rovati et al. (2012, fig. 7, digitalizada no pacote
# colorSpec do R); é dye imobilizado em sol-gel, então fica A CALIBRAR em solução.
ESPECTRO = {
    "acida": [(433.0, 36.0, 36.0, 1.0)],
    "basica": [(558.0, 41.6, 14.5, 2.734), (340.0, 71.9, 71.9, 0.179)],
    # Absorbância no pico da forma ácida = ε433(HI⁻) · c · l. ε433 ≈ 23 000 M⁻¹ cm⁻¹
    # (E1% de 600 a 700 em pH 6,5 e 430–435 nm, especificação de fornecedor);
    # c = 100 µM e l = 0,5 cm são A CALIBRAR (composição do mix não confirmada).
    "epsilon_433_acida_M_cm": 23000.0,
    "concentracao_uM": 100.0,
    "caminho_cm": 0.5,
    "fundo": 0.02,  # absorbância neutra do tubo e do mix, constante
    # turbidez do pirofosfato de magnésio: absorbância neutra proporcional ao
    # progresso da amplificação (0 sem produto, este valor com reação completa); A CALIBRAR
    "turbidez_max": 0.05,
    "temperatura_cor_iluminante_k": 6500.0,
}

# Cinética: do número de cópias e da eficiência ao tempo da virada (PROVISÓRIO).
# t0 = atraso + (log10_copias_ref - log10(cópias)) · min_por_log / eficiência
# t0 é o meio da queda de pH (minutos); k = k · eficiência (1/min).
CINETICA = {
    "atraso_min": 2.0,          # aquecimento e início da amplificação
    "log10_copias_ref": 9.0,    # cópias em que t0 = atraso (eficiência 1)
    "min_por_log": 1.5,         # minutos a mais por década a menos de cópias
    "k_por_min": (0.8, 1.2),    # velocidade da virada com eficiência 1
    "A": (0.90, 1.00),          # fração da queda de pH que uma reação completa atinge
    "sd_t0_min": 0.3,           # variação de t0 entre tubos iguais (desvio-padrão)
    "t0_max_min": 999.0,        # teto de t0 (reação que não vira na janela)
}

# Mix da reação (sorteados por tubo, uniformes na faixa)
MIX = {
    "ph_inicial": (8.5, 8.9),
    "queda_max": (2.2, 2.6),    # queda de pH de uma reação completa
    "deriva": (0.0, 0.003),     # queda lenta sem amplificação (pH/min)
}

# Leitura do resultado no tempo nominal, a partir da curva simulada
LEITURA = {
    "limiar_positivo": 0.5,     # fração ácida na leitura a partir da qual o tubo é positivo
    "limiar_falha": 0.5,        # fração ácida em t = 0 a partir da qual a reação é inválida
}

# Papéis de um tubo: amostra ou controle
PAPEIS = ("amostra", "controle_negativo", "controle_positivo")

# Classes ópticas: o que o tubo faz, sem dizer se é amostra ou controle.
# log10_copias: faixa do log10 da média de cópias (None = nenhuma cópia);
# eficiencia: faixa da eficiência; o mix vem de MIX, salvo as chaves de mix da classe.
CLASSES = {
    "positivo_forte": {"log10_copias": (4.0, 6.0), "eficiencia": (0.90, 1.00)},
    "positivo_perto_lod": {"log10_copias": (0.5, 2.5), "eficiencia": (0.75, 1.00)},
    # amplificação tardia de poucas cópias: típica de contaminação de um controle
    "positivo_tardio": {"log10_copias": (0.0, 1.5), "eficiencia": (0.85, 1.00)},
    "negativo": {"log10_copias": None, "eficiencia": (0.90, 1.00)},
    # falha de reação: a amostra acidifica o mix já no início (laranja ou amarelo em t = 0)
    "falha_reacao": {"log10_copias": None, "eficiencia": (0.90, 1.00),
                     "ph_inicial": (6.4, 7.2), "queda_max": (0.2, 0.5)},
}

# Composição de um lote: (nome do grupo, classe óptica, papel). O diagnóstico não
# vem daqui: sai da curva simulada no tempo de leitura (classes.diagnosticar).
LOTE = (
    ("positivo_forte", "positivo_forte", "amostra"),
    ("positivo_perto_lod", "positivo_perto_lod", "amostra"),
    ("negativo", "negativo", "amostra"),
    ("controle_negativo", "negativo", "controle_negativo"),
    ("controle_contaminado", "positivo_tardio", "controle_negativo"),
    ("falha_reacao", "falha_reacao", "amostra"),
)
