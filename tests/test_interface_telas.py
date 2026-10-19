"""Testes para as funções geradoras de contexto das telas da interface web."""

from pathlib import Path

from esteira.ponte import carregar_config as carregar_config_esteira
from interface.sessao import Sessao
from interface.telas import (
    AVISO_CURTO,
    tela_analise,
    tela_base,
    tela_desenho,
    tela_inicio,
    tela_parametros,
    tela_tubos,
)
from primers.entrada import carregar_config

RAIZ = Path(__file__).resolve().parents[1]
cfg = carregar_config(RAIZ / "config" / "primers.toml")
CFG_ESTEIRA = carregar_config_esteira(RAIZ / "config" / "esteira.toml")
cfg_web = {
    "entrada": {"tamanho_max_kb": 64, "alvo_desenho_max_pb": 2000},
    "analise": {"tm_como_calibracao": True},
    "heterodimeros": {"escala_kcal": [-10.0, -8.0, -6.0, -4.0]},
    "exemplos": [
        {
            "nome": "sintético",
            "conjunto": "tests/dados/conjunto_ok.csv",
            "alvo": "tests/dados/alvo_ok.fasta",
        },
        {
            "nome": "inexistente",
            "conjunto": "dados/nao_existe.csv",
            "alvo": "dados/nao_existe.fasta",
        },
    ],
}
ESCALA_TESTE = [-10.0, -8.0, -6.0, -4.0]

# conjunto_ok reprova em Tm e ponta na convenção do PrimerExplorer V5 (fase 3): itens comuns às telas
PONTA_F3 = {"item": "F3", "lado": "F", "texto": "ΔG da ponta 3′ −2,56 kcal/mol", "faixa": "máximo −4"}
PONTA_F1C = {"item": "F1c", "lado": "F", "texto": "ΔG da ponta 5′ −3,50 kcal/mol", "faixa": "máximo −4"}
TM_CONJUNTO_OK = [
    {"item": "F3", "lado": "F", "texto": "54,94 °C", "faixa": "faixa 59–61"},
    {"item": "B3", "lado": "B", "texto": "56,61 °C", "faixa": "faixa 59–61"},
    {"item": "F2", "lado": "F", "texto": "51,27 °C", "faixa": "faixa 59–61"},
    {"item": "F1c", "lado": "F", "texto": "58,90 °C", "faixa": "faixa 64–66"},
    {"item": "B2", "lado": "B", "texto": "52,89 °C", "faixa": "faixa 59–61"},
    {"item": "B1c", "lado": "B", "texto": "57,17 °C", "faixa": "faixa 64–66"},
]


def test_tela_base():
    sessao_vazia = Sessao(RAIZ, cfg, cfg_web)
    ctx_base_vazia = tela_base("inicio", sessao_vazia, ("feat/esteira", "abc1234"))
    assert ctx_base_vazia["projeto"] == "Nenhum arquivo carregado"
    assert ctx_base_vazia["aviso"] == AVISO_CURTO

    sessao_ex = Sessao(RAIZ, cfg, cfg_web)
    sessao_ex.carregar_exemplo("sintético")
    ctx_base_ex = tela_base("inicio", sessao_ex, ("feat/esteira", "abc1234"))
    assert ctx_base_ex["projeto"] == "conjunto_ok × alvo_ok"


def test_tela_inicio_e_sinonimo():
    sessao = Sessao(RAIZ, cfg, cfg_web)
    sessao.carregar_exemplo("sintético")
    sessao.analisar()

    dados_inicio = tela_inicio(sessao)
    conjunto = dados_inicio["conjunto"]
    assert conjunto is not None
    assert conjunto["resumo"] == "CSV · 6 primers · FIP e BIP montados a partir das partes"
    assert conjunto["reconhecidos"] == [
        {"texto": "F3", "lado": "F"},
        {"texto": "B3", "lado": "B"},
        {"texto": "F2", "lado": "F"},
        {"texto": "F1c", "lado": "F"},
        {"texto": "B2", "lado": "B"},
        {"texto": "B1c", "lado": "B"},
    ]
    alvo = dados_inicio["alvo"]
    assert alvo == {
        "nome_arquivo": "alvo_ok.fasta",
        "resumo": "FASTA · 1 registro · 216 pb",
        "cabecalho": "alvo_sintetico_ok",
    }
    assert dados_inicio["pode_analisar"] is True
    assert dados_inicio["pode_desenhar"] is True
    assert dados_inicio["motivo_sem_desenho"] == ""

    hist0 = {k: v for k, v in dados_inicio["historico"][0].items() if k != "data"}
    assert hist0 == {
        "tipo": "Análise",
        "conjunto": "conjunto_ok",
        "alvo": "alvo_ok",
        "resultado": "2 violações e 6 de Tm a calibrar · 6 de 6 regiões",  # fase 6
        "resultado_tipo": "viol",
        "href": "/analise/1",
    }

    # Teste de sinônimo
    sessao_sin = Sessao(RAIZ, cfg, cfg_web)
    conteudo_sin = "nome,sequencia\nF3,GTGACAGACTCTCCTTAGATA\nFLP,CGTATCCAGCGACAGGCA\n".encode()
    sessao_sin.carregar_conjunto("sin.csv", conteudo_sin)
    dados_sin = tela_inicio(sessao_sin)
    conj_sin = dados_sin["conjunto"]
    assert conj_sin is not None
    assert conj_sin["reconhecidos"] == [
        {"texto": "F3", "lado": "F"},
        {"texto": "LF ← FLP", "lado": "L"},
    ]
    assert conj_sin["resumo"] == "CSV · 2 primers"
    assert conj_sin["avisos"] == [
        "Conjunto incompleto: faltam primers essenciais (F3, B3, F2, F1c, B2, B1c)."
    ]


def test_tela_analise_exemplo():
    sessao = Sessao(RAIZ, cfg, cfg_web)
    sessao.carregar_exemplo("sintético")
    ex = sessao.analisar()

    dados = tela_analise(ex, cfg, tm_como_calibracao=True, escala_kcal=ESCALA_TESTE)
    assert dados["titulo"] == "conjunto_ok × alvo_ok"
    assert dados["subtitulo"] == (
        "8 primers · alvo alvo_ok, 216 pb · K⁺ 50 mM, NH₄⁺ 20 mM, Mg²⁺ 8 mM, dNTP 1,4 mM cada, 65 °C"
    )
    assert dados["kpis"] == [
        {"rotulo": "Primers lidos", "valor": "8", "detalhe": "", "tipo": "neutro"},
        # fase 6: com "a calibrar", Tm fora da faixa sai da contagem de violações
        {"rotulo": "Violações", "valor": "2", "detalhe": "e 6 de Tm a calibrar", "tipo": "viol"},
        {"rotulo": "Regiões localizadas", "valor": "6 de 6", "detalhe": "", "tipo": "neutro"},
        {"rotulo": "Distâncias", "valor": "7 de 7", "detalhe": "", "tipo": "neutro"},
    ]
    assert dados["aviso_calibracao"] == (
        "As faixas de Tm ainda não estão calibradas. 6 primers ficam fora da faixa de Tm, "
        "a até 7,7 °C dela. Leia as violações de Tm como comparação, não como reprovação."
    )

    mapa = dados["mapa"]
    assert mapa is not None
    assert mapa["segmentos"][1] == {
        "nome": "F2",
        "regiao": "F2",
        "ini": 31,
        "fim": 49,
        "esquerda": "14.352",
        "largura": "8.333",
        "fantasma": False,
    }
    assert [t["valor"] for t in mapa["ticks"]] == ["1", "40", "80", "120", "160", "200", "216"]
    assert mapa["ticks"][1]["esquerda"] == "18.519"
    assert mapa["ticks"][0]["alinhamento"] == "inicio"
    assert mapa["ticks"][-1]["alinhamento"] == "fim"
    assert mapa["faixas_fora"] == []

    assert dados["distancias"][5] == {
        "nome": "F2 → B2",
        "medida": "146",
        "faixa": "120–180",
        "dentro": True,
        "texto": "dentro",
        # fase 6: régua de -11,25 a 236,25 nt (0 a 1,25 · 180, com 5 % de folga de cada lado)
        "regua": {
            "banda_esquerda": "53.030",
            "banda_largura": "24.242",
            "marca": "63.535",
            "zero": "",
        },
    }

    primers = dados["primers"]
    assert primers[0] == {
        "nome": "F3",
        "lado": "F",
        "seq": "GTGACAGACTCTCCTTAGATA",
        "seq2": "",
        "regiao": "1–21",
        "regiao_viol": False,
        "nt": "21",
        "nt_viol": False,
        "gc": "42,86",
        "gc_viol": False,
        "tm": "54,94",
        "tm_estado": "calibrar",
        "tm_reacao": "58,75",
        "dg": "−2,56",
        "ponta": "3′",
        "grampo": "—",
        "homodimero": "−0,61",
        "situacao": "1 + Tm",
        "situacao_tipo": "viol",
    }
    assert primers[3]["grampo"] == "—"
    assert primers[6] == {
        "nome": "FIP",
        "lado": "X",
        "seq": "GTAAACCCCCAAGTGTTATTGG",
        "seq2": "CGCATATGGATTCTGTAC",
        "regiao": "—",
        "regiao_viol": False,
        "nt": "40",
        "nt_viol": False,
        "gc": "45,00",
        "gc_viol": False,
        "tm": "—",
        "tm_estado": "ok",
        "tm_reacao": "—",
        "dg": "—",
        "ponta": "",
        "grampo": "—",
        "homodimero": "−2,93",
        "situacao": "só métricas",
        "situacao_tipo": "neutro",
    }

    assert dados["grupos_violacoes"] == [
        {"titulo": "Composição e alvo · 2", "tipo": "viol", "itens": [PONTA_F3, PONTA_F1C]},
        {"titulo": "Tm — faixa a calibrar · 6", "tipo": "calibrar", "itens": TM_CONJUNTO_OK},
    ]
    assert dados["total_violacoes"] == 8

    het = dados["heterodimeros"]
    assert het["nomes"] == ["F3", "B3", "FIP", "BIP"]
    assert het["linhas"][0]["celulas"][0] == {"valor": "1,37", "faixa": 4, "diagonal": True}
    assert het["linhas"][1]["celulas"][3] == {"valor": "1,53", "faixa": 4, "diagonal": False}
    assert het["resumo"] == "Par mais estável: FIP × BIP, −1,03."
    assert het["legenda"] == [
        {"texto": "> −4", "faixa": 4},
        {"texto": "−6 a −4", "faixa": 3},
        {"texto": "−8 a −6", "faixa": 2},
        {"texto": "−10 a −8", "faixa": 1},
        {"texto": "≤ −10", "faixa": 0},
    ]

    assert dados["nota_fip_bip"] == (
        "FIP e BIP são montados das partes, sem ligante. "
        "Tm e ΔG de ponta não se aplicam a um primer quimérico."
    )
    assert dados["href_alternar_calibracao"] == "/analise/1?tm_calibracao=0"


def test_tela_analise_reguas():
    sessao = Sessao(RAIZ, cfg, cfg_web)
    sessao.carregar_exemplo("sintético")
    ex = sessao.analisar()

    dados = tela_analise(ex, cfg, tm_como_calibracao=True, escala_kcal=ESCALA_TESTE)
    assert len(dados["distancias"]) == 7

    for item in dados["distancias"]:
        regua = item["regua"]
        marca = float(regua["marca"])
        banda_esquerda = float(regua["banda_esquerda"])
        banda_largura = float(regua["banda_largura"])

        assert 0.0 <= marca <= 100.0
        assert 0.0 <= banda_esquerda
        assert banda_esquerda + banda_largura <= 100.001
        if item["dentro"]:
            assert banda_esquerda - 0.001 <= marca <= banda_esquerda + banda_largura + 0.001


def test_tela_analise_conjunto_lf():
    sessao = Sessao(RAIZ, cfg, cfg_web)
    conteudo_ok = (RAIZ / "tests/dados/conjunto_ok.csv").read_bytes()
    conteudo_lf = conteudo_ok + b"FLP,GTGCTGACGTGTCTGTTA\n"
    sessao.carregar_conjunto("conjunto_lf.csv", conteudo_lf)

    alvo_bytes = (RAIZ / "tests/dados/alvo_ok.fasta").read_bytes()
    sessao.carregar_alvo("alvo_ok.fasta", alvo_bytes)

    ex = sessao.analisar()
    dados = tela_analise(ex, cfg, tm_como_calibracao=True, escala_kcal=ESCALA_TESTE)

    assert dados["kpis"][1] == {
        "rotulo": "Violações",
        "valor": "4",
        "detalhe": "e 7 de Tm a calibrar",
        "tipo": "viol",
    }
    assert dados["kpis"][2] == {
        "rotulo": "Regiões localizadas",
        "valor": "6 de 7",
        "detalhe": "LF com 1 dif.",
        "tipo": "neutro",
    }

    assert dados["mapa"]["segmentos"][2] == {
        "nome": "LF",
        "regiao": "LF",
        "ini": 55,
        "fim": 73,
        "esquerda": "25.463",
        "largura": "8.333",
        "fantasma": True,
    }

    primer_lf = dados["primers"][8]
    assert primer_lf["regiao"] == "não localizado"
    assert primer_lf["regiao_viol"] is True
    assert primer_lf["tm"] == "55,98"
    assert primer_lf["tm_estado"] == "calibrar"
    assert primer_lf["situacao"] == "2 + Tm"
    assert primer_lf["lado"] == "L"

    ponta_lf = {"item": "LF", "lado": "L", "texto": "ΔG da ponta 3′ −3,74 kcal/mol", "faixa": "máximo −4"}
    nao_localizado_lf = {
        "item": "LF",
        "lado": "L",
        # fase 6: LF é lido pelo complemento reverso: a base 61 da fita + é a 13 do primer
        "texto": (
            "Não localizado: melhor posição 56–73 (fita +), 1 diferença "
            "(base 13 de 18 do LF: C no primer, o alvo pede A, a 5 nt da ponta 3′; posição 61)"
        ),
        "faixa": "",
    }
    tm_lf = {"item": "LF", "lado": "L", "texto": "55,98 °C", "faixa": "faixa 64–66"}
    assert dados["grupos_violacoes"] == [
        {
            "titulo": "Composição e alvo · 4",
            "tipo": "viol",
            "itens": [PONTA_F3, PONTA_F1C, ponta_lf, nao_localizado_lf],
        },
        {"titulo": "Tm — faixa a calibrar · 7", "tipo": "calibrar", "itens": [*TM_CONJUNTO_OK, tm_lf]},
    ]

    assert dados["aviso_calibracao"] == (
        "As faixas de Tm ainda não estão calibradas. 7 primers ficam fora da faixa de Tm, "
        "a até 8,0 °C dela. Leia as violações de Tm como comparação, não como reprovação."
    )


def test_tela_analise_conjunto_violacoes_sem_alvo():
    sessao = Sessao(RAIZ, cfg, cfg_web)
    conteudo_viol = (RAIZ / "tests/dados/conjunto_violacoes.csv").read_bytes()
    sessao.carregar_conjunto("conjunto_violacoes.csv", conteudo_viol)

    ex = sessao.analisar()

    # Com cal = True
    dados_cal = tela_analise(ex, cfg, tm_como_calibracao=True, escala_kcal=ESCALA_TESTE)
    assert dados_cal["titulo"] == "conjunto_violacoes"
    assert dados_cal["subtitulo"] == (
        "8 primers · sem alvo: regiões e distâncias não calculadas · "
        "K⁺ 50 mM, NH₄⁺ 20 mM, Mg²⁺ 8 mM, dNTP 1,4 mM cada, 65 °C"
    )
    assert dados_cal["kpis"][1] == {
        "rotulo": "Violações",
        "valor": "5",
        "detalhe": "e 6 de Tm a calibrar",
        "tipo": "viol",
    }
    assert dados_cal["kpis"][2] == {
        "rotulo": "Regiões localizadas",
        "valor": "—",
        "detalhe": "sem alvo",
        "tipo": "neutro",
    }
    assert dados_cal["mapa"] is None
    assert dados_cal["distancias"] == []
    assert dados_cal["aviso_calibracao"] == (
        "As faixas de Tm ainda não estão calibradas. 6 primers ficam fora da faixa de Tm, "
        "a até 24,0 °C dela. Leia as violações de Tm como comparação, não como reprovação."
    )
    assert dados_cal["grupos_violacoes"] == [
        {
            "titulo": "Composição e alvo · 5",
            "tipo": "viol",
            "itens": [
                {"item": "F3", "lado": "F", "texto": "GC 5,56 %", "faixa": "faixa 40–65"},
                {"item": "F3", "lado": "F", "texto": "ΔG da ponta 3′ −2,63 kcal/mol", "faixa": "máximo −4"},
                PONTA_F1C,
                {"item": "B1c", "lado": "B", "texto": "GC 31,58 %", "faixa": "faixa 40–65"},
                {"item": "B1c", "lado": "B", "texto": "ΔG da ponta 5′ −1,01 kcal/mol", "faixa": "máximo −4"},
            ],
        },
        {
            "titulo": "Tm — faixa a calibrar · 6",
            "tipo": "calibrar",
            "itens": [
                {"item": "F3", "lado": "F", "texto": "35,02 °C", "faixa": "faixa 59–61"},
                *TM_CONJUNTO_OK[1:5],
                {"item": "B1c", "lado": "B", "texto": "48,88 °C", "faixa": "faixa 64–66"},
            ],
        },
    ]
    assert dados_cal["primers"][0]["situacao"] == "2 + Tm"
    assert dados_cal["primers"][0]["regiao"] == "—"
    assert dados_cal["primers"][5]["situacao"] == "2 + Tm"

    # Com cal = False
    dados_sem_cal = tela_analise(ex, cfg, tm_como_calibracao=False, escala_kcal=ESCALA_TESTE)
    assert dados_sem_cal["kpis"][1]["detalhe"] == "6 de Tm"
    assert dados_sem_cal["grupos_violacoes"][1]["titulo"] == "Tm fora da faixa · 6"
    assert dados_sem_cal["grupos_violacoes"][1]["tipo"] == "viol"
    assert dados_sem_cal["primers"][0]["situacao"] == "3 violações"
    assert dados_sem_cal["aviso_calibracao"] == ""


def test_tela_desenho():
    sessao = Sessao(RAIZ, cfg, cfg_web)
    sessao.carregar_exemplo("sintético")
    d = sessao.desenhar()

    # Com sel=1, agrupar=False
    dados1 = tela_desenho(d, cfg, sel=1, agrupar=False)
    # fase 4: teto por janela (2 a cada 10 nt) e propostas de núcleos F2–B2c distintos
    assert dados1["candidatos"] == [
        {"texto": "F3 35", "lado": "F"},
        {"texto": "F2 35", "lado": "F"},
        {"texto": "F1 22", "lado": "F"},
        {"texto": "B1c 20", "lado": "B"},
        {"texto": "B2c 32", "lado": "B"},
        {"texto": "B3c 32", "lado": "B"},
    ]
    from interface.formato import milhar
    assert dados1["examinados"] == milhar(d.resultado.estatisticas["conjuntos_examinados"])
    assert d.resultado.estatisticas["conjuntos_examinados"] == 6
    assert dados1["relaxado"] is False
    assert dados1["pode_relaxar"] is False
    assert dados1["vazio"] is False

    lista = dados1["lista"]
    assert [x["num"] for x in lista] == ["1", "2", "3", "4"]
    assert lista[0] == {
        "num": "1",
        "extra": "",
        "selecionada": True,
        "href": "/desenho/1?sel=1&agrupar=0",
        "pontuacao": "0,29",
        "desvio_tm": "0,28 °C",
        "pior_dimero": "−0,03",
        "f1": "76–99",
        "f1_dif": False,
        "b1c": "106–128",
        "b1c_dif": False,
        "lb": "—",
        "lb_dif": False,
        "avisos": "sem LF, sem LB",
    }
    assert lista[1]["b1c"] == "119–140 (+13)"
    assert lista[1]["b1c_dif"] is True
    assert lista[2]["f1"] == "92–113 (+16)"
    assert lista[2]["f1_dif"] is True
    assert lista[3]["pontuacao"] == "5,97"
    assert lista[3]["b1c"] == "122–145 (+16)"

    assert dados1["nota_empate"] == ""
    assert dados1["href_alternar_agrupar"] == "/desenho/1?sel=1&agrupar=1"

    prop = dados1["proposta"]
    assert prop is not None
    assert prop["primers"][2] == {
        "nome": "FIP",
        "lado": "X",
        "seq": "ACCCCCAAGTGTTATTGGCAAGTT",
        "seq2": "AGACACGCATATGGATTCTGTAC",
        "regiao": "F1 76–99 · F2 27–49",
        "nt": "47",
        "copiar": "ACCCCCAAGTGTTATTGGCAAGTTAGACACGCATATGGATTCTGTAC",
    }
    assert prop["primers"][4] == {
        "nome": "LF",
        "lado": "L",
        "seq": "—",
        "seq2": "",
        "regiao": "sem candidato",
        "nt": "—",
        "copiar": "",
    }
    assert prop["primers"][5] == {
        "nome": "LB",
        "lado": "L",
        "seq": "—",
        "seq2": "",
        "regiao": "sem candidato",
        "nt": "—",
        "copiar": "",
    }
    assert prop["nota"] == (
        "Sem LF: nenhuma janela de 15–25 nt entre F2 e F1 (posições 50–75) passou nas regras de loop. "
        "Sem LB: nenhuma janela de 15–25 nt entre B1c e B2c (posições 129–145) passou nas regras de loop. "
        "FIP = F1c + F2 e BIP = B1c + B2, sem ligante."
    )
    assert [s["nome"] for s in prop["mapa"]["segmentos"]] == [
        "F3",
        "F2",
        "sem LF",
        "F1",
        "B1c",
        "sem LB",
        "B2c",
        "B3c",
    ]
    assert prop["mapa"]["segmentos"][2] == {
        "nome": "sem LF",
        "regiao": "LF",
        "ini": 49,
        "fim": 75,
        "esquerda": "22.685",
        "largura": "12.037",
        "fantasma": True,
    }

    # Com sel=3, agrupar=True: núcleos distintos não formam grupos
    dados_grp = tela_desenho(d, cfg, sel=3, agrupar=True)
    assert dados_grp["sel"] == 3
    assert len(dados_grp["lista"]) == 4
    assert dados_grp["lista"][0]["extra"] == ""
    assert dados_grp["lista"][0]["href"] == "/desenho/1?sel=1&agrupar=1"
    assert dados_grp["href_alternar_agrupar"] == "/desenho/1?sel=1&agrupar=0"

    # Com sel=2, agrupar=False
    dados_sel2 = tela_desenho(d, cfg, sel=2, agrupar=False)
    assert dados_sel2["proposta"]["num"] == "2"
    assert dados_sel2["lista"][1]["selecionada"] is True
    assert dados_sel2["proposta"]["primers"][3]["seq"] == "CTCCGTCCACCAGTGTGAATCG"

    # Com sel=99
    dados_sel99 = tela_desenho(d, cfg, sel=99, agrupar=False)
    assert dados_sel99["sel"] == 1


def test_tela_parametros():
    sessao_nova = Sessao(RAIZ, cfg, cfg_web)
    dados_nova = tela_parametros(sessao_nova, cfg)

    assert dados_nova["referencia"] is None
    assert dados_nova["calibracao"] is None
    assert dados_nova["distancias"][0] == {
        "nome": "F2 → B2",
        "faixa": "120–180",
        "valor": "—",
        "dentro": None,
    }
    assert dados_nova["reacao"][0] == {
        "nome": "Monovalentes",
        "valor": "70 mM",
        "nota": (
            "Na⁺ 0 + K⁺ 50 + NH₄⁺ 20 mM. Conferir: a mistura WarmStart é proprietária; "
            "o tampão isotérmico 1X da NEB tem 50 mM de KCl e 10 mM de (NH₄)₂SO₄."
        ),
        "destaque": False,
    }
    assert dados_nova["reacao"][2] == {
        "nome": "dNTP",
        "valor": "1,4 mM cada",
        "nota": "Total de 5,6 mM, o valor enviado ao primer3.",
        "destaque": False,
    }
    assert dados_nova["regras"][1] == {"nome": "GC", "valor": "40–65 %", "marca": "conferir"}
    assert dados_nova["desenho"][3] == {
        "nome": "Peso do dímero",
        "valor": "0,5 por kcal/mol",
        "marca": "provisório",
    }
    assert dados_nova["desenho"][5] == {"nome": "Ligante FIP/BIP", "valor": "nenhum", "marca": ""}

    # Após analisar o exemplo
    sessao = Sessao(RAIZ, cfg, cfg_web)
    sessao.carregar_exemplo("sintético")
    sessao.analisar()

    dados_analisado = tela_parametros(sessao, cfg)
    assert dados_analisado["referencia"] == {"rotulo": "conjunto_ok × alvo_ok", "id": 1}

    cal = dados_analisado["calibracao"]
    assert cal is not None
    assert [linha["nome"] for linha in cal["linhas"]] == [
        "Externos",
        "Internos, ponta 3′",
        "Internos, ponta 5′",
    ]
    assert cal["linhas"][0]["faixa"] == "59–61"
    assert cal["linhas"][0]["conc"] == "200 nM"
    assert cal["linhas"][0]["banda_esquerda"] == "39.286"
    assert cal["linhas"][0]["banda_largura"] == "7.143"
    assert cal["linhas"][0]["pontos"][0] == {
        "nome": "F3",
        "lado": "F",
        "cheio": True,
        "esquerda": "24.792",
        "rotulo": True,
    }

    assert [t["valor"] for t in cal["ticks"]] == [
        "48",
        "50",
        "52",
        "54",
        "56",
        "58",
        "60",
        "62",
        "64",
        "66",
        "68",
        "70",
        "72",
        "74",
        "76 °C",
    ]
    assert sorted(cal) == ["descricao", "linhas", "texto", "ticks"]
    assert cal["texto"] == (
        "Nenhum dos 6 primers medidos fica dentro da faixa. "
        "Na reação, a Tm dos primers em duplex fica em média 5,6 °C acima da Tm de referência. "
        "F1c e B1c entram na reação como dobramento, fora da média."
    )

    assert dados_analisado["distancias"][0]["valor"] == "146"
    assert dados_analisado["distancias"][0]["dentro"] is True
    assert dados_analisado["reacao"][2]["nota"] == "Total de 5,6 mM, o valor enviado ao primer3."


def test_tela_base_etapas_e_espera() -> None:
    """Verifica estados das etapas e valores de espera padrão na tela base."""
    sessao_vazia = Sessao(RAIZ, cfg, cfg_web)
    ctx_inicio = tela_base("inicio", sessao_vazia, ("feat/esteira", "abc1234"))
    assert [e["estado"] for e in ctx_inicio["etapas"]] == ["atual", "pendente", "pendente"]
    assert [e["nome"] for e in ctx_inicio["etapas"]] == ["Dados", "Análise", "Desenho"]
    assert [e["href"] for e in ctx_inicio["etapas"]] == ["/", "/analise", "/desenho"]
    assert [e["numero"] for e in ctx_inicio["etapas"]] == [1, 2, 3]
    assert ctx_inicio["espera"] == {"analise_ms": 2600, "desenho_ms": 3200, "envio_ms": 900}

    sessao = Sessao(RAIZ, cfg, cfg_web)
    sessao.carregar_exemplo("sintético")
    sessao.analisar()

    ctx_analise = tela_base("analise", sessao, ("feat/esteira", "abc1234"))
    assert [e["estado"] for e in ctx_analise["etapas"]] == ["feita", "atual", "pendente"]

    ctx_param = tela_base("parametros", sessao, ("feat/esteira", "abc1234"))
    assert [e["estado"] for e in ctx_param["etapas"]] == ["feita", "feita", "pendente"]

    sessao.desenhar()
    ctx_desenho = tela_base("desenho", sessao, ("feat/esteira", "abc1234"))
    assert [e["estado"] for e in ctx_desenho["etapas"]] == ["feita", "feita", "atual"]


def test_tela_base_espera_do_config() -> None:
    """Verifica carregamento de tempos de espera a partir da configuração web."""
    cfg_com_espera = dict(cfg_web)
    cfg_com_espera["carregamento"] = {"analise_ms": 10, "desenho_ms": 20, "envio_ms": 0}
    sessao = Sessao(RAIZ, cfg, cfg_com_espera)
    ctx = tela_base("inicio", sessao, ("feat/esteira", "abc1234"))
    assert ctx["espera"] == {"analise_ms": 10, "desenho_ms": 20, "envio_ms": 0}


def test_tela_tubos_padrao() -> None:
    """Sem análise: curvas manuais em minutos, com aviso, e cores do modelo do lampsim."""
    dados = tela_tubos({}, CFG_ESTEIRA)
    assert dados["modo"] == "manual"
    assert dados["aviso"].startswith("Nenhuma análise nesta sessão")
    assert dados["duracao_min"] == 30.0
    assert dados["duracao_max_min"] == 60.0
    assert dados["passo_min"] == 0.5
    assert dados["leitura_min"] == 15.0
    assert dados["velocidade"] == 60.0
    assert dados["velocidades"] == [30.0, 60.0, 120.0, 300.0]
    assert len(dados["paleta"]) == 101
    assert dados["cor_inicial"] == "#FE00B1"
    assert dados["cor_final"] == "#F3DF48"
    assert dados["t_min"] == []
    assert len(dados["tubos"]) == 9
    assert dados["tubos"][0] == {
        "nome": "Amostra 1",
        "curva": "linear",
        "inicio": 0.0,
        "duracao": 30.0,
        "alcance": 100.0,
    }
    assert dados["tubos"][8] == {
        "nome": "Referência",
        "curva": "linear",
        "inicio": 0.0,
        "duracao": 30.0,
        "alcance": 0.0,
    }


def test_tela_tubos_config_proprio_e_validacoes() -> None:
    """Verifica respeito a config próprio e sanitização de valores inválidos."""
    cfg_custom = {
        "tubos": {
            "duracao_min": 40,
            "velocidade": "rápido",
            "manual": {"ph_inicial": -1},
            "exemplares": [
                {
                    "nome": "Amostra Custom",
                    "curva": "x",
                    "inicio": 0,
                    "duracao": 0,
                    "alcance": 150,
                }
            ],
        }
    }
    dados = tela_tubos(cfg_custom, CFG_ESTEIRA)
    assert dados["duracao_min"] == 40.0
    assert dados["velocidade"] == 60.0
    assert dados["cor_inicial"] == "#FE00B1"
    assert len(dados["tubos"]) == 1
    assert dados["tubos"][0] == {
        "nome": "Amostra Custom",
        "curva": "linear",
        "inicio": 0.0,
        "duracao": 1.0,
        "alcance": 100.0,
    }


def test_tela_tubos_exemplares_vazio_e_limites_minimos() -> None:
    """Verifica que exemplares vazio cai nos 9 padrão e que limites inferiores são respeitados."""
    dados_vazio = tela_tubos({"tubos": {"exemplares": []}}, CFG_ESTEIRA)
    assert len(dados_vazio["tubos"]) == 9

    cfg_inferior = {
        "tubos": {
            "duracao_min": -10,
            "duracao_max_min": 20,
            "exemplares": [
                {
                    "nome": "Amostra Min",
                    "curva": "sigmoide",
                    "inicio": -5,
                    "duracao": -2,
                    "alcance": -50,
                }
            ],
        }
    }
    dados = tela_tubos(cfg_inferior, CFG_ESTEIRA)
    assert dados["duracao_min"] == 20.0
    assert dados["duracao_max_min"] == 20.0
    assert len(dados["tubos"]) == 1
    assert dados["tubos"][0] == {
        "nome": "Amostra Min",
        "curva": "sigmoide",
        "inicio": 0.0,
        "duracao": 1.0,
        "alcance": 0.0,
    }

