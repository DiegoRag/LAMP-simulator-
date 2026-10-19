"""Testes da Fase 1 - Sessão 3: coordenadas base 1, ponta 3', distâncias com aproximação e fita -."""

from pathlib import Path

from Bio import SeqIO

from interface.formato import coord
from interface.sessao import Sessao
from interface.telas import tela_analise
from primers.entrada import carregar_config
from primers.modelo import rc

RAIZ = Path(__file__).resolve().parents[1]
cfg = carregar_config(RAIZ / "config" / "primers.toml")
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
    ],
}
ESCALA_TESTE = [-10.0, -8.0, -6.0, -4.0]


def test_coord_base_1_genbank():
    assert coord(0, 21) == "1–21"
    assert coord(31, 49) == "32–49"
    assert coord(187, 206) == "188–206"


def test_tela_analise_regiao_com_um_erro_mantem_7_distancias_e_mostra_ponta_3p():
    sessao = Sessao(RAIZ, cfg, cfg_web)
    conteudo_conj = (RAIZ / "tests/dados/conjunto_ok.csv").read_bytes()
    sessao.carregar_conjunto("conjunto_ok.csv", conteudo_conj)

    registros = list(SeqIO.parse(str(RAIZ / "tests/dados/alvo_ok.fasta"), "fasta"))
    seq_orig = str(registros[0].seq)
    assert seq_orig[48] == "C"
    seq_mut = seq_orig[:48] + "A" + seq_orig[49:]
    alvo_bytes = f">{registros[0].id}\n{seq_mut}\n".encode("utf-8")
    sessao.carregar_alvo("alvo_mut.fasta", alvo_bytes)

    ex = sessao.analisar()
    dados = tela_analise(ex, cfg, tm_como_calibracao=True, escala_kcal=ESCALA_TESTE)

    assert dados["kpis"][2] == {
        "rotulo": "Regiões localizadas",
        "valor": "5 de 6",
        "detalhe": "F2 com 1 dif., 1 na ponta",  # fase 6: a base trocada é a 3' terminal de F2
        "tipo": "neutro",
    }
    assert dados["kpis"][3]["valor"] == "7 de 7"
    assert len(dados["distancias"]) == 7

    nomes_dist = [d["nome"] for d in dados["distancias"]]
    for esperado in ("F2 → B2", "F2 → F1", "F3 → F2", "F3 → B3"):
        assert esperado in nomes_dist

    itens_viol = dados["grupos_violacoes"][0]["itens"]
    viol_f2 = next(item for item in itens_viol if item["item"] == "F2")
    # fase 6: o erro é contado no primer (base 18 de 18 do F2) e situado em relação à ponta 3′
    assert viol_f2["texto"] == (
        "Não localizado: melhor posição 32–49 (fita +), 1 diferença (base 18 de 18 do F2: C no primer, "
        "o alvo pede A, a 0 nt da ponta 3′, dentro das 5 bases da ponta; posição 49)"
    )


def test_tela_analise_alvo_fita_menos():
    sessao = Sessao(RAIZ, cfg, cfg_web)
    conteudo_conj = (RAIZ / "tests/dados/conjunto_ok.csv").read_bytes()
    sessao.carregar_conjunto("conjunto_ok.csv", conteudo_conj)

    registros = list(SeqIO.parse(str(RAIZ / "tests/dados/alvo_ok.fasta"), "fasta"))
    seq_rc = rc(str(registros[0].seq))
    alvo_bytes = f">{registros[0].id}_rc\n{seq_rc}\n".encode("utf-8")
    sessao.carregar_alvo("alvo_rc.fasta", alvo_bytes)

    ex = sessao.analisar()
    dados = tela_analise(ex, cfg, tm_como_calibracao=True, escala_kcal=ESCALA_TESTE)

    assert "fita −" in dados["subtitulo"]
    assert dados["kpis"][2]["valor"] == "6 de 6"
    assert dados["kpis"][3]["valor"] == "7 de 7"
    assert any("fita complementar (−)" in a for a in dados["avisos"])
