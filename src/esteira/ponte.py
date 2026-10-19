"""Ponte entre a análise de primers e o simulador: incompatibilidades -> eficiência -> tubos.

É a única parte que conhece os dois lados: lê o que primers.analise mediu (primers.incompat)
e entrega ao lampsim a eficiência e o número de cópias de cada tubo. Os parâmetros do modelo
de eficiência ficam em config/esteira.toml, todos PROVISÓRIOS.

O contrato fecha a esteira (fase 6): selecao_da_analise escreve o selecao.json (selecao@2)
de uma análise, com incompat_* e f_previsto, e tubos_da_selecao lê esse arquivo e simula.
"""

import hashlib
import json
import re
import tomllib
from pathlib import Path

import numpy as np

from esteira.contrato import ESQUEMA_V2, validar_selecao
from lampsim.classes import amostra_de_reacao, simular
from lampsim.config import CINETICA, ESPECTRO, PAPEIS, QUIMICA
from lampsim.cor import cor_do_indicador, rgb_para_hex
from lampsim.quimica import fracao_acida
from primers.incompat import LACOS, Incompat, regioes_para_json

CONFIG_PADRAO = Path(__file__).resolve().parents[2] / "config" / "esteira.toml"
T0_MAX = float(CINETICA["t0_max_min"])
# número de acesso com versão (GenBank, RefSeq): AP018671.1, NZ_CP012345.1, NC_000913.3
_ACESSO = re.compile(r"^(?:[A-Z]{2}_)?[A-Z]{1,6}\d{5,9}\.\d+$")


def carregar_config(caminho: str | Path | None = None) -> dict:
    """Lê config/esteira.toml (ou o caminho dado)."""
    return tomllib.loads(Path(caminho or CONFIG_PADRAO).read_text(encoding="utf-8"))


def eficiencia_prevista(incompat: Incompat | None, cfg: dict) -> float | None:
    """Eficiência de 0 a 1 prevista pelas incompatibilidades; None se não houve medida."""
    if incompat is None:
        return None
    par = cfg["eficiencia"]
    f3, fo, piso = float(par["fator_ponta3"]), float(par["fator_outra"]), float(par["piso_laco"])
    eficiencia = 1.0
    for item in incompat.por_regiao:
        if item.encontrada:
            n3 = int(item.ponta3 or 0)
            fator = f3 ** n3 * fo ** (int(item.total or 0) - n3)
        else:
            fator = 0.0
        if item.regiao in LACOS:
            fator = max(piso, fator)
        eficiencia *= fator
    return float(min(1.0, max(0.0, eficiencia)))


def medida_por_alelo(incompat: Incompat | None, alelo_id: str, acesso: str | None, cfg: dict) -> dict:
    """Uma entrada de por_alelo do contrato selecao.json (esteira.contrato)."""
    return {
        "alelo_id": alelo_id,
        "acesso": acesso,
        "incompat_total": None if incompat is None else incompat.total,
        "incompat_3p": None if incompat is None else incompat.ponta3,
        "f_previsto": eficiencia_prevista(incompat, cfg),
    }


def acesso_do_id(id_registro: str | None) -> str | None:
    """O id do FASTA, se tiver a forma de um número de acesso com versão; senão None."""
    if id_registro and _ACESSO.match(id_registro):
        return id_registro
    return None


def assinatura(versao_codigo: str, experimento: str, configs: dict) -> dict:
    """Campo gerado_por do contrato: versão do código, experimento e SHA-256 dos configs."""
    texto = json.dumps(configs, sort_keys=True, ensure_ascii=False)
    return {
        "versao_codigo": versao_codigo,
        "experimento": experimento,
        "config_sha256": hashlib.sha256(texto.encode("utf-8")).hexdigest(),
    }


def selecao_da_analise(
    *,
    conjunto_id: str,
    origem: str,
    primers: dict[str, str],
    incompat: Incompat | None,
    alelo_id: str,
    acesso: str | None,
    cfg: dict,
    gerado_por: dict,
    aviso: str | None = None,
) -> dict:
    """selecao.json (selecao@2) de uma análise: um conjunto medido contra um alelo.

    A medida traz incompat_total, incompat_3p, f_previsto e, com medida, incompat_regioes.
    selecionado = eficiência prevista > 0; o motivo diz o porquê. Levanta ValueError se o
    resultado não passar em esteira.contrato.validar_selecao.
    """
    medida = medida_por_alelo(incompat, alelo_id, acesso, cfg)
    if incompat is not None:
        medida["incompat_regioes"] = regioes_para_json(incompat)
    f = medida["f_previsto"]
    if incompat is None:
        motivo = "sem medida contra o alvo: nada a prever"
    else:
        partes = [f"eficiência prevista {f:.3g} (modelo provisório de config/esteira.toml)"]
        if incompat.total is not None:
            partes.append(f"{incompat.total} incompatibilidades, {incompat.ponta3} na ponta que estende")
        if incompat.nao_encontradas:
            partes.append(f"sem posição no alvo: {', '.join(incompat.nao_encontradas)}")
        motivo = "; ".join(partes)
    conjunto = {
        "conjunto_id": conjunto_id,
        "origem": origem,
        "primers": dict(primers),
        "por_alelo": [medida],
        "selecionado": f is not None and f > 0.0,
        "motivo": motivo,
    }
    if aviso is not None:
        conjunto["aviso"] = aviso
    dados = {"esquema": ESQUEMA_V2, "gerado_por": dict(gerado_por), "conjuntos": [conjunto]}
    problemas = validar_selecao(dados)
    if problemas:
        raise ValueError("selecao.json inválido:\n- " + "\n- ".join(problemas))
    return dados


def medida_da_selecao(dados: dict, conjunto_id: str | None = None, alelo_id: str | None = None) -> dict:
    """A medida por_alelo de um conjunto (o primeiro, se não disser qual) contra um alelo."""
    problemas = validar_selecao(dados)
    if problemas:
        raise ValueError("selecao.json inválido:\n- " + "\n- ".join(problemas))
    conjuntos = dados["conjuntos"]
    if conjunto_id is None:
        conjunto = conjuntos[0]
    else:
        conjunto = next((c for c in conjuntos if c["conjunto_id"] == conjunto_id), None)
        if conjunto is None:
            raise KeyError(f"conjunto {conjunto_id!r} não está na seleção")
    medidas = conjunto["por_alelo"]
    if alelo_id is None:
        return medidas[0]
    medida = next((m for m in medidas if m["alelo_id"] == alelo_id), None)
    if medida is None:
        raise KeyError(f"alelo {alelo_id!r} não está no conjunto {conjunto['conjunto_id']!r}")
    return medida


def tubos_da_selecao(
    dados: dict,
    tubos: list[dict],
    t_min,
    semente: int,
    conjunto_id: str | None = None,
    alelo_id: str | None = None,
) -> list[dict]:
    """Lê o f_previsto do selecao.json e simula os tubos com ele (o lado leitor do contrato)."""
    medida = medida_da_selecao(dados, conjunto_id, alelo_id)
    f = medida.get("f_previsto")
    if f is None:
        raise ValueError(f"alelo {medida['alelo_id']!r} sem f_previsto: não há o que simular")
    return simular_tubos(float(f), tubos, t_min, semente)


def eixo_minutos(duracao_min: float, passo_min: float) -> np.ndarray:
    """Minutos de 0 até a duração, inclusive, no passo dado."""
    n = int(round(duracao_min / passo_min)) + 1
    return np.linspace(0.0, float(duracao_min), n)


def sementes_dos_tubos(semente: int, n: int) -> list[int]:
    """Uma semente por tubo, derivada da semente da série (numpy SeedSequence)."""
    return [int(s.generate_state(1)[0]) for s in np.random.SeedSequence(semente).spawn(n)]


def simular_tubos(eficiencia: float, tubos: list[dict], t_min, semente: int) -> list[dict]:
    """Simula cada tubo da série com o lampsim.

    tubos: dicionários com "nome", "copias" (média por tubo), "papel" (lampsim PAPEIS) e,
    opcional, "eficiencia" (para um controle positivo com molde próprio). A eficiência do
    conjunto vale para os tubos sem eficiência própria.
    """
    t = np.asarray(t_min, dtype=float)
    saida = []
    for tubo, s in zip(tubos, sementes_dos_tubos(semente, len(tubos)), strict=True):
        papel = tubo.get("papel", "amostra")
        if papel not in PAPEIS:
            raise ValueError(f"papel desconhecido: {papel!r}; opções: {list(PAPEIS)}")
        ef = tubo.get("eficiencia")
        ef = float(eficiencia if ef is None else ef)
        amostra = amostra_de_reacao(ef, float(tubo["copias"]), s, papel=papel, id_amostra=str(tubo["nome"]))
        traj = simular(amostra, t=t)
        saida.append({
            "nome": str(tubo["nome"]),
            "papel": papel,
            "copias_media": float(tubo["copias"]),
            "copias": amostra.copias,
            "eficiencia": ef,
            "t0_min": None if amostra.A == 0.0 or amostra.t0 >= T0_MAX else round(amostra.t0, 2),
            "diagnostico": amostra.diagnostico,
            "ph": [round(float(v), 3) for v in traj["ph"]],
            "cores": rgb_para_hex(traj["rgb"]),
        })
    return saida


def paleta_progresso(ph_inicial: float, queda_ph: float, n: int = 101) -> list[str]:
    """Cores do indicador do início ao fim de uma queda de pH, em n passos de progresso.

    Serve às curvas manuais da página Tubos: o progresso p (0 a 1) vira pH = pH_inicial −
    queda · p e a cor sai do mesmo modelo do lampsim, com a turbidez do produto.
    """
    p = np.linspace(0.0, 1.0, n)
    f = fracao_acida(ph_inicial - queda_ph * p, QUIMICA["pka"])
    return rgb_para_hex(cor_do_indicador(f, turbidez=ESPECTRO["turbidez_max"] * p))
