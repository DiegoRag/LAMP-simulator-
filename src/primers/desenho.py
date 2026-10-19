"""Desenho de conjuntos de primers LAMP no alvo inteiro (orquestra candidatos, núcleos e ranking)."""

from dataclasses import dataclass

from primers.candidatos import INFO_TIPO, TIPOS_BASE, Candidato, gerar_candidatos, reprova_homodimero
from primers.diagnostico import AVISO_RELAXADO, NOMES_REGRA, config_relaxada, diagnosticar
from primers.janelas import IndiceAlvo
from primers.modelo import ConjuntoLAMP, Regiao, normalizar, primers_de_regioes
from primers.nucleos import Nucleo, escolher_distintos, melhores_nucleos
from primers.restricoes import verificar
from primers.termo import MAX_NT_ESTRUTURA, Condicoes, dimero_ponta3

__all__ = [
    "AVISO",
    "AVISO_ASCII",
    "AVISO_RELAXADO",
    "INFO_TIPO",
    "NOMES_REGRA",
    "Candidato",
    "ConjuntoDesenhado",
    "ResultadoDesenho",
    "config_relaxada",
    "desenhar",
    "diagnosticar",
    "gerar_candidatos",
    "pior_dimero_tubo",
]

AVISO = (
    "Candidatos computacionais, não validados experimentalmente. A definição do conjunto "
    "que vai ao tubo é da equipe biomédica (ADR-0013)."
)

AVISO_ASCII = (
    "Candidatos computacionais, nao validados experimentalmente. A definicao do conjunto "
    "que vai ao tubo e da equipe biomedica (ADR-0013)."
)


@dataclass
class ConjuntoDesenhado:
    """Conjunto LAMP desenhado e ranqueado."""

    posicao: int
    pontuacao: float
    desvio_tm: float
    pior_heterodimero_kcal: float
    conjunto: ConjuntoLAMP
    regioes: dict[str, Regiao]
    avisos: list[str]


@dataclass
class ResultadoDesenho:
    """Resultado completo da etapa de desenho de primers LAMP."""

    conjuntos: list[ConjuntoDesenhado]
    estatisticas: dict


def pior_dimero_tubo(tubo: dict[str, str], cond: Condicoes, max_nt: int) -> tuple[float, list[str]]:
    """Menor ΔG (kcal/mol) de dimero_ponta3 entre os pares de primers do tubo, e os avisos.

    Primer com mais de max_nt nt (limite do primer3) entra pelos seus max_nt nt da ponta 3′, a que a
    polimerase estende, com um aviso. Sem pares, devolve 0.0.
    """
    avaliados: dict[str, str] = {}
    avisos: list[str] = []
    for nome, seq in tubo.items():
        if len(seq) > max_nt:
            avaliados[nome] = seq[-max_nt:]
            avisos.append(
                f"{nome} com {len(seq)} nt: dímeros avaliados nos {max_nt} nt "
                "da ponta 3′ (limite do primer3)"
            )
        else:
            avaliados[nome] = seq
    nomes = list(avaliados)
    dgs = [
        dimero_ponta3(avaliados[a], avaliados[b], cond).dg_kcal
        for i, a in enumerate(nomes)
        for b in nomes[i + 1 :]
    ]
    return (min(dgs) if dgs else 0.0), avisos


def desenhar(alvo: str, cfg: dict, n: int | None = None, *, relaxado: bool = False) -> ResultadoDesenho:
    """Desenha conjuntos de primers LAMP para o alvo inteiro.

    1. Candidatos com teto por janela (gerar_candidatos), o homodímero deixado para o passo 3.
    2. melhores_nucleos: todos os núcleos (F2, B2c) do alvo, os nucleos_guardados melhores na etapa 1.
    3. escolher_distintos: até pre_selecao núcleos que não repetem F2 e B2c e cujos seis candidatos passam
       no homodímero.
    4. Etapa 2: LF e LB (melhor desvio_tm que passa no homodímero), pior dímero com ponta 3′ entre todos
       os primers do tubo e pontuação = desvio_tm + pen_f2_b2 + peso_dimero × max(0, −pior).
    5. Os n melhores, autoverificados com verificar(). relaxado=True usa config_relaxada(cfg) em tudo.
    """
    alvo = normalizar(alvo)
    cfg_uso = config_relaxada(cfg) if relaxado else cfg
    d = cfg_uso["desenho"]
    cond = Condicoes.de_config(cfg_uso)
    max_nt = min(int(cfg_uso.get("dimeros", {}).get("max_nt_estrutura", MAX_NT_ESTRUTURA)), MAX_NT_ESTRUTURA)
    indice = IndiceAlvo(alvo, cfg_uso)
    eliminados: dict[str, dict[str, int]] = {}
    cands = gerar_candidatos(alvo, cfg_uso, indice=indice, eliminados=eliminados, testar_homodimero=False)

    homo: dict[tuple[str, Regiao], bool] = {}

    def homo_ok(c: Candidato) -> bool:
        chave = (c.tipo, c.regiao)
        if chave not in homo:
            homo[chave] = not reprova_homodimero(c.primer, c.tipo, cfg_uso, cond)
            if not homo[chave] and c.tipo in eliminados:
                eliminados[c.tipo]["homodimero"] += 1
        return homo[chave]

    guardados, total_nucleos = melhores_nucleos(cands, cfg_uso)
    distintos = escolher_distintos(guardados, cfg_uso, lambda nu: all(homo_ok(c) for c in nu.combo))

    def melhor_laco(tipo: str, ini: int, fim: int) -> Regiao | None:
        if fim <= ini:
            return None
        lista = gerar_candidatos(
            alvo, cfg_uso, tipos=(tipo,), dentro=Regiao(ini, fim), indice=indice, testar_homodimero=False
        )[tipo]
        for c in sorted(lista, key=lambda c: (c.desvio_tm, c.regiao.ini, c.regiao.fim)):
            if homo_ok(c):
                return c.regiao
        return None

    etapa2: list[tuple[float, Nucleo, float, ConjuntoLAMP, dict[str, Regiao], list[str]]] = []
    for nu in distintos:
        _, c_f2, c_f1, c_b1c, c_b2c, _ = nu.combo
        regioes = {t: c.regiao for t, c in zip(TIPOS_BASE, nu.combo, strict=True)}
        avisos: list[str] = []
        for tipo, ini, fim in (
            ("LF", c_f2.regiao.fim, c_f1.regiao.ini),
            ("LB", c_b1c.regiao.fim, c_b2c.regiao.ini),
        ):
            r_laco = melhor_laco(tipo, ini, fim)
            if r_laco is None:
                avisos.append(f"sem {tipo} candidato")
            else:
                regioes[tipo] = r_laco
        conjunto = primers_de_regioes(alvo, regioes, ligante=d["ligante"])
        tubo = {"F3": conjunto.F3, "B3": conjunto.B3, "FIP": conjunto.FIP, "BIP": conjunto.BIP}
        if conjunto.LF is not None:
            tubo["LF"] = conjunto.LF
        if conjunto.LB is not None:
            tubo["LB"] = conjunto.LB
        pior, avisos_dim = pior_dimero_tubo(tubo, cond, max_nt)
        avisos.extend(avisos_dim)
        if relaxado:
            avisos.append(AVISO_RELAXADO)
        pontuacao = nu.desvio_tm + nu.pen_f2_b2 + float(d["peso_dimero"]) * max(0.0, -pior)
        etapa2.append((pontuacao, nu, pior, conjunto, regioes, avisos))
    etapa2.sort(key=lambda x: (x[0], x[1].chave))

    limite_n = int(d["n_resultados"]) if n is None else n
    conjuntos: list[ConjuntoDesenhado] = []
    for posicao, (pontuacao, nu, pior, conjunto, regioes, avisos) in enumerate(etapa2[:limite_n], start=1):
        violacoes = verificar(conjunto, alvo, cfg_uso)
        if violacoes:
            raise RuntimeError(f"bug no desenho: {violacoes}")
        conjuntos.append(
            ConjuntoDesenhado(
                posicao=posicao,
                pontuacao=pontuacao,
                desvio_tm=nu.desvio_tm,
                pior_heterodimero_kcal=pior,
                conjunto=conjunto,
                regioes=regioes,
                avisos=avisos,
            )
        )

    stats = {
        "candidatos_por_tipo": {t: len(cands[t]) for t in TIPOS_BASE},
        "eliminados_por_tipo": {t: dict(eliminados[t]) for t in TIPOS_BASE},
        "conjuntos_examinados": total_nucleos,
        "nucleos_guardados": len(guardados),
        "nucleos_distintos": len(distintos),
        "relaxado": relaxado,
        "diagnostico": "",
    }
    if not conjuntos:
        stats["diagnostico"] = diagnosticar(stats, relaxado)
    return ResultadoDesenho(conjuntos=conjuntos, estatisticas=stats)
