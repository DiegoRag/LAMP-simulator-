"""Contrato entre a etapa 3 (primers) e a etapa 4 (lampsim): o arquivo selecao.json.

A etapa 3 escreve o arquivo e a esteira o lê. É a única coisa que um lado sabe do outro: o
primers não importa o lampsim e o lampsim não importa o primers (tests/test_fronteiras.py).

Versão 1, congelada em 16/10/2026 (roadmap da esteira). Mudar um campo é subir a versão do
esquema, nunca alterar a v1 em silêncio: o leitor recusa versão que não conhece.

    esquema           "lamp-sim/selecao@1"
    gerado_por        versao_codigo, experimento, config_sha256 (texto)
    conjuntos         lista não vazia; em cada conjunto:
      conjunto_id       nome do conjunto, com ponto ("fosA.1")
      origem            "artigo" (validado em bancada) ou "desenho" (candidato do desenhar)
      aviso             obrigatório quando a origem é "desenho" (ADR-0013, regra 2)
      primers           F3, B3, FIP e BIP obrigatórios; LF e LB opcionais; só ACGT, 5'->3'
      por_alelo         lista não vazia; em cada medida contra um alelo:
        alelo_id          nome do alelo, sem ponto ("fosA3")
        acesso            número de acesso público, ou null (alvo sintético)
        incompat_total    bases do conjunto que não casam com o alelo, ou null (não medido)
        incompat_3p       dessas, quantas caem nas 5 bases da ponta 3', ou null
        f_previsto        opcional: eficiência prevista por um modelo, de 0 a 1, ou null
      selecionado       true ou false
      motivo            texto

Versão 2 ("lamp-sim/selecao@2"), o que muda em relação à v1:
    conjunto_id       qualquer texto não vazio (com ou sem ponto)
    origem            "artigo", "desenho" ou "usuario" (conjunto carregado, sem afirmação de bancada)
    alelo_id          qualquer texto não vazio (permitindo ponto)
    primers           F3, B3, FIP, BIP ou 6 partes F3, B3, F2, F1c, B2, B1c; bases IUPAC/U
    incompat_3p       incompatibilidades a menos de [casamento] bases_ponta_3p nt da ponta que
                      ESTENDE: a ponta 3' do primer em F3, F2, B2, B3, LF e LB; a ponta 5' de F1c
                      e de B1c, que no haltere pareia com a ponta 3' da fita e fecha o laço
    incompat_regioes  opcional: lista, uma entrada por região do conjunto, com
                        regiao      "F3", "F2", "F1", "B1c", "B2c", "B3c", "LF" ou "LB", sem repetir
                        encontrada  true ou false
                        total       incompatibilidades da região (null se não encontrada)
                        ponta       dessas, quantas na ponta que estende (null se não encontrada)
                      Com todas encontradas, incompat_total e incompat_3p são as somas; com
                      alguma não encontrada, os dois são null. Assim um erro no meio de F1c e um
                      erro na ponta 5' de F1c, que impede o laço de fechar, ficam distintos.
"""

import json
from pathlib import Path

ESQUEMA = "lamp-sim/selecao@1"
ESQUEMA_V2 = "lamp-sim/selecao@2"
ESQUEMAS_ACEITOS = (ESQUEMA, ESQUEMA_V2)
ORIGENS = ("artigo", "desenho")
ORIGENS_V2 = ("artigo", "desenho", "usuario")
REGIOES_INCOMPAT = ("F3", "F2", "F1", "B1c", "B2c", "B3c", "LF", "LB")
PRIMERS_OBRIGATORIOS = ("F3", "B3", "FIP", "BIP")
PRIMERS_OPCIONAIS = ("LF", "LB")
PRIMERS_PARTES_V2 = ("F3", "B3", "F2", "F1c", "B2", "B1c")
PRIMERS_OPCIONAIS_V2 = ("LF", "LB", "FLP", "BLP", "LoopF", "LoopB", "ligante")
_BASES_V2 = frozenset("ACGTURYKMSWBDHVN")

_CAMPOS_RAIZ = {"esquema", "gerado_por", "conjuntos"}
_CAMPOS_GERADO_POR = {"versao_codigo", "experimento", "config_sha256"}
_CAMPOS_CONJUNTO = {"conjunto_id", "origem", "primers", "por_alelo", "selecionado", "motivo"}
_CAMPOS_CONJUNTO_OPCIONAIS = {"aviso"}
_CAMPOS_ALELO = {"alelo_id", "acesso", "incompat_total", "incompat_3p"}
_CAMPOS_ALELO_OPCIONAIS = {"f_previsto"}
_CAMPOS_ALELO_OPCIONAIS_V2 = {"f_previsto", "incompat_regioes"}
_CAMPOS_REGIAO_V2 = {"regiao", "encontrada", "total", "ponta"}


def validar_selecao(dados: object) -> list[str]:
    """Confere um selecao.json já lido com json.load. Devolve os problemas; vazia = válido.

    Campo desconhecido também é problema: um nome digitado errado ("mismatches_3p") que
    passasse calado faria o leitor cair num valor padrão e gerar dado errado sem erro nenhum.
    """
    if not isinstance(dados, dict):
        return ["a raiz deve ser um objeto JSON"]
    esquema = dados.get("esquema")
    if esquema not in ESQUEMAS_ACEITOS:
        return [f"esquema {esquema!r} não reconhecido; esperado {ESQUEMA!r} ou {ESQUEMA_V2!r}"]

    versao = 2 if esquema == ESQUEMA_V2 else 1
    problemas = _campos(dados, _CAMPOS_RAIZ, set(), "raiz")

    gerado_por = dados.get("gerado_por")
    if isinstance(gerado_por, dict):
        problemas += _campos(gerado_por, _CAMPOS_GERADO_POR, set(), "gerado_por")
        for campo in sorted(_CAMPOS_GERADO_POR & gerado_por.keys()):
            if not isinstance(gerado_por[campo], str):
                problemas.append(f"gerado_por.{campo} deve ser texto")
    elif "gerado_por" in dados:
        problemas.append("gerado_por deve ser um objeto")

    conjuntos = dados.get("conjuntos")
    if not isinstance(conjuntos, list) or not conjuntos:
        if "conjuntos" in dados:
            problemas.append("conjuntos deve ser uma lista não vazia")
        return problemas

    vistos: set[str] = set()
    for i, conjunto in enumerate(conjuntos):
        problemas += _validar_conjunto(conjunto, f"conjuntos[{i}]", versao=versao)
        cid = conjunto.get("conjunto_id") if isinstance(conjunto, dict) else None
        if isinstance(cid, str):
            if cid in vistos:
                problemas.append(f"conjuntos[{i}]: conjunto {cid!r} repetido")
            vistos.add(cid)
    return problemas


def ler_selecao(caminho: str | Path) -> dict:
    """Lê e confere um selecao.json. Levanta ValueError listando todos os problemas."""
    dados = json.loads(Path(caminho).read_text(encoding="utf-8"))
    problemas = validar_selecao(dados)
    if problemas:
        esquema = dados.get("esquema") if isinstance(dados, dict) else None
        nome = esquema if esquema in ESQUEMAS_ACEITOS else ESQUEMA
        raise ValueError(f"{caminho} não segue {nome}:\n- " + "\n- ".join(problemas))
    return dados


def _campos(obj: dict, obrigatorios: set[str], opcionais: set[str], onde: str) -> list[str]:
    faltando = sorted(obrigatorios - obj.keys())
    desconhecidos = sorted(obj.keys() - obrigatorios - opcionais)
    return [f"{onde}: falta o campo {c!r}" for c in faltando] + [
        f"{onde}: campo desconhecido {c!r}" for c in desconhecidos
    ]


def _contagem(valor: object) -> bool:
    """Inteiro >= 0, ou null (ainda não medido). bool fica de fora: em Python, True é um int."""
    return valor is None or (isinstance(valor, int) and not isinstance(valor, bool) and valor >= 0)


def _fracao(valor: object) -> bool:
    """Número de 0 a 1, ou null. NaN não passa, porque nenhuma comparação com ele é verdadeira."""
    if valor is None:
        return True
    return isinstance(valor, int | float) and not isinstance(valor, bool) and 0 <= valor <= 1


def _validar_conjunto(conjunto: object, onde: str, versao: int = 1) -> list[str]:
    if not isinstance(conjunto, dict):
        return [f"{onde} deve ser um objeto"]
    problemas = _campos(conjunto, _CAMPOS_CONJUNTO, _CAMPOS_CONJUNTO_OPCIONAIS, onde)

    cid = conjunto.get("conjunto_id")
    if "conjunto_id" in conjunto:
        if versao == 1:
            if not (isinstance(cid, str) and "." in cid):
                problemas.append(f"{onde}.conjunto_id {cid!r} deve ser texto com ponto, como 'fosA.1'")
        else:
            if not (isinstance(cid, str) and bool(cid.strip())):
                problemas.append(f"{onde}.conjunto_id {cid!r} deve ser texto não vazio")

    origem = conjunto.get("origem")
    if "origem" in conjunto and versao == 1 and origem not in ORIGENS:
        problemas.append(f"{onde}.origem {origem!r} deve ser 'artigo' ou 'desenho'")
    elif "origem" in conjunto and versao == 2 and origem not in ORIGENS_V2:
        problemas.append(f"{onde}.origem {origem!r} deve ser 'artigo', 'desenho' ou 'usuario'")
    aviso = conjunto.get("aviso")
    if "aviso" in conjunto and not isinstance(aviso, str):
        problemas.append(f"{onde}.aviso deve ser texto")
    elif origem == "desenho" and not (isinstance(aviso, str) and aviso.strip()):
        problemas.append(f"{onde}: conjunto de desenho sem aviso (ADR-0013, regra 2)")

    if "primers" in conjunto:
        problemas += _validar_primers(conjunto["primers"], f"{onde}.primers", versao=versao)
    if "por_alelo" in conjunto:
        problemas += _validar_alelos(conjunto["por_alelo"], f"{onde}.por_alelo", versao=versao)
    if "selecionado" in conjunto and not isinstance(conjunto["selecionado"], bool):
        problemas.append(f"{onde}.selecionado deve ser true ou false")
    if "motivo" in conjunto and not isinstance(conjunto["motivo"], str):
        problemas.append(f"{onde}.motivo deve ser texto")
    return problemas


def _validar_primers(primers: object, onde: str, versao: int = 1) -> list[str]:
    if not isinstance(primers, dict):
        return [f"{onde} deve ser um objeto"]

    if versao == 1:
        conhecidos = set(PRIMERS_OBRIGATORIOS) | set(PRIMERS_OPCIONAIS)
        problemas = _campos(primers, set(PRIMERS_OBRIGATORIOS), set(PRIMERS_OPCIONAIS), onde)
        for nome in sorted(primers.keys() & conhecidos):
            seq = primers[nome]
            if not (isinstance(seq, str) and seq and set(seq) <= set("ACGT")):
                problemas.append(f"{onde}.{nome} deve ser uma sequência só de A, C, G e T")
        return problemas

    conhecidos = set(PRIMERS_OBRIGATORIOS) | set(PRIMERS_PARTES_V2) | set(PRIMERS_OPCIONAIS_V2)
    desconhecidos = sorted(primers.keys() - conhecidos)
    problemas = [f"{onde}: campo desconhecido {c!r}" for c in desconhecidos]

    tem_tubo = {"FIP", "BIP"}.issubset(primers.keys())
    tem_partes = {"F2", "F1c", "B2", "B1c"}.issubset(primers.keys())

    faltando: set[str] = set()
    for obr in ("F3", "B3"):
        if obr not in primers:
            faltando.add(obr)

    if not tem_tubo and not tem_partes:
        tem_algum_tubo = bool({"FIP", "BIP"} & primers.keys())
        tem_alguma_parte = bool({"F2", "F1c", "B2", "B1c"} & primers.keys())
        if tem_alguma_parte and not tem_algum_tubo:
            faltando.update({"F2", "F1c", "B2", "B1c"} - primers.keys())
        else:
            faltando.update({"FIP", "BIP"} - primers.keys())

    for c in sorted(faltando):
        problemas.append(f"{onde}: falta o campo {c!r}")

    for nome in sorted(primers.keys() & conhecidos):
        seq = primers[nome]
        if nome == "ligante":
            valido = isinstance(seq, str) and set(seq.strip().upper()) <= _BASES_V2
        else:
            valido = isinstance(seq, str) and bool(seq.strip()) and set(seq.strip().upper()) <= _BASES_V2
        if not valido:
            problemas.append(
                f"{onde}.{nome} deve ser uma sequência IUPAC válida (A, C, G, T, U ou códigos IUPAC)"
            )
    return problemas


def _validar_alelos(por_alelo: object, onde: str, versao: int = 1) -> list[str]:
    if not isinstance(por_alelo, list) or not por_alelo:
        return [f"{onde} deve ser uma lista não vazia"]
    problemas: list[str] = []
    vistos: set[str] = set()
    for i, medida in enumerate(por_alelo):
        aqui = f"{onde}[{i}]"
        if not isinstance(medida, dict):
            problemas.append(f"{aqui} deve ser um objeto")
            continue
        opcionais = _CAMPOS_ALELO_OPCIONAIS if versao == 1 else _CAMPOS_ALELO_OPCIONAIS_V2
        problemas += _campos(medida, _CAMPOS_ALELO, opcionais, aqui)

        aid = medida.get("alelo_id")
        if "alelo_id" in medida:
            if versao == 1:
                if not (isinstance(aid, str) and aid and "." not in aid):
                    problemas.append(f"{aqui}.alelo_id {aid!r} deve ser texto sem ponto, como 'fosA3'")
                elif isinstance(aid, str):
                    if aid in vistos:
                        problemas.append(f"{aqui}: alelo {aid!r} repetido no conjunto")
                    vistos.add(aid)
            else:
                if not (isinstance(aid, str) and bool(aid.strip())):
                    problemas.append(f"{aqui}.alelo_id {aid!r} deve ser texto não vazio")
                elif isinstance(aid, str):
                    if aid in vistos:
                        problemas.append(f"{aqui}: alelo {aid!r} repetido no conjunto")
                    vistos.add(aid)

        if "acesso" in medida and not (medida["acesso"] is None or isinstance(medida["acesso"], str)):
            problemas.append(f"{aqui}.acesso deve ser texto ou null")
        for campo in ("incompat_total", "incompat_3p"):
            if campo in medida and not _contagem(medida[campo]):
                problemas.append(f"{aqui}.{campo} deve ser inteiro >= 0 ou null")
        total, ponta = medida.get("incompat_total"), medida.get("incompat_3p")
        if _contagem(total) and _contagem(ponta) and total is not None and ponta is not None:
            if ponta > total:
                problemas.append(f"{aqui}: incompat_3p ({ponta}) maior que incompat_total ({total})")
        if "f_previsto" in medida and not _fracao(medida["f_previsto"]):
            problemas.append(f"{aqui}.f_previsto deve ser um número de 0 a 1, ou null")
        if versao == 2 and "incompat_regioes" in medida:
            problemas += _validar_regioes(medida["incompat_regioes"], aqui, total, ponta)
    return problemas


def _validar_regioes(regioes: object, aqui: str, total: object, ponta: object) -> list[str]:
    """Confere incompat_regioes (v2) e a coerência com incompat_total e incompat_3p."""
    onde = f"{aqui}.incompat_regioes"
    if not isinstance(regioes, list) or not regioes:
        return [f"{onde} deve ser uma lista não vazia"]
    problemas: list[str] = []
    vistas: set[str] = set()
    for j, item in enumerate(regioes):
        aqui_j = f"{onde}[{j}]"
        if not isinstance(item, dict):
            problemas.append(f"{aqui_j} deve ser um objeto")
            continue
        problemas += _campos(item, _CAMPOS_REGIAO_V2, set(), aqui_j)
        reg = item.get("regiao")
        if "regiao" in item:
            if reg not in REGIOES_INCOMPAT:
                problemas.append(f"{aqui_j}.regiao {reg!r} deve ser uma de {list(REGIOES_INCOMPAT)}")
            elif reg in vistas:
                problemas.append(f"{aqui_j}: região {reg!r} repetida")
            else:
                vistas.add(reg)
        encontrada = item.get("encontrada")
        if "encontrada" in item and not isinstance(encontrada, bool):
            problemas.append(f"{aqui_j}.encontrada deve ser true ou false")
        for campo in ("total", "ponta"):
            if campo in item and not _contagem(item[campo]):
                problemas.append(f"{aqui_j}.{campo} deve ser inteiro >= 0 ou null")
        t_reg, p_reg = item.get("total"), item.get("ponta")
        if encontrada is True and (t_reg is None or p_reg is None):
            problemas.append(f"{aqui_j}: região encontrada precisa de total e ponta")
        if encontrada is False and (t_reg is not None or p_reg is not None):
            problemas.append(f"{aqui_j}: região não encontrada leva total e ponta null")
        if _contagem(t_reg) and _contagem(p_reg) and t_reg is not None and p_reg is not None:
            if p_reg > t_reg:
                problemas.append(f"{aqui_j}: ponta ({p_reg}) maior que total ({t_reg})")
    if problemas:
        return problemas

    if all(item["encontrada"] for item in regioes):
        soma_total = sum(item["total"] for item in regioes)
        soma_ponta = sum(item["ponta"] for item in regioes)
        if total != soma_total:
            problemas.append(f"{aqui}: incompat_total ({total}) difere da soma das regiões ({soma_total})")
        if ponta != soma_ponta:
            problemas.append(f"{aqui}: incompat_3p ({ponta}) difere da soma das regiões ({soma_ponta})")
    elif total is not None or ponta is not None:
        problemas.append(f"{aqui}: com região não encontrada, incompat_total e incompat_3p devem ser null")
    return problemas
