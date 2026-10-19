"""Interface web local da ferramenta de primers LAMP.

Escuta só no endereço configurado em config/web.toml (127.0.0.1); todo o estado fica
em memória e nenhuma sequência sai da máquina (ADR-0013, regra 3).
"""

import io
import json
import subprocess
import tempfile
import tomllib
import zipfile
from pathlib import Path
from typing import Annotated

from fastapi import FastAPI, File, Form, HTTPException, Request, UploadFile
from fastapi.exception_handlers import http_exception_handler
from fastapi.responses import JSONResponse, RedirectResponse, Response
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from starlette.exceptions import HTTPException as StarletteHTTPException
from starlette.types import Scope

import esteira.ponte
import primers.desenho
import primers.entrada
import primers.saida
from interface.sessao import Execucao, Sessao
from interface.telas import (
    tela_analise,
    tela_base,
    tela_desenho,
    tela_inicio,
    tela_parametros,
    tela_tubos,
)


def _obter_versao(raiz: Path) -> tuple[str, str]:
    """Obtém ramo e commit do Git ou valores padrão em caso de falha."""
    try:
        res_ramo = subprocess.run(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"],
            cwd=raiz,
            capture_output=True,
            text=True,
            check=False,
        )
        res_commit = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=raiz,
            capture_output=True,
            text=True,
            check=False,
        )
        if res_ramo.returncode == 0 and res_commit.returncode == 0:
            ramo = res_ramo.stdout.strip()
            commit = res_commit.stdout.strip()
            if ramo and commit:
                return (ramo, commit)
    except Exception:
        pass
    return ("desconhecido", "—")


def descrever_alvo_da_execucao(execucao: Execucao, cfg_primers: dict) -> dict | None:
    """O alvo da execução para o resultado.json (primers.saida.descrever_alvo), ou None."""
    if execucao.alvo is None:
        return None
    partes = execucao.alvo.descricao.split()
    return primers.saida.descrever_alvo(
        execucao.alvo.sequencia,
        id_registro=partes[0] if partes else None,
        arquivo=execucao.alvo.nome,
        circular=bool(cfg_primers.get("alvo", {}).get("circular", False)),
    )


def selecao_da_execucao(execucao: Execucao, cfg_primers: dict, cfg_esteira: dict) -> dict | None:
    """selecao.json (selecao@2) de uma análise com conjunto completo e alvo; None sem os dois."""
    leitura = execucao.leitura
    if execucao.tipo != "analise" or execucao.alvo is None or leitura is None or leitura.conjunto is None:
        return None
    c = leitura.conjunto
    primers_sel = {
        "F3": c.F3,
        "B3": c.B3,
        "F2": c.F2,
        "F1c": c.F1c,
        "B2": c.B2,
        "B1c": c.B1c,
        "FIP": c.FIP,
        "BIP": c.BIP,
        "ligante": c.ligante,
    }
    if c.LF is not None:
        primers_sel["LF"] = c.LF
    if c.LB is not None:
        primers_sel["LB"] = c.LB
    partes = execucao.alvo.descricao.split()
    alelo_id = partes[0] if partes else (Path(execucao.alvo.nome).stem or "alvo")
    gerado_por = esteira.ponte.assinatura(
        primers.saida.versao_do_codigo(),
        f"interface/analise/{execucao.id}",
        {"primers": cfg_primers, "esteira": cfg_esteira},
    )
    return esteira.ponte.selecao_da_analise(
        conjunto_id=execucao.rotulo_conjunto,
        origem=execucao.origem,
        primers=primers_sel,
        incompat=execucao.resultado.incompat,
        alelo_id=alelo_id,
        acesso=esteira.ponte.acesso_do_id(alelo_id),
        cfg=cfg_esteira,
        gerado_por=gerado_por,
        aviso=primers.desenho.AVISO if execucao.origem == "desenho" else None,
    )


class EstaticosRevalidados(StaticFiles):
    """Arquivos estáticos que o navegador revalida a cada carga."""

    async def get_response(self, path: str, scope: Scope) -> Response:
        resposta = await super().get_response(path, scope)
        resposta.headers["Cache-Control"] = "no-cache"
        return resposta


def criar_app(raiz: Path | None = None) -> FastAPI:
    """Cria e configura o aplicativo FastAPI da interface web."""
    if raiz is None:
        raiz = Path(__file__).resolve().parents[2]

    cfg_primers = primers.entrada.carregar_config(raiz / "config" / "primers.toml")
    cfg_web = tomllib.loads((raiz / "config" / "web.toml").read_text(encoding="utf-8"))
    cfg_esteira = esteira.ponte.carregar_config(raiz / "config" / "esteira.toml")

    sessao = Sessao(raiz, cfg_primers, cfg_web)
    versao = _obter_versao(raiz)

    app = FastAPI(title="LAMP Primers", docs_url=None, redoc_url=None, openapi_url=None)
    app.state.sessao = sessao

    templates = Jinja2Templates(directory=str(Path(__file__).parent / "templates"))
    app.mount(
        "/static",
        EstaticosRevalidados(directory=str(Path(__file__).parent / "static")),
        name="static",
    )

    @app.exception_handler(StarletteHTTPException)
    async def pagina_de_erro(request: Request, exc: StarletteHTTPException) -> Response:
        """Mostra o 404 no visual do site; os outros códigos seguem o tratamento padrão."""
        if exc.status_code != 404:
            return await http_exception_handler(request, exc)
        contexto = {
            "base": tela_base("erro", sessao, versao),
            "tela": {"detalhe": str(exc.detail)},
            "mensagem": None,
        }
        return templates.TemplateResponse(request, "erro.html", contexto, status_code=404)

    mensagem_pendente: dict[str, dict[str, str] | None] = {"item": None}

    def ir_para(request: Request, url: str) -> Response:
        """Redireciona com 303, ou devolve o destino em JSON para a tela de espera do navegador."""
        if request.headers.get("x-requested-with") == "lamp-carregamento":
            return JSONResponse({"destino": url})
        return RedirectResponse(url=url, status_code=303)

    def renderizar(request: Request, nome_template: str, pagina: str, tela: dict) -> Response:
        mensagem = mensagem_pendente["item"]
        mensagem_pendente["item"] = None
        contexto = {
            "base": tela_base(pagina, sessao, versao),
            "tela": tela,
            "mensagem": mensagem,
        }
        return templates.TemplateResponse(request, nome_template, contexto)

    @app.get("/")
    def inicio(request: Request) -> Response:
        return renderizar(request, "inicio.html", "inicio", tela_inicio(sessao))

    @app.post("/conjunto")
    async def carregar_conjunto(request: Request, arquivo: Annotated[UploadFile, File()]) -> Response:
        tamanho_max_kb = int(cfg_web["entrada"]["tamanho_max_kb"])
        conteudo = await arquivo.read(tamanho_max_kb * 1024 + 1)
        nome = Path(arquivo.filename or "arquivo").name
        try:
            sessao.carregar_conjunto(nome, conteudo)
        except (ValueError, RuntimeError) as e:
            mensagem_pendente["item"] = {"tipo": "erro", "texto": str(e)}
        return ir_para(request, "/")

    @app.post("/alvo")
    async def carregar_alvo(request: Request, arquivo: Annotated[UploadFile, File()]) -> Response:
        tamanho_max_kb = int(cfg_web["entrada"]["tamanho_max_kb"])
        conteudo = await arquivo.read(tamanho_max_kb * 1024 + 1)
        nome = Path(arquivo.filename or "arquivo").name
        try:
            sessao.carregar_alvo(nome, conteudo)
        except (ValueError, RuntimeError) as e:
            mensagem_pendente["item"] = {"tipo": "erro", "texto": str(e)}
        return ir_para(request, "/")

    @app.post("/exemplo")
    def carregar_exemplo(request: Request, nome: str = Form(...)) -> Response:
        try:
            sessao.carregar_exemplo(nome)
        except (ValueError, RuntimeError) as e:
            mensagem_pendente["item"] = {"tipo": "erro", "texto": str(e)}
        return ir_para(request, "/")

    @app.post("/analisar")
    def analisar(request: Request) -> Response:
        try:
            execucao = sessao.analisar()
            return ir_para(request, f"/analise/{execucao.id}")
        except (ValueError, RuntimeError) as e:
            mensagem_pendente["item"] = {"tipo": "erro", "texto": str(e)}
            return ir_para(request, "/")

    @app.get("/analise")
    def ver_ultima_analise() -> RedirectResponse:
        ultima = sessao.ultima("analise")
        if ultima is not None:
            return RedirectResponse(url=f"/analise/{ultima.id}", status_code=303)
        return RedirectResponse(url="/", status_code=303)

    @app.get("/analise/{id_execucao}")
    def ver_analise(
        request: Request,
        id_execucao: int,
        tm_calibracao: int | None = None,
    ) -> Response:
        try:
            execucao = sessao.execucao(id_execucao)
        except KeyError:
            raise HTTPException(status_code=404, detail="Execução não encontrada") from None
        if execucao.tipo != "analise":
            raise HTTPException(status_code=404, detail="Execução não encontrada")

        if tm_calibracao is None:
            tm_como_calibracao = bool(cfg_web["analise"]["tm_como_calibracao"])
        else:
            tm_como_calibracao = bool(tm_calibracao)

        tela = tela_analise(
            execucao,
            cfg_primers,
            tm_como_calibracao,
            cfg_web["heterodimeros"]["escala_kcal"],
        )
        return renderizar(request, "analise.html", "analise", tela)

    @app.get("/analise/{id_execucao}/exportar")
    def exportar_analise(id_execucao: int) -> Response:
        try:
            execucao = sessao.execucao(id_execucao)
        except KeyError:
            raise HTTPException(status_code=404, detail="Execução não encontrada") from None
        if execucao.tipo != "analise":
            raise HTTPException(status_code=404, detail="Execução não encontrada")

        resultado = execucao.resultado
        eficiencia = None
        if resultado.incompat is not None:
            eficiencia = {
                "f_previsto": esteira.ponte.eficiencia_prevista(resultado.incompat, cfg_esteira),
                "modelo": "esteira.ponte.eficiencia_prevista com config/esteira.toml (provisório)",
                "parametros": dict(cfg_esteira["eficiencia"]),
            }
        selecao = selecao_da_execucao(execucao, cfg_primers, cfg_esteira)
        with tempfile.TemporaryDirectory() as pasta:
            pasta_p = Path(pasta)
            primers.saida.exportar(
                resultado,
                pasta_p,
                cfg_primers,
                leitura=execucao.leitura,
                alvo=descrever_alvo_da_execucao(execucao, cfg_primers),
                eficiencia=eficiencia,
            )
            nomes = ["primers.csv", "heterodimeros.csv", "violacoes.csv", "resultado.json"]
            if selecao is not None:
                texto = json.dumps(selecao, ensure_ascii=False, indent=2)
                (pasta_p / "selecao.json").write_text(texto, encoding="utf-8")
                nomes.append("selecao.json")
            buf = io.BytesIO()
            with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
                for nome in nomes:
                    zf.write(pasta_p / nome, arcname=nome)
            conteudo = buf.getvalue()

        return Response(
            content=conteudo,
            media_type="application/zip",
            headers={"Content-Disposition": f'attachment; filename="analise_{id_execucao}.zip"'},
        )

    @app.post("/desenhar")
    def desenhar(request: Request, relaxado: int = Form(0)) -> Response:
        try:
            execucao = sessao.desenhar(relaxado=bool(relaxado))
            return ir_para(request, f"/desenho/{execucao.id}")
        except (ValueError, RuntimeError) as e:
            mensagem_pendente["item"] = {"tipo": "erro", "texto": str(e)}
            return ir_para(request, "/")

    @app.get("/desenho")
    def ver_ultimo_desenho() -> RedirectResponse:
        ultima = sessao.ultima("desenho")
        if ultima is not None:
            return RedirectResponse(url=f"/desenho/{ultima.id}", status_code=303)
        return RedirectResponse(url="/", status_code=303)

    @app.get("/desenho/{id_execucao}")
    def ver_desenho(
        request: Request,
        id_execucao: int,
        sel: int = 1,
        agrupar: int = 0,
    ) -> Response:
        try:
            execucao = sessao.execucao(id_execucao)
        except KeyError:
            raise HTTPException(status_code=404, detail="Execução não encontrada") from None
        if execucao.tipo != "desenho":
            raise HTTPException(status_code=404, detail="Execução não encontrada")

        tm_a_calibrar = bool(cfg_web["analise"]["tm_como_calibracao"])
        tela = tela_desenho(execucao, cfg_primers, sel, bool(agrupar), tm_a_calibrar=tm_a_calibrar)
        return renderizar(request, "desenho.html", "desenho", tela)

    @app.get("/desenho/{id_execucao}/exportar")
    def exportar_desenho(id_execucao: int) -> Response:
        try:
            execucao = sessao.execucao(id_execucao)
        except KeyError:
            raise HTTPException(status_code=404, detail="Execução não encontrada") from None
        if execucao.tipo != "desenho":
            raise HTTPException(status_code=404, detail="Execução não encontrada")

        sequencia = execucao.alvo.sequencia if execucao.alvo is not None else ""
        with tempfile.TemporaryDirectory() as pasta:
            pasta_p = Path(pasta)
            primers.saida.exportar_desenho(execucao.resultado, sequencia, pasta_p, cfg_primers)
            buf = io.BytesIO()
            with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
                for nome in ("conjuntos.csv", "resultado.json"):
                    zf.write(pasta_p / nome, arcname=nome)
            conteudo = buf.getvalue()

        return Response(
            content=conteudo,
            media_type="application/zip",
            headers={"Content-Disposition": f'attachment; filename="desenho_{id_execucao}.zip"'},
        )

    @app.post("/desenho/{id_execucao}/analisar")
    def analisar_proposta(request: Request, id_execucao: int, posicao: int = Form(...)) -> Response:
        try:
            execucao = sessao.analisar_proposta(id_execucao, posicao)
        except KeyError:
            raise HTTPException(status_code=404, detail="Execução não encontrada") from None
        except (ValueError, RuntimeError) as e:
            mensagem_pendente["item"] = {"tipo": "erro", "texto": str(e)}
            return ir_para(request, "/")
        return ir_para(request, f"/analise/{execucao.id}")

    @app.get("/parametros")
    def parametros(request: Request) -> Response:
        return renderizar(request, "parametros.html", "parametros", tela_parametros(sessao, cfg_primers))

    @app.get("/tubos")
    def tubos(request: Request, execucao: int | None = None, modo: str = "simulacao") -> Response:
        if execucao is None:
            escolhida = sessao.ultima("analise")
        else:
            try:
                escolhida = sessao.execucao(execucao)
            except KeyError:
                raise HTTPException(status_code=404, detail="Execução não encontrada") from None
            if escolhida.tipo != "analise":
                raise HTTPException(status_code=404, detail="Execução não encontrada")
        tela = tela_tubos(sessao.cfg_web, cfg_esteira, escolhida, modo)
        return renderizar(request, "tubos.html", "tubos", tela)

    return app
