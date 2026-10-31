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
from starlette.concurrency import run_in_threadpool
from starlette.datastructures import UploadFile as StarletteUploadFile
from starlette.exceptions import HTTPException as StarletteHTTPException
from starlette.types import Scope

import esteira.ponte
import primers.desenho
import primers.entrada
import primers.saida
from interface.analise_isolada import AnalisadorIsolado, tempo_max_s
from interface.sessao import Execucao, Sessao
from interface.sessoes import Sessoes, config_sessoes
from interface.telas import (
    tela_base,
    tela_desenho,
    tela_inicio,
    tela_tubos,
)
from interface.telas_fluxo import (
    ABAS_DETALHE,
    ABAS_PARAMETROS,
    ERRO_SEM_ARQUIVO,
    tela_desenho_vazio,
    tela_detalhe,
    tela_detalhe_vazia,
    tela_landing,
    tela_marca,
    tela_parametros_marca,
    tela_passo_alvo,
    tela_passo_primers,
    tela_resultado,
)

ABAS_EXTRAS = {chave for chave, _, sufixo in ABAS_DETALHE if sufixo}
ABAS_PARAMETROS_EXTRAS = {chave for chave, _, sufixo in ABAS_PARAMETROS if sufixo}


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

    versao = _obter_versao(raiz)

    app = FastAPI(title="LampSim", docs_url=None, redoc_url=None, openapi_url=None)

    cfg_sessoes = config_sessoes(cfg_web)
    sessoes = Sessoes(
        lambda: Sessao(raiz, cfg_primers, cfg_web),
        cfg_sessoes["validade_min"],
        cfg_sessoes["max_sessoes"],
    )
    app.state.sessoes = sessoes
    analisador = AnalisadorIsolado(raiz, cfg_primers, cfg_web, tempo_max_s(cfg_web))
    app.state.analisador = analisador

    @app.middleware("http")
    async def sessao_do_navegador(request: Request, call_next):
        if request.url.path.startswith("/static/"):
            return await call_next(request)
        entrada, nova = sessoes.obter(request.cookies.get(cfg_sessoes["cookie"]))
        request.state.sessao = entrada.sessao
        request.state.token = entrada.token
        resposta = await call_next(request)
        if nova:
            resposta.set_cookie(
                cfg_sessoes["cookie"],
                entrada.token,
                max_age=int(cfg_sessoes["validade_min"] * 60),
                path="/",
                httponly=True,
                samesite="lax",
                secure=cfg_sessoes["cookie_seguro"],
            )
        return resposta

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
        sessao_erro = getattr(request.state, "sessao", None)
        if sessao_erro is None:
            sessao_erro = Sessao(raiz, cfg_primers, cfg_web)
        contexto = {
            "base": tela_base("erro", sessao_erro, versao),
            "tela": {"detalhe": str(exc.detail)},
            "mensagem": None,
        }
        return templates.TemplateResponse(request, "erro.html", contexto, status_code=404)

    def ir_para(request: Request, url: str) -> Response:
        """Redireciona com 303, ou devolve o destino em JSON para a tela de espera do navegador."""
        if request.headers.get("x-requested-with") == "lamp-carregamento":
            return JSONResponse({"destino": url})
        return RedirectResponse(url=url, status_code=303)

    def renderizar(request: Request, nome_template: str, pagina: str, tela: dict) -> Response:
        token = getattr(request.state, "token", None)
        mensagem = sessoes.tomar_mensagem(token) if token else None
        sessao_atual = getattr(request.state, "sessao", None)
        if sessao_atual is None:
            sessao_atual = Sessao(raiz, cfg_primers, cfg_web)
        contexto = {
            "base": tela_base(pagina, sessao_atual, versao),
            "tela": tela,
            "mensagem": mensagem,
        }
        return templates.TemplateResponse(request, nome_template, contexto)

    def renderizar_marca(
        request: Request,
        nome_template: str,
        pagina: str,
        tela: dict,
        erro: str | None = None,
        status: int = 200,
    ) -> Response:
        contexto = {
            "marca": tela_marca(pagina),
            "tela": tela,
            "erro": erro,
        }
        return templates.TemplateResponse(request, nome_template, contexto, status_code=status)

    async def ler_envio(arquivo: object) -> tuple[str, bytes] | None:
        if not isinstance(arquivo, StarletteUploadFile) or not arquivo.filename:
            return None
        tamanho_max_kb = int(cfg_web["entrada"]["tamanho_max_kb"])
        conteudo = await arquivo.read(tamanho_max_kb * 1024 + 1)
        nome = Path(arquivo.filename).name
        return (nome, conteudo)

    @app.get("/")
    def landing(request: Request) -> Response:
        return renderizar_marca(request, "landing.html", "landing", tela_landing(cfg_web))

    @app.get("/dados")
    def inicio(request: Request) -> Response:
        return renderizar(request, "inicio.html", "inicio", tela_inicio(request.state.sessao))

    @app.get("/testar")
    def passo_primers(request: Request) -> Response:
        tela = tela_passo_primers(request.state.sessao)
        return renderizar_marca(request, "passo_primers.html", "testar", tela, erro=None, status=200)

    @app.post("/testar/primers")
    async def testar_primers(request: Request) -> Response:
        sessao = request.state.sessao
        formulario = await request.form()
        arquivo = formulario.get("arquivo")
        envio = await ler_envio(arquivo)
        if envio is None:
            tela = tela_passo_primers(sessao)
            return renderizar_marca(
                request, "passo_primers.html", "testar", tela, erro=ERRO_SEM_ARQUIVO, status=400
            )
        nome, conteudo = envio
        try:
            await run_in_threadpool(sessao.carregar_conjunto, nome, conteudo)
        except (ValueError, RuntimeError) as e:
            tela = tela_passo_primers(sessao)
            return renderizar_marca(request, "passo_primers.html", "testar", tela, erro=str(e), status=400)
        return RedirectResponse(url="/testar/alvo", status_code=303)

    @app.get("/testar/alvo")
    def passo_alvo(request: Request) -> Response:
        sessao = request.state.sessao
        if sessao.conjunto is None:
            return RedirectResponse(url="/testar", status_code=303)
        tela = tela_passo_alvo(sessao)
        return renderizar_marca(request, "passo_alvo.html", "testar", tela, erro=None, status=200)

    @app.post("/testar/alvo")
    async def testar_alvo(request: Request) -> Response:
        sessao = request.state.sessao
        if sessao.conjunto is None:
            return RedirectResponse(url="/testar", status_code=303)
        formulario = await request.form()
        arquivo = formulario.get("arquivo")
        envio = await ler_envio(arquivo)
        if envio is None:
            tela = tela_passo_alvo(sessao)
            return renderizar_marca(
                request, "passo_alvo.html", "testar", tela, erro=ERRO_SEM_ARQUIVO, status=400
            )
        nome, conteudo = envio
        try:
            await run_in_threadpool(sessao.carregar_alvo, nome, conteudo)
            conjunto = sessao.conjunto
            execucao = await run_in_threadpool(
                app.state.analisador.analisar, conjunto.nome, conjunto.conteudo, nome, conteudo
            )
            execucao = sessao.registrar(execucao)
        except (ValueError, RuntimeError) as e:
            tela = tela_passo_alvo(sessao)
            return renderizar_marca(request, "passo_alvo.html", "testar", tela, erro=str(e), status=400)
        return RedirectResponse(url=f"/testar/resultado/{execucao.id}", status_code=303)

    @app.get("/testar/resultado/{id_execucao}")
    def testar_resultado(request: Request, id_execucao: int) -> Response:
        sessao = request.state.sessao
        try:
            execucao = sessao.execucao(id_execucao)
        except KeyError:
            raise HTTPException(status_code=404, detail="Execução não encontrada") from None
        if execucao.tipo != "analise":
            raise HTTPException(status_code=404, detail="Execução não encontrada")

        tubos_ctx = tela_tubos(sessao.cfg_web, cfg_esteira, execucao, "simulacao")
        tela = tela_resultado(execucao, tubos_ctx, cfg_primers, cfg_web)
        return renderizar_marca(request, "resultado.html", "resultado", tela, erro=None, status=200)

    @app.post("/conjunto")
    async def carregar_conjunto(request: Request, arquivo: Annotated[UploadFile, File()]) -> Response:
        sessao = request.state.sessao
        tamanho_max_kb = int(cfg_web["entrada"]["tamanho_max_kb"])
        conteudo = await arquivo.read(tamanho_max_kb * 1024 + 1)
        nome = Path(arquivo.filename or "arquivo").name
        try:
            sessao.carregar_conjunto(nome, conteudo)
        except (ValueError, RuntimeError) as e:
            sessoes.definir_mensagem(request.state.token, {"tipo": "erro", "texto": str(e)})
        return ir_para(request, "/dados")

    @app.post("/alvo")
    async def carregar_alvo(request: Request, arquivo: Annotated[UploadFile, File()]) -> Response:
        sessao = request.state.sessao
        tamanho_max_kb = int(cfg_web["entrada"]["tamanho_max_kb"])
        conteudo = await arquivo.read(tamanho_max_kb * 1024 + 1)
        nome = Path(arquivo.filename or "arquivo").name
        try:
            sessao.carregar_alvo(nome, conteudo)
        except (ValueError, RuntimeError) as e:
            sessoes.definir_mensagem(request.state.token, {"tipo": "erro", "texto": str(e)})
        return ir_para(request, "/dados")

    @app.post("/exemplo")
    def carregar_exemplo(request: Request, nome: str = Form(...)) -> Response:
        sessao = request.state.sessao
        try:
            sessao.carregar_exemplo(nome)
        except (ValueError, RuntimeError) as e:
            sessoes.definir_mensagem(request.state.token, {"tipo": "erro", "texto": str(e)})
        return ir_para(request, "/dados")

    @app.post("/analisar")
    def analisar(request: Request) -> Response:
        sessao = request.state.sessao
        try:
            execucao = sessao.analisar()
            return ir_para(request, f"/analise/{execucao.id}")
        except (ValueError, RuntimeError) as e:
            sessoes.definir_mensagem(request.state.token, {"tipo": "erro", "texto": str(e)})
            return ir_para(request, "/dados")

    @app.get("/analise")
    def ver_ultima_analise(request: Request) -> Response:
        sessao = request.state.sessao
        ultima = sessao.ultima("analise")
        if ultima is not None:
            return RedirectResponse(url=f"/analise/{ultima.id}", status_code=303)
        return renderizar_marca(
            request, "analise_vazia.html", "analise", tela_detalhe_vazia(), status=200
        )

    @app.get("/analise/{id_execucao}")
    def ver_analise(
        request: Request,
        id_execucao: int,
        tm_calibracao: int | None = None,
    ) -> Response:
        sessao = request.state.sessao
        try:
            execucao = sessao.execucao(id_execucao)
        except KeyError:
            raise HTTPException(status_code=404, detail="Execução não encontrada") from None
        if execucao.tipo != "analise":
            raise HTTPException(status_code=404, detail="Execução não encontrada")

        tela = tela_detalhe(
            execucao,
            "resumo",
            cfg_primers,
            cfg_web,
            tm_calibracao=tm_calibracao,
        )
        return renderizar_marca(request, "analise.html", "analise", tela)

    @app.get("/analise/{id_execucao}/exportar")
    def exportar_analise(request: Request, id_execucao: int) -> Response:
        sessao = request.state.sessao
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

    @app.get("/analise/{id_execucao}/{aba}")
    def ver_analise_aba(
        request: Request,
        id_execucao: int,
        aba: str,
        tm_calibracao: int | None = None,
    ) -> Response:
        if aba not in ABAS_EXTRAS:
            raise HTTPException(status_code=404, detail="Aba não encontrada")
        sessao = request.state.sessao
        try:
            execucao = sessao.execucao(id_execucao)
        except KeyError:
            raise HTTPException(status_code=404, detail="Execução não encontrada") from None
        if execucao.tipo != "analise":
            raise HTTPException(status_code=404, detail="Execução não encontrada")

        tela = tela_detalhe(
            execucao,
            aba,
            cfg_primers,
            cfg_web,
            tm_calibracao=tm_calibracao,
        )
        return renderizar_marca(request, "analise.html", "analise", tela)

    @app.post("/desenhar")
    def desenhar(request: Request, relaxado: int = Form(0)) -> Response:
        sessao = request.state.sessao
        try:
            execucao = sessao.desenhar(relaxado=bool(relaxado))
            return ir_para(request, f"/desenho/{execucao.id}")
        except (ValueError, RuntimeError) as e:
            sessoes.definir_mensagem(request.state.token, {"tipo": "erro", "texto": str(e)})
            return ir_para(request, "/dados")

    @app.get("/desenho")
    def ver_ultimo_desenho(request: Request) -> Response:
        sessao = request.state.sessao
        ultima = sessao.ultima("desenho")
        if ultima is not None:
            return RedirectResponse(url=f"/desenho/{ultima.id}", status_code=303)
        return renderizar_marca(
            request, "analise_vazia.html", "desenho", tela_desenho_vazio(), status=200
        )

    @app.get("/desenho/{id_execucao}")
    def ver_desenho(
        request: Request,
        id_execucao: int,
        sel: int = 1,
        agrupar: int = 0,
    ) -> Response:
        sessao = request.state.sessao
        try:
            execucao = sessao.execucao(id_execucao)
        except KeyError:
            raise HTTPException(status_code=404, detail="Execução não encontrada") from None
        if execucao.tipo != "desenho":
            raise HTTPException(status_code=404, detail="Execução não encontrada")

        tm_a_calibrar = bool(cfg_web["analise"]["tm_como_calibracao"])
        tela = tela_desenho(execucao, cfg_primers, sel, bool(agrupar), tm_a_calibrar=tm_a_calibrar)
        analises_anteriores = [
            ex for ex in sessao.execucoes
            if ex.tipo == "analise" and ex.id < execucao.id
        ]
        tela["url_analise"] = f"/analise/{analises_anteriores[-1].id}" if analises_anteriores else None
        return renderizar_marca(request, "desenho.html", "desenho", tela)

    @app.get("/desenho/{id_execucao}/exportar")
    def exportar_desenho(request: Request, id_execucao: int) -> Response:
        sessao = request.state.sessao
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
        sessao = request.state.sessao
        try:
            execucao = sessao.analisar_proposta(id_execucao, posicao)
        except KeyError:
            raise HTTPException(status_code=404, detail="Execução não encontrada") from None
        except (ValueError, RuntimeError) as e:
            sessoes.definir_mensagem(request.state.token, {"tipo": "erro", "texto": str(e)})
            return ir_para(request, "/dados")
        return ir_para(request, f"/analise/{execucao.id}")

    @app.get("/parametros")
    def parametros(request: Request) -> Response:
        tela = tela_parametros_marca(request.state.sessao, cfg_primers, "calibracao")
        return renderizar_marca(request, "parametros.html", "parametros", tela)

    @app.get("/parametros/{aba}")
    def parametros_aba(request: Request, aba: str) -> Response:
        if aba not in ABAS_PARAMETROS_EXTRAS:
            raise HTTPException(status_code=404, detail="Aba não encontrada")
        tela = tela_parametros_marca(request.state.sessao, cfg_primers, aba)
        return renderizar_marca(request, "parametros.html", "parametros", tela)

    @app.get("/tubos")
    def tubos(request: Request, execucao: int | None = None, modo: str = "simulacao") -> Response:
        sessao = request.state.sessao
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
