# Como rodar o lamp-sim: guia para um Claude

Guia para um Claude (ou qualquer agente) que nunca trabalhou neste projeto. Ele ensina a
instalar, testar, rodar e mexer no código sem quebrar as regras do projeto. Tudo aqui foi
conferido em 18/10/2026 no Windows, com Python 3.14, no branch `feat/esteira`.

---

## 1. O que é o projeto

Simulador do Projeto LAMP (SENAC · FMABC). São quatro pacotes em `src/`:

| Pacote | O que faz | Autor |
|---|---|---|
| `primers` | Analisa um conjunto de primers LAMP contra um alvo e desenha conjuntos novos | Diego Ragozini |
| `lampsim` | Simula a reação colorimétrica: pH, cor do vermelho de fenol, classes, lotes | Henrique Prado |
| `esteira` | Liga os dois pelo contrato `selecao.json` (v2) e pela ponte análise → simulação | a dupla |
| `interface` | Site local (FastAPI) com Análise, Desenho, Parâmetros e Tubos | — |

Os números ficam em `config/` (`primers.toml`, `esteira.toml`, `web.toml`) e em
`src/lampsim/config.py`. **Todo resultado é provisório até a validação em bancada.**

---

## 2. Regras que você NÃO pode quebrar

1. **Nenhuma sequência sai da máquina (ADR-0013, regra 3).** O servidor escuta só em
   `127.0.0.1`. Nunca troque para `0.0.0.0`, nunca suba o site num serviço externo e nunca
   envie um FASTA para uma API.
2. **Sequências de isolados clínicos nunca entram no Git.** A pasta `dados/` está no
   `.gitignore`. Os testes usam só as sequências sintéticas de `tests/dados/`.
3. **`primers` e `lampsim` nunca se importam** (`tests/test_fronteiras.py`). A única ligação
   entre eles é a `esteira`.
4. **Nenhum número solto no código.** Os parâmetros vão no `config/`, e o código lê com um
   fallback via `.get(...)`.
5. **Aleatoriedade sempre com semente.** Mesma entrada e mesma semente dão o mesmo resultado.
6. **Versione a receita, não o resultado (LAMP-GIT-002).** `saida/`, `datasets/`, imagens,
   `.npy` e `_dados_local/` ficam fora do Git.
7. **Nada vai direto para o `main` (ADR-0010).** O gancho `.githooks/pre-push` recusa. Use um
   branch e abra um pedido de integração. **Nunca** use `--no-verify`.
8. **Commit, push e merge só com o "sim" do usuário.**

---

## 3. Instalar

Requer **Python 3.12 ou mais novo** (numpy 2.5.2 e scipy 1.18.1 não instalam em 3.11).

```bash
cd "C:/Users/ragoz/Documents/lamp project/lamp-sim-esteira"
python --version                      # tem de ser >= 3.12
pip install -r requirements.txt       # todas as versões são fixas (==)
git config core.hooksPath .githooks   # ativa o gancho que protege o main
```

---

## 4. Conferir que está tudo verde

Rode isto sempre antes e depois de mexer:

```bash
python -m ruff check .                # estilo; tem de dar "All checks passed!"
python -m pytest -q -p no:cacheprovider
```

Resultado esperado hoje: **354 testes passando** (~35 s) e o ruff limpo.

- O `pythonpath = ["src"]` do `pyproject.toml` já resolve os imports nos testes.
- Sem `primer3-py` e `biopython` (no Colab, por exemplo), os testes do `primers` ficam fora da
  coleta e os do `lampsim` continuam rodando (`tests/conftest.py`). Isso não é falha.
- Se um teste antigo quebrar, descubra antes se o comportamento mudou de propósito. Não apague
  nem enfraqueça o teste para fazê-lo passar.

---

## 5. Rodar o site (interface web)

No PowerShell:

```powershell
$env:PYTHONPATH = "src"
python -m interface
```

No Git Bash:

```bash
PYTHONPATH=src python -m interface
```

- O site abre em **http://127.0.0.1:8000**. A porta está em `config/web.toml`, seção `[servidor]`.
- O comando ocupa o terminal. Rode-o em segundo plano (ou numa aba de terminal própria) e pare
  com Ctrl+C.
- O estado fica só na memória: ao reiniciar, o histórico some.
- Para conferir sem navegador: `curl -s -o /dev/null -w "%{http_code}" http://127.0.0.1:8000/`
  precisa devolver `200`.
- Os exemplos da página Início:
  - `conjunto_ok — sintético`: usa `tests/dados/`, sempre disponível.
  - `fosA.1 × fosA6`: usa `dados/`, que existe só nesta máquina e não está no Git.

---

## 6. Rodar pela linha de comando

Rode da raiz, com `PYTHONPATH=src`. **A pasta de saída tem de ficar dentro de `saida/`**: fora
dela, o comando sai com o erro "deve ficar dentro de 'saida/'".

```bash
# analisar um conjunto existente
PYTHONPATH=src python -m primers analisar --entrada tests/dados/conjunto_ok.csv \
  --alvo tests/dados/alvo_ok.fasta --saida saida/analise

# desenhar conjuntos novos (alvo de até 2.000 pb na interface)
PYTHONPATH=src python -m primers desenhar --alvo tests/dados/alvo_ok.fasta --n 5 --saida saida/desenho
# --relaxado tenta com as faixas de [desenho_relaxado] quando o desenho sai vazio
```

A análise grava `primers.csv`, `violacoes.csv`, `heterodimeros.csv` e `resultado.json`.
No exemplo sintético, o esperado é: *"8 primers, 8 violacoes, 0 avisos"*. As coordenadas dos
arquivos exportados estão em **base 0 semiaberta**; a tela mostra base 1.

A `esteira` ainda não tem linha de comando própria: use-a como biblioteca
(`esteira.ponte`, `esteira.contrato`) ou pela página Tubos.

---

## 7. Simulador no Google Colab

1. Coloque a pasta do projeto em `Meu Drive/LAMP/lamp-sim`.
2. Abra `notebooks/00_setup.ipynb` no Colab e execute tudo. No fim deve aparecer "Estilo e
   testes: tudo verde".
3. Depois abra `notebooks/01_simulador.ipynb` e execute tudo. Os dados vão para
   `Meu Drive/LAMP/lamp-sim-dados`.

Fora do Colab, os notebooks gravam em `_dados_local/` (ignorado pelo Git).

---

## 8. Onde está cada coisa

| Precisa de | Onde |
|---|---|
| Regras de primer (Tm, GC, distâncias, ponta, desenho) | `config/primers.toml` |
| Servidor, limites de envio, exemplos, tubos | `config/web.toml` |
| Ponte primers → simulação (penalidades, eficiência) | `config/esteira.toml` |
| Números do simulador (cinética, cor) | `src/lampsim/config.py` |
| O que já foi corrigido | `MELHORIAS.md` |
| O que está em aberto | `PENDENCIAS.md` |
| Revisão original com os 88 problemas, specs e diffs das fases | `../correcao-revisao/` (fora do repositório) |
| Manual de uso para pessoas | Claude Docs: "Manual do LAMP Primers" |

---

## 9. Fluxo para mudar o código

1. `git status` limpo e branch próprio: `git switch -c tipo/descricao-curta` (`feat/`, `fix/`,
   `docs/`, `chore/`).
2. Leia o trecho de código e o config que a mudança toca. Siga o estilo do arquivo: português,
   linhas de até 110 caracteres, nenhuma variável chamada `l`.
3. Escreva ou ajuste os testes com sequências sintéticas e valores esperados calculados.
4. Rode `ruff` e `pytest` (seção 4). Se mexeu na interface, suba o site e confira a página.
5. Mostre o diff ao usuário. Faça o commit só com o "sim" dele, terminando a mensagem com:
   `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
6. Remotos: `origin` é o repositório da organização (`LAMP-Senac-FMABC/lamp-sim`); `pessoal`
   é o do Diego (`DiegoRag/LAMP-simulator-`, privado). Pergunte para qual remoto enviar.

Para delegar código ao Antigravity (`agy`), use a skill `/delegar` e o
`guia-delegacao-agentes.md`.

---

## 10. Problemas comuns

| Sintoma | Causa e solução |
|---|---|
| `pip install` falha no numpy ou no scipy | Python 3.11 ou mais antigo: use 3.12 ou mais novo |
| `ModuleNotFoundError: primers` | Faltou `PYTHONPATH=src` (ou `pip install -e .`) |
| "pasta de saida … deve ficar dentro de 'saida/'" | Use `--saida saida/<nome>` |
| Exemplo fosA não aparece no Início | A pasta `dados/` não existe nesta máquina; é normal |
| Push recusado: "enviar commits direto para o main" | Gancho da ADR-0010: envie um branch |
| Porta 8000 ocupada | Outro servidor ainda rodando: pare-o, ou mude a `porta` no `web.toml` |
| Desenho recusa o alvo | Acima de 2.000 pb (`alvo_desenho_max_pb`); a análise aceita genomas inteiros |
| Testes do `primers` "somem" | Faltam primer3-py e biopython: `pip install -r requirements.txt` |
