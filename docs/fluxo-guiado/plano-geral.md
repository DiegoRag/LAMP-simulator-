# Plano geral — fluxo guiado do LampSim (fase 1)

Repositório: `/home/diego/Documents/Lamp/lamp-sim`, branch `feat/fluxo-guiado` (a partir de
`feat/esteira`). Pacote afetado: `src/interface` (FastAPI + Jinja2 + JS sem framework).
Leia antes: `src/interface/app.py`, `src/interface/sessao.py`, `src/interface/telas.py`
(só `tela_base`, `tela_analise`, `resumo_violacoes`, `tela_tubos` e as funções que elas usam),
`src/interface/templates/base.html`, `templates/tubos.html`, `static/tubos.js`,
`config/web.toml`, `tests/test_interface_app.py`.

## 1. Objetivo

Trocar a porta de entrada do site por um fluxo linear, de tela cheia, para quem já sabe o que
quer: **landing → passo 1 (primers) → passo 2 (gene-alvo) → resultado (tubos)**. As telas novas
nascem na identidade visual LampSim (seção 6). As páginas antigas continuam existindo.

## 2. Papéis

| Papel | Modelo | Faz |
|---|---|---|
| Planejador | claude-opus-5-5-high | Escreve `fase-1/contrato.md`, `fase-1/s1.md … s9.md` e `fase-1/resumo.md` |
| Executores s1–s9 | claude-sonnet-5-5-medium (reserva gemini-3.8-flash-medium) | Um por spec, em `lamp-sim-w1 … w9`, arquivos disjuntos |
| Registrador (w10) | claude-sonnet-5-5-medium | Acrescenta a seção da fase ao `MELHORIAS.md` com base no diff real |
| Orquestrador | Claude Code | Sincroniza, roda testes e ruff, revisa diffs, integra |

## 3. Decisões adotadas (o executor não decide nada disto)

### 3.1 Rotas

| Método e rota | Comportamento |
|---|---|
| `GET /` | **Landing** nova (`landing.html`). Substitui a antiga página Início |
| `GET /dados` | A antiga página Início (`inicio.html`, `tela_inicio`), sem mudança de conteúdo. Fica fora da navbar; dá acesso ao Desenho, que volta depois |
| `GET /testar` | **Passo 1 de 2**: enviar o conjunto de primers (`passo_primers.html`) |
| `POST /testar/primers` | Recebe o arquivo (`arquivo`). Sucesso → 303 para `/testar/alvo`. Erro → **re-renderiza o passo 1 na mesma URL do POST com status 400** e a mensagem de erro em vermelho. Nunca redireciona em caso de erro |
| `GET /testar/alvo` | **Passo 2 de 2**: enviar o gene-alvo (`passo_alvo.html`). Sem conjunto carregado na sessão → 303 para `/testar` |
| `POST /testar/alvo` | Recebe o FASTA (`arquivo`), carrega o alvo e roda a análise (`Sessao.analisar()`). Sucesso → 303 para `/testar/resultado/{id}`. Erro no alvo **ou** na análise → re-renderiza o passo 2 com status 400 e a mensagem. Sem conjunto → 303 para `/testar` |
| `GET /testar/resultado/{id}` | **Resultado** (`resultado.html`): resumo + tubos simulados da análise `{id}` da sessão. Id inexistente ou que não é análise → 404 |

As rotas antigas continuam: `/conjunto`, `/alvo`, `/exemplo`, `/analisar`, `/analise…`,
`/desenho…`, `/desenhar`, `/parametros`, `/tubos`. Onde elas hoje voltam para `/` (sucesso de
envio e erros), passam a voltar para `/dados`.

### 3.2 Sessão por navegador (cookie)

Hoje existe **um** objeto `Sessao` para todo o servidor, e as mensagens também são globais. Num
site público, um visitante vê o estado do outro. Decisão:

- Módulo novo `src/interface/sessoes.py` com a classe `Sessoes`: um registro de `Sessao` por
  token. O token é `secrets.token_urlsafe(32)`.
- Cookie: nome, validade e limite vêm de `config/web.toml`, seção nova `[sessoes]`:
  `cookie = "lampsim_sessao"`, `validade_min = 120`, `max_sessoes = 500`,
  `cookie_seguro = false`. O código lê com `.get(...)` e o mesmo fallback.
- Cookie `HttpOnly`, `SameSite=Lax`, `Path=/`, `Max-Age = validade_min*60`; `Secure` só se
  `cookie_seguro` for verdadeiro.
- Token ausente, desconhecido ou expirado → nova sessão e novo cookie na resposta.
- Expiração: sessão sem uso há mais de `validade_min` é descartada. Ao passar de `max_sessoes`,
  descarta a menos usada recentemente. O relógio é injetável (`relogio=time.monotonic`) para
  testar sem esperar.
- As **mensagens pendentes** (`mensagem_pendente` em `app.py`) passam a ser por sessão.
- Thread-safe (`threading.Lock`).
- `app.py` usa um middleware HTTP que resolve a sessão, guarda em `request.state.sessao` e
  `request.state.token`, e grava o cookie quando a sessão é nova. Todas as rotas passam a usar a
  sessão da requisição (as que hoje não recebem `request` passam a receber). `app.state.sessoes`
  guarda o registro; `app.state.sessao` deixa de existir.
- `sessao.py` **não muda**.

### 3.3 Conteúdo das telas novas (textos em PT-BR, exatamente estes)

**Navbar das telas novas** (`base_marca.html`): logo horizontal colorida à esquerda (link para
`/`, com `alt="LampSim"`); à direita os links **Início** (`/`), **Testar primers** (`/testar`) e
**Parâmetros** (`/parametros`), com `aria-current="page"` na página atual. Não há trilho de
etapas, selo "Só nesta máquina" nem chave de animações nas telas novas.

**Landing** (`/`), uma tela cheia, direta:
- Logo completa (versão negativa) sobre fundo azul-marinho.
- Título: "Exames genéticos em qualquer lugar, a qualquer hora."
- Texto: "O LampSim substitui as máquinas LAMP convencionais por um aparelho muito menor,
  portátil e barato."
- Bloco curto "O que dá para fazer aqui": "Teste um conjunto de primers LAMP contra o gene-alvo e
  veja, em tubos simulados, se a reação vira de rosa (negativo) para amarelo (positivo)."
- Botão principal: "Testar primers" → `/testar`.
- Aviso discreto: "Resultados simulados e provisórios até a validação em bancada."
- Rodapé: "Projeto LAMP · Senac + FMABC".

**Passo 1** (`/testar`), tela cheia, conteúdo centralizado:
- Indicador "Passo 1 de 2".
- Título: "Envie o conjunto de primers".
- Explicação: "Um arquivo CSV com as colunas nome e sequencia, ou um FASTA com um registro por
  primer. Sequências sempre de 5′ para 3′." e "Nomes aceitos: F3, B3, F2, F1c, B2, B1c, FIP, BIP,
  LF e LB."
- Um exemplo curto de CSV num bloco `<pre>` (cabeçalho + 2 linhas de `tests/dados/conjunto_ok.csv`).
- **Um botão central grande**: "Escolher arquivo de primers". É um `<label>` estilizado para um
  `<input type="file" accept=".csv,.fasta,.fa,.fna">` visualmente oculto (mas focável). Ao
  escolher o arquivo, o formulário é enviado sozinho (JS). Não existe botão "Enviar" separado; um
  botão "Enviar" só aparece dentro de `<noscript>`.
- O botão fica dentro de uma moldura com os **cantos de visor** da logo (quatro cantos
  arredondados em azul-marinho, sem linha completa), o motivo da câmera do LampSim.
- Erro: mensagem em vermelho logo abaixo do botão, com `role="alert"`, ícone de alerta e o texto
  da exceção (já vem no formato "arquivo.csv: …").

**Passo 2** (`/testar/alvo`), mesmo layout:
- Indicador "Passo 2 de 2".
- Título: "Envie o gene-alvo".
- Explicação: "Um arquivo FASTA (.fasta, .fa ou .fna) com a sequência do gene que os primers devem
  reconhecer. Pode ter um ou vários registros; usamos o que mais combina com os primers."
- Resumo do passo 1: "Primers: {nome do arquivo} · {n} primers" e o link "Trocar primers" →
  `/testar`.
- Botão central: "Escolher arquivo do gene-alvo" (`accept=".fasta,.fa,.fna"`), com a mesma
  moldura de visor e o mesmo envio automático.
- Ao enviar, uma camada de espera cobre a tela: "Analisando os primers contra o gene-alvo…" (sem
  barra de progresso falsa; um indicador simples e o texto). Respeita `prefers-reduced-motion`.
- Erro: igual ao passo 1.

**Resultado** (`/testar/resultado/{id}`):
- Título "Resultado" com "{conjunto} × {alvo}".
- **Resumo** no topo, com um veredito em uma frase (regras na seção 3.4) e três números:
  "Regiões localizadas" (a de b), "Violações" (n, com "e m de Tm a calibrar" quando houver, a
  mesma conta do KPI da página de análise) e "Eficiência prevista" (texto de `texto_eficiencia`).
- Ações: "Ver análise detalhada" → `/analise/{id}`; "Exportar resultados" →
  `/analise/{id}/exportar`; "Testar outro conjunto" → `/testar`.
- Abaixo, os **tubos simulados**: bandeja, controles de tempo, faixas de cor, tabela "Leitura aos
  15 min" e gráfico de pH, reaproveitando `static/tubos.js` sem alterá-lo. Para isso o template usa
  exatamente os mesmos ganchos `data-*`, as classes que o JS consulta (`.tubo`, `.tubo-liquido`,
  `.tubos-faixa-linha`) e o `<script type="application/json" id="tubos-dados">` com o contexto de
  `tela_tubos`. Só o modo simulação: sem editor manual e sem botão "Usar curvas manuais".
- Se a simulação não for possível (`tela_tubos` devolveu `modo == "manual"`), mostrar o aviso de
  `tela_tubos` e o link para a análise detalhada, sem a bandeja.

### 3.4 Veredito do resumo (função pura, testável)

Entradas: a execução de análise e o contexto de `tela_tubos` (modo simulação).
Em ordem:
1. Eficiência prevista 0 → tom `"negativo"`: "Estes primers não devem amplificar este alvo."
   e, se houver regiões sem posição, "Regiões não encontradas no alvo: {lista}.".
2. Algum controle com diagnóstico `invalido` → tom `"alerta"`: "Simulação inválida: um controle
   não se comportou como esperado.".
3. Violações (sem as de Tm a calibrar) > 0 → tom `"alerta"`: "Estes primers devem amplificar este
   alvo, com ressalvas: {n} regra(s) do guia não cumprida(s).".
4. Senão → tom `"positivo"`: "Estes primers devem amplificar este alvo.".

Em todos os casos com simulação, acrescentar a frase de sensibilidade: "Na simulação, detectam
a partir de {menor nº de cópias de amostra com diagnóstico positivo aos 15 min} cópia(s)." ou,
se nenhuma amostra positivar, "Na simulação, nenhuma amostra positivou aos 15 min.".
Números formatados com `interface.formato.milhar`.

A contagem de violações usa `telas.resumo_violacoes(violacoes, tm_como_calibracao)` com
`tm_como_calibracao = cfg_web["analise"]["tm_como_calibracao"]` (chaves `n_violacoes`, `n_tm`,
`detalhe`). Os códigos de diagnóstico dos tubos são os de `telas.DIAGNOSTICO_TEXTO`
(`positivo`, `negativo`, `controle_ok`, `invalido`) e os papéis os de `telas.PAPEL_TEXTO`
(`amostra`, `controle_negativo`, `controle_positivo`).

### 3.5 Páginas antigas

- `base.html` (layout antigo): a navbar troca o trilho de etapas e os links "Tubos" e
  "Parâmetros" pelos mesmos três links das telas novas (Início `/`, Testar primers `/testar`,
  Parâmetros `/parametros`). Ficam a marca antiga, a chave de animações, o selo e a tela de
  espera, que as páginas antigas ainda usam. `telas.py` **não muda**.
- Os testes antigos que esperam o Início em `/`, o trilho ou redirecionamento para `/` são
  ajustados para `/dados` e para a navbar nova. Nenhum teste é apagado ou enfraquecido.

## 4. Divisão de arquivos (disjunta, um executor por linha)

| Spec | Pasta | PODE editar (só estes) |
|---|---|---|
| s1 | w1 | `src/interface/sessoes.py` (criar), `tests/test_fase1_s1.py` (criar) |
| s2 | w2 | `src/interface/app.py`, `config/web.toml` (só a seção `[sessoes]` nova), `tests/test_interface_app.py` (ajustes da seção 3.5), `tests/test_fase1_s2.py` (criar) |
| s3 | w3 | `src/interface/telas_fluxo.py` (criar), `tests/test_fase1_s3.py` (criar) |
| s4 | w4 | `src/interface/templates/base_marca.html` (criar), `src/interface/templates/landing.html` (criar) |
| s5 | w5 | `src/interface/templates/passo_primers.html` (criar), `src/interface/templates/passo_alvo.html` (criar) |
| s6 | w6 | `src/interface/templates/resultado.html` (criar) |
| s7 | w7 | `src/interface/static/marca.css` (criar) |
| s8 | w8 | `src/interface/static/marca-resultado.css` (criar), `src/interface/static/fluxo.js` (criar) |
| s9 | w9 | `src/interface/templates/base.html`, `tests/test_fase1_s9.py` (criar: teste de ponta a ponta do fluxo com `TestClient`, usando `tests/dados/conjunto_ok.csv` e `tests/dados/alvo_ok.fasta`) |

Já existem e **não são tocados por ninguém**: `static/marca/*.svg` (logos oficiais do kit:
`lampsim_horizontal_colorido.svg`, `lampsim_horizontal_negativo.svg`,
`lampsim_completo_colorido.svg`, `lampsim_completo_negativo.svg`,
`lampsim_simbolo_colorido.svg`, `lampsim_simbolo_negativo.svg`) e
`static/fontes/montserrat-latin.woff2` (Montserrat variável, pesos 400–700, subconjunto latino).

## 5. O contrato precisa fixar (fonte única para as sessões)

- Assinaturas e docstrings de `Sessoes` (s1) e de cada função de `telas_fluxo.py` (s3), com as
  chaves exatas dos dicionários de contexto que os templates usam.
- Para cada template: nome do arquivo, bloco(s) Jinja, o contexto recebido (`base`, `tela`,
  `erro`…) e **todas as classes CSS e atributos `data-*`** que ele usa. s7 e s8 estilizam só esse
  vocabulário; s4–s6 usam só esse vocabulário.
- O `base_marca.html` define os blocos que os outros templates estendem, os `<link>` de
  `marca.css` (sempre) e de `marca-resultado.css` (só no resultado, por bloco), o favicon
  (`/static/marca/lampsim_simbolo_colorido.svg`) e o `<script src="/static/fluxo.js" defer>`.
- O contrato de `fluxo.js`: quais seletores ele procura (envio automático ao escolher o arquivo,
  camada de espera no passo 2, bloqueio de duplo envio, reabrir o estado normal no `pageshow` do
  bfcache).
- As funções de `app.py` que os testes de s2 e s9 exercitam e o que cada rota devolve (status,
  `Location`, textos-chave que aparecem no HTML).

## 6. Identidade visual (o Manual de identidade visual LampSim manda)

- Cores: azul-marinho `#15344A` (base, textos, botões), rosa `#E83A73`, amarelo `#F4C44B`,
  branco `#FFFFFF`; fundo claro alternativo off-white `#F6F4EF`. Texto rosa pequeno sobre fundo
  claro usa `#C2275E`. Dosagem: branco ~55 %, azul ~30 %, rosa ~8 %, amarelo ~7 %; rosa e amarelo
  como destaque, de preferência juntos. **Nada de degradê** nos elementos de marca.
- Vermelho de erro (não é cor de marca, só para mensagens): texto `#B42318`, fundo `#FEF3F2`,
  borda esquerda `#B42318`. Nunca usar o rosa da marca para erro.
- Tipografia: Montserrat (`@font-face` local, `font-weight: 400 700`, `font-display: swap`).
  Bold 700 títulos; SemiBold 600 rótulos e botões; Regular 400 texto. Fallback
  `system-ui, sans-serif`.
- Logo: sempre os SVGs do kit, sem recolorir, distorcer, girar ou aplicar sombra. Horizontal
  colorida na navbar com **largura mínima de 160 px**; completa negativa na landing (fundo azul),
  largura mínima 120 px. Área de proteção livre em volta.
- As cores dos tubos e das faixas são científicas (vermelho de fenol) e continuam vindo de
  `tela_tubos`/`tubos.js`; não trocar pela paleta da marca.
- Acessibilidade: contraste AA, foco visível em tudo que é clicável, `prefers-reduced-motion`
  respeitado, rótulos reais nos inputs, `lang="pt-BR"`, link "Pular para o conteúdo".
- Responsivo: funciona a partir de 360 px de largura; as telas dos passos ocupam
  `min-height: 100svh` menos a navbar, com o conteúdo centralizado na vertical.

## 7. Regras fixas do projeto

1. Nenhum número solto no código Python: parâmetros em `config/*.toml`, lidos com `.get(...)` e
   fallback.
2. Código, comentários e textos em português; linhas de até 110 caracteres; nenhuma variável
   chamada `l`. `ruff check .` tem de passar.
3. Testes só com as sequências sintéticas de `tests/dados/`. Nenhuma sequência real.
4. Os 354 testes antigos continuam passando (com os ajustes da seção 3.5).
5. `primers` e `lampsim` nunca se importam (`tests/test_fronteiras.py`).
6. Templates com autoescape (padrão do Jinja2Templates); nada de `|safe` em texto de usuário.
7. Nenhuma dependência nova (nem Python, nem JS, nem CSS de CDN).
