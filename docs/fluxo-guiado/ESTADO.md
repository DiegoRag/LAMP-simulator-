# Fluxo guiado do LampSim: o que foi feito e onde paramos

Atualizado em 26/10/2026, por volta das 21h30, no notebook Linux do Diego. Branch
`feat/fluxo-guiado`, criado a partir de `feat/esteira`.

Para continuar em outra máquina ou sessão, leia este arquivo e depois o
[`plano-geral.md`](plano-geral.md), que tem todas as decisões. Os scripts de orquestração estão
em [`orquestracao/`](orquestracao/).

---

## 1. Linha do tempo de 26/10/2026

1. **Manual de identidade visual lido.** O arquivo `Manual de identidade visual LampSim.pdf`
   (v1.0) e o kit `lampsim-logo/` estão fora do repositório, em `~/Documents/Lamp/`. Resumo na
   seção 3.
2. **Auditoria de usabilidade do site atual**, com Playwright no Edge: percurso completo, erros
   de envio, telas de 1280 px e 375 px. Achados na seção 4. As capturas ficaram em
   `~/Documents/Lamp/.playwright-mcp/`, fora do Git.
3. **Mudança de prioridade do Diego.** Antes de corrigir a auditoria, o site ganha uma landing e
   um fluxo guiado de tela cheia: landing → passo 1 (primers) → passo 2 (gene-alvo) → resultado
   nos tubos. Perguntas e respostas na seção 2.
4. **Antigravity CLI (`agy`) instalado** e configurado para o fluxo de 10 agentes do
   `guia-delegacao-agentes.md`. Detalhes na seção 5.
5. **Commit local `fbc92ea` em `feat/esteira`**, sem push: o `.gitignore` passa a ignorar
   `graphify-out/` e `.delegar/`.
6. **Branch `feat/fluxo-guiado`**, com:
   - os logos oficiais do kit em `src/interface/static/marca/` (6 SVGs, sem alteração);
   - a Montserrat local em `src/interface/static/fontes/montserrat-latin.woff2` (variável,
     400–700, subconjunto latino, licença OFL; linha de copyright adicionada ao `OFL.txt`).
7. **Plano geral escrito** ([`plano-geral.md`](plano-geral.md)) e scripts de orquestração
   criados.
8. **O planejador (Opus 5.5 high no agy) falhou duas vezes** (seção 6). Paramos aqui.

---

## 2. Decisões do Diego sobre o fluxo guiado

| # | Pergunta | Decisão |
|---|---|---|
| 1 | Os tubos mostram um conjunto ou vários? | Sem resposta explícita; seguimos com **um conjunto por vez** (o que o formato de arquivo atual suporta) |
| 2 | Pular a Análise esconde o porquê do resultado | **Resumo simples no resultado**, com link para a análise detalhada |
| 3 | Estado global no site público | **Cookie de sessão** por navegador |
| 4 | Exemplo para quem não tem arquivo | **Depois.** O fluxo é para quem já sabe o que quer |
| 5 | Indicador de passo | **"Passo 1 de 2" e "Passo 2 de 2"** |
| 6 | Erros no fluxo | **Mensagem vermelha** indicando o problema, **sem mudar de página** |
| 7 | Navbar | **Início · Testar primers · Parâmetros**; o Desenho sai da navbar por enquanto |
| 8 | Exportar resultados | **Botão no resultado** |
| 9 | Histórico "Execuções desta sessão" | **Sai** por enquanto |
| 10 | O que é o LampSim (landing) | "Substitui as máquinas LAMP convencionais por uma versão muito menor, portátil e barata, que permite fazer exames genéticos em qualquer lugar, a qualquer hora." Sem as etapas físicas por enquanto; foco no online |
| 11 | Identidade visual | As telas novas **já nascem com a identidade LampSim** |
| 12 | Quem escreve o código | **10 agentes do agy**, do mesmo jeito que no PC de casa |

---

## 3. Identidade visual (resumo do manual)

- Cores: azul-marinho `#15344A` (base), rosa `#E83A73`, amarelo `#F4C44B`, branco. Dosagem de
  cerca de 55 % branco, 30 % azul, 8 % rosa e 7 % amarelo; rosa e amarelo aparecem pouco e, de
  preferência, juntos. Texto rosa pequeno sobre fundo claro usa `#C2275E`. Nada de degradê na
  marca.
- Fonte: Montserrat (Bold 700 títulos, SemiBold 600 rótulos e botões, Regular 400 texto). O nome
  da logo é um desenho: sempre usar os arquivos do kit.
- Logo: horizontal em cabeçalhos (mínimo 160 px de largura), completa em capas e na página
  principal (mínimo 120 px), símbolo no favicon (mínimo 32 px). Área de proteção igual à largura
  da cápsula. Não distorcer, recolorir, inverter, girar nem aplicar efeitos.
- As cores do vermelho de fenol nos tubos são científicas e não mudam para a paleta da marca.

---

## 4. Auditoria de usabilidade do site atual

Ainda vale para as páginas antigas. Parte dela é resolvida pelo fluxo guiado; o resto fica para
depois.

**Alta**
1. O erro de envio aparece no topo, mas a página volta rolada para o meio; o usuário não vê a
   mensagem, e o cartão continua mostrando o arquivo anterior como "Carregado". *(O fluxo novo
   resolve: o erro fica embaixo do botão.)*
2. Estado único para todos os visitantes: um erro enviado por fora apareceu no navegador. *(O
   cookie de sessão resolve.)*
3. Escolher o arquivo já envia sozinho; o botão "Enviar" é redundante. *(O fluxo novo resolve.)*

**Média**
4. Escolher uma proposta no ranking do Desenho recarrega a página e volta ao topo.
5. Falta um veredito em uma frase na Análise e nos Tubos. *(O resultado do fluxo novo resolve.)*
6. Na Análise, o botão de destaque é "Desenhar alternativas"; os Tubos ficam fora do trilho.
7. Números que não batem entre telas: "6 primers" no Início contra "8 primers lidos" na Análise;
   2 violações no resumo contra 8 na seção (já estava no `PENDENCIAS.md`); o Início diz "FASTA
   com um único registro", mas o código aceita vários (`primers/entrada.py:218`).
8. O cabeçalho mostra `conjunto_ok × alvo_ok` mesmo ao analisar uma proposta do desenho.
9. Alvo acima de 2.000 pb trava o Desenho sem oferecer saída.
10. O gráfico de Parâmetros tem rótulos sobrepostos (F3/B3, F2/B2, F1c/B1c).

**Baixa**
11. No celular, o trilho de etapas é cortado, o selo some e a chave de animações fica sem rótulo.
12. Os dois botões de escolher arquivo têm o mesmo nome acessível.
13. O link "Mostrar Tm fora da faixa como violação" aparece mesmo sem Tm fora da faixa.
14. A etapa marcada e o texto da tela de espera ficam dessincronizados.

**Visual (para a fase de identidade das páginas antigas):** o site se chama "LAMP Primers", usa
Open Sans e destaques laranja fora da paleta; as regiões F2/F1c (rosa) e B1c/B2c
(amarelo/laranja) podem ser confundidas com negativo/positivo.

**Bom:** textos de erro claros, foco visível pelo teclado, link "Pular para o conteúdo" e tubos
bem organizados no celular.

**Limite:** auditoria feita por IA. Falta testar com 3 a 5 pessoas da equipe biomédica.

---

## 5. Antigravity CLI e orquestração nesta máquina

- O que estava instalado era o **aplicativo de desktop** (`/usr/local/bin/antigravity`). O
  **CLI** foi instalado pelo script oficial (`https://antigravity.google/cli/install.sh`) em
  `~/.local/bin/agy`, versão 1.4.7, já autenticado.
- O settings fica em `~/.gemini/antigravity-cli/settings.json`, com
  `{"permissions": {"allow": [...]}}` e `write_file(...)` só para
  `~/Documents/Lamp/lamp-sim`, `lamp-sim-w1` … `lamp-sim-w10` e `~/Documents/Lamp/orquestracao`.
  Nada de `command(...)`.
- Os scripts do PC de casa não estavam aqui. As versões enxutas estão copiadas em
  [`orquestracao/`](orquestracao/) (os originais ficam em `~/Documents/Lamp/orquestracao/`):
  - `sincronizar.sh DEST…`: copia o repositório, sem `.git`, venv, dados e caches;
  - `executar-fase.sh N [K…]`: um executor por spec `sK.md`, em paralelo, cada um em
    `lamp-sim-wK`. `PRINCIPAL` e `RESERVA` podem ser trocados por variável de ambiente;
  - `integrar.py N BASE`: junta em `int-N` o que cada executor mudou e acusa conflitos e arquivos
    fora da lista;
  - `agy_retry.py`: troca para o modelo reserva quando a cota acaba ("RESOURCE_EXHAUSTED" ou
    "quota reached").
- As pastas `lamp-sim-w1` … `w10` e `orquestracao/base1` foram sincronizadas antes da linha nova
  do `OFL.txt`. Ao aplicar o `int-1`, copie só os arquivos alterados.

---

## 6. Onde paramos

- **Planejador, 1ª tentativa:** tentou rodar um comando, a ação foi negada e ele parou sem
  escrever nada (armadilha nº 4 do guia).
- **Planejador, 2ª tentativa:** retomado com `--conversation`, gastou cerca de 370 mil tokens e
  esgotou a **cota dos modelos Claude no agy** (Opus e Sonnet, cota da conta) sem gravar os
  arquivos. A cota volta por volta das **2h de 27/10/2026**.
- Os modelos **Gemini** (`gemini-4-argon-*`, `gemini-3.8-flash-*`) continuavam com cota.
- **Nenhum executor rodou** e nenhum código da fase foi escrito.

## 7. Próximos passos

1. Pelo guia, sem cota no planejador, o **Claude Code escreve `fase-1/contrato.md` e
   `s1.md … s9.md`** a partir do `plano-geral.md` e dos achados da seção 8.
2. Rodar os executores: `PRINCIPAL=gemini-3.8-flash-medium ./executar-fase.sh 1`, ou com Sonnet
   quando a cota voltar. O registrador (w10) roda depois.
3. `python3 integrar.py 1 base1` → `ruff` e `pytest` em `int-1` → revisar o diff → aplicar no
   repositório → registrador atualiza o `MELHORIAS.md`.
4. Subir o site, percorrer o fluxo novo no Edge (1280 e 375 px) e mostrar ao Diego.
5. Commit só com o "sim" do Diego.

## 8. Achados para o contrato (conferidos no código)

- **Passo 2, resumo dos primers:** usar `tela_inicio(sessao)["conjunto"]` →
  `"nome_arquivo"` e `"resumo"` (ex.: "CSV · 6 primers · FIP e BIP montados a partir das
  partes").
- **Resultado, números:** usar `tela_analise(execucao, cfg_primers, tm_como_calibracao,
  cfg_web["heterodimeros"]["escala_kcal"])["kpis"]`. Pegar os itens de rótulo "Violações"
  (`valor`, `detalhe`) e "Regiões localizadas" (o valor já vem como "6 de 6"). Eficiência:
  `texto_eficiencia`.
- **Veredito esperado para `conjunto_ok × alvo_ok`** (conferir no teste):
  - tom "alerta";
  - 2 violações, escrito como "2 regras do guia não cumpridas" (com `formato.plural`);
  - eficiência de 100 %;
  - sensibilidade "a partir de 1 cópia" (o tubo de 1 cópia vira em 15,0 min);
  - números de cópias formatados com `formato.milhar`.
- **Códigos:** `DIAGNOSTICO_TEXTO` (`positivo`, `negativo`, `controle_ok`, `invalido`) e
  `PAPEL_TEXTO` (`amostra`, `controle_negativo`, `controle_positivo`), em `telas.py`.
- **Middleware de sessão:** não criar sessão para caminhos `/static/`. O handler de 404 usa
  `getattr(request.state, "sessao", None)` ou um `Sessao` vazio novo (404 de estático).
- **Upload opcional nas rotas novas** (`UploadFile | None = None`): sem arquivo, mostrar o erro
  "Escolha um arquivo antes de enviar." re-renderizado com status 400. Isso evita um 422 em JSON.
- **Testes antigos de `tests/test_interface_app.py` que mudam:**
  - `test_inicio`: `/` vira a landing; o Início antigo vai para `/dados`;
  - `test_fluxo_analise`, `test_fluxo_desenho`, `test_mensagens_uma_vez`,
    `test_destino_em_json`, `test_destino_em_json_com_erro` e `test_envio_de_arquivo_em_json`:
    Location e destino `/` passam a `/dados`, e `GET /` passa a `GET /dados`;
  - `test_links_navegacao_tubos`: a navbar nova não tem "Tubos"; conferir o link "Testar
    primers" e o `aria-current` em `/parametros`.
- **Tubos no resultado:**
  - `tubos.js` não muda. Ele usa os atributos `data-*`, as classes `.tubo`, `.tubo-liquido`
    (variável `--cor`) e `.tubos-faixa-linha` (`.selecionada`), e a variável `--x` no marcador,
    no eixo e na marca de leitura;
  - s8 porta de `estilo.css` as regras das linhas ~2964–3393, ~3589–3620 e ~3766 em diante,
    sem o editor manual, trocando `var(--azul/--cinza/--linha/--laranja/--mono/--instrumento)`
    pelos tokens da marca;
  - pendência fora da fase 1: o gráfico de pH em `tubos.js` (linhas 313–378) desenha o canvas com
    Open Sans e as cores antigas.
