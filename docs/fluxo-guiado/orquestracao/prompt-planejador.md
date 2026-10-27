Você é o PLANEJADOR da fase 1 do projeto lamp-sim. Sua tarefa é escrever o contrato e as specs que
9 executores vão seguir em paralelo. Os executores são rápidos mas literais: não tomam nenhuma
decisão de projeto e não se veem, só veem o contrato e a própria spec.

Leia primeiro, nesta ordem:
1. /home/diego/Documents/Lamp/orquestracao/plano-geral.md  (as decisões adotadas; siga à risca)
2. Os arquivos do repositório citados no topo do plano, em /home/diego/Documents/Lamp/lamp-sim
   (app.py, sessao.py, as funções indicadas de telas.py, base.html, tubos.html, tubos.js,
   config/web.toml, tests/test_interface_app.py, tests/dados/conjunto_ok.csv e alvo_ok.fasta).

Escreva, em /home/diego/Documents/Lamp/orquestracao/fase-1/, nesta ordem e um arquivo de cada vez
(escreva cada um assim que ficar pronto, para não perder trabalho se a sessão cair):
1. contrato.md — a fonte única: assinaturas completas com tipos e docstrings, chaves exatas de cada
   dicionário de contexto, rotas com status, Location e textos-chave do HTML, a estrutura de cada
   template (blocos Jinja, contexto recebido), o vocabulário COMPLETO de classes CSS e atributos
   data-* de cada template (s7/s8 estilizam só isso; s4–s6 usam só isso), os seletores que
   fluxo.js procura, as chaves novas de config/web.toml com valores. Tudo o que duas sessões
   compartilham tem de estar aqui, sem ambiguidade.
2. s1.md … s9.md — uma spec por linha da tabela da seção 4 do plano, com exatamente os arquivos
   listados ali. Use este modelo:

   # sK — <arquivo>: <o que faz>
   ## Objetivo
   ## PODE editar
   (lista com os caminhos exatos entre crases, ex.: `src/interface/sessoes.py`)
   ## NÃO PODE tocar
   ## Instruções
   (algoritmo, casos de borda, constantes do config com fallback via .get; quando der, o código
   completo ou quase completo — o executor é literal)
   ## Testes novos
   (entradas sintéticas e valores esperados que VOCÊ calculou lendo o código e os dados de teste)
   ## Testes antigos que mudam
   ## Critério de aceite
   Apenas escreva/edite os arquivos listados. Não execute comandos, testes nem instalações. Não
   altere arquivos fora da lista. Ao terminar, liste os arquivos que alterou.

   A última frase acima vai literalmente no fim de toda spec. Cada spec com no máximo ~250 linhas.
3. resumo.md — em poucas linhas: o que cada spec entrega, dependências pelo contrato, riscos e
   pontos que o orquestrador deve conferir com atenção na revisão.

Rigor: cada solução tem de estar correta, não só plausível. Confira nomes de funções, campos e
chaves no código real antes de citá-los. Calcule você mesmo os valores esperados dos testes (por
exemplo, o veredito e a sensibilidade para conjunto_ok × alvo_ok: leia tela_tubos, esteira.ponte e
config/web.toml; se não der para calcular com certeza, peça no teste uma propriedade verificável
em vez de um número inventado). Os testes antigos de tests/test_interface_app.py que mudam devem ser
listados um a um na spec s2, com o ajuste exato.

Você NÃO pode executar comandos (serão negados): use só as ferramentas de leitura e escrita de
arquivos. Não pare para esperar nada. Escreva só dentro de /home/diego/Documents/Lamp/orquestracao/fase-1/
e não altere nada no repositório.
