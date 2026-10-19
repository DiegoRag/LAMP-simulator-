# lamp-sim

Frente de simulação do Projeto LAMP (Centro Universitário SENAC · FMABC), com simulação da reação LAMP colorimétrica, dados sintéticos, clusterização temporal e a ferramenta de primers.

Documentos de referência: LAMP-ORI-03 (orientação da área), LAMP-S02-ESP e LAMP-S02-PROT (simulador) e ADR-0013 (primers). A documentação do projeto vive no repositório `lamp-docs`.

**Todo resultado é provisório até a validação física.**

## O que tem aqui

```
lamp-sim/
├── src/
│   ├── primers/    análise e desenho de primers LAMP (Diego Ragozini)
│   ├── lampsim/    simulador: química, cor, classes e lotes por semente (Henrique Prado)
│   ├── esteira/    o que liga os dois: o contrato selecao.json e a ponte análise -> simulação (a dupla)
│   └── interface/  interface web local: análise e desenho no navegador
├── config/         primers.toml, esteira.toml e web.toml (os números do lampsim ficam em src/lampsim/config.py)
├── notebooks/      o simulador no Google Colab
├── tests/          testes automáticos (pytest)
└── .github/        CI: versões fixadas, nada de resultado versionado, estilo, testes e determinismo
```

Os dois pacotes nunca se importam (`tests/test_fronteiras.py`); a única ligação entre eles é o arquivo `selecao.json`, conferido por `esteira.contrato`.

## Como rodar no computador

Python 3.12 ou mais novo (numpy 2.5.2 e scipy 1.18.1, fixados em `requirements.txt`, não instalam em 3.11).

```bash
pip install -r requirements.txt
ruff check .
pytest -q
```

Sem primer3-py e biopython instalados (como no Colab), os testes do `primers` ficam fora da coleta e os do `lampsim` continuam rodando (`tests/conftest.py`).

## Simulador — `lampsim` (Henrique Prado)

**Situação: v0, fases 0 e 1.** O simulador gera a química (pH ao longo do tempo e forma do indicador) e a cor que sai dela. As imagens e as imperfeições de câmera entram na fase 2. Todo resultado é **provisório** até a validação física.

Módulos em `src/lampsim/`:
- `config.py` — todos os números do simulador (único lugar a editar)
- `quimica.py` — pH(t) logístico (normalizado para partir de zero) e Henderson-Hasselbalch
- `cinetica.py` — tempo da virada a partir da eficiência e do número de cópias (Poisson)
- `cor.py` — espectro de absorção do vermelho de fenol (Yao e Byrne, 2001) -> cor RGB (rota CIE e rota 3 bandas)
- `classes.py` — 5 classes ópticas, papel do tubo (amostra ou controle), sorteio, simulação e diagnóstico pela curva
- `manifesto.py` — lotes reprodutíveis por semente + manifesto com hash

### Como usar no Google Colab

1. Extraia o zip e envie a pasta `lamp-sim` inteira para o Google Drive, dentro de uma pasta `LAMP`. O caminho final deve ser `Meu Drive/LAMP/lamp-sim`.
2. No Drive, abra `lamp-sim/notebooks/00_setup.ipynb` com o Google Colaboratory (clique com o botão direito → *Abrir com*).
3. *Ambiente de execução → Executar tudo*. Autorize o acesso ao Drive quando pedir. No fim deve aparecer "Estilo e testes: tudo verde".
4. Abra `01_simulador.ipynb` do mesmo jeito e execute tudo.

Os dados gerados vão para `Meu Drive/LAMP/lamp-sim-dados` (lotes, resultados e logs), fora da pasta do código, para não entrarem no git. Não é preciso GPU.

Para editar o código, altere os arquivos `.py` de `src/lampsim/` direto no Drive. No Colab, reinicie a sessão (*Ambiente de execução → Reiniciar sessão*) para o notebook carregar a versão nova.

Fora do Colab, os notebooks gravam em `_dados_local/` (ignorado pelo git).

### Premissas e o que está em aberto

- **Indicador:** vermelho de fenol (pKa ≈ 8,0). É premissa: a composição do mix ainda não foi confirmada pela FMABC.
- **Faixas de t₀, k, A e pH:** provisórias. Serão calibradas com os tempos do artigo, por conjunto de primers, na fase 5.
- **Espectros:** bandas gaussianas ajustadas para dar a cor visual correta (rosa e amarelo). A banda secundária da forma básica é um ajuste empírico. Devem ser substituídas por um espectro medido quando houver.
- **Rotas de cor:** as duas estão implementadas (`rota="cie"` e `rota="3bandas"`) para a decisão de faseamento da LAMP-S02-ESP, Seção 2.4.

## Primers — `primers` (Diego Ragozini)

Lê e analisa um conjunto existente de primers LAMP e desenha conjuntos novos a partir da sequência-alvo.

Comandos rodados da raiz do repositório com `PYTHONPATH=src`, ou depois de `pip install -e .`:

```bash
python -m primers analisar --entrada conjunto.csv [--alvo alvo.fasta] --saida saida/analise
python -m primers desenhar --alvo alvo.fasta [--n 5] [--relaxado] --saida saida/desenho
```

Os parâmetros ficam em `config/primers.toml`. O desenho percorre o alvo inteiro (teto de candidatos por janela de
`[desenho] janela_nt` nt) e cada proposta vem de um núcleo F2–B2c diferente. Se nenhuma proposta sair, a CLI diz
qual regra eliminou os candidatos; `--relaxado` repete o desenho com as faixas de `[desenho_relaxado]` (GC e Tm
alargadas, fora do guia da Eiken). A saída fica sempre dentro de `saida/`, que não vai para o git.

Toda saída de `desenhar` é candidata computacional, não validada experimentalmente: a definição do conjunto que vai ao tubo é da equipe biomédica (ADR-0013). Dados reais ficam em `dados/`, fora do git, e sequência de isolado clínico nunca entra no repositório.

## A esteira e o contrato

A esteira é a única parte que importa os dois pacotes.

O contrato `lamp-sim/selecao@1`, descrito em `src/esteira/contrato.py`, é o arquivo que a etapa de seleção de primers escreve e a etapa de simulação lê; o exemplo versionado está em `tests/dados/selecao_v1_exemplo.json`. Mudar um campo é subir a versão do esquema, nunca alterar a v1 em silêncio.

A versão 2 (`lamp-sim/selecao@2`) aceita IDs livres, bases IUPAC e a origem `usuario` (conjunto carregado, sem afirmação de bancada). Por alelo, `incompat_3p` conta os erros na ponta que estende (3′ dos primers; 5′ de F1c e B1c, onde o laço fecha) e o campo opcional `incompat_regioes` dá a conta de cada região. Quem escreve é `esteira.ponte.selecao_da_analise`: a exportação da Análise, na interface, põe o `selecao.json` no zip, com `f_previsto`. Quem lê é `esteira.ponte.tubos_da_selecao`, que simula os tubos com o `f_previsto` do arquivo.

As demais etapas da esteira (escolher os dados de entrada, traduzir casamento em cinética, salvar a execução, a linha de comando) ainda vão entrar.

## Interface web

Interface web local com as telas Início, Análise, Desenho e Parâmetros da ferramenta de primers, no navegador.

Para abrir a partir da raiz do repositório:

```bash
PYTHONPATH=src python -m interface
```

(ou `python -m interface` depois de `pip install -e .`), acessando em seguida o endereço `http://127.0.0.1:8000`. O host e a porta ficam em `config/web.toml`.

A interface só escuta nesta máquina; arquivos e execuções ficam em memória e somem quando o servidor para; a exportação baixa um zip com os mesmos arquivos da linha de comando (a da Análise traz também o `selecao.json`); as fontes vêm do Google Fonts; sem internet, a página usa a fonte do sistema. Nenhuma sequência é enviada.

## Versionamento

- O repositório já existe no GitHub (`LAMP-Senac-FMABC/lamp-sim`); não se cria outro.
- Não há ramo `main`; o ramo padrão é `docs/kit-protecao`.
- Trabalho novo entra num ramo `feat/`, `fix/` ou `exp/` e sobe por pedido de integração, com a CI verde (ADR-0010); nunca `git push --force`.
- O gancho local se ativa uma vez por clone, com `git config core.hooksPath .githooks`.
- Versione a receita, não o resultado (LAMP-GIT-002): os manifestos de lote (`manifesto_<hash>.json`) podem ser versionados; trajetórias e imagens geradas, não.
