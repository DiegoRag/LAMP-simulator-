# Pendências do simulador LAMP

Questões que ficaram em aberto depois da correção da revisão de 17/10/2026 (6 fases,
commit `c9191a8`). O que já foi feito está em `MELHORIAS.md`.

## Dependem de dados de bancada

| Questão | Situação hoje | O que destrava |
|---|---|---|
| Faixas de Tm | Opção A adotada (convenção PrimerExplorer V5); Tm fora da faixa aparece como "a calibrar", não como erro (`config/web.toml`, `tm_como_calibracao`) | Comparar com um conjunto que funcionou em bancada |
| Cinética e eficiência do lampsim | Parâmetros provisórios, marcados "a calibrar" (`config/esteira.toml`, `src/lampsim/config.py`) | Curvas reais de virada por número de cópias |
| Concentrações dos primers | Mantidas em 200 / 1.600 / 400 nM (1:8:2, NEB); o artigo-base cita 1:4:1. Só a Tm da reação usa (`config/primers.toml`, `[concentracao_nM]`) | Conferir no artigo se os valores são da mistura ou finais |
| Limite de heterodímero | Só ordena propostas, não reprova (`peso_dimero = 0,5`, provisório) | Definir um limite de ΔG com a equipe |
| Valores "conferir" do guia da Eiken | GC, distâncias e condições do tampão ainda não checados no documento original | Ler o guia original |
| Faixas relaxadas do desenho | GC 30–70 % e Tm ±1 °C, provisórias, fora do guia (`[desenho_relaxado]`) | Validar com alvos ricos em AT ou GC |

## Código

| Questão | Detalhe | Esforço |
|---|---|---|
| Análise localiza o conjunto duas vezes | 6,6 s num genoma de 5 Mb; reaproveitar a busca deve baixar para ~3 s | Pequeno |
| Esteira sem linha de comando própria | `src/esteira/` só é usada pela interface e pelos testes | Pequeno |
| Título da seção Violações | Soma as violações de Tm, e o número fica diferente do KPI do topo | Pequeno |
| Coordenadas nos arquivos exportados | `conjuntos.csv` e `resultado.json` seguem em base 0 semiaberta (agora rotuladas); a tela usa base 1 | Decisão de contrato |

## Ficaram fora do desenho

- Ajustes do PrimerExplorer para alvos ricos em AT ou GC.
- Grampo e dímero do FIP e do BIP inteiros (hoje avaliados pela ponta 3′ ativa).
- Unicidade com erros de base (hoje a unicidade é exata, nas duas fitas).
- Alvo circular no desenho (a análise já aceita).

## Documentação

- A revisão no Claude Docs ainda descreve o site antigo. O manual foi atualizado em 18/10/2026.
- O merge de `feat/esteira` no `main` precisa de um pedido de integração com a esteira verde (ADR-0010).
