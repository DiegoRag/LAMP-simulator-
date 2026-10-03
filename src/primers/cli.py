"""Interface de linha de comando (CLI) para primers LAMP."""

import argparse
from pathlib import Path
import sys
from primers.analise import analisar
from primers.desenho import AVISO_ASCII, desenhar
from primers.entrada import carregar_config, ler_alvo, ler_conjunto
from primers.saida import exportar, exportar_desenho


def main(argv: list[str] | None = None) -> int:
    """Ponto de entrada da CLI."""
    parser = argparse.ArgumentParser(
        prog="primers",
        description="Ferramenta de primers LAMP: leitura, analise e verificacao.",
    )
    subparsers = parser.add_subparsers(dest="subcomando")

    p_analisar = subparsers.add_parser("analisar", help="Analisa um conjunto existente de primers.")
    p_analisar.add_argument("--entrada", required=True, help="Caminho para arquivo CSV ou FASTA com primers.")
    p_analisar.add_argument("--alvo", default=None, help="Caminho para arquivo FASTA com o alvo.")
    p_analisar.add_argument("--config", default="config/primers.toml", help="Caminho para primers.toml.")
    p_analisar.add_argument("--saida", required=True, help="Pasta de saida (dentro de saida/).")

    p_desenhar = subparsers.add_parser("desenhar", help="Desenha conjuntos de primers LAMP.")
    p_desenhar.add_argument("--alvo", required=True, help="Caminho para arquivo FASTA com o alvo.")
    p_desenhar.add_argument("--config", default="config/primers.toml", help="Caminho para primers.toml.")
    p_desenhar.add_argument("--n", type=int, default=None, help="Numero maximo de conjuntos a devolver.")
    p_desenhar.add_argument("--saida", required=True, help="Pasta de saida (dentro de saida/).")

    args = parser.parse_args(argv)

    if args.subcomando not in ("analisar", "desenhar"):
        parser.print_help(sys.stderr)
        return 2

    saida_p = Path(args.saida)
    saida_base = Path("saida").resolve()
    try:
        valido = saida_p.resolve().is_relative_to(saida_base)
    except (ValueError, AttributeError):
        valido = False

    if not valido:
        sys.stderr.write(f"erro: pasta de saida '{args.saida}' deve ficar dentro de 'saida/'.\n")
        return 2

    if args.subcomando == "analisar":
        try:
            cfg = carregar_config(args.config)
            alvo = ler_alvo(args.alvo) if args.alvo else None
            leitura = ler_conjunto(args.entrada, alvo=alvo)
        except ValueError as e:
            sys.stderr.write(f"erro: {e}\n")
            return 2

        resultado = analisar(leitura, alvo, cfg)
        exportar(resultado, args.saida, cfg)

        print(
            f"Analise concluida: {len(resultado.linhas)} primers, "
            f"{len(resultado.violacoes)} violacoes, "
            f"{len(resultado.avisos)} avisos. Pasta de saida: {args.saida}"
        )
        return 0

    try:
        cfg = carregar_config(args.config)
        alvo = ler_alvo(args.alvo)
    except ValueError as e:
        sys.stderr.write(f"erro: {e}\n")
        return 2

    resultado_desenho = desenhar(alvo, cfg, n=args.n)
    exportar_desenho(resultado_desenho, alvo, args.saida, cfg)

    print(AVISO_ASCII)
    stats = resultado_desenho.estatisticas
    cands_str = ", ".join(f"{k}: {v}" for k, v in stats["candidatos_por_tipo"].items())
    print(
        f"Conjuntos: {len(resultado_desenho.conjuntos)}, "
        f"candidatos por tipo: [{cands_str}], "
        f"examinados: {stats['conjuntos_examinados']}, "
        f"busca truncada: {stats['busca_truncada']}."
    )
    return 0
