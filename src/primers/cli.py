"""Interface de linha de comando (CLI) para primers LAMP."""

import argparse
import sys
import unicodedata
from pathlib import Path

from primers.analise import analisar
from primers.desenho import AVISO_ASCII, desenhar
from primers.entrada import carregar_config, escolher_registro_alvo, ler_alvo, ler_alvos, ler_conjunto
from primers.saida import descrever_alvo, exportar, exportar_desenho


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
    p_desenhar.add_argument(
        "--relaxado",
        action="store_true",
        help="Usa as faixas de [desenho_relaxado] do config (GC e Tm alargadas, fora do guia).",
    )

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
            registro_alvo = escolher_registro_alvo(ler_alvos(args.alvo)) if args.alvo else None
            alvo = registro_alvo.sequencia if registro_alvo is not None else None
            ligante_cfg = cfg.get("desenho", {}).get("ligante", "")
            max_erros_cfg = int(cfg.get("casamento", {}).get("max_erros", 3))
            circular_cfg = bool(cfg.get("alvo", {}).get("circular", False))
            leitura = ler_conjunto(
                args.entrada,
                alvo=alvo,
                ligante=ligante_cfg,
                max_erros=max_erros_cfg,
                circular=circular_cfg,
            )
        except ValueError as e:
            sys.stderr.write(f"erro: {e}\n")
            return 2

        resultado = analisar(leitura, alvo, cfg)
        alvo_json = None
        if registro_alvo is not None:
            alvo_json = descrever_alvo(
                registro_alvo.sequencia,
                id_registro=registro_alvo.id,
                arquivo=Path(args.alvo).name,
                circular=circular_cfg,
            )
        exportar(resultado, args.saida, cfg, leitura=leitura, alvo=alvo_json)

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

    resultado_desenho = desenhar(alvo, cfg, n=args.n, relaxado=args.relaxado)
    exportar_desenho(resultado_desenho, alvo, args.saida, cfg)

    print(AVISO_ASCII)
    stats = resultado_desenho.estatisticas
    cands_str = ", ".join(f"{k}: {v}" for k, v in stats["candidatos_por_tipo"].items())
    print(
        f"Conjuntos: {len(resultado_desenho.conjuntos)}, "
        f"candidatos por tipo: [{cands_str}], "
        f"nucleos examinados: {stats['conjuntos_examinados']}, "
        f"nucleos distintos: {stats['nucleos_distintos']}, "
        f"faixas relaxadas: {stats['relaxado']}."
    )
    if stats["diagnostico"]:
        print(f"Diagnostico: {_ascii(stats['diagnostico'])}")
    return 0


def _ascii(texto: str) -> str:
    """Texto sem acentos, para consoles que não mostram UTF-8 (como o AVISO_ASCII)."""
    return unicodedata.normalize("NFKD", texto).encode("ascii", "ignore").decode("ascii")
