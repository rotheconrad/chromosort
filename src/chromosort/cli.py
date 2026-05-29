"""Top-level ChromoSort command dispatcher."""

import argparse
import sys

from . import clean, cut, evaluate, fix_contigs, gapfill, manual, plot, reference_order, scaffold


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser(
        prog="chromo",
        description=(
            "Reference-order, clean, fix, cut, manually edit, scaffold, gapfill, "
            "plot, and review genome assembly contigs."
        ),
    )
    parser.add_argument(
        "command",
        nargs="?",
        choices=[
            "sort",
            "clean",
            "eval",
            "fix",
            "cut",
            "manual",
            "scaffold",
            "gapfill",
            "plot",
        ],
        help="Subcommand to run.",
    )

    if not argv or argv[0] in {"-h", "--help"}:
        parser.print_help()
        return

    command = argv[0]
    remaining = argv[1:]
    if command == "sort":
        reference_order.main(remaining, prog="chromo sort")
    elif command == "clean":
        clean.main(remaining, prog="chromo clean")
    elif command == "eval":
        evaluate.main(remaining, prog="chromo eval")
    elif command == "fix":
        fix_contigs.main(remaining, prog="chromo fix")
    elif command == "cut":
        cut.main(remaining, prog="chromo cut")
    elif command == "manual":
        manual.main(remaining, prog="chromo manual")
    elif command == "scaffold":
        scaffold.main(remaining, prog="chromo scaffold")
    elif command == "gapfill":
        gapfill.main(remaining, prog="chromo gapfill")
    elif command == "plot":
        plot.main(remaining, prog="chromo plot")
    else:
        parser.error(f"unknown command: {command}")


if __name__ == "__main__":
    main()
