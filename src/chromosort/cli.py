"""Top-level ChromoSort command dispatcher."""

import argparse
import sys

from . import cut, fix_contigs, manual, plot, reference_order, scaffold


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser(
        prog="chromo",
        description=(
            "Reference-order, fix, cut, manually edit, scaffold, plot, and "
            "review genome assembly contigs."
        ),
    )
    parser.add_argument(
        "command",
        nargs="?",
        choices=["sort", "fix", "cut", "manual", "scaffold", "plot"],
        help="Subcommand to run.",
    )

    if not argv or argv[0] in {"-h", "--help"}:
        parser.print_help()
        return

    command = argv[0]
    remaining = argv[1:]
    if command == "sort":
        reference_order.main(remaining, prog="chromo sort")
    elif command == "fix":
        fix_contigs.main(remaining, prog="chromo fix")
    elif command == "cut":
        cut.main(remaining, prog="chromo cut")
    elif command == "manual":
        manual.main(remaining, prog="chromo manual")
    elif command == "scaffold":
        scaffold.main(remaining, prog="chromo scaffold")
    elif command == "plot":
        plot.main(remaining, prog="chromo plot")
    else:
        parser.error(f"unknown command: {command}")


if __name__ == "__main__":
    main()
