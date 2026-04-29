"""Top-level ChromoSort command dispatcher."""

import argparse
import sys

from . import fix_contigs, reference_order


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser(
        prog="chromo",
        description="Reference-order and fix genome assembly contigs.",
    )
    parser.add_argument(
        "command",
        nargs="?",
        choices=["sort", "fix"],
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
    else:
        parser.error(f"unknown command: {command}")


if __name__ == "__main__":
    main()
