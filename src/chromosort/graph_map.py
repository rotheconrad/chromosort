#!/usr/bin/env python3
"""Project GFA unitig/path coordinates onto contig FASTA coordinates."""

import argparse
import sys
from pathlib import Path
from typing import Optional, Sequence

from .graph import (
    build_direct_projection,
    build_path_projection,
    summarize_projections,
    write_path_summary_tsv,
    write_projection_tsv,
    write_projection_warnings_tsv,
)
from .graph import GraphProjectionWarning, read_gfa
from .paths import ensure_parent_dir
from .reference_order import read_fasta_lengths


def parse_args(argv: Optional[Sequence[str]] = None, prog: Optional[str] = None):
    ap = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Build a GFA path/walk projection table from unitig graph coordinates "
            "to matching contig FASTA coordinates."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument(
        "--ctg-fasta",
        required=True,
        help="Contig FASTA whose record names should match GFA path/walk names.",
    )
    ap.add_argument(
        "--ctg-fai",
        default=None,
        help="Contig FASTA index. Defaults to <ctg-fasta>.fai when present.",
    )
    ap.add_argument(
        "--utg-gfa",
        required=True,
        help="Unitig or contig GFA, preferably a hifiasm .noseq.gfa when sequence is not needed.",
    )
    ap.add_argument(
        "--path-name",
        action="append",
        default=None,
        help=(
            "Specific GFA path/walk name to project. May be repeated. Defaults "
            "to all contig names from --ctg-fasta."
        ),
    )
    ap.add_argument(
        "--output-prefix",
        required=True,
        help=(
            "Output prefix. Writes <prefix>.utg_to_ctg.tsv, "
            "<prefix>.path_summary.tsv, and <prefix>.warnings.tsv."
        ),
    )
    ap.add_argument(
        "--trim-overlaps",
        action="store_true",
        help="Subtract each left overlap from the following path step span.",
    )
    ap.add_argument(
        "--length-tolerance-bp",
        type=int,
        default=1000,
        help="Absolute path-vs-FASTA length tolerance before warning.",
    )
    ap.add_argument(
        "--length-tolerance-frac",
        type=float,
        default=0.01,
        help="Fractional path-vs-FASTA length tolerance before warning.",
    )
    return ap.parse_args(argv)


def _output_paths(prefix):
    prefix = Path(prefix)
    return {
        "projection": Path(f"{prefix}.utg_to_ctg.tsv"),
        "summary": Path(f"{prefix}.path_summary.tsv"),
        "warnings": Path(f"{prefix}.warnings.tsv"),
    }


def _warn(warnings, severity, code, message):
    warnings.append(
        GraphProjectionWarning(
            severity=severity,
            code=code,
            path_name=".",
            segment=".",
            line_number=0,
            message=message,
        )
    )


def run(args):
    ctg_records, ctg_by_name = read_fasta_lengths(args.ctg_fasta, args.ctg_fai)
    requested_names = list(dict.fromkeys(args.path_name or [rec.name for rec in ctg_records]))
    fasta_lengths = {name: rec.length for name, rec in ctg_by_name.items()}
    warnings = []

    graph = read_gfa(args.utg_gfa)
    projections = build_path_projection(
        graph,
        path_names=requested_names,
        trim_overlaps=args.trim_overlaps,
        warnings=warnings,
    )

    if not projections:
        direct = build_direct_projection(graph, requested_names, warnings=warnings)
        if direct:
            _warn(
                warnings,
                "info",
                "direct_contig_projection",
                (
                    "No matching GFA P/W path projection was available; using direct "
                    "GFA segment-to-contig coordinates for matching FASTA names."
                ),
            )
            projections = direct

    summaries = summarize_projections(
        projections,
        fasta_lengths=fasta_lengths,
        requested_path_names=requested_names,
        tolerance_bp=args.length_tolerance_bp,
        tolerance_frac=args.length_tolerance_frac,
        warnings=warnings,
    )
    paths = _output_paths(args.output_prefix)
    for output_path in paths.values():
        ensure_parent_dir(output_path)
    write_projection_tsv(paths["projection"], projections)
    write_path_summary_tsv(paths["summary"], summaries)
    write_projection_warnings_tsv(paths["warnings"], warnings)
    return paths, projections, warnings


def main(argv: Optional[Sequence[str]] = None, prog: Optional[str] = None):
    try:
        args = parse_args(argv, prog=prog)
        paths, projections, warnings = run(args)
    except ValueError as exc:
        sys.stderr.write(f"ERROR: {exc}\n")
        raise SystemExit(2)

    sys.stderr.write(f"Wrote projection table: {paths['projection']}\n")
    sys.stderr.write(f"Wrote path summary: {paths['summary']}\n")
    sys.stderr.write(f"Wrote warnings: {paths['warnings']}\n")
    sys.stderr.write(f"Projected {len(projections)} GFA segment step(s).\n")
    if warnings:
        sys.stderr.write(f"Recorded {len(warnings)} warning/info row(s).\n")


if __name__ == "__main__":
    main()
