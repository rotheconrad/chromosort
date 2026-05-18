#!/usr/bin/env python3
"""Split FASTA contigs at user-specified cut positions."""

import argparse
import sys
from collections import defaultdict
from pathlib import Path
from typing import Optional, Sequence

from .reference_order import iter_fasta_records, read_fasta_lengths, write_wrapped


def parse_args(argv: Optional[Sequence[str]] = None, prog: Optional[str] = None):
    ap = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Cut assembly contigs at reviewed 1-based positions and write a new FASTA."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument(
        "-f",
        "--assembly-fasta",
        required=True,
        help="Assembly FASTA containing the contigs to cut.",
    )
    ap.add_argument(
        "--assembly-fai",
        default=None,
        help="Assembly FASTA index. Defaults to <assembly-fasta>.fai when present.",
    )
    ap.add_argument(
        "--cut",
        action="append",
        default=[],
        metavar="CONTIG:POS[,POS...]",
        help=(
            "Cut after one or more 1-based positions in a contig. May be repeated "
            "for multiple contigs."
        ),
    )
    ap.add_argument(
        "--cuts-file",
        default=None,
        help=(
            "Optional text/TSV/CSV file of reviewed cuts. Each non-comment line "
            "should contain contig and one or more positions."
        ),
    )
    ap.add_argument(
        "--contig",
        default=None,
        help=(
            "Convenience form for one contig, used with one or more --pos values. "
            "For multiple contigs, prefer repeated --cut entries or --cuts-file."
        ),
    )
    ap.add_argument(
        "--pos",
        action="append",
        default=[],
        nargs="+",
        type=int,
        help="Cut position(s) for --contig. May be repeated.",
    )
    ap.add_argument(
        "-o",
        "--output-fasta",
        required=True,
        help="Output FASTA with requested contigs replaced by cut pieces.",
    )
    ap.add_argument(
        "--report",
        required=True,
        help="TSV report describing every emitted cut piece.",
    )
    ap.add_argument(
        "--min-piece-bp",
        type=int,
        default=1,
        help="Minimum allowed length for every emitted piece.",
    )
    ap.add_argument(
        "--name-separator",
        default="_cut",
        help="Text inserted before the numeric cut-piece suffix.",
    )
    ap.add_argument(
        "--simple-headers",
        action="store_true",
        help="Write cut-piece FASTA headers containing only the new sequence ID.",
    )
    return ap.parse_args(argv)


def split_cut_file_line(line):
    clean = line.strip()
    if "\t" in clean:
        return [item.strip() for item in clean.split("\t") if item.strip()]
    if "," in clean:
        return [item.strip() for item in clean.split(",") if item.strip()]
    return clean.split()


def parse_position_token(token, source):
    try:
        return int(str(token).strip())
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid cut position {token!r} in {source}.") from exc


def parse_position_values(values, source):
    positions = []
    for value in values:
        for token in str(value).split(","):
            token = token.strip()
            if token:
                positions.append(parse_position_token(token, source))
    if not positions:
        raise ValueError(f"No cut positions found in {source}.")
    return positions


def parse_cut_spec(spec):
    if ":" not in spec:
        raise ValueError(
            f"Invalid --cut value {spec!r}; expected CONTIG:POS[,POS...]."
        )
    contig, raw_positions = spec.rsplit(":", 1)
    contig = contig.strip()
    if not contig:
        raise ValueError(f"Invalid --cut value {spec!r}; contig name is empty.")
    return contig, parse_position_values([raw_positions], f"--cut {spec!r}")


def read_cuts_file(path):
    rows = []
    with open(path) as fh:
        for line_number, line in enumerate(fh, start=1):
            clean = line.strip()
            if not clean or clean.startswith("#"):
                continue
            fields = split_cut_file_line(clean)
            if len(fields) < 2:
                raise ValueError(
                    f"Expected contig and position columns in {path}:{line_number}."
                )
            first = fields[0].lower()
            second = fields[1].lower()
            if first in {"contig", "sequence", "query"} and second.startswith("pos"):
                continue
            contig = fields[0]
            positions = parse_position_values(
                fields[1:],
                f"{path}:{line_number}",
            )
            rows.append((contig, positions))
    return rows


def collect_cut_positions(args):
    cuts = defaultdict(list)

    for spec in args.cut:
        contig, positions = parse_cut_spec(spec)
        cuts[contig].extend(positions)

    shorthand_positions = [pos for group in args.pos for pos in group]
    if args.contig or shorthand_positions:
        if not args.contig or not shorthand_positions:
            raise ValueError("Use --contig together with one or more --pos values.")
        cuts[args.contig].extend(shorthand_positions)

    if args.cuts_file:
        for contig, positions in read_cuts_file(args.cuts_file):
            cuts[contig].extend(positions)

    if not cuts:
        raise ValueError("Provide at least one cut via --cut, --cuts-file, or --contig/--pos.")

    return dict(cuts)


def duplicate_positions(positions):
    seen = set()
    duplicates = []
    duplicate_seen = set()
    for pos in positions:
        if pos in seen and pos not in duplicate_seen:
            duplicates.append(pos)
            duplicate_seen.add(pos)
        seen.add(pos)
    return duplicates


def normalize_cuts(raw_cuts, records, min_piece_bp):
    if min_piece_bp < 1:
        raise ValueError("--min-piece-bp must be at least 1.")

    lengths = {rec.name: rec.length for rec in records}
    normalized = {}
    for contig, positions in raw_cuts.items():
        if contig not in lengths:
            raise ValueError(f"Requested cut contig {contig!r} was not found in the assembly FASTA.")
        duplicates = duplicate_positions(positions)
        if duplicates:
            joined = ", ".join(str(pos) for pos in duplicates)
            raise ValueError(f"Duplicate cut position(s) for {contig!r}: {joined}.")

        length = lengths[contig]
        invalid = [pos for pos in positions if pos <= 0 or pos >= length]
        if invalid:
            joined = ", ".join(str(pos) for pos in invalid)
            raise ValueError(
                f"Cut position(s) for {contig!r} must be between 1 and {length - 1} "
                f"because cuts are made after the named base: {joined}."
            )

        sorted_positions = sorted(positions)
        starts = [0] + sorted_positions
        ends = sorted_positions + [length]
        short_pieces = [
            (start + 1, end, end - start)
            for start, end in zip(starts, ends)
            if end - start < min_piece_bp
        ]
        if short_pieces:
            preview = ", ".join(
                f"{start}-{end} ({size} bp)"
                for start, end, size in short_pieces[:5]
            )
            suffix = "..." if len(short_pieces) > 5 else ""
            raise ValueError(
                f"Cut plan for {contig!r} creates piece(s) shorter than "
                f"--min-piece-bp {min_piece_bp}: {preview}{suffix}."
            )
        normalized[contig] = sorted_positions

    return normalized


def piece_name(contig, part_index, separator):
    return f"{contig}{separator}{part_index:03d}"


def ensure_unique_output_names(records, cuts_by_contig, separator):
    seen = set()
    duplicates = []
    for rec in records:
        positions = cuts_by_contig.get(rec.name)
        if positions:
            names = [
                piece_name(rec.name, index, separator)
                for index in range(1, len(positions) + 2)
            ]
        else:
            names = [rec.name]
        for name in names:
            if name in seen and name not in duplicates:
                duplicates.append(name)
            seen.add(name)
    if duplicates:
        preview = ", ".join(duplicates[:5])
        suffix = "..." if len(duplicates) > 5 else ""
        raise ValueError(
            "Cut output would contain duplicate FASTA IDs: "
            f"{preview}{suffix}. Change --name-separator or input names."
        )


def cut_header(original_name, new_name, start, end, cut_positions, simple_headers):
    if simple_headers:
        return new_name
    return " ".join(
        [
            new_name,
            f"original={original_name}",
            f"slice={start + 1}-{end}",
            f"cut_after={','.join(str(pos) for pos in cut_positions)}",
        ]
    )


def write_cut_fasta(path, fasta_path, cuts_by_contig, separator, simple_headers):
    with open(path, "w") as out:
        for name, header, seq in iter_fasta_records(fasta_path):
            cut_positions = cuts_by_contig.get(name)
            if not cut_positions:
                out.write(header + "\n")
                write_wrapped(out, seq)
                continue

            starts = [0] + cut_positions
            ends = cut_positions + [len(seq)]
            for index, (start, end) in enumerate(zip(starts, ends), start=1):
                new_name = piece_name(name, index, separator)
                out.write(
                    f">{cut_header(name, new_name, start, end, cut_positions, simple_headers)}\n"
                )
                write_wrapped(out, seq[start:end])


def write_report(path, records, cuts_by_contig, separator):
    header = [
        "original_contig",
        "new_contig",
        "part_index",
        "slice_start",
        "slice_end",
        "piece_bp",
        "cut_after_positions",
    ]
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for rec in records:
            cut_positions = cuts_by_contig.get(rec.name)
            if not cut_positions:
                continue
            starts = [0] + cut_positions
            ends = cut_positions + [rec.length]
            for index, (start, end) in enumerate(zip(starts, ends), start=1):
                row = [
                    rec.name,
                    piece_name(rec.name, index, separator),
                    index,
                    start + 1,
                    end,
                    end - start,
                    ",".join(str(pos) for pos in cut_positions),
                ]
                out.write("\t".join(str(item) for item in row) + "\n")


def ensure_parent(path):
    parent = Path(path).parent
    if parent and str(parent) != ".":
        parent.mkdir(parents=True, exist_ok=True)


def run(args):
    raw_cuts = collect_cut_positions(args)
    records, _ = read_fasta_lengths(args.assembly_fasta, args.assembly_fai)
    cuts_by_contig = normalize_cuts(raw_cuts, records, args.min_piece_bp)
    ensure_unique_output_names(records, cuts_by_contig, args.name_separator)

    ensure_parent(args.output_fasta)
    ensure_parent(args.report)
    write_cut_fasta(
        args.output_fasta,
        args.assembly_fasta,
        cuts_by_contig,
        args.name_separator,
        args.simple_headers,
    )
    write_report(args.report, records, cuts_by_contig, args.name_separator)

    cut_contigs = len(cuts_by_contig)
    cut_positions = sum(len(positions) for positions in cuts_by_contig.values())
    emitted_pieces = cut_positions + cut_contigs
    sys.stderr.write(
        f"Cut {cut_contigs} contig(s) at {cut_positions} position(s), "
        f"emitting {emitted_pieces} piece(s).\n"
    )
    sys.stderr.write(f"Wrote cut FASTA: {args.output_fasta}\n")
    sys.stderr.write(f"Wrote cut report: {args.report}\n")


def main(argv: Optional[Sequence[str]] = None, prog: Optional[str] = None):
    args = parse_args(argv, prog=prog)
    try:
        run(args)
    except ValueError as exc:
        sys.stderr.write(f"ERROR: {exc}\n")
        sys.exit(2)


if __name__ == "__main__":
    main()
