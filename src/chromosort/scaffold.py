#!/usr/bin/env python3
"""
Build chromosome-level scaffold FASTA records from ChromoSort-ordered contigs.

The default gap model infers N runs from adjacent reference coordinates in the
chromo sort contig assignment report. Users can instead provide a fixed number
of Ns between neighboring contigs.
"""

import argparse
import csv
import sys
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

from .reference_order import iter_fasta_records, write_wrapped


REQUIRED_ASSIGNMENT_COLUMNS = {
    "contig",
    "kept",
    "new_name",
    "assigned_ref",
    "order_in_ref",
    "ref_start",
    "ref_end",
}


@dataclass
class AssignmentRow:
    contig: str
    new_name: str
    ref: str
    order_in_ref: int
    ref_start: int
    ref_end: int


@dataclass
class FastaRecord:
    name: str
    header: str
    seq: str


@dataclass
class ScaffoldMember:
    assignment: AssignmentRow
    record: FastaRecord


@dataclass
class GapRecord:
    scaffold: str
    left_contig: str
    right_contig: str
    left_ref_end: int
    right_ref_start: int
    raw_inferred_gap_bp: int
    gap_bp: int
    gap_mode: str


@dataclass
class ScaffoldRecord:
    name: str
    seq: str
    members: list
    gaps: list

    @property
    def sequence_bp(self):
        return sum(len(member.record.seq) for member in self.members)

    @property
    def gap_bp(self):
        return sum(gap.gap_bp for gap in self.gaps)


def parse_args(argv: Optional[Sequence[str]] = None, prog: Optional[str] = None):
    ap = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Scaffold a ChromoSort ordered FASTA into per-reference records "
            "using inferred or fixed N gaps."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument(
        "-f",
        "--ordered-fasta",
        required=True,
        help="Final ordered FASTA from chromo sort, optionally after chromo fix and re-sorting.",
    )
    ap.add_argument(
        "-a",
        "--assignments",
        required=True,
        help="Corresponding <prefix>.contig_assignments.tsv from chromo sort.",
    )
    ap.add_argument(
        "-o",
        "--output-prefix",
        required=True,
        help=(
            "Output prefix. Writes <prefix>.scaffold.fa, <prefix>.scaffold_gaps.tsv, "
            "<prefix>.scaffold_summary.tsv, and <prefix>.run_summary.txt."
        ),
    )
    ap.add_argument(
        "--fixed-gap-bp",
        "--fixed-gap",
        dest="fixed_gap_bp",
        type=int,
        default=None,
        help=(
            "Use this many Ns between neighboring contigs instead of inferring "
            "gap length from reference coordinates."
        ),
    )
    ap.add_argument(
        "--simple-headers",
        action="store_true",
        help="Write FASTA headers containing only the scaffold sequence ID.",
    )
    return ap.parse_args(argv)


def parse_int(value, field, row_name):
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Expected integer {field} for {row_name!r}, found {value!r}") from exc


def read_assignments(path):
    assignments = OrderedDict()
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = REQUIRED_ASSIGNMENT_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Assignment report {path} is missing required columns: {', '.join(sorted(missing))}"
            )
        for row in reader:
            kept = row["kept"].strip().lower() in {"yes", "true", "1"}
            if not kept:
                continue
            new_name = row["new_name"].strip()
            if not new_name or new_name == ".":
                continue
            if new_name in assignments:
                raise ValueError(f"Duplicate kept new_name {new_name!r} in {path}")
            assignments[new_name] = AssignmentRow(
                contig=row["contig"].strip(),
                new_name=new_name,
                ref=row["assigned_ref"].strip(),
                order_in_ref=parse_int(row["order_in_ref"], "order_in_ref", new_name),
                ref_start=parse_int(row["ref_start"], "ref_start", new_name),
                ref_end=parse_int(row["ref_end"], "ref_end", new_name),
            )
    return assignments


def read_ordered_fasta(path):
    records = []
    seen = set()
    for name, header, seq in iter_fasta_records(path):
        if name in seen:
            raise ValueError(f"Duplicate FASTA record {name!r} in {path}")
        seen.add(name)
        records.append(FastaRecord(name=name, header=header, seq=seq))
    return records


def group_scaffold_members(records, assignments):
    groups = OrderedDict()
    unassigned = []
    seen_assigned = set()

    for record in records:
        assignment = assignments.get(record.name)
        if assignment is None:
            unassigned.append(record)
            continue
        groups.setdefault(assignment.ref, []).append(
            ScaffoldMember(assignment=assignment, record=record)
        )
        seen_assigned.add(record.name)

    missing = [name for name in assignments if name not in seen_assigned]
    if missing:
        preview = ", ".join(missing[:5])
        suffix = "..." if len(missing) > 5 else ""
        raise ValueError(
            "Ordered FASTA is missing kept assignment sequence(s): "
            f"{preview}{suffix}. Use the ordered FASTA from the same chromo sort run."
        )

    return groups, unassigned


def inferred_gap(left, right):
    return right.assignment.ref_start - left.assignment.ref_end - 1


def build_scaffold(ref, members, fixed_gap_bp):
    pieces = []
    gaps = []
    gap_mode = "fixed" if fixed_gap_bp is not None else "inferred"

    for index, member in enumerate(members):
        if index:
            left = members[index - 1]
            raw_gap = inferred_gap(left, member)
            gap_bp = fixed_gap_bp if fixed_gap_bp is not None else max(0, raw_gap)
            gaps.append(
                GapRecord(
                    scaffold=ref,
                    left_contig=left.assignment.new_name,
                    right_contig=member.assignment.new_name,
                    left_ref_end=left.assignment.ref_end,
                    right_ref_start=member.assignment.ref_start,
                    raw_inferred_gap_bp=raw_gap,
                    gap_bp=gap_bp,
                    gap_mode=gap_mode,
                )
            )
            pieces.append("N" * gap_bp)
        pieces.append(member.record.seq)

    return ScaffoldRecord(
        name=ref,
        seq="".join(pieces),
        members=members,
        gaps=gaps,
    )


def build_scaffolds(groups, fixed_gap_bp):
    return [
        build_scaffold(ref, members, fixed_gap_bp)
        for ref, members in groups.items()
    ]


def scaffold_header(scaffold, simple_headers, gap_mode):
    if simple_headers:
        return scaffold.name
    fields = [
        scaffold.name,
        f"contigs={len(scaffold.members)}",
        f"sequence_bp={scaffold.sequence_bp}",
        f"gap_bp={scaffold.gap_bp}",
        f"gap_mode={gap_mode}",
    ]
    return " ".join(str(field) for field in fields)


def write_scaffold_fasta(path, scaffolds, unassigned, simple_headers, gap_mode):
    with open(path, "w") as out:
        for scaffold in scaffolds:
            out.write(f">{scaffold_header(scaffold, simple_headers, gap_mode)}\n")
            write_wrapped(out, scaffold.seq)
        for record in unassigned:
            out.write(record.header + "\n")
            write_wrapped(out, record.seq)


def write_gap_report(path, scaffolds):
    header = [
        "scaffold",
        "left_contig",
        "right_contig",
        "left_ref_end",
        "right_ref_start",
        "raw_inferred_gap_bp",
        "gap_bp",
        "gap_mode",
    ]
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for scaffold in scaffolds:
            for gap in scaffold.gaps:
                row = [
                    gap.scaffold,
                    gap.left_contig,
                    gap.right_contig,
                    gap.left_ref_end,
                    gap.right_ref_start,
                    gap.raw_inferred_gap_bp,
                    gap.gap_bp,
                    gap.gap_mode,
                ]
                out.write("\t".join(str(item) for item in row) + "\n")


def write_summary(path, scaffolds, unassigned):
    header = [
        "scaffold",
        "contigs",
        "scaffold_bp",
        "sequence_bp",
        "gap_bp",
        "gaps",
        "first_ref_start",
        "last_ref_end",
        "ordered_contigs",
    ]
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for scaffold in scaffolds:
            first_ref_start = scaffold.members[0].assignment.ref_start if scaffold.members else "."
            last_ref_end = scaffold.members[-1].assignment.ref_end if scaffold.members else "."
            row = [
                scaffold.name,
                len(scaffold.members),
                len(scaffold.seq),
                scaffold.sequence_bp,
                scaffold.gap_bp,
                len(scaffold.gaps),
                first_ref_start,
                last_ref_end,
                ",".join(member.assignment.new_name for member in scaffold.members),
            ]
            out.write("\t".join(str(item) for item in row) + "\n")
        for record in unassigned:
            row = [
                record.name,
                1,
                len(record.seq),
                len(record.seq),
                0,
                0,
                ".",
                ".",
                record.name,
            ]
            out.write("\t".join(str(item) for item in row) + "\n")


def write_run_summary(path, args, output_paths, scaffolds, unassigned):
    gap_mode = "fixed" if args.fixed_gap_bp is not None else "inferred"
    with open(path, "w") as out:
        out.write("chromo scaffold\n")
        out.write("\nInputs\n")
        out.write(f"ordered_fasta\t{args.ordered_fasta}\n")
        out.write(f"assignments\t{args.assignments}\n")
        out.write("\nGap model\n")
        out.write(f"gap_mode\t{gap_mode}\n")
        out.write(f"fixed_gap_bp\t{args.fixed_gap_bp if args.fixed_gap_bp is not None else '.'}\n")
        out.write("\nOutputs\n")
        for label, value in output_paths.items():
            out.write(f"{label}\t{value}\n")
        out.write("\nSummary\n")
        out.write(f"scaffolds\t{len(scaffolds)}\n")
        out.write(f"unassigned_records\t{len(unassigned)}\n")
        out.write(f"input_contigs_scaffolded\t{sum(len(scaffold.members) for scaffold in scaffolds)}\n")
        out.write(f"total_gap_bp\t{sum(scaffold.gap_bp for scaffold in scaffolds)}\n")


def run(args):
    if args.fixed_gap_bp is not None and args.fixed_gap_bp < 0:
        raise ValueError("--fixed-gap-bp must be zero or greater")

    prefix = Path(args.output_prefix)
    if prefix.parent and str(prefix.parent) != ".":
        prefix.parent.mkdir(parents=True, exist_ok=True)

    output_paths = {
        "scaffold_fasta": Path(str(prefix) + ".scaffold.fa"),
        "gap_report": Path(str(prefix) + ".scaffold_gaps.tsv"),
        "scaffold_summary": Path(str(prefix) + ".scaffold_summary.tsv"),
        "run_summary": Path(str(prefix) + ".run_summary.txt"),
    }

    assignments = read_assignments(args.assignments)
    records = read_ordered_fasta(args.ordered_fasta)
    groups, unassigned = group_scaffold_members(records, assignments)
    scaffolds = build_scaffolds(groups, args.fixed_gap_bp)
    gap_mode = "fixed" if args.fixed_gap_bp is not None else "inferred"

    write_scaffold_fasta(output_paths["scaffold_fasta"], scaffolds, unassigned, args.simple_headers, gap_mode)
    write_gap_report(output_paths["gap_report"], scaffolds)
    write_summary(output_paths["scaffold_summary"], scaffolds, unassigned)
    write_run_summary(output_paths["run_summary"], args, output_paths, scaffolds, unassigned)

    sys.stderr.write(f"Wrote scaffold FASTA: {output_paths['scaffold_fasta']}\n")
    sys.stderr.write(f"Wrote gap report: {output_paths['gap_report']}\n")
    sys.stderr.write(f"Wrote scaffold summary: {output_paths['scaffold_summary']}\n")


def main(argv: Optional[Sequence[str]] = None, prog: Optional[str] = None):
    args = parse_args(argv, prog=prog)
    try:
        run(args)
    except ValueError as exc:
        sys.stderr.write(f"ERROR: {exc}\n")
        sys.exit(2)


if __name__ == "__main__":
    main()
