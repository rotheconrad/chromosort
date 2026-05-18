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
    trim_left_bp: int = 0
    trim_right_bp: int = 0

    @property
    def seq(self):
        end = len(self.record.seq) - self.trim_right_bp
        return self.record.seq[self.trim_left_bp : end]

    @property
    def trimmed_bp(self):
        return self.trim_left_bp + self.trim_right_bp


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
    overlap_bp: int
    overlap_class: str
    overlap_left_ref_frac: float
    overlap_right_ref_frac: float
    overlap_policy: str
    overlap_action: str
    trimmed_bp: int
    sequence_overlap_identity: Optional[float] = None


@dataclass
class ScaffoldRecord:
    name: str
    seq: str
    members: list
    gaps: list

    @property
    def sequence_bp(self):
        return sum(len(member.seq) for member in self.members)

    @property
    def gap_bp(self):
        return sum(gap.gap_bp for gap in self.gaps)

    @property
    def overlap_bp(self):
        return sum(gap.overlap_bp for gap in self.gaps)

    @property
    def trimmed_bp(self):
        return sum(member.trimmed_bp for member in self.members)


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
        "--overlap-policy",
        choices=["zero-gap", "warn", "trim-reference", "trim-sequence"],
        default="zero-gap",
        help=(
            "How to handle adjacent contigs with overlapping reference spans. "
            "zero-gap writes no Ns and keeps both sequences unchanged; warn does "
            "the same but emits stderr warnings; trim-reference removes the "
            "reference-inferred overlap from the right contig; trim-sequence "
            "removes the right-contig prefix only when the terminal sequence "
            "overlap has high identity."
        ),
    )
    ap.add_argument(
        "--trim-sequence-min-identity",
        type=float,
        default=0.98,
        help=(
            "Minimum suffix/prefix identity required by "
            "--overlap-policy trim-sequence."
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


def ref_span_bp(assignment):
    return max(0, assignment.ref_end - assignment.ref_start + 1)


def fmt(value, digits=4):
    if value is None:
        return "."
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def classify_adjacent_overlap(left, right, overlap_bp):
    if overlap_bp <= 0:
        return "no_overlap"
    left_start = left.assignment.ref_start
    left_end = left.assignment.ref_end
    right_start = right.assignment.ref_start
    right_end = right.assignment.ref_end
    if right_start > left_start and right_end > left_end:
        return "terminal_overlap"
    if right_start >= left_start and right_end <= left_end:
        return "contained_overlap"
    if right_start <= left_start and right_end >= left_end:
        return "spanning_overlap"
    return "internal_overlap"


def sequence_overlap_identity(left_seq, right_seq, overlap_bp):
    size = min(overlap_bp, len(left_seq), len(right_seq))
    if size <= 0:
        return 0, None
    left_part = left_seq[-size:].upper()
    right_part = right_seq[:size].upper()
    matches = sum(1 for left_base, right_base in zip(left_part, right_part) if left_base == right_base)
    return size, matches / size


def apply_overlap_policy(left, right, overlap_bp, overlap_class, args):
    if overlap_bp <= 0:
        return "none", 0, None
    if args.overlap_policy in {"zero-gap", "warn"}:
        return "zero_gap", 0, None
    if overlap_class != "terminal_overlap":
        return "trim_skipped_nonterminal", 0, None

    if args.overlap_policy == "trim-reference":
        trim_bp = min(overlap_bp, len(right.seq))
        right.trim_left_bp += trim_bp
        return "trimmed_reference", trim_bp, None

    trim_bp, identity = sequence_overlap_identity(left.seq, right.seq, overlap_bp)
    if trim_bp > 0 and identity is not None and identity >= args.trim_sequence_min_identity:
        right.trim_left_bp += trim_bp
        return "trimmed_sequence", trim_bp, identity
    return "trim_skipped_sequence_identity", 0, identity


def build_scaffold(ref, members, fixed_gap_bp, args):
    pieces = []
    gaps = []
    gap_mode = "fixed" if fixed_gap_bp is not None else "inferred"

    for index, member in enumerate(members):
        if index:
            left = members[index - 1]
            raw_gap = inferred_gap(left, member)
            overlap_bp = max(0, -raw_gap)
            overlap_class = classify_adjacent_overlap(left, member, overlap_bp)
            overlap_action, trimmed_bp, sequence_identity = apply_overlap_policy(
                left,
                member,
                overlap_bp,
                overlap_class,
                args,
            )
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
                    overlap_bp=overlap_bp,
                    overlap_class=overlap_class,
                    overlap_left_ref_frac=(
                        overlap_bp / ref_span_bp(left.assignment)
                        if ref_span_bp(left.assignment)
                        else 0.0
                    ),
                    overlap_right_ref_frac=(
                        overlap_bp / ref_span_bp(member.assignment)
                        if ref_span_bp(member.assignment)
                        else 0.0
                    ),
                    overlap_policy=args.overlap_policy,
                    overlap_action=overlap_action,
                    trimmed_bp=trimmed_bp,
                    sequence_overlap_identity=sequence_identity,
                )
            )
            if overlap_bp > 0 and (
                args.overlap_policy == "warn" or overlap_action.startswith("trim")
            ):
                sys.stderr.write(
                    "WARNING: "
                    f"{ref} overlap between {left.assignment.new_name} and "
                    f"{member.assignment.new_name}: raw_gap={raw_gap}, "
                    f"overlap_bp={overlap_bp}, class={overlap_class}, "
                    f"action={overlap_action}, trimmed_bp={trimmed_bp}\n"
                )
            pieces.append("N" * gap_bp)
        pieces.append(member.seq)

    return ScaffoldRecord(
        name=ref,
        seq="".join(pieces),
        members=members,
        gaps=gaps,
    )


def build_scaffolds(groups, fixed_gap_bp, args):
    return [
        build_scaffold(ref, members, fixed_gap_bp, args)
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
        "overlap_bp",
        "overlap_class",
        "overlap_left_ref_frac",
        "overlap_right_ref_frac",
        "overlap_policy",
        "overlap_action",
        "trimmed_bp",
        "sequence_overlap_identity",
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
                    gap.overlap_bp,
                    gap.overlap_class,
                    fmt(gap.overlap_left_ref_frac),
                    fmt(gap.overlap_right_ref_frac),
                    gap.overlap_policy,
                    gap.overlap_action,
                    gap.trimmed_bp,
                    fmt(gap.sequence_overlap_identity, 3),
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
        "overlap_gaps",
        "overlap_bp",
        "trimmed_bp",
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
                sum(1 for gap in scaffold.gaps if gap.overlap_bp > 0),
                scaffold.overlap_bp,
                scaffold.trimmed_bp,
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
                0,
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
        out.write(f"overlap_policy\t{args.overlap_policy}\n")
        out.write(f"trim_sequence_min_identity\t{args.trim_sequence_min_identity}\n")
        out.write("\nOutputs\n")
        for label, value in output_paths.items():
            out.write(f"{label}\t{value}\n")
        out.write("\nSummary\n")
        out.write(f"scaffolds\t{len(scaffolds)}\n")
        out.write(f"unassigned_records\t{len(unassigned)}\n")
        out.write(f"input_contigs_scaffolded\t{sum(len(scaffold.members) for scaffold in scaffolds)}\n")
        out.write(f"total_gap_bp\t{sum(scaffold.gap_bp for scaffold in scaffolds)}\n")
        out.write(f"total_overlap_bp\t{sum(scaffold.overlap_bp for scaffold in scaffolds)}\n")
        out.write(f"total_trimmed_bp\t{sum(scaffold.trimmed_bp for scaffold in scaffolds)}\n")


def run(args):
    if args.fixed_gap_bp is not None and args.fixed_gap_bp < 0:
        raise ValueError("--fixed-gap-bp must be zero or greater")
    if not 0.0 <= args.trim_sequence_min_identity <= 1.0:
        raise ValueError("--trim-sequence-min-identity must be between 0 and 1")

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
    scaffolds = build_scaffolds(groups, args.fixed_gap_bp, args)
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
