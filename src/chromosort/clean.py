#!/usr/bin/env python3
"""Conservatively clean mostly-correct assemblies with sort and fix logic."""

import argparse
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

from . import fix_contigs, reference_order
from .agp import AgpFastaRecord, component_part, write_agp, write_component_tsv
from .fix_contigs import MODE_CHOICES
from .paths import ensure_output_dirs, ensure_parent_dir
from .reference_order import (
    FastaReader,
    alignment_source_from_args,
    build_match_metrics,
    choose_assignments,
    fmt,
    order_assignments,
    read_fasta_lengths,
    resolve_duplicate_overlaps,
    reverse_complement,
    write_assignment_report,
    write_chromosome_summary,
    write_discarded_fasta,
    write_match_report,
    write_wrapped,
)
from .submission import write_submission_checklist


FIX_SCOPE_CHOICES = (
    "kept",
    "split-candidates",
    "kept-and-split-candidates",
    "file",
)


@dataclass
class CleanRecord:
    source_contig: str
    source_length: int
    clean_status: str
    clean_name: str
    ref: str
    ref_start: int
    ref_end: int
    orientation: str
    reverse_complemented: bool
    assignment: reference_order.Assignment
    fix_selected: bool = False
    fix_status: str = "."
    part_index: str = "."
    slice_start: int = 0
    slice_end: int = 0
    align_start: str = "."
    align_end: str = "."
    avg_identity: Optional[float] = None
    segment_count: str = "."
    reason: str = "."

    @property
    def piece_bp(self):
        return self.slice_end - self.slice_start


def parse_args(argv: Optional[Sequence[str]] = None, prog: Optional[str] = None):
    ap = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Conservatively clean a mostly-correct assembly by first applying "
            "chromo sort placement/filtering to raw contigs, then applying "
            "chromo fix planning to retained raw contigs, and finally writing "
            "one cleaned, reference-ordered FASTA."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--ref-fasta", required=True, help="Reference FASTA.")
    ap.add_argument(
        "--ref-fai",
        default=None,
        help="Reference FASTA index. Defaults to <ref-fasta>.fai when present.",
    )
    ap.add_argument(
        "-f",
        "--assembly-fasta",
        required=True,
        help="Raw assembly FASTA to clean.",
    )
    ap.add_argument(
        "--assembly-fai",
        default=None,
        help="Assembly FASTA index. Defaults to <assembly-fasta>.fai when present.",
    )
    alignment_group = ap.add_mutually_exclusive_group(required=True)
    alignment_group.add_argument(
        "-c",
        "--coords",
        help="MUMmer show-coords file generated from the raw assembly FASTA.",
    )
    alignment_group.add_argument(
        "--paf",
        help="minimap2 PAF file generated from the raw assembly FASTA.",
    )
    ap.add_argument(
        "-o",
        "--output-prefix",
        required=True,
        help=(
            "Output prefix. Writes <prefix>.clean.fa, initial-sort reports, "
            "<prefix>.clean.agp, <prefix>.clean_components.tsv, "
            "<prefix>.fix_report.tsv, <prefix>.clean_contigs.tsv, "
            "<prefix>.clean_chromosome_summary.tsv, "
            "<prefix>.submission_checklist.tsv, and <prefix>.run_summary.txt."
        ),
    )
    ap.add_argument(
        "--discarded-fasta",
        default=None,
        help="Optional FASTA path for raw contigs discarded by the sort step.",
    )
    ap.add_argument(
        "--reports-only",
        action="store_true",
        help="Write reports but skip clean/discarded FASTA output.",
    )
    ap.add_argument(
        "--name-separator",
        default="_",
        help="Separator used in cleaned FASTA IDs.",
    )
    ap.add_argument(
        "--simple-headers",
        action="store_true",
        help="Write FASTA headers containing only the cleaned sequence ID.",
    )
    ap.add_argument(
        "--orient-to-reference",
        action="store_true",
        help="Reverse-complement emitted records whose dominant alignment is reverse-strand.",
    )

    sort = ap.add_argument_group("sort/filter options")
    sort.add_argument("--min-aligned-bp", type=int, default=100_000)
    sort.add_argument("--min-query-cov", type=float, default=0.50)
    sort.add_argument("--min-best-ref-share", type=float, default=0.50)
    sort.add_argument("--large-alignment-min-bp", type=int, default=10_000_000)
    sort.add_argument("--large-alignment-min-query-cov", type=float, default=0.45)
    sort.add_argument("--min-novel-ref-bp", type=int, default=50_000)
    sort.add_argument("--min-novel-ref-frac", type=float, default=0.20)
    sort.add_argument("--overlap-mode", choices=["span", "alignment"], default="span")
    sort.add_argument("--novel-ref-criteria", choices=["both", "either"], default="both")
    sort.add_argument("--no-overlap-filter", action="store_true")
    sort.add_argument("--min-terminal-extension-bp", type=int, default=100_000)
    sort.add_argument("--min-terminal-extension-frac", type=float, default=0.02)
    sort.add_argument(
        "--no-terminal-overlap-rescue",
        dest="terminal_overlap_rescue",
        action="store_false",
    )
    sort.set_defaults(terminal_overlap_rescue=True)
    sort.add_argument(
        "--no-protect-split-candidates",
        dest="protect_split_candidates",
        action="store_false",
    )
    sort.set_defaults(protect_split_candidates=True)
    sort.add_argument("--split-candidate-min-aligned-bp", type=int, default=100_000)
    sort.add_argument("--split-candidate-min-query-frac", type=float, default=0.05)
    sort.add_argument("--split-candidate-max-best-share", type=float, default=0.95)

    shared = ap.add_argument_group("alignment row filters")
    shared.add_argument(
        "--min-segment-idy",
        type=float,
        default=0.0,
        help="Ignore individual alignment rows below this percent identity.",
    )
    shared.add_argument(
        "--min-mapq",
        type=int,
        default=0,
        help="Ignore PAF rows below this MAPQ. Ignored for MUMmer coords.",
    )
    shared.add_argument(
        "--include-secondary-paf",
        action="store_true",
        help="Include minimap2 PAF rows marked with tp:A:S.",
    )

    fix = ap.add_argument_group("fix options")
    fix.add_argument(
        "--fix-scope",
        choices=FIX_SCOPE_CHOICES,
        default="kept",
        help=(
            "Which original raw contigs to inspect with the fix planner. "
            "'kept' inspects all sort-retained raw contigs; 'split-candidates' "
            "inspects retained contigs flagged by sort; 'kept-and-split-candidates' "
            "also includes split candidates that were not retained; 'file' uses "
            "--fix-targets-file intersected with retained contigs."
        ),
    )
    fix.add_argument(
        "--fix-targets-file",
        default=None,
        help="Text file with one original raw contig ID per line; requires --fix-scope file.",
    )
    fix.add_argument("--fix-mode", choices=MODE_CHOICES, default="conservative")
    fix.add_argument("--min-segment-bp", type=int, default=10_000)
    fix.add_argument("--max-merge-gap", type=int, default=1_000)
    fix.add_argument("--min-piece-bp", type=int, default=1)
    fix.add_argument("--breakpoint-penalty-bp", type=float, default=50_000.0)
    fix.add_argument("--min-piece-aligned-bp", type=int, default=50_000)
    fix.add_argument("--min-piece-query-frac", type=float, default=0.05)
    fix.add_argument("--complex-inversion-min-piece-aligned-bp", type=int, default=1_000_000)
    fix.add_argument("--complex-inversion-min-overlap-frac", type=float, default=0.50)
    fix.add_argument("--max-breakpoints-per-contig", type=int, default=4)

    args = ap.parse_args(argv)
    if args.fix_scope == "file" and not args.fix_targets_file:
        ap.error("--fix-scope file requires --fix-targets-file")
    if args.fix_scope != "file" and args.fix_targets_file:
        ap.error("--fix-targets-file can only be used with --fix-scope file")
    return args


def output_paths(prefix, discarded_fasta=None):
    prefix = Path(prefix)
    return {
        "clean_fasta": Path(str(prefix) + ".clean.fa"),
        "clean_agp": Path(str(prefix) + ".clean.agp"),
        "clean_components": Path(str(prefix) + ".clean_components.tsv"),
        "initial_assignments": Path(str(prefix) + ".initial_sort.contig_assignments.tsv"),
        "initial_matches": Path(str(prefix) + ".initial_sort.contig_ref_matches.tsv"),
        "initial_chromosome_summary": Path(
            str(prefix) + ".initial_sort.chromosome_summary.tsv"
        ),
        "fix_targets": Path(str(prefix) + ".fix_targets.txt"),
        "fix_report": Path(str(prefix) + ".fix_report.tsv"),
        "clean_contigs": Path(str(prefix) + ".clean_contigs.tsv"),
        "clean_chromosome_summary": Path(str(prefix) + ".clean_chromosome_summary.tsv"),
        "submission_checklist": Path(str(prefix) + ".submission_checklist.tsv"),
        "run_summary": Path(str(prefix) + ".run_summary.txt"),
        **({"discarded_fasta": Path(discarded_fasta)} if discarded_fasta else {}),
    }


def read_name_file(path):
    names = []
    seen = set()
    with open(path) as fh:
        for line in fh:
            name = line.strip()
            if not name or name.startswith("#") or name in seen:
                continue
            names.append(name)
            seen.add(name)
    return names


def select_fix_targets(query_records, assignments, args):
    retained = {name for name, assignment in assignments.items() if assignment.kept}
    if args.fix_scope == "kept":
        selected = [rec.name for rec in query_records if rec.name in retained]
    elif args.fix_scope == "split-candidates":
        selected = [
            rec.name
            for rec in query_records
            if rec.name in retained and assignments[rec.name].split_candidate
        ]
    elif args.fix_scope == "kept-and-split-candidates":
        selected = [
            rec.name
            for rec in query_records
            if rec.name in retained or assignments[rec.name].split_candidate
        ]
    else:
        requested = read_name_file(args.fix_targets_file)
        selected = [name for name in requested if name in retained]
    return selected


def write_fix_targets(path, targets):
    ensure_parent_dir(path)
    with open(path, "w") as out:
        for target in targets:
            out.write(f"{target}\n")


def make_fix_args(args):
    return argparse.Namespace(
        mode=args.fix_mode,
        min_segment_bp=args.min_segment_bp,
        min_segment_idy=args.min_segment_idy,
        min_mapq=args.min_mapq,
        include_secondary_paf=args.include_secondary_paf,
        max_merge_gap=args.max_merge_gap,
        min_piece_bp=args.min_piece_bp,
        breakpoint_penalty_bp=args.breakpoint_penalty_bp,
        min_piece_aligned_bp=args.min_piece_aligned_bp,
        min_piece_query_frac=args.min_piece_query_frac,
        complex_inversion_min_piece_aligned_bp=args.complex_inversion_min_piece_aligned_bp,
        complex_inversion_min_overlap_frac=args.complex_inversion_min_overlap_frac,
        max_breakpoints_per_contig=args.max_breakpoints_per_contig,
        name_separator=args.name_separator,
        orient_to_reference=args.orient_to_reference,
        simple_headers=args.simple_headers,
        pieces_only=False,
    )


def plan_fixes(args, alignment_path, alignment_format, targets):
    fix_args = make_fix_args(args)
    if not targets:
        return {}, fix_args
    blocks = fix_contigs.collect_blocks(
        alignment_path,
        alignment_format,
        targets,
        args.min_segment_bp,
        args.min_segment_idy,
        args.max_merge_gap,
        min_mapq=args.min_mapq,
        include_secondary_paf=args.include_secondary_paf,
    )
    plans = fix_contigs.build_plans(args.assembly_fasta, targets, blocks, fix_args)
    fix_contigs.apply_breakpoint_guard(plans, fix_args)
    return plans, fix_args


def clean_base_name(ref, source_contig, part_index, separator):
    if part_index == ".":
        return f"{ref}{separator}{source_contig}"
    return f"{ref}{separator}{source_contig}{separator}{fix_contigs.alpha_label(part_index - 1)}"


def uniquify_records(records, separator):
    used = Counter()
    for record in records:
        base = record.clean_name
        used[base] += 1
        if used[base] > 1:
            record.clean_name = f"{base}{separator}{used[base]}"


def build_clean_records(query_records, assignments, plans, args, ref_order):
    ref_rank = {name: idx for idx, name in enumerate(ref_order)}
    records = []
    report_rows = []

    for rec in query_records:
        assignment = assignments[rec.name]
        if not assignment.kept or assignment.best is None:
            report_rows.append(
                CleanRecord(
                    source_contig=rec.name,
                    source_length=rec.length,
                    clean_status=f"discarded_{assignment.status}",
                    clean_name=".",
                    ref=assignment.best.ref if assignment.best else ".",
                    ref_start=assignment.best.ref_start if assignment.best else 0,
                    ref_end=assignment.best.ref_end if assignment.best else 0,
                    orientation=assignment.best.orientation if assignment.best else ".",
                    reverse_complemented=False,
                    assignment=assignment,
                    reason=f"Sort status: {assignment.status}.",
                )
            )
            continue

        plan = plans.get(rec.name)
        if plan is not None and plan.status == "split":
            for piece in plan.pieces:
                clean_name = clean_base_name(piece.ref, rec.name, piece.part_index, args.name_separator)
                record = CleanRecord(
                    source_contig=rec.name,
                    source_length=rec.length,
                    clean_status="kept_split_piece",
                    clean_name=clean_name,
                    ref=piece.ref,
                    ref_start=piece.ref_start,
                    ref_end=piece.ref_end,
                    orientation=piece.orientation,
                    reverse_complemented=piece.reverse_complemented,
                    assignment=assignment,
                    fix_selected=True,
                    fix_status=plan.status,
                    part_index=piece.part_index,
                    slice_start=piece.slice_start,
                    slice_end=piece.slice_end,
                    align_start=piece.align_start,
                    align_end=piece.align_end,
                    avg_identity=piece.avg_identity,
                    segment_count=piece.segment_count,
                    reason=plan.reason,
                )
                records.append(record)
                report_rows.append(record)
            continue

        fix_selected = plan is not None
        clean_status = plan.status if plan is not None else "kept_unsplit"
        reason = plan.reason if plan is not None else "Retained by sort; not selected for fixing."
        reverse = args.orient_to_reference and assignment.best.orientation == "-"
        record = CleanRecord(
            source_contig=rec.name,
            source_length=rec.length,
            clean_status=clean_status,
            clean_name=clean_base_name(assignment.best.ref, rec.name, ".", args.name_separator),
            ref=assignment.best.ref,
            ref_start=assignment.best.ref_start,
            ref_end=assignment.best.ref_end,
            orientation=assignment.best.orientation,
            reverse_complemented=reverse,
            assignment=assignment,
            fix_selected=fix_selected,
            fix_status=plan.status if plan is not None else ".",
            slice_start=0,
            slice_end=rec.length,
            avg_identity=assignment.best.avg_identity,
            segment_count=assignment.best.segment_count,
            reason=reason,
        )
        records.append(record)
        report_rows.append(record)

    records.sort(
        key=lambda record: (
            ref_rank.get(record.ref, float("inf")),
            record.ref_start,
            record.ref_end,
            record.source_contig,
            record.part_index if isinstance(record.part_index, int) else 0,
        )
    )
    uniquify_records(records, args.name_separator)
    report_name_by_key = {
        (record.source_contig, record.part_index, record.slice_start, record.slice_end): record.clean_name
        for record in records
    }
    for row in report_rows:
        key = (row.source_contig, row.part_index, row.slice_start, row.slice_end)
        if key in report_name_by_key:
            row.clean_name = report_name_by_key[key]
    return records, report_rows


def clean_fasta_header(record, simple_headers):
    if simple_headers:
        return record.clean_name
    fields = [
        record.clean_name,
        f"original={record.source_contig}",
        f"ref={record.ref}",
        f"ref_start={record.ref_start}",
        f"ref_end={record.ref_end}",
        f"slice={record.slice_start + 1}-{record.slice_end}",
        f"orientation={record.orientation}",
        f"reverse_complemented={'yes' if record.reverse_complemented else 'no'}",
        f"clean_status={record.clean_status}",
    ]
    return " ".join(fields)


def clean_record_sequence(reader, record):
    seq = reader.fetch(record.source_contig)
    seq = seq[record.slice_start : record.slice_end]
    if record.reverse_complemented:
        seq = reverse_complement(seq)
    return seq


def write_clean_fasta(path, fasta_path, records, assembly_fai, simple_headers):
    ensure_parent_dir(path)
    reader = FastaReader(fasta_path, assembly_fai)
    try:
        with open(path, "w") as out:
            for record in records:
                seq = clean_record_sequence(reader, record)
                out.write(f">{clean_fasta_header(record, simple_headers)}\n")
                write_wrapped(out, seq)
    finally:
        reader.close()


def build_clean_fasta_records(fasta_path, records, assembly_fai):
    fasta_records = []
    reader = FastaReader(fasta_path, assembly_fai)
    try:
        for record in records:
            fasta_records.append(
                AgpFastaRecord(record.clean_name, clean_record_sequence(reader, record))
            )
    finally:
        reader.close()
    return fasta_records


def build_clean_agp_parts(records):
    parts = []
    for record in records:
        notes = f"{record.ref}:{record.ref_start}-{record.ref_end}"
        parts.append(
            component_part(
                object_name=record.clean_name,
                object_start=1,
                object_end=record.piece_bp,
                part_number=1,
                component_id=record.source_contig,
                component_start=record.slice_start + 1,
                component_end=record.slice_end,
                orientation="-" if record.reverse_complemented else "+",
                source="clean_split_piece"
                if record.clean_status == "kept_split_piece"
                else "clean_contig",
                status=record.clean_status,
                notes=notes,
            )
        )
    return parts


def write_clean_report(path, rows):
    header = [
        "source_contig",
        "source_length",
        "clean_status",
        "clean_name",
        "kept_by_sort",
        "sort_status",
        "sort_assigned_ref",
        "sort_query_cov",
        "sort_best_ref_share",
        "sort_overlap_class",
        "fix_selected",
        "fix_status",
        "part_index",
        "slice_start",
        "slice_end",
        "piece_bp",
        "dominant_ref",
        "ref_start",
        "ref_end",
        "orientation",
        "reverse_complemented",
        "avg_identity",
        "segment_count",
        "reason",
    ]
    ensure_parent_dir(path)
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for row in rows:
            assignment = row.assignment
            best = assignment.best
            values = [
                row.source_contig,
                row.source_length,
                row.clean_status,
                row.clean_name,
                "yes" if assignment.kept else "no",
                assignment.status,
                best.ref if best else ".",
                fmt(best.query_cov) if best else ".",
                fmt(assignment.best_ref_share),
                assignment.overlap_class,
                "yes" if row.fix_selected else "no",
                row.fix_status,
                row.part_index,
                row.slice_start + 1 if row.clean_name != "." else ".",
                row.slice_end if row.clean_name != "." else ".",
                row.piece_bp if row.clean_name != "." else ".",
                row.ref,
                row.ref_start,
                row.ref_end,
                row.orientation,
                "yes" if row.reverse_complemented else "no",
                fmt(row.avg_identity) if row.avg_identity is not None else ".",
                row.segment_count,
                row.reason,
            ]
            out.write("\t".join(str(value) for value in values) + "\n")


def write_clean_chromosome_summary(path, ref_records, clean_records):
    by_ref = defaultdict(list)
    for record in clean_records:
        by_ref[record.ref].append(record)
    ensure_parent_dir(path)
    with open(path, "w") as out:
        out.write(
            "\t".join(
                [
                    "ref",
                    "ref_length",
                    "clean_records",
                    "covered_ref_bp",
                    "ref_cov",
                    "ordered_clean_names",
                ]
            )
            + "\n"
        )
        for ref in ref_records:
            records = by_ref.get(ref.name, [])
            intervals = reference_order.merge_intervals(
                [(record.ref_start, record.ref_end + 1) for record in records]
            )
            covered = reference_order.interval_bp(intervals)
            ref_cov = covered / ref.length if ref.length else 0.0
            out.write(
                "\t".join(
                    [
                        ref.name,
                        str(ref.length),
                        str(len(records)),
                        str(covered),
                        f"{ref_cov:.4f}",
                        ",".join(record.clean_name for record in records) if records else ".",
                    ]
                )
                + "\n"
            )


def write_run_summary(
    path,
    args,
    paths,
    ref_records,
    query_records,
    assignments,
    clean_records,
    fix_targets,
    plans,
    skipped_unknown_query,
    alignment_format,
):
    sort_counts = Counter(assignment.status for assignment in assignments.values())
    clean_counts = Counter(record.clean_status for record in clean_records)
    fix_counts = Counter(plan.status for plan in plans.values())
    ensure_parent_dir(path)
    with open(path, "w") as out:
        out.write("chromo clean\n")
        out.write("\nInputs\n")
        out.write(f"ref_fasta\t{args.ref_fasta}\n")
        out.write(f"assembly_fasta\t{args.assembly_fasta}\n")
        out.write(f"coords\t{args.coords or '.'}\n")
        out.write(f"paf\t{args.paf or '.'}\n")
        out.write(f"alignment_format\t{alignment_format}\n")
        out.write("\nOutputs\n")
        for label, value in paths.items():
            out.write(f"{label}\t{value}\n")
        out.write("\nSummary\n")
        out.write(f"reference_sequences\t{len(ref_records)}\n")
        out.write(f"raw_assembly_contigs\t{len(query_records)}\n")
        out.write(f"sort_retained_contigs\t{sum(1 for a in assignments.values() if a.kept)}\n")
        out.write(f"sort_discarded_contigs\t{sum(1 for a in assignments.values() if not a.kept)}\n")
        out.write(f"fix_scope\t{args.fix_scope}\n")
        out.write(f"fix_mode\t{args.fix_mode}\n")
        out.write(f"fix_targets\t{len(fix_targets)}\n")
        out.write(f"clean_records\t{len(clean_records)}\n")
        out.write(f"alignment_rows_skipped_unknown_query\t{skipped_unknown_query}\n")
        for status, count in sorted(sort_counts.items()):
            out.write(f"sort_status_{status}\t{count}\n")
        for status, count in sorted(fix_counts.items()):
            out.write(f"fix_status_{status}\t{count}\n")
        for status, count in sorted(clean_counts.items()):
            out.write(f"clean_status_{status}\t{count}\n")
        out.write(
            "\nNext steps\n"
            "The cleaned FASTA was generated from raw alignment evidence. "
            "For final validation, re-run MUMmer or minimap2 against the "
            "clean FASTA, then inspect chromo plot or mummerplot from that "
            "clean-FASTA alignment.\n"
        )


def main(argv: Optional[Sequence[str]] = None, prog: Optional[str] = None):
    args = parse_args(argv, prog=prog)
    alignment_path, alignment_format = alignment_source_from_args(args)
    paths = output_paths(args.output_prefix, args.discarded_fasta)
    ensure_output_dirs(paths)

    ref_records, ref_by_name = read_fasta_lengths(args.ref_fasta, args.ref_fai)
    query_records, query_by_name = read_fasta_lengths(args.assembly_fasta, args.assembly_fai)
    matches, by_query, skipped_unknown_query = build_match_metrics(
        alignment_path,
        alignment_format,
        ref_by_name,
        query_by_name,
        args.min_segment_idy,
        min_mapq=args.min_mapq,
        include_secondary_paf=args.include_secondary_paf,
    )

    assignments = choose_assignments(query_records, by_query, args)
    resolve_duplicate_overlaps(assignments, args)
    initial_kept = order_assignments(
        assignments,
        [rec.name for rec in ref_records],
        args.name_separator,
        args.orient_to_reference,
    )

    fix_targets = select_fix_targets(query_records, assignments, args)
    plans, fix_args = plan_fixes(args, alignment_path, alignment_format, fix_targets)
    clean_records, clean_report_rows = build_clean_records(
        query_records,
        assignments,
        plans,
        args,
        [rec.name for rec in ref_records],
    )

    write_assignment_report(paths["initial_assignments"], query_records, assignments)
    write_match_report(paths["initial_matches"], matches)
    write_chromosome_summary(paths["initial_chromosome_summary"], ref_records, initial_kept)
    write_fix_targets(paths["fix_targets"], fix_targets)
    fix_contigs.write_report(paths["fix_report"], fix_targets, plans)
    write_clean_report(paths["clean_contigs"], clean_report_rows)
    write_clean_chromosome_summary(paths["clean_chromosome_summary"], ref_records, clean_records)
    write_run_summary(
        paths["run_summary"],
        args,
        paths,
        ref_records,
        query_records,
        assignments,
        clean_records,
        fix_targets,
        plans,
        skipped_unknown_query,
        alignment_format,
    )

    agp_parts = build_clean_agp_parts(clean_records)
    write_agp(paths["clean_agp"], agp_parts)
    write_component_tsv(paths["clean_components"], agp_parts)

    fasta_records = None
    if not args.reports_only:
        write_clean_fasta(
            paths["clean_fasta"],
            args.assembly_fasta,
            clean_records,
            args.assembly_fai,
            args.simple_headers,
        )
        fasta_records = build_clean_fasta_records(
            args.assembly_fasta,
            clean_records,
            args.assembly_fai,
        )
        if args.discarded_fasta:
            write_discarded_fasta(
                paths["discarded_fasta"],
                args.assembly_fasta,
                assignments,
                args.assembly_fai,
            )
    write_submission_checklist(
        paths["submission_checklist"],
        "chromo clean",
        paths,
        fasta_records,
        agp_parts,
    )

    sys.stderr.write(f"Retained {len(clean_records)} cleaned record(s).\n")
    sys.stderr.write(f"Selected {len(fix_targets)} retained raw contig(s) for fixing.\n")
    for status, count in sorted(Counter(plan.status for plan in plans.values()).items()):
        sys.stderr.write(f"  fix {status}: {count}\n")
    if not args.reports_only:
        sys.stderr.write(f"Wrote clean FASTA: {paths['clean_fasta']}\n")
    sys.stderr.write(f"Wrote clean AGP: {paths['clean_agp']}\n")
    sys.stderr.write(f"Wrote clean components: {paths['clean_components']}\n")
    sys.stderr.write(f"Wrote submission checklist: {paths['submission_checklist']}\n")
    sys.stderr.write(f"Wrote clean report: {paths['clean_contigs']}\n")
    if fix_args.mode != args.fix_mode:
        sys.stderr.write(f"Fix mode normalized to: {fix_args.mode}\n")


if __name__ == "__main__":
    main()
