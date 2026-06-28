"""Submission-oriented checklist helpers for FASTA-changing outputs."""

from collections import Counter

from .agp import GAP_COMPONENT_TYPES, group_parts_by_object
from .paths import ensure_parent_dir


ALLOWED_FASTA_BASES = set("ACGTNacgtn")


def record_name(record):
    return getattr(record, "name")


def record_sequence(record):
    return getattr(record, "seq")


def add_row(rows, section, item, status, value, detail="."):
    rows.append(
        {
            "section": section,
            "item": item,
            "status": status,
            "value": value,
            "detail": detail,
        }
    )


def status_for_count(count):
    return "ok" if count == 0 else "warn"


def object_span_length(part):
    return part.object_end - part.object_start + 1


def component_span_length(part):
    if part.component_start is None or part.component_end is None:
        return 0
    return part.component_end - part.component_start + 1


def summarize_invalid_bases(records):
    summaries = []
    total = 0
    for record in records:
        seq = record_sequence(record)
        invalid = sorted({base for base in seq if base not in ALLOWED_FASTA_BASES})
        if invalid:
            count = sum(1 for base in seq if base not in ALLOWED_FASTA_BASES)
            total += count
            summaries.append(f"{record_name(record)}:{''.join(invalid)}({count})")
    return total, ";".join(summaries) if summaries else "."


def summarize_terminal_ns(records):
    names = []
    for record in records:
        seq = record_sequence(record).upper()
        if seq and (seq[0] == "N" or seq[-1] == "N"):
            names.append(record_name(record))
    return len(names), ",".join(names) if names else "."


def fasta_record_map(records):
    return {record_name(record): record_sequence(record) for record in records}


def check_agp_contiguity(parts):
    problems = []
    for object_name, object_parts in group_parts_by_object(parts).items():
        expected_start = 1
        expected_part_number = 1
        for part in object_parts:
            if part.part_number != expected_part_number:
                problems.append(
                    f"{object_name}:part_number {part.part_number} != {expected_part_number}"
                )
            if part.object_start != expected_start:
                problems.append(
                    f"{object_name}:object_start {part.object_start} != {expected_start}"
                )
            if part.object_end < part.object_start:
                problems.append(f"{object_name}:invalid span {part.object_start}-{part.object_end}")
            expected_start = part.object_end + 1
            expected_part_number += 1
    return problems


def check_agp_lengths(parts):
    problems = []
    for part in parts:
        span = object_span_length(part)
        if part.part_type == "gap":
            if part.gap_length != span:
                problems.append(
                    f"{part.object_name}:{part.part_number}:gap {part.gap_length} != span {span}"
                )
        else:
            component_span = component_span_length(part)
            if component_span != span:
                problems.append(
                    f"{part.object_name}:{part.part_number}:component {component_span} != span {span}"
                )
    return problems


def check_fasta_agp_objects(records, parts):
    if records is None:
        return []
    fasta_names = set(fasta_record_map(records))
    agp_names = set(group_parts_by_object(parts))
    problems = []
    missing_agp = sorted(fasta_names - agp_names)
    missing_fasta = sorted(agp_names - fasta_names)
    if missing_agp:
        problems.append(f"FASTA records without AGP rows: {','.join(missing_agp)}")
    if missing_fasta:
        problems.append(f"AGP objects without FASTA records: {','.join(missing_fasta)}")
    return problems


def check_fasta_agp_lengths(records, parts):
    if records is None:
        return []
    seqs = fasta_record_map(records)
    problems = []
    for object_name, object_parts in group_parts_by_object(parts).items():
        seq = seqs.get(object_name)
        if seq is None or not object_parts:
            continue
        agp_length = max(part.object_end for part in object_parts)
        if len(seq) != agp_length:
            problems.append(f"{object_name}:FASTA {len(seq)} bp != AGP {agp_length} bp")
    return problems


def check_gap_spans_are_ns(records, parts):
    if records is None:
        return []
    seqs = fasta_record_map(records)
    problems = []
    for part in parts:
        if part.part_type != "gap":
            continue
        seq = seqs.get(part.object_name)
        if seq is None:
            continue
        span = seq[part.object_start - 1 : part.object_end]
        if span.upper() != "N" * len(span):
            problems.append(f"{part.object_name}:{part.object_start}-{part.object_end}")
    return problems


def summarize_sources(parts):
    counts = Counter(part.source for part in parts)
    if not counts:
        return "."
    return ",".join(f"{source}:{count}" for source, count in sorted(counts.items()))


def summarize_statuses(parts):
    counts = Counter(part.status for part in parts)
    if not counts:
        return "."
    return ",".join(f"{status}:{count}" for status, count in sorted(counts.items()))


def problem_status(problems):
    return "ok" if not problems else "fail"


def problem_detail(problems):
    return ";".join(problems[:20]) if problems else "."


def write_submission_checklist(path, command_name, output_paths, fasta_records, agp_parts):
    rows = []
    records = None if fasta_records is None else list(fasta_records)
    parts = list(agp_parts)
    gap_parts = [part for part in parts if part.part_type == "gap"]
    component_parts = [part for part in parts if part.part_type != "gap"]
    graph_fill_parts = [
        part for part in component_parts
        if part.source == "graph_fill"
    ]
    fallback_gap_parts = [
        part for part in gap_parts
        if part.source in {"fallback_gap", "scaffold_gap", "gap"}
    ]

    add_row(rows, "metadata", "command", "info", command_name)
    for label, value in output_paths.items():
        add_row(rows, "files", label, "info", value, "keep_with_submission")

    if records is None:
        add_row(
            rows,
            "fasta",
            "records",
            "warn",
            ".",
            "No final FASTA output was written; run with --apply before submission.",
        )
        add_row(rows, "fasta", "total_bp", "warn", ".")
        add_row(rows, "fasta", "non_acgtn_bases", "warn", ".")
        add_row(rows, "fasta", "terminal_n_records", "warn", ".")
    else:
        invalid_count, invalid_detail = summarize_invalid_bases(records)
        terminal_n_count, terminal_n_detail = summarize_terminal_ns(records)
        add_row(rows, "fasta", "records", "ok", len(records))
        add_row(rows, "fasta", "total_bp", "ok", sum(len(record_sequence(record)) for record in records))
        add_row(
            rows,
            "fasta",
            "non_acgtn_bases",
            status_for_count(invalid_count),
            invalid_count,
            invalid_detail,
        )
        add_row(
            rows,
            "fasta",
            "terminal_n_records",
            status_for_count(terminal_n_count),
            terminal_n_count,
            terminal_n_detail,
        )

    add_row(rows, "agp", "objects", "ok", len(group_parts_by_object(parts)))
    add_row(rows, "agp", "parts", "ok", len(parts))
    add_row(rows, "agp", "component_parts", "ok", len(component_parts))
    add_row(rows, "agp", "gap_parts", "ok", len(gap_parts))
    add_row(rows, "agp", "gap_bp", "ok", sum(part.gap_length or 0 for part in gap_parts))
    add_row(
        rows,
        "agp",
        "gap_component_types",
        "ok",
        ",".join(sorted({part.component_type for part in gap_parts} & GAP_COMPONENT_TYPES)) or ".",
    )

    contiguity_problems = check_agp_contiguity(parts)
    length_problems = check_agp_lengths(parts)
    object_problems = check_fasta_agp_objects(records, parts)
    fasta_length_problems = check_fasta_agp_lengths(records, parts)
    gap_span_problems = check_gap_spans_are_ns(records, parts)
    fasta_check_status = "warn" if records is None else None
    fasta_check_detail = "not_checked_no_final_fasta" if records is None else None
    add_row(
        rows,
        "checks",
        "agp_contiguous_objects",
        problem_status(contiguity_problems),
        len(contiguity_problems),
        problem_detail(contiguity_problems),
    )
    add_row(
        rows,
        "checks",
        "agp_part_lengths",
        problem_status(length_problems),
        len(length_problems),
        problem_detail(length_problems),
    )
    add_row(
        rows,
        "checks",
        "fasta_agp_object_match",
        fasta_check_status or problem_status(object_problems),
        len(object_problems),
        fasta_check_detail or problem_detail(object_problems),
    )
    add_row(
        rows,
        "checks",
        "fasta_agp_length_match",
        fasta_check_status or problem_status(fasta_length_problems),
        len(fasta_length_problems),
        fasta_check_detail or problem_detail(fasta_length_problems),
    )
    add_row(
        rows,
        "checks",
        "agp_gap_spans_are_ns",
        fasta_check_status or problem_status(gap_span_problems),
        len(gap_span_problems),
        fasta_check_detail or problem_detail(gap_span_problems),
    )

    add_row(rows, "provenance", "component_sources", "info", summarize_sources(parts))
    add_row(rows, "provenance", "part_statuses", "info", summarize_statuses(parts))
    add_row(rows, "provenance", "unresolved_gap_parts", "info", len(fallback_gap_parts))
    add_row(rows, "provenance", "unresolved_gap_bp", "info", sum(part.gap_length or 0 for part in fallback_gap_parts))
    add_row(rows, "provenance", "graph_fill_parts", "info", len(graph_fill_parts))
    add_row(rows, "provenance", "graph_fill_bp", "info", sum(object_span_length(part) for part in graph_fill_parts))
    add_row(
        rows,
        "submission",
        "external_validation",
        "info",
        "required",
        "Run the NCBI submission validator or table2asn workflow on the final FASTA and AGP.",
    )

    header = ["section", "item", "status", "value", "detail"]
    ensure_parent_dir(path)
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for row in rows:
            out.write("\t".join(str(row[column]) for column in header) + "\n")
