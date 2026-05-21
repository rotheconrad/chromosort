#!/usr/bin/env python3
"""Plan and optionally apply reviewed graph-based gap fills."""

import argparse
import sys
from collections import defaultdict, deque
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

from .graph import read_gfa
from .reference_order import reverse_complement, write_wrapped
from .scaffold import (
    classify_adjacent_overlap,
    graph_intermediate_nodes,
    graph_node_name,
    graph_orientation,
    graph_path_nodes,
    graph_start_keys,
    graph_target_keys,
    group_scaffold_members,
    inferred_gap,
    read_assignments,
    read_ordered_fasta,
)


@dataclass
class FillPlan:
    scaffold: str
    left_contig: str
    right_contig: str
    left_graph_node: str
    right_graph_node: str
    left_orientation: str
    right_orientation: str
    raw_inferred_gap_bp: int
    fallback_gap_bp: int
    overlap_bp: int
    overlap_class: str
    graph_status: str
    path_edges: Optional[int]
    path_nodes: str
    intermediate_nodes: str
    candidate_paths: int
    fill_status: str
    gaf_path_support: Optional[int] = None
    gaf_best_alt_support: Optional[int] = None
    hic_path_support: Optional[int] = None
    hic_best_alt_support: Optional[int] = None
    fill_sequence: str = ""
    fill_bp: int = 0
    right_trim_bp: int = 0
    applied: bool = False
    reason: str = "."


@dataclass
class FilledScaffold:
    name: str
    seq: str
    contigs: int
    filled_gaps: int
    fallback_gaps: int
    fill_bp: int
    fallback_gap_bp: int
    trimmed_bp: int


@dataclass(frozen=True)
class GafPathRecord:
    query: str
    path: list
    mapq: int


def parse_args(argv: Optional[Sequence[str]] = None, prog: Optional[str] = None):
    ap = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Plan graph-supported fills between adjacent ChromoSort sorted "
            "contigs, and optionally apply only unambiguous sequence paths."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument(
        "-f",
        "--ordered-fasta",
        required=True,
        help="Final ordered FASTA from chromo sort, optionally after chromo fix.",
    )
    ap.add_argument(
        "-a",
        "--assignments",
        required=True,
        help="Corresponding <prefix>.contig_assignments.tsv from chromo sort.",
    )
    ap.add_argument(
        "--gfa",
        required=True,
        help="Assembly graph GFA containing segment sequences and links.",
    )
    ap.add_argument(
        "--gaf",
        default=None,
        help=(
            "Optional GAF graph alignments. When provided, read paths can "
            "resolve otherwise-ambiguous GFA branches if one candidate path "
            "has unique support."
        ),
    )
    ap.add_argument(
        "--hic-pairs",
        default=None,
        help=(
            "Optional TSV of graph-node Hi-C contact counts. Expected columns "
            "are node_a, node_b, and count; a header row is allowed."
        ),
    )
    ap.add_argument(
        "-o",
        "--output-prefix",
        required=True,
        help=(
            "Output prefix. Writes <prefix>.fill_plan.tsv and "
            "<prefix>.run_summary.txt; with --apply also writes "
            "<prefix>.filled.fa."
        ),
    )
    ap.add_argument(
        "--apply",
        action="store_true",
        help="Write <prefix>.filled.fa using only fillable graph paths.",
    )
    ap.add_argument(
        "--fixed-gap-bp",
        type=int,
        default=None,
        help=(
            "Use this many Ns for unresolved gaps in --apply output instead "
            "of inferring gap length from reference coordinates."
        ),
    )
    ap.add_argument(
        "--max-path-edges",
        type=int,
        default=4,
        help="Maximum GFA link depth searched between adjacent sorted contigs.",
    )
    ap.add_argument(
        "--max-candidate-paths",
        type=int,
        default=2,
        help=(
            "Stop path enumeration after this many candidate paths. The "
            "default is enough to distinguish unique from ambiguous paths."
        ),
    )
    ap.add_argument(
        "--min-gaf-mapq",
        type=int,
        default=20,
        help="Minimum GAF MAPQ for a read path to support a candidate fill.",
    )
    ap.add_argument(
        "--min-gaf-path-support",
        type=int,
        default=1,
        help=(
            "Minimum supporting GAF read paths required to resolve an "
            "ambiguous GFA branch."
        ),
    )
    ap.add_argument(
        "--min-hic-path-support",
        type=int,
        default=1,
        help=(
            "Minimum summed Hi-C contact support required to resolve an "
            "ambiguous GFA branch."
        ),
    )
    ap.add_argument(
        "--max-fill-bp",
        type=int,
        default=1_000_000,
        help="Maximum inserted graph sequence allowed for one fill. Set negative to disable.",
    )
    ap.add_argument(
        "--include-fill-sequences",
        action="store_true",
        help="Include candidate fill sequences in the TSV plan.",
    )
    ap.add_argument(
        "--simple-headers",
        action="store_true",
        help="Write filled FASTA headers containing only the scaffold ID.",
    )
    return ap.parse_args(argv)


def graph_status_for_nodes(left_node, right_node):
    if left_node == "." and right_node == ".":
        return "missing_both_nodes"
    if left_node == ".":
        return "missing_left_node"
    if right_node == ".":
        return "missing_right_node"
    return None


def parse_gaf_path(path_text):
    nodes = []
    index = 0
    while index < len(path_text):
        orientation = path_text[index]
        if orientation not in {">", "<"}:
            raise ValueError(f"Malformed GAF path component near {path_text[index:]!r}")
        next_index = index + 1
        while next_index < len(path_text) and path_text[next_index] not in {">", "<"}:
            next_index += 1
        name = path_text[index + 1 : next_index]
        if not name:
            raise ValueError(f"Malformed empty GAF path component in {path_text!r}")
        nodes.append((name, "+" if orientation == ">" else "-"))
        index = next_index
    return nodes


def read_gaf(path, min_mapq=0):
    records = []
    with open(path) as fh:
        for line_number, line in enumerate(fh, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 12:
                raise ValueError(f"Malformed GAF row at line {line_number}: expected at least 12 columns")
            try:
                mapq = int(cols[11])
            except ValueError as exc:
                raise ValueError(f"Malformed GAF MAPQ at line {line_number}: {cols[11]!r}") from exc
            if mapq < min_mapq:
                continue
            records.append(
                GafPathRecord(
                    query=cols[0],
                    path=parse_gaf_path(cols[5]),
                    mapq=mapq,
                )
            )
    return records


def reverse_oriented_nodes(nodes):
    flipped = {"+": "-", "-": "+"}
    return [(name, flipped[orientation]) for name, orientation in reversed(nodes)]


def contains_oriented_subpath(read_path, candidate):
    if not candidate or len(candidate) > len(read_path):
        return False
    for index in range(0, len(read_path) - len(candidate) + 1):
        if read_path[index : index + len(candidate)] == candidate:
            return True
    return False


def path_oriented_nodes(path):
    if not path:
        return []
    nodes = [path[0].source_key]
    nodes.extend(edge.target_key for edge in path)
    return nodes


def gaf_support_for_path(gaf_records, path):
    candidate = path_oriented_nodes(path)
    reverse_candidate = reverse_oriented_nodes(candidate)
    support = 0
    for record in gaf_records:
        if contains_oriented_subpath(record.path, candidate) or contains_oriented_subpath(
            record.path,
            reverse_candidate,
        ):
            support += 1
    return support


def gaf_path_supports(gaf_records, paths):
    return [gaf_support_for_path(gaf_records, path) for path in paths]


def choose_gaf_supported_path(paths, supports, min_support):
    if not paths or not supports:
        return None
    best_support = max(supports)
    if best_support < min_support:
        return None
    if supports.count(best_support) != 1:
        return None
    return supports.index(best_support)


def read_hic_pairs(path):
    contacts = defaultdict(int)
    seen_data_row = False
    with open(path) as fh:
        for line_number, line in enumerate(fh, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                raise ValueError(
                    f"Malformed Hi-C pair row at line {line_number}: "
                    "expected 3 columns"
                )
            is_first_data_row = not seen_data_row
            seen_data_row = True
            try:
                count = int(cols[2])
            except ValueError as exc:
                if is_first_data_row:
                    continue
                raise ValueError(
                    f"Malformed Hi-C count at line {line_number}: {cols[2]!r}"
                ) from exc
            if count < 0:
                raise ValueError(
                    f"Malformed Hi-C count at line {line_number}: "
                    "count must be non-negative"
                )
            node_a, node_b = cols[0], cols[1]
            if not node_a or not node_b:
                raise ValueError(
                    f"Malformed Hi-C pair row at line {line_number}: "
                    "empty node name"
                )
            contacts[tuple(sorted((node_a, node_b)))] += count
    return dict(contacts)


def hic_contact_count(hic_contacts, node_a, node_b):
    if not hic_contacts:
        return 0
    return hic_contacts.get(tuple(sorted((node_a, node_b))), 0)


def hic_support_for_path(hic_contacts, path):
    nodes = [name for name, _ in path_oriented_nodes(path)]
    return sum(
        hic_contact_count(hic_contacts, left, right)
        for left, right in zip(nodes, nodes[1:])
    )


def hic_path_supports(hic_contacts, paths):
    return [hic_support_for_path(hic_contacts, path) for path in paths]


def choose_unique_supported_path(supports, min_support):
    if not supports:
        return None
    best_support = max(supports)
    if best_support < min_support:
        return None
    if supports.count(best_support) != 1:
        return None
    return supports.index(best_support)


def enumerate_graph_paths(
    graph,
    left_node,
    right_node,
    left_orientation,
    right_orientation,
    max_edges,
    max_paths,
):
    if max_edges < 1 or max_paths < 1:
        return []

    targets = graph_target_keys(right_node, right_orientation)
    if not targets:
        return []

    queue = deque()
    for start in graph_start_keys(left_node, left_orientation):
        queue.append((start, [], {start}))

    paths = []
    while queue and len(paths) < max_paths:
        (node, orientation), path, seen = queue.popleft()
        if len(path) >= max_edges:
            continue
        for edge in graph.outgoing(node, orientation):
            next_key = edge.target_key
            if next_key in seen:
                continue
            next_path = path + [edge]
            if next_key in targets:
                paths.append(next_path)
                if len(paths) >= max_paths:
                    break
                continue
            queue.append((next_key, next_path, seen | {next_key}))
    return paths


def oriented_sequence(seq, orientation):
    if seq is None:
        return None
    if orientation == "+":
        return seq
    if orientation == "-":
        return reverse_complement(seq)
    return None


def sequence_matches(observed, expected):
    return observed.upper() == expected.upper()


def reconstruct_fill_sequence(graph, path, left, right, max_fill_bp):
    if not path:
        return None, 0, "no_graph_path"

    left_seq = oriented_sequence(graph.nodes[path[0].source].sequence, path[0].source_orientation)
    right_seq = oriented_sequence(graph.nodes[path[-1].target].sequence, path[-1].target_orientation)
    if left_seq is None or right_seq is None:
        return None, 0, "unsequenced_flank_node"
    if not sequence_matches(left.record.seq, left_seq):
        return None, 0, "left_flank_sequence_mismatch"
    if not sequence_matches(right.record.seq, right_seq):
        return None, 0, "right_flank_sequence_mismatch"

    parts = []
    for index, edge in enumerate(path[:-1]):
        overlap_bp = edge.overlap_bp
        if overlap_bp is None:
            return None, 0, "unknown_overlap"
        node = graph.nodes[edge.target]
        seq = oriented_sequence(node.sequence, edge.target_orientation)
        if seq is None:
            return None, 0, "unsequenced_intermediate_node"
        if overlap_bp > len(seq):
            return None, 0, "invalid_overlap"
        parts.append(seq[overlap_bp:])

        next_edge = path[index + 1]
        if next_edge.overlap_bp is None:
            return None, 0, "unknown_overlap"

    right_trim_bp = path[-1].overlap_bp
    if right_trim_bp is None:
        return None, 0, "unknown_overlap"
    if right_trim_bp > len(right.record.seq):
        return None, 0, "invalid_overlap"

    fill_sequence = "".join(parts)
    if max_fill_bp >= 0 and len(fill_sequence) > max_fill_bp:
        return None, 0, "fill_too_long"
    return fill_sequence, right_trim_bp, "."


def plan_gap(
    scaffold_name,
    left,
    right,
    graph,
    args,
    gaf_records=None,
    hic_contacts=None,
):
    gaf_records = gaf_records or []
    hic_contacts = hic_contacts or {}
    raw_gap = inferred_gap(left, right)
    overlap_bp = max(0, -raw_gap)
    fallback_gap_bp = args.fixed_gap_bp if args.fixed_gap_bp is not None else max(0, raw_gap)
    overlap_class = classify_adjacent_overlap(left, right, overlap_bp)
    left_node = graph_node_name(left, graph)
    right_node = graph_node_name(right, graph)
    left_orientation = graph_orientation(left)
    right_orientation = graph_orientation(right)

    missing_status = graph_status_for_nodes(left_node, right_node)
    if missing_status is not None:
        return FillPlan(
            scaffold=scaffold_name,
            left_contig=left.assignment.new_name,
            right_contig=right.assignment.new_name,
            left_graph_node=left_node,
            right_graph_node=right_node,
            left_orientation=left_orientation,
            right_orientation=right_orientation,
            raw_inferred_gap_bp=raw_gap,
            fallback_gap_bp=fallback_gap_bp,
            overlap_bp=overlap_bp,
            overlap_class=overlap_class,
            graph_status=missing_status,
            path_edges=None,
            path_nodes=".",
            intermediate_nodes=".",
            candidate_paths=0,
            fill_status="missing_node",
            reason=missing_status,
        )

    paths = enumerate_graph_paths(
        graph,
        left_node,
        right_node,
        left_orientation,
        right_orientation,
        args.max_path_edges,
        args.max_candidate_paths,
    )
    if not paths:
        return FillPlan(
            scaffold=scaffold_name,
            left_contig=left.assignment.new_name,
            right_contig=right.assignment.new_name,
            left_graph_node=left_node,
            right_graph_node=right_node,
            left_orientation=left_orientation,
            right_orientation=right_orientation,
            raw_inferred_gap_bp=raw_gap,
            fallback_gap_bp=fallback_gap_bp,
            overlap_bp=overlap_bp,
            overlap_class=overlap_class,
            graph_status="no_path",
            path_edges=None,
            path_nodes=".",
            intermediate_nodes=".",
            candidate_paths=0,
            fill_status="no_graph_path",
            reason="no_graph_path",
        )

    supports = gaf_path_supports(gaf_records, paths) if gaf_records else []
    hic_supports = hic_path_supports(hic_contacts, paths) if hic_contacts else []
    selected_index = 0
    graph_status = "direct_edge" if len(paths[0]) == 1 else "short_path"

    if len(paths) > 1:
        gaf_choice = choose_gaf_supported_path(
            paths,
            supports,
            args.min_gaf_path_support,
        )
        hic_choice = choose_unique_supported_path(
            hic_supports,
            args.min_hic_path_support,
        )
        if gaf_choice is not None and hic_choice is not None and gaf_choice != hic_choice:
            path = paths[0]
            selected_support = supports[0] if supports else None
            alt_support = max(supports[1:], default=0) if supports else None
            selected_hic_support = hic_supports[0] if hic_supports else None
            alt_hic_support = max(hic_supports[1:], default=0) if hic_supports else None
            return FillPlan(
                scaffold=scaffold_name,
                left_contig=left.assignment.new_name,
                right_contig=right.assignment.new_name,
                left_graph_node=left_node,
                right_graph_node=right_node,
                left_orientation=left_orientation,
                right_orientation=right_orientation,
                raw_inferred_gap_bp=raw_gap,
                fallback_gap_bp=fallback_gap_bp,
                overlap_bp=overlap_bp,
                overlap_class=overlap_class,
                graph_status="ambiguous_paths",
                path_edges=len(path),
                path_nodes=graph_path_nodes(path),
                intermediate_nodes=graph_intermediate_nodes(path),
                candidate_paths=len(paths),
                fill_status="ambiguous_paths",
                gaf_path_support=selected_support,
                gaf_best_alt_support=alt_support,
                hic_path_support=selected_hic_support,
                hic_best_alt_support=alt_hic_support,
                reason="conflicting_gaf_hic_support",
            )
        if gaf_choice is None and hic_choice is None:
            path = paths[0]
            selected_support = supports[0] if supports else None
            alt_support = max(supports[1:], default=0) if supports else None
            selected_hic_support = hic_supports[0] if hic_supports else None
            alt_hic_support = max(hic_supports[1:], default=0) if hic_supports else None
            return FillPlan(
                scaffold=scaffold_name,
                left_contig=left.assignment.new_name,
                right_contig=right.assignment.new_name,
                left_graph_node=left_node,
                right_graph_node=right_node,
                left_orientation=left_orientation,
                right_orientation=right_orientation,
                raw_inferred_gap_bp=raw_gap,
                fallback_gap_bp=fallback_gap_bp,
                overlap_bp=overlap_bp,
                overlap_class=overlap_class,
                graph_status="ambiguous_paths",
                path_edges=len(path),
                path_nodes=graph_path_nodes(path),
                intermediate_nodes=graph_intermediate_nodes(path),
                candidate_paths=len(paths),
                fill_status="ambiguous_paths",
                gaf_path_support=selected_support,
                gaf_best_alt_support=alt_support,
                hic_path_support=selected_hic_support,
                hic_best_alt_support=alt_hic_support,
                reason="multiple_candidate_paths",
            )
        if gaf_choice is not None:
            selected_index = gaf_choice
            graph_status = "gaf_resolved_paths"
        else:
            selected_index = hic_choice
            graph_status = "hic_resolved_paths"

    path = paths[selected_index]
    selected_support = supports[selected_index] if supports else None
    alt_supports = [
        support for index, support in enumerate(supports)
        if index != selected_index
    ]
    alt_support = max(alt_supports, default=0) if supports else None
    selected_hic_support = hic_supports[selected_index] if hic_supports else None
    alt_hic_supports = [
        support for index, support in enumerate(hic_supports)
        if index != selected_index
    ]
    alt_hic_support = max(alt_hic_supports, default=0) if hic_supports else None
    path_nodes = graph_path_nodes(path)
    intermediate_nodes = graph_intermediate_nodes(path)

    fill_sequence, right_trim_bp, reason = reconstruct_fill_sequence(
        graph,
        path,
        left,
        right,
        args.max_fill_bp,
    )
    fillable = fill_sequence is not None
    return FillPlan(
        scaffold=scaffold_name,
        left_contig=left.assignment.new_name,
        right_contig=right.assignment.new_name,
        left_graph_node=left_node,
        right_graph_node=right_node,
        left_orientation=left_orientation,
        right_orientation=right_orientation,
        raw_inferred_gap_bp=raw_gap,
        fallback_gap_bp=fallback_gap_bp,
        overlap_bp=overlap_bp,
        overlap_class=overlap_class,
        graph_status=graph_status,
        path_edges=len(path),
        path_nodes=path_nodes,
        intermediate_nodes=intermediate_nodes,
        candidate_paths=len(paths),
        fill_status="fillable" if fillable else reason,
        gaf_path_support=selected_support,
        gaf_best_alt_support=alt_support,
        hic_path_support=selected_hic_support,
        hic_best_alt_support=alt_hic_support,
        fill_sequence=fill_sequence or "",
        fill_bp=len(fill_sequence or ""),
        right_trim_bp=right_trim_bp,
        reason=reason,
    )


def build_fill_plans(groups, graph, args, gaf_records=None, hic_contacts=None):
    plans = []
    for scaffold_name, members in groups.items():
        for left, right in zip(members, members[1:]):
            plans.append(
                plan_gap(
                    scaffold_name,
                    left,
                    right,
                    graph,
                    args,
                    gaf_records,
                    hic_contacts,
                )
            )
    return plans


def plans_by_scaffold(plans):
    grouped = defaultdict(list)
    for plan in plans:
        grouped[plan.scaffold].append(plan)
    return dict(grouped)


def filled_header(record, simple_headers):
    if simple_headers:
        return record.name
    return (
        f"{record.name} contigs={record.contigs} filled_gaps={record.filled_gaps} "
        f"fallback_gaps={record.fallback_gaps} fill_bp={record.fill_bp} "
        f"fallback_gap_bp={record.fallback_gap_bp} trimmed_bp={record.trimmed_bp}"
    )


def build_filled_scaffolds(groups, unassigned, plans, args):
    grouped_plans = plans_by_scaffold(plans)
    records = []
    for scaffold_name, members in groups.items():
        pieces = [members[0].record.seq] if members else []
        filled_gaps = 0
        fallback_gaps = 0
        fill_bp = 0
        fallback_gap_bp = 0
        trimmed_bp = 0
        for plan, member in zip(grouped_plans.get(scaffold_name, []), members[1:]):
            if plan.fill_status == "fillable":
                plan.applied = True
                pieces.append(plan.fill_sequence)
                pieces.append(member.record.seq[plan.right_trim_bp:])
                filled_gaps += 1
                fill_bp += plan.fill_bp
                trimmed_bp += plan.right_trim_bp
            else:
                pieces.append("N" * plan.fallback_gap_bp)
                pieces.append(member.record.seq)
                fallback_gaps += 1
                fallback_gap_bp += plan.fallback_gap_bp

        records.append(
            FilledScaffold(
                name=scaffold_name,
                seq="".join(pieces),
                contigs=len(members),
                filled_gaps=filled_gaps,
                fallback_gaps=fallback_gaps,
                fill_bp=fill_bp,
                fallback_gap_bp=fallback_gap_bp,
                trimmed_bp=trimmed_bp,
            )
        )

    for record in unassigned:
        records.append(
            FilledScaffold(
                name=record.name,
                seq=record.seq,
                contigs=1,
                filled_gaps=0,
                fallback_gaps=0,
                fill_bp=0,
                fallback_gap_bp=0,
                trimmed_bp=0,
            )
        )
    return records


def write_fill_plan(path, plans, include_sequences):
    header = [
        "scaffold",
        "left_contig",
        "right_contig",
        "left_graph_node",
        "right_graph_node",
        "left_orientation",
        "right_orientation",
        "raw_inferred_gap_bp",
        "fallback_gap_bp",
        "overlap_bp",
        "overlap_class",
        "graph_status",
        "path_edges",
        "path_nodes",
        "intermediate_nodes",
        "candidate_paths",
        "gaf_path_support",
        "gaf_best_alt_support",
        "hic_path_support",
        "hic_best_alt_support",
        "fill_status",
        "fill_bp",
        "right_trim_bp",
        "applied",
        "reason",
        "fill_sequence",
    ]
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for plan in plans:
            row = [
                plan.scaffold,
                plan.left_contig,
                plan.right_contig,
                plan.left_graph_node,
                plan.right_graph_node,
                plan.left_orientation,
                plan.right_orientation,
                plan.raw_inferred_gap_bp,
                plan.fallback_gap_bp,
                plan.overlap_bp,
                plan.overlap_class,
                plan.graph_status,
                plan.path_edges if plan.path_edges is not None else ".",
                plan.path_nodes,
                plan.intermediate_nodes,
                plan.candidate_paths,
                plan.gaf_path_support if plan.gaf_path_support is not None else ".",
                plan.gaf_best_alt_support if plan.gaf_best_alt_support is not None else ".",
                plan.hic_path_support if plan.hic_path_support is not None else ".",
                plan.hic_best_alt_support if plan.hic_best_alt_support is not None else ".",
                plan.fill_status,
                plan.fill_bp,
                plan.right_trim_bp,
                "yes" if plan.applied else "no",
                plan.reason,
                plan.fill_sequence if include_sequences and plan.fill_sequence else ".",
            ]
            out.write("\t".join(str(item) for item in row) + "\n")


def write_filled_fasta(path, records, simple_headers):
    with open(path, "w") as out:
        for record in records:
            out.write(f">{filled_header(record, simple_headers)}\n")
            write_wrapped(out, record.seq)


def write_run_summary(path, args, output_paths, plans, filled_records):
    fillable = sum(1 for plan in plans if plan.fill_status == "fillable")
    applied = sum(1 for plan in plans if plan.applied)
    with open(path, "w") as out:
        out.write("chromo fill\n")
        out.write("\nInputs\n")
        out.write(f"ordered_fasta\t{args.ordered_fasta}\n")
        out.write(f"assignments\t{args.assignments}\n")
        out.write(f"gfa\t{args.gfa}\n")
        out.write(f"gaf\t{args.gaf if args.gaf else '.'}\n")
        out.write(f"hic_pairs\t{args.hic_pairs if args.hic_pairs else '.'}\n")
        out.write("\nParameters\n")
        out.write(f"apply\t{args.apply}\n")
        out.write(f"fixed_gap_bp\t{args.fixed_gap_bp if args.fixed_gap_bp is not None else '.'}\n")
        out.write(f"max_path_edges\t{args.max_path_edges}\n")
        out.write(f"max_candidate_paths\t{args.max_candidate_paths}\n")
        out.write(f"max_fill_bp\t{args.max_fill_bp}\n")
        out.write(f"min_gaf_mapq\t{args.min_gaf_mapq}\n")
        out.write(f"min_gaf_path_support\t{args.min_gaf_path_support}\n")
        out.write(f"min_hic_path_support\t{args.min_hic_path_support}\n")
        out.write("\nOutputs\n")
        for label, output_path in output_paths.items():
            out.write(f"{label}\t{output_path}\n")
        out.write("\nSummary\n")
        out.write(f"planned_gaps\t{len(plans)}\n")
        out.write(f"fillable_gaps\t{fillable}\n")
        out.write(f"applied_gaps\t{applied}\n")
        for status, count in sorted(status_counts(plans).items()):
            out.write(f"status_{status}\t{count}\n")
        if filled_records is not None:
            out.write(f"filled_records\t{len(filled_records)}\n")
            out.write(f"filled_sequence_bp\t{sum(record.fill_bp for record in filled_records)}\n")
            out.write(f"fallback_gap_bp\t{sum(record.fallback_gap_bp for record in filled_records)}\n")
            out.write(f"trimmed_bp\t{sum(record.trimmed_bp for record in filled_records)}\n")


def status_counts(plans):
    counts = defaultdict(int)
    for plan in plans:
        counts[plan.fill_status] += 1
    return dict(counts)


def run(args):
    if args.fixed_gap_bp is not None and args.fixed_gap_bp < 0:
        raise ValueError("--fixed-gap-bp must be zero or greater")
    if args.max_path_edges < 1:
        raise ValueError("--max-path-edges must be at least 1")
    if args.max_candidate_paths < 1:
        raise ValueError("--max-candidate-paths must be at least 1")
    if args.min_gaf_mapq < 0:
        raise ValueError("--min-gaf-mapq must be zero or greater")
    if args.min_gaf_path_support < 1:
        raise ValueError("--min-gaf-path-support must be at least 1")
    if args.min_hic_path_support < 1:
        raise ValueError("--min-hic-path-support must be at least 1")

    prefix = Path(args.output_prefix)
    if prefix.parent and str(prefix.parent) != ".":
        prefix.parent.mkdir(parents=True, exist_ok=True)

    output_paths = {
        "fill_plan": Path(str(prefix) + ".fill_plan.tsv"),
        "run_summary": Path(str(prefix) + ".run_summary.txt"),
    }
    if args.apply:
        output_paths["filled_fasta"] = Path(str(prefix) + ".filled.fa")
    for output_path in output_paths.values():
        if output_path.parent and str(output_path.parent) != ".":
            output_path.parent.mkdir(parents=True, exist_ok=True)

    graph = read_gfa(args.gfa)
    gaf_records = read_gaf(args.gaf, args.min_gaf_mapq) if args.gaf else []
    hic_contacts = read_hic_pairs(args.hic_pairs) if args.hic_pairs else {}
    assignments = read_assignments(args.assignments)
    records = read_ordered_fasta(args.ordered_fasta)
    groups, unassigned = group_scaffold_members(records, assignments)
    plans = build_fill_plans(groups, graph, args, gaf_records, hic_contacts)

    filled_records = None
    if args.apply:
        filled_records = build_filled_scaffolds(groups, unassigned, plans, args)
        write_filled_fasta(output_paths["filled_fasta"], filled_records, args.simple_headers)

    write_fill_plan(output_paths["fill_plan"], plans, args.include_fill_sequences)
    write_run_summary(output_paths["run_summary"], args, output_paths, plans, filled_records)

    sys.stderr.write(f"Planned {len(plans)} graph gap fill(s).\n")
    for status, count in sorted(status_counts(plans).items()):
        sys.stderr.write(f"  {status}: {count}\n")
    sys.stderr.write(f"Wrote fill plan: {output_paths['fill_plan']}\n")
    if args.apply:
        sys.stderr.write(f"Wrote filled FASTA: {output_paths['filled_fasta']}\n")


def main(argv: Optional[Sequence[str]] = None, prog: Optional[str] = None):
    args = parse_args(argv, prog=prog)
    try:
        run(args)
    except ValueError as exc:
        sys.stderr.write(f"ERROR: {exc}\n")
        sys.exit(2)


if __name__ == "__main__":
    main()
