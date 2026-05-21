#!/usr/bin/env python3
"""Plan and optionally apply reviewed graph-based gap fills."""

import argparse
import csv
import json
import sys
from collections import defaultdict, deque
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

from .graph import read_gfa
from .reference_order import iter_paf, reverse_complement, write_wrapped
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
    ref_path_support: Optional[int] = None
    ref_best_alt_support: Optional[int] = None
    risk_flags: str = "."
    branch_complexity_score: int = 0
    high_degree_nodes: str = "."
    self_loop_nodes: str = "."
    unsequenced_nodes: str = "."
    cycles_avoided: int = 0
    fill_sequence: str = ""
    fill_bp: int = 0
    right_trim_bp: int = 0
    accept_fill: bool = False
    applied: bool = False
    reason: str = "."
    candidate_details: list = field(default_factory=list)


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


@dataclass(frozen=True)
class GraphPathSearchResult:
    paths: list
    cycles_avoided: int = 0


def parse_args(argv: Optional[Sequence[str]] = None, prog: Optional[str] = None):
    ap = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Plan graph-supported gap fills between adjacent ChromoSort sorted "
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
        "--ref-paf",
        default=None,
        help=(
            "Optional reference-to-assembly PAF. Candidate graph paths can use "
            "intermediate node reference placement to resolve an ambiguous "
            "branch when one path has unique expected-gap support."
        ),
    )
    ap.add_argument(
        "-o",
        "--output-prefix",
        required=True,
        help=(
            "Output prefix. Writes <prefix>.gapfill_plan.tsv and "
            "<prefix>.run_summary.txt; with --apply also writes "
            "<prefix>.gapfilled.fa."
        ),
    )
    ap.add_argument(
        "--apply",
        action="store_true",
        help="Write <prefix>.gapfilled.fa using only fillable graph paths.",
    )
    ap.add_argument(
        "--reviewed-plan",
        default=None,
        help=(
            "Optional edited gapfill plan TSV. When provided with --apply, only "
            "rows with accept_fill=yes are applied after the current graph "
            "path is rechecked."
        ),
    )
    ap.add_argument(
        "--review-html",
        default=None,
        help=(
            "Optional self-contained HTML review dashboard for the gapfill plan. "
            "The dashboard can export an edited reviewed-plan TSV."
        ),
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
        "--min-ref-path-support",
        type=int,
        default=1,
        help=(
            "Minimum expected-gap reference-placement support required to "
            "resolve an ambiguous graph branch from --ref-paf."
        ),
    )
    ap.add_argument(
        "--min-ref-paf-mapq",
        type=int,
        default=0,
        help="Minimum MAPQ for reference-placement PAF rows used by --ref-paf.",
    )
    ap.add_argument(
        "--min-ref-paf-idy",
        type=float,
        default=0.0,
        help="Minimum percent identity for reference-placement PAF rows used by --ref-paf.",
    )
    ap.add_argument(
        "--include-secondary-ref-paf",
        action="store_true",
        help="Include secondary PAF rows marked tp:A:S when reading --ref-paf.",
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
        help="Write gapfilled FASTA headers containing only the scaffold ID.",
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


def selected_and_alt_support(supports, selected_index):
    if not supports:
        return None, None
    selected = supports[selected_index]
    alternatives = [
        support for index, support in enumerate(supports)
        if index != selected_index
    ]
    return selected, max(alternatives, default=0)


def support_conflict_reason(choices):
    labels = {label for label, _index in choices}
    if labels == {"gaf", "hic"}:
        return "conflicting_gaf_hic_support"
    if "ref_paf" in labels:
        return "conflicting_ref_paf_support"
    return "conflicting_path_support"


def support_graph_status(choices):
    labels = [label for label, _index in choices]
    if "gaf" in labels:
        return "gaf_resolved_paths"
    if "hic" in labels:
        return "hic_resolved_paths"
    return "ref_paf_resolved_paths"


def read_ref_paf(path, min_identity=0.0, min_mapq=0, include_secondary=False):
    placements = defaultdict(list)
    for segment in iter_paf(
        path,
        min_identity=min_identity,
        min_mapq=min_mapq,
        include_secondary=include_secondary,
    ):
        placements[segment.query].append(segment)
    return dict(placements)


def expected_ref_gap_window(left, right):
    left_edge = left.assignment.ref_end
    right_edge = right.assignment.ref_start
    return min(left_edge, right_edge), max(left_edge, right_edge)


def segment_midpoint(segment):
    return (segment.ref_start + segment.ref_end) / 2.0


def ref_segment_support(segment):
    return int(round(segment.len_query * (segment.identity / 100.0)))


def ref_placement_support_for_node(placements, node, orientation, ref_name, ref_window):
    best = 0
    for segment in placements.get(node, []):
        if segment.ref != ref_name:
            continue
        if segment.orientation != orientation:
            continue
        if not ref_window[0] <= segment_midpoint(segment) <= ref_window[1]:
            continue
        best = max(best, ref_segment_support(segment))
    return best


def ref_support_for_path(ref_placements, path, ref_name, ref_window):
    support = 0
    for node, orientation in path_oriented_nodes(path)[1:-1]:
        support += ref_placement_support_for_node(
            ref_placements,
            node,
            orientation,
            ref_name,
            ref_window,
        )
    return support


def ref_path_supports(ref_placements, paths, left, right, scaffold_name):
    if not ref_placements:
        return []
    if left.assignment.ref != right.assignment.ref:
        return [0 for _path in paths]
    ref_window = expected_ref_gap_window(left, right)
    return [
        ref_support_for_path(ref_placements, path, scaffold_name, ref_window)
        for path in paths
    ]


def ordered_unique(values):
    seen = set()
    ordered = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return ordered


def join_values(values):
    values = list(values)
    return ",".join(values) if values else "."


def graph_node_degree(graph, node):
    return len(graph.incoming(node)) + len(graph.outgoing(node))


def graph_node_has_self_loop(graph, node):
    return any(edge.source == edge.target == node for edge in graph.outgoing(node))


def risk_nodes_for_paths(paths):
    return ordered_unique(
        node
        for path in paths
        for node, _orientation in path_oriented_nodes(path)
    )


def support_risk_flags(gaf_supports, hic_supports, ref_supports, args):
    flags = []
    for label, supports, threshold in (
        ("gaf", gaf_supports, args.min_gaf_path_support),
        ("hic", hic_supports, args.min_hic_path_support),
        ("ref_paf", ref_supports, args.min_ref_path_support),
    ):
        if not supports:
            continue
        best = max(supports)
        if best < threshold:
            flags.append(f"weak_{label}_support")
        elif supports.count(best) > 1:
            flags.append(f"tied_{label}_support")
    return flags


def path_risk_annotations(
    graph,
    paths,
    candidate_paths,
    cycles_avoided=0,
    extra_flags=None,
):
    nodes = risk_nodes_for_paths(paths)
    high_degree_nodes = [
        node for node in nodes
        if graph_node_degree(graph, node) > 2
    ]
    self_loop_nodes = [
        node for node in nodes
        if graph_node_has_self_loop(graph, node)
    ]
    unsequenced_nodes = [
        node for node in nodes
        if graph.nodes[node].sequence is None
    ]

    flags = []
    if candidate_paths > 1:
        flags.append("branching")
    if high_degree_nodes:
        flags.append("high_degree")
    if self_loop_nodes:
        flags.append("self_loop")
    if unsequenced_nodes:
        flags.append("unsequenced")
    if cycles_avoided:
        flags.append("cycles_avoided")
    flags.extend(extra_flags or [])
    flags = ordered_unique(flags)

    branch_score = max(0, candidate_paths - 1)
    branch_score += len(high_degree_nodes)
    branch_score += 2 * len(self_loop_nodes)
    branch_score += 2 * len(unsequenced_nodes)
    branch_score += cycles_avoided
    branch_score += len(extra_flags or [])

    return {
        "risk_flags": join_values(flags),
        "branch_complexity_score": branch_score,
        "high_degree_nodes": join_values(high_degree_nodes),
        "self_loop_nodes": join_values(self_loop_nodes),
        "unsequenced_nodes": join_values(unsequenced_nodes),
        "cycles_avoided": cycles_avoided,
    }


def support_value(supports, index):
    if not supports:
        return "."
    return supports[index]


def candidate_path_details(
    graph,
    paths,
    selected_index,
    left,
    right,
    args,
    gaf_supports=None,
    hic_supports=None,
    ref_supports=None,
    cycles_avoided=0,
):
    details = []
    gaf_supports = gaf_supports or []
    hic_supports = hic_supports or []
    ref_supports = ref_supports or []
    for index, path in enumerate(paths):
        fill_sequence, right_trim_bp, reason = reconstruct_fill_sequence(
            graph,
            path,
            left,
            right,
            args.max_fill_bp,
        )
        risk = path_risk_annotations(
            graph,
            [path],
            len(paths),
            cycles_avoided,
            extra_flags=[] if fill_sequence is not None else ["sequence_validation_failed"],
        )
        details.append(
            {
                "candidate_index": index + 1,
                "reported": "yes" if index == selected_index else "no",
                "path_nodes": graph_path_nodes(path),
                "intermediate_nodes": graph_intermediate_nodes(path),
                "path_edges": len(path),
                "gaf_support": support_value(gaf_supports, index),
                "hic_support": support_value(hic_supports, index),
                "ref_support": support_value(ref_supports, index),
                "candidate_status": "fillable" if fill_sequence is not None else reason,
                "fill_bp": len(fill_sequence or ""),
                "right_trim_bp": right_trim_bp,
                "risk_flags": risk["risk_flags"],
                "branch_complexity_score": risk["branch_complexity_score"],
                "fill_sequence": fill_sequence or ".",
            }
        )
    return details


def parse_review_accept(value, line_number):
    normalized = (value or "").strip().lower()
    if normalized in {"yes", "y", "true", "1", "accept", "accepted", "apply"}:
        return True
    if normalized in {"", ".", "no", "n", "false", "0", "reject", "skip"}:
        return False
    raise ValueError(
        f"Malformed accept_fill value at reviewed plan line {line_number}: "
        f"{value!r}"
    )


def fill_plan_key(plan):
    return (plan.scaffold, plan.left_contig, plan.right_contig)


def read_reviewed_plan(path):
    required = {"scaffold", "left_contig", "right_contig", "path_nodes", "accept_fill"}
    decisions = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = set(reader.fieldnames or [])
        missing = sorted(required - fieldnames)
        if missing:
            raise ValueError(
                "Reviewed gapfill plan is missing required column(s): "
                + ", ".join(missing)
            )
        for line_number, row in enumerate(reader, start=2):
            key = (
                row.get("scaffold", ""),
                row.get("left_contig", ""),
                row.get("right_contig", ""),
            )
            if not all(key):
                raise ValueError(
                    f"Malformed reviewed plan row at line {line_number}: "
                    "missing scaffold, left_contig, or right_contig"
                )
            if key in decisions:
                raise ValueError(
                    f"Duplicate reviewed plan row at line {line_number}: "
                    f"{key[0]} {key[1]} {key[2]}"
                )
            decisions[key] = {
                "accept": parse_review_accept(row.get("accept_fill", ""), line_number),
                "path_nodes": row.get("path_nodes", ""),
                "line_number": line_number,
            }
    return decisions


def apply_reviewed_plan(plans, decisions):
    plan_by_key = {fill_plan_key(plan): plan for plan in plans}
    for key, decision in decisions.items():
        if key not in plan_by_key and decision["accept"]:
            raise ValueError(
                "Reviewed plan accepts a gap that is not in the current plan "
                f"at line {decision['line_number']}: {key[0]} {key[1]} {key[2]}"
            )

    for plan in plans:
        decision = decisions.get(fill_plan_key(plan))
        if not decision or not decision["accept"]:
            continue
        if decision["path_nodes"] != plan.path_nodes:
            raise ValueError(
                "Reviewed plan accepted path no longer matches current graph "
                f"path for {plan.scaffold} {plan.left_contig}->{plan.right_contig}: "
                f"{decision['path_nodes']} != {plan.path_nodes}"
            )
        if plan.fill_status != "fillable":
            raise ValueError(
                "Reviewed plan accepts a gap that is not currently fillable "
                f"for {plan.scaffold} {plan.left_contig}->{plan.right_contig}: "
                f"{plan.fill_status}"
            )
        plan.accept_fill = True


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
        return GraphPathSearchResult([])

    targets = graph_target_keys(right_node, right_orientation)
    if not targets:
        return GraphPathSearchResult([])

    queue = deque()
    for start in graph_start_keys(left_node, left_orientation):
        queue.append((start, [], {start}))

    paths = []
    cycles_avoided = 0
    while queue and len(paths) < max_paths:
        (node, orientation), path, seen = queue.popleft()
        if len(path) >= max_edges:
            continue
        for edge in graph.outgoing(node, orientation):
            next_key = edge.target_key
            if next_key in seen:
                cycles_avoided += 1
                continue
            next_path = path + [edge]
            if next_key in targets:
                paths.append(next_path)
                if len(paths) >= max_paths:
                    break
                continue
            queue.append((next_key, next_path, seen | {next_key}))
    return GraphPathSearchResult(paths, cycles_avoided)


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
    ref_placements=None,
):
    gaf_records = gaf_records or []
    hic_contacts = hic_contacts or {}
    ref_placements = ref_placements or {}
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

    path_search = enumerate_graph_paths(
        graph,
        left_node,
        right_node,
        left_orientation,
        right_orientation,
        args.max_path_edges,
        args.max_candidate_paths,
    )
    paths = path_search.paths
    if not paths:
        risk = path_risk_annotations(
            graph,
            [],
            0,
            path_search.cycles_avoided,
        )
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
            **risk,
        )

    supports = gaf_path_supports(gaf_records, paths) if gaf_records else []
    hic_supports = hic_path_supports(hic_contacts, paths) if hic_contacts else []
    ref_supports = ref_path_supports(
        ref_placements,
        paths,
        left,
        right,
        scaffold_name,
    )
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
        ref_choice = choose_unique_supported_path(
            ref_supports,
            args.min_ref_path_support,
        )
        support_choices = [
            (label, choice)
            for label, choice in (
                ("gaf", gaf_choice),
                ("hic", hic_choice),
                ("ref_paf", ref_choice),
            )
            if choice is not None
        ]
        if len({choice for _label, choice in support_choices}) > 1:
            path = paths[0]
            risk = path_risk_annotations(
                graph,
                paths,
                len(paths),
                path_search.cycles_avoided,
                extra_flags=["conflicting_support"],
            )
            selected_support, alt_support = selected_and_alt_support(supports, 0)
            selected_hic_support, alt_hic_support = selected_and_alt_support(
                hic_supports,
                0,
            )
            selected_ref_support, alt_ref_support = selected_and_alt_support(
                ref_supports,
                0,
            )
            candidates = candidate_path_details(
                graph,
                paths,
                0,
                left,
                right,
                args,
                supports,
                hic_supports,
                ref_supports,
                path_search.cycles_avoided,
            )
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
                ref_path_support=selected_ref_support,
                ref_best_alt_support=alt_ref_support,
                reason=support_conflict_reason(support_choices),
                candidate_details=candidates,
                **risk,
            )
        if not support_choices:
            path = paths[0]
            risk = path_risk_annotations(
                graph,
                paths,
                len(paths),
                path_search.cycles_avoided,
                extra_flags=support_risk_flags(
                    supports,
                    hic_supports,
                    ref_supports,
                    args,
                ),
            )
            selected_support, alt_support = selected_and_alt_support(supports, 0)
            selected_hic_support, alt_hic_support = selected_and_alt_support(
                hic_supports,
                0,
            )
            selected_ref_support, alt_ref_support = selected_and_alt_support(
                ref_supports,
                0,
            )
            candidates = candidate_path_details(
                graph,
                paths,
                0,
                left,
                right,
                args,
                supports,
                hic_supports,
                ref_supports,
                path_search.cycles_avoided,
            )
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
                ref_path_support=selected_ref_support,
                ref_best_alt_support=alt_ref_support,
                reason="multiple_candidate_paths",
                candidate_details=candidates,
                **risk,
            )
        selected_index = support_choices[0][1]
        graph_status = support_graph_status(support_choices)

    path = paths[selected_index]
    selected_support, alt_support = selected_and_alt_support(supports, selected_index)
    selected_hic_support, alt_hic_support = selected_and_alt_support(
        hic_supports,
        selected_index,
    )
    selected_ref_support, alt_ref_support = selected_and_alt_support(
        ref_supports,
        selected_index,
    )
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
    risk = path_risk_annotations(
        graph,
        [path],
        len(paths),
        path_search.cycles_avoided,
        extra_flags=[] if fillable else ["sequence_validation_failed"],
    )
    candidates = candidate_path_details(
        graph,
        paths,
        selected_index,
        left,
        right,
        args,
        supports,
        hic_supports,
        ref_supports,
        path_search.cycles_avoided,
    )
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
        ref_path_support=selected_ref_support,
        ref_best_alt_support=alt_ref_support,
        fill_sequence=fill_sequence or "",
        fill_bp=len(fill_sequence or ""),
        right_trim_bp=right_trim_bp,
        reason=reason,
        candidate_details=candidates,
        **risk,
    )


def build_fill_plans(
    groups,
    graph,
    args,
    gaf_records=None,
    hic_contacts=None,
    ref_placements=None,
):
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
                    ref_placements,
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
            if plan.fill_status == "fillable" and (
                not args.reviewed_plan or plan.accept_fill
            ):
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


def fill_plan_header():
    return [
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
        "ref_path_support",
        "ref_best_alt_support",
        "risk_flags",
        "branch_complexity_score",
        "high_degree_nodes",
        "self_loop_nodes",
        "unsequenced_nodes",
        "cycles_avoided",
        "fill_status",
        "fill_bp",
        "right_trim_bp",
        "accept_fill",
        "applied",
        "reason",
        "fill_sequence",
    ]


def fill_plan_row(plan, include_sequences):
    return [
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
        plan.ref_path_support if plan.ref_path_support is not None else ".",
        plan.ref_best_alt_support if plan.ref_best_alt_support is not None else ".",
        plan.risk_flags,
        plan.branch_complexity_score,
        plan.high_degree_nodes,
        plan.self_loop_nodes,
        plan.unsequenced_nodes,
        plan.cycles_avoided,
        plan.fill_status,
        plan.fill_bp,
        plan.right_trim_bp,
        "yes" if plan.accept_fill else "no",
        "yes" if plan.applied else "no",
        plan.reason,
        plan.fill_sequence if include_sequences and plan.fill_sequence else ".",
    ]


def fill_plan_dict(plan, include_sequences):
    return {
        column: str(value)
        for column, value in zip(fill_plan_header(), fill_plan_row(plan, include_sequences))
    }


def candidate_detail_columns():
    return [
        "candidate_index",
        "reported",
        "path_nodes",
        "intermediate_nodes",
        "path_edges",
        "gaf_support",
        "hic_support",
        "ref_support",
        "candidate_status",
        "fill_bp",
        "right_trim_bp",
        "risk_flags",
        "branch_complexity_score",
        "fill_sequence",
    ]


def candidate_detail_for_review(detail, include_sequences):
    row = {
        column: str(detail.get(column, "."))
        for column in candidate_detail_columns()
    }
    if not include_sequences:
        row["fill_sequence"] = "."
    return row


def fill_plan_review_dict(plan, include_sequences):
    row = fill_plan_dict(plan, include_sequences)
    row["_candidate_paths"] = [
        candidate_detail_for_review(detail, include_sequences)
        for detail in plan.candidate_details
    ]
    return row


def write_fill_plan(path, plans, include_sequences):
    header = fill_plan_header()
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for plan in plans:
            row = fill_plan_row(plan, include_sequences)
            out.write("\t".join(str(item) for item in row) + "\n")


def json_for_script(value):
    return (
        json.dumps(value, sort_keys=True)
        .replace("&", "\\u0026")
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
    )


def write_review_html(path, plans, include_sequences):
    data = {
        "schema": "chromosort-gapfill-review-v1",
        "columns": fill_plan_header(),
        "candidateColumns": candidate_detail_columns(),
        "rows": [fill_plan_review_dict(plan, include_sequences) for plan in plans],
    }
    html_text = GAPFILL_REVIEW_HTML.replace(
        "__CHROMOSORT_GAPFILL_REVIEW_DATA__",
        json_for_script(data),
    )
    with open(path, "w") as out:
        out.write(html_text)


def write_filled_fasta(path, records, simple_headers):
    with open(path, "w") as out:
        for record in records:
            out.write(f">{filled_header(record, simple_headers)}\n")
            write_wrapped(out, record.seq)


def write_run_summary(path, args, output_paths, plans, filled_records):
    fillable = sum(1 for plan in plans if plan.fill_status == "fillable")
    applied = sum(1 for plan in plans if plan.applied)
    with open(path, "w") as out:
        out.write("chromo gapfill\n")
        out.write("\nInputs\n")
        out.write(f"ordered_fasta\t{args.ordered_fasta}\n")
        out.write(f"assignments\t{args.assignments}\n")
        out.write(f"gfa\t{args.gfa}\n")
        out.write(f"gaf\t{args.gaf if args.gaf else '.'}\n")
        out.write(f"hic_pairs\t{args.hic_pairs if args.hic_pairs else '.'}\n")
        out.write(f"ref_paf\t{args.ref_paf if args.ref_paf else '.'}\n")
        out.write(f"reviewed_plan\t{args.reviewed_plan if args.reviewed_plan else '.'}\n")
        out.write("\nParameters\n")
        out.write(f"apply\t{args.apply}\n")
        out.write(f"fixed_gap_bp\t{args.fixed_gap_bp if args.fixed_gap_bp is not None else '.'}\n")
        out.write(f"max_path_edges\t{args.max_path_edges}\n")
        out.write(f"max_candidate_paths\t{args.max_candidate_paths}\n")
        out.write(f"max_fill_bp\t{args.max_fill_bp}\n")
        out.write(f"min_gaf_mapq\t{args.min_gaf_mapq}\n")
        out.write(f"min_gaf_path_support\t{args.min_gaf_path_support}\n")
        out.write(f"min_hic_path_support\t{args.min_hic_path_support}\n")
        out.write(f"min_ref_path_support\t{args.min_ref_path_support}\n")
        out.write(f"min_ref_paf_mapq\t{args.min_ref_paf_mapq}\n")
        out.write(f"min_ref_paf_idy\t{args.min_ref_paf_idy}\n")
        out.write(f"include_secondary_ref_paf\t{args.include_secondary_ref_paf}\n")
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
    if args.min_ref_path_support < 1:
        raise ValueError("--min-ref-path-support must be at least 1")
    if args.min_ref_paf_mapq < 0:
        raise ValueError("--min-ref-paf-mapq must be zero or greater")
    if args.min_ref_paf_idy < 0:
        raise ValueError("--min-ref-paf-idy must be zero or greater")

    prefix = Path(args.output_prefix)
    if prefix.parent and str(prefix.parent) != ".":
        prefix.parent.mkdir(parents=True, exist_ok=True)

    output_paths = {
        "fill_plan": Path(str(prefix) + ".gapfill_plan.tsv"),
        "run_summary": Path(str(prefix) + ".run_summary.txt"),
    }
    if args.apply:
        output_paths["filled_fasta"] = Path(str(prefix) + ".gapfilled.fa")
    if args.review_html:
        output_paths["review_html"] = Path(args.review_html)
    for output_path in output_paths.values():
        if output_path.parent and str(output_path.parent) != ".":
            output_path.parent.mkdir(parents=True, exist_ok=True)

    graph = read_gfa(args.gfa)
    gaf_records = read_gaf(args.gaf, args.min_gaf_mapq) if args.gaf else []
    hic_contacts = read_hic_pairs(args.hic_pairs) if args.hic_pairs else {}
    ref_placements = (
        read_ref_paf(
            args.ref_paf,
            min_identity=args.min_ref_paf_idy,
            min_mapq=args.min_ref_paf_mapq,
            include_secondary=args.include_secondary_ref_paf,
        )
        if args.ref_paf
        else {}
    )
    assignments = read_assignments(args.assignments)
    records = read_ordered_fasta(args.ordered_fasta)
    groups, unassigned = group_scaffold_members(records, assignments)
    plans = build_fill_plans(
        groups,
        graph,
        args,
        gaf_records,
        hic_contacts,
        ref_placements,
    )
    if args.reviewed_plan:
        apply_reviewed_plan(plans, read_reviewed_plan(args.reviewed_plan))

    filled_records = None
    if args.apply:
        filled_records = build_filled_scaffolds(groups, unassigned, plans, args)
        write_filled_fasta(output_paths["filled_fasta"], filled_records, args.simple_headers)

    write_fill_plan(output_paths["fill_plan"], plans, args.include_fill_sequences)
    if args.review_html:
        write_review_html(output_paths["review_html"], plans, args.include_fill_sequences)
    write_run_summary(output_paths["run_summary"], args, output_paths, plans, filled_records)

    sys.stderr.write(f"Planned {len(plans)} graph gap fill(s).\n")
    for status, count in sorted(status_counts(plans).items()):
        sys.stderr.write(f"  {status}: {count}\n")
    sys.stderr.write(f"Wrote gapfill plan: {output_paths['fill_plan']}\n")
    if args.review_html:
        sys.stderr.write(f"Wrote review HTML: {output_paths['review_html']}\n")
    if args.apply:
        sys.stderr.write(f"Wrote gapfilled FASTA: {output_paths['filled_fasta']}\n")


GAPFILL_REVIEW_HTML = r"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>ChromoSort Gapfill Review</title>
  <style>
    :root {
      color-scheme: light;
      --border: #cbd5e1;
      --ink: #111827;
      --muted: #64748b;
      --panel: #f8fafc;
      --accent: #0f766e;
      --warn: #b45309;
      --bad: #b91c1c;
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      color: var(--ink);
      background: #ffffff;
      font: 14px/1.45 -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
    }
    header {
      padding: 18px 24px;
      border-bottom: 1px solid var(--border);
      background: var(--panel);
    }
    h1 {
      margin: 0;
      font-size: 22px;
      font-weight: 650;
      letter-spacing: 0;
    }
    main { padding: 18px 24px 28px; }
    .toolbar {
      display: flex;
      flex-wrap: wrap;
      gap: 10px;
      align-items: end;
      margin-bottom: 14px;
    }
    label {
      display: grid;
      gap: 4px;
      color: var(--muted);
      font-size: 12px;
      font-weight: 600;
    }
    input[type="search"], select {
      min-width: 190px;
      border: 1px solid var(--border);
      border-radius: 6px;
      padding: 8px 10px;
      color: var(--ink);
      background: #ffffff;
      font: inherit;
    }
    button {
      border: 1px solid #0f766e;
      border-radius: 6px;
      padding: 8px 12px;
      color: #ffffff;
      background: var(--accent);
      font: inherit;
      font-weight: 650;
      cursor: pointer;
    }
    button:hover { background: #115e59; }
    .summary {
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      margin-bottom: 14px;
    }
    .summary span {
      border: 1px solid var(--border);
      border-radius: 999px;
      padding: 4px 9px;
      background: #ffffff;
      color: var(--muted);
      font-size: 12px;
      font-weight: 650;
    }
    .table-wrap {
      overflow: auto;
      border: 1px solid var(--border);
      border-radius: 8px;
    }
    table {
      width: 100%;
      min-width: 1320px;
      border-collapse: collapse;
      background: #ffffff;
    }
    th, td {
      border-bottom: 1px solid #e5e7eb;
      padding: 7px 8px;
      text-align: left;
      vertical-align: top;
      white-space: nowrap;
    }
    th {
      position: sticky;
      top: 0;
      z-index: 1;
      background: #f1f5f9;
      color: #334155;
      font-size: 12px;
      font-weight: 700;
    }
    tr:last-child td { border-bottom: 0; }
    td.path, td.sequence {
      max-width: 360px;
      overflow: hidden;
      text-overflow: ellipsis;
    }
    .status-fillable { color: var(--accent); font-weight: 700; }
    .status-ambiguous_paths { color: var(--warn); font-weight: 700; }
    .status-no_graph_path,
    .status-missing_node,
    .status-left_flank_sequence_mismatch,
    .status-right_flank_sequence_mismatch { color: var(--bad); font-weight: 700; }
    .candidate-cell {
      background: #f8fafc;
      white-space: normal;
    }
    .candidate-title {
      margin-bottom: 6px;
      color: var(--muted);
      font-size: 12px;
      font-weight: 700;
    }
    .candidate-table {
      width: 100%;
      min-width: 980px;
      border-collapse: collapse;
      background: #ffffff;
    }
    .candidate-table th {
      position: static;
      background: #e2e8f0;
    }
    .candidate-table th,
    .candidate-table td {
      padding: 5px 7px;
      font-size: 12px;
      white-space: nowrap;
    }
  </style>
</head>
<body>
  <header>
    <h1>ChromoSort Gapfill Review</h1>
  </header>
  <main>
    <div class="toolbar">
      <label>Search<input id="search" type="search"></label>
      <label>Status<select id="status"></select></label>
      <button id="export">Export reviewed TSV</button>
    </div>
    <div class="summary" id="summary"></div>
    <div class="table-wrap">
      <table>
        <thead id="thead"></thead>
        <tbody id="tbody"></tbody>
      </table>
    </div>
  </main>
  <script id="chromosort-gapfill-review-data" type="application/json">__CHROMOSORT_GAPFILL_REVIEW_DATA__</script>
  <script>
    const data = JSON.parse(document.getElementById("chromosort-gapfill-review-data").textContent);
    const rows = data.rows.map(row => ({...row}));
    const columns = data.columns;
    const candidateColumns = data.candidateColumns || [];
    const displayColumns = [
      "accept_fill", "scaffold", "left_contig", "right_contig", "fill_status",
      "graph_status", "path_nodes", "gaf_path_support", "gaf_best_alt_support",
      "hic_path_support", "hic_best_alt_support", "ref_path_support",
      "ref_best_alt_support", "risk_flags", "branch_complexity_score",
      "high_degree_nodes", "self_loop_nodes", "unsequenced_nodes",
      "cycles_avoided", "fill_bp", "right_trim_bp", "fallback_gap_bp", "reason",
      "fill_sequence"
    ];
    const els = {
      search: document.getElementById("search"),
      status: document.getElementById("status"),
      export: document.getElementById("export"),
      summary: document.getElementById("summary"),
      thead: document.getElementById("thead"),
      tbody: document.getElementById("tbody")
    };

    function cls(value) {
      return String(value || "").replace(/[^a-zA-Z0-9_-]/g, "_");
    }

    function tsvCell(value) {
      return String(value ?? ".").replace(/[\t\r\n]+/g, " ");
    }

    function statusOptions() {
      const statuses = Array.from(new Set(rows.map(row => row.fill_status))).sort();
      els.status.replaceChildren();
      const all = document.createElement("option");
      all.value = "";
      all.textContent = "All";
      els.status.appendChild(all);
      for (const status of statuses) {
        const option = document.createElement("option");
        option.value = status;
        option.textContent = status;
        els.status.appendChild(option);
      }
    }

    function visibleRows() {
      const needle = els.search.value.trim().toLowerCase();
      const status = els.status.value;
      return rows.filter(row => {
        if (status && row.fill_status !== status) return false;
        if (!needle) return true;
        const rowHit = displayColumns.some(column => String(row[column] || "").toLowerCase().includes(needle));
        const candidateHit = (row._candidate_paths || []).some(candidate =>
          candidateColumns.some(column => String(candidate[column] || "").toLowerCase().includes(needle))
        );
        return rowHit || candidateHit;
      });
    }

    function renderSummary() {
      const fillable = rows.filter(row => row.fill_status === "fillable").length;
      const accepted = rows.filter(row => row.accept_fill === "yes").length;
      const applied = rows.filter(row => row.applied === "yes").length;
      const visible = visibleRows().length;
      const items = [
        `${rows.length} planned`,
        `${visible} visible`,
        `${fillable} fillable`,
        `${accepted} accepted`,
        `${applied} applied`
      ];
      els.summary.replaceChildren();
      for (const item of items) {
        const node = document.createElement("span");
        node.textContent = item;
        els.summary.appendChild(node);
      }
    }

    function renderHeader() {
      const tr = document.createElement("tr");
      for (const column of displayColumns) {
        const th = document.createElement("th");
        th.textContent = column;
        tr.appendChild(th);
      }
      els.thead.replaceChildren(tr);
    }

    function renderTable() {
      els.tbody.replaceChildren();
      for (const row of visibleRows()) {
        const tr = document.createElement("tr");
        for (const column of displayColumns) {
          const td = document.createElement("td");
          if (column === "accept_fill") {
            const checkbox = document.createElement("input");
            checkbox.type = "checkbox";
            checkbox.checked = row.accept_fill === "yes";
            checkbox.disabled = row.fill_status !== "fillable";
            checkbox.addEventListener("change", () => {
              row.accept_fill = checkbox.checked ? "yes" : "no";
              renderSummary();
            });
            td.appendChild(checkbox);
          } else {
            td.textContent = row[column] ?? ".";
          }
          if (column === "fill_status") td.className = `status-${cls(row.fill_status)}`;
          if (column === "path_nodes") td.className = "path";
          if (column === "fill_sequence") td.className = "sequence";
          tr.appendChild(td);
        }
        els.tbody.appendChild(tr);
        if ((row._candidate_paths || []).length > 1) {
          renderCandidateRows(row);
        }
      }
      renderSummary();
    }

    function renderCandidateRows(row) {
      const tr = document.createElement("tr");
      const td = document.createElement("td");
      td.colSpan = displayColumns.length;
      td.className = "candidate-cell";

      const title = document.createElement("div");
      title.className = "candidate-title";
      title.textContent = `Candidate path comparison for ${row.left_contig} -> ${row.right_contig}`;
      td.appendChild(title);

      const table = document.createElement("table");
      table.className = "candidate-table";
      const thead = document.createElement("thead");
      const headRow = document.createElement("tr");
      for (const column of candidateColumns) {
        const th = document.createElement("th");
        th.textContent = column;
        headRow.appendChild(th);
      }
      thead.appendChild(headRow);
      table.appendChild(thead);

      const tbody = document.createElement("tbody");
      for (const candidate of row._candidate_paths || []) {
        const candidateRow = document.createElement("tr");
        for (const column of candidateColumns) {
          const cell = document.createElement("td");
          cell.textContent = candidate[column] ?? ".";
          if (column === "path_nodes") cell.className = "path";
          if (column === "fill_sequence") cell.className = "sequence";
          candidateRow.appendChild(cell);
        }
        tbody.appendChild(candidateRow);
      }
      table.appendChild(tbody);
      td.appendChild(table);
      tr.appendChild(td);
      els.tbody.appendChild(tr);
    }

    function exportTsv() {
      const lines = [columns.join("\t")];
      for (const row of rows) {
        lines.push(columns.map(column => tsvCell(row[column])).join("\t"));
      }
      const blob = new Blob([lines.join("\n") + "\n"], {type: "text/tab-separated-values"});
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = "chromosort.gapfill.reviewed_plan.tsv";
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
    }

    statusOptions();
    renderHeader();
    renderTable();
    els.search.addEventListener("input", renderTable);
    els.status.addEventListener("change", renderTable);
    els.export.addEventListener("click", exportTsv);
  </script>
</body>
</html>
"""


def main(argv: Optional[Sequence[str]] = None, prog: Optional[str] = None):
    args = parse_args(argv, prog=prog)
    try:
        run(args)
    except ValueError as exc:
        sys.stderr.write(f"ERROR: {exc}\n")
        sys.exit(2)


if __name__ == "__main__":
    main()
