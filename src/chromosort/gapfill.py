#!/usr/bin/env python3
"""Plan and optionally apply reviewed graph-based gap fills."""

import argparse
import csv
import json
import sys
from collections import OrderedDict, defaultdict, deque
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

from .agp import (
    component_part,
    gap_part,
    group_parts_by_object,
    read_agp,
    write_agp,
    write_component_tsv,
)
from .gaf import (
    choose_gaf_supported_path,
    path_oriented_nodes,
    read_gaf,
    summarize_gaf_traversal,
)
from .graph import build_path_projection, projections_by_contig, read_gfa
from .longreads import read_long_read_paf, summarize_contig_bridge
from .paths import ensure_output_dirs, ensure_parent_dir
from .reference_order import iter_paf, reverse_complement, write_wrapped
from .review import read_review_events
from .scaffold import (
    AssignmentRow,
    FastaRecord,
    ScaffoldMember,
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
from .submission import write_submission_checklist


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
    gaf_support_status: str = "."
    gaf_selected_reads: str = "."
    hic_path_support: Optional[int] = None
    hic_best_alt_support: Optional[int] = None
    ref_path_support: Optional[int] = None
    ref_best_alt_support: Optional[int] = None
    longread_bridge_reads: Optional[int] = None
    longread_orientation_summary: str = "."
    longread_read_order_summary: str = "."
    longread_median_read_gap_bp: str = "."
    patch_candidate_count: int = 0
    patch_best_id: str = "."
    patch_best_source: str = "."
    patch_best_bp: Optional[int] = None
    patch_graph_status: str = "."
    patch_best_notes: str = "."
    patch_best_sequence: str = ""
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
    graph_path: list = field(default_factory=list)


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
class GraphPathSearchResult:
    paths: list
    cycles_avoided: int = 0


@dataclass(frozen=True)
class ProjectedTerminal:
    status: str
    node: str = "."
    orientation: str = "."
    path_name: str = "."
    step_index: Optional[int] = None
    side: str = "."


@dataclass(frozen=True)
class PatchCandidate:
    scaffold: str
    left_contig: str
    right_contig: str
    patch_id: str
    source: str
    sequence: str
    notes: str = "."


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
        default=None,
        help=(
            "Final ordered FASTA from chromo sort, optionally after chromo fix. "
            "Use together with --assignments."
        ),
    )
    ap.add_argument(
        "-a",
        "--assignments",
        default=None,
        help=(
            "Corresponding <prefix>.contig_assignments.tsv from chromo sort. "
            "Use together with --ordered-fasta."
        ),
    )
    ap.add_argument(
        "--scaffold-fasta",
        default=None,
        help=(
            "Existing scaffold FASTA whose component and N-gap layout is "
            "described by --agp. This is the post-scaffold input mode."
        ),
    )
    ap.add_argument(
        "--agp",
        default=None,
        help=(
            "AGP 2.1 component/gap map for --scaffold-fasta. Required when "
            "using post-scaffold gapfill mode."
        ),
    )
    ap.add_argument(
        "--gfa",
        required=True,
        help="Assembly graph GFA containing segment sequences and links.",
    )
    ap.add_argument(
        "--project-gfa-paths",
        action="store_true",
        help=(
            "Project ordered or AGP component names through matching GFA P/W "
            "path records, then plan from terminal unitigs. Sequence-bearing "
            "projected paths can be applied; no-sequence GFAs remain topology-only."
        ),
    )
    ap.add_argument(
        "--projection-trim-overlaps",
        action="store_true",
        help=(
            "When --project-gfa-paths is used, subtract path-step overlaps "
            "while building the GFA path projection."
        ),
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
        "--read-paf",
        default=None,
        help=(
            "Optional long-read-to-contig/ordered-FASTA PAF. Used as bridge "
            "evidence across adjacent contigs; it does not supply inserted "
            "gapfill sequence."
        ),
    )
    ap.add_argument(
        "--min-read-mapq",
        type=int,
        default=0,
        help="Minimum MAPQ for --read-paf rows.",
    )
    ap.add_argument(
        "--read-terminal-window-bp",
        type=int,
        default=5_000,
        help="Terminal window used when counting long-read bridges.",
    )
    ap.add_argument(
        "--read-min-anchor-bp",
        type=int,
        default=500,
        help="Minimum terminal anchor length required for a long-read bridge.",
    )
    ap.add_argument(
        "--patch-table",
        default=None,
        help=(
            "Optional TSV of external patch candidates keyed by scaffold, "
            "left_contig, and right_contig. Rows may include patch_sequence "
            "directly or patch_id plus --patch-fasta."
        ),
    )
    ap.add_argument(
        "--patch-fasta",
        default=None,
        help=(
            "Optional FASTA containing external patch sequences referenced by "
            "--patch-table patch_id values."
        ),
    )
    ap.add_argument(
        "-o",
        "--output-prefix",
        required=True,
        help=(
            "Output prefix. Writes <prefix>.gapfill_plan.tsv, "
            "<prefix>.gapfilled.agp, <prefix>.gapfilled_components.tsv, "
            "<prefix>.submission_checklist.tsv, and <prefix>.run_summary.txt; "
            "with --apply also writes <prefix>.gapfilled.fa."
        ),
    )
    ap.add_argument(
        "--apply",
        action="store_true",
        help=(
            "Write <prefix>.gapfilled.fa. Requires either --reviewed-plan or "
            "--apply-all-fillable."
        ),
    )
    ap.add_argument(
        "--apply-all-fillable",
        action="store_true",
        help=(
            "With --apply, apply every currently fillable graph path without a "
            "reviewed plan. Intended for deliberate exploratory or benchmark runs."
        ),
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


def gaf_summary_for(gaf_records, paths, selected_index, args):
    if not gaf_records:
        return None
    return summarize_gaf_traversal(
        gaf_records,
        paths,
        selected_index=selected_index,
        min_support=args.min_gaf_path_support,
    )


def gaf_summary_supports(summary):
    if summary is None:
        return []
    return [support.support for support in summary.path_supports]


def gaf_selected_reads(summary):
    if summary is None or not summary.selected_reads:
        return "."
    return ",".join(summary.selected_reads)


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
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if "schema" in set(reader.fieldnames or []):
            return read_review_event_plan(path)

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


def read_review_event_plan(path):
    decisions = {}
    for event in read_review_events(path, expected_task="gapfill"):
        if event.action != "fill_path":
            raise ValueError(
                f"Review event {event.event_id!r} has unsupported gapfill action "
                f"{event.action!r}."
            )
        key = (
            event.fields.get("scaffold", ""),
            event.fields.get("left_contig", ""),
            event.fields.get("right_contig", ""),
        )
        if not all(key):
            raise ValueError(
                f"Malformed gapfill review event {event.event_id!r}: "
                "missing scaffold, left_contig, or right_contig"
            )
        if key in decisions:
            raise ValueError(
                f"Duplicate gapfill review event {event.event_id!r}: "
                f"{key[0]} {key[1]} {key[2]}"
            )
        decisions[key] = {
            "accept": event.accept,
            "path_nodes": event.fields.get("path_nodes", ""),
            "line_number": event.event_id,
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


def sequence_has_suffix(observed, expected):
    return observed.upper().endswith(expected.upper())


def sequence_has_prefix(observed, expected):
    return observed.upper().startswith(expected.upper())


def projected_fill_status(fill_sequence, reason):
    if fill_sequence is not None:
        return "fillable"
    if reason and reason.startswith("unsequenced"):
        return "projected_path_planning_only"
    return reason


def reconstruct_projected_fill_sequence(graph, path, left, right, max_fill_bp):
    if not path:
        return None, 0, "no_graph_path"

    left_seq = oriented_sequence(
        graph.nodes[path[0].source].sequence,
        path[0].source_orientation,
    )
    right_seq = oriented_sequence(
        graph.nodes[path[-1].target].sequence,
        path[-1].target_orientation,
    )
    if left_seq is None or right_seq is None:
        return None, 0, "unsequenced_flank_node"
    if not sequence_has_suffix(left.record.seq, left_seq):
        return None, 0, "left_terminal_sequence_mismatch"
    if not sequence_has_prefix(right.record.seq, right_seq):
        return None, 0, "right_terminal_sequence_mismatch"

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
    if right_trim_bp > len(right_seq) or right_trim_bp > len(right.record.seq):
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

    gaf_summary = gaf_summary_for(gaf_records, paths, 0, args)
    supports = gaf_summary_supports(gaf_summary)
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
                gaf_support_status=gaf_summary.support_status if gaf_summary else ".",
                gaf_selected_reads=gaf_selected_reads(gaf_summary),
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
                gaf_support_status=gaf_summary.support_status if gaf_summary else ".",
                gaf_selected_reads=gaf_selected_reads(gaf_summary),
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
    selected_gaf_summary = gaf_summary_for(gaf_records, paths, selected_index, args)
    selected_support, alt_support = selected_and_alt_support(
        gaf_summary_supports(selected_gaf_summary),
        selected_index,
    )
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
        gaf_support_status=selected_gaf_summary.support_status if selected_gaf_summary else ".",
        gaf_selected_reads=gaf_selected_reads(selected_gaf_summary),
        hic_path_support=selected_hic_support,
        hic_best_alt_support=alt_hic_support,
        ref_path_support=selected_ref_support,
        ref_best_alt_support=alt_ref_support,
        fill_sequence=fill_sequence or "",
        fill_bp=len(fill_sequence or ""),
        right_trim_bp=right_trim_bp,
        reason=reason,
        candidate_details=candidates,
        graph_path=path,
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


AGP_ORIENTATIONS = {"+", "-", "?", "0", "na", "."}


def part_length(part):
    return part.object_end - part.object_start + 1


def records_by_name(records, path):
    by_name = OrderedDict()
    for record in records:
        if record.name in by_name:
            raise ValueError(f"Duplicate FASTA record {record.name!r} in {path}")
        by_name[record.name] = record
    return by_name


def validate_agp_against_scaffold(records, agp_parts, scaffold_fasta):
    if not agp_parts:
        raise ValueError("--agp did not contain any component or gap rows")

    record_map = records_by_name(records, scaffold_fasta)
    grouped_parts = group_parts_by_object(agp_parts)

    missing_records = [name for name in record_map if name not in grouped_parts]
    if missing_records:
        preview = ", ".join(missing_records[:5])
        suffix = "..." if len(missing_records) > 5 else ""
        raise ValueError(
            "Scaffold FASTA record(s) missing from AGP: "
            f"{preview}{suffix}."
        )

    unknown_objects = [name for name in grouped_parts if name not in record_map]
    if unknown_objects:
        preview = ", ".join(unknown_objects[:5])
        suffix = "..." if len(unknown_objects) > 5 else ""
        raise ValueError(
            "AGP object(s) missing from scaffold FASTA: "
            f"{preview}{suffix}."
        )

    for object_name, parts in grouped_parts.items():
        if not parts:
            continue
        if parts[0].part_type != "component" or parts[-1].part_type != "component":
            raise ValueError(
                f"AGP object {object_name!r} must start and end with component rows "
                "for scaffold/AGP gapfill mode."
            )

        record = record_map[object_name]
        cursor = 1
        seen_components = set()
        for expected_part_number, part in enumerate(parts, start=1):
            if part.part_number != expected_part_number:
                raise ValueError(
                    f"AGP object {object_name!r} has non-sequential part number "
                    f"{part.part_number}; expected {expected_part_number}."
                )
            if part.object_start != cursor:
                raise ValueError(
                    f"AGP object {object_name!r} is not contiguous at part "
                    f"{part.part_number}: expected object_start {cursor}, "
                    f"found {part.object_start}."
                )

            length = part_length(part)
            if part.part_type == "gap":
                if part.gap_length != length:
                    raise ValueError(
                        f"AGP object {object_name!r} gap part {part.part_number} "
                        f"has gap_length {part.gap_length}, but object span is {length}."
                    )
                observed = record.seq[part.object_start - 1 : part.object_end]
                if observed.upper() != "N" * length:
                    raise ValueError(
                        f"AGP object {object_name!r} gap part {part.part_number} "
                        "does not correspond to an N-only span in --scaffold-fasta."
                    )
            else:
                if not part.component_id or part.component_id == ".":
                    raise ValueError(
                        f"AGP object {object_name!r} component part {part.part_number} "
                        "has an empty component ID."
                    )
                if part.component_id in seen_components:
                    raise ValueError(
                        f"AGP object {object_name!r} uses component ID "
                        f"{part.component_id!r} more than once; direct scaffold/AGP "
                        "gapfill currently requires unique component IDs per object."
                    )
                seen_components.add(part.component_id)
                component_length = part.component_end - part.component_start + 1
                if component_length != length:
                    raise ValueError(
                        f"AGP object {object_name!r} component part {part.part_number} "
                        f"has component span {component_length}, but object span is {length}."
                    )
                if (part.orientation or ".").lower() not in AGP_ORIENTATIONS:
                    raise ValueError(
                        f"AGP object {object_name!r} component part {part.part_number} "
                        f"has unsupported orientation {part.orientation!r}."
                    )
            cursor = part.object_end + 1

        if cursor != len(record.seq) + 1:
            raise ValueError(
                f"AGP object {object_name!r} covers {cursor - 1} bp, but scaffold "
                f"FASTA record length is {len(record.seq)} bp."
            )

    return grouped_parts, record_map


def agp_orientation_for_graph(part):
    return part.orientation if part.orientation in {"+", "-"} else "."


def member_from_agp_part(part, scaffold_seq):
    seq = scaffold_seq[part.object_start - 1 : part.object_end]
    assignment = AssignmentRow(
        contig=part.component_id,
        new_name=part.component_id,
        ref=part.object_name,
        order_in_ref=part.part_number,
        ref_start=part.object_start,
        ref_end=part.object_end,
        orientation=agp_orientation_for_graph(part),
    )
    member = ScaffoldMember(
        assignment=assignment,
        record=FastaRecord(
            name=part.component_id,
            header=f">{part.component_id}",
            seq=seq,
        ),
    )
    member.agp_component_id = part.component_id
    member.agp_component_start = part.component_start
    member.agp_component_end = part.component_end
    member.agp_component_orientation = part.orientation
    member.agp_component_type = part.component_type
    member.agp_source = "agp_component"
    member.agp_notes = f"{part.object_name}:{part.object_start}-{part.object_end}"
    return member


def build_scaffold_agp_groups(records, agp_parts, scaffold_fasta):
    grouped_parts, record_map = validate_agp_against_scaffold(
        records,
        agp_parts,
        scaffold_fasta,
    )
    groups = OrderedDict()
    gap_parts_by_pair = {}

    for object_name, parts in grouped_parts.items():
        members = []
        pending_gaps = []
        previous_member = None
        scaffold_seq = record_map[object_name].seq

        for part in parts:
            if part.part_type == "gap":
                pending_gaps.append(part)
                continue

            member = member_from_agp_part(part, scaffold_seq)
            if previous_member is not None:
                key = (
                    object_name,
                    previous_member.assignment.new_name,
                    member.assignment.new_name,
                )
                gap_parts_by_pair[key] = tuple(pending_gaps)
            pending_gaps = []
            members.append(member)
            previous_member = member

        groups[object_name] = members

    return groups, gap_parts_by_pair


def attach_agp_gap_metadata(plan, gap_parts):
    if not gap_parts:
        return plan
    first = gap_parts[0]
    plan.agp_gap_type = first.gap_type
    plan.agp_linkage = first.linkage
    plan.agp_linkage_evidence = first.linkage_evidence
    plan.agp_gap_part_count = len(gap_parts)
    return plan


def blocked_agp_junction_plan(
    scaffold_name,
    left,
    right,
    graph,
    reason,
    fallback_gap_bp=0,
):
    raw_gap = inferred_gap(left, right)
    overlap_bp = max(0, -raw_gap)
    left_node = graph_node_name(left, graph)
    right_node = graph_node_name(right, graph)
    return FillPlan(
        scaffold=scaffold_name,
        left_contig=left.assignment.new_name,
        right_contig=right.assignment.new_name,
        left_graph_node=left_node,
        right_graph_node=right_node,
        left_orientation=graph_orientation(left),
        right_orientation=graph_orientation(right),
        raw_inferred_gap_bp=raw_gap,
        fallback_gap_bp=fallback_gap_bp,
        overlap_bp=overlap_bp,
        overlap_class=classify_adjacent_overlap(left, right, overlap_bp),
        graph_status=reason,
        path_edges=None,
        path_nodes=".",
        intermediate_nodes=".",
        candidate_paths=0,
        fill_status=reason,
        reason=reason,
    )


def build_agp_fill_plans(
    groups,
    gap_parts_by_pair,
    graph,
    args,
    gaf_records=None,
    hic_contacts=None,
    ref_placements=None,
):
    plans = []
    for scaffold_name, members in groups.items():
        for left, right in zip(members, members[1:]):
            key = (scaffold_name, left.assignment.new_name, right.assignment.new_name)
            gap_parts = gap_parts_by_pair.get(key, ())
            if not gap_parts:
                plans.append(
                    blocked_agp_junction_plan(
                        scaffold_name,
                        left,
                        right,
                        graph,
                        "not_agp_gap",
                        fallback_gap_bp=0,
                    )
                )
                continue
            plan = plan_gap(
                scaffold_name,
                left,
                right,
                graph,
                args,
                gaf_records,
                hic_contacts,
                ref_placements,
            )
            plans.append(attach_agp_gap_metadata(plan, gap_parts))
    return plans


def member_name_candidates(member):
    candidates = [
        member.assignment.new_name,
        member.assignment.contig,
        member.record.name,
    ]
    ordered = []
    for candidate in candidates:
        if candidate and candidate != "." and candidate not in ordered:
            ordered.append(candidate)
    return ordered


def requested_projection_names(groups):
    names = []
    for members in groups.values():
        for member in members:
            for candidate in member_name_candidates(member):
                if candidate not in names:
                    names.append(candidate)
    return names


def flip_orientation(orientation):
    if orientation == "+":
        return "-"
    if orientation == "-":
        return "+"
    return orientation


def member_projection_rows(member, projections):
    for candidate in member_name_candidates(member):
        rows = projections.get(candidate)
        if rows:
            return candidate, sorted(rows, key=lambda row: row.step_index)
    return ".", []


def projection_terminal_for_member(member, projections, side):
    path_name, rows = member_projection_rows(member, projections)
    if not rows:
        return ProjectedTerminal(status="missing_projection", side=side)

    line_numbers = {
        row.path_line_number for row in rows
        if row.path_line_number is not None
    }
    if len(line_numbers) > 1:
        return ProjectedTerminal(
            status="ambiguous_projection",
            path_name=path_name,
            side=side,
        )

    member_orientation = graph_orientation(member)
    use_forward_order = member_orientation != "-"
    if side == "right":
        row = rows[-1] if use_forward_order else rows[0]
    else:
        row = rows[0] if use_forward_order else rows[-1]

    orientation = row.segment_orientation
    if member_orientation == "-":
        orientation = flip_orientation(orientation)
    return ProjectedTerminal(
        status="present",
        node=row.segment,
        orientation=orientation,
        path_name=path_name,
        step_index=row.step_index,
        side=side,
    )


def projected_terminal_status(left_terminal, right_terminal):
    left_missing = left_terminal.status != "present"
    right_missing = right_terminal.status != "present"
    if left_missing and right_missing:
        if left_terminal.status == right_terminal.status:
            return f"both_{left_terminal.status}s"
        return "missing_or_ambiguous_projections"
    if left_missing:
        return f"left_{left_terminal.status}"
    if right_missing:
        return f"right_{right_terminal.status}"
    return None


def projected_candidate_path_details(
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
        fill_sequence, right_trim_bp, reason = reconstruct_projected_fill_sequence(
            graph,
            path,
            left,
            right,
            args.max_fill_bp,
        )
        candidate_status = projected_fill_status(fill_sequence, reason)
        extra_flags = ["projected_unitig"]
        if fill_sequence is None and candidate_status != "projected_path_planning_only":
            extra_flags.append("sequence_validation_failed")
        risk = path_risk_annotations(
            graph,
            [path],
            len(paths),
            cycles_avoided,
            extra_flags=extra_flags,
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
                "candidate_status": candidate_status,
                "fill_bp": len(fill_sequence or ""),
                "right_trim_bp": right_trim_bp,
                "risk_flags": risk["risk_flags"],
                "branch_complexity_score": risk["branch_complexity_score"],
                "fill_sequence": fill_sequence or ".",
            }
        )
    return details


def projected_projection_reason(left_terminal, right_terminal):
    parts = []
    for label, terminal in (("left", left_terminal), ("right", right_terminal)):
        if terminal.status != "present":
            parts.append(f"{label}:{terminal.status}")
        else:
            parts.append(
                f"{label}:{terminal.path_name}:{terminal.step_index}:{terminal.node}{terminal.orientation}"
            )
    return ";".join(parts)


def plan_projected_gap(
    scaffold_name,
    left,
    right,
    graph,
    projections,
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
    left_terminal = projection_terminal_for_member(left, projections, "right")
    right_terminal = projection_terminal_for_member(right, projections, "left")
    projection_status = projected_terminal_status(left_terminal, right_terminal)

    if projection_status is not None:
        return FillPlan(
            scaffold=scaffold_name,
            left_contig=left.assignment.new_name,
            right_contig=right.assignment.new_name,
            left_graph_node=left_terminal.node,
            right_graph_node=right_terminal.node,
            left_orientation=left_terminal.orientation,
            right_orientation=right_terminal.orientation,
            raw_inferred_gap_bp=raw_gap,
            fallback_gap_bp=fallback_gap_bp,
            overlap_bp=overlap_bp,
            overlap_class=overlap_class,
            graph_status=projection_status,
            path_edges=None,
            path_nodes=".",
            intermediate_nodes=".",
            candidate_paths=0,
            fill_status="missing_projection",
            risk_flags="projected_unitig",
            branch_complexity_score=1,
            reason=projected_projection_reason(left_terminal, right_terminal),
        )

    missing_status = graph_status_for_nodes(left_terminal.node, right_terminal.node)
    if missing_status is not None:
        return FillPlan(
            scaffold=scaffold_name,
            left_contig=left.assignment.new_name,
            right_contig=right.assignment.new_name,
            left_graph_node=left_terminal.node,
            right_graph_node=right_terminal.node,
            left_orientation=left_terminal.orientation,
            right_orientation=right_terminal.orientation,
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
            risk_flags="projected_unitig",
            branch_complexity_score=1,
            reason=missing_status,
        )

    path_search = enumerate_graph_paths(
        graph,
        left_terminal.node,
        right_terminal.node,
        left_terminal.orientation,
        right_terminal.orientation,
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
            extra_flags=["projected_unitig"],
        )
        return FillPlan(
            scaffold=scaffold_name,
            left_contig=left.assignment.new_name,
            right_contig=right.assignment.new_name,
            left_graph_node=left_terminal.node,
            right_graph_node=right_terminal.node,
            left_orientation=left_terminal.orientation,
            right_orientation=right_terminal.orientation,
            raw_inferred_gap_bp=raw_gap,
            fallback_gap_bp=fallback_gap_bp,
            overlap_bp=overlap_bp,
            overlap_class=overlap_class,
            graph_status="projected_no_path",
            path_edges=None,
            path_nodes=".",
            intermediate_nodes=".",
            candidate_paths=0,
            fill_status="no_graph_path",
            reason="projected_no_graph_path",
            **risk,
        )

    gaf_summary = gaf_summary_for(gaf_records, paths, 0, args)
    supports = gaf_summary_supports(gaf_summary)
    hic_supports = hic_path_supports(hic_contacts, paths) if hic_contacts else []
    ref_supports = ref_path_supports(
        ref_placements,
        paths,
        left,
        right,
        scaffold_name,
    )
    selected_index = 0
    graph_status = "projected_direct_edge" if len(paths[0]) == 1 else "projected_short_path"

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
            risk = path_risk_annotations(
                graph,
                paths,
                len(paths),
                path_search.cycles_avoided,
                extra_flags=["projected_unitig", "conflicting_support"],
            )
            selected_support, alt_support = selected_and_alt_support(supports, 0)
            selected_hic_support, alt_hic_support = selected_and_alt_support(hic_supports, 0)
            selected_ref_support, alt_ref_support = selected_and_alt_support(ref_supports, 0)
            candidates = projected_candidate_path_details(
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
                left_graph_node=left_terminal.node,
                right_graph_node=right_terminal.node,
                left_orientation=left_terminal.orientation,
                right_orientation=right_terminal.orientation,
                raw_inferred_gap_bp=raw_gap,
                fallback_gap_bp=fallback_gap_bp,
                overlap_bp=overlap_bp,
                overlap_class=overlap_class,
                graph_status="projected_ambiguous_paths",
                path_edges=len(paths[0]),
                path_nodes=graph_path_nodes(paths[0]),
                intermediate_nodes=graph_intermediate_nodes(paths[0]),
                candidate_paths=len(paths),
                fill_status="ambiguous_paths",
                gaf_path_support=selected_support,
                gaf_best_alt_support=alt_support,
                gaf_support_status=gaf_summary.support_status if gaf_summary else ".",
                gaf_selected_reads=gaf_selected_reads(gaf_summary),
                hic_path_support=selected_hic_support,
                hic_best_alt_support=alt_hic_support,
                ref_path_support=selected_ref_support,
                ref_best_alt_support=alt_ref_support,
                reason=support_conflict_reason(support_choices),
                candidate_details=candidates,
                **risk,
            )
        if not support_choices:
            risk = path_risk_annotations(
                graph,
                paths,
                len(paths),
                path_search.cycles_avoided,
                extra_flags=["projected_unitig"]
                + support_risk_flags(supports, hic_supports, ref_supports, args),
            )
            selected_support, alt_support = selected_and_alt_support(supports, 0)
            selected_hic_support, alt_hic_support = selected_and_alt_support(hic_supports, 0)
            selected_ref_support, alt_ref_support = selected_and_alt_support(ref_supports, 0)
            candidates = projected_candidate_path_details(
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
                left_graph_node=left_terminal.node,
                right_graph_node=right_terminal.node,
                left_orientation=left_terminal.orientation,
                right_orientation=right_terminal.orientation,
                raw_inferred_gap_bp=raw_gap,
                fallback_gap_bp=fallback_gap_bp,
                overlap_bp=overlap_bp,
                overlap_class=overlap_class,
                graph_status="projected_ambiguous_paths",
                path_edges=len(paths[0]),
                path_nodes=graph_path_nodes(paths[0]),
                intermediate_nodes=graph_intermediate_nodes(paths[0]),
                candidate_paths=len(paths),
                fill_status="ambiguous_paths",
                gaf_path_support=selected_support,
                gaf_best_alt_support=alt_support,
                gaf_support_status=gaf_summary.support_status if gaf_summary else ".",
                gaf_selected_reads=gaf_selected_reads(gaf_summary),
                hic_path_support=selected_hic_support,
                hic_best_alt_support=alt_hic_support,
                ref_path_support=selected_ref_support,
                ref_best_alt_support=alt_ref_support,
                reason="multiple_projected_candidate_paths",
                candidate_details=candidates,
                **risk,
            )
        selected_index = support_choices[0][1]
        graph_status = "projected_" + support_graph_status(support_choices)

    path = paths[selected_index]
    selected_gaf_summary = gaf_summary_for(gaf_records, paths, selected_index, args)
    selected_support, alt_support = selected_and_alt_support(
        gaf_summary_supports(selected_gaf_summary),
        selected_index,
    )
    selected_hic_support, alt_hic_support = selected_and_alt_support(
        hic_supports,
        selected_index,
    )
    selected_ref_support, alt_ref_support = selected_and_alt_support(
        ref_supports,
        selected_index,
    )
    fill_sequence, right_trim_bp, reason = reconstruct_projected_fill_sequence(
        graph,
        path,
        left,
        right,
        args.max_fill_bp,
    )
    fill_status = projected_fill_status(fill_sequence, reason)
    extra_flags = ["projected_unitig"]
    if fill_sequence is None and fill_status != "projected_path_planning_only":
        extra_flags.append("sequence_validation_failed")
    risk = path_risk_annotations(
        graph,
        [path],
        len(paths),
        path_search.cycles_avoided,
        extra_flags=extra_flags,
    )
    candidates = projected_candidate_path_details(
        graph,
        paths,
        selected_index,
        left,
        right,
        args,
        gaf_summary_supports(selected_gaf_summary),
        hic_supports,
        ref_supports,
        path_search.cycles_avoided,
    )
    return FillPlan(
        scaffold=scaffold_name,
        left_contig=left.assignment.new_name,
        right_contig=right.assignment.new_name,
        left_graph_node=left_terminal.node,
        right_graph_node=right_terminal.node,
        left_orientation=left_terminal.orientation,
        right_orientation=right_terminal.orientation,
        raw_inferred_gap_bp=raw_gap,
        fallback_gap_bp=fallback_gap_bp,
        overlap_bp=overlap_bp,
        overlap_class=overlap_class,
        graph_status=graph_status,
        path_edges=len(path),
        path_nodes=graph_path_nodes(path),
        intermediate_nodes=graph_intermediate_nodes(path),
        candidate_paths=len(paths),
        fill_status=fill_status,
        gaf_path_support=selected_support,
        gaf_best_alt_support=alt_support,
        gaf_support_status=selected_gaf_summary.support_status if selected_gaf_summary else ".",
        gaf_selected_reads=gaf_selected_reads(selected_gaf_summary),
        hic_path_support=selected_hic_support,
        hic_best_alt_support=alt_hic_support,
        ref_path_support=selected_ref_support,
        ref_best_alt_support=alt_ref_support,
        fill_sequence=fill_sequence or "",
        fill_bp=len(fill_sequence or ""),
        right_trim_bp=right_trim_bp,
        reason=reason,
        candidate_details=candidates,
        graph_path=path,
        **risk,
    )


def build_projected_fill_plans(
    groups,
    graph,
    projections,
    args,
    gaf_records=None,
    hic_contacts=None,
    ref_placements=None,
):
    plans = []
    for scaffold_name, members in groups.items():
        for left, right in zip(members, members[1:]):
            plans.append(
                plan_projected_gap(
                    scaffold_name,
                    left,
                    right,
                    graph,
                    projections,
                    args,
                    gaf_records,
                    hic_contacts,
                    ref_placements,
                )
            )
    return plans


def build_projected_agp_fill_plans(
    groups,
    gap_parts_by_pair,
    graph,
    projections,
    args,
    gaf_records=None,
    hic_contacts=None,
    ref_placements=None,
):
    plans = []
    for scaffold_name, members in groups.items():
        for left, right in zip(members, members[1:]):
            key = (scaffold_name, left.assignment.new_name, right.assignment.new_name)
            gap_parts = gap_parts_by_pair.get(key, ())
            if not gap_parts:
                plans.append(
                    blocked_agp_junction_plan(
                        scaffold_name,
                        left,
                        right,
                        graph,
                        "not_agp_gap",
                        fallback_gap_bp=0,
                    )
                )
                continue
            plan = plan_projected_gap(
                scaffold_name,
                left,
                right,
                graph,
                projections,
                args,
                gaf_records,
                hic_contacts,
                ref_placements,
            )
            plans.append(attach_agp_gap_metadata(plan, gap_parts))
    return plans


def best_read_bridge(read_evidence, left_member, right_member, args):
    supports = []
    for left_name in member_name_candidates(left_member):
        for right_name in member_name_candidates(right_member):
            supports.append(
                summarize_contig_bridge(
                    read_evidence,
                    left_name,
                    right_name,
                    terminal_window_bp=args.read_terminal_window_bp,
                    min_anchor_bp=args.read_min_anchor_bp,
                )
            )
    return max(supports, key=lambda item: item.bridge_count) if supports else None


def annotate_longread_bridge_support(groups, plans, read_evidence, args):
    if read_evidence is None:
        return
    pair_by_key = {}
    for scaffold_name, members in groups.items():
        for left, right in zip(members, members[1:]):
            pair_by_key[(scaffold_name, left.assignment.new_name, right.assignment.new_name)] = (
                left,
                right,
            )

    for plan in plans:
        pair = pair_by_key.get(fill_plan_key(plan))
        if pair is None:
            plan.longread_bridge_reads = 0
            continue
        support = best_read_bridge(read_evidence, pair[0], pair[1], args)
        if support is None:
            plan.longread_bridge_reads = 0
            continue
        plan.longread_bridge_reads = support.bridge_count
        plan.longread_orientation_summary = support.orientation_summary
        plan.longread_read_order_summary = support.read_order_summary
        median_gap = support.median_read_gap_bp
        plan.longread_median_read_gap_bp = "." if median_gap is None else str(median_gap)


PATCH_FIELD_ALIASES = {
    "scaffold": ("scaffold", "object", "object_name", "chrom", "chromosome"),
    "left_contig": ("left_contig", "left_component", "left", "component_left"),
    "right_contig": ("right_contig", "right_component", "right", "component_right"),
    "patch_id": ("patch_id", "id", "name", "record", "patch"),
    "source": ("source", "tool", "method", "caller"),
    "sequence": ("patch_sequence", "sequence", "seq", "fill_sequence"),
    "notes": ("notes", "note", "status", "description"),
}


def row_value(row, aliases, default="."):
    for alias in aliases:
        value = row.get(alias)
        if value is not None and str(value).strip() not in {"", "."}:
            return str(value).strip()
    return default


def read_patch_fasta(path):
    if not path:
        return {}
    records = read_ordered_fasta(path)
    return {record.name: record.seq for record in records}


def normalize_patch_sequence(sequence, row_label):
    sequence = (sequence or "").strip().upper()
    if not sequence or sequence == ".":
        raise ValueError(f"External patch {row_label} is missing a sequence")
    return sequence


def read_patch_candidates(path, fasta_path=None):
    fasta_records = read_patch_fasta(fasta_path)
    candidates = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Patch table {path} is empty or missing a header")
        for line_number, row in enumerate(reader, start=2):
            scaffold = row_value(row, PATCH_FIELD_ALIASES["scaffold"], "")
            left = row_value(row, PATCH_FIELD_ALIASES["left_contig"], "")
            right = row_value(row, PATCH_FIELD_ALIASES["right_contig"], "")
            if not (scaffold and left and right):
                raise ValueError(
                    f"Patch table row {line_number} is missing scaffold, "
                    "left_contig, or right_contig"
                )
            patch_id = row_value(
                row,
                PATCH_FIELD_ALIASES["patch_id"],
                f"patch_row_{line_number}",
            )
            sequence = row_value(row, PATCH_FIELD_ALIASES["sequence"], "")
            if not sequence and patch_id in fasta_records:
                sequence = fasta_records[patch_id]
            sequence = normalize_patch_sequence(sequence, patch_id)
            candidates.append(
                PatchCandidate(
                    scaffold=scaffold,
                    left_contig=left,
                    right_contig=right,
                    patch_id=patch_id,
                    source=row_value(row, PATCH_FIELD_ALIASES["source"], "."),
                    sequence=sequence,
                    notes=row_value(row, PATCH_FIELD_ALIASES["notes"], "."),
                )
            )
    return candidates


def patch_candidate_key(candidate):
    return candidate.scaffold, candidate.left_contig, candidate.right_contig


def patch_candidates_by_key(candidates):
    grouped = defaultdict(list)
    for candidate in candidates:
        grouped[patch_candidate_key(candidate)].append(candidate)
    return dict(grouped)


def patch_compare_status(plan, candidate):
    if not candidate:
        return "."
    if not plan.fill_sequence:
        return "patch_only_no_graph_sequence"
    if candidate.sequence.upper() == plan.fill_sequence.upper():
        return "exact_graph_match"
    if len(candidate.sequence) == len(plan.fill_sequence):
        return "same_length_graph_mismatch"
    return "graph_mismatch"


def patch_sort_key(plan, candidate):
    status = patch_compare_status(plan, candidate)
    status_rank = {
        "exact_graph_match": 0,
        "same_length_graph_mismatch": 1,
        "graph_mismatch": 2,
        "patch_only_no_graph_sequence": 3,
    }.get(status, 4)
    graph_length = len(plan.fill_sequence or "")
    length_delta = abs(len(candidate.sequence) - graph_length) if graph_length else 0
    return status_rank, length_delta, candidate.patch_id


def add_risk_flag(plan, flag):
    if not flag:
        return
    flags = [] if plan.risk_flags in {"", "."} else plan.risk_flags.split(",")
    if flag not in flags:
        flags.append(flag)
        plan.risk_flags = ",".join(flags) if flags else "."
        plan.branch_complexity_score += 1


def annotate_patch_candidates(plans, patch_candidates):
    if not patch_candidates:
        return
    grouped = patch_candidates_by_key(patch_candidates)
    for plan in plans:
        candidates = grouped.get(fill_plan_key(plan), [])
        if not candidates:
            continue
        best = min(candidates, key=lambda candidate: patch_sort_key(plan, candidate))
        status = patch_compare_status(plan, best)
        plan.patch_candidate_count = len(candidates)
        plan.patch_best_id = best.patch_id
        plan.patch_best_source = best.source
        plan.patch_best_bp = len(best.sequence)
        plan.patch_graph_status = status
        plan.patch_best_notes = best.notes
        plan.patch_best_sequence = best.sequence
        if len(candidates) > 1:
            add_risk_flag(plan, "multiple_external_patches")
        if status in {"same_length_graph_mismatch", "graph_mismatch"}:
            add_risk_flag(plan, "external_patch_graph_mismatch")
        elif status == "patch_only_no_graph_sequence":
            add_risk_flag(plan, "external_patch_only")


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
                plan.accept_fill or getattr(args, "apply_all_fillable", False)
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


def add_component_part(parts, object_name, cursor, part_number, component_id, start, end, orientation, source, status=".", notes="."):
    length = end - start + 1
    if length <= 0:
        return cursor, part_number
    parts.append(
        component_part(
            object_name=object_name,
            object_start=cursor,
            object_end=cursor + length - 1,
            part_number=part_number,
            component_id=component_id,
            component_start=start,
            component_end=end,
            orientation=orientation,
            source=source,
            status=status,
            notes=notes,
        )
    )
    return cursor + length, part_number + 1


def add_gap_part(
    parts,
    object_name,
    cursor,
    part_number,
    gap_bp,
    source,
    status=".",
    notes=".",
    gap_type="scaffold",
    linkage="yes",
    linkage_evidence="align_genus",
):
    if gap_bp <= 0:
        return cursor, part_number
    parts.append(
        gap_part(
            object_name=object_name,
            object_start=cursor,
            object_end=cursor + gap_bp - 1,
            part_number=part_number,
            gap_length=gap_bp,
            gap_type=gap_type,
            linkage=linkage,
            linkage_evidence=linkage_evidence,
            source=source,
            status=status,
            notes=notes,
        )
    )
    return cursor + gap_bp, part_number + 1


def component_fields_for_member(member, trim_left_bp=0):
    component_id = getattr(member, "agp_component_id", member.assignment.new_name)
    start = getattr(member, "agp_component_start", 1)
    end = getattr(member, "agp_component_end", len(member.record.seq))
    orientation = getattr(member, "agp_component_orientation", "+")
    source = getattr(member, "agp_source", "ordered_contig")
    notes = getattr(member, "agp_notes", member.assignment.contig)

    trim = min(max(trim_left_bp or 0, 0), len(member.record.seq))
    if trim:
        if hasattr(member, "agp_component_orientation") and orientation not in {"+", "-"}:
            raise ValueError(
                "Cannot emit trimmed AGP coordinates for component "
                f"{component_id!r} with orientation {orientation!r}."
            )
        if orientation == "-":
            end -= trim
        else:
            start += trim

    return component_id, start, end, orientation, source, notes


def oriented_component_span(node_length, orientation, trim_left_bp):
    trim = min(max(trim_left_bp or 0, 0), node_length)
    if orientation == "-":
        return 1, node_length - trim
    return trim + 1, node_length


def graph_fill_parts(graph, plan):
    parts = []
    for edge in plan.graph_path[:-1]:
        overlap_bp = edge.overlap_bp or 0
        node = graph.nodes[edge.target]
        start, end = oriented_component_span(node.length, edge.target_orientation, overlap_bp)
        if end >= start:
            parts.append((edge.target, start, end, edge.target_orientation))
    return parts


def build_gapfill_agp_parts(groups, unassigned, plans, graph):
    grouped_plans = plans_by_scaffold(plans)
    parts = []
    for scaffold_name, members in groups.items():
        if not members:
            continue
        cursor = 1
        part_number = 1
        first = members[0]
        component_id, start, end, orientation, source, notes = component_fields_for_member(first)
        cursor, part_number = add_component_part(
            parts,
            scaffold_name,
            cursor,
            part_number,
            component_id,
            start,
            end,
            orientation,
            source=source,
            status="unchanged",
            notes=notes,
        )
        for plan, member in zip(grouped_plans.get(scaffold_name, []), members[1:]):
            trim_left_bp = 0
            if plan.applied:
                for node, start, end, orientation in graph_fill_parts(graph, plan):
                    cursor, part_number = add_component_part(
                        parts,
                        scaffold_name,
                        cursor,
                        part_number,
                        node,
                        start,
                        end,
                        orientation,
                        source="graph_fill",
                        status=plan.graph_status,
                        notes=plan.path_nodes,
                    )
                trim_left_bp = plan.right_trim_bp
            else:
                cursor, part_number = add_gap_part(
                    parts,
                    scaffold_name,
                    cursor,
                    part_number,
                    plan.fallback_gap_bp,
                    source="fallback_gap",
                    status=plan.fill_status,
                    notes=f"{plan.left_contig}|{plan.right_contig}",
                    gap_type=getattr(plan, "agp_gap_type", "scaffold"),
                    linkage=getattr(plan, "agp_linkage", "yes"),
                    linkage_evidence=getattr(
                        plan,
                        "agp_linkage_evidence",
                        "align_genus",
                    ),
                )
            component_id, start, end, orientation, source, notes = component_fields_for_member(
                member,
                trim_left_bp,
            )
            cursor, part_number = add_component_part(
                parts,
                scaffold_name,
                cursor,
                part_number,
                component_id,
                start,
                end,
                orientation,
                source=source,
                status="trimmed" if trim_left_bp else "unchanged",
                notes=notes,
            )

    for record in unassigned:
        add_component_part(
            parts,
            record.name,
            1,
            1,
            record.name,
            1,
            len(record.seq),
            "+",
            source="unassigned_record",
            status="unchanged",
        )
    return parts


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
        "gaf_support_status",
        "gaf_selected_reads",
        "hic_path_support",
        "hic_best_alt_support",
        "ref_path_support",
        "ref_best_alt_support",
        "longread_bridge_reads",
        "longread_orientation_summary",
        "longread_read_order_summary",
        "longread_median_read_gap_bp",
        "patch_candidate_count",
        "patch_best_id",
        "patch_best_source",
        "patch_best_bp",
        "patch_graph_status",
        "patch_best_notes",
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
        "patch_best_sequence",
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
        plan.gaf_support_status,
        plan.gaf_selected_reads,
        plan.hic_path_support if plan.hic_path_support is not None else ".",
        plan.hic_best_alt_support if plan.hic_best_alt_support is not None else ".",
        plan.ref_path_support if plan.ref_path_support is not None else ".",
        plan.ref_best_alt_support if plan.ref_best_alt_support is not None else ".",
        plan.longread_bridge_reads if plan.longread_bridge_reads is not None else ".",
        plan.longread_orientation_summary,
        plan.longread_read_order_summary,
        plan.longread_median_read_gap_bp,
        plan.patch_candidate_count,
        plan.patch_best_id,
        plan.patch_best_source,
        plan.patch_best_bp if plan.patch_best_bp is not None else ".",
        plan.patch_graph_status,
        plan.patch_best_notes,
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
        plan.patch_best_sequence if include_sequences and plan.patch_best_sequence else ".",
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
    ensure_parent_dir(path)
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
    ensure_parent_dir(path)
    with open(path, "w") as out:
        out.write(html_text)


def write_filled_fasta(path, records, simple_headers):
    ensure_parent_dir(path)
    with open(path, "w") as out:
        for record in records:
            out.write(f">{filled_header(record, simple_headers)}\n")
            write_wrapped(out, record.seq)


def write_run_summary(path, args, output_paths, plans, filled_records):
    fillable = sum(1 for plan in plans if plan.fill_status == "fillable")
    applied = sum(1 for plan in plans if plan.applied)
    ensure_parent_dir(path)
    with open(path, "w") as out:
        out.write("chromo gapfill\n")
        out.write("\nInputs\n")
        out.write(
            "input_mode\t"
            f"{'scaffold_agp' if args.scaffold_fasta else 'ordered_fasta'}\n"
        )
        out.write(f"ordered_fasta\t{args.ordered_fasta if args.ordered_fasta else '.'}\n")
        out.write(f"assignments\t{args.assignments if args.assignments else '.'}\n")
        out.write(f"scaffold_fasta\t{args.scaffold_fasta if args.scaffold_fasta else '.'}\n")
        out.write(f"agp\t{args.agp if args.agp else '.'}\n")
        out.write(f"gfa\t{args.gfa}\n")
        out.write(
            "graph_mode\t"
            f"{'projected_gfa_paths' if getattr(args, 'project_gfa_paths', False) else 'direct_gfa_nodes'}\n"
        )
        out.write(f"gaf\t{args.gaf if args.gaf else '.'}\n")
        out.write(f"hic_pairs\t{args.hic_pairs if args.hic_pairs else '.'}\n")
        out.write(f"ref_paf\t{args.ref_paf if args.ref_paf else '.'}\n")
        out.write(f"read_paf\t{args.read_paf if args.read_paf else '.'}\n")
        out.write(f"patch_table\t{args.patch_table if args.patch_table else '.'}\n")
        out.write(f"patch_fasta\t{args.patch_fasta if args.patch_fasta else '.'}\n")
        out.write(f"reviewed_plan\t{args.reviewed_plan if args.reviewed_plan else '.'}\n")
        out.write("\nParameters\n")
        out.write(f"apply\t{args.apply}\n")
        out.write(f"apply_all_fillable\t{getattr(args, 'apply_all_fillable', False)}\n")
        out.write(f"project_gfa_paths\t{getattr(args, 'project_gfa_paths', False)}\n")
        out.write(f"projection_trim_overlaps\t{getattr(args, 'projection_trim_overlaps', False)}\n")
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
        out.write(f"min_read_mapq\t{args.min_read_mapq}\n")
        out.write(f"read_terminal_window_bp\t{args.read_terminal_window_bp}\n")
        out.write(f"read_min_anchor_bp\t{args.read_min_anchor_bp}\n")
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


def validate_plan_args(args):
    ordered_mode = bool(args.ordered_fasta or args.assignments)
    scaffold_agp_mode = bool(args.scaffold_fasta or args.agp)
    apply_all_fillable = getattr(args, "apply_all_fillable", False)
    apply_output = getattr(args, "apply", False)
    reviewed_plan = getattr(args, "reviewed_plan", None)
    if ordered_mode and scaffold_agp_mode:
        raise ValueError(
            "Choose either --ordered-fasta/--assignments or --scaffold-fasta/--agp, "
            "not both."
        )
    if not ordered_mode and not scaffold_agp_mode:
        raise ValueError(
            "Provide either --ordered-fasta with --assignments, or "
            "--scaffold-fasta with --agp."
        )
    if ordered_mode and not (args.ordered_fasta and args.assignments):
        raise ValueError("--ordered-fasta and --assignments must be provided together")
    if scaffold_agp_mode and not (args.scaffold_fasta and args.agp):
        raise ValueError("--scaffold-fasta and --agp must be provided together")
    if getattr(args, "patch_fasta", None) and not getattr(args, "patch_table", None):
        raise ValueError("--patch-fasta requires --patch-table")
    if apply_all_fillable and not apply_output:
        raise ValueError("--apply-all-fillable requires --apply")
    if apply_output and not reviewed_plan and not apply_all_fillable:
        raise ValueError(
            "--apply requires --reviewed-plan or --apply-all-fillable. "
            "Use --reviewed-plan for production runs, or --apply-all-fillable "
            "when you intentionally want every currently fillable path applied."
        )
    if reviewed_plan and apply_all_fillable:
        raise ValueError("Use either --reviewed-plan or --apply-all-fillable, not both")
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
    if args.min_read_mapq < 0:
        raise ValueError("--min-read-mapq must be zero or greater")
    if args.read_terminal_window_bp < 1:
        raise ValueError("--read-terminal-window-bp must be at least 1")
    if args.read_min_anchor_bp < 1:
        raise ValueError("--read-min-anchor-bp must be at least 1")


def build_plan_context(args):
    validate_plan_args(args)
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
    read_evidence = (
        read_long_read_paf(args.read_paf, min_mapq=args.min_read_mapq)
        if args.read_paf
        else None
    )
    patch_candidates = (
        read_patch_candidates(args.patch_table, getattr(args, "patch_fasta", None))
        if getattr(args, "patch_table", None)
        else []
    )
    if args.scaffold_fasta:
        records = read_ordered_fasta(args.scaffold_fasta)
        groups, gap_parts_by_pair = build_scaffold_agp_groups(
            records,
            read_agp(args.agp),
            args.scaffold_fasta,
        )
        unassigned = []
        if getattr(args, "project_gfa_paths", False):
            projections = projections_by_contig(
                build_path_projection(
                    graph,
                    path_names=requested_projection_names(groups),
                    trim_overlaps=getattr(args, "projection_trim_overlaps", False),
                )
            )
            plans = build_projected_agp_fill_plans(
                groups,
                gap_parts_by_pair,
                graph,
                projections,
                args,
                gaf_records,
                hic_contacts,
                ref_placements,
            )
        else:
            plans = build_agp_fill_plans(
                groups,
                gap_parts_by_pair,
                graph,
                args,
                gaf_records,
                hic_contacts,
                ref_placements,
            )
    else:
        assignments = read_assignments(args.assignments)
        records = read_ordered_fasta(args.ordered_fasta)
        groups, unassigned = group_scaffold_members(records, assignments)
        if getattr(args, "project_gfa_paths", False):
            projections = projections_by_contig(
                build_path_projection(
                    graph,
                    path_names=requested_projection_names(groups),
                    trim_overlaps=getattr(args, "projection_trim_overlaps", False),
                )
            )
            plans = build_projected_fill_plans(
                groups,
                graph,
                projections,
                args,
                gaf_records,
                hic_contacts,
                ref_placements,
            )
        else:
            plans = build_fill_plans(
                groups,
                graph,
                args,
                gaf_records,
                hic_contacts,
                ref_placements,
            )
    annotate_longread_bridge_support(groups, plans, read_evidence, args)
    annotate_patch_candidates(plans, patch_candidates)
    return groups, unassigned, plans


def run(args):
    prefix = Path(args.output_prefix)

    output_paths = {
        "fill_plan": Path(str(prefix) + ".gapfill_plan.tsv"),
        "agp": Path(str(prefix) + ".gapfilled.agp"),
        "components": Path(str(prefix) + ".gapfilled_components.tsv"),
        "submission_checklist": Path(str(prefix) + ".submission_checklist.tsv"),
        "run_summary": Path(str(prefix) + ".run_summary.txt"),
    }
    if args.apply:
        output_paths["filled_fasta"] = Path(str(prefix) + ".gapfilled.fa")
    if args.review_html:
        output_paths["review_html"] = Path(args.review_html)
    ensure_output_dirs(output_paths)

    groups, unassigned, plans = build_plan_context(args)
    if args.reviewed_plan:
        apply_reviewed_plan(plans, read_reviewed_plan(args.reviewed_plan))

    filled_records = None
    if args.apply:
        filled_records = build_filled_scaffolds(groups, unassigned, plans, args)
        write_filled_fasta(output_paths["filled_fasta"], filled_records, args.simple_headers)

    graph = read_gfa(args.gfa)
    agp_parts = build_gapfill_agp_parts(groups, unassigned, plans, graph)
    write_agp(output_paths["agp"], agp_parts)
    write_component_tsv(output_paths["components"], agp_parts)
    write_fill_plan(output_paths["fill_plan"], plans, args.include_fill_sequences)
    if args.review_html:
        write_review_html(output_paths["review_html"], plans, args.include_fill_sequences)
    write_submission_checklist(
        output_paths["submission_checklist"],
        "chromo gapfill",
        output_paths,
        filled_records,
        agp_parts,
    )
    write_run_summary(output_paths["run_summary"], args, output_paths, plans, filled_records)

    sys.stderr.write(f"Planned {len(plans)} graph gap fill(s).\n")
    for status, count in sorted(status_counts(plans).items()):
        sys.stderr.write(f"  {status}: {count}\n")
    sys.stderr.write(f"Wrote gapfill plan: {output_paths['fill_plan']}\n")
    sys.stderr.write(f"Wrote gapfilled AGP: {output_paths['agp']}\n")
    sys.stderr.write(f"Wrote gapfilled components: {output_paths['components']}\n")
    sys.stderr.write(f"Wrote submission checklist: {output_paths['submission_checklist']}\n")
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
