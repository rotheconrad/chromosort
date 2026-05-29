"""Prepare task-specific review tables for ChromoSort commands."""

import argparse
import sys
from collections import defaultdict
from pathlib import Path
from typing import Optional, Sequence

from . import fix_contigs
from . import gapfill as gapfill_mod
from . import scaffold as scaffold_mod
from .gaf import read_gaf, summarize_gaf_node, summarize_gaf_traversal
from .graph import graph_node_evidence, read_gfa
from .longreads import read_long_read_paf, summarize_breakpoint, summarize_contig_bridge
from .paths import ensure_output_dirs
from .reference_order import alignment_source_from_args
from .review import ReviewEvent, event_table_path, write_review_events


FIX_REVIEW_COLUMNS = [
    "source_contig",
    "new_contig",
    "part_index",
    "dominant_ref",
    "slice_start",
    "slice_end",
    "piece_bp",
    "alignment_query_start",
    "alignment_query_end",
    "orientation",
    "reverse_complemented",
    "avg_identity",
    "segment_count",
    "planner_status",
    "planner_breakpoints",
    "planner_ref_transition",
    "planner_score",
    "graph_node",
    "graph_node_status",
    "graph_neighbor_count",
    "graph_self_loop",
    "gaf_node",
    "gaf_node_reads",
    "gaf_node_traversals",
    "gaf_node_orientations",
    "gaf_path_examples",
    "longread_breakpoint_position",
    "longread_spanning_reads",
    "longread_split_reads",
    "longread_left_edge_reads",
    "longread_right_edge_reads",
    "longread_nearby_reads",
]

SCAFFOLD_REVIEW_COLUMNS = [
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
    "overlap_action",
    "graph_status",
    "graph_direct_edge",
    "graph_path_nodes",
    "gaf_path_nodes",
    "gaf_path_support",
    "gaf_best_alt_path_nodes",
    "gaf_best_alt_support",
    "gaf_support_status",
    "gaf_selected_reads",
    "longread_bridge_reads",
    "longread_orientation_summary",
    "longread_read_order_summary",
    "longread_median_read_gap_bp",
]

GAPFILL_REVIEW_COLUMNS = [
    column
    for column in gapfill_mod.fill_plan_header()
    if column not in {"accept_fill", "applied"}
] + [
    "longread_bridge_reads",
    "longread_orientation_summary",
    "longread_read_order_summary",
    "longread_median_read_gap_bp",
]


def add_fix_planner_args(ap):
    ap.add_argument(
        "--mode",
        choices=fix_contigs.MODE_CHOICES,
        default="conservative",
        help="Breakpoint planner mode used to prepare candidate fix rows.",
    )
    ap.add_argument("--min-segment-bp", type=int, default=10_000)
    ap.add_argument("--min-segment-idy", type=float, default=0.0)
    ap.add_argument("--min-mapq", type=int, default=0)
    ap.add_argument("--include-secondary-paf", action="store_true")
    ap.add_argument("--max-merge-gap", type=int, default=1_000)
    ap.add_argument("--min-piece-bp", type=int, default=1)
    ap.add_argument("--breakpoint-penalty-bp", type=float, default=50_000.0)
    ap.add_argument("--min-piece-aligned-bp", type=int, default=50_000)
    ap.add_argument("--min-piece-query-frac", type=float, default=0.05)
    ap.add_argument("--complex-inversion-min-piece-aligned-bp", type=int, default=1_000_000)
    ap.add_argument("--complex-inversion-min-overlap-frac", type=float, default=0.50)
    ap.add_argument("--max-breakpoints-per-contig", type=int, default=4)
    ap.add_argument("--name-separator", default="-")
    ap.add_argument("--orient-to-reference", action="store_true")
    ap.set_defaults(simple_headers=False, pieces_only=False)


def parse_fix_args(argv=None, prog=None):
    ap = argparse.ArgumentParser(
        prog=prog,
        description="Prepare an editable TSV of candidate chromo fix decisions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument(
        "-f",
        "--assembly-fasta",
        required=True,
        help="Assembly FASTA containing contigs to evaluate for fixing.",
    )
    alignment_group = ap.add_mutually_exclusive_group(required=True)
    alignment_group.add_argument("-c", "--coords", help="MUMmer show-coords file.")
    alignment_group.add_argument("--paf", help="minimap2 PAF file.")
    ap.add_argument("--contigs", nargs="+", default=[], help="Contigs to evaluate.")
    ap.add_argument("--contigs-file", default=None, help="One contig name per line.")
    ap.add_argument("--all", action="store_true", dest="all_contigs", help="Evaluate all split-signal contigs.")
    ap.add_argument("--gfa", default=None, help="Optional assembly graph GFA for node context.")
    ap.add_argument("--gaf", default=None, help="Optional long-read-to-graph GAF for advisory graph traversal context.")
    ap.add_argument("--min-gaf-mapq", type=int, default=20, help="Minimum GAF MAPQ for advisory graph traversal context.")
    ap.add_argument("--read-paf", default=None, help="Optional long-read-to-assembly PAF.")
    ap.add_argument("--min-read-mapq", type=int, default=0, help="Minimum MAPQ for --read-paf rows.")
    ap.add_argument("--read-window-bp", type=int, default=5_000, help="Breakpoint window for read evidence.")
    ap.add_argument("--read-min-anchor-bp", type=int, default=1_000, help="Minimum read anchor on each side of a breakpoint.")
    ap.add_argument(
        "-o",
        "--output-prefix",
        required=True,
        help="Output prefix. Writes <prefix>.fix_review.tsv.",
    )
    add_fix_planner_args(ap)
    return ap.parse_args(argv)


def parse_scaffold_args(argv=None, prog=None):
    ap = argparse.ArgumentParser(
        prog=prog,
        description="Prepare an editable TSV of candidate chromo scaffold junction decisions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument(
        "-f",
        "--ordered-fasta",
        required=True,
        help="Final ordered FASTA from chromo sort.",
    )
    ap.add_argument(
        "-a",
        "--assignments",
        required=True,
        help="Matching <prefix>.contig_assignments.tsv from chromo sort.",
    )
    ap.add_argument(
        "-o",
        "--output-prefix",
        required=True,
        help="Output prefix. Writes <prefix>.scaffold_review.tsv.",
    )
    ap.add_argument("--gfa", default=None, help="Optional assembly graph GFA for junction context.")
    ap.add_argument("--gaf", default=None, help="Optional long-read-to-graph GAF for junction path support.")
    ap.add_argument("--min-gaf-mapq", type=int, default=20, help="Minimum GAF MAPQ for read-path support.")
    ap.add_argument("--min-gaf-path-support", type=int, default=1)
    ap.add_argument("--gaf-max-candidate-paths", type=int, default=2)
    ap.add_argument("--read-paf", default=None, help="Optional long-read-to-assembly PAF.")
    ap.add_argument("--min-read-mapq", type=int, default=0, help="Minimum MAPQ for --read-paf rows.")
    ap.add_argument("--read-terminal-window-bp", type=int, default=5_000)
    ap.add_argument("--read-min-anchor-bp", type=int, default=500)
    ap.add_argument("--fixed-gap-bp", type=int, default=None)
    ap.add_argument(
        "--overlap-policy",
        choices=["zero-gap", "warn", "trim-reference", "trim-sequence"],
        default="zero-gap",
    )
    ap.add_argument("--trim-sequence-min-identity", type=float, default=0.98)
    ap.add_argument(
        "--graph-overlap-policy",
        choices=["report", "warn", "confirm"],
        default="report",
    )
    ap.add_argument("--graph-max-path-edges", type=int, default=4)
    ap.set_defaults(simple_headers=False)
    return ap.parse_args(argv)


def parse_gapfill_args(argv=None, prog=None):
    ap = argparse.ArgumentParser(
        prog=prog,
        description="Prepare an editable TSV of candidate chromo gapfill decisions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("-f", "--ordered-fasta", required=True, help="Final ordered FASTA from chromo sort.")
    ap.add_argument(
        "-a",
        "--assignments",
        required=True,
        help="Matching <prefix>.contig_assignments.tsv from chromo sort.",
    )
    ap.add_argument("--gfa", required=True, help="Assembly graph GFA containing segment sequences and links.")
    ap.add_argument("--gaf", default=None, help="Optional GAF graph alignments for read-path support.")
    ap.add_argument("--hic-pairs", default=None, help="Optional graph-node contact count TSV.")
    ap.add_argument("--ref-paf", default=None, help="Optional reference-to-assembly PAF for graph-node placement.")
    ap.add_argument("--read-paf", default=None, help="Optional long-read-to-assembly PAF for contig-end bridge support.")
    ap.add_argument("--min-read-mapq", type=int, default=0, help="Minimum MAPQ for --read-paf rows.")
    ap.add_argument("--read-terminal-window-bp", type=int, default=5_000)
    ap.add_argument("--read-min-anchor-bp", type=int, default=500)
    ap.add_argument(
        "-o",
        "--output-prefix",
        required=True,
        help="Output prefix. Writes <prefix>.gapfill_review.tsv.",
    )
    ap.add_argument("--fixed-gap-bp", type=int, default=None)
    ap.add_argument("--max-path-edges", type=int, default=4)
    ap.add_argument("--max-candidate-paths", type=int, default=2)
    ap.add_argument("--min-gaf-mapq", type=int, default=20)
    ap.add_argument("--min-gaf-path-support", type=int, default=1)
    ap.add_argument("--min-hic-path-support", type=int, default=1)
    ap.add_argument("--min-ref-path-support", type=int, default=1)
    ap.add_argument("--min-ref-paf-mapq", type=int, default=0)
    ap.add_argument("--min-ref-paf-idy", type=float, default=0.0)
    ap.add_argument("--include-secondary-ref-paf", action="store_true")
    ap.add_argument("--max-fill-bp", type=int, default=1_000_000)
    ap.add_argument(
        "--include-fill-sequences",
        action="store_true",
        help="Include candidate fill sequences in the review table.",
    )
    return ap.parse_args(argv)


def fmt(value, digits=3):
    if value is None:
        return "."
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def bool_text(value):
    if value is None:
        return "."
    return "yes" if value else "no"


def graph_fields(graph, contig):
    if graph is None:
        return {
            "graph_node": ".",
            "graph_node_status": ".",
            "graph_neighbor_count": ".",
            "graph_self_loop": ".",
        }
    evidence = graph_node_evidence(graph, [contig])
    return {
        "graph_node": evidence.graph_node,
        "graph_node_status": evidence.graph_node_status,
        "graph_neighbor_count": fmt(evidence.graph_neighbor_count, 0),
        "graph_self_loop": bool_text(evidence.graph_self_loop),
    }


def empty_gaf_node_fields():
    return {
        "gaf_node": ".",
        "gaf_node_reads": ".",
        "gaf_node_traversals": ".",
        "gaf_node_orientations": ".",
        "gaf_path_examples": ".",
    }


def gaf_node_fields(gaf_records, graph, contig):
    if not gaf_records:
        return empty_gaf_node_fields()

    candidates = []
    if graph is not None:
        evidence = graph_node_evidence(graph, [contig])
        if evidence.graph_node and evidence.graph_node != ".":
            candidates.append(evidence.graph_node)
    candidates.append(contig)

    seen = []
    for candidate in candidates:
        if candidate and candidate != "." and candidate not in seen:
            seen.append(candidate)
    summaries = [summarize_gaf_node(gaf_records, candidate) for candidate in seen]
    summary = max(summaries, key=lambda item: item.traversal_count) if summaries else None
    if summary is None:
        return empty_gaf_node_fields()

    return {
        "gaf_node": summary.node,
        "gaf_node_reads": str(summary.read_count),
        "gaf_node_traversals": str(summary.traversal_count),
        "gaf_node_orientations": summary.orientation_summary,
        "gaf_path_examples": ";".join(summary.path_examples) if summary.path_examples else ".",
    }


def read_support_fields(read_evidence, contig, position, args):
    if read_evidence is None or position is None:
        return {
            "longread_breakpoint_position": ".",
            "longread_spanning_reads": ".",
            "longread_split_reads": ".",
            "longread_left_edge_reads": ".",
            "longread_right_edge_reads": ".",
            "longread_nearby_reads": ".",
        }
    support = summarize_breakpoint(
        read_evidence,
        contig,
        position,
        window_bp=args.read_window_bp,
        min_anchor_bp=args.read_min_anchor_bp,
    )
    return {
        "longread_breakpoint_position": str(position),
        "longread_spanning_reads": str(support.spanning_count),
        "longread_split_reads": str(support.split_count),
        "longread_left_edge_reads": str(support.left_edge_count),
        "longread_right_edge_reads": str(support.right_edge_count),
        "longread_nearby_reads": str(support.nearby_count),
    }


def split_piece_event(plan, piece, index, graph, gaf_records, read_evidence, args):
    is_last = index == len(plan.pieces) - 1
    breakpoint_position = None if is_last else piece.slice_end
    fields = {
        "source_contig": piece.original_contig,
        "new_contig": piece.new_name,
        "part_index": piece.part_index,
        "dominant_ref": piece.ref,
        "slice_start": piece.slice_start + 1,
        "slice_end": piece.slice_end,
        "piece_bp": piece.length,
        "alignment_query_start": piece.align_start + 1,
        "alignment_query_end": piece.align_end,
        "orientation": piece.orientation,
        "reverse_complemented": bool_text(piece.reverse_complemented),
        "avg_identity": fmt(piece.avg_identity),
        "segment_count": piece.segment_count,
        "planner_status": plan.status,
        "planner_breakpoints": plan.planner_breakpoints,
        "planner_ref_transition": bool_text(plan.planner_ref_transition),
        "planner_score": fmt(plan.planner_score, 1),
    }
    fields.update(graph_fields(graph, piece.original_contig))
    fields.update(gaf_node_fields(gaf_records, graph, piece.original_contig))
    fields.update(read_support_fields(read_evidence, piece.original_contig, breakpoint_position, args))
    return ReviewEvent(
        event_id=f"fix:{piece.original_contig}:{piece.part_index}",
        task="fix",
        action="split_piece",
        target=piece.original_contig,
        accept=True,
        status="candidate",
        confidence=".",
        reason=plan.reason,
        fields=fields,
    )


def no_split_event(contig, plan, graph, gaf_records):
    fields = {
        "source_contig": contig,
        "new_contig": ".",
        "part_index": ".",
        "dominant_ref": ".",
        "slice_start": ".",
        "slice_end": ".",
        "piece_bp": ".",
        "alignment_query_start": ".",
        "alignment_query_end": ".",
        "orientation": ".",
        "reverse_complemented": ".",
        "avg_identity": ".",
        "segment_count": ".",
        "planner_status": plan.status,
        "planner_breakpoints": plan.planner_breakpoints,
        "planner_ref_transition": bool_text(plan.planner_ref_transition),
        "planner_score": fmt(plan.planner_score, 1),
        "longread_breakpoint_position": ".",
        "longread_spanning_reads": ".",
        "longread_split_reads": ".",
        "longread_left_edge_reads": ".",
        "longread_right_edge_reads": ".",
        "longread_nearby_reads": ".",
    }
    fields.update(graph_fields(graph, contig))
    fields.update(gaf_node_fields(gaf_records, graph, contig))
    return ReviewEvent(
        event_id=f"fix:{contig}:no_split",
        task="fix",
        action="no_split",
        target=contig,
        accept=False,
        status=plan.status,
        confidence=".",
        reason=plan.reason,
        fields=fields,
    )


def build_fix_events(args):
    alignment_path, alignment_format = alignment_source_from_args(args)
    explicit_requested = fix_contigs.read_requested_contigs(args.contigs, args.contigs_file)
    if args.all_contigs and explicit_requested:
        raise ValueError("use either --all or --contigs/--contigs-file, not both")
    if not explicit_requested and not args.all_contigs:
        raise ValueError("provide at least one contig via --contigs/--contigs-file or use --all")
    if args.min_gaf_mapq < 0:
        raise ValueError("--min-gaf-mapq must be zero or greater")

    collect_for = None if args.all_contigs else explicit_requested
    blocks_by_contig = fix_contigs.collect_blocks(
        alignment_path,
        alignment_format,
        collect_for,
        args.min_segment_bp,
        args.min_segment_idy,
        args.max_merge_gap,
        min_mapq=args.min_mapq,
        include_secondary_paf=args.include_secondary_paf,
    )
    requested = (
        fix_contigs.all_requested_contigs(args.assembly_fasta, blocks_by_contig, args)
        if args.all_contigs
        else explicit_requested
    )
    plans = fix_contigs.build_plans(args.assembly_fasta, requested, blocks_by_contig, args)
    fix_contigs.apply_breakpoint_guard(plans, args)

    graph = read_gfa(args.gfa) if args.gfa else None
    gaf_records = read_gaf(args.gaf, args.min_gaf_mapq) if args.gaf else []
    read_evidence = (
        read_long_read_paf(args.read_paf, min_mapq=args.min_read_mapq)
        if args.read_paf
        else None
    )

    events = []
    for contig in requested:
        plan = plans[contig]
        if plan.status == "split":
            for index, piece in enumerate(plan.pieces):
                events.append(split_piece_event(plan, piece, index, graph, gaf_records, read_evidence, args))
        else:
            events.append(no_split_event(contig, plan, graph, gaf_records))
    return events


def graph_gap_fields(graph_records, gap):
    record = graph_records.get((gap.scaffold, gap.left_contig, gap.right_contig))
    if record is None:
        return {
            "graph_status": ".",
            "graph_direct_edge": ".",
            "graph_path_nodes": ".",
        }
    return {
        "graph_status": record.graph_status,
        "graph_direct_edge": bool_text(record.direct_edge),
        "graph_path_nodes": record.shortest_path_nodes,
    }


def empty_gaf_fields():
    return {
        "gaf_path_nodes": ".",
        "gaf_path_support": ".",
        "gaf_best_alt_path_nodes": ".",
        "gaf_best_alt_support": ".",
        "gaf_support_status": ".",
        "gaf_selected_reads": ".",
    }


def gaf_scaffold_fields(graph, gaf_records, left_member, right_member, args):
    if not gaf_records:
        return empty_gaf_fields()
    if graph is None:
        fields = empty_gaf_fields()
        fields["gaf_support_status"] = "missing_gfa"
        return fields

    left_node = scaffold_mod.graph_node_name(left_member, graph)
    right_node = scaffold_mod.graph_node_name(right_member, graph)
    missing_status = scaffold_mod.missing_graph_status(left_node, right_node)
    if missing_status is not None:
        fields = empty_gaf_fields()
        fields["gaf_support_status"] = missing_status
        return fields

    paths = gapfill_mod.enumerate_graph_paths(
        graph,
        left_node,
        right_node,
        scaffold_mod.graph_orientation(left_member),
        scaffold_mod.graph_orientation(right_member),
        args.graph_max_path_edges,
        args.gaf_max_candidate_paths,
    ).paths
    if not paths:
        fields = empty_gaf_fields()
        fields["gaf_support_status"] = "no_graph_path"
        return fields

    summary = summarize_gaf_traversal(
        gaf_records,
        paths,
        selected_index=0,
        min_support=args.min_gaf_path_support,
    )
    return {
        "gaf_path_nodes": summary.selected_path_nodes,
        "gaf_path_support": str(summary.selected_support),
        "gaf_best_alt_path_nodes": summary.best_alt_path_nodes,
        "gaf_best_alt_support": str(summary.best_alt_support),
        "gaf_support_status": summary.support_status,
        "gaf_selected_reads": ",".join(summary.selected_reads) if summary.selected_reads else ".",
    }


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


def read_bridge_fields(read_evidence, left_member, right_member, args):
    if read_evidence is None:
        return {
            "longread_bridge_reads": ".",
            "longread_orientation_summary": ".",
            "longread_read_order_summary": ".",
            "longread_median_read_gap_bp": ".",
        }
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
    support = max(supports, key=lambda item: item.bridge_count) if supports else None
    if support is None:
        return {
            "longread_bridge_reads": "0",
            "longread_orientation_summary": ".",
            "longread_read_order_summary": ".",
            "longread_median_read_gap_bp": ".",
        }
    median_gap = support.median_read_gap_bp
    return {
        "longread_bridge_reads": str(support.bridge_count),
        "longread_orientation_summary": support.orientation_summary,
        "longread_read_order_summary": support.read_order_summary,
        "longread_median_read_gap_bp": "." if median_gap is None else str(median_gap),
    }


def build_scaffold_events(args):
    if args.fixed_gap_bp is not None and args.fixed_gap_bp < 0:
        raise ValueError("--fixed-gap-bp must be zero or greater")
    if args.graph_max_path_edges < 1:
        raise ValueError("--graph-max-path-edges must be at least 1")
    if args.gaf_max_candidate_paths < 1:
        raise ValueError("--gaf-max-candidate-paths must be at least 1")
    if args.min_gaf_mapq < 0:
        raise ValueError("--min-gaf-mapq must be zero or greater")
    if args.min_gaf_path_support < 1:
        raise ValueError("--min-gaf-path-support must be at least 1")
    if args.graph_overlap_policy != "report" and not args.gfa:
        raise ValueError("--graph-overlap-policy warn/confirm requires --gfa")
    if args.gaf and not args.gfa:
        raise ValueError("--gaf requires --gfa")
    assignments = scaffold_mod.read_assignments(args.assignments)
    records = scaffold_mod.read_ordered_fasta(args.ordered_fasta)
    groups, unassigned = scaffold_mod.group_scaffold_members(records, assignments)
    del unassigned
    graph = read_gfa(args.gfa) if args.gfa else None
    scaffolds = scaffold_mod.build_scaffolds(groups, args.fixed_gap_bp, args, graph)
    graph_records = {}
    if graph is not None:
        graph_records = {
            (record.scaffold, record.left_contig, record.right_contig): record
            for record in scaffold_mod.build_graph_gap_records(
                scaffolds,
                graph,
                args.graph_max_path_edges,
            )
        }
    read_evidence = (
        read_long_read_paf(args.read_paf, min_mapq=args.min_read_mapq)
        if args.read_paf
        else None
    )
    gaf_records = read_gaf(args.gaf, args.min_gaf_mapq) if args.gaf else []

    events = []
    for scaffold in scaffolds:
        for gap, left, right in zip(scaffold.gaps, scaffold.members, scaffold.members[1:]):
            fields = {
                "scaffold": gap.scaffold,
                "left_contig": gap.left_contig,
                "right_contig": gap.right_contig,
                "left_ref_end": gap.left_ref_end,
                "right_ref_start": gap.right_ref_start,
                "raw_inferred_gap_bp": gap.raw_inferred_gap_bp,
                "gap_bp": gap.gap_bp,
                "gap_mode": gap.gap_mode,
                "overlap_bp": gap.overlap_bp,
                "overlap_class": gap.overlap_class,
                "overlap_action": gap.overlap_action,
            }
            fields.update(graph_gap_fields(graph_records, gap))
            fields.update(gaf_scaffold_fields(graph, gaf_records, left, right, args))
            fields.update(read_bridge_fields(read_evidence, left, right, args))
            events.append(
                ReviewEvent(
                    event_id=f"scaffold:{gap.scaffold}:{gap.left_contig}:{gap.right_contig}",
                    task="scaffold",
                    action="scaffold_gap",
                    target=f"{gap.scaffold}:{gap.left_contig}|{gap.right_contig}",
                    accept=True,
                    status="candidate",
                    confidence=".",
                    reason="scaffold_junction",
                    fields=fields,
                )
            )
    return events


def adjacent_member_pairs(groups):
    pairs = {}
    for scaffold_name, members in groups.items():
        for left, right in zip(members, members[1:]):
            key = (scaffold_name, left.assignment.new_name, right.assignment.new_name)
            pairs[key] = (left, right)
    return pairs


def gapfill_read_bridge_fields(read_evidence, member_pairs, plan, args):
    if read_evidence is None:
        return {
            "longread_bridge_reads": ".",
            "longread_orientation_summary": ".",
            "longread_read_order_summary": ".",
            "longread_median_read_gap_bp": ".",
        }
    pair = member_pairs.get(gapfill_mod.fill_plan_key(plan))
    if pair is None:
        return {
            "longread_bridge_reads": "0",
            "longread_orientation_summary": ".",
            "longread_read_order_summary": ".",
            "longread_median_read_gap_bp": ".",
        }
    left, right = pair
    return read_bridge_fields(read_evidence, left, right, args)


def build_gapfill_events(args):
    groups, _unassigned, plans = gapfill_mod.build_plan_context(args)
    member_pairs = adjacent_member_pairs(groups)
    read_evidence = (
        read_long_read_paf(args.read_paf, min_mapq=args.min_read_mapq)
        if args.read_paf
        else None
    )

    events = []
    for plan in plans:
        fields = gapfill_mod.fill_plan_dict(plan, args.include_fill_sequences)
        fields.pop("accept_fill", None)
        fields.pop("applied", None)
        fields.update(gapfill_read_bridge_fields(read_evidence, member_pairs, plan, args))
        events.append(
            ReviewEvent(
                event_id=f"gapfill:{plan.scaffold}:{plan.left_contig}:{plan.right_contig}",
                task="gapfill",
                action="fill_path",
                target=f"{plan.scaffold}:{plan.left_contig}|{plan.right_contig}",
                accept=plan.fill_status == "fillable",
                status=plan.fill_status,
                confidence=".",
                reason=plan.reason,
                fields=fields,
            )
        )
    return events


def run_fix(args):
    prefix = Path(args.output_prefix)
    output_paths = {
        "fix_review": event_table_path(prefix, ".fix_review.tsv"),
    }
    ensure_output_dirs(output_paths)
    events = build_fix_events(args)
    write_review_events(output_paths["fix_review"], events, extra_columns=FIX_REVIEW_COLUMNS)
    status_counts = defaultdict(int)
    action_counts = defaultdict(int)
    for event in events:
        status_counts[event.status] += 1
        action_counts[event.action] += 1
    sys.stderr.write(f"Wrote fix review table: {output_paths['fix_review']}\n")
    for action, count in sorted(action_counts.items()):
        sys.stderr.write(f"  {action}: {count}\n")
    for status, count in sorted(status_counts.items()):
        sys.stderr.write(f"  {status}: {count}\n")


def run_scaffold(args):
    prefix = Path(args.output_prefix)
    output_paths = {
        "scaffold_review": event_table_path(prefix, ".scaffold_review.tsv"),
    }
    ensure_output_dirs(output_paths)
    events = build_scaffold_events(args)
    write_review_events(
        output_paths["scaffold_review"],
        events,
        extra_columns=SCAFFOLD_REVIEW_COLUMNS,
    )
    sys.stderr.write(f"Wrote scaffold review table: {output_paths['scaffold_review']}\n")
    sys.stderr.write(f"  scaffold_gap: {len(events)}\n")


def run_gapfill(args):
    prefix = Path(args.output_prefix)
    output_paths = {
        "gapfill_review": event_table_path(prefix, ".gapfill_review.tsv"),
    }
    ensure_output_dirs(output_paths)
    events = build_gapfill_events(args)
    write_review_events(
        output_paths["gapfill_review"],
        events,
        extra_columns=GAPFILL_REVIEW_COLUMNS,
    )
    status_counts = defaultdict(int)
    for event in events:
        status_counts[event.status] += 1
    sys.stderr.write(f"Wrote gapfill review table: {output_paths['gapfill_review']}\n")
    for status, count in sorted(status_counts.items()):
        sys.stderr.write(f"  {status}: {count}\n")


def main(argv: Optional[Sequence[str]] = None, prog: Optional[str] = None):
    if argv is None:
        argv = sys.argv[1:]
    parser = argparse.ArgumentParser(
        prog=prog,
        description="Prepare TSV review tables for ChromoSort command decisions.",
    )
    parser.add_argument(
        "mode",
        nargs="?",
        choices=["fix", "scaffold", "gapfill"],
        help="Review table mode to run.",
    )
    if not argv or argv[0] in {"-h", "--help"}:
        parser.print_help()
        return

    mode = argv[0]
    remaining = argv[1:]
    try:
        if mode == "fix":
            run_fix(parse_fix_args(remaining, prog=f"{prog} fix" if prog else None))
        elif mode == "scaffold":
            run_scaffold(parse_scaffold_args(remaining, prog=f"{prog} scaffold" if prog else None))
        elif mode == "gapfill":
            run_gapfill(parse_gapfill_args(remaining, prog=f"{prog} gapfill" if prog else None))
        else:
            parser.error(f"unknown eval mode: {mode}")
    except ValueError as exc:
        sys.stderr.write(f"ERROR: {exc}\n")
        sys.exit(2)


if __name__ == "__main__":
    main()
