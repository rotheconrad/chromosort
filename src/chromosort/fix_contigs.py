#!/usr/bin/env python3
"""
Split chimeric contigs using coords or PAF alignments.

This module is intentionally conservative by default. --contigs and
--contigs-file choose a user-reviewed subset to inspect, while --all scans every
contig with a split signal. The selected contigs then go through the same
mode-specific planner so small discordant blocks, INDEL-sized gaps, and local
SV-like noise can be smoothed over before breakpoints are accepted.
"""

import argparse
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

from .graph import graph_node_evidence, read_gfa
from .reference_order import (
    alignment_source_from_args,
    iter_alignments,
    iter_fasta_records,
    reverse_complement,
    write_wrapped,
)


@dataclass
class QueryBlock:
    contig: str
    ref: str
    ref_start: int
    ref_end: int
    query_start: int
    query_end: int
    orientation: str
    aligned_bp: int
    identity_bp: float
    segment_count: int

    @property
    def avg_identity(self):
        return self.identity_bp / self.aligned_bp if self.aligned_bp else 0.0

    @property
    def weighted_bp(self):
        return self.identity_bp / 100.0

    @property
    def query_span(self):
        return self.query_end - self.query_start


@dataclass
class SplitPiece:
    original_contig: str
    new_name: str
    part_index: int
    ref: str
    ref_start: int
    ref_end: int
    slice_start: int
    slice_end: int
    align_start: int
    align_end: int
    orientation: str
    avg_identity: float
    segment_count: int
    reverse_complemented: bool = False

    @property
    def length(self):
        return self.slice_end - self.slice_start


@dataclass
class ContigPlan:
    contig: str
    status: str
    pieces: List[SplitPiece]
    reason: str
    planner_score: float = 0.0
    planner_breakpoints: int = 0
    planner_ref_transition: bool = False


@dataclass
class BlockGroup:
    blocks: List[QueryBlock]
    summary: QueryBlock
    discordant_bp: float
    dominant_weighted_bp: float
    total_weighted_bp: float


MODE_CHOICES = ("conservative", "chromosome", "comprehensive", "sensitive")


def parse_args(argv: Optional[Sequence[str]] = None, prog: Optional[str] = None):
    ap = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Split selected chimeric contigs into reference-labeled pieces using "
            "MUMmer coords or minimap2 PAF alignments."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument(
        "-f",
        "--assembly-fasta",
        required=True,
        help="Assembly FASTA containing the contigs to fix.",
    )
    alignment_group = ap.add_mutually_exclusive_group(required=True)
    alignment_group.add_argument(
        "-c",
        "--coords",
        help="MUMmer show-coords file for reference-vs-assembly alignment.",
    )
    alignment_group.add_argument(
        "--paf",
        help="minimap2 PAF file for reference-vs-assembly alignment.",
    )
    ap.add_argument(
        "--contigs",
        nargs="+",
        default=[],
        help="Names of contigs to inspect and split when chimeric.",
    )
    ap.add_argument(
        "--contigs-file",
        default=None,
        help="Optional text file with one contig name per line.",
    )
    ap.add_argument(
        "--all",
        action="store_true",
        dest="all_contigs",
        help=(
            "Inspect all contigs with a passing split signal. Use --contigs or "
            "--contigs-file instead to inspect only a reviewed subset."
        ),
    )
    ap.add_argument(
        "--mode",
        choices=MODE_CHOICES,
        default="conservative",
        help=(
            "Breakpoint planner to use for selected contigs. conservative "
            "smooths weak discordance and cuts reference transitions plus "
            "complex same-reference orientation events; chromosome cuts only "
            "reference transitions; comprehensive also considers all "
            "same-reference orientation changes; sensitive disables smoothing "
            "and cuts every passing reference/orientation transition."
        ),
    )
    ap.add_argument(
        "--auto",
        action="store_true",
        dest="legacy_auto",
        help=argparse.SUPPRESS,
    )
    ap.add_argument(
        "--auto-sensitive",
        action="store_true",
        dest="legacy_auto_sensitive",
        help=argparse.SUPPRESS,
    )
    ap.add_argument(
        "--auto-split-inversions",
        action="store_true",
        dest="legacy_auto_split_inversions",
        help=argparse.SUPPRESS,
    )
    ap.add_argument(
        "--auto-complex-inversions",
        action=argparse.BooleanOptionalAction,
        default=None,
        dest="legacy_auto_complex_inversions",
        help=argparse.SUPPRESS,
    )
    ap.add_argument(
        "-o",
        "--output-fasta",
        required=True,
        help="Output FASTA with selected chimeras replaced by split pieces.",
    )
    ap.add_argument(
        "--report",
        required=True,
        help="TSV report describing split pieces and unsplit requested contigs.",
    )
    ap.add_argument(
        "--gfa",
        default=None,
        help=(
            "Optional assembly graph GFA. When provided, writes a report-only "
            "graph context table for requested contigs."
        ),
    )
    ap.add_argument(
        "--graph-report",
        default=None,
        help=(
            "Optional path for --gfa graph context output. Defaults to the "
            "--report path with a .graph.tsv suffix."
        ),
    )
    ap.add_argument(
        "--graph-guard",
        action="store_true",
        help=(
            "Requires --gfa. Emit conservative warnings when planned splits "
            "intersect simple graph neighborhoods or when unsplit reviewed "
            "contigs sit in complex graph neighborhoods. This does not change "
            "breakpoints or sequence output."
        ),
    )
    ap.add_argument(
        "--min-segment-bp",
        type=int,
        default=10_000,
        help="Minimum query-aligned bp for an alignment segment to inform splitting.",
    )
    ap.add_argument(
        "--min-segment-idy",
        type=float,
        default=0.0,
        help="Minimum percent identity for an alignment segment to inform splitting.",
    )
    ap.add_argument(
        "--min-mapq",
        type=int,
        default=0,
        help="Ignore PAF rows below this MAPQ. Ignored for MUMmer coords.",
    )
    ap.add_argument(
        "--include-secondary-paf",
        action="store_true",
        help="Include minimap2 PAF rows marked with tp:A:S. By default they are skipped.",
    )
    ap.add_argument(
        "--max-merge-gap",
        type=int,
        default=1_000,
        help=(
            "Merge neighboring alignment rows for the same reference and "
            "orientation when separated by no more than this many query bp."
        ),
    )
    ap.add_argument(
        "--min-piece-bp",
        type=int,
        default=1,
        help="Minimum emitted split-piece length.",
    )
    ap.add_argument(
        "--breakpoint-penalty-bp",
        type=float,
        default=50_000.0,
        dest="breakpoint_penalty_bp",
        help=(
            "Breakpoint penalty for smoothed modes, measured as identity-weighted "
            "aligned bp. A split is kept only when doing so explains more "
            "discordant support than this penalty."
        ),
    )
    ap.add_argument(
        "--auto-breakpoint-penalty-bp",
        type=float,
        dest="breakpoint_penalty_bp",
        help=argparse.SUPPRESS,
    )
    ap.add_argument(
        "--min-piece-aligned-bp",
        type=int,
        default=50_000,
        dest="min_piece_aligned_bp",
        help=(
            "Minimum dominant aligned bp required in each smoothed split piece. "
            "This prevents cutting off weak local SV or repeat hits."
        ),
    )
    ap.add_argument(
        "--auto-min-piece-aligned-bp",
        type=int,
        dest="min_piece_aligned_bp",
        help=argparse.SUPPRESS,
    )
    ap.add_argument(
        "--min-piece-query-frac",
        type=float,
        default=0.05,
        dest="min_piece_query_frac",
        help=(
            "Minimum query-span fraction required for each smoothed split piece. "
            "Set to 0 to disable this fraction check."
        ),
    )
    ap.add_argument(
        "--auto-min-piece-query-frac",
        type=float,
        dest="min_piece_query_frac",
        help=argparse.SUPPRESS,
    )
    ap.add_argument(
        "--complex-inversion-min-piece-aligned-bp",
        type=int,
        default=1_000_000,
        dest="complex_inversion_min_piece_aligned_bp",
        help=(
            "Minimum dominant aligned bp for pieces used to classify a "
            "same-reference orientation event as complex rather than a simple "
            "contiguous inversion."
        ),
    )
    ap.add_argument(
        "--auto-complex-inversion-min-piece-aligned-bp",
        type=int,
        dest="complex_inversion_min_piece_aligned_bp",
        help=argparse.SUPPRESS,
    )
    ap.add_argument(
        "--complex-inversion-min-overlap-frac",
        type=float,
        default=0.50,
        dest="complex_inversion_min_overlap_frac",
        help=(
            "Minimum reference-span overlap fraction, relative to the smaller "
            "piece, for a same-reference orientation event to be considered "
            "complex."
        ),
    )
    ap.add_argument(
        "--auto-complex-inversion-min-overlap-frac",
        type=float,
        dest="complex_inversion_min_overlap_frac",
        help=argparse.SUPPRESS,
    )
    ap.add_argument(
        "--max-breakpoints-per-contig",
        type=int,
        default=4,
        dest="max_breakpoints_per_contig",
        help=(
            "Maximum accepted breakpoints per contig. Set to a negative value "
            "to disable this guardrail."
        ),
    )
    ap.add_argument(
        "--auto-max-breakpoints",
        type=int,
        dest="max_breakpoints_per_contig",
        help=argparse.SUPPRESS,
    )
    ap.add_argument(
        "--name-separator",
        default="-",
        help="Separator used in split FASTA IDs: REF<sep>CONTIG<sep>PART.",
    )
    ap.add_argument(
        "--orient-to-reference",
        action="store_true",
        help="Reverse-complement split pieces whose alignment block is reverse-strand.",
    )
    ap.add_argument(
        "--simple-headers",
        action="store_true",
        help="Write FASTA headers containing only the new sequence ID.",
    )
    ap.add_argument(
        "--pieces-only",
        action="store_true",
        help="Write only split pieces instead of a full fixed assembly FASTA.",
    )
    args = ap.parse_args(argv)
    if args.legacy_auto:
        args.all_contigs = True
    if args.legacy_auto_sensitive:
        args.mode = "sensitive"
    elif args.legacy_auto_split_inversions:
        args.mode = "comprehensive"
    elif args.legacy_auto_complex_inversions is False:
        args.mode = "chromosome"
    return args


def read_requested_contigs(names: Iterable[str], path: Optional[str]) -> List[str]:
    requested = []
    seen = set()

    def add(name):
        clean = name.strip()
        if not clean or clean.startswith("#") or clean in seen:
            return
        requested.append(clean)
        seen.add(clean)

    for name in names:
        add(name)
    if path:
        with open(path) as fh:
            for line in fh:
                add(line)
    return requested


def query_interval(segment):
    return min(segment.query_start, segment.query_end) - 1, max(segment.query_start, segment.query_end)


def ref_interval(segment):
    return min(segment.ref_start, segment.ref_end) - 1, max(segment.ref_start, segment.ref_end)


def collect_blocks(
    alignment_path,
    alignment_format,
    requested,
    min_segment_bp,
    min_segment_idy,
    max_merge_gap,
    min_mapq=0,
    include_secondary_paf=False,
):
    raw_blocks = defaultdict(list)
    requested_set = set(requested) if requested is not None else None
    for segment in iter_alignments(
        alignment_path,
        alignment_format,
        min_identity=min_segment_idy,
        min_mapq=min_mapq,
        include_secondary_paf=include_secondary_paf,
    ):
        if requested_set is not None and segment.query not in requested_set:
            continue
        if segment.len_query < min_segment_bp:
            continue
        start, end = query_interval(segment)
        ref_start, ref_end = ref_interval(segment)
        raw_blocks[segment.query].append(
            QueryBlock(
                contig=segment.query,
                ref=segment.ref,
                ref_start=ref_start,
                ref_end=ref_end,
                query_start=start,
                query_end=end,
                orientation=segment.orientation,
                aligned_bp=end - start,
                identity_bp=segment.identity * (end - start),
                segment_count=1,
            )
        )

    merged = {}
    for contig, blocks in raw_blocks.items():
        merged[contig] = merge_query_blocks(blocks, max_merge_gap)
    return merged


def merge_query_blocks(blocks, max_merge_gap):
    merged = []
    for block in sorted(blocks, key=lambda item: (item.query_start, item.query_end, item.ref)):
        if not merged:
            merged.append(block)
            continue

        last = merged[-1]
        same_target = last.ref == block.ref and last.orientation == block.orientation
        close_enough = block.query_start <= last.query_end + max_merge_gap
        if same_target and close_enough:
            last.ref_start = min(last.ref_start, block.ref_start)
            last.ref_end = max(last.ref_end, block.ref_end)
            last.query_start = min(last.query_start, block.query_start)
            last.query_end = max(last.query_end, block.query_end)
            last.aligned_bp += block.aligned_bp
            last.identity_bp += block.identity_bp
            last.segment_count += block.segment_count
        else:
            merged.append(block)
    return merged


def split_signature(block):
    return block.ref, block.orientation


def auto_signature(block, include_orientation):
    if include_orientation:
        return split_signature(block)
    return block.ref, "."


def has_split_signal(blocks, include_orientation=True):
    if len(blocks) < 2:
        return False
    return len({auto_signature(block, include_orientation) for block in blocks}) > 1


def group_support(blocks, include_orientation=True):
    support = defaultdict(float)
    aligned_bp = defaultdict(int)
    identity_bp = defaultdict(float)
    ref_start = {}
    ref_end = {}
    segment_count = defaultdict(int)
    first_seen = {}

    for index, block in enumerate(blocks):
        signature = auto_signature(block, include_orientation)
        support[signature] += block.weighted_bp
        aligned_bp[signature] += block.aligned_bp
        identity_bp[signature] += block.identity_bp
        segment_count[signature] += block.segment_count
        ref_start[signature] = min(ref_start.get(signature, block.ref_start), block.ref_start)
        ref_end[signature] = max(ref_end.get(signature, block.ref_end), block.ref_end)
        first_seen.setdefault(signature, index)

    if not support:
        return None

    dominant = max(
        support,
        key=lambda signature: (
            support[signature],
            aligned_bp[signature],
            -first_seen[signature],
        ),
    )
    total_weighted_bp = sum(support.values())
    return {
        "signature": dominant,
        "dominant_weighted_bp": support[dominant],
        "total_weighted_bp": total_weighted_bp,
        "discordant_bp": total_weighted_bp - support[dominant],
        "aligned_bp": aligned_bp[dominant],
        "identity_bp": identity_bp[dominant],
        "ref_start": ref_start[dominant],
        "ref_end": ref_end[dominant],
        "segment_count": sum(segment_count.values()),
    }


def dominant_orientation(blocks, ref=None):
    relevant = [block for block in blocks if ref is None or block.ref == ref]
    plus_bp = sum(block.aligned_bp for block in relevant if block.orientation == "+")
    minus_bp = sum(block.aligned_bp for block in relevant if block.orientation == "-")
    return "-" if minus_bp > plus_bp else "+"


def summarize_group(blocks, include_orientation=True):
    support = group_support(blocks, include_orientation)
    ref, orientation = support["signature"]
    if not include_orientation:
        orientation = dominant_orientation(blocks, ref)
    query_start = min(block.query_start for block in blocks)
    query_end = max(block.query_end for block in blocks)
    summary = QueryBlock(
        contig=blocks[0].contig,
        ref=ref,
        ref_start=support["ref_start"],
        ref_end=support["ref_end"],
        query_start=query_start,
        query_end=query_end,
        orientation=orientation,
        aligned_bp=support["aligned_bp"],
        identity_bp=support["identity_bp"],
        segment_count=support["segment_count"],
    )
    return BlockGroup(
        blocks=blocks,
        summary=summary,
        discordant_bp=support["discordant_bp"],
        dominant_weighted_bp=support["dominant_weighted_bp"],
        total_weighted_bp=support["total_weighted_bp"],
    )


def collapse_adjacent_blocks_by_signature(blocks, include_orientation=True):
    collapsed = []
    current = []
    current_signature = None

    for block in blocks:
        signature = auto_signature(block, include_orientation)
        if current and signature != current_signature:
            collapsed.append(summarize_group(current, include_orientation).summary)
            current = []
        current.append(block)
        current_signature = signature

    if current:
        collapsed.append(summarize_group(current, include_orientation).summary)
    return collapsed


def auto_group_piece_is_supported(group, seq_len, args):
    span_frac = group.summary.query_span / seq_len if seq_len else 0.0
    return (
        group.summary.aligned_bp >= args.min_piece_aligned_bp
        and span_frac >= args.min_piece_query_frac
    )


def segment_blocks_for_auto(blocks, args, include_orientation):
    nblocks = len(blocks)
    if nblocks == 0:
        return []

    cost_cache = {}

    def interval_cost(start, end):
        key = (start, end)
        if key not in cost_cache:
            cost_cache[key] = group_support(
                blocks[start:end],
                include_orientation,
            )["discordant_bp"]
        return cost_cache[key]

    # Dynamic programming: retaining a discordant block inside a segment costs
    # its identity-weighted aligned bp; adding a breakpoint pays a fixed penalty.
    dp = [(float("inf"), float("inf"), None) for _ in range(nblocks + 1)]
    dp[0] = (0.0, 0, None)

    for end in range(1, nblocks + 1):
        best = dp[end]
        for start in range(0, end):
            previous_cost, previous_breaks, _ = dp[start]
            if previous_cost == float("inf"):
                continue
            breakpoint_cost = args.breakpoint_penalty_bp if start else 0.0
            breakpoint_count = previous_breaks + (1 if start else 0)
            candidate = (
                previous_cost + interval_cost(start, end) + breakpoint_cost,
                breakpoint_count,
                start,
            )
            if (candidate[0], candidate[1]) < (best[0], best[1]):
                best = candidate
        dp[end] = best

    groups = []
    end = nblocks
    while end > 0:
        start = dp[end][2]
        groups.append(summarize_group(blocks[start:end], include_orientation))
        end = start
    groups.reverse()
    return groups


def auto_plan_score(blocks, groups, args, include_orientation):
    if len(groups) < 2:
        return 0.0
    whole_discordant = group_support(blocks, include_orientation)["discordant_bp"]
    split_discordant = sum(group.discordant_bp for group in groups)
    breakpoints = len(groups) - 1
    return whole_discordant - split_discordant - (
        args.breakpoint_penalty_bp * breakpoints
    )


def ref_span_overlap_fraction(left, right):
    left_start, left_end = left.summary.ref_start, left.summary.ref_end
    right_start, right_end = right.summary.ref_start, right.summary.ref_end
    overlap = min(left_end, right_end) - max(left_start, right_start)
    if overlap <= 0:
        return 0.0
    smaller = min(left_end - left_start, right_end - right_start)
    return overlap / smaller if smaller > 0 else 0.0


def same_ref_orientation_plan_is_complex(groups, args):
    if len(groups) < 2:
        return False
    if len({group.summary.ref for group in groups}) != 1:
        return False
    if len({group.summary.orientation for group in groups}) < 2:
        return False

    for left, right in zip(groups, groups[1:]):
        if (
            left.summary.aligned_bp
            < args.complex_inversion_min_piece_aligned_bp
            or right.summary.aligned_bp
            < args.complex_inversion_min_piece_aligned_bp
        ):
            continue
        if (
            ref_span_overlap_fraction(left, right)
            >= args.complex_inversion_min_overlap_frac
        ):
            return True
    return False


def auto_orientation_groups(blocks, args):
    return merge_adjacent_auto_groups_by_signature(
        segment_blocks_for_auto(blocks, args, True),
        True,
    )


def auto_include_orientation(blocks, args):
    if not has_split_signal(blocks, True):
        return False
    if args.mode in {"comprehensive", "sensitive"}:
        return True
    if args.mode == "chromosome":
        return False
    groups = auto_orientation_groups(blocks, args)
    return same_ref_orientation_plan_is_complex(groups, args)


def merge_adjacent_auto_groups_by_signature(groups, include_orientation=True):
    merged = []
    for group in groups:
        if (
            merged
            and auto_signature(merged[-1].summary, include_orientation)
            == auto_signature(group.summary, include_orientation)
        ):
            combined = summarize_group(
                merged[-1].blocks + group.blocks,
                include_orientation,
            )
            if auto_signature(combined.summary, include_orientation) == auto_signature(
                merged[-1].summary,
                include_orientation,
            ):
                merged[-1] = combined
                continue
        merged.append(group)
    return merged


def merge_unsupported_terminal_groups(groups, seq_len, args, include_orientation):
    groups = list(groups)
    changed = False
    while len(groups) > 1 and not auto_group_piece_is_supported(groups[0], seq_len, args):
        groups[0:2] = [
            summarize_group(
                groups[0].blocks + groups[1].blocks,
                include_orientation,
            )
        ]
        changed = True
    while len(groups) > 1 and not auto_group_piece_is_supported(groups[-1], seq_len, args):
        groups[-2:] = [
            summarize_group(
                groups[-2].blocks + groups[-1].blocks,
                include_orientation,
            )
        ]
        changed = True
    if changed:
        groups = merge_adjacent_auto_groups_by_signature(groups, include_orientation)
    return groups, changed


def count_smoothed_transitions(groups, include_orientation=True):
    count = 0
    for group in groups:
        count += sum(
            1
            for left, right in zip(group.blocks, group.blocks[1:])
            if auto_signature(left, include_orientation)
            != auto_signature(right, include_orientation)
        )
    return count


def all_requested_contigs(fasta_path, blocks_by_contig, args):
    requested = []
    for name, _, _ in iter_fasta_records(fasta_path):
        blocks = blocks_by_contig.get(name, [])
        if args.mode == "sensitive":
            include = has_split_signal(blocks, True)
        else:
            include = has_split_signal(blocks, False) or auto_include_orientation(
                blocks,
                args,
            )
        if include:
            requested.append(name)
    return requested


def alpha_label(index):
    label = ""
    value = index
    while True:
        label = chr(ord("a") + (value % 26)) + label
        value = value // 26 - 1
        if value < 0:
            return label


def piece_name(ref, contig, part_index, separator):
    return f"{ref}{separator}{contig}{separator}{alpha_label(part_index)}"


def boundary_between(left, right):
    return max(0, (left.query_end + right.query_start) // 2)


def pieces_from_ordered_blocks(contig, seq_len, blocks, boundaries, args):
    starts = [0] + boundaries
    ends = boundaries + [seq_len]

    pieces = []
    for index, (block, start, end) in enumerate(zip(blocks, starts, ends)):
        start = max(0, min(start, seq_len))
        end = max(0, min(end, seq_len))
        if end - start < args.min_piece_bp:
            continue
        reverse_piece = args.orient_to_reference and block.orientation == "-"
        pieces.append(
            SplitPiece(
                original_contig=contig,
                new_name=piece_name(block.ref, contig, index, args.name_separator),
                part_index=index + 1,
                ref=block.ref,
                ref_start=block.ref_start,
                ref_end=block.ref_end,
                slice_start=start,
                slice_end=end,
                align_start=block.query_start,
                align_end=block.query_end,
                orientation=block.orientation,
                avg_identity=block.avg_identity,
                segment_count=block.segment_count,
                reverse_complemented=reverse_piece,
            )
        )
    return pieces


def build_split_plan(contig, seq_len, blocks, args):
    if not blocks:
        return ContigPlan(
            contig=contig,
            status="not_split_no_alignment",
            pieces=[],
            reason="No passing alignment segments were found for this requested contig.",
        )
    if len(blocks) == 1:
        return ContigPlan(
            contig=contig,
            status="not_split_single_block",
            pieces=[],
            reason="Only one passing alignment block was found.",
        )
    if not has_split_signal(blocks):
        return ContigPlan(
            contig=contig,
            status="not_split_single_target",
            pieces=[],
            reason="All passing alignment blocks map to the same reference sequence and orientation.",
        )

    split_blocks = collapse_adjacent_blocks_by_signature(blocks)
    if len(split_blocks) < 2:
        return ContigPlan(
            contig=contig,
            status="not_split_single_target",
            pieces=[],
            reason=(
                "After collapsing adjacent same-reference/orientation alignment "
                "runs, only one split target remained."
            ),
        )

    boundaries = [
        boundary_between(left, right)
        for left, right in zip(split_blocks, split_blocks[1:])
    ]
    pieces = pieces_from_ordered_blocks(contig, seq_len, split_blocks, boundaries, args)

    if len(pieces) < 2:
        return ContigPlan(
            contig=contig,
            status="not_split_short_piece",
            pieces=[],
            reason="Fewer than two split pieces passed the minimum piece length.",
        )

    return ContigPlan(
        contig=contig,
        status="split",
        pieces=pieces,
        reason=(
            f"{len(pieces)} query-ordered pieces inferred from "
            f"{len(split_blocks)} reference/orientation runs collapsed from "
            f"{len(blocks)} alignment blocks."
        ),
        planner_breakpoints=max(0, len(pieces) - 1),
        planner_ref_transition=len({piece.ref for piece in pieces}) > 1,
    )


def build_smoothed_split_plan(contig, seq_len, blocks, args):
    if not blocks:
        return ContigPlan(
            contig=contig,
            status="not_split_no_alignment",
            pieces=[],
            reason="No passing alignment segments were found for this candidate.",
        )
    if len(blocks) == 1:
        return ContigPlan(
            contig=contig,
            status="not_split_single_block",
            pieces=[],
            reason="Only one passing alignment block was found.",
        )

    include_orientation = auto_include_orientation(blocks, args)
    if not has_split_signal(blocks, False) and not include_orientation:
        return ContigPlan(
            contig=contig,
            status="not_split_single_target",
            pieces=[],
            reason=(
                "All passing alignment blocks map to the same reference sequence"
                + (
                    ". Same-reference orientation changes are ignored in this "
                    "mode unless they look complex enough to separate from "
                    "simple reference-relative inversions. Use --mode "
                    "comprehensive to consider all orientation changes or "
                    "--mode sensitive to disable smoothing."
                )
            ),
        )

    groups = merge_adjacent_auto_groups_by_signature(
        segment_blocks_for_auto(blocks, args, include_orientation),
        include_orientation,
    )
    initial_unsupported = [
        group
        for group in groups
        if not auto_group_piece_is_supported(group, seq_len, args)
    ]
    groups, merged_unsupported_terminal = merge_unsupported_terminal_groups(
        groups,
        seq_len,
        args,
        include_orientation,
    )
    if len(groups) < 2:
        discordant_bp = groups[0].discordant_bp if groups else 0.0
        if initial_unsupported and merged_unsupported_terminal:
            weakest = min(group.summary.aligned_bp for group in initial_unsupported)
            weakest_frac = min(
                group.summary.query_span / seq_len if seq_len else 0.0
                for group in initial_unsupported
            )
            return ContigPlan(
                contig=contig,
                status="not_split_smooth",
                pieces=[],
                reason=(
                    "Breakpoint smoothing rejected the best split because at least one "
                    "terminal piece failed support thresholds and was merged into "
                    "the neighboring piece "
                    f"({weakest} aligned bp; required {args.min_piece_aligned_bp}; "
                    f"query span fraction {weakest_frac:.4f}; "
                    f"required {args.min_piece_query_frac:.4f})."
                ),
            )
        return ContigPlan(
            contig=contig,
            status="not_split_smooth",
            pieces=[],
            reason=(
                "Breakpoint smoothing retained one piece: discordant support "
                f"({discordant_bp:.1f} identity-weighted bp) did not overcome "
                f"the breakpoint penalty ({args.breakpoint_penalty_bp:.1f})."
            ),
        )

    unsupported = [
        group
        for group in groups
        if not auto_group_piece_is_supported(group, seq_len, args)
    ]
    if unsupported:
        weakest = min(group.summary.aligned_bp for group in unsupported)
        weakest_frac = min(
            group.summary.query_span / seq_len if seq_len else 0.0
            for group in unsupported
        )
        return ContigPlan(
            contig=contig,
            status="not_split_smooth",
            pieces=[],
            reason=(
                "Breakpoint smoothing rejected the best split because at least one "
                "piece failed support thresholds "
                f"({weakest} aligned bp; required {args.min_piece_aligned_bp}; "
                f"query span fraction {weakest_frac:.4f}; "
                f"required {args.min_piece_query_frac:.4f})."
            ),
        )

    boundaries = [
        boundary_between(left.blocks[-1], right.blocks[0])
        for left, right in zip(groups, groups[1:])
    ]
    summary_blocks = [group.summary for group in groups]
    pieces = pieces_from_ordered_blocks(contig, seq_len, summary_blocks, boundaries, args)

    if len(pieces) < 2:
        return ContigPlan(
            contig=contig,
            status="not_split_short_piece",
            pieces=[],
            reason="Fewer than two smoothed split pieces passed the minimum piece length.",
        )

    breakpoints = len(pieces) - 1
    score = auto_plan_score(blocks, groups, args, include_orientation)
    ref_transition = len({group.summary.ref for group in groups}) > 1
    smoothed = count_smoothed_transitions(groups, include_orientation)
    return ContigPlan(
        contig=contig,
        status="split",
        pieces=pieces,
        reason=(
            f"{len(pieces)} smoothed pieces inferred from {len(blocks)} alignment blocks "
            f"after smoothing {smoothed} weak transition(s) with breakpoint penalty "
            f"{args.breakpoint_penalty_bp:.1f}; planner_score={score:.1f}."
        ),
        planner_score=score,
        planner_breakpoints=breakpoints,
        planner_ref_transition=ref_transition,
    )


def build_plans(fasta_path, requested, blocks_by_contig, args):
    seq_lengths = {name: len(seq) for name, _, seq in iter_fasta_records(fasta_path)}
    plans: Dict[str, ContigPlan] = {}
    for contig in requested:
        if contig not in seq_lengths:
            plans[contig] = ContigPlan(
                contig=contig,
                status="target_missing",
                pieces=[],
                reason="Requested contig was not found in the assembly FASTA.",
            )
            continue
        if args.mode == "sensitive":
            plans[contig] = build_split_plan(
                contig,
                seq_lengths[contig],
                blocks_by_contig.get(contig, []),
                args,
            )
        else:
            plans[contig] = build_smoothed_split_plan(
                contig,
                seq_lengths[contig],
                blocks_by_contig.get(contig, []),
                args,
            )
    return plans


def apply_breakpoint_guard(plans, args):
    if args.max_breakpoints_per_contig < 0:
        return

    for plan in plans.values():
        if plan.status != "split":
            continue
        breakpoints = plan.planner_breakpoints or max(0, len(plan.pieces) - 1)
        if breakpoints <= args.max_breakpoints_per_contig:
            continue

        plan.status = "not_split_too_many_breakpoints"
        plan.pieces = []
        plan.reason = (
            "Split plan was rejected by the per-contig breakpoint guardrail: "
            f"this plan required {breakpoints} breakpoint(s), which is above "
            f"--max-breakpoints-per-contig {args.max_breakpoints_per_contig}. "
            f"planner_score={plan.planner_score:.1f}."
        )


def fasta_header(piece, simple_headers):
    if simple_headers:
        return piece.new_name
    fields = [
        piece.new_name,
        f"original={piece.original_contig}",
        f"ref={piece.ref}",
        f"slice={piece.slice_start + 1}-{piece.slice_end}",
        f"alignment={piece.align_start + 1}-{piece.align_end}",
        f"orientation={piece.orientation}",
        f"reverse_complemented={'yes' if piece.reverse_complemented else 'no'}",
        f"avg_identity={piece.avg_identity:.3f}",
    ]
    return " ".join(fields)


def write_fixed_fasta(path, fasta_path, plans, args):
    with open(path, "w") as out:
        for name, header, seq in iter_fasta_records(fasta_path):
            plan = plans.get(name)
            if plan is None:
                if not args.pieces_only:
                    out.write(header + "\n")
                    write_wrapped(out, seq)
                continue

            if plan.status != "split":
                if not args.pieces_only:
                    out.write(header + "\n")
                    write_wrapped(out, seq)
                continue

            for piece in plan.pieces:
                piece_seq = seq[piece.slice_start : piece.slice_end]
                if piece.reverse_complemented:
                    piece_seq = reverse_complement(piece_seq)
                out.write(f">{fasta_header(piece, args.simple_headers)}\n")
                write_wrapped(out, piece_seq)


def fmt(value, digits=3):
    if value is None:
        return "."
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def fmt_bool(value):
    if value is None:
        return "."
    return "yes" if value else "no"


def write_report(path, requested, plans):
    header = [
        "original_contig",
        "status",
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
        "reason",
    ]
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for contig in requested:
            plan = plans[contig]
            if plan.status != "split":
                row = [
                    contig,
                    plan.status,
                    ".",
                    ".",
                    ".",
                    ".",
                    ".",
                    ".",
                    ".",
                    ".",
                    ".",
                    ".",
                    ".",
                    ".",
                    plan.reason,
                ]
                out.write("\t".join(str(item) for item in row) + "\n")
                continue
            for piece in plan.pieces:
                row = [
                    contig,
                    plan.status,
                    piece.new_name,
                    piece.part_index,
                    piece.ref,
                    piece.slice_start + 1,
                    piece.slice_end,
                    piece.length,
                    piece.align_start + 1,
                    piece.align_end,
                    piece.orientation,
                    "yes" if piece.reverse_complemented else "no",
                    fmt(piece.avg_identity),
                    piece.segment_count,
                    plan.reason,
                ]
                out.write("\t".join(str(item) for item in row) + "\n")


def default_graph_report_path(report_path):
    return Path(report_path).with_suffix(".graph.tsv")


def graph_fix_note(plan, node_evidence):
    present = node_evidence.graph_node_status == "present"
    if plan.status == "split":
        return (
            "split_source_present_in_graph_review"
            if present
            else "split_source_missing_from_graph"
        )
    return "graph_context_only" if present else "source_missing_from_graph"


def write_graph_report(path, requested, plans, graph):
    header = [
        "original_contig",
        "status",
        "graph_node",
        "graph_node_status",
        "graph_node_length",
        "graph_node_has_sequence",
        "graph_in_degree",
        "graph_out_degree",
        "graph_neighbor_count",
        "graph_self_loop",
        "split_pieces",
        "planner_breakpoints",
        "planner_ref_transition",
        "planner_score",
        "graph_note",
    ]
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for contig in requested:
            plan = plans[contig]
            node_evidence = graph_node_evidence(graph, [contig])
            split_pieces = (
                ",".join(piece.new_name for piece in plan.pieces)
                if plan.status == "split"
                else "."
            )
            row = [
                contig,
                plan.status,
                node_evidence.graph_node,
                node_evidence.graph_node_status,
                fmt(node_evidence.graph_node_length, 0),
                fmt_bool(node_evidence.graph_node_has_sequence),
                fmt(node_evidence.graph_in_degree, 0),
                fmt(node_evidence.graph_out_degree, 0),
                fmt(node_evidence.graph_neighbor_count, 0),
                fmt_bool(node_evidence.graph_self_loop),
                split_pieces,
                plan.planner_breakpoints,
                fmt_bool(plan.planner_ref_transition),
                fmt(plan.planner_score, 1),
                graph_fix_note(plan, node_evidence),
            ]
            out.write("\t".join(str(item) for item in row) + "\n")


def graph_guard_fix_warnings(requested, plans, graph, stream=None):
    stream = stream or sys.stderr
    for contig in requested:
        plan = plans[contig]
        node_evidence = graph_node_evidence(graph, [contig])
        if node_evidence.graph_node_status != "present":
            continue
        neighbor_count = node_evidence.graph_neighbor_count or 0
        self_loop = bool(node_evidence.graph_self_loop)
        is_simple = neighbor_count <= 2 and not self_loop
        is_complex = neighbor_count > 2 or self_loop
        if plan.status == "split" and is_simple:
            stream.write(
                "WARNING: graph guard: "
                f"{contig} is planned for splitting, but its GFA node "
                f"{node_evidence.graph_node} has a simple neighborhood "
                f"(neighbors={neighbor_count}, self_loop=no). Review the "
                "breakpoints before applying this edit.\n"
            )
        elif plan.status != "split" and is_complex:
            stream.write(
                "WARNING: graph guard: "
                f"{contig} was not split, but its GFA node "
                f"{node_evidence.graph_node} is graph-complex "
                f"(neighbors={neighbor_count}, self_loop={'yes' if self_loop else 'no'}). "
                "Manual graph review may be useful.\n"
            )


def main(argv: Optional[Sequence[str]] = None, prog: Optional[str] = None):
    args = parse_args(argv, prog=prog)
    if args.graph_report and not args.gfa:
        sys.stderr.write("ERROR: --graph-report requires --gfa\n")
        sys.exit(2)
    if args.graph_guard and not args.gfa:
        sys.stderr.write("ERROR: --graph-guard requires --gfa\n")
        sys.exit(2)
    alignment_path, alignment_format = alignment_source_from_args(args)
    explicit_requested = read_requested_contigs(args.contigs, args.contigs_file)
    if args.all_contigs and explicit_requested:
        sys.stderr.write("ERROR: use either --all or --contigs/--contigs-file, not both\n")
        sys.exit(2)
    if not explicit_requested and not args.all_contigs:
        sys.stderr.write("ERROR: provide at least one contig via --contigs/--contigs-file or use --all\n")
        sys.exit(2)

    graph_report_path = None
    if args.gfa:
        graph_report_path = (
            Path(args.graph_report)
            if args.graph_report
            else default_graph_report_path(args.report)
        )

    output_paths = [Path(args.output_fasta), Path(args.report)]
    if graph_report_path is not None:
        output_paths.append(graph_report_path)
    for output_path in output_paths:
        if output_path.parent and str(output_path.parent) != ".":
            output_path.parent.mkdir(parents=True, exist_ok=True)

    collect_for = None if args.all_contigs else explicit_requested
    blocks_by_contig = collect_blocks(
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
        all_requested_contigs(args.assembly_fasta, blocks_by_contig, args)
        if args.all_contigs
        else explicit_requested
    )
    plans = build_plans(
        args.assembly_fasta,
        requested,
        blocks_by_contig,
        args,
    )
    apply_breakpoint_guard(plans, args)
    write_fixed_fasta(args.output_fasta, args.assembly_fasta, plans, args)
    write_report(args.report, requested, plans)
    if args.gfa:
        graph = read_gfa(args.gfa)
        write_graph_report(graph_report_path, requested, plans, graph)
        if args.graph_guard:
            graph_guard_fix_warnings(requested, plans, graph)

    status_counts = defaultdict(int)
    for plan in plans.values():
        status_counts[plan.status] += 1
    sys.stderr.write(f"Processed {len(requested)} requested contigs.\n")
    for status, count in sorted(status_counts.items()):
        sys.stderr.write(f"  {status}: {count}\n")
    sys.stderr.write(f"Wrote fixed FASTA: {args.output_fasta}\n")
    sys.stderr.write(f"Wrote report: {args.report}\n")
    if args.gfa:
        sys.stderr.write(f"Wrote graph report: {graph_report_path}\n")


if __name__ == "__main__":
    main()
