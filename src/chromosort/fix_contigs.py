#!/usr/bin/env python3
"""
Split user-nominated chimeric contigs using MUMmer show-coords alignments.

This module is intentionally conservative by default. User-nominated contigs
are split directly from query-ordered alignment blocks. With --auto, candidate
blocks are first passed through a breakpoint-penalty segmentation step so small
discordant blocks, INDEL-sized gaps, and local SV-like noise are smoothed over.
"""

import argparse
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

from .reference_order import (
    iter_coords,
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


@dataclass
class BlockGroup:
    blocks: List[QueryBlock]
    summary: QueryBlock
    discordant_bp: float
    dominant_weighted_bp: float
    total_weighted_bp: float


def parse_args(argv: Optional[Sequence[str]] = None, prog: Optional[str] = None):
    ap = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Split user-nominated chimeric contigs into reference-labeled pieces "
            "using MUMmer show-coords alignments."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument(
        "-f",
        "--assembly-fasta",
        required=True,
        help="Assembly FASTA containing the contigs to fix.",
    )
    ap.add_argument(
        "-c",
        "--coords",
        required=True,
        help="MUMmer show-coords file for reference-vs-assembly alignment.",
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
        "--auto",
        action="store_true",
        help=(
            "Automatically inspect contigs with passing alignment blocks that "
            "change reference sequence or orientation. Auto mode uses "
            "breakpoint-penalty smoothing unless --auto-sensitive is set."
        ),
    )
    ap.add_argument(
        "--auto-sensitive",
        action="store_true",
        help=(
            "In auto mode, split every passing reference/orientation transition "
            "without breakpoint-penalty smoothing. This is useful for debugging "
            "or intentionally sensitive scans."
        ),
    )
    ap.add_argument(
        "-o",
        "--output-fasta",
        required=True,
        help="Output FASTA with nominated chimeras replaced by split pieces.",
    )
    ap.add_argument(
        "--report",
        required=True,
        help="TSV report describing split pieces and unsplit requested contigs.",
    )
    ap.add_argument(
        "--min-segment-bp",
        type=int,
        default=10_000,
        help="Minimum show-coords LEN2 for an alignment segment to inform splitting.",
    )
    ap.add_argument(
        "--min-segment-idy",
        type=float,
        default=0.0,
        help="Minimum percent identity for an alignment segment to inform splitting.",
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
        "--auto-breakpoint-penalty-bp",
        type=float,
        default=50_000.0,
        help=(
            "Auto-mode breakpoint penalty, measured as identity-weighted aligned bp. "
            "A split is kept only when doing so explains more discordant support "
            "than this penalty."
        ),
    )
    ap.add_argument(
        "--auto-min-piece-aligned-bp",
        type=int,
        default=50_000,
        help=(
            "Minimum dominant aligned bp required in each auto-split piece. "
            "This prevents auto mode from cutting off weak local SV or repeat hits."
        ),
    )
    ap.add_argument(
        "--auto-min-piece-query-frac",
        type=float,
        default=0.01,
        help=(
            "Minimum query-span fraction required for each auto-split piece. "
            "Set to 0 to disable this fraction check."
        ),
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
    return ap.parse_args(argv)


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


def collect_blocks(coords_path, requested, min_segment_bp, min_segment_idy, max_merge_gap):
    raw_blocks = defaultdict(list)
    requested_set = set(requested) if requested is not None else None
    for segment in iter_coords(coords_path, min_segment_idy):
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


def has_split_signal(blocks):
    if len(blocks) < 2:
        return False
    return len({split_signature(block) for block in blocks}) > 1


def group_support(blocks):
    support = defaultdict(float)
    aligned_bp = defaultdict(int)
    identity_bp = defaultdict(float)
    ref_start = {}
    ref_end = {}
    segment_count = defaultdict(int)
    first_seen = {}

    for index, block in enumerate(blocks):
        signature = split_signature(block)
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


def summarize_group(blocks):
    support = group_support(blocks)
    ref, orientation = support["signature"]
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


def auto_group_piece_is_supported(group, seq_len, args):
    span_frac = group.summary.query_span / seq_len if seq_len else 0.0
    return (
        group.summary.aligned_bp >= args.auto_min_piece_aligned_bp
        and span_frac >= args.auto_min_piece_query_frac
    )


def segment_blocks_for_auto(blocks, args):
    nblocks = len(blocks)
    if nblocks == 0:
        return []

    cost_cache = {}

    def interval_cost(start, end):
        key = (start, end)
        if key not in cost_cache:
            cost_cache[key] = group_support(blocks[start:end])["discordant_bp"]
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
            breakpoint_cost = args.auto_breakpoint_penalty_bp if start else 0.0
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
        groups.append(summarize_group(blocks[start:end]))
        end = start
    groups.reverse()
    return groups


def count_smoothed_transitions(groups):
    count = 0
    for group in groups:
        count += sum(
            1
            for left, right in zip(group.blocks, group.blocks[1:])
            if split_signature(left) != split_signature(right)
        )
    return count


def auto_requested_contigs(fasta_path, blocks_by_contig, explicit_contigs):
    explicit_set = set(explicit_contigs)
    requested = list(explicit_contigs)
    for name, _, _ in iter_fasta_records(fasta_path):
        if name in explicit_set:
            continue
        blocks = blocks_by_contig.get(name, [])
        if has_split_signal(blocks):
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

    boundaries = [boundary_between(left, right) for left, right in zip(blocks, blocks[1:])]
    pieces = pieces_from_ordered_blocks(contig, seq_len, blocks, boundaries, args)

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
        reason=f"{len(pieces)} query-ordered pieces inferred from {len(blocks)} alignment blocks.",
    )


def build_auto_split_plan(contig, seq_len, blocks, args):
    if args.auto_sensitive:
        return build_split_plan(contig, seq_len, blocks, args)

    if not blocks:
        return ContigPlan(
            contig=contig,
            status="not_split_no_alignment",
            pieces=[],
            reason="No passing alignment segments were found for this auto candidate.",
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

    groups = segment_blocks_for_auto(blocks, args)
    if len(groups) < 2:
        discordant_bp = groups[0].discordant_bp if groups else 0.0
        return ContigPlan(
            contig=contig,
            status="not_split_auto_smooth",
            pieces=[],
            reason=(
                "Auto smoothing retained one piece: discordant support "
                f"({discordant_bp:.1f} identity-weighted bp) did not overcome "
                f"the breakpoint penalty ({args.auto_breakpoint_penalty_bp:.1f})."
            ),
        )

    unsupported = [
        group
        for group in groups
        if not auto_group_piece_is_supported(group, seq_len, args)
    ]
    if unsupported:
        weakest = min(group.summary.aligned_bp for group in unsupported)
        return ContigPlan(
            contig=contig,
            status="not_split_auto_smooth",
            pieces=[],
            reason=(
                "Auto smoothing rejected the best split because at least one "
                "piece had weak dominant support "
                f"({weakest} aligned bp; required {args.auto_min_piece_aligned_bp})."
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
            reason="Fewer than two auto split pieces passed the minimum piece length.",
        )

    smoothed = count_smoothed_transitions(groups)
    return ContigPlan(
        contig=contig,
        status="split",
        pieces=pieces,
        reason=(
            f"{len(pieces)} auto pieces inferred from {len(blocks)} alignment blocks "
            f"after smoothing {smoothed} weak transition(s) with breakpoint penalty "
            f"{args.auto_breakpoint_penalty_bp:.1f}."
        ),
    )


def build_plans(fasta_path, requested, explicit_requested, blocks_by_contig, args):
    seq_lengths = {name: len(seq) for name, _, seq in iter_fasta_records(fasta_path)}
    plans: Dict[str, ContigPlan] = {}
    explicit_set = set(explicit_requested)
    for contig in requested:
        if contig not in seq_lengths:
            plans[contig] = ContigPlan(
                contig=contig,
                status="target_missing",
                pieces=[],
                reason="Requested contig was not found in the assembly FASTA.",
            )
            continue
        if args.auto and contig not in explicit_set:
            plans[contig] = build_auto_split_plan(
                contig,
                seq_lengths[contig],
                blocks_by_contig.get(contig, []),
                args,
            )
        else:
            plans[contig] = build_split_plan(
                contig,
                seq_lengths[contig],
                blocks_by_contig.get(contig, []),
                args,
            )
    return plans


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


def main(argv: Optional[Sequence[str]] = None, prog: Optional[str] = None):
    args = parse_args(argv, prog=prog)
    explicit_requested = read_requested_contigs(args.contigs, args.contigs_file)
    if not explicit_requested and not args.auto:
        sys.stderr.write("ERROR: provide at least one contig via --contigs/--contigs-file or use --auto\n")
        sys.exit(2)

    for output_path in [Path(args.output_fasta), Path(args.report)]:
        if output_path.parent and str(output_path.parent) != ".":
            output_path.parent.mkdir(parents=True, exist_ok=True)

    collect_for = None if args.auto else explicit_requested
    blocks_by_contig = collect_blocks(
        args.coords,
        collect_for,
        args.min_segment_bp,
        args.min_segment_idy,
        args.max_merge_gap,
    )
    requested = (
        auto_requested_contigs(args.assembly_fasta, blocks_by_contig, explicit_requested)
        if args.auto
        else explicit_requested
    )
    plans = build_plans(
        args.assembly_fasta,
        requested,
        explicit_requested,
        blocks_by_contig,
        args,
    )
    write_fixed_fasta(args.output_fasta, args.assembly_fasta, plans, args)
    write_report(args.report, requested, plans)

    status_counts = defaultdict(int)
    for plan in plans.values():
        status_counts[plan.status] += 1
    sys.stderr.write(f"Processed {len(requested)} requested contigs.\n")
    for status, count in sorted(status_counts.items()):
        sys.stderr.write(f"  {status}: {count}\n")
    sys.stderr.write(f"Wrote fixed FASTA: {args.output_fasta}\n")
    sys.stderr.write(f"Wrote report: {args.report}\n")


if __name__ == "__main__":
    main()
