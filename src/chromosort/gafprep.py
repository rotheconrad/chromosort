"""Prepare targeted GraphAligner inputs from ChromoSort review evidence."""

import argparse
import csv
import gzip
import hashlib
import heapq
import os
import shlex
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .graph import parse_gfa_overlap_bp, parse_gfa_tags
from .paths import ensure_output_dirs, ensure_parent_dir
from .reference_order import parse_paf_line, parse_paf_tags, read_fasta_lengths


REVIEW_TYPES = {"fix", "scaffold", "gapfill"}
MISSING_VALUES = {"", "."}

TARGET_COLUMNS = [
    "target_id",
    "review_table",
    "review_type",
    "review_row_id",
    "contig",
    "start",
    "end",
    "target_type",
    "reason",
    "source_fields",
]

SELECTED_READ_COLUMNS = [
    "read_id",
    "read_length",
    "assembly_contig",
    "aln_start",
    "aln_end",
    "strand",
    "mapq",
    "aln_block_length",
    "matches",
    "identity_estimate",
    "selection_reasons",
    "target_ids",
    "review_row_ids",
    "is_decision_target",
    "is_background",
]

LINK_COLUMNS = [
    "read_id",
    "review_table",
    "review_type",
    "review_row_id",
    "event_type",
    "contig",
    "target_id",
    "selection_reason",
    "distance_to_target",
    "overlap_bp",
]

DROPPED_LINK_COLUMNS = [
    "line_number",
    "node_a",
    "orientation_a",
    "node_b",
    "orientation_b",
    "node_a_length",
    "node_b_length",
    "cigar",
    "overlap_len",
    "reason",
    "near_selected_target",
]


@dataclass(frozen=True)
class TargetInterval:
    """One decision-relevant assembly interval used to select reads."""

    target_id: str
    review_table: str
    review_type: str
    review_row_id: str
    event_type: str
    contig: str
    start: int
    end: int
    target_type: str
    reason: str
    source_fields: str


@dataclass(frozen=True)
class ReadPafAlignment:
    """One read-to-assembly PAF row with fields needed for audit output."""

    read_id: str
    read_length: int
    read_start: int
    read_end: int
    contig: str
    contig_length: int
    aln_start: int
    aln_end: int
    strand: str
    matches: int
    aln_block_length: int
    mapq: int
    identity_estimate: float
    is_secondary: bool
    line_number: int


@dataclass(frozen=True)
class ReadTargetHit:
    """Association between one selected read and one review-derived target."""

    read_id: str
    review_table: str
    review_type: str
    review_row_id: str
    event_type: str
    contig: str
    target_id: str
    selection_reason: str
    distance_to_target: int
    overlap_bp: int


@dataclass
class SelectedRead:
    """Deduplicated selected read plus all retained reasons and links."""

    read_id: str
    read_length: int = 0
    best_alignment: Optional[ReadPafAlignment] = None
    selection_reasons: set = field(default_factory=set)
    target_ids: set = field(default_factory=set)
    review_row_ids: set = field(default_factory=set)
    is_decision_target: bool = False
    is_background: bool = False
    links: Dict[Tuple[str, str, str], ReadTargetHit] = field(default_factory=dict)

    def consider_alignment(self, alignment: ReadPafAlignment):
        self.read_length = max(self.read_length, alignment.read_length)
        if self.best_alignment is None or alignment_rank(alignment) > alignment_rank(
            self.best_alignment
        ):
            self.best_alignment = alignment

    def add_reason(self, reason: str, alignment: Optional[ReadPafAlignment] = None):
        if reason:
            self.selection_reasons.add(reason)
        if reason == "background_bin":
            self.is_background = True
        else:
            self.is_decision_target = True
        if alignment is not None:
            self.consider_alignment(alignment)

    def add_link(self, hit: ReadTargetHit):
        key = (hit.target_id, hit.review_row_id, hit.contig)
        existing = self.links.get(key)
        if existing is None or hit_rank(hit) > hit_rank(existing):
            self.links[key] = hit
        self.target_ids.add(hit.target_id)
        self.review_row_ids.add(hit.review_row_id)


@dataclass
class PafScanStats:
    """Transparent counters for the read-to-assembly PAF scan."""

    data_rows: int = 0
    parseable_rows: int = 0
    malformed_rows: int = 0
    kept_rows: int = 0
    filtered_secondary_rows: int = 0
    filtered_mapq_rows: int = 0
    filtered_aligned_bp_rows: int = 0


@dataclass
class SelectionResult:
    """Read selection result plus scan counters."""

    selected: Dict[str, SelectedRead]
    paf_stats: PafScanStats
    target_hit_reads: int
    background_bins_seen: int
    background_reads_selected: int


def parse_args(argv=None, prog=None):
    ap = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Prepare targeted GraphAligner inputs from read-to-assembly PAF "
            "and ChromoSort review tables. For broad review runs, generate "
            "fix, scaffold, and gapfill tables together with chromo eval all, "
            "then pass all three tables here."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--assembly-fasta", required=True, help="Assembly FASTA used as the read-PAF target.")
    ap.add_argument("--assembly-gfa", required=True, help="Assembly graph GFA to sanitize for GraphAligner.")
    ap.add_argument("--read-paf", required=True, help="Long-read-to-assembly PAF.")
    ap.add_argument("--reads", required=True, help="Input FASTQ/FASTQ.GZ containing the reads.")
    ap.add_argument(
        "--eval-review-table",
        action="append",
        required=True,
        help=(
            "ChromoSort eval/review TSV. May be supplied more than once; "
            "the recommended broad workflow is the three tables from chromo eval all."
        ),
    )
    ap.add_argument(
        "--review-type",
        choices=["auto", "fix", "scaffold", "gapfill"],
        default="auto",
        help="Override review table type inference for all supplied tables.",
    )
    ap.add_argument("-o", "--output-prefix", required=True, help="Output prefix for audit tables and inputs.")
    ap.add_argument("--target-padding", type=int, default=50_000)
    ap.add_argument("--contig-end-window", type=int, default=100_000)
    ap.add_argument("--target-reads-per-interval", type=int, default=20)
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-aligned-bp", type=int, default=5_000)
    ap.add_argument("--include-secondary-paf", action="store_true")
    ap.add_argument("--background-bin-size", type=int, default=1_000_000)
    ap.add_argument("--background-reads-per-bin", type=int, default=1)
    ap.add_argument("--max-reads", type=int, default=None)
    ap.add_argument("--max-reads-per-contig", type=int, default=None)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--graphaligner-threads", type=int, default=16)
    ap.add_argument("--graphaligner-preset", default="vg")
    ap.add_argument("--precise-clipping", default="0.9")
    ap.add_argument("--graphaligner-bin", default="GraphAligner")
    return ap.parse_args(argv)


def validate_args(args):
    nonnegative = [
        ("--target-padding", args.target_padding),
        ("--contig-end-window", args.contig_end_window),
        ("--min-mapq", args.min_mapq),
        ("--min-aligned-bp", args.min_aligned_bp),
        ("--background-reads-per-bin", args.background_reads_per_bin),
    ]
    for name, value in nonnegative:
        if value < 0:
            raise ValueError(f"{name} must be zero or greater")
    positive = [
        ("--target-reads-per-interval", args.target_reads_per_interval),
        ("--background-bin-size", args.background_bin_size),
        ("--graphaligner-threads", args.graphaligner_threads),
    ]
    for name, value in positive:
        if value < 1:
            raise ValueError(f"{name} must be at least 1")
    if args.max_reads is not None and args.max_reads < 1:
        raise ValueError("--max-reads must be at least 1 when supplied")
    if args.max_reads_per_contig is not None and args.max_reads_per_contig < 1:
        raise ValueError("--max-reads-per-contig must be at least 1 when supplied")


def open_text(path):
    path = Path(path)
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def write_tsv(path, columns, rows):
    ensure_parent_dir(path)
    with open(path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, ".") for column in columns})


def present(value):
    return value is not None and str(value).strip() not in MISSING_VALUES


def first_present(row, names):
    for name in names:
        value = row.get(name)
        if present(value):
            return str(value).strip()
    return None


def parse_int_value(value):
    if not present(value):
        return None
    text = str(value).strip().replace(",", "")
    try:
        return int(text)
    except ValueError:
        try:
            return int(float(text))
        except ValueError:
            return None


def parse_position_list(value):
    if not present(value):
        return []
    positions = []
    for token in str(value).replace(";", ",").replace("|", ",").split(","):
        token = token.strip()
        if not token:
            continue
        parsed = parse_int_value(token)
        if parsed is not None:
            positions.append(parsed)
    return positions


def source_fields_text(row, names):
    fields = []
    for name in names:
        value = row.get(name)
        if present(value):
            fields.append(f"{name}={value}")
    return ";".join(fields) if fields else "."


def infer_review_type(path, fieldnames, row, override):
    if override != "auto":
        return override
    task = row.get("task")
    if task in REVIEW_TYPES:
        return task

    lowered_path = Path(path).name.lower()
    for review_type in sorted(REVIEW_TYPES):
        if review_type in lowered_path:
            return review_type

    columns = set(fieldnames or [])
    if {"source_contig", "slice_start", "slice_end"} & columns or "planner_breakpoints" in columns:
        return "fix"
    if {"fill_status", "path_nodes", "fill_sequence", "candidate_path_count"} & columns:
        return "gapfill"
    if {"left_contig", "right_contig", "gap_bp", "gap_mode"} & columns:
        return "scaffold"
    return "fix"


def review_row_id(row, row_number):
    return first_present(row, ["event_id", "review_row_id", "row_id", "id"]) or f"row{row_number}"


def review_event_type(row):
    return first_present(row, ["action", "event_type", "task", "status"]) or "."


class TargetBuilder:
    """Build target intervals while keeping IDs stable and table-scoped."""

    def __init__(self, contig_lengths, target_padding, contig_end_window):
        self.contig_lengths = contig_lengths
        self.target_padding = target_padding
        self.contig_end_window = contig_end_window
        self.targets = []
        self.next_index = 1

    def add_interval(
        self,
        review_table,
        review_type,
        row_id,
        event_type,
        contig,
        start,
        end,
        target_type,
        reason,
        source_fields,
    ):
        if not contig or contig == ".":
            return
        length = self.contig_lengths.get(contig)
        if start is None or end is None:
            return
        start = max(1, int(start))
        end = max(start, int(end))
        if length is not None:
            start = min(start, length)
            end = min(end, length)
            if end < start:
                return
        target_id = f"target_{self.next_index:06d}"
        self.next_index += 1
        self.targets.append(
            TargetInterval(
                target_id=target_id,
                review_table=review_table,
                review_type=review_type,
                review_row_id=row_id,
                event_type=event_type,
                contig=contig,
                start=start,
                end=end,
                target_type=target_type,
                reason=reason or ".",
                source_fields=source_fields or ".",
            )
        )

    def add_position(
        self,
        review_table,
        review_type,
        row_id,
        event_type,
        contig,
        position,
        target_type,
        reason,
        source_fields,
    ):
        if position is None:
            return
        self.add_interval(
            review_table,
            review_type,
            row_id,
            event_type,
            contig,
            position - self.target_padding,
            position + self.target_padding,
            target_type,
            reason,
            source_fields,
        )

    def add_contig_end(
        self,
        review_table,
        review_type,
        row_id,
        event_type,
        contig,
        side,
        target_type,
        reason,
        source_fields,
    ):
        length = self.contig_lengths.get(contig)
        if length is None:
            return
        window = min(self.contig_end_window, length)
        if side == "left":
            start, end = 1, window
        elif side == "right":
            start, end = length - window + 1, length
        else:
            raise ValueError(f"unknown contig side: {side!r}")
        self.add_interval(
            review_table,
            review_type,
            row_id,
            event_type,
            contig,
            start,
            end,
            target_type,
            reason,
            source_fields,
        )

    def add_both_ends(
        self,
        review_table,
        review_type,
        row_id,
        event_type,
        contig,
        target_type_prefix,
        reason,
        source_fields,
    ):
        self.add_contig_end(
            review_table,
            review_type,
            row_id,
            event_type,
            contig,
            "left",
            f"{target_type_prefix}_left_end",
            reason,
            source_fields,
        )
        self.add_contig_end(
            review_table,
            review_type,
            row_id,
            event_type,
            contig,
            "right",
            f"{target_type_prefix}_right_end",
            reason,
            source_fields,
        )


def build_fix_targets(builder, review_table, review_type, row_id, row):
    contig = first_present(row, ["source_contig", "target", "contig", "assembly_contig"])
    event_type = review_event_type(row)
    reason = first_present(row, ["reason", "status", "planner_status", "action"]) or "fix_review"
    source_names = [
        "source_contig",
        "target",
        "action",
        "status",
        "reason",
        "planner_status",
        "planner_breakpoints",
        "longread_breakpoint_position",
        "slice_start",
        "slice_end",
        "alignment_query_start",
        "alignment_query_end",
        "graph_unitig",
        "gaf_node",
        "longread_spanning_reads",
        "longread_split_reads",
    ]
    source_fields = source_fields_text(row, source_names)
    positions = []
    for name in ["longread_breakpoint_position", "breakpoint_position", "breakpoint", "position"]:
        value = parse_int_value(row.get(name))
        if value is not None:
            positions.append((value, name))
    for value in parse_position_list(row.get("planner_breakpoints")):
        positions.append((value, "planner_breakpoints"))

    length = builder.contig_lengths.get(contig)
    slice_start = parse_int_value(row.get("slice_start"))
    slice_end = parse_int_value(row.get("slice_end"))
    if slice_start is not None and slice_start > 1:
        positions.append((slice_start, "slice_start"))
    if slice_end is not None and (length is None or slice_end < length):
        positions.append((slice_end, "slice_end"))
    if not positions:
        for name in ["alignment_query_start", "alignment_query_end"]:
            value = parse_int_value(row.get(name))
            if value is not None:
                positions.append((value, name))

    seen_positions = set()
    for position, source_name in positions:
        if position in seen_positions:
            continue
        seen_positions.add(position)
        builder.add_position(
            review_table,
            review_type,
            row_id,
            event_type,
            contig,
            position,
            "fix_breakpoint",
            f"{reason};source={source_name}",
            source_fields,
        )

    if not seen_positions and contig:
        builder.add_both_ends(
            review_table,
            review_type,
            row_id,
            event_type,
            contig,
            "fix_fallback",
            f"{reason};fallback=contig_end",
            source_fields,
        )


def build_junction_targets(builder, review_table, review_type, row_id, row):
    left = first_present(row, ["left_contig", "left_component", "left", "upstream_contig"])
    right = first_present(row, ["right_contig", "right_component", "right", "downstream_contig"])
    event_type = review_event_type(row)
    reason = first_present(row, ["reason", "fill_status", "gaf_support_status", "graph_status", "action"])
    if not reason:
        reason = f"{review_type}_review"
    source_names = [
        "scaffold",
        "left_contig",
        "right_contig",
        "left_component",
        "right_component",
        "action",
        "status",
        "reason",
        "fill_status",
        "graph_status",
        "graph_direct_edge",
        "graph_path_nodes",
        "path_nodes",
        "gaf_path_nodes",
        "gaf_support_status",
        "longread_bridge_reads",
    ]
    source_fields = source_fields_text(row, source_names)
    prefix = "gapfill" if review_type == "gapfill" else "scaffold"

    added = False
    if left:
        builder.add_contig_end(
            review_table,
            review_type,
            row_id,
            event_type,
            left,
            "right",
            f"{prefix}_left_flank",
            f"{reason};flank=left",
            source_fields,
        )
        added = True
    if right:
        builder.add_contig_end(
            review_table,
            review_type,
            row_id,
            event_type,
            right,
            "left",
            f"{prefix}_right_flank",
            f"{reason};flank=right",
            source_fields,
        )
        added = True

    if not added:
        contig = first_present(row, ["target", "contig", "source_contig", "assembly_contig"])
        if contig:
            builder.add_both_ends(
                review_table,
                review_type,
                row_id,
                event_type,
                contig,
                f"{prefix}_fallback",
                f"{reason};fallback=contig_end",
                source_fields,
            )


def build_targets(review_tables, review_type_override, contig_lengths, args):
    builder = TargetBuilder(contig_lengths, args.target_padding, args.contig_end_window)
    for table in review_tables:
        with open(table, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if reader.fieldnames is None:
                raise ValueError(f"Review table is empty: {table}")
            for row_number, row in enumerate(reader, start=1):
                review_type = infer_review_type(table, reader.fieldnames, row, review_type_override)
                row_id = review_row_id(row, row_number)
                review_table = Path(table).name
                if review_type == "fix":
                    build_fix_targets(builder, review_table, review_type, row_id, row)
                elif review_type in {"scaffold", "gapfill"}:
                    build_junction_targets(builder, review_table, review_type, row_id, row)
                else:
                    raise ValueError(f"Unsupported review type for {table}: {review_type!r}")
    return builder.targets


def target_to_row(target):
    return {
        "target_id": target.target_id,
        "review_table": target.review_table,
        "review_type": target.review_type,
        "review_row_id": target.review_row_id,
        "contig": target.contig,
        "start": str(target.start),
        "end": str(target.end),
        "target_type": target.target_type,
        "reason": target.reason,
        "source_fields": target.source_fields,
    }


def parse_read_paf(path, min_mapq, min_aligned_bp, include_secondary, stats=None):
    if stats is None:
        stats = PafScanStats()
    with open_text(path) as fh:
        for line_number, line in enumerate(fh, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            stats.data_rows += 1
            segment = parse_paf_line(line)
            cols = line.rstrip("\n").split("\t")
            if segment is None or len(cols) < 12:
                stats.malformed_rows += 1
                continue
            stats.parseable_rows += 1
            tags = parse_paf_tags(cols)
            is_secondary = tags.get("tp") == ("A", "S")
            if is_secondary and not include_secondary:
                stats.filtered_secondary_rows += 1
                continue
            try:
                matches = int(cols[9])
                block_length = int(cols[10])
                mapq = int(cols[11])
            except ValueError:
                stats.malformed_rows += 1
                continue
            if mapq < min_mapq:
                stats.filtered_mapq_rows += 1
                continue
            if block_length < min_aligned_bp:
                stats.filtered_aligned_bp_rows += 1
                continue
            stats.kept_rows += 1
            yield ReadPafAlignment(
                read_id=segment.query,
                read_length=segment.query_length,
                read_start=segment.query_start,
                read_end=segment.query_end,
                contig=segment.ref,
                contig_length=segment.ref_length,
                aln_start=segment.ref_start,
                aln_end=segment.ref_end,
                strand=segment.orientation,
                matches=matches,
                aln_block_length=block_length,
                mapq=mapq,
                identity_estimate=segment.identity,
                is_secondary=is_secondary,
                line_number=line_number,
            )


def read_paf_target_lengths(path):
    """Return target-name lengths advertised by a read-to-assembly PAF."""

    lengths = {}
    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 7:
                continue
            try:
                length = int(cols[6])
            except ValueError:
                continue
            lengths.setdefault(cols[5], length)
    return lengths


def interval_overlap_bp(a_start, a_end, b_start, b_end):
    left = max(a_start, b_start)
    right = min(a_end, b_end)
    if right < left:
        return 0
    return right - left + 1


def interval_distance(a_start, a_end, b_start, b_end):
    if interval_overlap_bp(a_start, a_end, b_start, b_end):
        return 0
    if a_end < b_start:
        return b_start - a_end
    return a_start - b_end


def alignment_rank(alignment):
    return (
        alignment.mapq,
        alignment.aln_block_length,
        alignment.matches,
        alignment.identity_estimate,
        -alignment.line_number,
        alignment.read_id,
    )


def candidate_rank(alignment, hit):
    return (
        hit.overlap_bp > 0,
        hit.overlap_bp,
        -hit.distance_to_target,
        alignment.mapq,
        alignment.aln_block_length,
        alignment.matches,
        alignment.identity_estimate,
        -alignment.line_number,
        alignment.read_id,
    )


def hit_rank(hit):
    return (
        hit.overlap_bp > 0,
        hit.overlap_bp,
        -hit.distance_to_target,
        hit.target_id,
    )


def stable_hash_int(seed, *parts):
    text = "\t".join([str(seed), *(str(part) for part in parts)])
    digest = hashlib.sha256(text.encode("utf-8")).hexdigest()
    return int(digest[:16], 16)


def add_background_candidate(reservoir, bin_key, alignment, args):
    hash_value = stable_hash_int(args.seed, "background", bin_key[0], bin_key[1], alignment.read_id)
    item = (-hash_value, -alignment.line_number, alignment)
    heap = reservoir[bin_key]
    if len(heap) < args.background_reads_per_bin:
        heapq.heappush(heap, item)
    elif item > heap[0]:
        heapq.heapreplace(heap, item)


def build_read_hit(alignment, target):
    overlap = interval_overlap_bp(alignment.aln_start, alignment.aln_end, target.start, target.end)
    distance = interval_distance(alignment.aln_start, alignment.aln_end, target.start, target.end)
    reason = "target_overlap" if overlap else "target_nearby"
    return ReadTargetHit(
        read_id=alignment.read_id,
        review_table=target.review_table,
        review_type=target.review_type,
        review_row_id=target.review_row_id,
        event_type=target.event_type,
        contig=target.contig,
        target_id=target.target_id,
        selection_reason=reason,
        distance_to_target=distance,
        overlap_bp=overlap,
    )


def select_reads(read_paf, targets, args):
    targets_by_contig = defaultdict(list)
    for target in targets:
        targets_by_contig[target.contig].append(target)

    target_candidates: Dict[str, Dict[str, Tuple[tuple, ReadPafAlignment]]] = defaultdict(dict)
    read_target_hits: Dict[str, Dict[Tuple[str, str, str], ReadTargetHit]] = defaultdict(dict)
    background_reservoir = defaultdict(list)
    stats = PafScanStats()

    for alignment in parse_read_paf(
        read_paf,
        args.min_mapq,
        args.min_aligned_bp,
        args.include_secondary_paf,
        stats=stats,
    ):
        if args.background_reads_per_bin:
            midpoint = (alignment.aln_start + alignment.aln_end) // 2
            bin_index = max(0, (midpoint - 1) // args.background_bin_size)
            add_background_candidate(background_reservoir, (alignment.contig, bin_index), alignment, args)

        for target in targets_by_contig.get(alignment.contig, []):
            overlap = interval_overlap_bp(alignment.aln_start, alignment.aln_end, target.start, target.end)
            if overlap <= 0:
                continue
            hit = build_read_hit(alignment, target)
            read_hit_key = (hit.target_id, hit.review_row_id, hit.contig)
            existing_hit = read_target_hits[alignment.read_id].get(read_hit_key)
            if existing_hit is None or hit_rank(hit) > hit_rank(existing_hit):
                read_target_hits[alignment.read_id][read_hit_key] = hit

            rank = candidate_rank(alignment, hit)
            existing = target_candidates[target.target_id].get(alignment.read_id)
            if existing is None or rank > existing[0]:
                target_candidates[target.target_id][alignment.read_id] = (rank, alignment)

    selected: Dict[str, SelectedRead] = {}

    def ensure_selected(read_id):
        if read_id not in selected:
            selected[read_id] = SelectedRead(read_id=read_id)
        return selected[read_id]

    for target_id in sorted(target_candidates):
        candidates = sorted(
            target_candidates[target_id].values(),
            key=lambda item: item[0],
            reverse=True,
        )
        for _rank, alignment in candidates[: args.target_reads_per_interval]:
            summary = ensure_selected(alignment.read_id)
            summary.add_reason("target_interval", alignment)

    for read_id, hit_map in read_target_hits.items():
        target_ids = {hit.target_id for hit in hit_map.values()}
        contigs = {hit.contig for hit in hit_map.values()}
        if len(target_ids) >= 2 or len(contigs) >= 2:
            ensure_selected(read_id).add_reason("bridge_two_targets")

    background_selected = 0
    for bin_key in sorted(background_reservoir):
        for _hash_value, _line_number, alignment in sorted(background_reservoir[bin_key], reverse=True):
            summary = ensure_selected(alignment.read_id)
            summary.add_reason("background_bin", alignment)
            background_selected += 1

    for alignment in parse_read_paf(
        read_paf,
        args.min_mapq,
        args.min_aligned_bp,
        args.include_secondary_paf,
    ):
        summary = selected.get(alignment.read_id)
        if summary is not None:
            summary.consider_alignment(alignment)

    for read_id, summary in selected.items():
        for hit in read_target_hits.get(read_id, {}).values():
            summary.add_link(hit)

    selected = enforce_read_limits(selected, args)
    return SelectionResult(
        selected=selected,
        paf_stats=stats,
        target_hit_reads=len(read_target_hits),
        background_bins_seen=len(background_reservoir),
        background_reads_selected=background_selected,
    )


def selected_rank(summary):
    aln = summary.best_alignment
    aln_rank = alignment_rank(aln) if aln is not None else (-1, -1, -1, -1.0, 0, summary.read_id)
    return (
        summary.is_decision_target,
        len(summary.target_ids),
        not summary.is_background,
        *aln_rank,
        summary.read_id,
    )


def sorted_selected(selected):
    return sorted(selected.values(), key=selected_rank, reverse=True)


def enforce_read_limits(selected, args):
    if args.max_reads is None and args.max_reads_per_contig is None:
        return selected
    kept = {}
    per_contig = defaultdict(int)
    for summary in sorted_selected(selected):
        contig = summary.best_alignment.contig if summary.best_alignment is not None else "."
        if args.max_reads_per_contig is not None and per_contig[contig] >= args.max_reads_per_contig:
            continue
        kept[summary.read_id] = summary
        per_contig[contig] += 1
        if args.max_reads is not None and len(kept) >= args.max_reads:
            break
    return kept


def selected_read_row(summary):
    aln = summary.best_alignment
    reasons = ";".join(sorted(summary.selection_reasons)) if summary.selection_reasons else "."
    target_ids = ";".join(sorted(summary.target_ids)) if summary.target_ids else "."
    review_ids = ";".join(sorted(summary.review_row_ids)) if summary.review_row_ids else "."
    return {
        "read_id": summary.read_id,
        "read_length": str(summary.read_length or (aln.read_length if aln else 0)),
        "assembly_contig": aln.contig if aln else ".",
        "aln_start": str(aln.aln_start) if aln else ".",
        "aln_end": str(aln.aln_end) if aln else ".",
        "strand": aln.strand if aln else ".",
        "mapq": str(aln.mapq) if aln else ".",
        "aln_block_length": str(aln.aln_block_length) if aln else ".",
        "matches": str(aln.matches) if aln else ".",
        "identity_estimate": f"{aln.identity_estimate:.3f}" if aln else ".",
        "selection_reasons": reasons,
        "target_ids": target_ids,
        "review_row_ids": review_ids,
        "is_decision_target": "yes" if summary.is_decision_target else "no",
        "is_background": "yes" if summary.is_background else "no",
    }


def link_rows(summary):
    rows = []
    bridge = "bridge_two_targets" in summary.selection_reasons
    background = "background_bin" in summary.selection_reasons and not summary.is_decision_target
    for hit in summary.links.values():
        reasons = [hit.selection_reason]
        if bridge:
            reasons.append("bridge_two_targets")
        if background:
            reasons.append("background_bin")
        rows.append(
            {
                "read_id": summary.read_id,
                "review_table": hit.review_table,
                "review_type": hit.review_type,
                "review_row_id": hit.review_row_id,
                "event_type": hit.event_type,
                "contig": hit.contig,
                "target_id": hit.target_id,
                "selection_reason": ";".join(dict.fromkeys(reasons)),
                "distance_to_target": str(hit.distance_to_target),
                "overlap_bp": str(hit.overlap_bp),
            }
        )
    return sorted(rows, key=lambda row: (row["read_id"], row["target_id"], row["review_row_id"]))


def write_selected_outputs(paths, selected):
    ordered = sorted_selected(selected)
    write_tsv(paths["selected_reads"], SELECTED_READ_COLUMNS, [selected_read_row(row) for row in ordered])

    links = []
    for summary in ordered:
        links.extend(link_rows(summary))
    write_tsv(paths["selected_read_review_links"], LINK_COLUMNS, links)

    ensure_parent_dir(paths["selected_read_ids"])
    with open(paths["selected_read_ids"], "w") as out:
        for summary in ordered:
            out.write(f"{summary.read_id}\n")


def iter_fastq_records(path):
    with open_text(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            sequence = fh.readline()
            plus = fh.readline()
            quality = fh.readline()
            if not quality:
                raise ValueError(f"Malformed FASTQ record near header {header.rstrip()!r}")
            if not header.startswith("@"):
                raise ValueError(f"Malformed FASTQ header: {header.rstrip()!r}")
            read_id = header[1:].strip().split()[0]
            yield read_id, header, sequence, plus, quality


def extract_fastq(reads_path, output_path, selected_ids):
    ensure_parent_dir(output_path)
    selected_ids = set(selected_ids)
    found = set()
    with gzip.open(output_path, "wt") as out:
        for read_id, header, sequence, plus, quality in iter_fastq_records(reads_path):
            if read_id not in selected_ids:
                continue
            out.write(header)
            out.write(sequence)
            out.write(plus)
            out.write(quality)
            found.add(read_id)
    missing = sorted(selected_ids - found)
    return len(found), missing


def gfa_node_lengths(path):
    lengths = {}
    with open_text(path) as fh:
        for line_number, line in enumerate(fh, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] != "S" or len(fields) < 3:
                continue
            name = fields[1]
            sequence = fields[2]
            tags = parse_gfa_tags(fields[3:])
            tagged_length = tags.get("LN")
            if sequence != "*":
                lengths[name] = len(sequence)
            elif isinstance(tagged_length, int):
                lengths[name] = tagged_length
            else:
                lengths[name] = None
    return lengths


def sanitize_gfa(input_path, output_path, dropped_links_path, target_contigs):
    ensure_parent_dir(output_path)
    ensure_parent_dir(dropped_links_path)
    lengths = gfa_node_lengths(input_path)
    summary = defaultdict(int)
    target_contigs = set(target_contigs)
    dropped_rows = []

    with open_text(input_path) as inp, open(output_path, "w") as out:
        for line_number, line in enumerate(inp, start=1):
            stripped = line.rstrip("\n")
            if not stripped or stripped.startswith("#"):
                out.write(line)
                summary["kept_lines"] += 1
                continue
            fields = stripped.split("\t")
            record_type = fields[0]
            summary[f"input_{record_type}_records"] += 1

            if record_type == "A":
                summary["removed_A_records"] += 1
                continue

            if record_type == "L" and len(fields) >= 6:
                node_a = fields[1]
                orientation_a = fields[2]
                node_b = fields[3]
                orientation_b = fields[4]
                cigar = fields[5]
                try:
                    overlap = parse_gfa_overlap_bp(cigar)
                except ValueError:
                    overlap = None
                    summary["malformed_link_overlap_cigars_kept"] += 1
                len_a = lengths.get(node_a)
                len_b = lengths.get(node_b)
                consumes_a = overlap is not None and len_a is not None and overlap >= len_a
                consumes_b = overlap is not None and len_b is not None and overlap >= len_b
                if consumes_a or consumes_b:
                    reasons = []
                    if consumes_a:
                        reasons.append("overlap_consumes_node_a")
                    if consumes_b:
                        reasons.append("overlap_consumes_node_b")
                    near_target = node_a in target_contigs or node_b in target_contigs
                    if near_target:
                        summary["dropped_links_near_selected_target"] += 1
                    summary["dropped_L_records"] += 1
                    dropped_rows.append(
                        {
                            "line_number": str(line_number),
                            "node_a": node_a,
                            "orientation_a": orientation_a,
                            "node_b": node_b,
                            "orientation_b": orientation_b,
                            "node_a_length": "." if len_a is None else str(len_a),
                            "node_b_length": "." if len_b is None else str(len_b),
                            "cigar": cigar,
                            "overlap_len": str(overlap),
                            "reason": ";".join(reasons),
                            "near_selected_target": "yes" if near_target else "no",
                        }
                    )
                    continue

            out.write(line)
            summary["kept_lines"] += 1
            summary[f"kept_{record_type}_records"] += 1

    write_tsv(dropped_links_path, DROPPED_LINK_COLUMNS, dropped_rows)
    summary["dropped_links"] = summary["dropped_L_records"]
    summary["target_contigs"] = len(target_contigs)
    if summary["dropped_links_near_selected_target"]:
        summary["gaf_evidence_limited_by_graph_sanitization_near_target"] = 1
    else:
        summary["gaf_evidence_limited_by_graph_sanitization_near_target"] = 0
    if summary["dropped_links"] >= 10:
        summary["many_dropped_links"] = 1
    else:
        summary["many_dropped_links"] = 0
    return dict(summary)


def write_gfa_sanitize_summary(path, input_path, output_path, summary):
    rows = [
        {"metric": "input_gfa", "value": str(input_path)},
        {"metric": "output_gfa", "value": str(output_path)},
    ]
    for key in sorted(summary):
        rows.append({"metric": key, "value": str(summary[key])})
    if summary.get("gaf_evidence_limited_by_graph_sanitization_near_target"):
        rows.append(
            {
                "metric": "warning",
                "value": "GAF evidence limited by graph sanitization near target",
            }
        )
    elif summary.get("many_dropped_links"):
        rows.append(
            {
                "metric": "warning",
                "value": "GAF evidence may be limited by extensive graph sanitization",
            }
        )
    write_tsv(path, ["metric", "value"], rows)


def shell_quote(value):
    return shlex.quote(str(value))


def write_graphaligner_script(path, args, graph_gfa, reads_fastq, output_gaf):
    ensure_parent_dir(path)
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        "# Targeted GraphAligner run prepared by chromo gafprep.",
        "# This script aligns only the selected read subset; edit scheduler",
        "# directives, module loads, or thread counts as needed for your HPC.",
        f"{shell_quote(args.graphaligner_bin)} \\",
        f"  -g {shell_quote(graph_gfa)} \\",
        f"  -f {shell_quote(reads_fastq)} \\",
        f"  -a {shell_quote(output_gaf)} \\",
        f"  -t {shell_quote(args.graphaligner_threads)} \\",
        f"  -x {shell_quote(args.graphaligner_preset)} \\",
        f"  --precise-clipping {shell_quote(args.precise_clipping)}",
        "",
    ]
    with open(path, "w") as out:
        out.write("\n".join(lines))
    mode = os.stat(path).st_mode
    os.chmod(path, mode | 0o111)


def summary_rows(args, targets, selection, extracted_count, missing_reads, gfa_summary, paths):
    stats = selection.paf_stats
    rows = [
        {"metric": "assembly_fasta", "value": str(args.assembly_fasta)},
        {"metric": "assembly_gfa", "value": str(args.assembly_gfa)},
        {"metric": "read_paf", "value": str(args.read_paf)},
        {"metric": "reads", "value": str(args.reads)},
        {"metric": "review_tables", "value": ",".join(args.eval_review_table)},
        {"metric": "target_count", "value": str(len(targets))},
        {"metric": "selected_read_count", "value": str(len(selection.selected))},
        {"metric": "extracted_read_count", "value": str(extracted_count)},
        {"metric": "missing_selected_read_count", "value": str(len(missing_reads))},
        {"metric": "missing_selected_reads", "value": ",".join(missing_reads) if missing_reads else "."},
        {"metric": "paf_data_rows", "value": str(stats.data_rows)},
        {"metric": "paf_parseable_rows", "value": str(stats.parseable_rows)},
        {"metric": "paf_malformed_rows", "value": str(stats.malformed_rows)},
        {"metric": "paf_kept_rows", "value": str(stats.kept_rows)},
        {"metric": "paf_filtered_secondary_rows", "value": str(stats.filtered_secondary_rows)},
        {"metric": "paf_filtered_mapq_rows", "value": str(stats.filtered_mapq_rows)},
        {"metric": "paf_filtered_aligned_bp_rows", "value": str(stats.filtered_aligned_bp_rows)},
        {"metric": "target_hit_reads", "value": str(selection.target_hit_reads)},
        {"metric": "background_bins_seen", "value": str(selection.background_bins_seen)},
        {"metric": "background_reads_selected", "value": str(selection.background_reads_selected)},
        {"metric": "selected_fastq", "value": str(paths["selected_fastq"])},
        {"metric": "graphaligner_gfa", "value": str(paths["graphaligner_gfa"])},
        {"metric": "graphaligner_script", "value": str(paths["graphaligner_script"])},
    ]
    for key in sorted(gfa_summary):
        rows.append({"metric": f"gfa_{key}", "value": str(gfa_summary[key])})
    if gfa_summary.get("gaf_evidence_limited_by_graph_sanitization_near_target"):
        rows.append(
            {
                "metric": "warning",
                "value": "GAF evidence limited by graph sanitization near target",
            }
        )
    elif gfa_summary.get("many_dropped_links"):
        rows.append(
            {
                "metric": "warning",
                "value": "GAF evidence may be limited by extensive graph sanitization",
            }
        )
    return rows


def output_paths(prefix):
    prefix = Path(prefix)
    return {
        "targets": Path(str(prefix) + ".targets.tsv"),
        "selected_reads": Path(str(prefix) + ".selected_reads.tsv"),
        "selected_read_review_links": Path(str(prefix) + ".selected_read_review_links.tsv"),
        "selected_read_ids": Path(str(prefix) + ".selected_read_ids.txt"),
        "selected_fastq": Path(str(prefix) + ".selected.fastq.gz"),
        "graphaligner_gfa": Path(str(prefix) + ".graphaligner.gfa"),
        "graphaligner_script": Path(str(prefix) + ".graphaligner.sh"),
        "summary": Path(str(prefix) + ".summary.tsv"),
        "dropped_gfa_links": Path(str(prefix) + ".dropped_gfa_links.tsv"),
        "gfa_sanitize_summary": Path(str(prefix) + ".gfa_sanitize_summary.tsv"),
        "graphaligner_gaf": Path(str(prefix) + ".gaf"),
    }


def run(args):
    validate_args(args)
    paths = output_paths(args.output_prefix)
    ensure_output_dirs(paths)
    _records, contig_record_map = read_fasta_lengths(args.assembly_fasta)
    contig_lengths = {name: record.length for name, record in contig_record_map.items()}
    for contig, length in read_paf_target_lengths(args.read_paf).items():
        contig_lengths.setdefault(contig, length)

    targets = build_targets(args.eval_review_table, args.review_type, contig_lengths, args)
    write_tsv(paths["targets"], TARGET_COLUMNS, [target_to_row(target) for target in targets])

    selection = select_reads(args.read_paf, targets, args)
    write_selected_outputs(paths, selection.selected)
    selected_ids = list(selection.selected)
    extracted_count, missing_reads = extract_fastq(args.reads, paths["selected_fastq"], selected_ids)

    target_contigs = {target.contig for target in targets}
    gfa_summary = sanitize_gfa(
        args.assembly_gfa,
        paths["graphaligner_gfa"],
        paths["dropped_gfa_links"],
        target_contigs,
    )
    write_gfa_sanitize_summary(
        paths["gfa_sanitize_summary"],
        args.assembly_gfa,
        paths["graphaligner_gfa"],
        gfa_summary,
    )
    write_graphaligner_script(
        paths["graphaligner_script"],
        args,
        paths["graphaligner_gfa"],
        paths["selected_fastq"],
        paths["graphaligner_gaf"],
    )
    write_tsv(
        paths["summary"],
        ["metric", "value"],
        summary_rows(args, targets, selection, extracted_count, missing_reads, gfa_summary, paths),
    )

    sys.stderr.write(f"Wrote targeted GraphAligner prep outputs with prefix: {args.output_prefix}\n")
    sys.stderr.write(f"  targets: {len(targets)}\n")
    sys.stderr.write(f"  selected reads: {len(selection.selected)}\n")
    sys.stderr.write(f"  extracted reads: {extracted_count}\n")
    if missing_reads:
        sys.stderr.write(f"  missing selected reads in FASTQ: {len(missing_reads)}\n")
    if gfa_summary.get("gaf_evidence_limited_by_graph_sanitization_near_target"):
        sys.stderr.write("  WARNING: GAF evidence limited by graph sanitization near target\n")


def main(argv: Optional[Sequence[str]] = None, prog: Optional[str] = None):
    try:
        run(parse_args(argv, prog=prog))
    except ValueError as exc:
        sys.stderr.write(f"ERROR: {exc}\n")
        sys.exit(2)


if __name__ == "__main__":
    main()
