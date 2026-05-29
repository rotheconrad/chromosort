"""Long-read backmapping evidence helpers."""

from collections import Counter, defaultdict
from dataclasses import dataclass, field
from statistics import median
from typing import Dict, Iterable, List, Optional, Sequence

from .reference_order import iter_paf


@dataclass(frozen=True)
class LongReadAlignment:
    """One long-read alignment against an assembly contig."""

    read: str
    read_length: int
    read_start: int
    read_end: int
    contig: str
    contig_length: int
    contig_start: int
    contig_end: int
    orientation: str
    identity: float
    mapq: Optional[int]
    is_secondary: bool = False

    @property
    def read_span_bp(self):
        return self.read_end - self.read_start + 1

    @property
    def contig_span_bp(self):
        return self.contig_end - self.contig_start + 1


@dataclass
class LongReadEvidence:
    """Long-read alignments indexed by read and assembly contig."""

    alignments: List[LongReadAlignment]
    by_read: Dict[str, List[LongReadAlignment]] = field(default_factory=dict)
    by_contig: Dict[str, List[LongReadAlignment]] = field(default_factory=dict)

    @classmethod
    def from_alignments(cls, alignments: Iterable[LongReadAlignment]):
        ordered = list(alignments)
        by_read = defaultdict(list)
        by_contig = defaultdict(list)
        for aln in ordered:
            by_read[aln.read].append(aln)
            by_contig[aln.contig].append(aln)
        for grouped in by_read.values():
            grouped.sort(key=lambda item: (item.read_start, item.read_end, item.contig))
        for grouped in by_contig.values():
            grouped.sort(key=lambda item: (item.contig_start, item.contig_end, item.read))
        return cls(
            alignments=ordered,
            by_read=dict(by_read),
            by_contig=dict(by_contig),
        )


@dataclass(frozen=True)
class BreakpointSupport:
    """Read evidence around one candidate cut after ``position``."""

    contig: str
    position: int
    window_bp: int
    min_anchor_bp: int
    spanning_reads: Sequence[str]
    split_reads: Sequence[str]
    left_edge_reads: Sequence[str]
    right_edge_reads: Sequence[str]
    nearby_reads: Sequence[str]

    @property
    def spanning_count(self):
        return len(self.spanning_reads)

    @property
    def split_count(self):
        return len(self.split_reads)

    @property
    def left_edge_count(self):
        return len(self.left_edge_reads)

    @property
    def right_edge_count(self):
        return len(self.right_edge_reads)

    @property
    def nearby_count(self):
        return len(self.nearby_reads)

    def as_dict(self):
        return {
            "contig": self.contig,
            "position": str(self.position),
            "window_bp": str(self.window_bp),
            "min_anchor_bp": str(self.min_anchor_bp),
            "spanning_reads": str(self.spanning_count),
            "split_reads": str(self.split_count),
            "left_edge_reads": str(self.left_edge_count),
            "right_edge_reads": str(self.right_edge_count),
            "nearby_reads": str(self.nearby_count),
        }


@dataclass(frozen=True)
class ReadBridge:
    """One read linking terminal alignments on two contigs."""

    read: str
    left_orientation: str
    right_orientation: str
    read_order: str
    read_gap_bp: int

    @property
    def orientation_key(self):
        return f"{self.left_orientation}/{self.right_orientation}"


@dataclass(frozen=True)
class ContigBridgeSupport:
    """Read evidence connecting terminal regions of two contigs."""

    left_contig: str
    right_contig: str
    left_side: str
    right_side: str
    terminal_window_bp: int
    min_anchor_bp: int
    bridges: Sequence[ReadBridge]

    @property
    def bridge_count(self):
        return len(self.bridges)

    @property
    def orientation_summary(self):
        counts = Counter(bridge.orientation_key for bridge in self.bridges)
        if not counts:
            return "."
        return ",".join(f"{key}:{counts[key]}" for key in sorted(counts))

    @property
    def read_order_summary(self):
        counts = Counter(bridge.read_order for bridge in self.bridges)
        if not counts:
            return "."
        return ",".join(f"{key}:{counts[key]}" for key in sorted(counts))

    @property
    def median_read_gap_bp(self):
        if not self.bridges:
            return None
        value = median(bridge.read_gap_bp for bridge in self.bridges)
        return int(value) if float(value).is_integer() else value

    @property
    def bridge_reads(self):
        return tuple(bridge.read for bridge in self.bridges)

    def as_dict(self):
        median_gap = self.median_read_gap_bp
        return {
            "left_contig": self.left_contig,
            "right_contig": self.right_contig,
            "left_side": self.left_side,
            "right_side": self.right_side,
            "terminal_window_bp": str(self.terminal_window_bp),
            "min_anchor_bp": str(self.min_anchor_bp),
            "bridge_reads": str(self.bridge_count),
            "orientation_summary": self.orientation_summary,
            "read_order_summary": self.read_order_summary,
            "median_read_gap_bp": "." if median_gap is None else str(median_gap),
        }


def read_long_read_paf(
    path,
    min_identity=0.0,
    min_mapq=0,
    include_secondary=False,
):
    """Read long-read-to-assembly PAF into indexed evidence."""

    alignments = []
    for segment in iter_paf(
        path,
        min_identity=min_identity,
        min_mapq=min_mapq,
        include_secondary=include_secondary,
    ):
        alignments.append(
            LongReadAlignment(
                read=segment.query,
                read_length=segment.query_length,
                read_start=segment.query_start,
                read_end=segment.query_end,
                contig=segment.ref,
                contig_length=segment.ref_length,
                contig_start=segment.ref_start,
                contig_end=segment.ref_end,
                orientation=segment.orientation,
                identity=segment.identity,
                mapq=segment.mapq,
                is_secondary=segment.is_secondary,
            )
        )
    return LongReadEvidence.from_alignments(alignments)


def _overlap_bp(start, end, query_start, query_end):
    left = max(start, query_start)
    right = min(end, query_end)
    if right < left:
        return 0
    return right - left + 1


def _sorted_names(values):
    return tuple(sorted(values))


def summarize_breakpoint(
    evidence: LongReadEvidence,
    contig,
    position,
    window_bp=5_000,
    min_anchor_bp=1_000,
):
    """Summarize long-read support near a candidate cut after ``position``."""

    if position < 1:
        raise ValueError("position must be 1 or greater")
    if window_bp < 0:
        raise ValueError("window_bp must be zero or greater")
    if min_anchor_bp < 1:
        raise ValueError("min_anchor_bp must be at least 1")

    spanning = set()
    left_edge = set()
    right_edge = set()
    nearby = set()

    left_anchor_start = max(1, position - min_anchor_bp + 1)
    left_anchor_end = position
    right_anchor_start = position + 1
    right_anchor_end = position + min_anchor_bp

    for aln in evidence.by_contig.get(contig, []):
        left_anchor = _overlap_bp(
            aln.contig_start,
            aln.contig_end,
            left_anchor_start,
            left_anchor_end,
        )
        right_anchor = _overlap_bp(
            aln.contig_start,
            aln.contig_end,
            right_anchor_start,
            right_anchor_end,
        )
        if left_anchor >= min_anchor_bp and right_anchor >= min_anchor_bp:
            spanning.add(aln.read)

        if _overlap_bp(
            aln.contig_start,
            aln.contig_end,
            max(1, position - window_bp),
            position + window_bp,
        ):
            nearby.add(aln.read)
        if abs(aln.contig_end - position) <= window_bp:
            left_edge.add(aln.read)
        if abs(aln.contig_start - (position + 1)) <= window_bp:
            right_edge.add(aln.read)

    split = (left_edge & right_edge) - spanning
    return BreakpointSupport(
        contig=contig,
        position=position,
        window_bp=window_bp,
        min_anchor_bp=min_anchor_bp,
        spanning_reads=_sorted_names(spanning),
        split_reads=_sorted_names(split),
        left_edge_reads=_sorted_names(left_edge - spanning),
        right_edge_reads=_sorted_names(right_edge - spanning),
        nearby_reads=_sorted_names(nearby),
    )


def _terminal_anchor_bp(aln, side, terminal_window_bp):
    if side == "start":
        return _overlap_bp(aln.contig_start, aln.contig_end, 1, terminal_window_bp)
    if side == "end":
        start = max(1, aln.contig_length - terminal_window_bp + 1)
        return _overlap_bp(aln.contig_start, aln.contig_end, start, aln.contig_length)
    raise ValueError(f"invalid contig side: {side!r}")


def _read_gap(left, right):
    if left.read_start <= right.read_start:
        return right.read_start - left.read_end - 1
    return left.read_start - right.read_end - 1


def _read_order(left, right):
    if left.read_end < right.read_start:
        return "left_before_right"
    if right.read_end < left.read_start:
        return "right_before_left"
    return "overlap_on_read"


def _best_bridge_pair(left_alignments, right_alignments):
    best = None
    best_key = None
    for left in left_alignments:
        for right in right_alignments:
            key = (abs(_read_gap(left, right)), -min(left.mapq or 0, right.mapq or 0))
            if best is None or key < best_key:
                best = (left, right)
                best_key = key
    return best


def summarize_contig_bridge(
    evidence: LongReadEvidence,
    left_contig,
    right_contig,
    left_side="end",
    right_side="start",
    terminal_window_bp=5_000,
    min_anchor_bp=500,
):
    """Summarize reads connecting terminal regions of two contigs."""

    if terminal_window_bp < 1:
        raise ValueError("terminal_window_bp must be at least 1")
    if min_anchor_bp < 1:
        raise ValueError("min_anchor_bp must be at least 1")

    bridges = []
    for read, alignments in evidence.by_read.items():
        left_hits = [
            aln
            for aln in alignments
            if aln.contig == left_contig
            and _terminal_anchor_bp(aln, left_side, terminal_window_bp) >= min_anchor_bp
        ]
        right_hits = [
            aln
            for aln in alignments
            if aln.contig == right_contig
            and _terminal_anchor_bp(aln, right_side, terminal_window_bp) >= min_anchor_bp
        ]
        pair = _best_bridge_pair(left_hits, right_hits)
        if pair is None:
            continue
        left, right = pair
        bridges.append(
            ReadBridge(
                read=read,
                left_orientation=left.orientation,
                right_orientation=right.orientation,
                read_order=_read_order(left, right),
                read_gap_bp=_read_gap(left, right),
            )
        )

    bridges.sort(key=lambda bridge: bridge.read)
    return ContigBridgeSupport(
        left_contig=left_contig,
        right_contig=right_contig,
        left_side=left_side,
        right_side=right_side,
        terminal_window_bp=terminal_window_bp,
        min_anchor_bp=min_anchor_bp,
        bridges=tuple(bridges),
    )


def read_depth_at(evidence: LongReadEvidence, contig, position):
    """Count unique reads with an alignment covering a 1-based contig position."""

    reads = {
        aln.read
        for aln in evidence.by_contig.get(contig, [])
        if aln.contig_start <= position <= aln.contig_end
    }
    return len(reads)


def read_depth_window(evidence: LongReadEvidence, contig, start, end):
    """Count unique reads with any alignment overlapping a 1-based closed window."""

    if end < start:
        raise ValueError("end must be greater than or equal to start")
    reads = {
        aln.read
        for aln in evidence.by_contig.get(contig, [])
        if _overlap_bp(aln.contig_start, aln.contig_end, start, end) > 0
    }
    return len(reads)
