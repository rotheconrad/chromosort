"""Assembly graph parsing helpers for graph-aware ChromoSort workflows."""

import gzip
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


ORIENTATIONS = {"+", "-"}
OVERLAP_RE = re.compile(r"(\d+)([MIDNSHP=X])")
OVERLAP_CONSUMES_BOTH = {"M", "=", "X"}
COMPLEX_OVERLAP_OPS = {"I", "D", "N", "S", "H", "P"}


@dataclass(frozen=True)
class GraphNode:
    """One GFA segment record."""

    name: str
    length: int
    sequence: Optional[str] = None
    tags: Dict[str, object] = field(default_factory=dict)
    line_number: int = 0

    @property
    def has_sequence(self):
        return self.sequence is not None


@dataclass(frozen=True)
class GraphEdge:
    """One oriented GFA link record."""

    source: str
    source_orientation: str
    target: str
    target_orientation: str
    overlap: str
    overlap_bp: Optional[int]
    tags: Dict[str, object] = field(default_factory=dict)
    line_number: int = 0

    @property
    def source_key(self):
        return self.source, self.source_orientation

    @property
    def target_key(self):
        return self.target, self.target_orientation


@dataclass
class AssemblyGraph:
    """Parsed subset of an assembly graph needed by ChromoSort."""

    nodes: Dict[str, GraphNode]
    edges: List[GraphEdge]
    path: Optional[Path] = None
    record_counts: Dict[str, int] = field(default_factory=dict)

    def __post_init__(self):
        outgoing = defaultdict(list)
        incoming = defaultdict(list)
        for edge in self.edges:
            outgoing[edge.source_key].append(edge)
            incoming[edge.target_key].append(edge)
        self._outgoing = dict(outgoing)
        self._incoming = dict(incoming)

    def outgoing(self, node, orientation=None):
        """Return outgoing edges for a node, optionally in one orientation."""

        if orientation is not None:
            validate_orientation(orientation, "orientation")
            return list(self._outgoing.get((node, orientation), []))

        result = []
        for orient in sorted(ORIENTATIONS):
            result.extend(self._outgoing.get((node, orient), []))
        return result

    def incoming(self, node, orientation=None):
        """Return incoming edges for a node, optionally in one orientation."""

        if orientation is not None:
            validate_orientation(orientation, "orientation")
            return list(self._incoming.get((node, orientation), []))

        result = []
        for orient in sorted(ORIENTATIONS):
            result.extend(self._incoming.get((node, orient), []))
        return result

    def direct_edges(
        self,
        source,
        target,
        source_orientation=None,
        target_orientation=None,
    ):
        """Return oriented links from source to target matching optional orientations."""

        if source_orientation is not None:
            validate_orientation(source_orientation, "source_orientation")
        if target_orientation is not None:
            validate_orientation(target_orientation, "target_orientation")

        result = []
        candidate_orientations = (
            [source_orientation] if source_orientation is not None else sorted(ORIENTATIONS)
        )
        for orient in candidate_orientations:
            for edge in self._outgoing.get((source, orient), []):
                if edge.target != target:
                    continue
                if target_orientation is not None and edge.target_orientation != target_orientation:
                    continue
                result.append(edge)
        return result

    def has_direct_edge(
        self,
        source,
        target,
        source_orientation=None,
        target_orientation=None,
    ):
        return bool(
            self.direct_edges(
                source,
                target,
                source_orientation=source_orientation,
                target_orientation=target_orientation,
            )
        )


@dataclass(frozen=True)
class GraphNodeEvidence:
    graph_node: str
    graph_node_status: str
    graph_node_length: Optional[int]
    graph_node_has_sequence: Optional[bool]
    graph_in_degree: Optional[int]
    graph_out_degree: Optional[int]
    graph_neighbor_count: Optional[int]
    graph_self_loop: Optional[bool]


@dataclass(frozen=True)
class GraphLinkEvidence:
    graph_link_status: str
    graph_direct_edge: bool
    graph_direct_edge_orientations: str
    graph_direct_edge_overlap_bp: str


def open_text(path):
    path = Path(path)
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def validate_orientation(value, field):
    if value not in ORIENTATIONS:
        raise ValueError(f"Invalid GFA orientation for {field}: {value!r}.")


def parse_gfa_tags(fields: Iterable[str]):
    tags = {}
    for raw in fields:
        parts = raw.split(":", 2)
        if len(parts) != 3:
            continue
        tag, tag_type, value = parts
        if tag_type == "i":
            try:
                tags[tag] = int(value)
            except ValueError:
                tags[tag] = value
        elif tag_type == "f":
            try:
                tags[tag] = float(value)
            except ValueError:
                tags[tag] = value
        elif tag_type == "Z":
            tags[tag] = value
        else:
            tags[tag] = value
    return tags


def parse_gfa_overlap_bp(overlap):
    """
    Parse simple GFA link overlap CIGAR strings.

    Returns the number of overlap bases when the CIGAR is composed only of
    operations that consume both linked segments. More complex CIGARs are kept
    as raw strings and return None so later sequence-changing code cannot
    accidentally treat them as exact trim lengths.
    """

    if overlap == "*":
        return None

    total = 0
    position = 0
    saw_component = False
    for match in OVERLAP_RE.finditer(overlap):
        if match.start() != position:
            raise ValueError(f"Malformed GFA overlap CIGAR: {overlap!r}.")
        saw_component = True
        length = int(match.group(1))
        op = match.group(2)
        if op in OVERLAP_CONSUMES_BOTH:
            total += length
        elif op in COMPLEX_OVERLAP_OPS:
            return None
        else:
            return None
        position = match.end()

    if not saw_component or position != len(overlap):
        raise ValueError(f"Malformed GFA overlap CIGAR: {overlap!r}.")
    return total


def _parse_segment(fields, line_number, strict):
    if len(fields) < 3:
        raise ValueError(f"Malformed GFA S record at line {line_number}.")

    name = fields[1]
    sequence = None if fields[2] == "*" else fields[2]
    tags = parse_gfa_tags(fields[3:])
    tagged_length = tags.get("LN")
    if sequence is None:
        length = tagged_length if isinstance(tagged_length, int) else 0
    else:
        length = len(sequence)
        if strict and isinstance(tagged_length, int) and tagged_length != length:
            raise ValueError(
                f"GFA segment {name!r} line {line_number} has sequence length "
                f"{length} but LN:i:{tagged_length}."
            )

    return GraphNode(
        name=name,
        length=length,
        sequence=sequence,
        tags=tags,
        line_number=line_number,
    )


def _parse_link(fields, line_number):
    if len(fields) < 6:
        raise ValueError(f"Malformed GFA L record at line {line_number}.")

    source = fields[1]
    source_orientation = fields[2]
    target = fields[3]
    target_orientation = fields[4]
    validate_orientation(source_orientation, "source_orientation")
    validate_orientation(target_orientation, "target_orientation")
    overlap = fields[5]
    overlap_bp = parse_gfa_overlap_bp(overlap)
    tags = parse_gfa_tags(fields[6:])
    return GraphEdge(
        source=source,
        source_orientation=source_orientation,
        target=target,
        target_orientation=target_orientation,
        overlap=overlap,
        overlap_bp=overlap_bp,
        tags=tags,
        line_number=line_number,
    )


def read_gfa(path, strict=True):
    """Read GFA S and L records into an AssemblyGraph."""

    nodes = {}
    edges = []
    record_counts = defaultdict(int)
    path = Path(path)

    with open_text(path) as fh:
        for line_number, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            record_type = fields[0]
            record_counts[record_type] += 1

            if record_type == "S":
                node = _parse_segment(fields, line_number, strict)
                if node.name in nodes:
                    raise ValueError(f"Duplicate GFA segment {node.name!r} at line {line_number}.")
                nodes[node.name] = node
            elif record_type == "L":
                edges.append(_parse_link(fields, line_number))

    if strict:
        missing = []
        for edge in edges:
            if edge.source not in nodes:
                missing.append((edge.source, edge.line_number))
            if edge.target not in nodes:
                missing.append((edge.target, edge.line_number))
        if missing:
            preview = ", ".join(f"{name!r} on line {line}" for name, line in missing[:5])
            suffix = "..." if len(missing) > 5 else ""
            raise ValueError(f"GFA link references missing segment(s): {preview}{suffix}.")

    return AssemblyGraph(
        nodes=nodes,
        edges=edges,
        path=path,
        record_counts=dict(record_counts),
    )


def oriented_node(name, orientation):
    validate_orientation(orientation, "orientation")
    return name, orientation


def format_oriented_node(item):
    name, orientation = item
    return f"{name}{orientation}"


def resolve_graph_node(graph, candidates):
    """Resolve a sequence/report name to the first matching GFA segment."""

    for candidate in candidates:
        if candidate and candidate != "." and candidate in graph.nodes:
            return candidate
    return "."


def _unique_neighbor_count(edges, node):
    neighbors = set()
    for edge in edges:
        if edge.source != node:
            neighbors.add(edge.source)
        if edge.target != node:
            neighbors.add(edge.target)
    return len(neighbors)


def graph_node_evidence(graph, candidates):
    """Summarize whether and how a report row maps to a GFA segment."""

    graph_node = resolve_graph_node(graph, candidates)
    if graph_node == ".":
        return GraphNodeEvidence(
            graph_node=".",
            graph_node_status="missing_node",
            graph_node_length=None,
            graph_node_has_sequence=None,
            graph_in_degree=None,
            graph_out_degree=None,
            graph_neighbor_count=None,
            graph_self_loop=None,
        )

    node = graph.nodes[graph_node]
    incoming = graph.incoming(graph_node)
    outgoing = graph.outgoing(graph_node)
    self_loop = any(edge.source == edge.target == graph_node for edge in incoming + outgoing)
    return GraphNodeEvidence(
        graph_node=graph_node,
        graph_node_status="present",
        graph_node_length=node.length,
        graph_node_has_sequence=node.has_sequence,
        graph_in_degree=len(incoming),
        graph_out_degree=len(outgoing),
        graph_neighbor_count=_unique_neighbor_count(incoming + outgoing, graph_node),
        graph_self_loop=self_loop,
    )


def _edge_orientations(edges):
    if not edges:
        return "."
    return ",".join(
        f"{edge.source}{edge.source_orientation}>{edge.target}{edge.target_orientation}"
        for edge in edges
    )


def _edge_overlap_bps(edges):
    values = []
    for edge in edges:
        value = "." if edge.overlap_bp is None else str(edge.overlap_bp)
        if value not in values:
            values.append(value)
    return ",".join(values) if values else "."


def graph_link_evidence(graph, left_node, right_node):
    """Summarize direct GFA links between two resolved graph nodes."""

    if left_node == "." and right_node == ".":
        status = "missing_both_nodes"
        direct_edges = []
    elif left_node == ".":
        status = "missing_left_node"
        direct_edges = []
    elif right_node == ".":
        status = "missing_right_node"
        direct_edges = []
    else:
        forward = graph.direct_edges(left_node, right_node)
        reverse = graph.direct_edges(right_node, left_node)
        direct_edges = forward + reverse
        if forward and reverse:
            status = "bidirectional_direct_edge"
        elif forward:
            status = "direct_edge"
        elif reverse:
            status = "reverse_direct_edge"
        else:
            status = "no_direct_edge"

    return GraphLinkEvidence(
        graph_link_status=status,
        graph_direct_edge=bool(direct_edges),
        graph_direct_edge_orientations=_edge_orientations(direct_edges),
        graph_direct_edge_overlap_bp=_edge_overlap_bps(direct_edges),
    )
