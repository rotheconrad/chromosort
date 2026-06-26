"""Assembly graph parsing helpers for graph-aware ChromoSort workflows."""

import gzip
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional


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


@dataclass(frozen=True)
class GraphPathStep:
    """One oriented segment step in a GFA path or walk."""

    segment: str
    orientation: str
    overlap_to_next: Optional[str] = None
    overlap_bp_to_next: Optional[int] = None


@dataclass(frozen=True)
class GraphPath:
    """One GFA path or walk record."""

    name: str
    steps: List[GraphPathStep]
    line_number: int
    tags: Dict[str, object] = field(default_factory=dict)
    record_type: str = "P"


@dataclass(frozen=True)
class GraphProjection:
    """Projected graph segment interval on a contig/path coordinate system."""

    path_name: str
    contig: str
    segment: str
    unitig: str
    segment_orientation: str
    step_index: int
    contig_start_0: int
    contig_end_0: int
    contig_start_1: int
    contig_end_1: int
    segment_length: int
    segment_start_0: int
    segment_end_0: int
    source_gfa: str
    has_sequence: bool
    overlap_left_bp: Optional[int]
    overlap_right_bp: Optional[int]
    duplicated_segment_count: int
    is_reused_segment: bool
    path_line_number: int = 0


@dataclass(frozen=True)
class GraphProjectionWarning:
    """Structured warning emitted while building graph projections."""

    severity: str
    code: str
    path_name: str
    segment: str
    line_number: int
    message: str


@dataclass(frozen=True)
class GraphPathSummary:
    """Summary of one projected path for report writing."""

    path_name: str
    fasta_contig: str
    step_count: int
    projected_bp: int
    fasta_length: Optional[int]
    length_diff_bp: Optional[int]
    length_status: str
    missing_segments: int
    zero_length_segments: int
    reused_segments: int
    source_gfa: str


@dataclass(frozen=True)
class GraphCoordinateEvidence:
    """Projection-aware graph context for one contig coordinate."""

    status: str
    contig: str
    position_1: Optional[int]
    containing_unitig: str = "."
    unitig_orientation: str = "."
    unitig_offset_0: Optional[int] = None
    distance_to_nearest_unitig_boundary: Optional[int] = None
    nearest_graph_junction: str = "."
    graph_in_degree: Optional[int] = None
    graph_out_degree: Optional[int] = None
    near_path_step_boundary: Optional[bool] = None
    path_step_boundary_distance: Optional[int] = None


@dataclass
class AssemblyGraph:
    """Parsed subset of an assembly graph needed by ChromoSort."""

    nodes: Dict[str, GraphNode]
    edges: List[GraphEdge]
    paths: List[GraphPath] = field(default_factory=list)
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
        paths_by_name = defaultdict(list)
        for path in self.paths:
            paths_by_name[path.name].append(path)
        self._paths_by_name = dict(paths_by_name)

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

    def paths_named(self, name):
        """Return all path/walk records with a given name."""

        return list(self._paths_by_name.get(name, []))


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


def _parse_path_step_token(token, line_number):
    if len(token) < 2 or token[-1] not in ORIENTATIONS:
        raise ValueError(
            f"Malformed GFA path step {token!r} at line {line_number}; "
            "expected segment name followed by + or -."
        )
    segment = token[:-1]
    orientation = token[-1]
    validate_orientation(orientation, "path_step_orientation")
    return segment, orientation


def _parse_path(fields, line_number, strict):
    if len(fields) < 4:
        raise ValueError(f"Malformed GFA P record at line {line_number}.")

    name = fields[1]
    step_tokens = [] if fields[2] in {"", "*"} else fields[2].split(",")
    raw_overlaps = [] if fields[3] in {"", "*"} else fields[3].split(",")
    expected_overlaps = max(0, len(step_tokens) - 1)
    if strict and len(raw_overlaps) not in {0, expected_overlaps}:
        raise ValueError(
            f"GFA P record {name!r} line {line_number} has {len(step_tokens)} "
            f"step(s) but {len(raw_overlaps)} overlap CIGAR(s)."
        )

    steps = []
    for idx, token in enumerate(step_tokens):
        segment, orientation = _parse_path_step_token(token, line_number)
        overlap = raw_overlaps[idx] if idx < len(raw_overlaps) else None
        overlap_bp = parse_gfa_overlap_bp(overlap) if overlap else None
        steps.append(
            GraphPathStep(
                segment=segment,
                orientation=orientation,
                overlap_to_next=overlap,
                overlap_bp_to_next=overlap_bp,
            )
        )

    return GraphPath(
        name=name,
        steps=steps,
        line_number=line_number,
        tags=parse_gfa_tags(fields[4:]),
        record_type="P",
    )


def _parse_int_or_text(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return value


def _parse_walk_steps(walk, line_number):
    if walk in {"", "*"}:
        return []

    steps = []
    position = 0
    for match in re.finditer(r"([<>])([^<>]+)", walk):
        if match.start() != position:
            raise ValueError(f"Malformed GFA W walk at line {line_number}: {walk!r}.")
        direction = match.group(1)
        segment = match.group(2)
        orientation = "+" if direction == ">" else "-"
        steps.append(GraphPathStep(segment=segment, orientation=orientation))
        position = match.end()

    if position != len(walk):
        raise ValueError(f"Malformed GFA W walk at line {line_number}: {walk!r}.")
    return steps


def _parse_walk(fields, line_number):
    if len(fields) < 7:
        raise ValueError(f"Malformed GFA W record at line {line_number}.")

    sample = fields[1]
    hap_index = fields[2]
    seq_id = fields[3]
    path_name = seq_id if seq_id not in {"", "*"} else sample
    tags = {
        "W_sample": sample,
        "W_hap_index": _parse_int_or_text(hap_index),
        "W_sequence": seq_id,
        "W_sequence_start": _parse_int_or_text(fields[4]),
        "W_sequence_end": _parse_int_or_text(fields[5]),
    }
    tags.update(parse_gfa_tags(fields[7:]))
    return GraphPath(
        name=path_name,
        steps=_parse_walk_steps(fields[6], line_number),
        line_number=line_number,
        tags=tags,
        record_type="W",
    )


def read_gfa(path, strict=True):
    """Read GFA S, L, P, and W records into an AssemblyGraph."""

    nodes = {}
    edges = []
    paths = []
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
            elif record_type == "P":
                paths.append(_parse_path(fields, line_number, strict))
            elif record_type == "W":
                paths.append(_parse_walk(fields, line_number))

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
        paths=paths,
        path=path,
        record_counts=dict(record_counts),
    )


def _warn(warnings, severity, code, path_name, segment, line_number, message):
    if warnings is not None:
        warnings.append(
            GraphProjectionWarning(
                severity=severity,
                code=code,
                path_name=path_name or ".",
                segment=segment or ".",
                line_number=line_number or 0,
                message=message,
            )
        )


def _selected_paths(graph, path_names, warnings):
    if path_names is None:
        selected = list(graph.paths)
    else:
        if not graph.paths:
            _warn(
                warnings,
                "warning",
                "no_paths",
                ".",
                ".",
                0,
                "No GFA P/W path records were found.",
            )
            return []
        selected = []
        for name in dict.fromkeys(path_names):
            matches = graph.paths_named(name)
            if matches:
                selected.extend(matches)
            else:
                _warn(
                    warnings,
                    "warning",
                    "missing_path",
                    name,
                    ".",
                    0,
                    f"No GFA P/W path named {name!r} was found.",
                )

    if not selected:
        code = "no_paths" if path_names is None else "no_matching_paths"
        message = (
            "No GFA P/W path records were found."
            if path_names is None
            else "No requested GFA P/W path records were found."
        )
        _warn(warnings, "warning", code, ".", ".", 0, message)
    return selected


def _trimmed_segment_interval(length, orientation, overlap_left_bp):
    if not overlap_left_bp:
        return 0, length

    trim = min(max(0, overlap_left_bp), length)
    if orientation == "+":
        return trim, length
    return 0, length - trim


def build_path_projection(graph, path_names=None, trim_overlaps=False, warnings=None):
    """
    Project GFA path/walk segment steps onto path/contig coordinates.

    Coordinates in the returned ``GraphProjection`` rows are zero-based,
    half-open for ``*_0`` fields and one-based inclusive for ``*_1`` fields.
    By default segment overlaps are reported but not subtracted from path spans.
    """

    selected_paths = _selected_paths(graph, path_names, warnings)
    segment_counts = Counter(
        step.segment
        for graph_path in selected_paths
        for step in graph_path.steps
    )
    source_gfa = str(graph.path) if graph.path else "."
    projections = []

    for graph_path in selected_paths:
        if not graph_path.steps:
            _warn(
                warnings,
                "warning",
                "empty_path",
                graph_path.name,
                ".",
                graph_path.line_number,
                f"GFA {graph_path.record_type} path {graph_path.name!r} has no segment steps.",
            )
            continue

        cursor = 0
        previous_overlap = None
        for idx, step in enumerate(graph_path.steps):
            node = graph.nodes.get(step.segment)
            if node is None:
                _warn(
                    warnings,
                    "warning",
                    "missing_segment",
                    graph_path.name,
                    step.segment,
                    graph_path.line_number,
                    (
                        f"GFA path {graph_path.name!r} references missing segment "
                        f"{step.segment!r}; skipping this step."
                    ),
                )
                previous_overlap = step.overlap_bp_to_next
                continue

            if node.length == 0 and not node.has_sequence and "LN" not in node.tags:
                _warn(
                    warnings,
                    "warning",
                    "missing_segment_length",
                    graph_path.name,
                    step.segment,
                    node.line_number,
                    (
                        f"Segment {step.segment} has no sequence and no LN:i length; "
                        "cannot project coordinates reliably."
                    ),
                )

            segment_start_0, segment_end_0 = (
                _trimmed_segment_interval(node.length, step.orientation, previous_overlap)
                if trim_overlaps
                else (0, node.length)
            )
            span = max(0, segment_end_0 - segment_start_0)
            contig_start_0 = cursor
            contig_end_0 = cursor + span
            duplicated_segment_count = segment_counts.get(step.segment, 0)
            projections.append(
                GraphProjection(
                    path_name=graph_path.name,
                    contig=graph_path.name,
                    segment=step.segment,
                    unitig=step.segment,
                    segment_orientation=step.orientation,
                    step_index=idx + 1,
                    contig_start_0=contig_start_0,
                    contig_end_0=contig_end_0,
                    contig_start_1=contig_start_0 + 1,
                    contig_end_1=contig_end_0,
                    segment_length=node.length,
                    segment_start_0=segment_start_0,
                    segment_end_0=segment_end_0,
                    source_gfa=source_gfa,
                    has_sequence=node.has_sequence,
                    overlap_left_bp=previous_overlap,
                    overlap_right_bp=step.overlap_bp_to_next,
                    duplicated_segment_count=duplicated_segment_count,
                    is_reused_segment=duplicated_segment_count > 1,
                    path_line_number=graph_path.line_number,
                )
            )
            cursor = contig_end_0
            previous_overlap = step.overlap_bp_to_next

    return projections


def build_direct_projection(graph, contig_names=None, warnings=None):
    """
    Build direct node-to-contig projections for contig-level GFAs.

    This is useful when graph segment names already match FASTA contig names and
    no path/walk projection is needed or available.
    """

    names = list(contig_names) if contig_names is not None else list(graph.nodes)
    source_gfa = str(graph.path) if graph.path else "."
    projections = []
    for idx, name in enumerate(names, start=1):
        node = graph.nodes.get(name)
        if node is None:
            continue
        if node.length == 0 and not node.has_sequence and "LN" not in node.tags:
            _warn(
                warnings,
                "warning",
                "missing_segment_length",
                name,
                name,
                node.line_number,
                (
                    f"Segment {name} has no sequence and no LN:i length; "
                    "cannot project coordinates reliably."
                ),
            )
        projections.append(
            GraphProjection(
                path_name=name,
                contig=name,
                segment=name,
                unitig=name,
                segment_orientation="+",
                step_index=idx,
                contig_start_0=0,
                contig_end_0=node.length,
                contig_start_1=1,
                contig_end_1=node.length,
                segment_length=node.length,
                segment_start_0=0,
                segment_end_0=node.length,
                source_gfa=source_gfa,
                has_sequence=node.has_sequence,
                overlap_left_bp=None,
                overlap_right_bp=None,
                duplicated_segment_count=1,
                is_reused_segment=False,
                path_line_number=node.line_number,
            )
        )
    if names and not projections:
        _warn(
            warnings,
            "warning",
            "no_direct_contig_matches",
            ".",
            ".",
            0,
            "No GFA segment names matched the requested FASTA contig names.",
        )
    return projections


def _length_status(projected_bp, fasta_length, tolerance_bp, tolerance_frac):
    if fasta_length is None:
        return "no_fasta_length"
    diff = abs(projected_bp - fasta_length)
    tolerance = max(tolerance_bp, int(round(fasta_length * tolerance_frac)))
    if diff <= tolerance:
        return "ok"
    return "length_mismatch"


def summarize_projections(
    projections,
    fasta_lengths=None,
    requested_path_names=None,
    tolerance_bp=1000,
    tolerance_frac=0.01,
    warnings=None,
):
    """Summarize projected path lengths and compare to optional FASTA lengths."""

    fasta_lengths = fasta_lengths or {}
    by_path = defaultdict(list)
    for projection in projections:
        by_path[projection.path_name].append(projection)
    missing_segment_counts = Counter(
        warning.path_name
        for warning in (warnings or [])
        if warning.code == "missing_segment"
    )

    summary_names = list(dict.fromkeys(requested_path_names or by_path.keys()))
    for name in by_path:
        if name not in summary_names:
            summary_names.append(name)

    summaries = []
    for name in summary_names:
        rows = by_path.get(name, [])
        if not rows:
            fasta_length = fasta_lengths.get(name)
            summaries.append(
                GraphPathSummary(
                    path_name=name,
                    fasta_contig=name if fasta_length is not None else ".",
                    step_count=0,
                    projected_bp=0,
                    fasta_length=fasta_length,
                    length_diff_bp=None if fasta_length is None else -fasta_length,
                    length_status="missing_projection",
                    missing_segments=1,
                    zero_length_segments=0,
                    reused_segments=0,
                    source_gfa=".",
                )
            )
            continue
        projected_bp = sum(row.contig_end_0 - row.contig_start_0 for row in rows)
        fasta_length = fasta_lengths.get(name)
        status = _length_status(projected_bp, fasta_length, tolerance_bp, tolerance_frac)
        diff = None if fasta_length is None else projected_bp - fasta_length
        if status == "length_mismatch":
            _warn(
                warnings,
                "warning",
                "path_length_mismatch",
                name,
                ".",
                rows[0].path_line_number if rows else 0,
                (
                    f"Projected path {name!r} length {projected_bp} bp differs from "
                    f"matching FASTA contig length {fasta_length} bp by {diff} bp."
                ),
            )

        summaries.append(
            GraphPathSummary(
                path_name=name,
                fasta_contig=name if fasta_length is not None else ".",
                step_count=len(rows),
                projected_bp=projected_bp,
                fasta_length=fasta_length,
                length_diff_bp=diff,
                length_status=status,
                missing_segments=missing_segment_counts.get(name, 0),
                zero_length_segments=sum(1 for row in rows if row.segment_length == 0),
                reused_segments=sum(1 for row in rows if row.is_reused_segment),
                source_gfa=rows[0].source_gfa if rows else ".",
            )
        )
    return summaries


def _tsv_value(value):
    if value is None:
        return "."
    if isinstance(value, bool):
        return "yes" if value else "no"
    return str(value)


def write_projection_tsv(path, projections):
    header = [
        "path_name",
        "contig",
        "segment",
        "unitig",
        "segment_orientation",
        "step_index",
        "contig_start_0",
        "contig_end_0",
        "contig_start_1",
        "contig_end_1",
        "segment_length",
        "segment_start_0",
        "segment_end_0",
        "source_gfa",
        "has_sequence",
        "overlap_left_bp",
        "overlap_right_bp",
        "duplicated_segment_count",
        "is_reused_segment",
        "path_line_number",
    ]
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for row in projections:
            out.write("\t".join(_tsv_value(getattr(row, col)) for col in header) + "\n")


def write_path_summary_tsv(path, summaries):
    header = [
        "path_name",
        "fasta_contig",
        "step_count",
        "projected_bp",
        "fasta_length",
        "length_diff_bp",
        "length_status",
        "missing_segments",
        "zero_length_segments",
        "reused_segments",
        "source_gfa",
    ]
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for row in summaries:
            out.write("\t".join(_tsv_value(getattr(row, col)) for col in header) + "\n")


def write_projection_warnings_tsv(path, warnings):
    header = ["severity", "code", "path_name", "segment", "line_number", "message"]
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for warning in warnings:
            out.write("\t".join(_tsv_value(getattr(warning, col)) for col in header) + "\n")


def projections_by_contig(projections):
    by_contig = defaultdict(list)
    for projection in projections:
        by_contig[projection.contig].append(projection)
    for rows in by_contig.values():
        rows.sort(key=lambda row: (row.contig_start_0, row.contig_end_0, row.step_index))
    return dict(by_contig)


def graph_coordinate_evidence(
    graph,
    projections,
    contig,
    position_1,
    boundary_window_bp=5000,
):
    """Summarize graph context for a contig coordinate using projections."""

    if graph is None or not projections or position_1 is None:
        return GraphCoordinateEvidence("no_projection", contig, position_1)

    position_0 = max(0, int(position_1) - 1)
    candidates = [
        row
        for row in projections_by_contig(projections).get(contig, [])
        if row.contig_start_0 <= position_0 < row.contig_end_0
    ]
    if not candidates:
        return GraphCoordinateEvidence("outside_projection", contig, position_1)

    projection = candidates[0]
    if projection.segment_orientation == "+":
        unitig_offset = projection.segment_start_0 + (position_0 - projection.contig_start_0)
    else:
        unitig_offset = projection.segment_end_0 - (position_0 - projection.contig_start_0) - 1

    left_distance = position_0 - projection.contig_start_0
    right_distance = projection.contig_end_0 - position_0
    boundary_distance = min(left_distance, right_distance)
    node_evidence = graph_node_evidence(graph, [projection.segment])
    is_junction = (
        node_evidence.graph_self_loop
        or (node_evidence.graph_in_degree or 0) > 1
        or (node_evidence.graph_out_degree or 0) > 1
        or (node_evidence.graph_neighbor_count or 0) > 2
    )
    return GraphCoordinateEvidence(
        status="present",
        contig=contig,
        position_1=position_1,
        containing_unitig=projection.segment,
        unitig_orientation=projection.segment_orientation,
        unitig_offset_0=unitig_offset,
        distance_to_nearest_unitig_boundary=boundary_distance,
        nearest_graph_junction=projection.segment if is_junction else ".",
        graph_in_degree=node_evidence.graph_in_degree,
        graph_out_degree=node_evidence.graph_out_degree,
        near_path_step_boundary=boundary_distance <= boundary_window_bp,
        path_step_boundary_distance=boundary_distance,
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
