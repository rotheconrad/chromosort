"""GAF read-to-graph evidence helpers."""

from collections import Counter
from dataclasses import dataclass


@dataclass(frozen=True)
class GafPathRecord:
    """One read path from a GAF row."""

    query: str
    path: list
    mapq: int


@dataclass(frozen=True)
class GafPathSupport:
    """Read support for one oriented graph path."""

    path_nodes: str
    support: int
    reads: tuple

    @property
    def read_count(self):
        return len(self.reads)


@dataclass(frozen=True)
class GafTraversalSummary:
    """Support summary comparing a selected path against alternatives."""

    selected_path_nodes: str = "."
    selected_support: int = 0
    best_alt_path_nodes: str = "."
    best_alt_support: int = 0
    best_path_nodes: str = "."
    best_support: int = 0
    support_status: str = "no_paths"
    selected_reads: tuple = ()
    path_supports: tuple = ()

    @property
    def selected_read_count(self):
        return len(self.selected_reads)


@dataclass(frozen=True)
class GafNodeSupport:
    """GAF read-path evidence touching one graph node."""

    node: str
    traversal_count: int
    reads: tuple
    orientation_summary: str = "."
    path_examples: tuple = ()

    @property
    def read_count(self):
        return len(self.reads)


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


def oriented_nodes_text(nodes):
    if not nodes:
        return "."
    return ",".join(f"{name}{orientation}" for name, orientation in nodes)


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


def gaf_records_for_oriented_nodes(gaf_records, nodes):
    candidate = list(nodes)
    reverse_candidate = reverse_oriented_nodes(candidate)
    supporting = []
    for record in gaf_records:
        if contains_oriented_subpath(record.path, candidate) or contains_oriented_subpath(
            record.path,
            reverse_candidate,
        ):
            supporting.append(record)
    return supporting


def gaf_support_for_oriented_nodes(gaf_records, nodes):
    return len(gaf_records_for_oriented_nodes(gaf_records, nodes))


def gaf_support_for_path(gaf_records, path):
    return gaf_support_for_oriented_nodes(gaf_records, path_oriented_nodes(path))


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


def gaf_path_support_details(gaf_records, paths):
    details = []
    for path in paths:
        nodes = path_oriented_nodes(path)
        records = gaf_records_for_oriented_nodes(gaf_records, nodes)
        reads = tuple(sorted({record.query for record in records}))
        details.append(
            GafPathSupport(
                path_nodes=oriented_nodes_text(nodes),
                support=len(records),
                reads=reads,
            )
        )
    return tuple(details)


def summarize_gaf_traversal(gaf_records, paths, selected_index=0, min_support=1):
    paths = list(paths)
    if not paths:
        return GafTraversalSummary()
    if selected_index < 0 or selected_index >= len(paths):
        raise ValueError("selected_index is outside the path list")

    path_supports = gaf_path_support_details(gaf_records or [], paths)
    selected = path_supports[selected_index]
    alternatives = [
        support for index, support in enumerate(path_supports)
        if index != selected_index
    ]
    best = max(path_supports, key=lambda item: item.support)
    best_alternatives = [
        support for support in alternatives
        if support.support == max((item.support for item in alternatives), default=0)
    ]
    best_alt = best_alternatives[0] if best_alternatives else None

    support_values = [support.support for support in path_supports]
    if not any(support_values):
        status = "no_support"
    elif selected.support < min_support and best.support < min_support:
        status = "weak_support"
    elif support_values.count(best.support) > 1:
        status = "tied_support"
    elif best is selected:
        status = "supports_selected"
    else:
        status = "supports_alternate"

    return GafTraversalSummary(
        selected_path_nodes=selected.path_nodes,
        selected_support=selected.support,
        best_alt_path_nodes=best_alt.path_nodes if best_alt else ".",
        best_alt_support=best_alt.support if best_alt else 0,
        best_path_nodes=best.path_nodes,
        best_support=best.support,
        support_status=status,
        selected_reads=selected.reads,
        path_supports=path_supports,
    )


def summarize_gaf_node(gaf_records, node, max_examples=3):
    if not node or node == ".":
        return GafNodeSupport(node=".", traversal_count=0, reads=())

    records = []
    orientations = Counter()
    examples = []
    for record in gaf_records or []:
        matched = False
        for path_node, orientation in record.path:
            if path_node != node:
                continue
            orientations[orientation] += 1
            matched = True
        if not matched:
            continue
        records.append(record)
        example = oriented_nodes_text(record.path)
        if len(examples) < max_examples and example not in examples:
            examples.append(example)

    orientation_summary = (
        ",".join(f"{orientation}:{orientations[orientation]}" for orientation in sorted(orientations))
        if orientations
        else "."
    )
    return GafNodeSupport(
        node=node,
        traversal_count=len(records),
        reads=tuple(sorted({record.query for record in records})),
        orientation_summary=orientation_summary,
        path_examples=tuple(examples),
    )
