"""GAF read-to-graph evidence helpers."""

from dataclasses import dataclass


@dataclass(frozen=True)
class GafPathRecord:
    """One read path from a GAF row."""

    query: str
    path: list
    mapq: int


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
