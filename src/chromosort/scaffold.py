#!/usr/bin/env python3
"""
Build chromosome-level scaffold FASTA records from ChromoSort-ordered contigs.

The default gap model infers N runs from adjacent reference coordinates in the
chromo sort contig assignment report. Users can instead provide a fixed number
of Ns between neighboring contigs.
"""

import argparse
import csv
import sys
from collections import OrderedDict, deque
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

from .graph import ORIENTATIONS, format_oriented_node, read_gfa
from .paths import ensure_output_dirs, ensure_parent_dir
from .reference_order import iter_fasta_records, write_wrapped


REQUIRED_ASSIGNMENT_COLUMNS = {
    "contig",
    "kept",
    "new_name",
    "assigned_ref",
    "order_in_ref",
    "ref_start",
    "ref_end",
}


@dataclass
class AssignmentRow:
    contig: str
    new_name: str
    ref: str
    order_in_ref: int
    ref_start: int
    ref_end: int
    orientation: str = "."


@dataclass
class FastaRecord:
    name: str
    header: str
    seq: str


@dataclass
class ScaffoldMember:
    assignment: AssignmentRow
    record: FastaRecord
    trim_left_bp: int = 0
    trim_right_bp: int = 0

    @property
    def seq(self):
        end = len(self.record.seq) - self.trim_right_bp
        return self.record.seq[self.trim_left_bp : end]

    @property
    def trimmed_bp(self):
        return self.trim_left_bp + self.trim_right_bp


@dataclass
class GapRecord:
    scaffold: str
    left_contig: str
    right_contig: str
    left_ref_end: int
    right_ref_start: int
    raw_inferred_gap_bp: int
    gap_bp: int
    gap_mode: str
    overlap_bp: int
    overlap_class: str
    overlap_left_ref_frac: float
    overlap_right_ref_frac: float
    overlap_policy: str
    graph_overlap_policy: str
    overlap_action: str
    graph_overlap_action: str
    trimmed_bp: int
    sequence_overlap_identity: Optional[float] = None


@dataclass
class GraphGapRecord:
    scaffold: str
    left_contig: str
    right_contig: str
    left_graph_node: str
    right_graph_node: str
    left_orientation: str
    right_orientation: str
    graph_status: str
    direct_edge: bool
    direct_edge_orientations: str
    direct_edge_overlap_bp: str
    shortest_path_edges: Optional[int]
    shortest_path_nodes: str
    candidate_intermediate_nodes: str
    max_path_edges: int


@dataclass
class ScaffoldRecord:
    name: str
    seq: str
    members: list
    gaps: list

    @property
    def sequence_bp(self):
        return sum(len(member.seq) for member in self.members)

    @property
    def gap_bp(self):
        return sum(gap.gap_bp for gap in self.gaps)

    @property
    def overlap_bp(self):
        return sum(gap.overlap_bp for gap in self.gaps)

    @property
    def trimmed_bp(self):
        return sum(member.trimmed_bp for member in self.members)


def parse_args(argv: Optional[Sequence[str]] = None, prog: Optional[str] = None):
    ap = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Scaffold a ChromoSort ordered FASTA into per-reference records "
            "using inferred or fixed N gaps."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument(
        "-f",
        "--ordered-fasta",
        required=True,
        help="Final ordered FASTA from chromo sort, optionally after chromo fix and re-sorting.",
    )
    ap.add_argument(
        "-a",
        "--assignments",
        required=True,
        help="Corresponding <prefix>.contig_assignments.tsv from chromo sort.",
    )
    ap.add_argument(
        "-o",
        "--output-prefix",
        required=True,
        help=(
            "Output prefix. Writes <prefix>.scaffold.fa, <prefix>.scaffold_gaps.tsv, "
            "<prefix>.scaffold_summary.tsv, and <prefix>.run_summary.txt."
        ),
    )
    ap.add_argument(
        "--fixed-gap-bp",
        "--fixed-gap",
        dest="fixed_gap_bp",
        type=int,
        default=None,
        help=(
            "Use this many Ns between neighboring contigs instead of inferring "
            "gap length from reference coordinates."
        ),
    )
    ap.add_argument(
        "--overlap-policy",
        choices=["zero-gap", "warn", "trim-reference", "trim-sequence"],
        default="zero-gap",
        help=(
            "How to handle adjacent contigs with overlapping reference spans. "
            "zero-gap writes no Ns and keeps both sequences unchanged; warn does "
            "the same but emits stderr warnings; trim-reference removes the "
            "reference-inferred overlap from the right contig; trim-sequence "
            "removes the right-contig prefix only when the terminal sequence "
            "overlap has high identity."
        ),
    )
    ap.add_argument(
        "--trim-sequence-min-identity",
        type=float,
        default=0.98,
        help=(
            "Minimum suffix/prefix identity required by "
            "--overlap-policy trim-sequence."
        ),
    )
    ap.add_argument(
        "--simple-headers",
        action="store_true",
        help="Write FASTA headers containing only the scaffold sequence ID.",
    )
    ap.add_argument(
        "--gfa",
        default=None,
        help=(
            "Optional assembly graph GFA. When provided, writes "
            "<prefix>.graph_gaps.tsv with report-only graph evidence for "
            "each adjacent scaffold junction."
        ),
    )
    ap.add_argument(
        "--graph-overlap-policy",
        choices=["report", "warn", "confirm"],
        default="report",
        help=(
            "How graph evidence should affect negative reference-gap overlaps. "
            "report keeps graph evidence report-only; warn emits warnings for "
            "graph-confirmed overlaps; confirm allows a direct oriented GFA "
            "edge to trim only terminal overlaps that would otherwise use the "
            "zero-gap/warn overlap policy."
        ),
    )
    ap.add_argument(
        "--graph-max-path-edges",
        type=int,
        default=4,
        help=(
            "Maximum explicit GFA link depth searched when reporting short "
            "paths between adjacent scaffold contigs."
        ),
    )
    return ap.parse_args(argv)


def parse_int(value, field, row_name):
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Expected integer {field} for {row_name!r}, found {value!r}") from exc


def read_assignments(path):
    assignments = OrderedDict()
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = REQUIRED_ASSIGNMENT_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Assignment report {path} is missing required columns: {', '.join(sorted(missing))}"
            )
        for row in reader:
            kept = row["kept"].strip().lower() in {"yes", "true", "1"}
            if not kept:
                continue
            new_name = row["new_name"].strip()
            if not new_name or new_name == ".":
                continue
            if new_name in assignments:
                raise ValueError(f"Duplicate kept new_name {new_name!r} in {path}")
            assignments[new_name] = AssignmentRow(
                contig=row["contig"].strip(),
                new_name=new_name,
                ref=row["assigned_ref"].strip(),
                order_in_ref=parse_int(row["order_in_ref"], "order_in_ref", new_name),
                ref_start=parse_int(row["ref_start"], "ref_start", new_name),
                ref_end=parse_int(row["ref_end"], "ref_end", new_name),
                orientation=(row.get("orientation") or ".").strip() or ".",
            )
    return assignments


def read_ordered_fasta(path):
    records = []
    seen = set()
    for name, header, seq in iter_fasta_records(path):
        if name in seen:
            raise ValueError(f"Duplicate FASTA record {name!r} in {path}")
        seen.add(name)
        records.append(FastaRecord(name=name, header=header, seq=seq))
    return records


def group_scaffold_members(records, assignments):
    groups = OrderedDict()
    unassigned = []
    seen_assigned = set()

    for record in records:
        assignment = assignments.get(record.name)
        if assignment is None:
            unassigned.append(record)
            continue
        groups.setdefault(assignment.ref, []).append(
            ScaffoldMember(assignment=assignment, record=record)
        )
        seen_assigned.add(record.name)

    missing = [name for name in assignments if name not in seen_assigned]
    if missing:
        preview = ", ".join(missing[:5])
        suffix = "..." if len(missing) > 5 else ""
        raise ValueError(
            "Ordered FASTA is missing kept assignment sequence(s): "
            f"{preview}{suffix}. Use the ordered FASTA from the same chromo sort run."
        )

    return groups, unassigned


def inferred_gap(left, right):
    return right.assignment.ref_start - left.assignment.ref_end - 1


def ref_span_bp(assignment):
    return max(0, assignment.ref_end - assignment.ref_start + 1)


def fmt(value, digits=4):
    if value is None:
        return "."
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def classify_adjacent_overlap(left, right, overlap_bp):
    if overlap_bp <= 0:
        return "no_overlap"
    left_start = left.assignment.ref_start
    left_end = left.assignment.ref_end
    right_start = right.assignment.ref_start
    right_end = right.assignment.ref_end
    if right_start > left_start and right_end > left_end:
        return "terminal_overlap"
    if right_start >= left_start and right_end <= left_end:
        return "contained_overlap"
    if right_start <= left_start and right_end >= left_end:
        return "spanning_overlap"
    return "internal_overlap"


def sequence_overlap_identity(left_seq, right_seq, overlap_bp):
    size = min(overlap_bp, len(left_seq), len(right_seq))
    if size <= 0:
        return 0, None
    left_part = left_seq[-size:].upper()
    right_part = right_seq[:size].upper()
    matches = sum(1 for left_base, right_base in zip(left_part, right_part) if left_base == right_base)
    return size, matches / size


def graph_direct_edge_for_members(left, right, graph):
    if graph is None:
        return False, "."
    left_node = graph_node_name(left, graph)
    right_node = graph_node_name(right, graph)
    if missing_graph_status(left_node, right_node) is not None:
        return False, "."
    left_orientation = graph_orientation(left)
    right_orientation = graph_orientation(right)
    direct = graph.direct_edges(
        left_node,
        right_node,
        source_orientation=left_orientation if left_orientation in ORIENTATIONS else None,
        target_orientation=right_orientation if right_orientation in ORIENTATIONS else None,
    )
    return bool(direct), edge_orientations(direct)


def apply_overlap_policy(left, right, overlap_bp, overlap_class, args, graph=None):
    graph_overlap_action = "."
    if overlap_bp <= 0:
        return "none", graph_overlap_action, 0, None

    graph_direct, graph_edges = graph_direct_edge_for_members(left, right, graph)
    if graph_direct:
        graph_overlap_action = f"direct_edge:{graph_edges}"
        if (
            args.graph_overlap_policy == "confirm"
            and args.overlap_policy in {"zero-gap", "warn"}
            and overlap_class == "terminal_overlap"
        ):
            trim_bp = min(overlap_bp, len(right.seq))
            right.trim_left_bp += trim_bp
            return "graph_confirmed_trim_reference", graph_overlap_action, trim_bp, None

    if args.overlap_policy in {"zero-gap", "warn"}:
        return "zero_gap", graph_overlap_action, 0, None
    if overlap_class != "terminal_overlap":
        return "trim_skipped_nonterminal", graph_overlap_action, 0, None

    if args.overlap_policy == "trim-reference":
        trim_bp = min(overlap_bp, len(right.seq))
        right.trim_left_bp += trim_bp
        return "trimmed_reference", graph_overlap_action, trim_bp, None

    trim_bp, identity = sequence_overlap_identity(left.seq, right.seq, overlap_bp)
    if trim_bp > 0 and identity is not None and identity >= args.trim_sequence_min_identity:
        right.trim_left_bp += trim_bp
        return "trimmed_sequence", graph_overlap_action, trim_bp, identity
    return "trim_skipped_sequence_identity", graph_overlap_action, 0, identity


def build_scaffold(ref, members, fixed_gap_bp, args, graph=None):
    pieces = []
    gaps = []
    gap_mode = "fixed" if fixed_gap_bp is not None else "inferred"

    for index, member in enumerate(members):
        if index:
            left = members[index - 1]
            raw_gap = inferred_gap(left, member)
            overlap_bp = max(0, -raw_gap)
            overlap_class = classify_adjacent_overlap(left, member, overlap_bp)
            overlap_action, graph_overlap_action, trimmed_bp, sequence_identity = apply_overlap_policy(
                left,
                member,
                overlap_bp,
                overlap_class,
                args,
                graph,
            )
            gap_bp = fixed_gap_bp if fixed_gap_bp is not None else max(0, raw_gap)
            gaps.append(
                GapRecord(
                    scaffold=ref,
                    left_contig=left.assignment.new_name,
                    right_contig=member.assignment.new_name,
                    left_ref_end=left.assignment.ref_end,
                    right_ref_start=member.assignment.ref_start,
                    raw_inferred_gap_bp=raw_gap,
                    gap_bp=gap_bp,
                    gap_mode=gap_mode,
                    overlap_bp=overlap_bp,
                    overlap_class=overlap_class,
                    overlap_left_ref_frac=(
                        overlap_bp / ref_span_bp(left.assignment)
                        if ref_span_bp(left.assignment)
                        else 0.0
                    ),
                    overlap_right_ref_frac=(
                        overlap_bp / ref_span_bp(member.assignment)
                        if ref_span_bp(member.assignment)
                        else 0.0
                    ),
                    overlap_policy=args.overlap_policy,
                    graph_overlap_policy=args.graph_overlap_policy,
                    overlap_action=overlap_action,
                    graph_overlap_action=graph_overlap_action,
                    trimmed_bp=trimmed_bp,
                    sequence_overlap_identity=sequence_identity,
                )
            )
            if overlap_bp > 0 and (
                args.overlap_policy == "warn" or overlap_action.startswith("trim")
                or (
                    args.graph_overlap_policy in {"warn", "confirm"}
                    and graph_overlap_action != "."
                )
            ):
                sys.stderr.write(
                    "WARNING: "
                    f"{ref} overlap between {left.assignment.new_name} and "
                    f"{member.assignment.new_name}: raw_gap={raw_gap}, "
                    f"overlap_bp={overlap_bp}, class={overlap_class}, "
                    f"action={overlap_action}, graph_action={graph_overlap_action}, "
                    f"trimmed_bp={trimmed_bp}\n"
                )
            pieces.append("N" * gap_bp)
        pieces.append(member.seq)

    return ScaffoldRecord(
        name=ref,
        seq="".join(pieces),
        members=members,
        gaps=gaps,
    )


def build_scaffolds(groups, fixed_gap_bp, args, graph=None):
    return [
        build_scaffold(ref, members, fixed_gap_bp, args, graph)
        for ref, members in groups.items()
    ]


def graph_node_name(member, graph):
    candidates = [member.assignment.contig, member.assignment.new_name, member.record.name]
    for candidate in candidates:
        if candidate in graph.nodes:
            return candidate
    return "."


def graph_orientation(member):
    orientation = member.assignment.orientation
    return orientation if orientation in ORIENTATIONS else "."


def graph_start_keys(node, orientation):
    if node == ".":
        return []
    if orientation in ORIENTATIONS:
        return [(node, orientation)]
    return [(node, orient) for orient in sorted(ORIENTATIONS)]


def graph_target_keys(node, orientation):
    if node == ".":
        return set()
    if orientation in ORIENTATIONS:
        return {(node, orientation)}
    return {(node, orient) for orient in ORIENTATIONS}


def find_shortest_graph_path(graph, left_node, right_node, left_orientation, right_orientation, max_edges):
    if max_edges < 1:
        return []

    targets = graph_target_keys(right_node, right_orientation)
    if not targets:
        return []

    queue = deque((start, []) for start in graph_start_keys(left_node, left_orientation))
    seen = {start for start, _ in queue}
    while queue:
        (node, orientation), path = queue.popleft()
        if len(path) >= max_edges:
            continue
        for edge in graph.outgoing(node, orientation):
            next_key = edge.target_key
            next_path = path + [edge]
            if next_key in targets:
                return next_path
            if next_key in seen:
                continue
            seen.add(next_key)
            queue.append((next_key, next_path))
    return []


def graph_path_nodes(path):
    if not path:
        return "."
    nodes = [format_oriented_node(path[0].source_key)]
    nodes.extend(format_oriented_node(edge.target_key) for edge in path)
    return ",".join(nodes)


def graph_intermediate_nodes(path):
    if len(path) <= 1:
        return "."
    intermediates = []
    for edge in path[:-1]:
        node = format_oriented_node(edge.target_key)
        if node not in intermediates:
            intermediates.append(node)
    return ",".join(intermediates) if intermediates else "."


def edge_orientations(edges):
    if not edges:
        return "."
    return ",".join(
        f"{edge.source}{edge.source_orientation}>{edge.target}{edge.target_orientation}"
        for edge in edges
    )


def edge_overlap_bps(edges):
    values = []
    for edge in edges:
        value = fmt(edge.overlap_bp)
        if value not in values:
            values.append(value)
    return ",".join(values) if values else "."


def missing_graph_status(left_node, right_node):
    if left_node == "." and right_node == ".":
        return "missing_both_nodes"
    if left_node == ".":
        return "missing_left_node"
    if right_node == ".":
        return "missing_right_node"
    return None


def graph_gap_record(scaffold, left, right, graph, max_path_edges):
    left_node = graph_node_name(left, graph)
    right_node = graph_node_name(right, graph)
    left_orientation = graph_orientation(left)
    right_orientation = graph_orientation(right)

    status = missing_graph_status(left_node, right_node)
    if status is not None:
        return GraphGapRecord(
            scaffold=scaffold.name,
            left_contig=left.assignment.new_name,
            right_contig=right.assignment.new_name,
            left_graph_node=left_node,
            right_graph_node=right_node,
            left_orientation=left_orientation,
            right_orientation=right_orientation,
            graph_status=status,
            direct_edge=False,
            direct_edge_orientations=".",
            direct_edge_overlap_bp=".",
            shortest_path_edges=None,
            shortest_path_nodes=".",
            candidate_intermediate_nodes=".",
            max_path_edges=max_path_edges,
        )

    direct = graph.direct_edges(
        left_node,
        right_node,
        source_orientation=left_orientation if left_orientation in ORIENTATIONS else None,
        target_orientation=right_orientation if right_orientation in ORIENTATIONS else None,
    )
    any_orientation_direct = graph.direct_edges(left_node, right_node)
    path = find_shortest_graph_path(
        graph,
        left_node,
        right_node,
        left_orientation,
        right_orientation,
        max_path_edges,
    )

    if direct:
        status = "direct_edge"
    elif any_orientation_direct:
        status = "direct_edge_orientation_mismatch"
    elif path:
        status = "short_path"
    else:
        status = "no_path"

    return GraphGapRecord(
        scaffold=scaffold.name,
        left_contig=left.assignment.new_name,
        right_contig=right.assignment.new_name,
        left_graph_node=left_node,
        right_graph_node=right_node,
        left_orientation=left_orientation,
        right_orientation=right_orientation,
        graph_status=status,
        direct_edge=bool(direct),
        direct_edge_orientations=edge_orientations(direct or any_orientation_direct),
        direct_edge_overlap_bp=edge_overlap_bps(direct or any_orientation_direct),
        shortest_path_edges=len(path) if path else None,
        shortest_path_nodes=graph_path_nodes(path),
        candidate_intermediate_nodes=graph_intermediate_nodes(path),
        max_path_edges=max_path_edges,
    )


def build_graph_gap_records(scaffolds, graph, max_path_edges):
    records = []
    for scaffold in scaffolds:
        for left, right in zip(scaffold.members, scaffold.members[1:]):
            records.append(graph_gap_record(scaffold, left, right, graph, max_path_edges))
    return records


def scaffold_header(scaffold, simple_headers, gap_mode):
    if simple_headers:
        return scaffold.name
    fields = [
        scaffold.name,
        f"contigs={len(scaffold.members)}",
        f"sequence_bp={scaffold.sequence_bp}",
        f"gap_bp={scaffold.gap_bp}",
        f"gap_mode={gap_mode}",
    ]
    return " ".join(str(field) for field in fields)


def write_scaffold_fasta(path, scaffolds, unassigned, simple_headers, gap_mode):
    ensure_parent_dir(path)
    with open(path, "w") as out:
        for scaffold in scaffolds:
            out.write(f">{scaffold_header(scaffold, simple_headers, gap_mode)}\n")
            write_wrapped(out, scaffold.seq)
        for record in unassigned:
            out.write(record.header + "\n")
            write_wrapped(out, record.seq)


def write_gap_report(path, scaffolds):
    header = [
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
        "overlap_left_ref_frac",
        "overlap_right_ref_frac",
        "overlap_policy",
        "graph_overlap_policy",
        "overlap_action",
        "graph_overlap_action",
        "trimmed_bp",
        "sequence_overlap_identity",
    ]
    ensure_parent_dir(path)
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for scaffold in scaffolds:
            for gap in scaffold.gaps:
                row = [
                    gap.scaffold,
                    gap.left_contig,
                    gap.right_contig,
                    gap.left_ref_end,
                    gap.right_ref_start,
                    gap.raw_inferred_gap_bp,
                    gap.gap_bp,
                    gap.gap_mode,
                    gap.overlap_bp,
                    gap.overlap_class,
                    fmt(gap.overlap_left_ref_frac),
                    fmt(gap.overlap_right_ref_frac),
                    gap.overlap_policy,
                    gap.graph_overlap_policy,
                    gap.overlap_action,
                    gap.graph_overlap_action,
                    gap.trimmed_bp,
                    fmt(gap.sequence_overlap_identity, 3),
                ]
                out.write("\t".join(str(item) for item in row) + "\n")


def write_graph_gap_report(path, records):
    header = [
        "scaffold",
        "left_contig",
        "right_contig",
        "left_graph_node",
        "right_graph_node",
        "left_orientation",
        "right_orientation",
        "graph_status",
        "direct_edge",
        "direct_edge_orientations",
        "direct_edge_overlap_bp",
        "shortest_path_edges",
        "shortest_path_nodes",
        "candidate_intermediate_nodes",
        "max_path_edges",
    ]
    ensure_parent_dir(path)
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for record in records:
            row = [
                record.scaffold,
                record.left_contig,
                record.right_contig,
                record.left_graph_node,
                record.right_graph_node,
                record.left_orientation,
                record.right_orientation,
                record.graph_status,
                "yes" if record.direct_edge else "no",
                record.direct_edge_orientations,
                record.direct_edge_overlap_bp,
                record.shortest_path_edges if record.shortest_path_edges is not None else ".",
                record.shortest_path_nodes,
                record.candidate_intermediate_nodes,
                record.max_path_edges,
            ]
            out.write("\t".join(str(item) for item in row) + "\n")


def write_summary(path, scaffolds, unassigned):
    header = [
        "scaffold",
        "contigs",
        "scaffold_bp",
        "sequence_bp",
        "gap_bp",
        "gaps",
        "overlap_gaps",
        "overlap_bp",
        "trimmed_bp",
        "first_ref_start",
        "last_ref_end",
        "ordered_contigs",
    ]
    ensure_parent_dir(path)
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for scaffold in scaffolds:
            first_ref_start = scaffold.members[0].assignment.ref_start if scaffold.members else "."
            last_ref_end = scaffold.members[-1].assignment.ref_end if scaffold.members else "."
            row = [
                scaffold.name,
                len(scaffold.members),
                len(scaffold.seq),
                scaffold.sequence_bp,
                scaffold.gap_bp,
                len(scaffold.gaps),
                sum(1 for gap in scaffold.gaps if gap.overlap_bp > 0),
                scaffold.overlap_bp,
                scaffold.trimmed_bp,
                first_ref_start,
                last_ref_end,
                ",".join(member.assignment.new_name for member in scaffold.members),
            ]
            out.write("\t".join(str(item) for item in row) + "\n")
        for record in unassigned:
            row = [
                record.name,
                1,
                len(record.seq),
                len(record.seq),
                0,
                0,
                0,
                0,
                0,
                ".",
                ".",
                record.name,
            ]
            out.write("\t".join(str(item) for item in row) + "\n")


def write_run_summary(path, args, output_paths, scaffolds, unassigned):
    gap_mode = "fixed" if args.fixed_gap_bp is not None else "inferred"
    ensure_parent_dir(path)
    with open(path, "w") as out:
        out.write("chromo scaffold\n")
        out.write("\nInputs\n")
        out.write(f"ordered_fasta\t{args.ordered_fasta}\n")
        out.write(f"assignments\t{args.assignments}\n")
        out.write("\nGap model\n")
        out.write(f"gap_mode\t{gap_mode}\n")
        out.write(f"fixed_gap_bp\t{args.fixed_gap_bp if args.fixed_gap_bp is not None else '.'}\n")
        out.write(f"overlap_policy\t{args.overlap_policy}\n")
        out.write(f"graph_overlap_policy\t{args.graph_overlap_policy}\n")
        out.write(f"trim_sequence_min_identity\t{args.trim_sequence_min_identity}\n")
        out.write("\nGraph evidence\n")
        out.write(f"gfa\t{args.gfa if args.gfa else '.'}\n")
        out.write(f"graph_max_path_edges\t{args.graph_max_path_edges}\n")
        out.write("\nOutputs\n")
        for label, value in output_paths.items():
            out.write(f"{label}\t{value}\n")
        out.write("\nSummary\n")
        out.write(f"scaffolds\t{len(scaffolds)}\n")
        out.write(f"unassigned_records\t{len(unassigned)}\n")
        out.write(f"input_contigs_scaffolded\t{sum(len(scaffold.members) for scaffold in scaffolds)}\n")
        out.write(f"total_gap_bp\t{sum(scaffold.gap_bp for scaffold in scaffolds)}\n")
        out.write(f"total_overlap_bp\t{sum(scaffold.overlap_bp for scaffold in scaffolds)}\n")
        out.write(f"total_trimmed_bp\t{sum(scaffold.trimmed_bp for scaffold in scaffolds)}\n")


def run(args):
    if args.fixed_gap_bp is not None and args.fixed_gap_bp < 0:
        raise ValueError("--fixed-gap-bp must be zero or greater")
    if not 0.0 <= args.trim_sequence_min_identity <= 1.0:
        raise ValueError("--trim-sequence-min-identity must be between 0 and 1")
    if args.graph_max_path_edges < 1:
        raise ValueError("--graph-max-path-edges must be at least 1")
    if args.graph_overlap_policy != "report" and not args.gfa:
        raise ValueError("--graph-overlap-policy warn/confirm requires --gfa")

    prefix = Path(args.output_prefix)

    output_paths = {
        "scaffold_fasta": Path(str(prefix) + ".scaffold.fa"),
        "gap_report": Path(str(prefix) + ".scaffold_gaps.tsv"),
        "scaffold_summary": Path(str(prefix) + ".scaffold_summary.tsv"),
        "run_summary": Path(str(prefix) + ".run_summary.txt"),
    }
    if args.gfa:
        output_paths["graph_gap_report"] = Path(str(prefix) + ".graph_gaps.tsv")
    ensure_output_dirs(output_paths)

    assignments = read_assignments(args.assignments)
    records = read_ordered_fasta(args.ordered_fasta)
    groups, unassigned = group_scaffold_members(records, assignments)
    graph = read_gfa(args.gfa) if args.gfa else None
    scaffolds = build_scaffolds(groups, args.fixed_gap_bp, args, graph)
    gap_mode = "fixed" if args.fixed_gap_bp is not None else "inferred"
    graph_gap_records = []
    if graph is not None:
        graph_gap_records = build_graph_gap_records(
            scaffolds,
            graph,
            args.graph_max_path_edges,
        )

    write_scaffold_fasta(output_paths["scaffold_fasta"], scaffolds, unassigned, args.simple_headers, gap_mode)
    write_gap_report(output_paths["gap_report"], scaffolds)
    if graph is not None:
        write_graph_gap_report(output_paths["graph_gap_report"], graph_gap_records)
    write_summary(output_paths["scaffold_summary"], scaffolds, unassigned)
    write_run_summary(output_paths["run_summary"], args, output_paths, scaffolds, unassigned)

    sys.stderr.write(f"Wrote scaffold FASTA: {output_paths['scaffold_fasta']}\n")
    sys.stderr.write(f"Wrote gap report: {output_paths['gap_report']}\n")
    if graph is not None:
        sys.stderr.write(f"Wrote graph gap report: {output_paths['graph_gap_report']}\n")
    sys.stderr.write(f"Wrote scaffold summary: {output_paths['scaffold_summary']}\n")


def main(argv: Optional[Sequence[str]] = None, prog: Optional[str] = None):
    args = parse_args(argv, prog=prog)
    try:
        run(args)
    except ValueError as exc:
        sys.stderr.write(f"ERROR: {exc}\n")
        sys.exit(2)


if __name__ == "__main__":
    main()
