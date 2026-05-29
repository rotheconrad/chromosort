#!/usr/bin/env python3
"""Generate and apply ChromoSort manual-edit dashboards."""

import argparse
import json
import sys
from collections import OrderedDict, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

from . import __version__
from .graph import graph_node_evidence, read_gfa
from .paths import ensure_parent_dir
from .reference_order import (
    alignment_source_from_args,
    build_match_metrics,
    iter_alignments,
    iter_fasta_records,
    read_fasta_lengths,
    reverse_complement,
    write_wrapped,
)
from .review import read_review_events


SCHEMA = "chromosort-manual-v1"
MANUAL_TASKS = {"fix", "scaffold", "gapfill"}


@dataclass
class ManualPiece:
    id: str
    source: str
    name: str
    start: int
    end: int
    strand: str = "+"
    scaffold: str = "unplaced"
    removed: bool = False

    @property
    def length(self):
        return self.end - self.start + 1


def parse_generate_args(argv=None, prog=None, manual_task="general"):
    ap = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Generate a self-contained HTML dashboard for manually editing an "
            "assembly against reference dot plots."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--ref-fasta", required=True, help="Reference FASTA.")
    ap.add_argument(
        "--ref-fai",
        default=None,
        help="Reference FASTA index. Defaults to <ref-fasta>.fai when present.",
    )
    ap.add_argument(
        "-f",
        "--assembly-fasta",
        required=True,
        help="Assembly FASTA whose contigs should be edited.",
    )
    ap.add_argument(
        "--assembly-fai",
        default=None,
        help="Assembly FASTA index. Defaults to <assembly-fasta>.fai when present.",
    )
    alignment_group = ap.add_mutually_exclusive_group(required=True)
    alignment_group.add_argument("-c", "--coords", help="MUMmer show-coords file.")
    alignment_group.add_argument("--paf", help="minimap2 PAF file.")
    ap.add_argument(
        "-o",
        "--output-html",
        required=True,
        help="Output self-contained HTML dashboard.",
    )
    ap.add_argument(
        "--suggested-output-fasta",
        default=None,
        help="Suggested FASTA filename used by the dashboard export button.",
    )
    ap.add_argument(
        "--embed-sequences",
        action="store_true",
        help=(
            "Embed assembly sequences in the HTML so FASTA export works without "
            "loading the original FASTA in the browser. Best for small/demo assemblies."
        ),
    )
    ap.add_argument(
        "--gfa",
        default=None,
        help=(
            "Optional assembly graph GFA. When provided, the dashboard embeds "
            "per-contig graph node context for manual review."
        ),
    )
    ap.add_argument(
        "--read-paf",
        default=None,
        help=(
            "Optional long-read-to-assembly PAF. In task dashboards, matching "
            "longread_* review-table fields are shown as a separate evidence panel."
        ),
    )
    ap.add_argument(
        "--min-read-mapq",
        type=int,
        default=0,
        help="MAPQ threshold recorded for optional long-read PAF evidence panels.",
    )
    ap.add_argument(
        "--gaf",
        default=None,
        help=(
            "Optional long-read-to-graph GAF. In task dashboards, matching "
            "gaf_* review-table fields are shown as a separate evidence panel."
        ),
    )
    ap.add_argument(
        "--min-gaf-mapq",
        type=int,
        default=20,
        help="MAPQ threshold recorded for optional GAF evidence panels.",
    )
    ap.add_argument(
        "--review-table",
        default=None,
        help=(
            "Optional shared review-event TSV from chromo eval. In task modes, "
            "the dashboard opens around these candidate events."
        ),
    )
    ap.add_argument(
        "--min-segment-bp",
        type=int,
        default=0,
        help="Minimum query-aligned bp for an alignment row to appear in the dashboard.",
    )
    ap.add_argument(
        "--min-segment-idy",
        type=float,
        default=0.0,
        help="Ignore alignment rows below this percent identity.",
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
        "--max-segments",
        type=int,
        default=0,
        help="Maximum number of alignment rows embedded after filtering; 0 means no limit.",
    )
    args = ap.parse_args(argv)
    if args.min_read_mapq < 0:
        ap.error("--min-read-mapq must be >= 0.")
    if args.min_gaf_mapq < 0:
        ap.error("--min-gaf-mapq must be >= 0.")
    args.manual_task = manual_task
    return args


def parse_apply_args(argv=None, prog=None):
    ap = argparse.ArgumentParser(
        prog=prog,
        description="Apply a ChromoSort manual-edit recipe to an assembly FASTA.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument(
        "-f",
        "--assembly-fasta",
        required=True,
        help="Original assembly FASTA referenced by the manual recipe.",
    )
    ap.add_argument(
        "--recipe",
        required=True,
        help="Recipe JSON exported by the manual dashboard.",
    )
    ap.add_argument(
        "-o",
        "--output-fasta",
        required=True,
        help="Output FASTA with manual edits applied.",
    )
    ap.add_argument(
        "--report",
        default=None,
        help="Optional TSV report describing emitted pieces/scaffolds.",
    )
    ap.add_argument(
        "--scaffold",
        action=argparse.BooleanOptionalAction,
        default=None,
        help=(
            "Override the recipe scaffold setting. When enabled, active pieces "
            "are grouped by their scaffold labels and joined with N gaps."
        ),
    )
    ap.add_argument(
        "--gap-bp",
        type=int,
        default=None,
        help="Override the recipe scaffold gap size.",
    )
    return ap.parse_args(argv)


def fmt_float(value, digits=4):
    if value is None:
        return None
    return round(float(value), digits)


def load_dashboard_segments(args, ref_by_name, query_by_name):
    alignment_path, alignment_format = alignment_source_from_args(args)
    segments = []
    skipped = defaultdict(int)

    for seg in iter_alignments(
        alignment_path,
        alignment_format,
        min_identity=args.min_segment_idy,
        min_mapq=args.min_mapq,
        include_secondary_paf=args.include_secondary_paf,
    ):
        if seg.ref not in ref_by_name:
            skipped["unknown_ref"] += 1
            continue
        if seg.query not in query_by_name:
            skipped["unknown_query"] += 1
            continue
        if seg.len_query < args.min_segment_bp:
            skipped["short_segment"] += 1
            continue
        segments.append(
            {
                "ref": seg.ref,
                "query": seg.query,
                "refStart": min(seg.ref_start, seg.ref_end),
                "refEnd": max(seg.ref_start, seg.ref_end),
                "queryStart": min(seg.query_start, seg.query_end),
                "queryEnd": max(seg.query_start, seg.query_end),
                "orientation": seg.orientation,
                "identity": fmt_float(seg.identity, 3),
                "queryBp": seg.len_query,
                "refBp": seg.len_ref,
                "mapq": seg.mapq,
            }
        )
        if args.max_segments and len(segments) >= args.max_segments:
            break

    return segments, dict(skipped)


def sorted_best_matches(by_query):
    best = {}
    all_matches = {}
    for query, matches in by_query.items():
        ranked = sorted(
            matches,
            key=lambda m: (m.merged_query_bp, m.raw_query_bp, m.avg_identity),
            reverse=True,
        )
        all_matches[query] = ranked
        if ranked:
            best[query] = ranked[0]
    return best, all_matches


def match_to_dict(match):
    return {
        "ref": match.ref,
        "queryLength": match.query_length,
        "refLength": match.ref_length,
        "segmentCount": match.segment_count,
        "rawQueryBp": match.raw_query_bp,
        "mergedQueryBp": match.merged_query_bp,
        "queryCov": fmt_float(match.query_cov),
        "mergedRefBp": match.merged_ref_bp,
        "refCov": fmt_float(match.ref_cov),
        "refStart": match.ref_start,
        "refEnd": match.ref_end,
        "refMidpoint": fmt_float(match.ref_midpoint, 1),
        "orientation": match.orientation,
        "avgIdentity": fmt_float(match.avg_identity, 3),
    }


def graph_edge_to_dict(edge, direction):
    other_node = edge.source if direction == "incoming" else edge.target
    other_orientation = (
        edge.source_orientation if direction == "incoming" else edge.target_orientation
    )
    return {
        "direction": direction,
        "source": edge.source,
        "sourceOrientation": edge.source_orientation,
        "target": edge.target,
        "targetOrientation": edge.target_orientation,
        "otherNode": other_node,
        "otherOrientation": other_orientation,
        "overlap": edge.overlap,
        "overlapBp": edge.overlap_bp,
    }


def graph_complexity(evidence):
    if evidence.graph_node_status != "present":
        return "missing"
    if evidence.graph_self_loop:
        return "self_loop"
    if (
        (evidence.graph_in_degree or 0) > 1
        or (evidence.graph_out_degree or 0) > 1
        or (evidence.graph_neighbor_count or 0) > 2
    ):
        return "branching"
    return "simple"


def graph_context_to_dict(graph, evidence):
    data = {
        "graphNode": evidence.graph_node,
        "graphNodeStatus": evidence.graph_node_status,
        "graphNodeLength": evidence.graph_node_length,
        "graphNodeHasSequence": evidence.graph_node_has_sequence,
        "graphInDegree": evidence.graph_in_degree,
        "graphOutDegree": evidence.graph_out_degree,
        "graphNeighborCount": evidence.graph_neighbor_count,
        "graphSelfLoop": evidence.graph_self_loop,
        "graphComplexity": graph_complexity(evidence),
        "graphCoverage": None,
        "incoming": [],
        "outgoing": [],
    }
    if evidence.graph_node_status != "present":
        return data

    node = graph.nodes[evidence.graph_node]
    coverage = node.tags.get("RC", node.tags.get("KC", node.tags.get("dp")))
    if isinstance(coverage, (int, float, str)):
        data["graphCoverage"] = coverage
    data["incoming"] = [
        graph_edge_to_dict(edge, "incoming")
        for edge in graph.incoming(evidence.graph_node)
    ]
    data["outgoing"] = [
        graph_edge_to_dict(edge, "outgoing")
        for edge in graph.outgoing(evidence.graph_node)
    ]
    return data


def build_graph_context(query_records, graph):
    if graph is None:
        return {}
    return {
        rec.name: graph_context_to_dict(
            graph,
            graph_node_evidence(graph, [rec.name]),
        )
        for rec in query_records
    }


def build_query_items(query_records, best_matches, all_matches, ref_rank, graph_context=None):
    graph_context = graph_context or {}
    indexed_records = list(enumerate(query_records))

    def order_key(item):
        assembly_index, rec = item
        best = best_matches.get(rec.name)
        if best is None:
            return (1, assembly_index)
        return (
            0,
            ref_rank.get(best.ref, 10**12),
            best.ref_start,
            best.ref_end,
            rec.name,
        )

    ordered = []
    for initial_rank, (assembly_index, rec) in enumerate(sorted(indexed_records, key=order_key), start=1):
        best = best_matches.get(rec.name)
        matches = all_matches.get(rec.name, [])
        total_matched = sum(match.merged_query_bp for match in matches)
        best_share = best.merged_query_bp / total_matched if best and total_matched else 0.0
        ordered.append(
            {
                "name": rec.name,
                "length": rec.length,
                "assemblyOrder": assembly_index + 1,
                "initialRank": initial_rank,
                "aligned": best is not None,
                "bestRef": best.ref if best else None,
                "refStart": best.ref_start if best else None,
                "refEnd": best.ref_end if best else None,
                "refMidpoint": fmt_float(best.ref_midpoint, 1) if best else None,
                "orientation": best.orientation if best else None,
                "mergedQueryBp": best.merged_query_bp if best else 0,
                "queryCov": fmt_float(best.query_cov) if best else 0.0,
                "avgIdentity": fmt_float(best.avg_identity, 3) if best else None,
                "bestRefShare": fmt_float(best_share) if best else 0.0,
                "totalRefsMatched": len(matches),
                "matches": [match_to_dict(match) for match in matches],
                "graph": graph_context.get(rec.name),
            }
        )
    return ordered


def annotate_graph_neighbor_alignment(query_items):
    by_name = {item["name"]: item for item in query_items}
    for item in query_items:
        graph = item.get("graph")
        if not graph:
            continue
        item_ref = item.get("bestRef")
        for direction in ("incoming", "outgoing"):
            for edge in graph.get(direction, []):
                neighbor = by_name.get(edge.get("otherNode"))
                neighbor_ref = neighbor.get("bestRef") if neighbor else None
                edge["otherAligned"] = bool(neighbor and neighbor.get("aligned"))
                edge["otherBestRef"] = neighbor_ref
                edge["sameBestRef"] = bool(item_ref and neighbor_ref == item_ref)


def build_initial_pieces(query_items):
    pieces = []
    for index, item in enumerate(query_items, start=1):
        pieces.append(
            {
                "id": f"p{index}",
                "source": item["name"],
                "name": item["name"],
                "start": 1,
                "end": item["length"],
                "strand": "+",
                "scaffold": item["bestRef"] or "unplaced",
                "removed": False,
                "aligned": item["aligned"],
                "bestRef": item["bestRef"],
            }
        )
    return pieces


def review_event_sources(event):
    fields = event.fields
    candidates = [
        fields.get("source_contig"),
        fields.get("new_contig"),
        fields.get("left_contig"),
        fields.get("right_contig"),
        fields.get("scaffold"),
        event.target,
    ]
    sources = []
    for candidate in candidates:
        if candidate in {None, "", "."}:
            continue
        text = str(candidate)
        for delimiter in (":", "|", ","):
            text = text.replace(delimiter, " ")
        for part in text.split():
            if part and part != "." and part not in sources:
                sources.append(part)
    return sources


def review_event_to_dict(event):
    return {
        "event_id": event.event_id,
        "task": event.task,
        "action": event.action,
        "target": event.target,
        "accept": event.accept,
        "status": event.status,
        "confidence": event.confidence,
        "reason": event.reason,
        "notes": event.notes,
        "fields": event.fields,
        "sources": review_event_sources(event),
    }


def load_review_events(path, task):
    if not path:
        return []
    expected_task = task if task in MANUAL_TASKS else None
    return [
        review_event_to_dict(event)
        for event in read_review_events(path, expected_task=expected_task)
    ]


def build_dashboard_data(args):
    alignment_path, alignment_format = alignment_source_from_args(args)
    ref_records, ref_by_name = read_fasta_lengths(args.ref_fasta, args.ref_fai)
    query_records, query_by_name = read_fasta_lengths(args.assembly_fasta, args.assembly_fai)
    ref_lengths = {name: rec for name, rec in ref_by_name.items()}
    query_lengths = {name: rec for name, rec in query_by_name.items()}

    matches, by_query, skipped_unknown_query = build_match_metrics(
        alignment_path,
        alignment_format,
        ref_lengths,
        query_lengths,
        args.min_segment_idy,
        min_mapq=args.min_mapq,
        include_secondary_paf=args.include_secondary_paf,
    )
    best_matches, all_matches = sorted_best_matches(by_query)
    ref_rank = {rec.name: idx for idx, rec in enumerate(ref_records)}
    graph = read_gfa(args.gfa) if args.gfa else None
    graph_context = build_graph_context(query_records, graph)
    query_items = build_query_items(
        query_records,
        best_matches,
        all_matches,
        ref_rank,
        graph_context,
    )
    annotate_graph_neighbor_alignment(query_items)
    segments, skipped = load_dashboard_segments(args, ref_by_name, query_by_name)
    if skipped_unknown_query:
        skipped["unknown_query_in_metrics"] = skipped_unknown_query

    sequences = {}
    if args.embed_sequences:
        sequences = {name: seq for name, _, seq in iter_fasta_records(args.assembly_fasta)}

    suggested = args.suggested_output_fasta
    if not suggested:
        suggested = str(Path(args.assembly_fasta).with_suffix(".manual.fa").name)

    manual_task = getattr(args, "manual_task", "general")
    review_events = load_review_events(args.review_table, manual_task)

    return {
        "schema": SCHEMA,
        "version": __version__,
        "mode": f"manual-{manual_task}" if manual_task in MANUAL_TASKS else "manual-dashboard",
        "reviewTask": manual_task,
        "inputs": {
            "refFasta": str(args.ref_fasta),
            "assemblyFasta": str(args.assembly_fasta),
            "alignmentFormat": alignment_format,
            "alignmentPath": str(alignment_path),
            "gfa": str(args.gfa) if args.gfa else None,
            "readPaf": str(args.read_paf) if args.read_paf else None,
            "gaf": str(args.gaf) if args.gaf else None,
            "reviewTable": str(args.review_table) if args.review_table else None,
        },
        "settings": {
            "minSegmentBp": args.min_segment_bp,
            "minSegmentIdy": args.min_segment_idy,
            "minMapq": args.min_mapq,
            "minReadMapq": args.min_read_mapq,
            "minGafMapq": args.min_gaf_mapq,
            "includeSecondaryPaf": bool(args.include_secondary_paf),
            "maxSegments": args.max_segments,
            "embedSequences": bool(args.embed_sequences),
            "suggestedOutputFasta": suggested,
        },
        "stats": {
            "referenceSequences": len(ref_records),
            "assemblyContigs": len(query_records),
            "alignedContigs": sum(1 for item in query_items if item["aligned"]),
            "unalignedContigs": sum(1 for item in query_items if not item["aligned"]),
            "embeddedSegments": len(segments),
            "allMatchRows": len(matches),
            "embeddedSequenceCount": len(sequences),
            "embeddedSequenceBp": sum(len(seq) for seq in sequences.values()),
            "reviewEvents": len(review_events),
            "graphContigsPresent": sum(
                1 for item in query_items
                if item.get("graph")
                and item["graph"].get("graphNodeStatus") == "present"
            ),
            "graphContigsMissing": sum(
                1 for item in query_items
                if item.get("graph")
                and item["graph"].get("graphNodeStatus") != "present"
            ),
            "graphBranchingContigs": sum(
                1 for item in query_items
                if item.get("graph")
                and item["graph"].get("graphComplexity") in {"branching", "self_loop"}
            ),
            "skipped": skipped,
        },
        "refRecords": [
            {"name": rec.name, "length": rec.length, "order": idx + 1}
            for idx, rec in enumerate(ref_records)
        ],
        "queryRecords": query_items,
        "initialPieces": build_initial_pieces(query_items),
        "segments": segments,
        "sequences": sequences,
        "reviewEvents": review_events,
    }


def json_for_script(data):
    return json.dumps(data, ensure_ascii=True, separators=(",", ":")).replace("</", "<\\/")


def write_dashboard(path, data):
    ensure_parent_dir(path)
    html = DASHBOARD_HTML.replace("__CHROMOSORT_MANUAL_DATA__", json_for_script(data))
    with open(path, "w") as out:
        out.write(html)


def normalize_recipe_piece(raw, index):
    try:
        source = str(raw["source"])
        start = int(raw["start"])
        end = int(raw["end"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"Malformed recipe piece at index {index}.") from exc

    name = str(raw.get("name") or f"{source}_manual{index:03d}")
    piece_id = str(raw.get("id") or f"p{index}")
    strand = str(raw.get("strand") or "+")
    scaffold = str(raw.get("scaffold") or "unplaced")
    removed = bool(raw.get("removed", False))
    if strand not in {"+", "-"}:
        raise ValueError(f"Recipe piece {piece_id!r} has invalid strand {strand!r}.")
    if start < 1 or end < start:
        raise ValueError(f"Recipe piece {piece_id!r} has invalid slice {start}-{end}.")
    return ManualPiece(
        id=piece_id,
        source=source,
        name=name,
        start=start,
        end=end,
        strand=strand,
        scaffold=scaffold,
        removed=removed,
    )


def read_recipe(path):
    with open(path) as fh:
        recipe = json.load(fh)
    if recipe.get("schema") != SCHEMA:
        raise ValueError(f"Recipe {path} is not a {SCHEMA} recipe.")
    raw_pieces = recipe.get("pieces")
    if not isinstance(raw_pieces, list):
        raise ValueError(f"Recipe {path} is missing a pieces list.")
    pieces = [
        normalize_recipe_piece(raw, index)
        for index, raw in enumerate(raw_pieces, start=1)
    ]
    return recipe, pieces


def active_pieces(pieces):
    return [piece for piece in pieces if not piece.removed]


def fetch_piece_sequences(fasta_path, pieces):
    wanted = {piece.source for piece in pieces}
    sequences = {}
    for name, _, seq in iter_fasta_records(fasta_path):
        if name in wanted:
            sequences[name] = seq
    missing = sorted(wanted - set(sequences))
    if missing:
        preview = ", ".join(missing[:5])
        suffix = "..." if len(missing) > 5 else ""
        raise ValueError(f"Assembly FASTA is missing recipe source contig(s): {preview}{suffix}.")
    return sequences


def piece_sequence(piece, sequences):
    seq = sequences[piece.source]
    if piece.end > len(seq):
        raise ValueError(
            f"Recipe piece {piece.id!r} slice {piece.start}-{piece.end} exceeds "
            f"source contig {piece.source!r} length {len(seq)}."
        )
    out = seq[piece.start - 1 : piece.end]
    if piece.strand == "-":
        out = reverse_complement(out)
    return out


def unique_piece_names(pieces):
    counts = defaultdict(int)
    names = {}
    for index, piece in enumerate(pieces, start=1):
        base = piece.name or f"{piece.source}_manual{index:03d}"
        counts[base] += 1
        names[piece.id] = base if counts[base] == 1 else f"{base}_{counts[base]}"
    return names


def group_scaffold_pieces(pieces):
    groups = OrderedDict()
    for piece in pieces:
        scaffold = piece.scaffold or "unplaced"
        groups.setdefault(scaffold, []).append(piece)
    return groups


def write_apply_fasta(path, fasta_path, pieces, scaffold_enabled, gap_bp):
    if gap_bp < 0:
        raise ValueError("Manual scaffold gap bp must be non-negative.")
    ensure_parent_dir(path)
    sequences = fetch_piece_sequences(fasta_path, pieces)
    with open(path, "w") as out:
        if scaffold_enabled:
            for scaffold, members in group_scaffold_pieces(pieces).items():
                seq_parts = []
                for idx, piece in enumerate(members):
                    if idx:
                        seq_parts.append("N" * gap_bp)
                    seq_parts.append(piece_sequence(piece, sequences))
                out.write(f">{scaffold} manual_pieces={len(members)} gap_bp={gap_bp}\n")
                write_wrapped(out, "".join(seq_parts))
        else:
            names = unique_piece_names(pieces)
            for piece in pieces:
                out.write(
                    f">{names[piece.id]} original={piece.source} "
                    f"slice={piece.start}-{piece.end} strand={piece.strand} "
                    f"scaffold={piece.scaffold}\n"
                )
                write_wrapped(out, piece_sequence(piece, sequences))


def write_apply_report(path, pieces, scaffold_enabled):
    ensure_parent_dir(path)
    header = [
        "piece_id",
        "source",
        "output_name",
        "scaffold",
        "slice_start",
        "slice_end",
        "piece_bp",
        "strand",
        "removed",
        "export_mode",
    ]
    names = unique_piece_names(active_pieces(pieces))
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for piece in pieces:
            row = [
                piece.id,
                piece.source,
                "." if piece.removed else names.get(piece.id, piece.name),
                piece.scaffold,
                piece.start,
                piece.end,
                piece.length,
                piece.strand,
                "yes" if piece.removed else "no",
                "scaffold" if scaffold_enabled else "pieces",
            ]
            out.write("\t".join(str(item) for item in row) + "\n")


def run_generate(args):
    data = build_dashboard_data(args)
    write_dashboard(args.output_html, data)
    sys.stderr.write(f"Wrote manual dashboard: {args.output_html}\n")
    sys.stderr.write(
        f"Embedded {data['stats']['embeddedSegments']} alignment row(s) for "
        f"{data['stats']['assemblyContigs']} contig(s).\n"
    )
    if args.embed_sequences:
        sys.stderr.write(
            f"Embedded {data['stats']['embeddedSequenceCount']} sequence(s), "
            f"{data['stats']['embeddedSequenceBp']} bp.\n"
        )


def run_apply(args):
    recipe, pieces = read_recipe(args.recipe)
    active = active_pieces(pieces)
    scaffold_enabled = (
        bool(recipe.get("scaffoldEnabled", False))
        if args.scaffold is None
        else bool(args.scaffold)
    )
    gap_bp = args.gap_bp
    if gap_bp is None:
        gap_bp = int(recipe.get("gapBp", 100))
    write_apply_fasta(args.output_fasta, args.assembly_fasta, active, scaffold_enabled, gap_bp)
    if args.report:
        write_apply_report(args.report, pieces, scaffold_enabled)
    sys.stderr.write(
        f"Applied manual recipe with {len(active)} active piece(s) "
        f"and {len(pieces) - len(active)} removed piece(s).\n"
    )
    sys.stderr.write(f"Wrote manual FASTA: {args.output_fasta}\n")
    if args.report:
        sys.stderr.write(f"Wrote manual report: {args.report}\n")


def main(argv: Optional[Sequence[str]] = None, prog: Optional[str] = None):
    if argv is None:
        argv = sys.argv[1:]

    try:
        if argv and argv[0] == "apply":
            apply_prog = f"{prog} apply" if prog else None
            args = parse_apply_args(argv[1:], prog=apply_prog)
            run_apply(args)
        elif argv and argv[0] in MANUAL_TASKS:
            mode = argv[0]
            mode_prog = f"{prog} {mode}" if prog else None
            args = parse_generate_args(argv[1:], prog=mode_prog, manual_task=mode)
            run_generate(args)
        else:
            args = parse_generate_args(argv, prog=prog)
            run_generate(args)
    except ValueError as exc:
        sys.stderr.write(f"ERROR: {exc}\n")
        sys.exit(2)


DASHBOARD_HTML = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>ChromoSort Manual</title>
<style>
:root {
  --bg: #f8fafc;
  --panel: #ffffff;
  --panel-2: #f1f5f9;
  --ink: #111827;
  --muted: #64748b;
  --line: #cbd5e1;
  --blue: #2563eb;
  --red: #dc2626;
  --green: #15803d;
  --amber: #a16207;
  --focus: #0f766e;
}
* { box-sizing: border-box; }
body {
  margin: 0;
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Arial, sans-serif;
  background: var(--bg);
  color: var(--ink);
}
header {
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 16px;
  min-height: 58px;
  padding: 10px 16px;
  border-bottom: 1px solid var(--line);
  background: var(--panel);
}
h1 {
  margin: 0;
  font-size: 18px;
  font-weight: 700;
  letter-spacing: 0;
}
button, input, select {
  font: inherit;
}
button {
  min-height: 34px;
  border: 1px solid var(--line);
  border-radius: 6px;
  background: var(--panel);
  color: var(--ink);
  padding: 6px 10px;
  cursor: pointer;
}
button:hover { border-color: var(--focus); }
button.primary {
  color: #ffffff;
  background: var(--focus);
  border-color: var(--focus);
}
button.danger { color: #991b1b; }
button:disabled {
  cursor: not-allowed;
  opacity: 0.45;
}
input[type="text"], input[type="number"], select {
  min-height: 34px;
  border: 1px solid var(--line);
  border-radius: 6px;
  background: var(--panel);
  color: var(--ink);
  padding: 6px 8px;
}
input[type="file"] {
  max-width: 280px;
}
.app {
  display: grid;
  grid-template-columns: minmax(290px, 380px) minmax(0, 1fr);
  height: calc(100vh - 58px);
}
.sidebar {
  display: flex;
  flex-direction: column;
  min-width: 0;
  border-right: 1px solid var(--line);
  background: var(--panel);
}
.side-tools, .detail-tools, .export-tools {
  display: flex;
  flex-wrap: wrap;
  gap: 8px;
  align-items: center;
}
.side-tools {
  padding: 12px;
  border-bottom: 1px solid var(--line);
}
#search {
  flex: 1 1 140px;
  min-width: 120px;
}
#graphFilter {
  flex: 0 1 150px;
  min-width: 132px;
}
.stats {
  display: grid;
  grid-template-columns: repeat(2, minmax(0, 1fr));
  gap: 8px;
  padding: 12px;
  border-bottom: 1px solid var(--line);
  background: var(--panel-2);
}
.stat {
  min-width: 0;
}
.stat strong {
  display: block;
  font-size: 15px;
}
.stat span {
  display: block;
  color: var(--muted);
  font-size: 12px;
}
.piece-list {
  overflow: auto;
  min-height: 0;
}
.event-panel {
  display: grid;
  gap: 6px;
  padding: 10px 12px;
  border-bottom: 1px solid var(--line);
  background: #f8fafc;
}
.event-panel.hidden { display: none; }
.event-title {
  font-size: 12px;
  font-weight: 700;
  color: var(--muted);
  text-transform: uppercase;
}
.event-list {
  display: grid;
  gap: 6px;
  max-height: 190px;
  overflow: auto;
}
.event-row {
  display: grid;
  gap: 3px;
  padding: 7px;
  border: 1px solid #e5e7eb;
  border-radius: 6px;
  background: #ffffff;
  cursor: pointer;
  font-size: 12px;
}
.event-row:hover { border-color: var(--focus); }
.event-row.selected {
  outline: 2px solid var(--focus);
  outline-offset: -2px;
}
.event-row strong, .event-row span {
  min-width: 0;
  overflow: hidden;
  text-overflow: ellipsis;
  white-space: nowrap;
}
.piece-row {
  display: grid;
  grid-template-columns: minmax(0, 1fr) auto;
  gap: 8px;
  padding: 10px 12px;
  border-bottom: 1px solid #e5e7eb;
  cursor: pointer;
}
.piece-row:hover { background: #f8fafc; }
.piece-row.selected {
  outline: 2px solid var(--focus);
  outline-offset: -2px;
  background: #ecfdf5;
}
.piece-row.removed {
  opacity: 0.55;
}
.piece-name {
  overflow: hidden;
  text-overflow: ellipsis;
  white-space: nowrap;
  font-weight: 650;
}
.piece-meta {
  display: flex;
  flex-wrap: wrap;
  gap: 6px;
  margin-top: 4px;
  color: var(--muted);
  font-size: 12px;
}
.badge {
  display: inline-flex;
  align-items: center;
  min-height: 20px;
  border: 1px solid var(--line);
  border-radius: 999px;
  padding: 1px 6px;
  background: #ffffff;
}
.badge.red { color: var(--red); border-color: #fecaca; }
.badge.green { color: var(--green); border-color: #bbf7d0; }
.badge.amber { color: var(--amber); border-color: #fde68a; }
.badge.blue { color: var(--blue); border-color: #bfdbfe; }
.workspace {
  display: grid;
  grid-template-rows: auto minmax(260px, 1fr) auto;
  min-width: 0;
  min-height: 0;
}
.plot-head {
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 12px;
  padding: 12px 14px;
  border-bottom: 1px solid var(--line);
  background: var(--panel);
}
.plot-title {
  min-width: 0;
}
.plot-title strong, .plot-title span {
  display: block;
  overflow: hidden;
  text-overflow: ellipsis;
  white-space: nowrap;
}
.plot-title span {
  color: var(--muted);
  font-size: 12px;
  margin-top: 2px;
}
.plot-area {
  min-width: 0;
  min-height: 0;
  padding: 12px;
}
#dotplot {
  display: block;
  width: 100%;
  height: 100%;
  min-height: 260px;
  border: 1px solid var(--line);
  border-radius: 8px;
  background: #ffffff;
}
.details {
  display: grid;
  grid-template-columns: minmax(0, 1fr);
  gap: 10px;
  padding: 12px 14px;
  border-top: 1px solid var(--line);
  background: var(--panel);
}
.field-grid {
  display: grid;
  grid-template-columns: repeat(5, minmax(100px, 1fr));
  gap: 8px;
}
.graph-panel {
  display: grid;
  gap: 8px;
  padding: 10px;
  border: 1px solid #e5e7eb;
  border-radius: 8px;
  background: #f8fafc;
}
.graph-panel.hidden { display: none; }
.graph-panel-title {
  display: flex;
  flex-wrap: wrap;
  gap: 6px;
  align-items: center;
  font-size: 13px;
  font-weight: 700;
}
.graph-neighbor-list {
  display: grid;
  gap: 6px;
}
.graph-neighbor-row {
  display: grid;
  grid-template-columns: auto minmax(110px, 1fr) minmax(140px, 1.3fr) auto auto;
  gap: 6px;
  align-items: center;
  min-height: 30px;
  padding: 5px 7px;
  border: 1px solid #e5e7eb;
  border-left-width: 4px;
  border-radius: 6px;
  background: #ffffff;
  font-size: 12px;
}
.graph-neighbor-row.incoming { border-left-color: var(--blue); }
.graph-neighbor-row.outgoing { border-left-color: var(--focus); }
.graph-neighbor-row span, .graph-neighbor-row strong {
  min-width: 0;
  overflow: hidden;
  text-overflow: ellipsis;
  white-space: nowrap;
}
.evidence-panels {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(210px, 1fr));
  gap: 8px;
}
.evidence-panels.hidden { display: none; }
.evidence-panel {
  display: grid;
  gap: 7px;
  min-width: 0;
  padding: 10px;
  border: 1px solid #e5e7eb;
  border-radius: 8px;
  background: #ffffff;
}
.evidence-title {
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 6px;
  min-width: 0;
  font-size: 12px;
  font-weight: 750;
  color: var(--ink);
}
.evidence-title span {
  overflow: hidden;
  text-overflow: ellipsis;
  white-space: nowrap;
}
.evidence-row {
  display: grid;
  grid-template-columns: minmax(76px, 0.8fr) minmax(0, 1.2fr);
  gap: 6px;
  min-width: 0;
  font-size: 12px;
}
.evidence-row span:first-child {
  color: var(--muted);
  font-weight: 650;
}
.evidence-row span:last-child {
  min-width: 0;
  overflow: hidden;
  text-overflow: ellipsis;
  white-space: nowrap;
}
.graph-side {
  color: var(--muted);
  font-weight: 650;
  text-transform: uppercase;
  font-size: 10px;
}
.graph-note {
  color: var(--muted);
  font-size: 12px;
}
label {
  display: grid;
  gap: 4px;
  min-width: 0;
  color: var(--muted);
  font-size: 12px;
}
label input, label select {
  color: var(--ink);
  font-size: 14px;
}
.status {
  color: var(--muted);
  font-size: 12px;
  min-height: 16px;
}
.export-box {
  display: flex;
  flex-wrap: wrap;
  align-items: center;
  gap: 8px;
  padding-top: 8px;
  border-top: 1px solid #e5e7eb;
}
@media (max-width: 860px) {
  .app { grid-template-columns: 1fr; height: auto; }
  .sidebar { border-right: 0; border-bottom: 1px solid var(--line); max-height: 45vh; }
  .workspace { min-height: 650px; }
  .field-grid { grid-template-columns: repeat(2, minmax(120px, 1fr)); }
  .graph-neighbor-row { grid-template-columns: auto minmax(100px, 1fr); }
  .graph-neighbor-row span:nth-child(3),
  .graph-neighbor-row span:nth-child(4) { display: none; }
}
</style>
</head>
<body>
<header>
  <div>
    <h1>ChromoSort Manual</h1>
    <div class="status" id="datasetLabel"></div>
  </div>
  <div class="export-tools">
    <button id="exportRecipe">Export recipe</button>
    <button class="primary" id="exportFasta">Export FASTA</button>
  </div>
</header>
<main class="app">
  <aside class="sidebar">
    <div class="side-tools">
      <input id="search" type="text" placeholder="Search contigs">
      <select id="graphFilter" aria-label="Graph filter">
        <option value="">All graph</option>
        <option value="simple">simple</option>
        <option value="branching">branching</option>
        <option value="self_loop">self_loop</option>
        <option value="missing">missing</option>
      </select>
      <button id="removeUnaligned">Remove unaligned</button>
      <button id="restoreAll">Restore all</button>
    </div>
    <div class="event-panel hidden" id="eventPanel">
      <div class="event-title" id="eventTitle">Review events</div>
      <div class="event-list" id="eventList"></div>
    </div>
    <div class="stats" id="stats"></div>
    <div class="piece-list" id="pieceList"></div>
  </aside>
  <section class="workspace">
    <div class="plot-head">
      <div class="plot-title">
        <strong id="selectedName">No contig selected</strong>
        <span id="selectedMeta"></span>
        <span id="selectedGraph"></span>
      </div>
      <div class="detail-tools">
        <button id="moveUp">Up</button>
        <button id="moveDown">Down</button>
        <button id="invert">Invert</button>
        <button class="danger" id="toggleRemove">Remove</button>
      </div>
    </div>
    <div class="plot-area">
      <canvas id="dotplot"></canvas>
    </div>
    <div class="details">
      <div class="field-grid">
        <label>Cut after bp<input id="cutPos" type="number" min="1" step="1"></label>
        <label>Scaffold<input id="scaffoldName" type="text"></label>
        <label>Gap bp<input id="gapBp" type="number" min="0" step="1" value="100"></label>
        <label>FASTA file<input id="fastaFile" type="file" accept=".fa,.fasta,.fna,.txt"></label>
        <label>Recipe file<input id="recipeFile" type="file" accept=".json"></label>
      </div>
      <div class="graph-panel hidden" id="graphPanel"></div>
      <div class="evidence-panels hidden" id="evidencePanels"></div>
      <div class="detail-tools">
        <button id="addCut">Add breakpoint</button>
        <label><input id="scaffoldMode" type="checkbox"> Scaffold export</label>
      </div>
      <div class="status" id="status"></div>
    </div>
  </section>
</main>
<script>
window.CHROMOSORT_MANUAL_DATA = __CHROMOSORT_MANUAL_DATA__;
</script>
<script>
(() => {
  const data = window.CHROMOSORT_MANUAL_DATA;
  let pieces = data.initialPieces.map(p => ({...p}));
  const reviewEvents = data.reviewEvents || [];
  let selectedEventId = reviewEvents.length ? reviewEvents[0].event_id : null;
  let selectedId = reviewEvents.length
    ? (pieceIdForEvent(reviewEvents[0]) || (pieces.length ? pieces[0].id : null))
    : (pieces.length ? pieces[0].id : null);
  let nextId = pieces.length + 1;
  let sequences = {...(data.sequences || {})};
  let scaffoldEnabled = false;

  const refs = data.refRecords || [];
  const queriesByName = new Map((data.queryRecords || []).map(q => [q.name, q]));
  const segmentsByQuery = new Map();
  for (const seg of data.segments || []) {
    if (!segmentsByQuery.has(seg.query)) segmentsByQuery.set(seg.query, []);
    segmentsByQuery.get(seg.query).push(seg);
  }
  const els = {
    datasetLabel: document.getElementById("datasetLabel"),
    stats: document.getElementById("stats"),
    search: document.getElementById("search"),
    graphFilter: document.getElementById("graphFilter"),
    eventPanel: document.getElementById("eventPanel"),
    eventTitle: document.getElementById("eventTitle"),
    eventList: document.getElementById("eventList"),
    pieceList: document.getElementById("pieceList"),
    selectedName: document.getElementById("selectedName"),
    selectedMeta: document.getElementById("selectedMeta"),
    selectedGraph: document.getElementById("selectedGraph"),
    graphPanel: document.getElementById("graphPanel"),
    evidencePanels: document.getElementById("evidencePanels"),
    dotplot: document.getElementById("dotplot"),
    status: document.getElementById("status"),
    cutPos: document.getElementById("cutPos"),
    scaffoldName: document.getElementById("scaffoldName"),
    gapBp: document.getElementById("gapBp"),
    scaffoldMode: document.getElementById("scaffoldMode"),
    fastaFile: document.getElementById("fastaFile"),
    recipeFile: document.getElementById("recipeFile")
  };

  function fmtBp(value) {
    if (value >= 1000000000) return (value / 1000000000).toFixed(2) + " Gb";
    if (value >= 1000000) return (value / 1000000).toFixed(1) + " Mb";
    if (value >= 1000) return (value / 1000).toFixed(1) + " kb";
    return String(value) + " bp";
  }

  function setStatus(message) {
    els.status.textContent = message || "";
  }

  function selectedPiece() {
    return pieces.find(p => p.id === selectedId) || pieces[0] || null;
  }

  function selectedEvent() {
    return reviewEvents.find(evt => evt.event_id === selectedEventId) || null;
  }

  function eventSources(evt) {
    return (evt && evt.sources ? evt.sources : []).map(String);
  }

  function pieceIdForEvent(evt) {
    const sources = new Set(eventSources(evt));
    if (!sources.size) return null;
    const piece = pieces.find(p => sources.has(p.source) || sources.has(p.name));
    return piece ? piece.id : null;
  }

  function pieceEventCount(piece) {
    return reviewEvents.filter(evt => {
      const sources = new Set(eventSources(evt));
      return sources.has(piece.source) || sources.has(piece.name);
    }).length;
  }

  function pieceLength(piece) {
    return piece.end - piece.start + 1;
  }

  function pieceBadges(piece) {
    const q = queriesByName.get(piece.source) || {};
    const graph = q.graph || null;
    const badges = [];
    badges.push(`<span class="badge">${fmtBp(pieceLength(piece))}</span>`);
    badges.push(`<span class="badge ${piece.strand === "-" ? "red" : "green"}">${piece.strand}</span>`);
    badges.push(`<span class="badge ${q.aligned ? "green" : "amber"}">${q.bestRef || "unaligned"}</span>`);
    if (graph) {
      const graphClass = graphBadgeClass(graph);
      const graphLabel = graph.graphNodeStatus === "present"
        ? graph.graphNode
        : "graph missing";
      badges.push(`<span class="badge ${graphClass}">${escapeHtml(graphLabel)}</span>`);
      if (graph.graphComplexity && graph.graphComplexity !== "missing") {
        badges.push(`<span class="badge ${graphClass}">${escapeHtml(graph.graphComplexity)}</span>`);
      }
    }
    const eventCount = pieceEventCount(piece);
    if (eventCount) badges.push(`<span class="badge blue">${eventCount} event${eventCount === 1 ? "" : "s"}</span>`);
    if (piece.removed) badges.push('<span class="badge red">removed</span>');
    return badges.join("");
  }

  function graphBadgeClass(graph) {
    if (!graph || graph.graphNodeStatus !== "present") return "red";
    if (graph.graphComplexity === "simple") return "green";
    if (graph.graphComplexity === "branching" || graph.graphComplexity === "self_loop") return "amber";
    return "blue";
  }

  function graphNeighbors(edges) {
    return (edges || []).map(edge => {
      const overlap = edge.overlapBp === null || edge.overlapBp === undefined ? edge.overlap : edge.overlapBp + " bp";
      return edge.otherNode + edge.otherOrientation + " (" + overlap + ")";
    }).join(", ");
  }

  function graphDetail(graph) {
    if (!graph) return "";
    if (graph.graphNodeStatus !== "present") return "Graph: missing from GFA";
    const pieces = [
      "Graph: " + graph.graphNode,
      "complexity " + graph.graphComplexity,
      "in/out " + graph.graphInDegree + "/" + graph.graphOutDegree,
      "neighbors " + graph.graphNeighborCount
    ];
    if (graph.graphCoverage !== null && graph.graphCoverage !== undefined) {
      pieces.push("coverage " + graph.graphCoverage);
    }
    const incoming = graphNeighbors(graph.incoming);
    const outgoing = graphNeighbors(graph.outgoing);
    if (incoming) pieces.push("in " + incoming);
    if (outgoing) pieces.push("out " + outgoing);
    return pieces.join(" | ");
  }

  function graphFilterValue(graph) {
    if (!graph) return "missing";
    return graph.graphComplexity || "missing";
  }

  function graphNeighborBadge(edge) {
    if (edge.sameBestRef) return '<span class="badge green">same ref</span>';
    if (edge.otherAligned) {
      return `<span class="badge blue">${escapeHtml(edge.otherBestRef || "aligned")}</span>`;
    }
    return '<span class="badge amber">not placed</span>';
  }

  function graphNeighborRow(edge) {
    const isIncoming = edge.direction === "incoming";
    const side = isIncoming ? "upstream" : "downstream";
    const overlap = edge.overlapBp === null || edge.overlapBp === undefined
      ? (edge.overlap || ".")
      : edge.overlapBp + " bp";
    const orientedNode = String(edge.otherNode || ".") + String(edge.otherOrientation || "");
    const link = String(edge.source || ".") + String(edge.sourceOrientation || "") +
      " -> " + String(edge.target || ".") + String(edge.targetOrientation || "");
    return `<div class="graph-neighbor-row ${isIncoming ? "incoming" : "outgoing"}">
      <span class="graph-side">${side}</span>
      <strong title="${escapeHtml(orientedNode)}">${escapeHtml(orientedNode)}</strong>
      <span title="${escapeHtml(link)}">${escapeHtml(link)}</span>
      <span title="${escapeHtml(overlap)}">${escapeHtml(overlap)}</span>
      ${graphNeighborBadge(edge)}
    </div>`;
  }

  function eventFields(evt) {
    return (evt && evt.fields) || {};
  }

  function presentValue(value) {
    if (value === null || value === undefined) return false;
    const text = String(value).trim();
    return text !== "" && text !== ".";
  }

  function firstPresent(...values) {
    for (const value of values) {
      if (presentValue(value)) return value;
    }
    return ".";
  }

  function fieldValue(fields, names) {
    for (const name of names) {
      if (presentValue(fields[name])) return fields[name];
    }
    return ".";
  }

  function hasPrefixedEvidence(fields, prefixes) {
    const source = fields || {};
    return Object.keys(source).some(key => {
      return prefixes.some(prefix => key.startsWith(prefix)) && presentValue(source[key]);
    });
  }

  function percentValue(value) {
    if (value === null || value === undefined || value === "") return ".";
    const numeric = Number(value);
    if (Number.isNaN(numeric)) return String(value);
    return (numeric * 100).toFixed(1) + "%";
  }

  function evidencePanel(title, rows) {
    const visibleRows = rows.filter(row => presentValue(row[1]));
    if (!visibleRows.length) return "";
    const body = visibleRows.map(([label, value]) => {
      return `<div class="evidence-row">
        <span>${escapeHtml(label)}</span>
        <span title="${escapeHtml(value)}">${escapeHtml(value)}</span>
      </div>`;
    }).join("");
    return `<div class="evidence-panel">
      <div class="evidence-title"><span>${escapeHtml(title)}</span></div>
      ${body}
    </div>`;
  }

  function renderEvidencePanels(q, graph, evt) {
    if (!els.evidencePanels) return;
    const inputs = data.inputs || {};
    const settings = data.settings || {};
    const fields = eventFields(evt);
    const panels = [];

    panels.push(evidencePanel("Alignment", [
      ["file", inputs.alignmentPath],
      ["format", inputs.alignmentFormat],
      ["selected", firstPresent(q && q.name, evt && evt.target)],
      ["best ref", q && q.bestRef],
      ["query cov", q && q.queryCov ? percentValue(q.queryCov) : "."]
    ]));

    if (inputs.gfa || graph || hasPrefixedEvidence(fields, ["graph_", "left_graph_", "right_graph_"])) {
      panels.push(evidencePanel("GFA", [
        ["file", inputs.gfa],
        ["node", firstPresent(graph && graph.graphNode, fields.graph_node, fields.left_graph_node)],
        ["right node", fields.right_graph_node],
        ["status", firstPresent(fields.graph_status, fields.graph_node_status, graph && graph.graphNodeStatus)],
        ["complexity", graph && graph.graphComplexity],
        ["path", fieldValue(fields, ["graph_path_nodes", "path_nodes"])],
        ["direct edge", fields.graph_direct_edge]
      ]));
    }

    if (inputs.readPaf || hasPrefixedEvidence(fields, ["longread_"])) {
      panels.push(evidencePanel("Long-read PAF", [
        ["file", inputs.readPaf],
        ["min MAPQ", settings.minReadMapq],
        ["bridges", fields.longread_bridge_reads],
        ["spanning", fields.longread_spanning_reads],
        ["split", fields.longread_split_reads],
        ["left edge", fields.longread_left_edge_reads],
        ["right edge", fields.longread_right_edge_reads],
        ["nearby", fields.longread_nearby_reads],
        ["orientation", fields.longread_orientation_summary],
        ["read order", fields.longread_read_order_summary],
        ["median gap", fields.longread_median_read_gap_bp]
      ]));
    }

    if (inputs.gaf || hasPrefixedEvidence(fields, ["gaf_"])) {
      panels.push(evidencePanel("Long-read GAF", [
        ["file", inputs.gaf],
        ["min MAPQ", settings.minGafMapq],
        ["status", fields.gaf_support_status],
        ["node", fields.gaf_node],
        ["node reads", fields.gaf_node_reads],
        ["traversals", fields.gaf_node_traversals],
        ["orientations", fields.gaf_node_orientations],
        ["path", fieldValue(fields, ["gaf_path_nodes", "gaf_path_examples"])],
        ["support", fields.gaf_path_support],
        ["alt path", fields.gaf_best_alt_path_nodes],
        ["alt support", fields.gaf_best_alt_support],
        ["reads", fields.gaf_selected_reads]
      ]));
    }

    const html = panels.filter(Boolean).join("");
    if (!html) {
      els.evidencePanels.classList.add("hidden");
      els.evidencePanels.innerHTML = "";
      return;
    }
    els.evidencePanels.classList.remove("hidden");
    els.evidencePanels.innerHTML = html;
  }

  function renderGraphPanel(q, graph) {
    if (!els.graphPanel) return;
    if (!graph) {
      els.graphPanel.classList.add("hidden");
      els.graphPanel.innerHTML = "";
      return;
    }
    els.graphPanel.classList.remove("hidden");
    if (graph.graphNodeStatus !== "present") {
      els.graphPanel.innerHTML = `<div class="graph-panel-title">
        Graph neighborhood <span class="badge red">missing</span>
      </div>
      <div class="graph-note">This contig was not found as an S node in the supplied GFA.</div>`;
      return;
    }
    const graphClass = graphBadgeClass(graph);
    const neighborRows = [
      ...(graph.incoming || []),
      ...(graph.outgoing || [])
    ].map(graphNeighborRow);
    const rows = neighborRows.length
      ? neighborRows.join("")
      : '<div class="graph-note">No immediate incoming or outgoing GFA neighbors.</div>';
    const refLabel = q && q.bestRef ? "best ref " + q.bestRef : "unaligned";
    els.graphPanel.innerHTML = `<div class="graph-panel-title">
      Graph neighborhood
      <span class="badge ${graphClass}">${escapeHtml(graph.graphComplexity || "unknown")}</span>
      <span class="badge">${escapeHtml(refLabel)}</span>
      <span class="badge">in/out ${graph.graphInDegree}/${graph.graphOutDegree}</span>
    </div>
    <div class="graph-neighbor-list">${rows}</div>`;
  }

  function renderStats() {
    const active = pieces.filter(p => !p.removed);
    const removed = pieces.length - active.length;
    const cutPieces = pieces.filter(p => {
      const q = queriesByName.get(p.source);
      return q && (p.start !== 1 || p.end !== q.length);
    }).length;
    const stats = [
      ["Active", active.length],
      ["Removed", removed],
      ["Pieces", pieces.length],
      ["Cut pieces", cutPieces],
      ["Aligned", data.stats.alignedContigs],
      ["Unaligned", data.stats.unalignedContigs]
    ];
    if (reviewEvents.length) {
      stats.push(["Events", reviewEvents.length]);
    }
    if (data.inputs.gfa) {
      stats.push(["Graph present", data.stats.graphContigsPresent || 0]);
      stats.push(["Graph branching", data.stats.graphBranchingContigs || 0]);
    }
    els.stats.innerHTML = stats.map(([label, value]) => `<div class="stat"><strong>${value}</strong><span>${label}</span></div>`).join("");
  }

  function eventBadgeClass(evt) {
    if (!evt.accept) return "amber";
    if (evt.status === "fillable" || evt.status === "candidate" || evt.status === "split") return "green";
    return "blue";
  }

  function renderEventList() {
    if (!els.eventPanel || !els.eventList) return;
    if (!reviewEvents.length) {
      els.eventPanel.classList.add("hidden");
      els.eventList.innerHTML = "";
      return;
    }
    els.eventPanel.classList.remove("hidden");
    els.eventTitle.textContent = (data.reviewTask || "review") + " events";
    els.eventList.innerHTML = reviewEvents.map(evt => {
      const selected = evt.event_id === selectedEventId ? " selected" : "";
      const fields = eventFields(evt);
      const support = firstPresent(
        fields.gaf_support_status,
        fields.gaf_path_support,
        fields.longread_bridge_reads,
        fields.longread_spanning_reads,
        fields.graph_status,
        ""
      );
      return `<div class="event-row${selected}" data-event-id="${escapeHtml(evt.event_id)}">
        <strong title="${escapeHtml(evt.target)}">${escapeHtml(evt.action)} | ${escapeHtml(evt.target)}</strong>
        <span>${escapeHtml(evt.status)} | ${escapeHtml(evt.reason || ".")}</span>
        <span><span class="badge ${eventBadgeClass(evt)}">${evt.accept ? "accepted" : "review"}</span> ${escapeHtml(support)}</span>
      </div>`;
    }).join("");
    els.eventList.querySelectorAll(".event-row").forEach(row => {
      row.addEventListener("click", () => {
        selectedEventId = row.dataset.eventId;
        const evt = selectedEvent();
        const pieceId = pieceIdForEvent(evt);
        if (pieceId) selectedId = pieceId;
        render();
      });
    });
  }

  function renderList() {
    const needle = els.search.value.trim().toLowerCase();
    const graphFilter = els.graphFilter ? els.graphFilter.value : "";
    const rows = pieces.filter(piece => {
      const q = queriesByName.get(piece.source) || {};
      const graph = q.graph || {};
      if (graphFilter && graphFilterValue(graph) !== graphFilter) return false;
      if (!needle) return true;
      return [
        piece.name,
        piece.source,
        piece.scaffold,
        piece.bestRef || "",
        graph.graphNode || "",
        graph.graphNodeStatus || "",
        graph.graphComplexity || "",
        graphDetail(graph)
      ]
        .some(v => String(v).toLowerCase().includes(needle));
    });
    els.pieceList.innerHTML = rows.map(piece => {
      const selected = piece.id === selectedId ? " selected" : "";
      const removed = piece.removed ? " removed" : "";
      return `<div class="piece-row${selected}${removed}" data-id="${escapeHtml(piece.id)}">
        <div>
          <div class="piece-name">${escapeHtml(piece.name)}</div>
          <div class="piece-meta">${pieceBadges(piece)}</div>
          <div class="piece-meta">${escapeHtml(piece.source)}:${piece.start}-${piece.end} | ${escapeHtml(piece.scaffold)}</div>
        </div>
        <div>${piece.removed ? "off" : "on"}</div>
      </div>`;
    }).join("");
    els.pieceList.querySelectorAll(".piece-row").forEach(row => {
      row.addEventListener("click", () => {
        selectedId = row.dataset.id;
        render();
      });
    });
  }

  function renderDetails() {
    const piece = selectedPiece();
    if (!piece) {
      els.selectedName.textContent = "No contig selected";
      els.selectedMeta.textContent = "";
      els.selectedGraph.textContent = "";
      renderGraphPanel(null, null);
      renderEvidencePanels(null, null, null);
      return;
    }
    const q = queriesByName.get(piece.source) || {};
    const graph = q.graph || null;
    const evt = selectedEvent();
    els.selectedName.textContent = piece.name;
    els.selectedMeta.textContent = [
      piece.source + ":" + piece.start + "-" + piece.end,
      fmtBp(pieceLength(piece)),
      q.bestRef ? "best " + q.bestRef + ":" + q.refStart + "-" + q.refEnd : "unaligned",
      q.queryCov ? "query cov " + (q.queryCov * 100).toFixed(1) + "%" : ""
    ].filter(Boolean).join(" | ");
    els.selectedGraph.textContent = [
      graphDetail(graph),
      evt ? "Event: " + evt.action + " | " + evt.status + " | " + evt.target : ""
    ].filter(Boolean).join(" | ");
    renderGraphPanel(q, graph);
    renderEvidencePanels(q, graph, evt);
    els.cutPos.value = "";
    els.cutPos.min = String(piece.start);
    els.cutPos.max = String(piece.end - 1);
    els.scaffoldName.value = piece.scaffold || "unplaced";
    document.getElementById("toggleRemove").textContent = piece.removed ? "Restore" : "Remove";
  }

  function render() {
    renderEventList();
    renderStats();
    renderList();
    renderDetails();
    drawPlot();
  }

  function escapeHtml(value) {
    return String(value).replace(/[&<>"']/g, ch => ({
      "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;"
    }[ch]));
  }

  function refOffsets() {
    const offsets = new Map();
    let cursor = 0;
    for (const rec of refs) {
      offsets.set(rec.name, cursor);
      cursor += rec.length;
    }
    return {offsets, total: cursor};
  }

  function resizeCanvas(canvas) {
    const rect = canvas.getBoundingClientRect();
    const ratio = window.devicePixelRatio || 1;
    const width = Math.max(320, Math.floor(rect.width * ratio));
    const height = Math.max(260, Math.floor(rect.height * ratio));
    if (canvas.width !== width || canvas.height !== height) {
      canvas.width = width;
      canvas.height = height;
    }
    return {width, height, ratio};
  }

  function drawPlot() {
    const canvas = els.dotplot;
    const ctx = canvas.getContext("2d");
    const {width, height, ratio} = resizeCanvas(canvas);
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.clearRect(0, 0, width, height);
    ctx.scale(ratio, ratio);

    const cssWidth = width / ratio;
    const cssHeight = height / ratio;
    ctx.fillStyle = "#ffffff";
    ctx.fillRect(0, 0, cssWidth, cssHeight);

    const piece = selectedPiece();
    if (!piece) return;
    const q = queriesByName.get(piece.source);
    if (!q) return;

    const left = 64;
    const right = 16;
    const top = 28;
    const bottom = 48;
    const plotW = Math.max(40, cssWidth - left - right);
    const plotH = Math.max(40, cssHeight - top - bottom);
    const {offsets, total} = refOffsets();
    const queryLen = q.length || piece.end;

    ctx.fillStyle = "#f8fafc";
    ctx.strokeStyle = "#94a3b8";
    ctx.lineWidth = 1;
    ctx.fillRect(left, top, plotW, plotH);
    ctx.strokeRect(left, top, plotW, plotH);

    ctx.font = "11px Arial";
    ctx.fillStyle = "#64748b";
    ctx.textAlign = "left";
    for (const rec of refs) {
      const x = left + (offsets.get(rec.name) / Math.max(1, total)) * plotW;
      ctx.strokeStyle = "#cbd5e1";
      ctx.beginPath();
      ctx.moveTo(x, top);
      ctx.lineTo(x, top + plotH);
      ctx.stroke();
      ctx.save();
      ctx.translate(x + 3, top + plotH + 16);
      ctx.rotate(-Math.PI / 7);
      ctx.fillText(rec.name, 0, 0);
      ctx.restore();
    }

    const y0 = top + ((piece.start - 1) / queryLen) * plotH;
    const y1 = top + (piece.end / queryLen) * plotH;
    ctx.fillStyle = "rgba(15, 118, 110, 0.08)";
    ctx.fillRect(left, y0, plotW, Math.max(1, y1 - y0));
    ctx.strokeStyle = "#0f766e";
    ctx.strokeRect(left, y0, plotW, Math.max(1, y1 - y0));

    for (const seg of segmentsByQuery.get(piece.source) || []) {
      const refOffset = offsets.get(seg.ref);
      if (refOffset === undefined) continue;
      const x1 = left + ((refOffset + seg.refStart - 1) / Math.max(1, total)) * plotW;
      const x2 = left + ((refOffset + seg.refEnd) / Math.max(1, total)) * plotW;
      const qStart = top + ((seg.queryStart - 1) / queryLen) * plotH;
      const qEnd = top + (seg.queryEnd / queryLen) * plotH;
      ctx.strokeStyle = seg.orientation === "-" ? "#dc2626" : "#2563eb";
      ctx.globalAlpha = 0.72;
      ctx.lineWidth = 1.4;
      ctx.beginPath();
      if (seg.orientation === "-") {
        ctx.moveTo(x1, qEnd);
        ctx.lineTo(x2, qStart);
      } else {
        ctx.moveTo(x1, qStart);
        ctx.lineTo(x2, qEnd);
      }
      ctx.stroke();
      ctx.globalAlpha = 1;
    }

    ctx.fillStyle = "#111827";
    ctx.font = "12px Arial";
    ctx.textAlign = "center";
    ctx.fillText("Reference position", left + plotW / 2, cssHeight - 8);
    ctx.save();
    ctx.translate(16, top + plotH / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText("Contig position", 0, 0);
    ctx.restore();
  }

  function selectedIndex() {
    return pieces.findIndex(p => p.id === selectedId);
  }

  function moveSelected(delta) {
    const idx = selectedIndex();
    if (idx < 0) return;
    const next = idx + delta;
    if (next < 0 || next >= pieces.length) return;
    const [piece] = pieces.splice(idx, 1);
    pieces.splice(next, 0, piece);
    render();
  }

  function toggleRemove() {
    const piece = selectedPiece();
    if (!piece) return;
    piece.removed = !piece.removed;
    render();
  }

  function invertSelected() {
    const piece = selectedPiece();
    if (!piece) return;
    piece.strand = piece.strand === "-" ? "+" : "-";
    render();
  }

  function sourcePieceCount(source) {
    return pieces.filter(p => p.source === source).length;
  }

  function addBreakpoint(pos) {
    const piece = selectedPiece();
    if (!piece) return;
    const cut = Number(pos);
    if (!Number.isInteger(cut) || cut < piece.start || cut >= piece.end) {
      setStatus("Breakpoint must fall inside the selected piece.");
      return;
    }
    const idx = selectedIndex();
    const count = sourcePieceCount(piece.source);
    const left = {
      ...piece,
      id: "p" + nextId++,
      name: piece.source + "_manual" + String(count).padStart(3, "0"),
      end: cut
    };
    const right = {
      ...piece,
      id: "p" + nextId++,
      name: piece.source + "_manual" + String(count + 1).padStart(3, "0"),
      start: cut + 1
    };
    pieces.splice(idx, 1, left, right);
    selectedId = right.id;
    setStatus("Added breakpoint after bp " + cut + ".");
    render();
  }

  function removeUnaligned() {
    for (const piece of pieces) {
      const q = queriesByName.get(piece.source);
      if (!q || !q.aligned) piece.removed = true;
    }
    render();
  }

  function restoreAll() {
    for (const piece of pieces) piece.removed = false;
    render();
  }

  function updateScaffoldName() {
    const piece = selectedPiece();
    if (!piece) return;
    piece.scaffold = els.scaffoldName.value.trim() || "unplaced";
    render();
  }

  function recipe() {
    return {
      schema: data.schema,
      chromosortVersion: data.version,
      sourceAssembly: data.inputs.assemblyFasta,
      scaffoldEnabled,
      gapBp: Math.max(0, Number(els.gapBp.value) || 0),
      pieces: pieces.map(p => ({
        id: p.id,
        source: p.source,
        name: p.name,
        start: p.start,
        end: p.end,
        strand: p.strand,
        scaffold: p.scaffold || "unplaced",
        removed: Boolean(p.removed)
      }))
    };
  }

  function downloadText(filename, text, type) {
    const blob = new Blob([text], {type});
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
  }

  function parseFasta(text) {
    const parsed = {};
    let name = null;
    let parts = [];
    for (const raw of text.split(/\r?\n/)) {
      if (!raw) continue;
      if (raw.startsWith(">")) {
        if (name) parsed[name] = parts.join("");
        name = raw.slice(1).trim().split(/\s+/)[0];
        parts = [];
      } else {
        parts.push(raw.trim());
      }
    }
    if (name) parsed[name] = parts.join("");
    return parsed;
  }

  function reverseComplement(seq) {
    const map = {
      A: "T", C: "G", G: "C", T: "A", R: "Y", Y: "R", K: "M", M: "K",
      S: "S", W: "W", B: "V", D: "H", H: "D", V: "B", N: "N",
      a: "t", c: "g", g: "c", t: "a", r: "y", y: "r", k: "m", m: "k",
      s: "s", w: "w", b: "v", d: "h", h: "d", v: "b", n: "n"
    };
    let out = "";
    for (let i = seq.length - 1; i >= 0; i--) out += map[seq[i]] || "N";
    return out;
  }

  function wrapSeq(seq) {
    const rows = [];
    for (let i = 0; i < seq.length; i += 80) rows.push(seq.slice(i, i + 80));
    return rows.join("\n");
  }

  function activePieces() {
    return pieces.filter(p => !p.removed);
  }

  function missingSources(active) {
    const needed = new Set(active.map(p => p.source));
    return [...needed].filter(name => !sequences[name]);
  }

  function pieceSeq(piece) {
    const seq = sequences[piece.source];
    let out = seq.slice(piece.start - 1, piece.end);
    if (piece.strand === "-") out = reverseComplement(out);
    return out;
  }

  function uniqueNames(active) {
    const counts = new Map();
    const names = new Map();
    for (const piece of active) {
      const base = piece.name || piece.source;
      const next = (counts.get(base) || 0) + 1;
      counts.set(base, next);
      names.set(piece.id, next === 1 ? base : base + "_" + next);
    }
    return names;
  }

  function buildPieceFasta(active) {
    const names = uniqueNames(active);
    return active.map(piece => {
      const name = names.get(piece.id);
      const header = `>${name} original=${piece.source} slice=${piece.start}-${piece.end} strand=${piece.strand} scaffold=${piece.scaffold || "unplaced"}`;
      return header + "\n" + wrapSeq(pieceSeq(piece));
    }).join("\n") + "\n";
  }

  function buildScaffoldFasta(active) {
    const gap = "N".repeat(Math.max(0, Number(els.gapBp.value) || 0));
    const groups = new Map();
    for (const piece of active) {
      const scaffold = piece.scaffold || "unplaced";
      if (!groups.has(scaffold)) groups.set(scaffold, []);
      groups.get(scaffold).push(piece);
    }
    const records = [];
    for (const [scaffold, members] of groups.entries()) {
      const seq = members.map(pieceSeq).join(gap);
      records.push(`>${scaffold} manual_pieces=${members.length} gap_bp=${gap.length}\n${wrapSeq(seq)}`);
    }
    return records.join("\n") + "\n";
  }

  function exportFasta() {
    const active = activePieces();
    const missing = missingSources(active);
    if (missing.length) {
      setStatus("Load the assembly FASTA before exporting. Missing: " + missing.slice(0, 5).join(", "));
      return;
    }
    const fasta = scaffoldEnabled ? buildScaffoldFasta(active) : buildPieceFasta(active);
    downloadText(data.settings.suggestedOutputFasta || "chromosort.manual.fa", fasta, "text/plain");
    setStatus("Exported FASTA with " + active.length + " active piece(s).");
  }

  function loadFastaFile(file) {
    if (!file) return;
    const reader = new FileReader();
    reader.onload = () => {
      sequences = {...sequences, ...parseFasta(String(reader.result || ""))};
      setStatus("Loaded " + Object.keys(sequences).length + " FASTA sequence(s).");
    };
    reader.readAsText(file);
  }

  function loadRecipeFile(file) {
    if (!file) return;
    const reader = new FileReader();
    reader.onload = () => {
      try {
        const incoming = JSON.parse(String(reader.result || "{}"));
        if (incoming.schema !== data.schema || !Array.isArray(incoming.pieces)) {
          throw new Error("Not a ChromoSort manual recipe.");
        }
        pieces = incoming.pieces.map(p => ({...p}));
        scaffoldEnabled = Boolean(incoming.scaffoldEnabled);
        els.scaffoldMode.checked = scaffoldEnabled;
        els.gapBp.value = String(incoming.gapBp ?? 100);
        selectedId = pieces.length ? pieces[0].id : null;
        nextId = pieces.length + 1;
        setStatus("Loaded manual recipe.");
        render();
      } catch (err) {
        setStatus(err.message || String(err));
      }
    };
    reader.readAsText(file);
  }

  els.datasetLabel.textContent = `${data.inputs.assemblyFasta} against ${data.inputs.refFasta}`;
  if (!data.inputs.gfa && els.graphFilter) els.graphFilter.disabled = true;
  els.search.addEventListener("input", renderList);
  els.graphFilter.addEventListener("change", renderList);
  document.getElementById("moveUp").addEventListener("click", () => moveSelected(-1));
  document.getElementById("moveDown").addEventListener("click", () => moveSelected(1));
  document.getElementById("invert").addEventListener("click", invertSelected);
  document.getElementById("toggleRemove").addEventListener("click", toggleRemove);
  document.getElementById("removeUnaligned").addEventListener("click", removeUnaligned);
  document.getElementById("restoreAll").addEventListener("click", restoreAll);
  document.getElementById("addCut").addEventListener("click", () => addBreakpoint(els.cutPos.value));
  els.scaffoldName.addEventListener("change", updateScaffoldName);
  els.scaffoldMode.addEventListener("change", () => {
    scaffoldEnabled = els.scaffoldMode.checked;
    setStatus(scaffoldEnabled ? "Scaffold export enabled." : "Piece export enabled.");
  });
  els.fastaFile.addEventListener("change", () => loadFastaFile(els.fastaFile.files[0]));
  els.recipeFile.addEventListener("change", () => loadRecipeFile(els.recipeFile.files[0]));
  document.getElementById("exportRecipe").addEventListener("click", () => {
    downloadText("chromosort.manual.recipe.json", JSON.stringify(recipe(), null, 2) + "\n", "application/json");
  });
  document.getElementById("exportFasta").addEventListener("click", exportFasta);
  els.dotplot.addEventListener("click", evt => {
    const piece = selectedPiece();
    const q = piece ? queriesByName.get(piece.source) : null;
    if (!piece || !q) return;
    const rect = els.dotplot.getBoundingClientRect();
    const top = 28;
    const bottom = 48;
    const plotH = Math.max(40, rect.height - top - bottom);
    const y = evt.clientY - rect.top;
    const bp = Math.round(((y - top) / plotH) * q.length);
    if (bp >= piece.start && bp < piece.end) {
      els.cutPos.value = String(bp);
      setStatus("Breakpoint staged after bp " + bp + ".");
    }
  });
  window.addEventListener("resize", drawPlot);

  if (Object.keys(sequences).length) {
    setStatus("Sequences are embedded; FASTA export is ready.");
  } else {
    setStatus("Load the assembly FASTA to export edited sequence.");
  }
  render();
})();
</script>
</body>
</html>
"""


if __name__ == "__main__":
    main()
