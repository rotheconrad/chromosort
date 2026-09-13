"""Bounded, stage-bound inspection proposals; never a sequence-edit planner.

The correction planner remains authoritative for automatic output. This module
exposes frame-local orientation evidence and competing partitions for review.
"""

import html
import json
from collections import defaultdict
from dataclasses import asdict
from types import SimpleNamespace

from . import fix_contigs as fix
from .continuity import assess_plans
from .graph import (build_direct_projection, build_path_projection,
                    graph_coordinate_evidence, read_gfa)
from .multireference import write_tsv
from .provenance import sha256_file, stable_id, write_json
from .reference_order import merge_intervals


SCHEMA = "chromosort-review-proposals-v1"
ORIENTATION_LIMIT = 25
EVENT_COLUMNS = ["event_id", "kind", "target", "start0", "end0", "positions", "action",
                 "decision", "assembly_sha256", "reference_frames", "supporting_block_ids",
                 "opposing_block_ids", "overlapping_block_ids", "loss_stages", "reasons",
                 "objective", "coordinate_evidence", "figures", "implicit_edits"]


def add_arguments(parser):
    parser.add_argument("--inspect-orientation", action="store_true",
                        help="Review aligned internal orientation discrepancies before correction size filtering/smoothing; never reverse automatically.")
    parser.add_argument("--inspect-min-orientation-bp", type=int, default=1000,
                        help="Independent review interval floor, in native bases.")
    parser.add_argument("--alternative-plans", type=int, default=0,
                        help="Up to 20 review-only competing segmentation plans per contig; 0 disables.")
    parser.add_argument("--review-max-blocks", type=int, default=128,
                        help="Skip new inspection paths above this per-contig block cap (1–256); never truncate planner inputs.")


def validate_arguments(args):
    if args.inspect_min_orientation_bp < 1 or not 0 <= args.alternative_plans <= 20 or not 1 <= args.review_max_blocks <= 256:
        raise ValueError("Review requires a positive orientation floor, 0–20 alternative plans and 1–256 blocks")


def segment_key(segment):
    return tuple(asdict(segment).values())


def overlap(start, end, block):
    return max(0, min(end, block["end0"]) - max(start, block["start0"]))


def block_registry(bundle, args):
    """Coalesce identical mapping rows for inspection, preserving multiplicity.

    This does not change the additive weights used by automatic correction or
    the alternative objective. Origin labels are provenance, not independence.
    """
    eligible = {segment_key(s) for s in bundle.correction_segments(min_mapq=args.min_mapq)}
    registry = {}
    roles = {r["reference_id"]: r.get("role", "supporting") for r in bundle.references}
    for source, segment in bundle.iter_segments(min_mapq=args.min_mapq, with_evidence=True):
        if args.contigs and segment.query not in args.contigs:
            continue
        seq = bundle.sequences[segment.ref]
        key = segment_key(segment)
        bid = stable_id("alignment", bundle.assembly_digest, bundle.reference_digests[seq.reference_id], key)
        if bid not in registry:
            start, end = fix.query_interval(segment)
            rstart, rend = fix.ref_interval(segment)
            losses = ([] if seq.backbone else ["not_backbone"])
            if seq.backbone and key not in eligible:
                losses.append("withheld_cross_reference_overlap")
            if segment.len_query < args.min_segment_bp:
                losses.append("below_correction_segment_floor")
            registry[bid] = {
                "block_id": bid, "target": segment.query, "start0": start, "end0": end,
                "reference_start0": rstart, "reference_end0": rend, "orientation": segment.orientation,
                "frame": segment.ref, "reference": {**asdict(seq), "role": roles[seq.reference_id],
                    "sha256": bundle.reference_digests[seq.reference_id]},
                "identity_percent": segment.identity, "mapq": segment.mapq,
                "aligned_query_bp": segment.len_query, "backbone": seq.backbone,
                "correction_eligible": key in eligible, "loss_stages": losses,
                "row_count": 0, "sources": [],
            }
        block = registry[bid]
        block["row_count"] += 1
        provenance = {k: source.get(k) for k in ["evidence_id", "kind", "sha256", "reference_id",
            "reference_sha256", "assembly_sha256", "dependency_group", "origin_id", "evidence_class"]}
        if provenance not in block["sources"]:
            block["sources"].append(provenance)
    return sorted(registry.values(), key=lambda b: (b["target"], b["start0"], b["end0"], b["frame"], b["block_id"]))


def objective(groups, penalty):
    discordant = sum(g.discordant_bp for g in groups)
    penalty_cost = penalty * max(0, len(groups) - 1)
    return {"discordant_bp": discordant, "penalty_bp": penalty_cost,
            "total": discordant + penalty_cost, "breakpoints": max(0, len(groups) - 1)}


def group_records(groups, include_orientation, seq_len, args, raw):
    result = []
    for group in groups:
        summary = group.summary
        member_ids = sorted({b["block_id"] for b in raw if b["correction_eligible"]
            and "below_correction_segment_floor" not in b["loss_stages"]
            and any(b["frame"] == m.ref and b["orientation"] == m.orientation
                    and m.query_start <= b["start0"] < b["end0"] <= m.query_end
                    and m.ref_start <= b["reference_start0"] < b["reference_end0"] <= m.ref_end
                    for m in group.blocks)})
        supporting = sorted(b["block_id"] for b in raw if b["block_id"] in member_ids
            and b["frame"] == summary.ref and (not include_orientation or b["orientation"] == summary.orientation))
        result.append({"start0": summary.query_start, "end0": summary.query_end,
            "signature": list(fix.auto_signature(summary, include_orientation)),
            "block_ids": member_ids, "merged_block_count": len(group.blocks),
            "supporting_block_ids": supporting, "opposing_block_ids": sorted(set(member_ids) - set(supporting)),
            "discordant_bp": group.discordant_bp, "dominant_weighted_bp": group.dominant_weighted_bp,
            "total_weighted_bp": group.total_weighted_bp, "dominant_aligned_bp": summary.aligned_bp,
            "query_span_fraction": summary.query_span / seq_len,
            "passes_group_support": fix.auto_group_piece_is_supported(group, seq_len, args)})
    return result


def partition_trace(blocks, seq_len, args, raw):
    if not blocks:
        return {"include_orientation": False, "best_partition": None, "after_terminal_merge": []}
    include = fix.auto_include_orientation(blocks, args)
    groups = fix.segment_blocks_for_auto(blocks, args, include)
    best = {**objective(groups, args.breakpoint_penalty_bp),
            "groups": group_records(groups, include, seq_len, args, raw)}
    merged = fix.merge_adjacent_auto_groups_by_signature(groups, include)
    terminal, changed = fix.merge_unsupported_terminal_groups(merged, seq_len, args, include)
    return {"include_orientation": include, "best_partition": best,
        "after_terminal_merge": group_records(terminal, include, seq_len, args, raw),
        "terminal_groups_merged": changed, "objective_comparable_to_automatic": args.mode != "sensitive",
        "comparison_note": "Sensitive automatic mode does not optimize this objective." if args.mode == "sensitive"
                           else "Objective precedes terminal merging, support and breakpoint guards, and read continuity."}


def new_event(bundle, target, kind, start, end, positions):
    return {"event_id": stable_id("inspection", bundle.assembly_digest, target, kind, positions),
        "kind": kind, "target": target, "start0": start, "end0": end, "positions": positions,
        "action": "inspect", "decision": "pending", "assembly_sha256": bundle.assembly_digest,
        "supporting_block_ids": [], "opposing_block_ids": [], "overlapping_block_ids": [],
        "reference_frames": [], "frame_contexts": [], "loss_stages": [], "reasons": [],
        "objective": None, "coordinate_evidence": [], "figures": [], "implicit_edits": False}


def orientation_events(bundle, target, raw, trace, args):
    frames = defaultdict(list)
    for b in raw:
        if b["backbone"]:
            frames[b["frame"]].append(b)
    events, contexts, below_floor = {}, [], 0
    for frame, blocks in sorted(frames.items()):
        first, last = min(b["start0"] for b in blocks), max(b["end0"] for b in blocks)
        left = {b["orientation"] for b in blocks if b["start0"] == first}
        right = {b["orientation"] for b in blocks if b["end0"] == last}
        anchor = next(iter(left)) if len(left) == 1 and left == right else None
        context = {"frame": frame, "anchor_orientation": anchor,
                   "mapped_start0": first, "mapped_end0": last,
                   "status": "ambiguous_terminal_orientation" if anchor is None else "internal_orientation_context"}
        if len({b["orientation"] for b in blocks}) == 1:
            context["status"] = "whole_record_orientation" if anchor == "-" else "uniform_forward_orientation"
        contexts.append(context)
        if anchor is None or context["status"] != "internal_orientation_context":
            continue
        contrary = [b for b in blocks if b["orientation"] != anchor]
        anchors = [b for b in blocks if b["orientation"] == anchor]
        for start, end in merge_intervals([(b["start0"], b["end0"]) for b in contrary]):
            if end - start < args.inspect_min_orientation_bp:
                below_floor += 1
                continue
            left_blocks = [b for b in anchors if b["start0"] < start]
            right_blocks = [b for b in anchors if b["end0"] > end]
            if not left_blocks or not right_blocks:
                continue
            support = [b for b in contrary if overlap(start, end, b)]
            opposing = {b["block_id"] for b in left_blocks + right_blocks}
            overlapping = [b for b in raw if overlap(start, end, b)]
            event = events.setdefault((start, end), new_event(bundle, target, "orientation_discrepancy", start, end, [start, end]))
            losses = {stage for b in support for stage in b["loss_stages"]}
            retained = [b for b in support if b["correction_eligible"] and "below_correction_segment_floor" not in b["loss_stages"]]
            # Sensitive correction bypasses objective smoothing and its group
            # support/terminal merges. A diagnostic DP cannot explain its loss.
            if retained and args.mode != "sensitive":
                if not trace["include_orientation"]:
                    losses.add("orientation_not_used_in_mode")
                else:
                    best = trace.get("best_partition")
                    for b in retained:
                        matching = [g for g in (best or {}).get("groups", []) if b["block_id"] in g["block_ids"]]
                        if matching and all(g["signature"] != [frame, b["orientation"]] for g in matching):
                            losses.add("smoothed_by_objective")
                        elif not any(g["signature"] == [frame, b["orientation"]]
                                     and b["block_id"] in g["block_ids"] for g in trace["after_terminal_merge"]):
                            losses.add("merged_with_unsupported_terminal_group")
                        elif any(not g["passes_group_support"] for g in trace["after_terminal_merge"]):
                            losses.add("fails_group_support")
            if any(b["frame"] != frame for b in overlapping):
                losses.add("cross_reference_context_requires_review")
            own_overlaps = [b for b in anchors if overlap(start, end, b)]
            if own_overlaps:
                losses.add("overlapping_opposite_orientation_evidence")
            event["supporting_block_ids"].extend(b["block_id"] for b in support)
            event["opposing_block_ids"].extend(opposing)
            event["overlapping_block_ids"].extend(b["block_id"] for b in overlapping)
            event["reference_frames"].append(frame)
            event["loss_stages"].extend(losses)
            event["frame_contexts"].append({**context, "interval_orientation": support[0]["orientation"],
                "supporting_block_ids": sorted(b["block_id"] for b in support),
                "left_flank_block_ids": sorted(b["block_id"] for b in left_blocks),
                "right_flank_block_ids": sorted(b["block_id"] for b in right_blocks),
                "opposing_overlap_bp": sum(e - s for s, e in merge_intervals([
                    (max(start, b["start0"]), min(end, b["end0"])) for b in own_overlaps])),
                "supporting_union_bp": end - start,
                "supporting_row_count": sum(b["row_count"] for b in support),
                "supporting_additive_bp": sum((b["end0"] - b["start0"]) * b["row_count"] for b in support)})
    ordered = [events[key] for key in sorted(events)]
    for event in ordered:
        for key in ["supporting_block_ids", "opposing_block_ids", "overlapping_block_ids", "reference_frames", "loss_stages"]:
            event[key] = sorted(set(event[key]))
        event["reasons"] = [f"Aligned internal orientation discrepancy [{event['start0']}, {event['end0']}) in {len(event['reference_frames'])} reference frame(s); inspect both native edges. Biological variation and alignment ambiguity remain unresolved.",
            "Mapped evidence is covered by the alignment union, so this is distinct from unaligned-gap inspection."]
    return ordered[:ORIENTATION_LIMIT], {"frames": contexts, "eligible_events": len(ordered),
        "frame_intervals_below_review_floor": below_floor,
        "omitted_by_event_cap": max(0, len(ordered) - ORIENTATION_LIMIT)}


def kbest_partitions(blocks, args, include_orientation, limit):
    """K shortest paths through the interval DAG, with deterministic ties.

    Interval costs call the same additive function as automatic correction.
    Bound paths per prefix before canonical group merging; uniqueness filtering
    can therefore yield fewer than the requested number of displayed plans.
    """
    costs = {(start, end): fix.group_support(blocks[start:end], include_orientation)["discordant_bp"]
             for end in range(1, len(blocks) + 1) for start in range(end)}
    paths = [[(0.0, ())]]
    for end in range(1, len(blocks) + 1):
        candidates = []
        for start in range(end):
            for cost, boundaries in paths[start]:
                cuts = boundaries + ((start,) if start else ())
                candidates.append((cost + costs[start, end] + (args.breakpoint_penalty_bp if start else 0), cuts))
        paths.append(sorted(candidates, key=lambda p: (p[0], len(p[1]), p[1]))[:limit])
    return paths[-1]


def alternative_events(bundle, target, blocks, raw, trace, automatic_positions, args):
    if not blocks:
        return [], {"status": "no_correction_blocks", "raw_paths_examined": 0}
    include = trace["include_orientation"]
    paths = kbest_partitions(blocks, args, include, 4 * (args.alternative_plans + 1))
    unique, invalid = {}, 0
    length = bundle.query_by_name[target].length
    for _, boundaries in paths:
        edges = [0, *boundaries, len(blocks)]
        groups = [fix.summarize_group(blocks[s:e], include) for s, e in zip(edges, edges[1:])]
        groups = fix.merge_adjacent_auto_groups_by_signature(groups, include)
        positions = tuple(fix.boundary_between(a.blocks[-1], b.blocks[0]) for a, b in zip(groups, groups[1:]))
        if not positions or list(positions) == automatic_positions:
            continue
        if not all(0 < p < length for p in positions) or list(positions) != sorted(set(positions)):
            invalid += 1
            continue
        score = objective(groups, args.breakpoint_penalty_bp)
        signature = tuple(fix.auto_signature(g.summary, include) for g in groups)
        key = (score["total"], len(positions), positions, signature, boundaries)
        if positions not in unique or key < unique[positions][0]:
            unique[positions] = (key, groups, score)
    events = []
    for rank, (positions, (_, groups, score)) in enumerate(sorted(unique.items(), key=lambda p: p[1][0])[:args.alternative_plans], 1):
        event = new_event(bundle, target, "alternative_segmentation", positions[0], positions[-1], list(positions))
        records = group_records(groups, include, length, args, raw)
        event["objective"] = {**score, "delta_to_best_partition": score["total"] - trace["best_partition"]["total"],
            "rank": rank, "groups": records, "include_orientation": include,
            "objective_comparable_to_automatic": trace["objective_comparable_to_automatic"]}
        event["supporting_block_ids"] = sorted({bid for g in records for bid in g["supporting_block_ids"]})
        event["reference_frames"] = sorted({g["signature"][0] for g in records})
        event["opposing_block_ids"] = sorted({b["block_id"] for b in raw if b["block_id"] not in event["supporting_block_ids"]})
        event["overlapping_block_ids"] = sorted({b["block_id"] for b in raw if any(b["start0"] < p < b["end0"] for p in positions)})
        event["loss_stages"] = ["alternative_partition_for_review"]
        if score["total"] > trace["best_partition"]["total"]:
            event["loss_stages"].append("higher_objective_than_best_partition")
        if not all(g["passes_group_support"] for g in records):
            event["loss_stages"].append("fails_group_support")
        if any(e - s < args.min_piece_bp for s, e in zip([0, *positions], [*positions, length])):
            event["loss_stages"].append("fails_minimum_emitted_piece")
        event["reasons"] = [f"Alternative segmentation at native boundaries {list(positions)}; objective {score['total']:.6g}, delta {event['objective']['delta_to_best_partition']:.6g} to the best pre-guard partition. These units are identity-weighted aligned bases, not confidence.",
            "Bounded search precedes terminal merging and automatic guards. Select each desired cut explicitly; acknowledging this event makes no edit."]
        if args.mode == "sensitive":
            event["reasons"].append("Sensitive automatic correction bypasses this diagnostic objective and smoothing support checks; objective differences do not explain its automatic choices.")
        events.append(event)
    return events, {"status": "evaluated", "raw_paths_examined": len(paths),
        "unique_alternatives": len(unique), "invalid_native_partitions": invalid,
        "emitted": len(events), "paths_per_prefix": 4 * (args.alternative_plans + 1)}


def coordinate_evidence(events, bundle, args):
    positions = defaultdict(set)
    for event in events:
        positions[event["target"]].update(event["positions"])
    plans = {target: SimpleNamespace(status="split", pieces=[SimpleNamespace(slice_end=p)
             for p in sorted(points)] + [SimpleNamespace(slice_end=bundle.query_by_name[target].length)])
             for target, points in positions.items()}
    report_args = SimpleNamespace(**{**vars(args), "read_continuity_policy": "report"})
    reads = {(r["contig"], r["position"]): r for r in assess_plans(plans, bundle, report_args)}
    graphs = [e for e in bundle.evidence if e["kind"] == "gfa"]
    graph, projections, warnings = None, [], []
    if len(graphs) == 1:
        graph = read_gfa(bundle.resolve(graphs[0]["path"]))
        # Resolve each contig independently: a path for one must not prevent a
        # direct projection for another. Require exact total coordinate length.
        for target in sorted(positions):
            projected = build_path_projection(graph, path_names=[target], warnings=warnings)
            if not projected:
                projected = build_direct_projection(graph, [target], warnings=warnings)
            if projected and max(p.contig_end_0 for p in projected) == bundle.query_by_name[target].length:
                projections.extend(projected)
    for event in events:
        for point in event["positions"]:
            # A boundary after base p is between 1-based bases p and p+1. Show
            # both projections; never borrow graph context from a nearby cut.
            sides = {side: asdict(graph_coordinate_evidence(graph, projections, event["target"], base))
                     for side, base in [("left_base", point), ("right_base", point + 1)]}
            available = all(v["status"] == "present" for v in sides.values())
            graph_record = {"status": "coordinate_context" if available else
                "unavailable_no_graph" if not graphs else
                "unavailable_multiple_graph_sources" if len(graphs) > 1 else "unavailable_exact_projection",
                "evidence_id": graphs[0]["evidence_id"] if len(graphs) == 1 else None,
                "evidence_sha256": graphs[0]["sha256"] if len(graphs) == 1 else None,
                "assembly_sha256": bundle.assembly_digest, "sides": sides,
                "interpretation": "Graph coordinate context does not establish error or continuity."}
            event["coordinate_evidence"].append({"position": point,
                "reads": reads[event["target"], point], "graph": graph_record})
    return [asdict(w) for w in warnings]


def build(bundle, args, planner_args, automatic_rows):
    raw = block_registry(bundle, args)
    by_target = defaultdict(list)
    for block in raw:
        by_target[block["target"]].append(block)
    merged = fix.collect_blocks(bundle.path, "manifest", args.contigs or None, args.min_segment_bp,
                               planner_args.min_segment_idy, planner_args.max_merge_gap, min_mapq=args.min_mapq)
    cuts, automatic = defaultdict(list), defaultdict(list)
    for row in automatic_rows:
        automatic[row["target"]].append({k: row[k] for k in ["action", "position", "reason"]})
        if row["action"] == "cut":
            cuts[row["target"]].append(row["position"])
    result = {"schema": SCHEMA, "assembly_sha256": bundle.assembly_digest,
        "manifest_sha256": sha256_file(bundle.path), "inputs": bundle.input_records(),
        "parameters": {k: getattr(args, k) for k in ["inspect_orientation", "inspect_min_orientation_bp",
            "alternative_plans", "review_max_blocks", "mode", "min_mapq", "min_segment_bp"]},
        "objective_definition": {"signature": "reference frame, plus orientation when enabled by mode",
            "weight": "identity_percent * aligned_query_bp / 100; additive, including redundant overlaps",
            "discordance": "sum signature weights minus maximum signature weight within each group",
            "total": "sum group discordance + breakpoint_penalty_bp * (group_count - 1)",
            "breakpoint_penalty_bp": planner_args.breakpoint_penalty_bp, "calibrated_confidence": False},
        "limits": {"orientation_events_per_contig": ORIENTATION_LIMIT, "max_blocks": args.review_max_blocks,
            "policy": "Skip an inspection path above its block cap; orientation events sorted by native start/end. No automatic inputs are truncated."},
        "blocks": raw, "contigs": [], "events": []}
    for rec in bundle.query_records:
        if args.contigs and rec.name not in args.contigs:
            continue
        blocks, mappings = merged.get(rec.name, []), by_target[rec.name]
        trace = partition_trace(blocks, rec.length, planner_args, mappings) if len(blocks) <= args.review_max_blocks else {
            "include_orientation": False, "best_partition": None, "after_terminal_merge": [], "status": "skipped_block_cap"}
        ledger = {"target": rec.name, "canonical_alignment_blocks": len(mappings), "merged_correction_blocks": len(blocks),
            "automatic_scan_proposals": automatic[rec.name], "automatic_positions": sorted(cuts[rec.name]), "trace": trace}
        if args.inspect_orientation:
            if len(mappings) > args.review_max_blocks or len(blocks) > args.review_max_blocks:
                ledger["orientation"] = {"status": "skipped_block_cap"}
            else:
                events, status = orientation_events(bundle, rec.name, mappings, trace, args)
                ledger["orientation"] = {"status": "evaluated", **status}
                for event in events:
                    event["automatic_positions"] = sorted(cuts[rec.name])
                    event["loss_stages"].append("represented_at_both_automatic_boundaries" if all(
                        p in cuts[rec.name] for p in event["positions"]) else "not_represented_at_both_automatic_boundaries")
                result["events"].extend(events)
        if args.alternative_plans:
            if len(blocks) > args.review_max_blocks:
                ledger["alternatives"] = {"status": "skipped_block_cap"}
            else:
                settings = SimpleNamespace(**{**vars(planner_args), "alternative_plans": args.alternative_plans})
                events, status = alternative_events(bundle, rec.name, blocks, mappings, trace, sorted(cuts[rec.name]), settings)
                ledger["alternatives"] = status
                result["events"].extend(events)
        result["contigs"].append(ledger)
    result["graph_projection_warnings"] = coordinate_evidence(result["events"], bundle, args)
    return result


def write_report(out, record):
    write_json(out / "review_proposals.json", record)
    flat = [{k: json.dumps(e[k], sort_keys=True) if isinstance(e[k], (list, dict)) or e[k] is None else e[k]
             for k in EVENT_COLUMNS} for e in record["events"]]
    write_tsv(out / "review_proposals.tsv", flat, EVENT_COLUMNS)
    escape = lambda value: html.escape(str(value), quote=True)
    sections = []
    for e in record["events"]:
        links = " · ".join(f'<a href="{escape(f["path"])}.svg">Boundary {f["position"]}</a>' for f in e["figures"])
        evidence_rows = "".join(f'<tr><td>{c["position"]}</td><td>{escape(c["reads"]["continuity_status"])}</td>'
            f'<td>{c["reads"]["spanning_molecules"]}</td><td>{escape(c["reads"]["evidence_id"])}</td>'
            f'<td>{c["reads"]["anchor_bp"]} bp / MAPQ {c["reads"]["min_mapq"]}</td>'
            f'<td>{escape(c["graph"]["status"])}</td></tr>' for c in e["coordinate_evidence"])
        details = {k: e[k] for k in ["reference_frames", "frame_contexts", "loss_stages", "objective", "coordinate_evidence"]}
        selected = set(e["supporting_block_ids"] + e["opposing_block_ids"] + e["overlapping_block_ids"])
        block_rows = "".join("<tr>" + "".join(f"<td>{escape(b[k])}</td>" for k in
            ["block_id", "frame", "start0", "end0", "orientation", "row_count", "loss_stages"]) + "</tr>"
            for b in record["blocks"] if b["block_id"] in selected)
        sections.append(f'<article id="{escape(e["event_id"])}"><h2>{escape(e["kind"])} · {escape(e["target"])}</h2>'
            f'<p><code>{escape(e["event_id"])}</code> · pending · native positions {escape(e["positions"])}</p>'
            + "".join(f"<p>{escape(reason)}</p>" for reason in e["reasons"]) + f"<p>{links}</p>"
            + '<table><tr><th>Native boundary</th><th>Read evidence</th><th>Spanning read names</th><th>Source</th><th>Anchor / filter</th><th>Graph evidence</th></tr>'
            + evidence_rows + '</table>'
            + "<details><summary>Objective, frame context and coordinate evidence</summary><pre>"
            + escape(json.dumps(details, indent=2)) + "</pre></details><details><summary>Supporting, opposing and overlapping alignments</summary>"
            + "<table><tr><th>Block</th><th>Frame</th><th>Start</th><th>End</th><th>Strand</th><th>Rows</th><th>Correction stages</th></tr>"
            + block_rows + "</table></details></article>")
    (out / "review-proposals.html").write_text('<!doctype html><html lang="en"><meta charset="utf-8">'
        '<title>ChromoSort review proposals</title><style>body{font:15px system-ui;color:#243d43;max-width:1200px;margin:32px auto;padding:0 20px}article{border-top:1px solid #ccd;padding:16px 0}table{border-collapse:collapse;width:100%}td,th{padding:8px;border-bottom:1px solid #dde;text-align:left}pre{white-space:pre-wrap;overflow-wrap:anywhere}code,td{overflow-wrap:anywhere}summary{cursor:pointer;padding:10px}a{color:#176958}</style>'
        '<h1>Review proposals</h1><p>These events expose evidence for inspection. An aligned inversion can be genuine variation. '
        'Acknowledging an inspection changes no sequence. Use explicit reviewed interval edits or cuts to act on selected coordinates.</p>'
        '<p><a href="index.html">Decisions and read panels</a> · <a href="review_proposals.json">Full evidence JSON</a> · '
        '<a href="review_proposals.tsv">Event TSV</a></p><p>Native intervals are 0-based and half-open; boundary p lies between bases p and p+1 (1-based). '
        'Objective differences are additive aligned-base costs, not probabilities. Missing evidence does not confirm an error.</p>'
        + "".join(sections) + '<details><summary>Per-contig limits, exclusions and automatic-plan context</summary><pre>'
        + escape(json.dumps(record["contigs"], indent=2)) + '</pre></details></html>')
