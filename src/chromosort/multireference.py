"""Copy-preserving assignment and ordering over labeled reference genomes."""

import csv
from collections import defaultdict
from dataclasses import asdict, replace
from pathlib import Path

from . import reference_order as order
from .agp import write_agp, write_component_tsv
from .manifest import InputBundle, namespaced
from .provenance import audit_record, stable_id, write_json, sha256_file


def write_tsv(path, rows, columns=None):
    rows = list(rows)
    columns = columns or list(rows[0] if rows else [])
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def union_identity_score(segments):
    """Integral of maximum identity per aligned query base (no duplicated votes)."""
    events = defaultdict(list)
    for index, seg in enumerate(segments):
        start, end = order.closed_to_half_open(seg.query_start, seg.query_end)
        events[start].append((index, seg.identity / 100))
        events[end].append((index, None))
    active = {}
    score = 0.0
    previous = 0
    for position in sorted(events):
        if active:
            score += (position - previous) * max(active.values())
        for index, identity in events[position]:
            if identity is None:
                active.pop(index, None)
            else:
                active[index] = identity
        previous = position
    return score


def assign(bundle, args):
    metrics, by_query, _ = order.build_match_metrics(
        bundle.path, "manifest", bundle.ref_by_name, bundle.query_by_name,
        args.min_segment_idy, args.min_mapq, args.include_secondary_paf,
    )
    grouped_segments = defaultdict(list)
    for seg in bundle.iter_segments(args.min_segment_idy, args.min_mapq, args.include_secondary_paf):
        grouped_segments[(seg.query, seg.ref)].append(seg)
    scores = {key: union_identity_score(value) for key, value in grouped_segments.items()}
    assignments, labels, candidates = {}, {}, []
    for rec in bundle.query_records:
        matches = sorted(by_query.get(rec.name, []), key=lambda m: (-scores[(m.query, m.ref)], m.ref))
        group_best = {}
        for match in matches:
            group_best.setdefault(bundle.sequences[match.ref].group, match)
        ranked = list(group_best.values())
        best = ranked[0] if ranked else None
        second = ranked[1] if len(ranked) > 1 else None
        score = scores[(rec.name, best.ref)] if best else 0
        second_score = scores[(rec.name, second.ref)] if second else 0
        margin = (score - second_score) / score if score else 0
        best_seq = bundle.sequences[best.ref] if best else None
        backbone = bundle.backbones.get(best_seq.group) if best_seq else None
        placement = next((m for m in matches if m.ref == backbone), None)
        partition, reason = "placed", "assigned"
        if best is None:
            partition, reason = "unplaced", "no_alignment"
        elif best.merged_query_bp < args.min_aligned_bp or best.query_cov < args.min_query_cov:
            partition, reason = "unplaced", "insufficient_alignment"
        elif second and margin <= args.ambiguity_margin:
            partition, reason = "ambiguous", "competing_output_groups"
        elif not placement or placement.merged_query_bp < args.min_aligned_bp or placement.query_cov < args.min_query_cov:
            partition, reason = "unplaced", "insufficient_backbone_alignment"
        if partition != "placed":
            placement = None
        assignment = order.Assignment(
            rec.name, rec.length, reason, partition == "placed", placement, second,
            score / (score + second_score) if score else 0, len(matches),
        )
        assignments[rec.name] = assignment
        labels[rec.name] = {
            "event_id": stable_id("assignment", bundle.assembly_digest, rec.name),
            "partition": partition, "best_ref": best.ref if best else ".",
            "best_score_bp": round(score, 6), "alternate_ref": second.ref if second else ".",
            "alternate_score_bp": round(second_score, 6), "decision_margin": round(margin, 8),
            "ambiguity_margin": args.ambiguity_margin,
            "reference_id": best_seq.reference_id if best_seq else "",
            "chromosome": best_seq.chromosome if best_seq else "",
            "output_group": best_seq.output_group if best_seq else "",
            "subgenome": best_seq.subgenome if best_seq else "",
            "homoeolog_group": best_seq.homoeolog_group if best_seq else "",
            "copy_label": best_seq.copy_label if best_seq else "",
            "scaffold_backbone": backbone or ".",
            "scaffold_id": best_seq.output_id if placement else ".",
            "scaffold_eligible": "yes" if placement else "no",
            "scaffold_conflict": ".",
            "assembly_sha256": bundle.assembly_digest,
            "reference_disagreement": ";".join(m.ref for m in matches if placement and
                                                bundle.sequences[m.ref].group == best_seq.group and
                                                m.orientation != placement.orientation) or ".",
        }
        for rank, match in enumerate(matches, 1):
            seq = bundle.sequences[match.ref]
            candidates.append({
                "contig": rec.name, "rank": rank, "reference_key": seq.key,
                **asdict(seq), "score_bp": round(scores[(rec.name, seq.key)], 6),
                "merged_query_bp": match.merged_query_bp, "query_cov": match.query_cov,
                "ref_start": match.ref_start, "ref_end": match.ref_end + 1,
                "orientation": match.orientation, "coordinate_system": "0-based-half-open",
            })
    if args.filter_within_copy:
        if any(a.kept and not labels[a.query]["copy_label"] for a in assignments.values()):
            raise ValueError("--filter-within-copy requires explicit copy_label for every placed contig")
        order.resolve_duplicate_overlaps(assignments, args)
        for name, assignment in assignments.items():
            if not assignment.kept and labels[name]["partition"] == "placed":
                labels[name]["partition"] = "unplaced"
                labels[name]["scaffold_id"] = "."
                labels[name]["scaffold_eligible"] = "no"
                assignment.best = None
    by_backbone = defaultdict(list)
    for assignment in assignments.values():
        if assignment.kept and assignment.best:
            by_backbone[assignment.best.ref].append(assignment)
    for placed in by_backbone.values():
        placed.sort(key=lambda item: item.best.ref_start)
        for i, left in enumerate(placed):
            for right in placed[i + 1:]:
                if right.best.ref_start > left.best.ref_end:
                    break
                overlap = order.intersect_bp(left.best.merged_ref_intervals, right.best.merged_ref_intervals)
                shorter = min(left.best.merged_ref_bp, right.best.merged_ref_bp)
                if shorter and overlap / shorter >= 0.5:
                    for item in (left, right):
                        labels[item.query]["scaffold_eligible"] = "no"
                        labels[item.query]["scaffold_conflict"] = "overlapping_copy_or_repeat_requires_review"
    return assignments, labels, candidates, metrics


def run(args, bundle=None):
    bundle = bundle or InputBundle(args.manifest, args.assembly_fasta)
    if args.graph_guard and not args.gfa:
        raise ValueError("--graph-guard requires --gfa")
    if not 0 <= args.ambiguity_margin <= 1:
        raise ValueError("--ambiguity-margin must be between 0 and 1")
    assignments, labels, candidates, metrics = assign(bundle, args)
    ordered = order.order_assignments(assignments, [r.name for r in bundle.ref_records],
                                      args.name_separator, args.orient_to_reference)
    for assignment in ordered:
        assignment.new_name = namespaced("placed", assignment.best.ref, assignment.query)
    if args.retain_all:
        for rec in bundle.query_records:
            assignment = assignments[rec.name]
            if not assignment.kept:
                assignment.kept = True
                assignment.new_name = namespaced(labels[rec.name]["partition"], rec.name)
                ordered.append(assignment)
    prefix = str(args.output_prefix)
    paths = {"assignments": prefix + ".contig_assignments.tsv", "matches": prefix + ".contig_ref_matches.tsv",
             "candidates": prefix + ".reference_candidates.tsv", "fasta": prefix + ".ordered.fa",
             "agp": prefix + ".ordered.agp", "components": prefix + ".ordered_components.tsv",
             "accounting": prefix + ".accounting.tsv", "summary": prefix + ".run_summary.txt",
             "chromosome_summary": prefix + ".chromosome_summary.tsv",
             "submission_checklist": prefix + ".submission_checklist.tsv"}
    order.write_assignment_report(paths["assignments"], bundle.query_records, assignments)
    with open(paths["assignments"], newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        row.update(labels[row["contig"]])
    write_tsv(paths["assignments"], rows)
    write_tsv(paths["candidates"], candidates, list(candidates[0]) if candidates else ["contig", "rank", "reference_key"])
    order.write_match_report(paths["matches"], metrics)
    order.write_chromosome_summary(paths["chromosome_summary"], bundle.ref_records, ordered)
    if args.gfa:
        graph = order.read_gfa(args.gfa)
        paths["graph_assignments"] = prefix + ".graph_assignments.tsv"
        order.write_graph_assignment_report(paths["graph_assignments"], bundle.query_records, assignments, graph)
        if args.graph_guard:
            order.graph_guard_sort_warnings(bundle.query_records, assignments, graph)
    parts = order.build_ordered_agp_parts(ordered)
    parts = [replace(part, notes=part.notes + ";" +
              ";".join(f"{k}={labels[part.component_id][k]}" for k in
                       ("partition", "reference_id", "chromosome", "output_group", "subgenome", "copy_label", "scaffold_backbone")))
             for part in parts]
    write_agp(paths["agp"], parts)
    write_component_tsv(paths["components"], parts)
    write_tsv(paths["accounting"], [{
        "input_id": rec.name, "input_start": 0, "input_end": rec.length,
        "output_id": assignments[rec.name].new_name, "output_start": 0,
        "output_end": rec.length if assignments[rec.name].kept else 0,
        "strand": "-" if assignments[rec.name].reverse_complemented else "+",
        "status": "retained" if assignments[rec.name].kept else "excluded",
        "partition": labels[rec.name]["partition"], "reason": assignments[rec.name].status,
    } for rec in bundle.query_records])
    if not args.reports_only:
        order.write_ordered_fasta(paths["fasta"], bundle.assembly, ordered, None, args.simple_headers)
        digest = sha256_file(paths["fasta"])
        for row in rows:
            row["ordered_fasta_sha256"] = digest
        write_tsv(paths["assignments"], rows)
    else:
        paths.pop("fasta")
    order.write_submission_checklist(
        paths["submission_checklist"], "chromo sort --manifest", paths,
        None if args.reports_only else order.build_ordered_fasta_records(bundle.assembly, ordered, None), parts)
    if args.discarded_fasta:
        order.write_discarded_fasta(args.discarded_fasta, bundle.assembly, assignments, None)
        paths["discarded"] = str(args.discarded_fasta)
    Path(paths["summary"]).write_text(
        f"Manifest: {bundle.path}\nAssembly SHA-256: {bundle.assembly_digest}\n"
        f"Contigs retained: {len(ordered)}/{len(bundle.query_records)}\n"
        f"Retain all: {args.retain_all}\nFilter within declared copy: {args.filter_within_copy}\n"
        "Score: maximum identity integrated over the union of query bases per reference frame.\n"
        "Alternative reference group support: maximum per frame, never sum.\n"
        "Reference discordance is a review signal; supplied labels are not inferred phasing.\n")
    audit = audit_record("sort --manifest", bundle.input_paths() + ([args.gfa] if args.gfa else []), paths.values(),
                         {k: v for k, v in vars(args).items()})
    write_json(prefix + ".audit.json", audit)
    return paths
