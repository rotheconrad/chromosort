"""Read continuity can defer automatic cuts; it cannot establish biological truth."""

from .longreads import LongReadEvidence, summarize_breakpoint
from .readalign import load_reads


REPORT_COLUMNS = ["contig", "position", "evidence_id", "evidence_kind", "dependency_group",
                  "read_sample", "min_mapq", "anchor_bp", "min_spanning", "spanning_molecules",
                  "spanning_read_ids", "alignments_without_cigar", "continuity_status",
                  "automatic_action", "policy", "assembly_sha256", "evidence_sha256"]


def add_arguments(parser):
    parser.add_argument("--read-continuity-policy", choices=["review", "report"], default="review",
                        help="Manifest inputs: review defers a whole automatic contig plan when reads span a proposed cut; report leaves the algorithm unchanged. Explicit reviewed plans bypass this advisory guard.")
    parser.add_argument("--read-evidence-id", help="Manifest read source for continuity; correlated PAF/SAM/BAM representations are never pooled.")
    parser.add_argument("--continuity-min-mapq", type=int, default=20)
    parser.add_argument("--continuity-anchor-bp", type=int, default=1000)
    parser.add_argument("--continuity-min-spanning", type=int, default=2,
                        help="Distinct CIGAR-contiguous read molecules that trigger review, not a biological confidence threshold.")
    parser.add_argument("--samtools", default="samtools")


def select_source(bundle, evidence_id=None):
    entries = [e for e in bundle.evidence if e["kind"] in {"read_paf", "sam", "bam"}]
    if evidence_id:
        selected = [e for e in entries if e["evidence_id"] == evidence_id]
        if not selected:
            raise ValueError("Read evidence_id absent from manifest")
        return selected[0]
    if not entries:
        return None
    groups = {e.get("dependency_group") or e.get("readset_id") or e["evidence_id"] for e in entries}
    if len(groups) > 1:
        raise ValueError("Select one read source with --read-evidence-id; independent sources are not pooled")
    # Prefer portable CIGAR-bearing PAF. SAM/BAM flags are not needed to establish
    # that one aligned block spans both anchors. Missing CIGAR stays unavailable.
    return min(entries, key=lambda e: ({"read_paf": 0, "sam": 1, "bam": 2}[e["kind"]], e["evidence_id"]))


def assess_plans(plans, bundle, args, apply=False):
    """Return original boundary evidence; optionally defer the entire contig plan.

    Deferral preserves sequence and orientation. Rejoining a subset of oriented
    pieces would invent a new plan, so selective intervention remains reviewed.
    Lack of spanning reads never authorizes a cut or confirms an assembly error.
    """
    if bundle is None:
        return []
    mapq = getattr(args, "continuity_min_mapq", 20)
    anchor = getattr(args, "continuity_anchor_bp", 1000)
    threshold = getattr(args, "continuity_min_spanning", 2)
    policy = getattr(args, "read_continuity_policy", "review")
    if not 0 <= mapq <= 254 or min(anchor, threshold) < 1:
        raise ValueError("Continuity requires MAPQ 0–254 and positive anchor/molecule thresholds")
    if not any(plan.status == "split" for plan in plans.values()):
        return []
    source = select_source(bundle, getattr(args, "read_evidence_id", None))
    loaded = None
    if source and source["kind"] != "bam":
        key = ("continuity", source["evidence_id"], mapq)
        if key not in bundle.read_cache:
            bundle.read_cache[key] = load_reads(bundle.resolve(source["path"]),
                "paf" if source["kind"] == "read_paf" else "sam", min_mapq=mapq,
                targets=bundle.query_by_name)
        loaded = bundle.read_cache[key]
    rows = []
    for contig, plan in plans.items():
        if plan.status != "split":
            continue
        contig_rows = []
        for piece in plan.pieces[:-1]:
            position = piece.slice_end
            if source and source["kind"] == "bam":
                low, high = max(1, position - anchor), min(bundle.query_by_name[contig].length, position + anchor)
                loaded = load_reads(bundle.resolve(source["path"]), "bam", min_mapq=mapq,
                    samtools=getattr(args, "samtools", "samtools"), region=f"{contig}:{low}-{high}",
                    targets=bundle.query_by_name)
            alignments = loaded.evidence.by_contig.get(contig, []) if loaded else []
            usable = [a for a in alignments if a.aligned_blocks is not None]
            evidence = LongReadEvidence.from_alignments(usable)
            support = summarize_breakpoint(evidence, contig, position, anchor, anchor)
            status = ("unavailable_no_read_evidence" if source is None else
                      "review_read_continuity" if support.spanning_count >= threshold else
                      "unavailable_cigar" if alignments and not usable else "no_spanning_conflict")
            contig_rows.append({"contig": contig, "position": position,
                "evidence_id": source["evidence_id"] if source else ".",
                "evidence_kind": source["kind"] if source else ".",
                "dependency_group": (source.get("dependency_group") or source.get("readset_id") or ".") if source else ".",
                "read_sample": source.get("read_sample", ".") if source else ".",
                "min_mapq": mapq, "anchor_bp": anchor, "min_spanning": threshold,
                "spanning_molecules": support.spanning_count, "spanning_read_ids": ";".join(support.spanning_reads) or ".",
                "alignments_without_cigar": sum(a.aligned_blocks is None for a in alignments),
                "continuity_status": status, "automatic_action": "unchanged", "policy": policy,
                "assembly_sha256": bundle.assembly_digest, "evidence_sha256": source["sha256"] if source else "."})
        conflict = any(r["continuity_status"] == "review_read_continuity" for r in contig_rows)
        if conflict and policy == "review":
            for row in contig_rows:
                row["automatic_action"] = "defer_contig_plan_for_review"
            if apply:
                plan.status = "not_split_read_continuity"
                plan.pieces = []
                plan.reason += " Automatic plan deferred: CIGAR-contiguous read molecules span a proposed cut. Inspect the original proposals with workflow scan or eval fix; continuity is advisory, not biological truth."
        rows.extend(contig_rows)
    return rows
