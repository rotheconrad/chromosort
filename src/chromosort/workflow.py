"""Supervised analysis handoff, exact decision replay and fresh-stage validation."""

import argparse
import csv
import html
import json
import subprocess
from collections import defaultdict
from dataclasses import replace
from pathlib import Path

from . import __version__, evaluate, fix_contigs, manual, plot, readplot, repairs, review_proposals
from .agp import group_parts_by_object, read_agp
from .manifest import InputBundle, namespaced
from .multireference import write_tsv
from .provenance import check_digest, file_record, sha256_file, stable_id, write_json
from .reference_order import iter_fasta_records, merge_intervals, reverse_complement
from .review import ReviewEvent, write_review_events


DECISION_COLUMNS = ["schema", "event_id", "target", "action", "proposed_position", "position",
                    "decision", "reason", "reviewer", "notes", "assembly_sha256"]


def read_decisions(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not set(DECISION_COLUMNS).issubset(reader.fieldnames or []):
            raise ValueError("Decision table is missing required columns")
        rows = list(reader)
    if len({r["event_id"] for r in rows}) != len(rows):
        raise ValueError("Duplicate decision event_id")
    return rows


def verify_outputs(assembly, fasta, agp, require_complete=False):
    """Independently reconstruct output strings and account for source intervals."""
    source = {name: seq for name, _, seq in iter_fasta_records(assembly)}
    target = {name: seq for name, _, seq in iter_fasta_records(fasta)}
    if len(source) != sum(1 for _ in iter_fasta_records(assembly)) or len(target) != sum(1 for _ in iter_fasta_records(fasta)):
        raise ValueError("Duplicate FASTA names prevent exact AGP validation")
    parts = group_parts_by_object(read_agp(agp))
    if set(parts) != set(target):
        raise ValueError("FASTA and AGP object names differ")
    used = defaultdict(list)
    for name, entries in parts.items():
        sequence = []
        cursor = 1
        for index, part in enumerate(entries, 1):
            if part.object_start != cursor or part.part_number != index:
                raise ValueError(f"Noncontiguous AGP object: {name}")
            if part.part_type == "gap":
                piece = "N" * part.gap_length
            else:
                if part.component_id not in source:
                    raise ValueError(f"AGP component missing from source: {part.component_id}")
                start, end = part.component_start - 1, part.component_end
                if not 0 <= start < end <= len(source[part.component_id]):
                    raise ValueError(f"AGP component outside source: {part.component_id}")
                piece = source[part.component_id][start:end]
                if part.orientation == "-":
                    piece = reverse_complement(piece)
                elif part.orientation != "+":
                    raise ValueError("Exact reconstruction requires known component orientation")
                used[part.component_id].append((start, end))
            if len(piece) != part.object_end - part.object_start + 1:
                raise ValueError(f"AGP component/object lengths disagree: {name}")
            sequence.append(piece)
            cursor = part.object_end + 1
        if "".join(sequence) != target[name]:
            raise ValueError(f"AGP reconstruction differs from output FASTA: {name}")
    accounting = []
    for name, seq in source.items():
        spans = sorted(used[name])
        union = merge_intervals(spans)
        retained = sum(end - start for start, end in union)
        duplicated = sum(end - start for start, end in spans) - retained
        missing = len(seq) - retained
        if require_complete and (duplicated or missing):
            raise ValueError(f"Source coverage is incomplete or duplicated: {name}; missing={missing}, duplicated={duplicated}")
        accounting.append({"input_id": name, "input_bp": len(seq), "retained_bp": retained,
                           "omitted_bp": missing, "duplicated_bp": duplicated})
    return accounting


def write_dashboard(path, rows, figures, manifest_digest, has_review_proposals=False):
    data = json.dumps(rows).replace("<", "\\u003c")
    figure_html = "".join(f'<article id="{html.escape(eid)}"><h3>{html.escape(eid)}</h3>'
                          f'<a href="{html.escape(name)}.pdf">PDF</a> · <a href="{html.escape(name)}.svg">SVG</a>'
                          f'<img src="{html.escape(name)}.svg" alt="Read evidence for {html.escape(eid)}"></article>'
                          for eid, name in figures)
    page = r"""<!doctype html><html lang="en"><meta charset="utf-8"><title>ChromoSort evidence review</title>
<style>body{font:15px system-ui;color:#243d43;background:#f5f8fa;margin:36px auto;max-width:1400px;padding:0 24px}h1{font-size:30px}a{color:#176958}table{border-collapse:collapse;background:white;width:100%}td,th{padding:12px;border-bottom:1px solid #dbe4e7;text-align:left;vertical-align:top}input,select,button{font:inherit;padding:7px;border:1px solid #b8c7ca;border-radius:4px}input[type=number]{width:110px}input.note{width:160px}button{background:#176958;color:white;cursor:pointer}img{width:100%;background:white}article{margin-top:24px}small{color:#587277}code{word-break:break-all}</style>
<h1>ChromoSort evidence review</h1>
<p>Review each proposed cut against the alignments and read/graph evidence. Accept, reject, defer, or refine the native cut position. Rejected cuts retain the original sequence. Discordance alone does not establish an assembly error.</p>
<p><a href="manual.html">Interactive alignment and graph dashboard</a> · <a href="dotplot.svg">Dot plot</a> · <a href="decisions.tsv">Initial decision table</a>__PROPOSAL_LINK__</p>
<p><label>Reviewer <input id="reviewer" placeholder="Your name or review label"></label> <button id="export">Download decisions.tsv</button> <span id="status"></span></p>
<table><thead><tr><th>Event / contig</th><th>Proposal</th><th>Decision</th><th>Cut after base</th><th>Notes</th></tr></thead><tbody id="rows"></tbody></table>
<p><small>Input manifest SHA-256: __DIGEST__. Applying a decision requires the saved input identities. Realign changed FASTA before further alignment-dependent analysis.</small></p>
__FIGURES__
<script>const rows=__ROWS__;const esc=v=>String(v).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
document.getElementById('rows').innerHTML=rows.map((r,i)=>`<tr><td><a href="#${esc(r.event_id)}">${esc(r.event_id)}</a><br>${esc(r.target)}</td><td>${esc(r.reason)}</td><td><select data-i="${i}" data-k="decision"><option>pending</option><option>accept</option><option>reject</option><option>defer</option><option>refine</option></select></td><td><input type="number" min="1" value="${esc(r.position)}" data-i="${i}" data-k="position" ${r.action!=='cut'?'disabled':''}></td><td><input class="note" data-i="${i}" data-k="notes"></td></tr>`).join('');
document.getElementById('rows').addEventListener('change',e=>{if(e.target.dataset.k)rows[Number(e.target.dataset.i)][e.target.dataset.k]=e.target.value});
document.getElementById('export').onclick=()=>{let reviewer=document.getElementById('reviewer').value.trim();if(!reviewer){document.getElementById('status').textContent='Enter a reviewer label.';return;}rows.forEach(r=>r.reviewer=reviewer);const columns=__COLUMNS__;const clean=v=>String(v??'').replace(/[\t\r\n]/g,' ');const tsv=[columns.join('\t'),...rows.map(r=>columns.map(c=>clean(r[c])).join('\t'))].join('\n')+'\n';let a=document.createElement('a');a.href=URL.createObjectURL(new Blob([tsv],{type:'text/tab-separated-values'}));a.download='decisions.tsv';a.click();URL.revokeObjectURL(a.href);document.getElementById('status').textContent=rows.some(r=>r.decision==='pending')?'Downloaded with pending decisions.':'Decisions ready for reviewed apply.';};</script></html>"""
    Path(path).write_text(page.replace("__ROWS__", data).replace("__COLUMNS__", json.dumps(DECISION_COLUMNS))
                         .replace("__FIGURES__", figure_html).replace("__DIGEST__", manifest_digest)
                         .replace("__PROPOSAL_LINK__", ' · <a href="review-proposals.html">Orientation and alternative-plan evidence</a>' if has_review_proposals else ""))


def scan(args):
    review_proposals.validate_arguments(args)
    if min(args.inspect_min_gap_bp, args.min_segment_bp, args.min_piece_bp, args.read_window_bp) < 1:
        raise ValueError("Inspection, segment, piece and read-window thresholds must be positive")
    bundle = InputBundle(args.manifest)
    out = Path(args.output_dir).resolve()
    if (out / "scan.json").exists():
        raise ValueError("Choose a new scan directory to preserve the existing review record")
    out.mkdir(parents=True, exist_ok=True)
    planner = ["--manifest", str(bundle.path), "-o", str(out / "candidates"), "--mode", args.mode,
               "--min-segment-bp", str(args.min_segment_bp), "--min-mapq", str(args.min_mapq),
               "--min-piece-bp", str(args.min_piece_bp), "--read-window-bp", str(args.read_window_bp),
               "--read-min-anchor-bp", str(args.read_min_anchor_bp)]
    planner += ["--contigs", *args.contigs] if args.contigs else ["--all"]
    planner += ["--read-continuity-policy", args.read_continuity_policy,
                "--continuity-min-mapq", str(args.continuity_min_mapq),
                "--continuity-anchor-bp", str(args.continuity_anchor_bp),
                "--continuity-min-spanning", str(args.continuity_min_spanning), "--samtools", args.samtools]
    if args.read_evidence_id:
        planner += ["--read-evidence-id", args.read_evidence_id]
    for kind, flag in [("gfa", "--gfa"), ("gaf", "--gaf"), ("read_paf", "--read-paf")]:
        entries = [e for e in bundle.evidence if e["kind"] == kind]
        if len(entries) == 1:
            planner += [flag, str(bundle.resolve(entries[0]["path"]))]
    planner_args = evaluate.parse_fix_args(planner)
    events = evaluate.build_fix_events(planner_args)
    rows, dashboard_events, inspection_intervals = [], [], {}
    for event in events:
        # The planner marks the last piece with no breakpoint.
        if event.action == "split_piece":
            position = int(event.fields["slice_end"])
            if position >= bundle.query_by_name[event.target].length:
                continue
            action = "cut"
        else:
            position, action = 0, "keep"
        eid = stable_id("decision", bundle.assembly_digest, event.target, action, position)
        rows.append({"schema": "chromosort-decision-v1", "event_id": eid, "target": event.target,
                     "action": action, "proposed_position": position, "position": position,
                     "decision": "pending", "reason": event.reason, "reviewer": "", "notes": "",
                     "assembly_sha256": bundle.assembly_digest})
        dashboard_events.append(replace(event, event_id=eid, accept=False))
    automatic_rows = list(rows)
    aligned = defaultdict(list)
    for segment in bundle.correction_segments(min_mapq=args.min_mapq):
        aligned[segment.query].append((min(segment.query_start, segment.query_end) - 1,
                                       max(segment.query_start, segment.query_end)))
    for contig, intervals in aligned.items():
        if args.contigs and contig not in args.contigs:
            continue
        merged = merge_intervals(intervals)
        for left, right in zip(merged, merged[1:]):
            start, end = left[1], right[0]
            if end - start < args.inspect_min_gap_bp or min(left[1] - left[0], right[1] - right[0]) < args.min_segment_bp:
                continue
            position = (start + end) // 2
            eid = stable_id("decision", bundle.assembly_digest, contig, "inspect", start, end)
            reason = f"No passing backbone alignment covers internal interval [{start}, {end}); inspect both interval edges. Interior read continuity does not establish continuity across the boundaries. No cut is proposed."
            inspection_intervals[eid] = (start, end)
            rows.append({"schema": "chromosort-decision-v1", "event_id": eid, "target": contig,
                         "action": "inspect", "proposed_position": position, "position": position,
                         "decision": "pending", "reason": reason, "reviewer": "", "notes": "",
                         "assembly_sha256": bundle.assembly_digest})
            dashboard_events.append(ReviewEvent(eid, "fix", "inspect", contig, reason=reason,
                                    fields={"source_contig": contig, "slice_start": start + 1, "slice_end": end,
                                            "assembly_sha256": bundle.assembly_digest}))
    proposal_record = None
    detailed_events = {}
    if args.inspect_orientation or args.alternative_plans:
        proposal_record = review_proposals.build(bundle, args, planner_args, automatic_rows)
        detailed_events = {event["event_id"]: event for event in proposal_record["events"]}
        for event in proposal_record["events"]:
            position = (event["start0"] + event["end0"]) // 2
            reason = " ".join(event["reasons"])
            inspection_intervals[event["event_id"]] = event["positions"]
            rows.append({"schema": "chromosort-decision-v1", "event_id": event["event_id"], "target": event["target"],
                "action": "inspect", "proposed_position": position, "position": position,
                "decision": "pending", "reason": reason, "reviewer": "", "notes": "",
                "assembly_sha256": bundle.assembly_digest})
            dashboard_events.append(ReviewEvent(event["event_id"], "fix", "inspect", event["target"], reason=reason,
                fields={"source_contig": event["target"], "slice_start": event["start0"] + 1,
                        "slice_end": max(event["start0"] + 1, event["end0"]),
                        "inspection_kind": event["kind"], "inspection_positions": json.dumps(event["positions"]),
                        "review_proposal_html": "review-proposals.html#" + event["event_id"],
                        "assembly_sha256": bundle.assembly_digest}))
    write_tsv(out / "decisions.tsv", rows, DECISION_COLUMNS)
    write_review_events(out / "candidates.fix_review.tsv", dashboard_events)
    plots = []
    for row in rows:
        if row["action"] not in {"cut", "inspect"}:
            continue
        positions = list(inspection_intervals[row["event_id"]]) if row["action"] == "inspect" else [row["position"]]
        position = positions[0]
        length = bundle.query_by_name[row["target"]].length
        relative = "figures/" + row["event_id"].replace(":", "_")
        options = ["--manifest", str(bundle.path), "--contig", row["target"], "--event-id", row["event_id"],
                   "--breakpoint", str(position), "--start", str(max(0, position - args.read_window_bp)),
                   "--end", str(min(length, position + args.read_window_bp)), "--action", "proposed cut" if row["action"] == "cut" else "inspect alignment gap: start boundary",
                   "--min-anchor-bp", str(args.read_min_anchor_bp), "--read-window-bp", str(args.read_window_bp),
                   "--bin-bp", str(max(1, args.read_window_bp // 100)), "--samtools", args.samtools, "-o", str(out / relative)]
        if args.read_evidence_id:
            options += ["--evidence-id", args.read_evidence_id]
        for index, point in enumerate(positions):
            suffix = "" if index == 0 else "_end" if index == 1 else f"_boundary{index + 1}"
            key = row["event_id"] + ("" if index == 0 else "-end" if index == 1 else f"-boundary{index + 1}")
            detail = detailed_events.get(row["event_id"])
            action_label = ("inspect " + detail["kind"].replace("_", " ") + f": boundary {index + 1}" if detail else
                            "proposed cut" if row["action"] == "cut" else
                            "inspect alignment gap: " + ("start" if index == 0 else "end") + " boundary")
            for flag, value in (("--breakpoint", point), ("--start", max(0, point - args.read_window_bp)),
                                ("--end", min(length, point + args.read_window_bp)),
                                ("--action", action_label), ("-o", out / (relative + suffix))):
                options[options.index(flag) + 1] = str(value)
            readplot.run(readplot.parse_args(options), bundle=bundle)
            plots.append((key, relative + suffix))
            if detail is not None:
                detail["figures"].append({"position": point, "path": relative + suffix})
    figure_links = dict(plots)
    dashboard_events = [replace(event, fields={**event.fields, "read_figure_svg": figure_links[event.event_id] + ".svg",
                        **({"read_figure_end_svg": figure_links[event.event_id + "-end"] + ".svg"} if event.event_id + "-end" in figure_links else {})})
                        if event.event_id in figure_links else event for event in dashboard_events]
    write_review_events(out / "candidates.fix_review.tsv", dashboard_events)
    from .reference_order import main as sort_main
    sort_main(["--manifest", str(bundle.path), "-o", str(out / "placement"), "--retain-all", "--reports-only",
               "--min-aligned-bp", str(args.min_segment_bp), "--min-mapq", str(args.min_mapq)])
    dashboard_args = ["--manifest", str(bundle.path), "--output-html", str(out / "manual.html"),
                      "--review-table", str(out / "candidates.fix_review.tsv"),
                      "--assignments", str(out / "placement.contig_assignments.tsv")]
    graphs = [e for e in bundle.evidence if e["kind"] == "gfa"]
    if len(graphs) == 1:
        dashboard_args += ["--gfa", str(bundle.resolve(graphs[0]["path"]))]
    manual.main(["fix", *dashboard_args])
    plot.main(["--manifest", str(bundle.path), "-o", str(out / "dotplot"), "--formats", "svg", "pdf", "png"])
    if proposal_record is not None:
        review_proposals.write_report(out, proposal_record)
    write_dashboard(out / "index.html", rows, plots, sha256_file(bundle.path), proposal_record is not None)
    write_json(out / "scan.json", {
        "schema": "chromosort-scan-v1", "chromosort_version": __version__,
        "manifest": file_record(bundle.path), "assembly": file_record(bundle.assembly),
        "inputs": bundle.input_records(),
        "parameters": vars(args), "proposals": rows,
        "decision_columns": DECISION_COLUMNS, "state": "awaiting_human_decisions",
        "figures": plots,
        **({"review_proposals": file_record(out / "review_proposals.json")} if proposal_record is not None else {}),
    })
    return out


def load_scan(scan_dir):
    scan_path = Path(scan_dir).resolve() / "scan.json"
    record = json.loads(scan_path.read_text())
    if record.get("schema") != "chromosort-scan-v1":
        raise ValueError("Unsupported scan record")
    check_digest(record["manifest"]["path"], record["manifest"]["sha256"], "review manifest")
    for entry in record.get("inputs", []):
        check_digest(entry["path"], entry["sha256"], "review input")
    if record.get("review_proposals"):
        check_digest(scan_path.parent / "review_proposals.json", record["review_proposals"]["sha256"], "review proposals")
    bundle = InputBundle(record["manifest"]["path"])
    return scan_path, record, bundle


def propose_edit(args):
    scan_path, record, bundle = load_scan(args.scan_dir)
    proposals = {r["event_id"]: r for r in record["proposals"]}
    if args.event_id not in proposals:
        raise ValueError("Edit must link to an existing scan event")
    row = {"schema": "chromosort-edit-v1", "event_id": args.event_id,
           "target": proposals[args.event_id]["target"], "action": args.action,
           "start": args.start, "end": args.end, "decision": "pending",
           "reviewer": args.reviewer, "notes": args.notes,
           "assembly_sha256": bundle.assembly_digest, "scan_sha256": sha256_file(scan_path)}
    row["edit_id"] = repairs.edit_id(row)
    repairs.validate_edits([row], proposals, bundle, sha256_file(scan_path), require_decisions=False)
    out = Path(args.output_dir).resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError("Choose a new edit directory to preserve the existing review record")
    out.mkdir(parents=True, exist_ok=True)
    figures = repairs.write_edit_review(out, [row], bundle, record["parameters"])
    write_json(out / "edit.json", {"schema": "chromosort-edit-review-v1", "chromosort_version": __version__,
        "scan": file_record(scan_path), "inputs": bundle.input_records(), "edit": row,
        "figures": figures, "state": "pending_explicit_review_no_sequence_change"})
    return out / "reviewed_edits.tsv"


def apply(args):
    scan_path, record, bundle = load_scan(args.scan_dir)
    rows = read_decisions(args.decisions)
    proposals = {r["event_id"]: r for r in record["proposals"]}
    if {r["event_id"] for r in rows} != set(proposals):
        raise ValueError("Every proposed event must appear exactly once in the decision table")
    cuts, outcomes = defaultdict(set), []
    for row in rows:
        proposal = proposals[row["event_id"]]
        for key in ("schema", "event_id", "target", "action", "proposed_position", "assembly_sha256"):
            if str(row[key]) != str(proposal[key]):
                raise ValueError(f"Decision changed immutable proposal field {key}: {row['event_id']}")
        if row["decision"] not in {"accept", "reject", "refine", "defer"}:
            raise ValueError(f"Explicit decision required for {row['event_id']}")
        if not row["reviewer"].strip():
            raise ValueError(f"Reviewer label required for {row['event_id']}")
        position = int(row["position"])
        if row["decision"] != "refine" and position != int(proposal["position"]):
            raise ValueError("Changed cut position requires decision=refine")
        outcome = ("retained_" + row["decision"] if row["decision"] in {"reject", "defer"} else
                   "cut_requested" if row["action"] == "cut" else "acknowledged_no_sequence_edit")
        outcomes.append({"event_id": row["event_id"], "source": "scan", "target": row["target"],
            "action": row["action"], "start": position, "end": position, "decision": row["decision"],
            "outcome": outcome, "reviewer": row["reviewer"], "notes": row["notes"]})
        if row["action"] in {"keep", "inspect"}:
            if row["decision"] == "refine":
                raise ValueError("Refine is supported for proposed cuts; use a saved manual recipe to add other edits")
            continue
        if row["decision"] in {"accept", "refine"}:
            if not 0 < position < bundle.query_by_name[row["target"]].length:
                raise ValueError("Reviewed cut must be inside the native contig")
            if position in cuts[row["target"]]:
                raise ValueError("Two accepted decisions resolve to the same cut")
            cuts[row["target"]].add(position)
    edits = repairs.read_edits(args.reviewed_edits) if args.reviewed_edits else []
    repairs.validate_edits(edits, proposals, bundle, sha256_file(scan_path))
    retained_cuts = {(r["target"], int(r["position"])) for r in rows
                     if r["action"] == "cut" and r["decision"] in {"reject", "defer"}}
    if any(r["action"] == "cut" and r["decision"] == "accept" and
           (r["target"], int(r["start"])) in retained_cuts for r in edits):
        raise ValueError("Explicit cut contradicts a rejected/deferred scan cut at the same boundary")
    pieces, joined = repairs.compose_edits(bundle, cuts, edits)
    outcomes.extend({"event_id": row["edit_id"], "source": "explicit_edit", "target": row["target"],
        "action": row["action"], "start": int(row["start"]), "end": int(row["end"]),
        "decision": row["decision"], "outcome": row["action"] + "_requested" if row["decision"] == "accept" else "retained_" + row["decision"],
        "reviewer": row["reviewer"], "notes": row["notes"]} for row in edits)
    out = Path(args.output_dir).resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError("Choose a new apply directory to preserve the existing audit record")
    out.mkdir(parents=True, exist_ok=True)
    recipe = {"schema": manual.SCHEMA, "chromosortVersion": __version__, "sourceAssembly": str(bundle.assembly),
              "sourceAssemblySha256": bundle.assembly_digest, "scaffoldEnabled": joined, "gapBp": 0,
              "pieces": pieces, "decisions": rows, "reviewedEdits": edits, "requireCompleteCoverage": True}
    write_json(out / "recipe.json", recipe)
    write_tsv(out / "decisions.tsv", rows, DECISION_COLUMNS)
    write_tsv(out / "decision_outcomes.tsv", outcomes,
              ["event_id", "source", "target", "action", "start", "end", "decision", "outcome", "reviewer", "notes"])
    figures = repairs.write_edit_review(out, edits, bundle, record["parameters"]) if edits else []
    plan = {"schema": "chromosort-reviewed-plan-v1", "chromosort_version": __version__,
        "scan": file_record(scan_path), "decisions": file_record(args.decisions),
        "reviewed_edits": file_record(args.reviewed_edits) if args.reviewed_edits else None,
        "inputs": bundle.input_records(), "recipe": file_record(out / "recipe.json"),
        "outcomes": outcomes, "figures": figures,
        "state": "explicit_plan_only_no_sequence_output" if args.plan_only else "explicit_plan_ready_for_apply"}
    write_json(out / "plan.json", plan)
    if args.plan_only:
        return out / "recipe.json"
    fasta = out / "reviewed.fa"
    manual.main(["apply", "--recipe", str(out / "recipe.json"), "--assembly-fasta", str(bundle.assembly),
                 "-o", str(fasta), "--report", str(out / "edits.tsv")])
    accounting = verify_outputs(bundle.assembly, fasta, str(fasta) + ".agp", require_complete=True)
    write_tsv(out / "accounting.tsv", accounting)
    for outcome in outcomes:
        if outcome["outcome"].endswith("_requested"):
            outcome["outcome"] = "applied_" + outcome["action"]
    write_tsv(out / "decision_outcomes.tsv", outcomes,
              ["event_id", "source", "target", "action", "start", "end", "decision", "outcome", "reviewer", "notes"])
    write_json(out / "apply.json", {
        "schema": "chromosort-apply-v1", "chromosort_version": __version__,
        "scan": file_record(scan_path), "decisions": file_record(args.decisions),
        "assembly": file_record(bundle.assembly), "recipe": file_record(out / "recipe.json"),
        "plan": file_record(out / "plan.json"), "reviewed_edits": plan["reviewed_edits"],
        "decision_outcomes": outcomes,
        "output": file_record(fasta), "agp": file_record(str(fasta) + ".agp"),
        "accounting": accounting, "state": "requires_fresh_alignments_and_validation",
    })
    return fasta


def align(args):
    if args.threads < 1:
        raise ValueError("--threads must be positive")
    bundle = InputBundle(args.manifest)
    out = Path(args.output_dir).resolve()
    if (out / "manifest.json").exists():
        raise ValueError("Choose a new alignment directory")
    out.mkdir(parents=True, exist_ok=True)
    assembly = Path(args.assembly_fasta).resolve()
    digest = sha256_file(assembly)
    version = subprocess.run([args.minimap2, "--version"], check=True, text=True, capture_output=True).stdout.strip()
    evidence = []
    for ref_id, reference in bundle.reference_files.items():
        path = out / (namespaced(ref_id) + ".paf")
        command = [args.minimap2, "-x", args.preset, "-c", "--secondary=no", "-t", str(args.threads), str(reference), str(assembly)]
        with open(path, "w") as output, open(str(path) + ".log", "w") as error:
            subprocess.run(command, stdout=output, stderr=error, check=True)
        evidence.append({"evidence_id": ref_id + "_alignment", "kind": "assembly_paf", "path": path.name,
                         "sha256": sha256_file(path), "assembly_sha256": digest, "reference_id": ref_id,
                         "reference_sha256": bundle.reference_digests[ref_id], "command": command,
                         "software": {"minimap2": version}, "coordinate_system": "PAF",
                         "evidence_class": "mapped_assembly", "dependency_group": "reference_alignments"})
    if args.reads:
        path = out / "reads.paf"
        command = [args.minimap2, "-x", args.read_preset, "-c", "--secondary=no", "-t", str(args.threads), str(assembly), str(Path(args.reads).resolve())]
        with open(path, "w") as output, open(str(path) + ".log", "w") as error:
            subprocess.run(command, stdout=output, stderr=error, check=True)
        evidence.append({"evidence_id": "reads", "kind": "read_paf", "path": path.name, "sha256": sha256_file(path),
                         "assembly_sha256": digest, "readset_id": args.read_sample, "read_sample": args.read_sample,
                         "command": command, "software": {"minimap2": version}, "reads": file_record(args.reads),
                         "coordinate_system": "PAF", "evidence_class": "mapped_reads", "dependency_group": args.read_sample})
    refs = [{**row, "fasta": str(bundle.reference_files[row["reference_id"]])} for row in bundle.references]
    manifest = {"schema": "chromosort-input-v1", "assembly": {**bundle.data["assembly"], "path": str(assembly),
                "sha256": digest, "stage": args.stage}, "references": refs,
                "reference_sequences": bundle.rows("reference_sequences"), "evidence": evidence}
    write_json(out / "manifest.json", manifest)
    InputBundle(out / "manifest.json")
    return out / "manifest.json"


def validate(args):
    bundle = InputBundle(args.manifest)
    audit = json.loads(Path(args.apply_audit).read_text())
    if audit.get("schema") != "chromosort-apply-v1":
        raise ValueError("Unsupported apply audit")
    if bundle.assembly_digest != audit["output"]["sha256"]:
        raise ValueError("Validation manifest does not describe the reviewed output FASTA")
    check_digest(audit["assembly"]["path"], audit["assembly"]["sha256"], "original assembly")
    check_digest(audit["agp"]["path"], audit["agp"]["sha256"], "output AGP")
    accounting = verify_outputs(audit["assembly"]["path"], bundle.assembly, audit["agp"]["path"], True)
    if not any(e["kind"] in {"assembly_paf", "assembly_coords"} for e in bundle.evidence):
        raise ValueError("Validation requires fresh assembly alignments")
    segments = list(bundle.iter_segments())
    output = {"schema": "chromosort-validation-v1", "manifest": file_record(bundle.path),
              "apply_audit": file_record(args.apply_audit), "assembly": file_record(bundle.assembly),
              "accounting": accounting, "alignment_rows": len(segments),
              "state": "identity_and_component_reconstruction_validated",
              "biological_validation": "unresolved; requires independent evidence and reviewed interpretation"}
    write_json(args.output, output)
    return output


def main(argv=None, prog=None):
    parser = argparse.ArgumentParser(prog=prog, description="Supervised scan, decision replay, fresh alignment and validation.")
    commands = parser.add_subparsers(dest="task", required=True)
    scan_parser = commands.add_parser("scan")
    scan_parser.add_argument("--manifest", required=True)
    scan_parser.add_argument("--output-dir", required=True)
    scan_parser.add_argument("--mode", choices=fix_contigs.MODE_CHOICES, default="conservative")
    scan_parser.add_argument("--contigs", nargs="*")
    scan_parser.add_argument("--min-segment-bp", type=int, default=10000)
    scan_parser.add_argument("--min-piece-bp", type=int, default=10000)
    scan_parser.add_argument("--inspect-min-gap-bp", type=int, default=1000,
                            help="Review-only unaligned-gap floor between long aligned flanks; does not alter correction smoothing or propose cuts.")
    scan_parser.add_argument("--min-mapq", type=int, default=20)
    review_proposals.add_arguments(scan_parser)
    from .continuity import add_arguments
    add_arguments(scan_parser)
    scan_parser.add_argument("--read-window-bp", type=int, default=10000)
    scan_parser.add_argument("--read-min-anchor-bp", type=int, default=1000)
    apply_parser = commands.add_parser("apply")
    apply_parser.add_argument("--scan-dir", required=True)
    apply_parser.add_argument("--decisions", required=True)
    apply_parser.add_argument("--output-dir", required=True)
    apply_parser.add_argument("--reviewed-edits", help="Explicit native reverse/cut edits linked to this scan; pending edits are rejected.")
    apply_parser.add_argument("--plan-only", action="store_true", help="Write the complete reviewed recipe and evidence preview without a FASTA.")
    edit_parser = commands.add_parser("edit", help="Prepare one pending explicit interval edit and native boundary panels; makes no sequence change.")
    edit_parser.add_argument("--scan-dir", required=True)
    edit_parser.add_argument("--event-id", required=True)
    edit_parser.add_argument("--action", choices=["reverse", "cut"], required=True)
    edit_parser.add_argument("--start", type=int, required=True, help="Native 0-based interval start; cuts use start=end.")
    edit_parser.add_argument("--end", type=int, required=True, help="Native exclusive interval end; reverse complements [start,end).")
    edit_parser.add_argument("--reviewer", default="")
    edit_parser.add_argument("--notes", required=True, help="Rationale for the proposed edit; not an inferred confidence score.")
    edit_parser.add_argument("--output-dir", required=True)
    align_parser = commands.add_parser("align")
    align_parser.add_argument("--manifest", required=True)
    align_parser.add_argument("--assembly-fasta", required=True)
    align_parser.add_argument("--output-dir", required=True)
    align_parser.add_argument("--stage", default="reviewed")
    align_parser.add_argument("--minimap2", default="minimap2")
    align_parser.add_argument("--preset", choices=["asm5", "asm10", "asm20"], default="asm5")
    align_parser.add_argument("--threads", type=int, default=1)
    align_parser.add_argument("--reads")
    align_parser.add_argument("--read-sample", default="unspecified")
    align_parser.add_argument("--read-preset", choices=["map-hifi", "map-ont", "map-pb"], default="map-hifi")
    validate_parser = commands.add_parser("validate")
    validate_parser.add_argument("--manifest", required=True)
    validate_parser.add_argument("--apply-audit", required=True)
    validate_parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    return {"scan": scan, "edit": propose_edit, "apply": apply, "align": align, "validate": validate}[args.task](args)
