"""Explicit native-coordinate edits; no inference of biological correctness."""

import csv
import html
from collections import defaultdict

from . import readplot
from .manifest import namespaced
from .multireference import write_tsv
from .provenance import stable_id


EDIT_COLUMNS = ["schema", "edit_id", "event_id", "target", "action", "start", "end",
                "decision", "reviewer", "notes", "assembly_sha256", "scan_sha256"]


def edit_id(row):
    return stable_id("edit", row["scan_sha256"], row["event_id"], row["target"],
                     row["action"], int(row["start"]), int(row["end"]))


def read_edits(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not set(EDIT_COLUMNS).issubset(reader.fieldnames or []):
            raise ValueError("Reviewed edits table is missing required columns")
        rows = list(reader)
    if len({r["edit_id"] for r in rows}) != len(rows):
        raise ValueError("Duplicate edit_id")
    return rows


def validate_edits(rows, proposals, bundle, scan_digest, require_decisions=True):
    for row in rows:
        if row["schema"] != "chromosort-edit-v1":
            raise ValueError("Unsupported interval edit schema")
        if row["assembly_sha256"] != bundle.assembly_digest or row["scan_sha256"] != scan_digest:
            raise ValueError("Reviewed edit belongs to a different assembly or scan")
        event = proposals.get(row["event_id"])
        if event is None or row["target"] != event["target"]:
            raise ValueError("Reviewed edit must link to a scan event on the same target")
        start, end = int(row["start"]), int(row["end"])
        length = bundle.query_by_name[row["target"]].length
        if row["action"] == "reverse":
            valid = 0 <= start < end <= length
        elif row["action"] == "cut":
            valid = 0 < start == end < length
        else:
            raise ValueError("Explicit edits support reverse or cut")
        if not valid:
            raise ValueError("Edit requires native 0-based interval bounds; cuts use equal internal boundaries")
        if row["edit_id"] != edit_id(row):
            raise ValueError("Edit geometry changed; generate a new edit record")
        allowed = {"accept", "reject", "defer"} if require_decisions else {"pending", "accept", "reject", "defer"}
        if row["decision"] not in allowed:
            raise ValueError(f"Explicit edit decision required for {row['edit_id']}")
        if require_decisions and not row["reviewer"].strip():
            raise ValueError("Reviewer label required for interval edits")
        if not row["notes"].strip():
            raise ValueError("Interval edits require a rationale in notes")


def compose_edits(bundle, cuts, rows):
    """Partition native sources; reverse disjoint slices within output groups.

    An accepted cut cannot lie inside a reversed interval. This avoids silently
    remapping a cut to another output boundary or reordering separated records.
    Rejected/deferred edits never contribute boundaries or orientation changes.
    """
    reversals = defaultdict(list)
    for row in rows:
        if row["decision"] != "accept":
            continue
        start, end = int(row["start"]), int(row["end"])
        if row["action"] == "cut":
            if start in cuts[row["target"]]:
                raise ValueError("Two accepted decisions resolve to the same cut")
            cuts[row["target"]].add(start)
        else:
            reversals[row["target"]].append((start, end))
    for target, intervals in reversals.items():
        intervals.sort()
        if any(a[1] > b[0] for a, b in zip(intervals, intervals[1:])):
            raise ValueError("Accepted reverse intervals overlap; use disjoint edits")
        if any(start < cut < end for start, end in intervals for cut in cuts[target]):
            raise ValueError("Accepted cut lies inside a reversed interval; resolve conflicting edits explicitly")
    pieces, output_names = [], set()
    joined = any(reversals.values())
    for rec in bundle.query_records:
        boundaries = [0, *sorted(cuts[rec.name]), rec.length]
        for index, (start, end) in enumerate(zip(boundaries, boundaries[1:]), 1):
            name = namespaced("piece", rec.name, index) if cuts[rec.name] else rec.name
            if name in output_names:
                raise ValueError("Reviewed output names collide; rename the input before preparing a new scan")
            output_names.add(name)
            flips = [(lo, hi) for lo, hi in reversals[rec.name] if start <= lo < hi <= end]
            inner = sorted({start, end, *(x for interval in flips for x in interval)})
            for part, (lo, hi) in enumerate(zip(inner, inner[1:]), 1):
                pieces.append({"id": namespaced("slice", rec.name, index, part) if joined else name, "source": rec.name,
                    "name": name, "start": lo + 1, "end": hi,
                    "strand": "-" if (lo, hi) in flips else "+",
                    "scaffold": name if joined else "unplaced", "removed": False})
    return pieces, joined


def write_edit_review(out, rows, bundle, parameters):
    """Render the exact requested edges, including terminal edits without a cut marker."""
    write_tsv(out / "reviewed_edits.tsv", rows, EDIT_COLUMNS)
    figures = []
    for row in rows:
        length = bundle.query_by_name[row["target"]].length
        window = int(parameters.get("read_window_bp", 10000))
        boundaries = [("start", int(row["start"]))]
        if row["action"] == "reverse":
            boundaries.append(("end", int(row["end"])))
        for edge, position in boundaries:
            relative = "figures/" + row["edit_id"].replace(":", "_") + "_" + edge
            options = ["--manifest", str(bundle.path), "--contig", row["target"],
                "--event-id", row["edit_id"], "--start", str(max(0, position - window)),
                "--end", str(min(length, position + window)), "--action", f"Explicit {row['action']} {edge}: {row['decision']}",
                "--min-anchor-bp", str(parameters.get("read_min_anchor_bp", 1000)),
                "--read-window-bp", str(window), "--bin-bp", str(max(1, window // 100)),
                "--samtools", parameters.get("samtools", "samtools"), "-o", str(out / relative)]
            if 0 < position < length:
                options += ["--breakpoint", str(position)]
            if row["action"] == "reverse":
                options += ["--boundary-label", f"reverse {edge} boundary"]
            if parameters.get("read_evidence_id"):
                options += ["--evidence-id", parameters["read_evidence_id"]]
            readplot.run(readplot.parse_args(options), bundle=bundle)
            figures.append({"edit_id": row["edit_id"], "edge": edge, "position": position, "path": relative})
    escape = lambda value: html.escape(str(value), quote=True)
    table = "".join("<tr>" + "".join(f"<td>{escape(row[key])}</td>" for key in
        ["edit_id", "event_id", "target", "action", "start", "end", "decision", "reviewer", "notes"]) + "</tr>" for row in rows)
    panels = "".join(f'<article><h2>{escape(f["edit_id"])}: {escape(f["edge"])} boundary {f["position"]}</h2>'
        f'<a href="{f["path"]}.pdf">PDF</a> · <a href="{f["path"]}.svg">SVG</a>'
        f'<img src="{f["path"]}.svg" alt="Exact requested edit boundary"></article>' for f in figures)
    (out / "edit-review.html").write_text('<!doctype html><html lang="en"><meta charset="utf-8">'
        '<title>Explicit edit review</title><style>body{font:15px system-ui;margin:32px;color:#243d43}table{border-collapse:collapse}td,th{padding:8px;border:1px solid #ccd}img{width:100%}</style>'
        '<h1>Explicit edit review</h1><p>These coordinates were supplied explicitly; reference discordance and absent read bridges do not establish an error. '
        'Edit the decision, reviewer and notes in <a href="reviewed_edits.tsv">reviewed_edits.tsv</a>. '
        'Accept applies the requested edit; reject and defer preserve its native sequence. Pending edits cannot be applied. '
        'Generate a new edit record to change coordinates or action.</p>'
        '<p>Intervals are 0-based and half-open. Both edges of each reverse interval are shown in native coordinates. '
        'Terminal edges have no internal-breakpoint count. Realign every changed FASTA before interpreting its evidence.</p>'
        '<table><tr><th>Edit</th><th>Linked event</th><th>Target</th><th>Action</th><th>Start</th><th>End</th><th>Decision</th><th>Reviewer</th><th>Rationale</th></tr>'
        + table + '</table>' + panels + '</html>')
    return figures
