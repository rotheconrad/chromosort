"""Reproducible native-coordinate read evidence panels.

The region transforms, display-coordinate packing, track order, and visual
design adapt the improved user-authored read-snapshot renderer documented in
docs/read-figure-port.md. No project data or sample constants are included.
"""

import argparse
import hashlib
import json
import math
from collections import defaultdict
from dataclasses import asdict, dataclass, fields
from pathlib import Path

from . import __version__

from .longreads import LongReadEvidence, summarize_breakpoint
from .manifest import InputBundle, table
from .multireference import write_tsv
from .plot import make_line, make_rect, make_text, nice_tick_step, write_pdf, write_png, write_svg
from .provenance import file_record, write_json, check_digest, sha256_file
from .readalign import ReadInput, load_reads
from .reference_order import merge_intervals, open_text


@dataclass(frozen=True)
class Region:
    contig: str
    start: int
    end: int
    reverse: bool = False

    @property
    def length(self):
        return self.end - self.start

    def axis(self, coordinate):
        return self.end - coordinate if self.reverse else coordinate - self.start

    def visible(self, start, end):
        return max(self.start, start) < min(self.end, end)

    def span(self, start, end):
        return sorted((self.axis(max(start, self.start)), self.axis(min(end, self.end))))


def coverage_bins(alignments, region, bin_bp):
    """Mean unique-molecule coverage per bin, using full filtered evidence."""
    molecules = defaultdict(list)
    for aln in alignments:
        blocks = aln.aligned_blocks if aln.aligned_blocks is not None else [(aln.contig_start, aln.contig_end)]
        molecules[aln.read].extend((max(start - 1, region.start), min(end, region.end))
                                 for start, end in blocks if region.visible(start - 1, end))
    sums = [0] * math.ceil(region.length / bin_bp)
    for intervals in molecules.values():
        for start, end in merge_intervals(intervals):
            for index in range((start - region.start) // bin_bp, (end - region.start - 1) // bin_bp + 1):
                low = region.start + index * bin_bp
                high = min(region.end, low + bin_bp)
                sums[index] += max(0, min(end, high) - max(start, low))
    return [value / min(bin_bp, region.length - i * bin_bp) for i, value in enumerate(sums)]


def pack_reads(alignments, region, max_reads, max_rows, priority=()):
    molecules = defaultdict(list)
    for aln in alignments:
        molecules[aln.read].append(aln)
    preferred = set(priority)
    selected = sorted(molecules, key=lambda name: (name not in preferred, hashlib.sha256(name.encode()).hexdigest()))[:max_reads]
    geometries = []
    for name in selected:
        segments = sorted(molecules[name], key=lambda a: region.span(a.contig_start - 1, a.contig_end))
        spans = [region.span(a.contig_start - 1, a.contig_end) for a in segments]
        geometries.append((min(s[0] for s in spans), max(s[1] for s in spans), name, segments))
    row_ends = [-float("inf")] * max_rows
    packed = []
    for start, end, name, segments in sorted(geometries):
        for row, previous in enumerate(row_ends):
            if start > previous + region.length / 150:
                row_ends[row] = end
                packed.append((row, name, segments))
                break
    return packed


def load_overlays(bundle, kind, region):
    rows = []
    for entry in bundle.evidence:
        if entry["kind"] != kind:
            continue
        path = bundle.resolve(entry["path"])
        if kind == "variants" and (str(path).endswith(".vcf") or str(path).endswith(".vcf.gz")):
            with open_text(path) as handle:
                for line in handle:
                    if line.startswith("#"):
                        continue
                    fields = line.rstrip().split("\t")
                    if len(fields) < 8:
                        raise ValueError("Malformed VCF overlay")
                    if fields[0] != region.contig:
                        continue
                    start = int(fields[1]) - 1
                    info = dict(item.split("=", 1) for item in fields[7].split(";") if "=" in item)
                    for alt in fields[4].split(","):
                        if alt == ".":
                            continue
                        label = "SV" if alt.startswith("<") or "[" in alt or "]" in alt else (
                            "SNP" if len(fields[3]) == len(alt) == 1 else "INDEL")
                        rows.append({"contig": region.contig, "start": start,
                                     "end": int(info.get("END", start + len(fields[3]))),
                                     "label": label, "kind": label, "filter": fields[6]})
        else:
            rows.extend(table(path))
    result = []
    for row in rows:
        start, end = int(row["start"]), int(row["end"])
        if start < 0 or end <= start:
            raise ValueError("Overlay intervals must be 0-based half-open with positive length")
        if row["contig"] == region.contig and region.visible(start, end):
            result.append({**row, "start": start, "end": end})
    return result


def render_items(region, alignments, packed, depth, args, summary, features=(), variants=()):
    width, left, right = args.width, 150, args.width - 36
    plot_width = right - left
    feature_height = max(30, len(features) * 20)
    coverage_y = 185 + feature_height
    orientation_y = coverage_y + 100
    variant_y = orientation_y + 82
    reads_y = variant_y + (94 if variants else 30)
    rows = max((row for row, _, _ in packed), default=0) + 1
    axis_y = reads_y + max(90, rows * 15 + 34)
    height = axis_y + 116
    items = [make_rect(0, 0, width, height, "#ffffff")]

    def text(x, y, value, size=12, anchor="start", fill="#263a40"):
        items.append(make_text(x, y, str(value), size=size, anchor=anchor, fill=fill))

    def x(coordinate):
        return left + region.axis(coordinate) / region.length * plot_width

    text(28, 30, "Read evidence | " + args.event_id, 19)
    text(28, 54, f"{summary['read_sample']} reads vs {summary['assembly_sample']} assembly | stage: {summary['stage']}")
    text(28, 74, f"{region.contig}:{region.start:,}-{region.end:,} | {'reversed' if region.reverse else 'native'} display | action: {args.action}")
    text(28, 94, f"MAPQ >= {args.min_mapq}; secondary {'included' if args.include_secondary else 'excluded'}; supplementary retained when available")
    text(28, 114, f"{summary['molecules_in_window']} molecules / {summary['alignments_in_window']} alignments in window; displayed {summary['molecules_displayed']} / {summary['alignments_displayed']}")
    text(28, 134, f"Spanning: {summary['spanning_reads']} | split across boundary: {summary['split_reads']} | clipped near boundary: {summary['clipped_near_boundary']}")
    text(28, 154, f"Evidence role: {summary['evidence_role']}; clipping: {summary['clipping_detail']}; supplementary: {summary['supplementary_detail']}", 11)
    text(22, 186, "Features", 12)
    for index, feature in enumerate(features):
        lo, hi = sorted((x(max(region.start, feature["start"])), x(min(region.end, feature["end"]))))
        y = 181 + index * 20
        items.append(make_rect(lo, y, max(1, hi - lo), 6, "#8392ab"))
        text(min(right - 5, max(left + 5, (lo + hi) / 2)), y - 3, feature["label"][:75], 10,
             anchor="middle" if lo > left + 150 and hi < right - 150 else ("start" if lo <= left + 150 else "end"))
    if not features:
        text(left, 186, "No feature track supplied", 10, fill="#69777b")
    text(22, coverage_y + 22, "Coverage")
    text(22, coverage_y + 39, "molecules", 10)
    items.append(make_rect(left, coverage_y, plot_width, 75, "#f5f8fa", "#dce2e5"))
    maximum = max(1, max(depth or [0]))
    for i, value in enumerate(depth):
        start = region.start + i * args.bin_bp
        end = min(region.end, start + args.bin_bp)
        lo, hi = sorted((x(start), x(end)))
        items.append(make_rect(lo, coverage_y + 75 - value / maximum * 60, max(0.5, hi - lo), value / maximum * 60, "#93aaba"))
    text(left + 8, coverage_y + 13, f"{summary['coverage_detail']}; maximum bin mean {max(depth or [0]):.2f}x; bin {args.bin_bp:,} bp", 10)
    text(22, orientation_y + 22, "Orientation")
    items.append(make_rect(left, orientation_y, plot_width, 61, "#f5f8fa", "#dce2e5"))
    for orientation, color, y in [("+", "#29846e", orientation_y + 21), ("-", "#a65b89", orientation_y + 41)]:
        counts = coverage_bins([a for a in alignments if a.orientation == orientation], region, args.bin_bp)
        for i, value in enumerate(counts):
            start = region.start + i * args.bin_bp
            lo, hi = sorted((x(start), x(min(region.end, start + args.bin_bp))))
            items.append(make_rect(lo, y, max(.5, hi - lo), value / maximum * 17, color))
    text(left + 8, orientation_y + 12, "Native + (green) / - (purple); labels retain native strand on reversed axes", 10)
    text(22, variant_y + 14, "Variants")
    if not variants:
        text(left, variant_y + 14, "Unavailable: no overlapping supplied calls", 10, fill="#69777b")
    else:
        for index, lane in enumerate(("SNP", "INDEL", "SV")):
            y = variant_y + 14 + index * 24
            text(left, y, lane, 10)
            for variant in variants:
                if variant.get("kind", "SV").upper() == lane:
                    px = x(variant["start"])
                    items.append(make_rect(max(left, min(right - 3, px)), y - 8, 3, 10, "#bd743b"))
    text(22, reads_y + 16, "Read molecules")
    items.append(make_rect(left, reads_y, plot_width, axis_y - reads_y - 12, "#fafbfc", "#dce2e5"))
    for row, name, segments in packed:
        y = reads_y + 17 + row * 15
        if len(segments) > 1:
            spans = [region.span(a.contig_start - 1, a.contig_end) for a in segments]
            items.append(make_line(left + min(s[0] for s in spans) / region.length * plot_width, y,
                                   left + max(s[1] for s in spans) / region.length * plot_width, y,
                                   "#98a2a5", dash=[3, 3]))
        for aln in segments:
            color = "#29846e" if aln.orientation == "+" else "#a65b89"
            blocks = aln.aligned_blocks if aln.aligned_blocks is not None else [(aln.contig_start, aln.contig_end)]
            lo, hi = region.span(aln.contig_start - 1, aln.contig_end)
            items.append(make_line(left + lo / region.length * plot_width, y, left + hi / region.length * plot_width, y, color, width=.6))
            for start, end in blocks:
                if region.visible(start - 1, end):
                    a, b = sorted((x(max(region.start, start - 1)), x(min(region.end, end))))
                    item = make_line(a, y, b, y, color, width=2.7,
                                     dash=[5, 2] if aln.is_supplementary else None)
                    item["title"] = f"{name}; native {aln.orientation}; {aln.contig_start}-{aln.contig_end}; MAPQ {aln.mapq}; CIGAR {aln.cigar or 'unavailable'}"
                    items.append(item)
            arrow_coordinate = aln.contig_end if aln.orientation == "+" else aln.contig_start - 1
            if region.start < arrow_coordinate < region.end:
                px = x(arrow_coordinate)
                direction = -1 if (aln.orientation == "-") != region.reverse else 1
                items.extend([make_line(px, y, px - direction * 4, y - 3, color),
                              make_line(px, y, px - direction * 4, y + 3, color)])
            for position, clip in [(aln.contig_start - 1, aln.clip_left), (aln.contig_end, aln.clip_right)]:
                if clip and region.start <= position <= region.end:
                    items.append(make_line(x(position), y - 4, x(position), y + 4, "#e19532", width=2))
    if not alignments:
        text(left + 15, reads_y + 35, "No passing read evidence in this window", 13)
    if args.breakpoint is not None:
        px = x(args.breakpoint)
        items.append(make_line(px, 168, px, axis_y, "#cb4247", width=1.2, dash=[4, 3]))
        label = "inspection boundary" if args.action.startswith("inspect") else "cut after base"
        text(px, axis_y + 48, f"{label} {args.breakpoint:,}", 11, "middle", "#b73339")
    items.append(make_line(left, axis_y, right, axis_y, "#506369"))
    step = max(1, int(nice_tick_step(region.length)))
    for coordinate in range(math.ceil(region.start / step) * step, region.end + 1, step):
        px = x(coordinate)
        items.append(make_line(px, axis_y, px, axis_y + 5, "#506369"))
        text(px, axis_y + 21, f"{coordinate:,}", 10, "middle")
    text(left, axis_y + 70, "Native coordinates: 0-based boundaries. Dashed bars: supplementary; orange marks: terminal clipping.", 10)
    text(left, axis_y + 87, "Selection: spanning/split first, then read-ID hash; rows packed in display coordinates. Coverage uses all passing molecules.", 10)
    text(left, axis_y + 104, "Read continuity and reference discordance are evidence for review, not assembly-error truth.", 10)
    return items, width, height


def run(args, bundle=None):
    bundle = bundle or InputBundle(args.manifest)
    if args.contig not in bundle.query_by_name:
        raise ValueError(f"Unknown assembly contig: {args.contig}")
    length = bundle.query_by_name[args.contig].length
    region = Region(args.contig, args.start, args.end if args.end is not None else length, args.reverse)
    if not 0 <= region.start < region.end <= length:
        raise ValueError("Read window must be inside the exact manifest assembly")
    if args.breakpoint is not None and not region.start < args.breakpoint < region.end:
        raise ValueError("Breakpoint must be strictly inside the displayed window")
    if min(args.bin_bp, args.max_reads, args.max_rows, args.min_anchor_bp, args.min_clip_bp, args.read_window_bp) < 1 or args.width < 1000:
        raise ValueError("Positive bin/read/row/anchor/clip/window limits and width >= 1000 are required")
    available = [e for e in bundle.evidence if e["kind"] in {"read_paf", "sam", "bam"}]
    selected = [e for e in available if e["evidence_id"] == args.evidence_id] if args.evidence_id else available
    if len(selected) > 1:
        raise ValueError("Select one read evidence source with --evidence-id")
    if args.evidence_id and not selected:
        raise ValueError("Read evidence_id absent from manifest")
    entry = selected[0] if selected else None
    cache_key = (entry["evidence_id"], args.min_mapq, args.include_secondary) if entry else None
    if entry and entry["kind"] != "bam" and cache_key in bundle.read_cache:
        loaded = bundle.read_cache[cache_key]
    else:
        loaded = load_reads(bundle.resolve(entry["path"]), "paf" if entry["kind"] == "read_paf" else entry["kind"],
                            min_mapq=args.min_mapq, include_secondary=args.include_secondary,
                            samtools=args.samtools, region=f"{region.contig}:{region.start + 1}-{region.end}",
                            targets=bundle.query_by_name) if entry else ReadInput(
                                LongReadEvidence.from_alignments([]), {}, "unavailable", "unavailable")
        if entry and entry["kind"] != "bam":
            bundle.read_cache[cache_key] = loaded
    evidence = loaded.evidence
    alignments = [a for a in evidence.by_contig.get(region.contig, []) if region.visible(a.contig_start - 1, a.contig_end)]
    support = summarize_breakpoint(evidence, region.contig, args.breakpoint, args.read_window_bp, args.min_anchor_bp) if args.breakpoint else None
    priority = [*support.spanning_reads, *support.split_reads] if support else []
    packed = pack_reads(alignments, region, args.max_reads, args.max_rows, priority)
    summary = {
        "schema": "chromosort-read-figure-v1", "chromosort_version": __version__, "event_id": args.event_id,
        "manifest": file_record(bundle.path),
        "input_digests": sorted(record["sha256"] for record in bundle.input_records()),
        "assembly_sha256": bundle.assembly_digest, "stage": bundle.data["assembly"].get("stage", "unspecified"),
        "assembly_sample": bundle.data["assembly"].get("sample_id", "unspecified"),
        "read_sample": (entry.get("read_sample") or entry.get("readset_id", "unspecified")) if entry else "unavailable",
        "evidence_role": entry.get("evidence_class", "unspecified") if entry else "unavailable",
        "region": asdict(region), "settings": vars(args), "filters": loaded.counts,
        "coverage_detail": loaded.detail, "spanning_reads": support.spanning_count if support else ".",
        "split_reads": support.split_count if support else ".", "support": asdict(support) if support else None,
        "clipped_near_boundary": len({a.read for a in alignments if args.breakpoint is not None and
            (((a.clip_left or 0) >= args.min_clip_bp and abs(a.contig_start - 1 - args.breakpoint) <= args.read_window_bp)
             or ((a.clip_right or 0) >= args.min_clip_bp and abs(a.contig_end - args.breakpoint) <= args.read_window_bp))}) if loaded.format in {"sam", "bam"} else "unavailable",
        "molecules_in_window": len({a.read for a in alignments}), "alignments_in_window": len(alignments),
        "molecules_displayed": len(packed), "alignments_displayed": sum(len(a) for _, _, a in packed),
        "supplementary_detail": "available" if loaded.format in {"bam", "sam"} else "unavailable in PAF",
        "clipping_detail": "SAM terminal CIGAR clips" if loaded.format in {"bam", "sam"} else "unavailable in PAF",
        "bam_partner_scope": "overlapping alignments only; SA is annotation, not an extra alignment",
    }
    features, variants = load_overlays(bundle, "features", region), load_overlays(bundle, "variants", region)
    depth = coverage_bins(alignments, region, args.bin_bp)
    items, width, height = render_items(region, alignments, packed, depth, args, summary, features, variants)
    outputs = {}
    for fmt in args.formats:
        path = str(args.output_prefix) + "." + fmt
        {"svg": write_svg, "pdf": write_pdf, "png": write_png}[fmt](path, items, width, height)
        outputs[fmt] = file_record(path, Path(args.output_prefix).parent)
    summary["outputs"] = outputs
    summary["evidence"] = entry
    summary["depth_bins"] = depth
    write_json(str(args.output_prefix) + ".figure.json", summary)
    write_tsv(str(args.output_prefix) + ".reads.tsv", [
        {"row": row, **asdict(aln)} for row, _, segments in packed for aln in segments],
        ["row", *[field.name for field in fields(alignments[0])]]
        if alignments else ["row", "read"])
    return summary


def parse_args(argv=None, prog=None):
    parser = argparse.ArgumentParser(prog=prog, description="Render event-linked read evidence as reproducible SVG/PDF/PNG.")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--evidence-id")
    parser.add_argument("--contig", required=True)
    parser.add_argument("--start", type=int, default=0, help="0-based inclusive boundary")
    parser.add_argument("--end", type=int, help="0-based exclusive boundary; default contig length")
    parser.add_argument("--breakpoint", type=int, help="Cut after this many native bases")
    parser.add_argument("--event-id", default="overview")
    parser.add_argument("--action", default="review")
    parser.add_argument("--reverse", action="store_true")
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--include-secondary", action="store_true")
    parser.add_argument("--min-anchor-bp", type=int, default=1000)
    parser.add_argument("--min-clip-bp", type=int, default=100, help="Minimum terminal CIGAR clip for near-boundary molecule counts")
    parser.add_argument("--read-window-bp", type=int, default=5000)
    parser.add_argument("--bin-bp", type=int, default=100)
    parser.add_argument("--max-reads", type=int, default=100)
    parser.add_argument("--max-rows", type=int, default=40)
    parser.add_argument("--width", type=int, default=1400)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--formats", nargs="+", choices=["svg", "pdf", "png"], default=["svg", "pdf", "png"])
    parser.add_argument("-o", "--output-prefix", required=True)
    return parser.parse_args(argv)


def main(argv=None, prog=None):
    import sys
    from types import SimpleNamespace
    argv = sys.argv[1:] if argv is None else list(argv)
    if argv and argv[0] == "replay":
        parser = argparse.ArgumentParser(prog=f"{prog or 'chromo reads'} replay")
        parser.add_argument("--figure-json", required=True)
        parser.add_argument("--manifest", help="Relocated manifest with the same content identity")
        parser.add_argument("-o", "--output-prefix", required=True)
        options = parser.parse_args(argv[1:])
        saved = json.loads(Path(options.figure_json).read_text())
        if saved.get("schema") != "chromosort-read-figure-v1" or saved["chromosort_version"] != __version__:
            raise ValueError("Figure replay requires the saved schema and ChromoSort version")
        path = options.manifest or saved["manifest"]["path"]
        check_digest(path, saved["manifest"]["sha256"], "figure manifest")
        bundle = InputBundle(path)
        if sorted(sha256_file(p) for p in bundle.input_paths()) != saved["input_digests"]:
            raise ValueError("Figure inputs changed since rendering")
        settings = dict(saved["settings"], manifest=str(path), output_prefix=options.output_prefix)
        return run(SimpleNamespace(**settings))
    return run(parse_args(argv, prog))
