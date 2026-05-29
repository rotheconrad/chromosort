#!/usr/bin/env python3
"""Draw lightweight dot plots from MUMmer coords or minimap2 PAF."""

import argparse
import csv
import html
import math
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

from .paths import ensure_parent_dir
from .reference_order import (
    alignment_source_from_args,
    iter_alignments,
    read_fasta_lengths,
)


PLUS_COLOR = "#2563eb"
MINUS_COLOR = "#dc2626"
GRID_COLOR = "#d1d5db"
TEXT_COLOR = "#111827"
MUTED_TEXT = "#6b7280"
BACKGROUND = "#ffffff"


@dataclass(frozen=True)
class AxisRecord:
    name: str
    length: int
    axis_start: int = 1


def parse_args(argv: Optional[Sequence[str]] = None, prog: Optional[str] = None):
    ap = argparse.ArgumentParser(
        prog=prog,
        description="Draw dot plots from existing MUMmer coords or minimap2 PAF alignments.",
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
        help="Assembly/query FASTA.",
    )
    ap.add_argument(
        "--assembly-fai",
        default=None,
        help="Assembly FASTA index. Defaults to <assembly-fasta>.fai when present.",
    )
    alignment_group = ap.add_mutually_exclusive_group(required=True)
    alignment_group.add_argument(
        "-c",
        "--coords",
        help="MUMmer show-coords file.",
    )
    alignment_group.add_argument(
        "--paf",
        help="minimap2 PAF file.",
    )
    ap.add_argument(
        "-o",
        "--output-prefix",
        required=True,
        help=(
            "Output prefix. Writes <prefix>.<format> and, with --per-ref, "
            "<prefix>.<ref>.<format>."
        ),
    )
    ap.add_argument(
        "--formats",
        "--format",
        nargs="+",
        choices=["pdf", "svg", "png"],
        default=["pdf"],
        help="Output format(s). PDF is the default.",
    )
    ap.add_argument(
        "--per-ref",
        action="store_true",
        help="Also write one plot per reference sequence that has plotted alignments.",
    )
    ap.add_argument(
        "--assignments",
        default=None,
        help=(
            "Optional chromo sort <prefix>.contig_assignments.tsv. When provided, "
            "the query axis is ordered by the kept contigs in the assignment report "
            "instead of by assembly FASTA order."
        ),
    )
    ap.add_argument(
        "--per-ref-query-order",
        choices=["fasta", "ref"],
        default="fasta",
        help=(
            "Query order for per-reference plots. fasta preserves assembly FASTA order; "
            "ref orders query contigs by their first hit on that reference."
        ),
    )
    ap.add_argument(
        "--width",
        type=int,
        default=1600,
        help="Plot width in pixels.",
    )
    ap.add_argument(
        "--height",
        type=int,
        default=1200,
        help="Plot height in pixels.",
    )
    ap.add_argument(
        "--min-segment-bp",
        type=int,
        default=0,
        help="Minimum query-aligned bp for a row to be drawn.",
    )
    ap.add_argument(
        "--min-segment-idy",
        type=float,
        default=0.0,
        help="Ignore individual alignment rows below this percent identity.",
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
        help="Maximum number of alignment rows to draw after filtering; 0 means no limit.",
    )
    return ap.parse_args(argv)


def safe_name(name):
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in name)


def load_segments(args, ref_by_name, query_by_name):
    alignment_path, alignment_format = alignment_source_from_args(args)
    segments = []
    skipped_unknown_ref = 0
    skipped_unknown_query = 0
    skipped_short = 0

    for seg in iter_alignments(
        alignment_path,
        alignment_format,
        min_identity=args.min_segment_idy,
        min_mapq=args.min_mapq,
        include_secondary_paf=args.include_secondary_paf,
    ):
        if seg.ref not in ref_by_name:
            skipped_unknown_ref += 1
            continue
        if seg.query not in query_by_name:
            skipped_unknown_query += 1
            continue
        if seg.len_query < args.min_segment_bp:
            skipped_short += 1
            continue
        segments.append(seg)
        if args.max_segments and len(segments) >= args.max_segments:
            break

    return segments, {
        "skipped_unknown_ref": skipped_unknown_ref,
        "skipped_unknown_query": skipped_unknown_query,
        "skipped_short": skipped_short,
    }


def offsets_for_records(records):
    offsets = {}
    cursor = 0
    for rec in records:
        offsets[rec.name] = cursor
        cursor += rec.length
    return offsets, cursor


def read_assignment_order(path, query_by_name, ref_records):
    ref_rank = {rec.name: idx for idx, rec in enumerate(ref_records)}
    ordered = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("kept") != "yes":
                continue
            contig = row.get("contig")
            ref = row.get("assigned_ref", ".")
            if not contig or contig not in query_by_name:
                continue
            order_text = row.get("order_in_ref", ".")
            ref_start_text = row.get("ref_start", ".")
            try:
                order_in_ref = int(order_text)
            except (TypeError, ValueError):
                order_in_ref = 10**12
            try:
                ref_start = int(ref_start_text)
            except (TypeError, ValueError):
                ref_start = 10**12
            ordered.append(
                (
                    ref_rank.get(ref, 10**12),
                    order_in_ref,
                    ref_start,
                    contig,
                    query_by_name[contig],
                )
            )
    return [item[-1] for item in sorted(ordered)]


def axis_starts_for_records(records):
    return {rec.name: getattr(rec, "axis_start", 1) for rec in records}


def segment_coords(seg, x_offsets, y_offsets, x_starts=None, y_starts=None):
    x_offset = x_offsets[seg.ref]
    y_offset = y_offsets[seg.query]
    x_start = (x_starts or {}).get(seg.ref, 1)
    y_start = (y_starts or {}).get(seg.query, 1)
    r0 = x_offset + min(seg.ref_start, seg.ref_end) - x_start
    r1 = x_offset + max(seg.ref_start, seg.ref_end) - x_start + 1
    q0 = y_offset + min(seg.query_start, seg.query_end) - y_start
    q1 = y_offset + max(seg.query_start, seg.query_end) - y_start + 1
    if seg.orientation == "-":
        return r0, q1, r1, q0
    return r0, q0, r1, q1


def scale(value, total, screen_start, screen_size):
    if total <= 0:
        return screen_start
    return screen_start + (value / total) * screen_size


def fmt_bp(value):
    if value >= 1_000_000_000:
        return f"{value / 1_000_000_000:.2f} Gb"
    if value >= 1_000_000:
        return f"{value / 1_000_000:.1f} Mb"
    if value >= 1_000:
        return f"{value / 1_000:.1f} kb"
    return f"{value} bp"


def make_text(x, y, value, size=12, anchor="start", rotate=None, fill=TEXT_COLOR):
    return {
        "type": "text",
        "x": x,
        "y": y,
        "value": str(value),
        "size": size,
        "anchor": anchor,
        "rotate": rotate,
        "fill": fill,
    }


def make_line(x1, y1, x2, y2, stroke, width=1.0, opacity=1.0):
    return {
        "type": "line",
        "x1": x1,
        "y1": y1,
        "x2": x2,
        "y2": y2,
        "stroke": stroke,
        "width": width,
        "opacity": opacity,
    }


def make_rect(x, y, width, height, fill, stroke="none", opacity=1.0):
    return {
        "type": "rect",
        "x": x,
        "y": y,
        "width": width,
        "height": height,
        "fill": fill,
        "stroke": stroke,
        "opacity": opacity,
    }


def svg_text(item):
    x = item["x"]
    y = item["y"]
    rotate = item["rotate"]
    transform = f' transform="rotate({rotate} {x:.2f} {y:.2f})"' if rotate else ""
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" font-size="{item["size"]}" '
        f'text-anchor="{item["anchor"]}" fill="{item["fill"]}"{transform}>'
        f"{html.escape(item['value'])}</text>"
    )


def svg_line(item):
    return (
        f'<line x1="{item["x1"]:.2f}" y1="{item["y1"]:.2f}" '
        f'x2="{item["x2"]:.2f}" y2="{item["y2"]:.2f}" '
        f'stroke="{item["stroke"]}" stroke-width="{item["width"]}" '
        f'opacity="{item["opacity"]}" vector-effect="non-scaling-stroke"/>'
    )


def svg_rect(item):
    return (
        f'<rect x="{item["x"]:.2f}" y="{item["y"]:.2f}" '
        f'width="{item["width"]:.2f}" height="{item["height"]:.2f}" '
        f'fill="{item["fill"]}" stroke="{item["stroke"]}" opacity="{item["opacity"]}"/>'
    )


def draw_separators(elements, records, offsets, total, origin, size, is_x_axis, max_labels=60):
    label_every = max(1, len(records) // max_labels) if records else 1
    for idx, rec in enumerate(records):
        pos = scale(offsets[rec.name], total, origin[0 if is_x_axis else 1], size[0 if is_x_axis else 1])
        if is_x_axis:
            elements.append(make_line(pos, origin[1], pos, origin[1] + size[1], GRID_COLOR, 0.8, 0.65))
            if idx % label_every == 0:
                elements.append(make_text(pos + 3, origin[1] - 8, rec.name, size=10, rotate=-35, fill=MUTED_TEXT))
        else:
            elements.append(make_line(origin[0], pos, origin[0] + size[0], pos, GRID_COLOR, 0.8, 0.65))
            if idx % label_every == 0:
                elements.append(make_text(origin[0] - 8, pos + 3, rec.name, size=10, anchor="end", fill=MUTED_TEXT))


def build_plot_items(title, ref_records, query_records, segments, width, height):
    margin_left = 130
    margin_right = 40
    margin_top = 92
    margin_bottom = 92
    plot_width = max(100, width - margin_left - margin_right)
    plot_height = max(100, height - margin_top - margin_bottom)
    origin = (margin_left, margin_top)
    size = (plot_width, plot_height)

    ref_offsets, ref_total = offsets_for_records(ref_records)
    query_offsets, query_total = offsets_for_records(query_records)
    ref_starts = axis_starts_for_records(ref_records)
    query_starts = axis_starts_for_records(query_records)

    elements = [
        make_rect(0, 0, width, height, BACKGROUND),
        make_text(24, 32, title, size=18),
        make_text(
            24,
            54,
            f"{len(segments)} alignments | x: reference ({fmt_bp(ref_total)}) | y: query ({fmt_bp(query_total)})",
            size=12,
            fill=MUTED_TEXT,
        ),
        make_rect(origin[0], origin[1], plot_width, plot_height, "#f9fafb", stroke="#9ca3af"),
    ]

    draw_separators(elements, ref_records, ref_offsets, ref_total, origin, size, True)
    draw_separators(elements, query_records, query_offsets, query_total, origin, size, False)

    for seg in segments:
        r0, q0, r1, q1 = segment_coords(
            seg,
            ref_offsets,
            query_offsets,
            ref_starts,
            query_starts,
        )
        x1 = scale(r0, ref_total, origin[0], plot_width)
        x2 = scale(r1, ref_total, origin[0], plot_width)
        y1 = scale(q0, query_total, origin[1], plot_height)
        y2 = scale(q1, query_total, origin[1], plot_height)
        color = MINUS_COLOR if seg.orientation == "-" else PLUS_COLOR
        elements.append(make_line(x1, y1, x2, y2, color, width=1.1, opacity=0.62))

    elements.extend(
        [
            make_line(origin[0], origin[1] + plot_height, origin[0] + plot_width, origin[1] + plot_height, "#374151", 1.2),
            make_line(origin[0], origin[1], origin[0], origin[1] + plot_height, "#374151", 1.2),
            make_text(origin[0] + plot_width / 2, height - 24, "Reference position", size=13, anchor="middle"),
            make_text(22, origin[1] + plot_height / 2, "Query position", size=13, anchor="middle", rotate=-90),
            make_line(width - 210, 30, width - 170, 30, PLUS_COLOR, width=2.0),
            make_text(width - 160, 34, "forward", size=12),
            make_line(width - 210, 50, width - 170, 50, MINUS_COLOR, width=2.0),
            make_text(width - 160, 54, "reverse", size=12),
        ]
    )
    return elements


def write_svg(path, elements, width, height):
    ensure_parent_dir(path)
    with open(path, "w") as out:
        out.write(
            '<svg xmlns="http://www.w3.org/2000/svg" '
            f'width="{width}" height="{height}" viewBox="0 0 {width} {height}">\n'
        )
        out.write(
            "<style>text{font-family:Arial,Helvetica,sans-serif;}"
            "line{stroke-linecap:round;}</style>\n"
        )
        for item in elements:
            if item["type"] == "line":
                out.write(svg_line(item) + "\n")
            elif item["type"] == "rect":
                out.write(svg_rect(item) + "\n")
            elif item["type"] == "text":
                out.write(svg_text(item) + "\n")
        out.write("\n</svg>\n")


def hex_rgb(color):
    color = color.lstrip("#")
    return tuple(int(color[i : i + 2], 16) for i in (0, 2, 4))


def pdf_rgb(color):
    r, g, b = hex_rgb(color)
    return f"{r / 255:.4f} {g / 255:.4f} {b / 255:.4f}"


def pdf_escape(value):
    return value.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


def text_width(value, size):
    return len(value) * size * 0.55


def pdf_text(item, page_height):
    value = pdf_escape(item["value"])
    size = item["size"]
    x = item["x"]
    y = page_height - item["y"]
    anchor = item["anchor"]
    if anchor == "middle":
        x -= text_width(item["value"], size) / 2
    elif anchor == "end":
        x -= text_width(item["value"], size)
    fill = pdf_rgb(item["fill"])
    rotate = item["rotate"]
    if rotate:
        theta = math.radians(-rotate)
        cos_t = math.cos(theta)
        sin_t = math.sin(theta)
        return (
            "q\n"
            f"{fill} rg\n"
            f"{cos_t:.6f} {sin_t:.6f} {-sin_t:.6f} {cos_t:.6f} {x:.2f} {y:.2f} cm\n"
            f"BT /F1 {size} Tf 0 0 Td ({value}) Tj ET\n"
            "Q\n"
        )
    return f"{fill} rg\nBT /F1 {size} Tf {x:.2f} {y:.2f} Td ({value}) Tj ET\n"


def write_pdf(path, elements, width, height):
    stream = []
    for item in elements:
        if item["type"] == "rect":
            x = item["x"]
            y = height - item["y"] - item["height"]
            w = item["width"]
            h = item["height"]
            if item["fill"] != "none":
                stream.append(f"{pdf_rgb(item['fill'])} rg\n{x:.2f} {y:.2f} {w:.2f} {h:.2f} re f\n")
            if item["stroke"] != "none":
                stream.append(f"{pdf_rgb(item['stroke'])} RG\n0.8 w\n{x:.2f} {y:.2f} {w:.2f} {h:.2f} re S\n")
        elif item["type"] == "line":
            stream.append(
                f"{pdf_rgb(item['stroke'])} RG\n{item['width']:.2f} w\n"
                f"{item['x1']:.2f} {height - item['y1']:.2f} m "
                f"{item['x2']:.2f} {height - item['y2']:.2f} l S\n"
            )
        elif item["type"] == "text":
            stream.append(pdf_text(item, height))

    content = "".join(stream).encode("latin-1", errors="replace")
    objects = [
        b"<< /Type /Catalog /Pages 2 0 R >>",
        b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
        (
            f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {width} {height}] "
            "/Resources << /Font << /F1 5 0 R >> >> /Contents 4 0 R >>"
        ).encode("ascii"),
        b"<< /Length " + str(len(content)).encode("ascii") + b" >>\nstream\n" + content + b"endstream",
        b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>",
    ]
    ensure_parent_dir(path)
    with open(path, "wb") as out:
        out.write(b"%PDF-1.4\n")
        offsets = []
        for idx, obj in enumerate(objects, start=1):
            offsets.append(out.tell())
            out.write(f"{idx} 0 obj\n".encode("ascii"))
            out.write(obj)
            out.write(b"\nendobj\n")
        xref = out.tell()
        out.write(f"xref\n0 {len(objects) + 1}\n".encode("ascii"))
        out.write(b"0000000000 65535 f \n")
        for offset in offsets:
            out.write(f"{offset:010d} 00000 n \n".encode("ascii"))
        out.write(
            f"trailer << /Size {len(objects) + 1} /Root 1 0 R >>\n"
            f"startxref\n{xref}\n%%EOF\n".encode("ascii")
        )


def write_png(path, elements, width, height):
    ensure_parent_dir(path)
    try:
        from PIL import Image, ImageDraw, ImageFont
    except ImportError as exc:
        raise RuntimeError("PNG output requires Pillow to be installed.") from exc

    image = Image.new("RGB", (width, height), BACKGROUND)
    draw = ImageDraw.Draw(image, "RGBA")
    font_cache = {}

    def font(size):
        if size not in font_cache:
            try:
                font_cache[size] = ImageFont.truetype("Arial.ttf", size)
            except Exception:
                font_cache[size] = ImageFont.load_default()
        return font_cache[size]

    for item in elements:
        if item["type"] == "rect":
            fill = hex_rgb(item["fill"]) if item["fill"] != "none" else None
            stroke = hex_rgb(item["stroke"]) if item["stroke"] != "none" else None
            xy = [item["x"], item["y"], item["x"] + item["width"], item["y"] + item["height"]]
            draw.rectangle(xy, fill=fill, outline=stroke)
        elif item["type"] == "line":
            color = (*hex_rgb(item["stroke"]), int(255 * item["opacity"]))
            draw.line([item["x1"], item["y1"], item["x2"], item["y2"]], fill=color, width=max(1, round(item["width"])))
        elif item["type"] == "text":
            fnt = font(item["size"])
            x = item["x"]
            y = item["y"] - item["size"]
            bbox = draw.textbbox((0, 0), item["value"], font=fnt)
            tw = bbox[2] - bbox[0]
            if item["anchor"] == "middle":
                x -= tw / 2
            elif item["anchor"] == "end":
                x -= tw
            color = hex_rgb(item["fill"])
            if item["rotate"]:
                text_img = Image.new("RGBA", (max(1, tw + 4), max(1, item["size"] + 8)), (255, 255, 255, 0))
                text_draw = ImageDraw.Draw(text_img)
                text_draw.text((2, 2), item["value"], font=fnt, fill=(*color, 255))
                rotated = text_img.rotate(-item["rotate"], expand=True)
                image.paste(rotated, (round(x), round(y)), rotated)
            else:
                draw.text((x, y), item["value"], font=fnt, fill=(*color, 255))
    image.save(path)


def draw_plot(path, title, ref_records, query_records, segments, width, height, output_format):
    elements = build_plot_items(title, ref_records, query_records, segments, width, height)
    if output_format == "svg":
        write_svg(path, elements, width, height)
    elif output_format == "pdf":
        write_pdf(path, elements, width, height)
    elif output_format == "png":
        write_png(path, elements, width, height)
    else:
        raise ValueError(f"Unsupported plot format: {output_format}")


def per_ref_query_records(ref_name, segments, query_records, order):
    hits = [seg for seg in segments if seg.ref == ref_name]
    if not hits:
        return []

    query_names_to_plot = {rec.name for rec in query_records}
    hits = [seg for seg in hits if seg.query in query_names_to_plot]
    if not hits:
        return []

    span_bounds = {}
    for seg in hits:
        start = min(seg.query_start, seg.query_end)
        end = max(seg.query_start, seg.query_end)
        if seg.query not in span_bounds:
            span_bounds[seg.query] = [start, end]
        else:
            span_bounds[seg.query][0] = min(span_bounds[seg.query][0], start)
            span_bounds[seg.query][1] = max(span_bounds[seg.query][1], end)

    if order == "fasta":
        ordered_names = [rec.name for rec in query_records if rec.name in span_bounds]
    else:
        first_ref_pos = defaultdict(lambda: float("inf"))
        for seg in hits:
            first_ref_pos[seg.query] = min(
                first_ref_pos[seg.query],
                min(seg.ref_start, seg.ref_end),
            )
        ordered_names = sorted(
            span_bounds,
            key=lambda name: (first_ref_pos[name], name),
        )

    return [
        AxisRecord(
            name=name,
            length=span_bounds[name][1] - span_bounds[name][0] + 1,
            axis_start=span_bounds[name][0],
        )
        for name in ordered_names
    ]


def main(argv: Optional[Sequence[str]] = None, prog: Optional[str] = None):
    args = parse_args(argv, prog=prog)
    prefix = Path(args.output_prefix)

    ref_records, ref_by_name = read_fasta_lengths(args.ref_fasta, args.ref_fai)
    query_records, query_by_name = read_fasta_lengths(args.assembly_fasta, args.assembly_fai)
    if args.assignments:
        assignment_records = read_assignment_order(args.assignments, query_by_name, ref_records)
        if assignment_records:
            query_records = assignment_records
    segments, skipped = load_segments(args, ref_by_name, query_by_name)
    query_names_to_plot = {rec.name for rec in query_records}
    segments = [seg for seg in segments if seg.query in query_names_to_plot]

    written = []
    for output_format in args.formats:
        whole_plot = Path(f"{prefix}.{output_format}")
        draw_plot(
            whole_plot,
            "ChromoSort whole-genome dot plot",
            ref_records,
            query_records,
            segments,
            args.width,
            args.height,
            output_format,
        )
        written.append(whole_plot)

    if args.per_ref:
        by_ref = defaultdict(list)
        for seg in segments:
            by_ref[seg.ref].append(seg)
        for ref_rec in ref_records:
            ref_segments = by_ref.get(ref_rec.name, [])
            if not ref_segments:
                continue
            query_subset = per_ref_query_records(
                ref_rec.name,
                segments,
                query_records,
                args.per_ref_query_order,
            )
            if not query_subset:
                continue
            query_subset_names = {rec.name for rec in query_subset}
            ref_segments = [seg for seg in ref_segments if seg.query in query_subset_names]
            if not ref_segments:
                continue
            for output_format in args.formats:
                ref_plot = Path(f"{prefix}.{safe_name(ref_rec.name)}.{output_format}")
                draw_plot(
                    ref_plot,
                    f"ChromoSort dot plot: {ref_rec.name}",
                    [ref_rec],
                    query_subset,
                    ref_segments,
                    args.width,
                    args.height,
                    output_format,
                )
                written.append(ref_plot)

    sys.stderr.write(f"Plotted {len(segments)} alignment rows.\n")
    for label, count in skipped.items():
        if count:
            sys.stderr.write(f"  {label}: {count}\n")
    for path in written:
        sys.stderr.write(f"Wrote plot: {path}\n")


if __name__ == "__main__":
    main()
