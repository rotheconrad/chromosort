import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().parent / "data"


def run_plot(tmp_path, *alignment_args, formats=None):
    prefix = tmp_path / "sample_plot"
    cmd = [
        sys.executable,
        "-m",
        "chromosort.plot",
        "--ref-fasta",
        str(DATA / "ref.fa"),
        "--assembly-fasta",
        str(DATA / "assembly.fa"),
        *alignment_args,
        "--output-prefix",
        str(prefix),
        "--per-ref",
        "--width",
        "900",
        "--height",
        "700",
    ]
    if formats:
        cmd.extend(["--formats", *formats])
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    subprocess.run(cmd, cwd=ROOT, check=True, env=env)
    return prefix


def run_sort_for_assignments(tmp_path):
    prefix = tmp_path / "sample_sort"
    cmd = [
        sys.executable,
        "-m",
        "chromosort.reference_order",
        "--ref-fasta",
        str(DATA / "ref.fa"),
        "--assembly-fasta",
        str(DATA / "assembly.fa"),
        "--paf",
        str(DATA / "sample.paf"),
        "--output-prefix",
        str(prefix),
        "--min-aligned-bp",
        "10",
        "--min-query-cov",
        "0.10",
        "--split-candidate-min-aligned-bp",
        "20",
        "--split-candidate-min-query-frac",
        "0.20",
    ]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    subprocess.run(cmd, cwd=ROOT, check=True, env=env)
    return prefix.with_suffix(".contig_assignments.tsv")


class PlotTests(unittest.TestCase):
    def test_coords_plot_writes_pdf_by_default(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = run_plot(Path(tmp), "--coords", str(DATA / "sample.coords"))
            pdf = prefix.with_suffix(".pdf")
            self.assertTrue(pdf.exists())
            self.assertEqual(pdf.read_bytes()[:5], b"%PDF-")

    def test_coords_plot_can_write_svg(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = run_plot(
                Path(tmp),
                "--coords",
                str(DATA / "sample.coords"),
                formats=["svg"],
            )
            svg = prefix.with_suffix(".svg")
            self.assertTrue(svg.exists())
            text = svg.read_text()
            self.assertIn("<svg", text)
            self.assertIn("ChromoSort whole-genome dot plot", text)
            self.assertIn("Reference Position (bp)", text)
            self.assertIn("Query Position (bp)", text)
            self.assertIn("stroke-dasharray", text)
            self.assertIn("forward", text)
            self.assertIn("reverse", text)

    def test_query_axis_title_uses_same_gap_as_reference_axis_title(self):
        sys.path.insert(0, str(ROOT / "src"))
        from chromosort.plot import build_plot_items
        from chromosort.reference_order import FastaIndexRecord

        records = [FastaIndexRecord("chr1", 200, 0, 0, 0)]
        elements = build_plot_items("test", records, records, [], 900, 700)
        plot_rect = next(
            item for item in elements if item["type"] == "rect" and item["fill"] == "#f9fafb"
        )
        reference_label = next(
            item for item in elements if item["type"] == "text" and item["value"].startswith("Reference Position")
        )
        query_label = next(
            item for item in elements if item["type"] == "text" and item["value"].startswith("Query Position")
        )
        bottom_gap = reference_label["y"] - (plot_rect["y"] + plot_rect["height"])
        right_gap = query_label["x"] - (plot_rect["x"] + plot_rect["width"])

        self.assertAlmostEqual(right_gap, bottom_gap)

    def test_position_grid_is_visually_distinct_from_sequence_boundaries(self):
        sys.path.insert(0, str(ROOT / "src"))
        from chromosort.plot import (
            BOUNDARY_GRID_COLOR,
            MAJOR_POSITION_GRID_COLOR,
            MINOR_POSITION_GRID_COLOR,
            build_plot_items,
        )
        from chromosort.reference_order import FastaIndexRecord

        records = [
            FastaIndexRecord("chr1", 100, 0, 0, 0),
            FastaIndexRecord("chr2", 100, 0, 0, 0),
        ]
        elements = build_plot_items("test", records, records, [], 900, 700)
        major_grid_lines = [
            item for item in elements
            if item["type"] == "line" and item["stroke"] == MAJOR_POSITION_GRID_COLOR
        ]
        minor_grid_lines = [
            item for item in elements
            if item["type"] == "line" and item["stroke"] == MINOR_POSITION_GRID_COLOR
        ]
        boundary_lines = [
            item for item in elements
            if item["type"] == "line" and item["stroke"] == BOUNDARY_GRID_COLOR
        ]

        self.assertTrue(major_grid_lines)
        self.assertTrue(minor_grid_lines)
        self.assertTrue(boundary_lines)
        self.assertTrue(all(item["dash"] for item in major_grid_lines))
        self.assertTrue(all(item["dash"] for item in minor_grid_lines))
        self.assertTrue(all(item["dash"] is None for item in boundary_lines))
        self.assertGreater(boundary_lines[0]["width"], major_grid_lines[0]["width"])
        self.assertGreater(major_grid_lines[0]["width"], minor_grid_lines[0]["width"])
        self.assertGreater(major_grid_lines[0]["opacity"], minor_grid_lines[0]["opacity"])

    def test_major_and_minor_axis_ticks_are_drawn(self):
        sys.path.insert(0, str(ROOT / "src"))
        from chromosort.plot import (
            AXIS_COLOR,
            MAJOR_TICK_LENGTH,
            MAJOR_TICK_WIDTH,
            MINOR_TICK_LENGTH,
            MINOR_TICK_WIDTH,
            build_plot_items,
        )
        from chromosort.reference_order import FastaIndexRecord

        records = [FastaIndexRecord("chr1", 200, 0, 0, 0)]
        elements = build_plot_items("test", records, records, [], 900, 700)
        plot_rect = next(
            item for item in elements if item["type"] == "rect" and item["fill"] == "#f9fafb"
        )
        bottom = plot_rect["y"] + plot_rect["height"]
        right = plot_rect["x"] + plot_rect["width"]
        tick_lines = [
            item for item in elements
            if item["type"] == "line" and item["stroke"] == AXIS_COLOR
        ]
        x_tick_lengths = [
            item["y2"] - item["y1"]
            for item in tick_lines
            if item["x1"] == item["x2"] and item["y1"] == bottom and item["y2"] > bottom
        ]
        y_tick_lengths = [
            item["x2"] - item["x1"]
            for item in tick_lines
            if item["y1"] == item["y2"] and item["x1"] == right and item["x2"] > right
        ]

        self.assertIn(MAJOR_TICK_LENGTH, x_tick_lengths)
        self.assertIn(MINOR_TICK_LENGTH, x_tick_lengths)
        self.assertIn(MAJOR_TICK_LENGTH, y_tick_lengths)
        self.assertIn(MINOR_TICK_LENGTH, y_tick_lengths)
        self.assertGreater(MAJOR_TICK_WIDTH, MINOR_TICK_WIDTH)

    def test_paf_plot_writes_per_reference_svgs(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = run_plot(
                Path(tmp),
                "--paf",
                str(DATA / "sample.paf"),
                formats=["svg", "png"],
            )
            self.assertTrue(prefix.with_suffix(".svg").exists())
            self.assertTrue(prefix.with_suffix(".png").exists())
            self.assertEqual(prefix.with_suffix(".png").read_bytes()[:8], b"\x89PNG\r\n\x1a\n")
            self.assertTrue(Path(f"{prefix}.chr1.svg").exists())
            self.assertTrue(Path(f"{prefix}.chr2.png").exists())

    def test_sel_ref_limits_whole_and_per_reference_plots(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = run_plot(
                Path(tmp),
                "--paf",
                str(DATA / "sample.paf"),
                "--sel-ref",
                "chr2",
                formats=["svg"],
            )
            whole_svg = prefix.with_suffix(".svg").read_text()
            self.assertIn("2 alignments | x: reference (100 bp) | y: query (160 bp)", whole_svg)
            self.assertNotIn(">chr1<", whole_svg)
            self.assertIn(">chr2<", whole_svg)
            self.assertFalse(Path(f"{prefix}.chr1.svg").exists())
            self.assertTrue(Path(f"{prefix}.chr2.svg").exists())

    def test_sel_ref_reports_missing_reference(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / "sample_plot"
            cmd = [
                sys.executable,
                "-m",
                "chromosort.plot",
                "--ref-fasta",
                str(DATA / "ref.fa"),
                "--assembly-fasta",
                str(DATA / "assembly.fa"),
                "--paf",
                str(DATA / "sample.paf"),
                "--sel-ref",
                "chrMissing",
                "--output-prefix",
                str(prefix),
            ]
            env = os.environ.copy()
            env["PYTHONPATH"] = str(ROOT / "src")
            result = subprocess.run(cmd, cwd=ROOT, env=env, text=True, capture_output=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("--sel-ref reference sequence(s) not found", result.stderr)

    def test_per_ref_plot_crops_query_axis_to_current_reference_hits(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = run_plot(
                Path(tmp),
                "--paf",
                str(DATA / "sample.paf"),
                formats=["svg"],
            )
            svg = Path(f"{prefix}.chr2.svg").read_text()
            self.assertIn("2 alignments | x: reference (100 bp) | y: query (90 bp)", svg)
            self.assertNotIn("y: query (160 bp)", svg)

    def test_per_ref_plot_compacts_unaligned_query_gaps(self):
        sys.path.insert(0, str(ROOT / "src"))
        from chromosort.plot import (
            axis_intervals_for_records,
            offsets_for_records,
            per_ref_query_records,
            segment_coords,
        )
        from chromosort.reference_order import FastaIndexRecord, Segment

        ref_records = [FastaIndexRecord("chr1", 200, 0, 0, 0)]
        query_records = [FastaIndexRecord("contig1", 100, 0, 0, 0)]
        segments = [
            Segment("chr1", "contig1", 1, 10, 1, 10, 10, 10, 99.0, 200, 100, "+"),
            Segment("chr1", "contig1", 11, 20, 91, 100, 10, 10, 99.0, 200, 100, "+"),
        ]

        compact_records = per_ref_query_records("chr1", segments, query_records, "fasta")
        self.assertEqual(compact_records[0].length, 20)
        ref_offsets, _ = offsets_for_records(ref_records)
        query_offsets, _ = offsets_for_records(compact_records)
        _, q0, _, q1 = segment_coords(
            segments[1],
            ref_offsets,
            query_offsets,
            axis_intervals_for_records(ref_records),
            axis_intervals_for_records(compact_records),
        )
        self.assertEqual((q0, q1), (10, 20))

    def test_plot_can_use_assignment_order_without_realignment(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            assignments = run_sort_for_assignments(tmp_path)
            prefix = tmp_path / "ordered_plot"
            cmd = [
                sys.executable,
                "-m",
                "chromosort.plot",
                "--ref-fasta",
                str(DATA / "ref.fa"),
                "--assembly-fasta",
                str(DATA / "assembly.fa"),
                "--paf",
                str(DATA / "sample.paf"),
                "--assignments",
                str(assignments),
                "--output-prefix",
                str(prefix),
            ]
            env = os.environ.copy()
            env["PYTHONPATH"] = str(ROOT / "src")
            subprocess.run(cmd, cwd=ROOT, check=True, env=env)
            self.assertTrue(prefix.with_suffix(".pdf").exists())


if __name__ == "__main__":
    unittest.main()
