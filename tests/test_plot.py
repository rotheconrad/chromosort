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
            self.assertIn("Reference position (bp)", text)
            self.assertIn("Query position (bp)", text)
            self.assertIn("forward", text)
            self.assertIn("reverse", text)

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
