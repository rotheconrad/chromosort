import csv
import io
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = Path(__file__).resolve().parent / "data"
GRAPH_DATA = DATA / "graph_gotchas"


def run_chromosort(tmp_path, *extra_args, alignment_args=None):
    prefix = tmp_path / "sample"
    if alignment_args is None:
        alignment_args = ["--coords", str(DATA / "sample.coords")]
    cmd = [
        sys.executable,
        "-m",
        "chromosort.reference_order",
        "--ref-fasta",
        str(DATA / "ref.fa"),
        "--assembly-fasta",
        str(DATA / "assembly.fa"),
        *alignment_args,
        "--output-prefix",
        str(prefix),
        "--min-aligned-bp",
        "10",
        "--min-query-cov",
        "0.10",
        "--min-novel-ref-bp",
        "1",
        "--orient-to-reference",
        *extra_args,
    ]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    subprocess.run(cmd, cwd=ROOT, check=True, env=env)
    return prefix


def read_assignments(path):
    with open(path, newline="") as fh:
        return {
            row["contig"]: row
            for row in csv.DictReader(fh, delimiter="\t")
        }


def fasta_headers(path):
    headers = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                headers.append(line[1:].strip().split()[0])
    return headers


def write_coords(path, rows):
    lines = [
        "/tmp/ref.fa /tmp/assembly.fa",
        "NUCMER",
        "",
        "    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |  [LEN R]  [LEN Q]  |  [COV R]  [COV Q]  | [TAGS]",
        "===============================================================================================================================",
    ]
    for ref, query, rs, re, qs, qe, ref_len, query_len in rows:
        len_ref = abs(re - rs) + 1
        len_query = abs(qe - qs) + 1
        lines.append(
            f"{rs:8d} {re:8d}  | {qs:8d} {qe:8d}  | "
            f"{len_ref:8d} {len_query:8d}  |    99.00  | "
            f"{ref_len:8d} {query_len:8d}  |     1.00     1.00  | {ref}\t{query}"
        )
    path.write_text("\n".join(lines) + "\n")


def run_custom_sort(tmp_path, ref_fasta, assembly_fasta, coords, *extra_args):
    prefix = tmp_path / "custom"
    cmd = [
        sys.executable,
        "-m",
        "chromosort.reference_order",
        "--ref-fasta",
        str(ref_fasta),
        "--assembly-fasta",
        str(assembly_fasta),
        "--coords",
        str(coords),
        "--output-prefix",
        str(prefix),
        "--min-aligned-bp",
        "10",
        "--min-query-cov",
        "0.50",
        "--min-novel-ref-bp",
        "1",
        *extra_args,
    ]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    subprocess.run(cmd, cwd=ROOT, check=True, env=env)
    return prefix


def run_graph_sort(tmp_path, *extra_args):
    prefix = tmp_path / "graph"
    cmd = [
        sys.executable,
        "-m",
        "chromosort.reference_order",
        "--ref-fasta",
        str(GRAPH_DATA / "ref.fa"),
        "--assembly-fasta",
        str(GRAPH_DATA / "assembly.fa"),
        "--paf",
        str(GRAPH_DATA / "unitig_to_ref.paf"),
        "--output-prefix",
        str(prefix),
        "--min-aligned-bp",
        "4",
        "--min-query-cov",
        "0.50",
        "--min-novel-ref-bp",
        "1",
        "--gfa",
        str(GRAPH_DATA / "unitigs.gfa"),
        *extra_args,
    ]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    subprocess.run(cmd, cwd=ROOT, check=True, env=env)
    return prefix


class ReferenceOrderTests(unittest.TestCase):
    def test_graph_guard_warns_when_discarded_overlap_has_direct_graph_link(self):
        from chromosort.graph import read_gfa
        from chromosort.reference_order import (
            Assignment,
            graph_guard_sort_warnings,
        )

        graph = read_gfa(GRAPH_DATA / "unitigs.gfa")
        assignments = {
            "repeat_shared": Assignment(
                query="repeat_shared",
                query_length=8,
                status="duplicate_overlap",
                kept=False,
                best=None,
                second=None,
                best_ref_share=1.0,
                total_refs_matched=1,
                overlap_best_contig="bridge_good",
            ),
            "bridge_good": Assignment(
                query="bridge_good",
                query_length=8,
                status="kept",
                kept=True,
                best=None,
                second=None,
                best_ref_share=1.0,
                total_refs_matched=1,
                new_name="chr1_bridge_good",
            ),
        }
        stream = io.StringIO()

        graph_guard_sort_warnings(
            [SimpleNamespace(name="repeat_shared")],
            assignments,
            graph,
            stream,
        )

        warning = stream.getvalue()
        self.assertIn("graph guard", warning)
        self.assertIn("repeat_shared", warning)
        self.assertIn("repeat_shared+>bridge_good+", warning)

    def test_reference_order_and_duplicate_overlap_status(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = run_chromosort(Path(tmp))

            assignments = read_assignments(prefix.with_suffix(".contig_assignments.tsv"))
            self.assertEqual(assignments["contigA"]["status"], "kept")
            self.assertEqual(assignments["contigB"]["status"], "kept")
            self.assertEqual(assignments["contigC"]["status"], "kept")
            self.assertEqual(assignments["contigNo"]["status"], "no_alignment")

            duplicate = assignments["contigDup"]
            self.assertEqual(duplicate["status"], "duplicate_overlap")
            self.assertEqual(duplicate["kept"], "no")
            self.assertEqual(duplicate["overlap_best_contig"], "contigA")
            self.assertEqual(int(duplicate["novel_ref_bp"]), 0)

            self.assertEqual(
                fasta_headers(prefix.with_suffix(".ordered.fa")),
                [
                    "chr1_contigA",
                    "chr1_contigB",
                    "chr2_contigC",
                ],
            )

    def test_discarded_fasta_creates_missing_parent_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            discarded_fasta = tmp_path / "new_output_dir" / "discarded.fa"

            run_chromosort(
                tmp_path,
                "--discarded-fasta",
                str(discarded_fasta),
            )

            self.assertTrue(discarded_fasta.exists())
            self.assertIn(">contigNo status=no_alignment", discarded_fasta.read_text())

    def test_no_overlap_filter_keeps_otherwise_good_duplicate(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = run_chromosort(Path(tmp), "--no-overlap-filter")

            assignments = read_assignments(prefix.with_suffix(".contig_assignments.tsv"))
            self.assertEqual(assignments["contigDup"]["status"], "kept")

            headers = fasta_headers(prefix.with_suffix(".ordered.fa"))
            self.assertIn("chr1_contigDup", headers)

    def test_split_candidate_is_protected_from_duplicate_overlap_filter(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = run_chromosort(
                Path(tmp),
                "--split-candidate-min-aligned-bp",
                "20",
                "--split-candidate-min-query-frac",
                "0.20",
            )

            assignments = read_assignments(prefix.with_suffix(".contig_assignments.tsv"))
            candidate = assignments["contigChimeraDup"]
            self.assertEqual(candidate["status"], "kept_split_candidate")
            self.assertEqual(candidate["kept"], "yes")
            self.assertEqual(candidate["split_candidate"], "yes")
            self.assertEqual(candidate["split_candidate_refs"], "chr1,chr2")

            headers = fasta_headers(prefix.with_suffix(".ordered.fa"))
            self.assertIn("chr1_contigChimeraDup", headers)

    def test_paf_input_matches_reference_order_decisions(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = run_chromosort(
                Path(tmp),
                "--split-candidate-min-aligned-bp",
                "20",
                "--split-candidate-min-query-frac",
                "0.20",
                alignment_args=["--paf", str(DATA / "sample.paf")],
            )

            assignments = read_assignments(prefix.with_suffix(".contig_assignments.tsv"))
            self.assertEqual(assignments["contigA"]["status"], "kept")
            self.assertEqual(assignments["contigB"]["status"], "kept")
            self.assertEqual(assignments["contigC"]["status"], "kept_terminal_overlap")
            self.assertEqual(assignments["contigC"]["overlap_class"], "terminal_overlap")
            self.assertEqual(assignments["contigDup"]["status"], "duplicate_overlap")
            self.assertEqual(assignments["contigChimeraDup"]["status"], "kept_split_candidate")
            self.assertEqual(assignments["contigChimeraDup"]["split_candidate_refs"], "chr1,chr2")

            self.assertEqual(
                fasta_headers(prefix.with_suffix(".ordered.fa")),
                [
                    "chr1_contigA",
                    "chr1_contigChimeraDup",
                    "chr1_contigB",
                    "chr2_contigC",
                ],
            )

    def test_graph_assignment_report_summarizes_gfa_context(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = run_graph_sort(Path(tmp))

            graph_rows = read_assignments(prefix.with_suffix(".graph_assignments.tsv"))
            self.assertEqual(graph_rows["left"]["graph_node"], "left")
            self.assertEqual(graph_rows["left"]["graph_node_status"], "present")
            self.assertEqual(graph_rows["left"]["graph_node_has_sequence"], "yes")
            self.assertEqual(graph_rows["left"]["overlap_graph_status"], ".")
            self.assertEqual(graph_rows["left"]["graph_note"], "node_context_only")

            self.assertEqual(graph_rows["bridge_alt"]["graph_self_loop"], "yes")
            self.assertEqual(graph_rows["isolated"]["graph_neighbor_count"], "0")

    def test_low_coverage_split_candidate_is_protected(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ref = tmp_path / "ref.fa"
            assembly = tmp_path / "assembly.fa"
            coords = tmp_path / "sample.coords"
            ref.write_text(">chr1\n" + "A" * 100 + "\n>chr2\n" + "C" * 100 + "\n")
            assembly.write_text(">contigSplit\n" + "G" * 100 + "\n")
            write_coords(
                coords,
                [
                    ("chr1", "contigSplit", 1, 40, 1, 40, 100, 100),
                    ("chr2", "contigSplit", 1, 35, 41, 75, 100, 100),
                ],
            )

            prefix = run_custom_sort(
                tmp_path,
                ref,
                assembly,
                coords,
                "--split-candidate-min-aligned-bp",
                "30",
                "--split-candidate-min-query-frac",
                "0.30",
            )

            assignments = read_assignments(prefix.with_suffix(".contig_assignments.tsv"))
            self.assertEqual(assignments["contigSplit"]["status"], "kept_split_candidate")
            self.assertEqual(assignments["contigSplit"]["kept"], "yes")

    def test_large_alignment_rescues_near_threshold_contig(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ref = tmp_path / "ref.fa"
            assembly = tmp_path / "assembly.fa"
            coords = tmp_path / "sample.coords"
            ref.write_text(">chr1\n" + "A" * 200 + "\n")
            assembly.write_text(">contigLarge\n" + "G" * 100 + "\n")
            write_coords(
                coords,
                [("chr1", "contigLarge", 1, 48, 1, 48, 200, 100)],
            )

            prefix = run_custom_sort(
                tmp_path,
                ref,
                assembly,
                coords,
                "--large-alignment-min-bp",
                "40",
                "--large-alignment-min-query-cov",
                "0.45",
            )

            assignments = read_assignments(prefix.with_suffix(".contig_assignments.tsv"))
            self.assertEqual(assignments["contigLarge"]["status"], "kept_large_alignment")
            self.assertEqual(assignments["contigLarge"]["kept"], "yes")

    def test_span_overlap_and_split_claims_filter_duplicate_fragments(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ref = tmp_path / "ref.fa"
            assembly = tmp_path / "assembly.fa"
            coords = tmp_path / "sample.coords"
            ref.write_text(">chr1\n" + "A" * 220 + "\n>chr2\n" + "C" * 120 + "\n")
            assembly.write_text(
                ">contigAnchor\n" + "A" * 80 + "\n"
                ">contigGapDup\n" + "C" * 21 + "\n"
                ">contigSplit\n" + "G" * 100 + "\n"
                ">contigSplitDup\n" + "T" * 30 + "\n"
            )
            write_coords(
                coords,
                [
                    ("chr1", "contigAnchor", 1, 40, 1, 40, 220, 80),
                    ("chr1", "contigAnchor", 81, 120, 41, 80, 220, 80),
                    ("chr1", "contigGapDup", 50, 70, 1, 21, 220, 21),
                    ("chr2", "contigSplit", 1, 60, 1, 60, 120, 100),
                    ("chr1", "contigSplit", 121, 160, 61, 100, 220, 100),
                    ("chr1", "contigSplitDup", 130, 159, 1, 30, 220, 30),
                ],
            )

            prefix = run_custom_sort(
                tmp_path,
                ref,
                assembly,
                coords,
                "--split-candidate-min-aligned-bp",
                "30",
                "--split-candidate-min-query-frac",
                "0.30",
            )

            assignments = read_assignments(prefix.with_suffix(".contig_assignments.tsv"))
            self.assertEqual(assignments["contigAnchor"]["status"], "kept")
            self.assertEqual(assignments["contigSplit"]["status"], "kept_split_candidate")
            self.assertEqual(assignments["contigGapDup"]["status"], "duplicate_overlap")
            self.assertEqual(assignments["contigGapDup"]["overlap_best_contig"], "contigAnchor")
            self.assertEqual(assignments["contigSplitDup"]["status"], "duplicate_overlap")
            self.assertEqual(assignments["contigSplitDup"]["overlap_best_contig"], "contigSplit")

    def test_terminal_overlap_is_kept_and_reported_separately(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ref = tmp_path / "ref.fa"
            assembly = tmp_path / "assembly.fa"
            coords = tmp_path / "sample.coords"
            ref.write_text(">chr1\n" + "A" * 240 + "\n")
            assembly.write_text(
                ">contigLeft\n" + "A" * 120 + "\n"
                ">contigRight\n" + "C" * 100 + "\n"
            )
            write_coords(
                coords,
                [
                    ("chr1", "contigLeft", 1, 120, 1, 120, 240, 120),
                    ("chr1", "contigRight", 101, 200, 1, 100, 240, 100),
                ],
            )

            prefix = run_custom_sort(tmp_path, ref, assembly, coords)

            assignments = read_assignments(prefix.with_suffix(".contig_assignments.tsv"))
            self.assertEqual(assignments["contigLeft"]["status"], "kept")
            self.assertEqual(assignments["contigRight"]["status"], "kept_terminal_overlap")
            self.assertEqual(assignments["contigRight"]["overlap_class"], "terminal_overlap")
            self.assertEqual(assignments["contigRight"]["terminal_extension_side"], "right")
            self.assertEqual(assignments["contigRight"]["terminal_extension_bp"], "80")

    def test_terminal_extension_rescue_keeps_large_one_sided_extension(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ref = tmp_path / "ref.fa"
            assembly = tmp_path / "assembly.fa"
            coords = tmp_path / "sample.coords"
            ref.write_text(">chr1\n" + "A" * 1200 + "\n")
            assembly.write_text(
                ">contigAnchor\n" + "A" * 1000 + "\n"
                ">contigExtension\n" + "C" * 901 + "\n"
            )
            write_coords(
                coords,
                [
                    ("chr1", "contigAnchor", 1, 1000, 1, 1000, 1200, 1000),
                    ("chr1", "contigExtension", 200, 1100, 1, 901, 1200, 901),
                ],
            )

            default_prefix = run_custom_sort(tmp_path, ref, assembly, coords)
            default_assignments = read_assignments(default_prefix.with_suffix(".contig_assignments.tsv"))
            self.assertEqual(default_assignments["contigExtension"]["status"], "terminal_overlap")
            self.assertEqual(default_assignments["contigExtension"]["kept"], "no")

            rescued_prefix = run_custom_sort(
                tmp_path,
                ref,
                assembly,
                coords,
                "--min-terminal-extension-bp",
                "50",
                "--min-terminal-extension-frac",
                "0.05",
            )
            rescued_assignments = read_assignments(rescued_prefix.with_suffix(".contig_assignments.tsv"))
            self.assertEqual(rescued_assignments["contigExtension"]["status"], "kept_terminal_overlap")
            self.assertEqual(rescued_assignments["contigExtension"]["kept"], "yes")
            self.assertEqual(rescued_assignments["contigExtension"]["terminal_extension_bp"], "100")


if __name__ == "__main__":
    unittest.main()
