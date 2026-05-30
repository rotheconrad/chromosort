import csv
import io
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = Path(__file__).resolve().parent / "data" / "chimeric"
NOISY_DATA = Path(__file__).resolve().parent / "data" / "noisy_fix"


def run_fix_contigs(tmp_path, *extra_args, contigs=None, data=DATA, alignment_args=None):
    output_fasta = tmp_path / "fixed.fa"
    report = tmp_path / "fixed.tsv"
    if contigs is None:
        contigs = ["contig_04", "contig_12"]
    if alignment_args is None:
        alignment_args = ["--coords", str(data / "sample.coords")]
    cmd = [
        sys.executable,
        "-m",
        "chromosort.fix_contigs",
        "--assembly-fasta",
        str(data / "assembly.fa"),
        *alignment_args,
        "--output-fasta",
        str(output_fasta),
        "--report",
        str(report),
        "--min-segment-bp",
        "5",
        *extra_args,
    ]
    if contigs:
        cmd[cmd.index("--output-fasta"):cmd.index("--output-fasta")] = ["--contigs", *contigs]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    subprocess.run(cmd, cwd=ROOT, check=True, env=env)
    return output_fasta, report


def read_fasta(path):
    records = {}
    name = None
    parts = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
        if name is not None:
            records[name] = "".join(parts)
    return records


def read_report(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class FixContigsTests(unittest.TestCase):
    def test_cli_reports_progress_to_stderr(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            output_fasta = tmp_path / "fixed.fa"
            report = tmp_path / "fixed.tsv"
            env = os.environ.copy()
            env["PYTHONPATH"] = str(ROOT / "src")

            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "chromosort.cli",
                    "fix",
                    "--assembly-fasta",
                    str(DATA / "assembly.fa"),
                    "--coords",
                    str(DATA / "sample.coords"),
                    "--contigs",
                    "contig_04",
                    "--output-fasta",
                    str(output_fasta),
                    "--report",
                    str(report),
                    "--min-segment-bp",
                    "5",
                ],
                cwd=ROOT,
                check=True,
                env=env,
                text=True,
                capture_output=True,
            )

            self.assertIn("Starting chromo fix:", result.stderr)
            self.assertIn("Reading coords alignments:", result.stderr)
            self.assertIn("Finished alignment scan:", result.stderr)
            self.assertIn("Writing fixed FASTA:", result.stderr)
            self.assertIn("Processed 1 requested contigs.", result.stderr)

    def test_graph_guard_warns_on_simple_split_target(self):
        from chromosort.fix_contigs import ContigPlan, graph_guard_fix_warnings
        from chromosort.graph import read_gfa

        graph_data = Path(__file__).resolve().parent / "data" / "graph_gotchas"
        graph = read_gfa(graph_data / "unitigs.gfa")
        plans = {
            "isolated": ContigPlan(
                contig="isolated",
                status="split",
                pieces=[],
                reason="test plan",
            )
        }
        stream = io.StringIO()

        graph_guard_fix_warnings(["isolated"], plans, graph, stream)

        warning = stream.getvalue()
        self.assertIn("graph guard", warning)
        self.assertIn("isolated", warning)
        self.assertIn("simple neighborhood", warning)

    def test_user_split_collapses_same_target_alignment_gaps(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            data = tmp_path / "data"
            data.mkdir()
            (data / "assembly.fa").write_text(
                ">contig_gap_chimera\n"
                + "A" * 20
                + "N" * 10
                + "A" * 20
                + "C" * 20
                + "N" * 10
                + "G" * 20
                + "\n"
            )
            (data / "sample.coords").write_text(
                "/tmp/ref.fa /tmp/assembly.fa\n"
                "NUCMER\n\n"
                "    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |  [LEN R]  [LEN Q]  |  [COV R]  [COV Q]  | [TAGS]\n"
                "===============================================================================================================================\n"
                "       1       20  |       1       20  |       20       20  |   100.00  |      200      100  |    10.00    20.00  | chrom01\tcontig_gap_chimera\n"
                "      31       50  |      31       50  |       20       20  |   100.00  |      200      100  |    10.00    20.00  | chrom01\tcontig_gap_chimera\n"
                "       1       20  |      51       70  |       20       20  |   100.00  |      200      100  |    10.00    20.00  | chrom07\tcontig_gap_chimera\n"
                "      81      100  |      81      100  |       20       20  |   100.00  |      200      100  |    10.00    20.00  | chrom01\tcontig_gap_chimera\n"
            )

            output_fasta, report = run_fix_contigs(
                tmp_path,
                "--mode",
                "sensitive",
                "--simple-headers",
                contigs=["contig_gap_chimera"],
                data=data,
            )

            self.assertEqual(
                list(read_fasta(output_fasta)),
                [
                    "chrom01-contig_gap_chimera-a",
                    "chrom07-contig_gap_chimera-b",
                    "chrom01-contig_gap_chimera-c",
                ],
            )
            rows = read_report(report)
            self.assertEqual(
                [(row["new_contig"], row["slice_start"], row["slice_end"]) for row in rows],
                [
                    ("chrom01-contig_gap_chimera-a", "1", "50"),
                    ("chrom07-contig_gap_chimera-b", "51", "75"),
                    ("chrom01-contig_gap_chimera-c", "76", "100"),
                ],
            )

    def test_all_rejects_candidates_with_too_many_breakpoints_per_contig(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            data = tmp_path / "data"
            data.mkdir()
            (data / "assembly.fa").write_text(
                ">contig_many_transitions\n" + "A" * 60 + "\n"
            )
            rows = [
                (
                    ref,
                    start,
                    end,
                )
                for ref, start, end in [
                    ("chrom01", 1, 10),
                    ("chrom02", 11, 20),
                    ("chrom03", 21, 30),
                    ("chrom04", 31, 40),
                    ("chrom05", 41, 50),
                    ("chrom06", 51, 60),
                ]
            ]
            coords = [
                "/tmp/ref.fa /tmp/assembly.fa",
                "NUCMER",
                "",
                "    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |  [LEN R]  [LEN Q]  |  [COV R]  [COV Q]  | [TAGS]",
                "===============================================================================================================================",
            ]
            for ref, start, end in rows:
                coords.append(
                    f"       1       10  |      {start:2d}      {end:2d}  |       10       10  |   100.00  |      100       60  |    10.00    16.67  | {ref}\tcontig_many_transitions"
                )
            (data / "sample.coords").write_text("\n".join(coords) + "\n")

            output_fasta, report = run_fix_contigs(
                tmp_path,
                "--all",
                "--breakpoint-penalty-bp",
                "1",
                "--min-piece-aligned-bp",
                "5",
                "--min-piece-query-frac",
                "0",
                "--simple-headers",
                data=data,
                contigs=[],
            )

            self.assertEqual(list(read_fasta(output_fasta)), ["contig_many_transitions"])
            rows = read_report(report)
            self.assertEqual(rows[0]["status"], "not_split_too_many_breakpoints")

    def test_default_mode_rejects_tiny_reference_transition_piece(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            data = tmp_path / "data"
            data.mkdir()
            (data / "assembly.fa").write_text(
                ">contig_small_head\n" + "A" * 200 + "\n"
            )
            (data / "sample.coords").write_text(
                "/tmp/ref.fa /tmp/assembly.fa\n"
                "NUCMER\n\n"
                "    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |  [LEN R]  [LEN Q]  |  [COV R]  [COV Q]  | [TAGS]\n"
                "===============================================================================================================================\n"
                "       1        8  |       1        8  |        8        8  |   100.00  |      300      200  |     2.67     4.00  | chrom02\tcontig_small_head\n"
                "       1      192  |       9      200  |      192      192  |   100.00  |      300      200  |    64.00    96.00  | chrom19\tcontig_small_head\n"
            )

            output_fasta, report = run_fix_contigs(
                tmp_path,
                "--all",
                "--breakpoint-penalty-bp",
                "1",
                "--min-piece-aligned-bp",
                "5",
                "--simple-headers",
                data=data,
                contigs=[],
            )

            self.assertEqual(list(read_fasta(output_fasta)), ["contig_small_head"])
            rows = read_report(report)
            self.assertEqual(rows[0]["status"], "not_split_smooth")
            self.assertIn("failed support thresholds", rows[0]["reason"])
            self.assertIn("query span fraction 0.0400", rows[0]["reason"])

    def test_default_mode_merges_tiny_terminal_noise_before_large_split(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            data = tmp_path / "data"
            data.mkdir()
            (data / "assembly.fa").write_text(
                ">contig_tiny_head_then_split\n" + "A" * 300 + "\n"
            )
            (data / "sample.coords").write_text(
                "/tmp/ref.fa /tmp/assembly.fa\n"
                "NUCMER\n\n"
                "    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |  [LEN R]  [LEN Q]  |  [COV R]  [COV Q]  | [TAGS]\n"
                "===============================================================================================================================\n"
                "       1        8  |       1        8  |        8        8  |   100.00  |      300      300  |     2.67     2.67  | chrom02\tcontig_tiny_head_then_split\n"
                "       1      142  |       9      150  |      142      142  |   100.00  |      300      300  |    47.33    47.33  | chrom19\tcontig_tiny_head_then_split\n"
                "       1      150  |     151      300  |      150      150  |   100.00  |      300      300  |    50.00    50.00  | chrom20\tcontig_tiny_head_then_split\n"
            )

            output_fasta, report = run_fix_contigs(
                tmp_path,
                "--all",
                "--breakpoint-penalty-bp",
                "1",
                "--min-piece-aligned-bp",
                "5",
                "--simple-headers",
                data=data,
                contigs=[],
            )

            self.assertEqual(
                list(read_fasta(output_fasta)),
                [
                    "chrom19-contig_tiny_head_then_split-a",
                    "chrom20-contig_tiny_head_then_split-b",
                ],
            )
            rows = read_report(report)
            self.assertEqual([row["dominant_ref"] for row in rows], ["chrom19", "chrom20"])
            self.assertTrue(all(row["status"] == "split" for row in rows))

    def test_default_mode_ignores_same_reference_inversions(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            data = tmp_path / "data"
            data.mkdir()
            (data / "assembly.fa").write_text(
                ">contig_inversion_only\n" + "A" * 60 + "\n"
            )
            (data / "sample.coords").write_text(
                "/tmp/ref.fa /tmp/assembly.fa\n"
                "NUCMER\n\n"
                "    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |  [LEN R]  [LEN Q]  |  [COV R]  [COV Q]  | [TAGS]\n"
                "===============================================================================================================================\n"
                "       1       20  |       1       20  |       20       20  |   100.00  |      100       60  |    20.00    33.33  | chrom01\tcontig_inversion_only\n"
                "      21       40  |      40       21  |       20       20  |   100.00  |      100       60  |    20.00    33.33  | chrom01\tcontig_inversion_only\n"
                "      41       60  |      41       60  |       20       20  |   100.00  |      100       60  |    20.00    33.33  | chrom01\tcontig_inversion_only\n"
            )

            output_fasta, report = run_fix_contigs(
                tmp_path,
                "--all",
                "--simple-headers",
                data=data,
                contigs=[],
            )

            self.assertEqual(list(read_fasta(output_fasta)), ["contig_inversion_only"])
            self.assertEqual(read_report(report), [])

    def test_default_mode_splits_complex_same_reference_inversions(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            data = tmp_path / "data"
            data.mkdir()
            (data / "assembly.fa").write_text(
                ">contig_complex_inversion\n" + "A" * 90 + "\n"
            )
            (data / "sample.coords").write_text(
                "/tmp/ref.fa /tmp/assembly.fa\n"
                "NUCMER\n\n"
                "    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |  [LEN R]  [LEN Q]  |  [COV R]  [COV Q]  | [TAGS]\n"
                "===============================================================================================================================\n"
                "       1       80  |       1       30  |       80       30  |   100.00  |      100       90  |    80.00    33.33  | chrom01\tcontig_complex_inversion\n"
                "      30       60  |      60       31  |       31       30  |   100.00  |      100       90  |    31.00    33.33  | chrom01\tcontig_complex_inversion\n"
                "      50       90  |      61       90  |       41       30  |   100.00  |      100       90  |    41.00    33.33  | chrom01\tcontig_complex_inversion\n"
            )

            output_fasta, report = run_fix_contigs(
                tmp_path,
                "--all",
                "--simple-headers",
                "--breakpoint-penalty-bp",
                "5",
                "--min-piece-aligned-bp",
                "10",
                "--min-piece-query-frac",
                "0",
                "--complex-inversion-min-piece-aligned-bp",
                "10",
                data=data,
                contigs=[],
            )

            self.assertEqual(
                list(read_fasta(output_fasta)),
                [
                    "chrom01-contig_complex_inversion-a",
                    "chrom01-contig_complex_inversion-b",
                    "chrom01-contig_complex_inversion-c",
                ],
            )
            rows = read_report(report)
            self.assertEqual(len(rows), 3)
            self.assertTrue(all(row["status"] == "split" for row in rows))

    def test_sensitive_mode_splits_selected_chimeric_contigs(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_fasta, report = run_fix_contigs(Path(tmp), "--mode", "sensitive")

            records = read_fasta(output_fasta)
            self.assertEqual(
                list(records),
                [
                    "contig_01",
                    "chrom02-contig_04-a",
                    "chrom07-contig_04-b",
                    "chrom04-contig_12-a",
                    "chrom05-contig_12-b",
                    "chrom04-contig_12-c",
                    "contig_inv_mid",
                    "contig_inv_end",
                ],
            )
            self.assertEqual(records["contig_01"], "NNNNNNNNNN")
            self.assertEqual(records["chrom02-contig_04-a"], "A" * 20)
            self.assertEqual(records["chrom07-contig_04-b"], "C" * 20)
            self.assertEqual(records["chrom04-contig_12-a"], "G" * 5)
            self.assertEqual(records["chrom05-contig_12-b"], "T" * 30)
            self.assertEqual(records["chrom04-contig_12-c"], "G" * 5)

            rows = read_report(report)
            self.assertEqual(len(rows), 5)
            self.assertTrue(all(row["status"] == "split" for row in rows))
            self.assertEqual(
                [(row["new_contig"], row["slice_start"], row["slice_end"]) for row in rows],
                [
                    ("chrom02-contig_04-a", "1", "20"),
                    ("chrom07-contig_04-b", "21", "40"),
                    ("chrom04-contig_12-a", "1", "5"),
                    ("chrom05-contig_12-b", "6", "35"),
                    ("chrom04-contig_12-c", "36", "40"),
                ],
            )

    def test_paf_input_splits_selected_chimeric_contigs(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_fasta, report = run_fix_contigs(
                Path(tmp),
                "--mode",
                "sensitive",
                alignment_args=["--paf", str(DATA / "sample.paf")],
            )

            self.assertEqual(
                list(read_fasta(output_fasta)),
                [
                    "contig_01",
                    "chrom02-contig_04-a",
                    "chrom07-contig_04-b",
                    "chrom04-contig_12-a",
                    "chrom05-contig_12-b",
                    "chrom04-contig_12-c",
                    "contig_inv_mid",
                    "contig_inv_end",
                ],
            )
            rows = read_report(report)
            self.assertEqual(
                [(row["new_contig"], row["slice_start"], row["slice_end"]) for row in rows],
                [
                    ("chrom02-contig_04-a", "1", "20"),
                    ("chrom07-contig_04-b", "21", "40"),
                    ("chrom04-contig_12-a", "1", "5"),
                    ("chrom05-contig_12-b", "6", "35"),
                    ("chrom04-contig_12-c", "36", "40"),
                ],
            )

    def test_graph_report_summarizes_requested_contigs(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            gfa = tmp_path / "assembly.gfa"
            gfa.write_text(
                "H\tVN:Z:1.0\n"
                "S\tcontig_04\t*\tLN:i:40\n"
                "S\tcontig_12\t*\tLN:i:40\n"
                "L\tcontig_04\t+\tcontig_12\t+\t0M\n"
            )

            _, report = run_fix_contigs(
                tmp_path,
                "--mode",
                "sensitive",
                "--gfa",
                str(gfa),
            )

            rows = {
                row["original_contig"]: row
                for row in read_report(report.with_suffix(".graph.tsv"))
            }
            self.assertEqual(rows["contig_04"]["status"], "split")
            self.assertEqual(rows["contig_04"]["graph_node_status"], "present")
            self.assertEqual(rows["contig_04"]["graph_note"], "split_source_present_in_graph_review")
            self.assertIn("chrom02-contig_04-a", rows["contig_04"]["split_pieces"])
            self.assertEqual(rows["contig_12"]["graph_out_degree"], "0")

    def test_pieces_only_omits_untargeted_contigs(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_fasta, _ = run_fix_contigs(
                Path(tmp),
                "--mode",
                "sensitive",
                "--pieces-only",
                "--simple-headers",
            )

            records = read_fasta(output_fasta)
            self.assertNotIn("contig_01", records)
            self.assertEqual(len(records), 5)

    def test_splits_middle_and_terminal_inversion_blocks(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_fasta, report = run_fix_contigs(
                Path(tmp),
                "--mode",
                "sensitive",
                "--pieces-only",
                "--simple-headers",
                contigs=["contig_inv_mid", "contig_inv_end"],
            )

            records = read_fasta(output_fasta)
            self.assertEqual(
                list(records),
                [
                    "chrom06-contig_inv_mid-a",
                    "chrom06-contig_inv_mid-b",
                    "chrom06-contig_inv_mid-c",
                    "chrom08-contig_inv_end-a",
                    "chrom08-contig_inv_end-b",
                ],
            )
            rows = read_report(report)
            self.assertEqual(
                [(row["new_contig"], row["orientation"]) for row in rows],
                [
                    ("chrom06-contig_inv_mid-a", "+"),
                    ("chrom06-contig_inv_mid-b", "-"),
                    ("chrom06-contig_inv_mid-c", "+"),
                    ("chrom08-contig_inv_end-a", "+"),
                    ("chrom08-contig_inv_end-b", "-"),
                ],
            )

    def test_all_sensitive_mode_detects_chromosome_and_orientation_transitions(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_fasta, report = run_fix_contigs(
                Path(tmp),
                "--all",
                "--mode",
                "sensitive",
                "--max-breakpoints-per-contig",
                "-1",
                "--simple-headers",
                contigs=[],
            )

            records = read_fasta(output_fasta)
            self.assertEqual(
                list(records),
                [
                    "contig_01",
                    "chrom02-contig_04-a",
                    "chrom07-contig_04-b",
                    "chrom04-contig_12-a",
                    "chrom05-contig_12-b",
                    "chrom04-contig_12-c",
                    "chrom06-contig_inv_mid-a",
                    "chrom06-contig_inv_mid-b",
                    "chrom06-contig_inv_mid-c",
                    "chrom08-contig_inv_end-a",
                    "chrom08-contig_inv_end-b",
                ],
            )
            rows = read_report(report)
            self.assertEqual({row["original_contig"] for row in rows}, {
                "contig_04",
                "contig_12",
                "contig_inv_mid",
                "contig_inv_end",
            })

    def test_comprehensive_mode_splits_large_events_and_ignores_sv_noise(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_fasta, report = run_fix_contigs(
                Path(tmp),
                "--all",
                "--mode",
                "comprehensive",
                "--simple-headers",
                "--min-segment-bp",
                "10",
                "--breakpoint-penalty-bp",
                "50",
                "--min-piece-aligned-bp",
                "50",
                "--max-breakpoints-per-contig",
                "-1",
                data=NOISY_DATA,
                contigs=[],
            )

            records = read_fasta(output_fasta)
            self.assertEqual(
                list(records),
                [
                    "contig_indel_noise",
                    "contig_small_inv_sv",
                    "contig_repeat_noise",
                    "chrom02-contig_true_chimera-a",
                    "chrom07-contig_true_chimera-b",
                    "chrom06-contig_true_inv_mid-a",
                    "chrom06-contig_true_inv_mid-b",
                    "chrom06-contig_true_inv_mid-c",
                    "chrom10-contig_complex_chimera_indels-a",
                    "chrom11-contig_complex_chimera_indels-b",
                    "contig_terminal_inv_noise",
                    "chrom13-contig_true_inv_end-a",
                    "chrom13-contig_true_inv_end-b",
                ],
            )

            rows = read_report(report)
            rows_by_contig = {}
            for row in rows:
                rows_by_contig.setdefault(row["original_contig"], []).append(row)

            self.assertNotIn("contig_indel_noise", rows_by_contig)
            for contig in [
                "contig_small_inv_sv",
                "contig_repeat_noise",
                "contig_terminal_inv_noise",
            ]:
                self.assertEqual(rows_by_contig[contig][0]["status"], "not_split_smooth")

            self.assertEqual(
                [row["new_contig"] for row in rows_by_contig["contig_true_chimera"]],
                ["chrom02-contig_true_chimera-a", "chrom07-contig_true_chimera-b"],
            )
            self.assertEqual(
                [(row["new_contig"], row["orientation"]) for row in rows_by_contig["contig_true_inv_mid"]],
                [
                    ("chrom06-contig_true_inv_mid-a", "+"),
                    ("chrom06-contig_true_inv_mid-b", "-"),
                    ("chrom06-contig_true_inv_mid-c", "+"),
                ],
            )
            self.assertEqual(
                [row["new_contig"] for row in rows_by_contig["contig_complex_chimera_indels"]],
                [
                    "chrom10-contig_complex_chimera_indels-a",
                    "chrom11-contig_complex_chimera_indels-b",
                ],
            )
            self.assertEqual(
                [(row["new_contig"], row["orientation"]) for row in rows_by_contig["contig_true_inv_end"]],
                [
                    ("chrom13-contig_true_inv_end-a", "+"),
                    ("chrom13-contig_true_inv_end-b", "-"),
                ],
            )

    def test_contig_subset_uses_same_planner_as_all_scope(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            all_dir = tmp_path / "all"
            subset_dir = tmp_path / "subset"
            all_dir.mkdir()
            subset_dir.mkdir()
            common_args = [
                "--mode",
                "comprehensive",
                "--simple-headers",
                "--min-segment-bp",
                "10",
                "--breakpoint-penalty-bp",
                "50",
                "--min-piece-aligned-bp",
                "50",
                "--max-breakpoints-per-contig",
                "-1",
            ]
            targets = [
                "contig_true_chimera",
                "contig_true_inv_mid",
                "contig_small_inv_sv",
            ]

            _, all_report = run_fix_contigs(
                all_dir,
                "--all",
                *common_args,
                data=NOISY_DATA,
                contigs=[],
            )
            _, subset_report = run_fix_contigs(
                subset_dir,
                *common_args,
                data=NOISY_DATA,
                contigs=targets,
            )

            def selected_rows(path):
                return sorted(
                    [
                        (
                            row["original_contig"],
                            row["status"],
                            row["new_contig"],
                            row["slice_start"],
                            row["slice_end"],
                            row["orientation"],
                        )
                        for row in read_report(path)
                        if row["original_contig"] in targets
                    ]
                )

            self.assertEqual(selected_rows(subset_report), selected_rows(all_report))


if __name__ == "__main__":
    unittest.main()
