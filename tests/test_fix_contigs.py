import csv
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().parent / "data" / "chimeric"
NOISY_DATA = Path(__file__).resolve().parent / "data" / "noisy_fix"


def run_fix_contigs(tmp_path, *extra_args, contigs=None, data=DATA):
    output_fasta = tmp_path / "fixed.fa"
    report = tmp_path / "fixed.tsv"
    if contigs is None:
        contigs = ["contig_04", "contig_12"]
    cmd = [
        sys.executable,
        "-m",
        "chromosort.fix_contigs",
        "--assembly-fasta",
        str(data / "assembly.fa"),
        "--coords",
        str(data / "sample.coords"),
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

    def test_auto_rejects_candidates_with_too_many_breakpoints(self):
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
                "--auto",
                "--auto-breakpoint-penalty-bp",
                "1",
                "--auto-min-piece-aligned-bp",
                "5",
                "--auto-min-piece-query-frac",
                "0",
                "--simple-headers",
                data=data,
                contigs=[],
            )

            self.assertEqual(list(read_fasta(output_fasta)), ["contig_many_transitions"])
            rows = read_report(report)
            self.assertEqual(rows[0]["status"], "not_split_auto_too_many_breakpoints")

    def test_auto_default_rejects_tiny_reference_transition_piece(self):
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
                "--auto",
                "--auto-breakpoint-penalty-bp",
                "1",
                "--auto-min-piece-aligned-bp",
                "5",
                "--simple-headers",
                data=data,
                contigs=[],
            )

            self.assertEqual(list(read_fasta(output_fasta)), ["contig_small_head"])
            rows = read_report(report)
            self.assertEqual(rows[0]["status"], "not_split_auto_smooth")
            self.assertIn("failed support thresholds", rows[0]["reason"])
            self.assertIn("query span fraction 0.0400", rows[0]["reason"])

    def test_auto_default_ignores_same_reference_inversions(self):
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
                "--auto",
                "--simple-headers",
                data=data,
                contigs=[],
            )

            self.assertEqual(list(read_fasta(output_fasta)), ["contig_inversion_only"])
            self.assertEqual(read_report(report), [])

    def test_auto_default_splits_complex_same_reference_inversions(self):
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
                "--auto",
                "--simple-headers",
                "--auto-breakpoint-penalty-bp",
                "5",
                "--auto-min-piece-aligned-bp",
                "10",
                "--auto-min-piece-query-frac",
                "0",
                "--auto-complex-inversion-min-piece-aligned-bp",
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

    def test_splits_user_nominated_chimeric_contigs(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_fasta, report = run_fix_contigs(Path(tmp))

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

    def test_pieces_only_omits_untargeted_contigs(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_fasta, _ = run_fix_contigs(Path(tmp), "--pieces-only", "--simple-headers")

            records = read_fasta(output_fasta)
            self.assertNotIn("contig_01", records)
            self.assertEqual(len(records), 5)

    def test_splits_middle_and_terminal_inversion_blocks(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_fasta, report = run_fix_contigs(
                Path(tmp),
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

    def test_auto_detects_chromosome_and_orientation_transitions(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_fasta, report = run_fix_contigs(
                Path(tmp),
                "--auto",
                "--auto-sensitive",
                "--auto-max-breakpoints",
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

    def test_auto_smoothing_splits_large_events_and_ignores_sv_noise(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_fasta, report = run_fix_contigs(
                Path(tmp),
                "--auto",
                "--auto-split-inversions",
                "--simple-headers",
                "--min-segment-bp",
                "10",
                "--auto-breakpoint-penalty-bp",
                "50",
                "--auto-min-piece-aligned-bp",
                "50",
                "--auto-max-breakpoints",
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
                self.assertEqual(rows_by_contig[contig][0]["status"], "not_split_auto_smooth")

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


if __name__ == "__main__":
    unittest.main()
