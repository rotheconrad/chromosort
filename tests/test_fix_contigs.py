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
                "--simple-headers",
                "--min-segment-bp",
                "10",
                "--auto-breakpoint-penalty-bp",
                "50",
                "--auto-min-piece-aligned-bp",
                "50",
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
