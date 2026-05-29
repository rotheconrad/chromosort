import csv
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().parent / "data" / "chimeric"
SCAFFOLD_DATA = Path(__file__).resolve().parent / "data" / "scaffold"


def run_cli(*args):
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    return subprocess.run(
        [sys.executable, "-m", "chromosort.cli", *args],
        cwd=ROOT,
        check=True,
        env=env,
        text=True,
        capture_output=True,
    )


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fasta_headers(path):
    headers = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                headers.append(line[1:].split()[0])
    return headers


class EvalTests(unittest.TestCase):
    def test_eval_fix_writes_review_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / "sample"
            run_cli(
                "eval",
                "fix",
                "--assembly-fasta",
                str(DATA / "assembly.fa"),
                "--coords",
                str(DATA / "sample.coords"),
                "--contigs",
                "contig_04",
                "--output-prefix",
                str(prefix),
                "--min-segment-bp",
                "5",
                "--mode",
                "sensitive",
            )

            rows = read_tsv(Path(str(prefix) + ".fix_review.tsv"))

        self.assertEqual([row["action"] for row in rows], ["split_piece", "split_piece"])
        self.assertEqual([row["accept"] for row in rows], ["yes", "yes"])
        self.assertEqual(rows[0]["schema"], "chromosort-review-event-v1")
        self.assertEqual(rows[0]["source_contig"], "contig_04")
        self.assertEqual(rows[0]["new_contig"], "chrom02-contig_04-a")
        self.assertEqual(rows[1]["slice_start"], "21")

    def test_fix_can_apply_eval_fix_review_table_without_alignment(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / "sample"
            fixed = Path(tmp) / "fixed.fa"
            report = Path(tmp) / "fixed.tsv"
            run_cli(
                "eval",
                "fix",
                "--assembly-fasta",
                str(DATA / "assembly.fa"),
                "--coords",
                str(DATA / "sample.coords"),
                "--contigs",
                "contig_04",
                "--output-prefix",
                str(prefix),
                "--min-segment-bp",
                "5",
                "--mode",
                "sensitive",
            )
            run_cli(
                "fix",
                "--assembly-fasta",
                str(DATA / "assembly.fa"),
                "--reviewed-plan",
                str(Path(str(prefix) + ".fix_review.tsv")),
                "--output-fasta",
                str(fixed),
                "--report",
                str(report),
            )

            headers = fasta_headers(fixed)
            rows = read_tsv(report)

        self.assertIn("chrom02-contig_04-a", headers)
        self.assertIn("chrom07-contig_04-b", headers)
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["reason"], "Reviewed fix plan accepted 2 split_piece row(s).")

    def test_eval_scaffold_writes_review_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / "sample"
            read_paf = Path(tmp) / "reads.paf"
            read_paf.write_text(
                "read1\t20\t0\t12\t+\tcontigA\t4\t0\t4\t4\t4\t60\n"
                "read1\t20\t13\t20\t+\tcontigB\t3\t0\t3\t3\t3\t60\n"
            )
            run_cli(
                "eval",
                "scaffold",
                "--ordered-fasta",
                str(SCAFFOLD_DATA / "ordered.fa"),
                "--assignments",
                str(SCAFFOLD_DATA / "assignments.tsv"),
                "--read-paf",
                str(read_paf),
                "--read-min-anchor-bp",
                "2",
                "--output-prefix",
                str(prefix),
            )

            rows = read_tsv(Path(str(prefix) + ".scaffold_review.tsv"))

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["task"], "scaffold")
        self.assertEqual(rows[0]["action"], "scaffold_gap")
        self.assertEqual(rows[0]["accept"], "yes")
        self.assertEqual(rows[0]["left_contig"], "chr1_contigA")
        self.assertEqual(rows[0]["right_contig"], "chr1_contigB")
        self.assertEqual(rows[0]["gap_bp"], "5")
        self.assertEqual(rows[0]["longread_bridge_reads"], "1")
        self.assertEqual(rows[0]["longread_read_order_summary"], "left_before_right:1")


if __name__ == "__main__":
    unittest.main()
