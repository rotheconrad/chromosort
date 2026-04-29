import csv
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().parent / "data"


def run_chromosort(tmp_path, *extra_args):
    prefix = tmp_path / "sample"
    cmd = [
        sys.executable,
        "-m",
        "chromosort.reference_order",
        "--ref-fasta",
        str(DATA / "ref.fa"),
        "--assembly-fasta",
        str(DATA / "assembly.fa"),
        "--coords",
        str(DATA / "sample.coords"),
        "--output-prefix",
        str(prefix),
        "--min-aligned-bp",
        "10",
        "--min-query-cov",
        "0.10",
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


class ReferenceOrderTests(unittest.TestCase):
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

    def test_no_overlap_filter_keeps_otherwise_good_duplicate(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = run_chromosort(Path(tmp), "--no-overlap-filter")

            assignments = read_assignments(prefix.with_suffix(".contig_assignments.tsv"))
            self.assertEqual(assignments["contigDup"]["status"], "kept")

            headers = fasta_headers(prefix.with_suffix(".ordered.fa"))
            self.assertIn("chr1_contigDup", headers)


if __name__ == "__main__":
    unittest.main()
