import csv
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().parent / "data" / "chimeric"


def run_fix_contigs(tmp_path, *extra_args):
    output_fasta = tmp_path / "fixed.fa"
    report = tmp_path / "fixed.tsv"
    cmd = [
        sys.executable,
        "-m",
        "chromosort.fix_contigs",
        "--assembly-fasta",
        str(DATA / "assembly.fa"),
        "--coords",
        str(DATA / "sample.coords"),
        "--contigs",
        "contig_04",
        "contig_12",
        "--output-fasta",
        str(output_fasta),
        "--report",
        str(report),
        "--min-segment-bp",
        "5",
        *extra_args,
    ]
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


if __name__ == "__main__":
    unittest.main()

