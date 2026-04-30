import csv
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().parent / "data" / "scaffold"


def run_scaffold(tmp_path, *extra_args):
    prefix = tmp_path / "sample"
    cmd = [
        sys.executable,
        "-m",
        "chromosort.scaffold",
        "--ordered-fasta",
        str(DATA / "ordered.fa"),
        "--assignments",
        str(DATA / "assignments.tsv"),
        "--output-prefix",
        str(prefix),
        "--simple-headers",
        *extra_args,
    ]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    subprocess.run(cmd, cwd=ROOT, check=True, env=env)
    return prefix


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


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class ScaffoldTests(unittest.TestCase):
    def test_scaffold_uses_inferred_reference_gap_by_default(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = run_scaffold(Path(tmp))

            records = read_fasta(Path(str(prefix) + ".scaffold.fa"))
            self.assertEqual(records["chr1"], "AAAANNNNNTTT")
            self.assertEqual(records["chr2"], "GG")

            gaps = read_tsv(Path(str(prefix) + ".scaffold_gaps.tsv"))
            self.assertEqual(len(gaps), 1)
            self.assertEqual(gaps[0]["raw_inferred_gap_bp"], "5")
            self.assertEqual(gaps[0]["gap_bp"], "5")
            self.assertEqual(gaps[0]["gap_mode"], "inferred")

            summary = {
                row["scaffold"]: row
                for row in read_tsv(Path(str(prefix) + ".scaffold_summary.tsv"))
            }
            self.assertEqual(summary["chr1"]["contigs"], "2")
            self.assertEqual(summary["chr1"]["gap_bp"], "5")
            self.assertEqual(summary["chr1"]["scaffold_bp"], "12")

    def test_scaffold_can_use_fixed_gap(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = run_scaffold(Path(tmp), "--fixed-gap-bp", "2")

            records = read_fasta(Path(str(prefix) + ".scaffold.fa"))
            self.assertEqual(records["chr1"], "AAAANNTTT")

            gaps = read_tsv(Path(str(prefix) + ".scaffold_gaps.tsv"))
            self.assertEqual(gaps[0]["raw_inferred_gap_bp"], "5")
            self.assertEqual(gaps[0]["gap_bp"], "2")
            self.assertEqual(gaps[0]["gap_mode"], "fixed")


if __name__ == "__main__":
    unittest.main()
