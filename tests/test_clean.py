import csv
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data" / "chimeric"


def run_clean(tmp_path, *extra_args):
    prefix = tmp_path / "sample"
    cmd = [
        sys.executable,
        "-m",
        "chromosort.clean",
        "--ref-fasta",
        str(DATA / "ref.fa"),
        "--assembly-fasta",
        str(DATA / "assembly.fa"),
        "--coords",
        str(DATA / "sample.coords"),
        "--output-prefix",
        str(prefix),
        "--min-aligned-bp",
        "1",
        "--min-query-cov",
        "0.1",
        "--min-novel-ref-bp",
        "1",
        "--split-candidate-min-aligned-bp",
        "1",
        "--split-candidate-min-query-frac",
        "0.01",
        "--min-segment-bp",
        "1",
        "--breakpoint-penalty-bp",
        "0",
        "--min-piece-aligned-bp",
        "1",
        "--min-piece-query-frac",
        "0",
        "--complex-inversion-min-piece-aligned-bp",
        "1",
        "--orient-to-reference",
        *extra_args,
    ]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    subprocess.run(cmd, cwd=ROOT, check=True, env=env)
    return prefix


def read_fasta_headers(path):
    headers = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                headers.append(line[1:].strip().split()[0])
    return headers


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class CleanTests(unittest.TestCase):
    def test_clean_discards_unaligned_and_splits_retained_chimeras(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = run_clean(Path(tmp))

            headers = read_fasta_headers(Path(str(prefix) + ".clean.fa"))
            self.assertNotIn("contig_01", headers)
            self.assertIn("chrom02_contig_04_a", headers)
            self.assertIn("chrom07_contig_04_b", headers)
            self.assertIn("chrom04_contig_12_a", headers)
            self.assertIn("chrom05_contig_12_b", headers)
            self.assertEqual(headers.index("chrom02_contig_04_a"), 0)

            rows = read_tsv(Path(str(prefix) + ".clean_contigs.tsv"))
            by_contig = {}
            for row in rows:
                by_contig.setdefault(row["source_contig"], []).append(row)

            self.assertEqual(by_contig["contig_01"][0]["clean_status"], "discarded_no_alignment")
            self.assertEqual(
                [row["clean_status"] for row in by_contig["contig_04"]],
                ["kept_split_piece", "kept_split_piece"],
            )
            self.assertEqual(by_contig["contig_inv_mid"][0]["fix_status"], "not_split_single_target")

            summary = Path(str(prefix) + ".run_summary.txt").read_text()
            self.assertIn("clean_status_kept_split_piece\t5", summary)
            self.assertIn("re-run MUMmer or minimap2", summary)

    def test_fix_scope_split_candidates_leaves_other_retained_contigs_unsplit(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = run_clean(Path(tmp), "--fix-scope", "split-candidates")
            rows = read_tsv(Path(str(prefix) + ".clean_contigs.tsv"))
            inv_mid = [row for row in rows if row["source_contig"] == "contig_inv_mid"][0]

            self.assertEqual(inv_mid["clean_status"], "kept_unsplit")
            self.assertEqual(inv_mid["fix_selected"], "no")
            targets = Path(str(prefix) + ".fix_targets.txt").read_text().splitlines()
            self.assertEqual(targets, ["contig_04", "contig_12"])


if __name__ == "__main__":
    unittest.main()
