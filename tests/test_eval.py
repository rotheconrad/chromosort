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


def write_gapfill_fixture(tmp_path):
    ordered = tmp_path / "ordered.fa"
    assignments = tmp_path / "assignments.tsv"
    graph = tmp_path / "graph.gfa"

    ordered.write_text(">chr1_left\nAAAACC\n>chr1_right\nTTCCCC\n")
    with open(assignments, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "contig",
                "kept",
                "new_name",
                "assigned_ref",
                "order_in_ref",
                "ref_start",
                "ref_end",
                "orientation",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(
            [
                {
                    "contig": "left",
                    "kept": "yes",
                    "new_name": "chr1_left",
                    "assigned_ref": "chr1",
                    "order_in_ref": 1,
                    "ref_start": 1,
                    "ref_end": 6,
                    "orientation": "+",
                },
                {
                    "contig": "right",
                    "kept": "yes",
                    "new_name": "chr1_right",
                    "assigned_ref": "chr1",
                    "order_in_ref": 2,
                    "ref_start": 21,
                    "ref_end": 26,
                    "orientation": "+",
                },
            ]
        )
    graph.write_text(
        "H\tVN:Z:1.0\n"
        "S\tleft\tAAAACC\tLN:i:6\n"
        "S\tgapper\tCCGGGGTT\tLN:i:8\n"
        "S\tright\tTTCCCC\tLN:i:6\n"
        "L\tleft\t+\tgapper\t+\t2M\n"
        "L\tgapper\t+\tright\t+\t2M\n"
    )
    return ordered, assignments, graph


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
            graph = Path(tmp) / "graph.gfa"
            gaf = Path(tmp) / "reads.gaf"
            read_paf.write_text(
                "read1\t20\t0\t12\t+\tcontigA\t4\t0\t4\t4\t4\t60\n"
                "read1\t20\t13\t20\t+\tcontigB\t3\t0\t3\t3\t3\t60\n"
            )
            graph.write_text(
                "H\tVN:Z:1.0\n"
                "S\tcontigA\tAAAA\tLN:i:4\n"
                "S\tbridge\tNN\tLN:i:2\n"
                "S\tcontigB\tTTT\tLN:i:3\n"
                "L\tcontigA\t+\tbridge\t+\t0M\n"
                "L\tbridge\t+\tcontigB\t+\t0M\n"
            )
            gaf.write_text(
                "read_graph\t20\t0\t20\t+\t>contigA>bridge>contigB\t20\t0\t20\t20\t20\t60\n"
            )
            run_cli(
                "eval",
                "scaffold",
                "--ordered-fasta",
                str(SCAFFOLD_DATA / "ordered.fa"),
                "--assignments",
                str(SCAFFOLD_DATA / "assignments.tsv"),
                "--gfa",
                str(graph),
                "--gaf",
                str(gaf),
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
        self.assertEqual(rows[0]["gaf_path_nodes"], "contigA+,bridge+,contigB+")
        self.assertEqual(rows[0]["gaf_path_support"], "1")
        self.assertEqual(rows[0]["gaf_support_status"], "supports_selected")
        self.assertEqual(rows[0]["longread_bridge_reads"], "1")
        self.assertEqual(rows[0]["longread_read_order_summary"], "left_before_right:1")

    def test_eval_gapfill_review_table_can_drive_gapfill_apply(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ordered, assignments, graph = write_gapfill_fixture(tmp_path)
            read_paf = tmp_path / "reads.paf"
            eval_prefix = tmp_path / "eval"
            applied_prefix = tmp_path / "applied"
            read_paf.write_text(
                "read1\t20\t0\t12\t+\tleft\t6\t0\t6\t6\t6\t60\n"
                "read1\t20\t13\t20\t+\tright\t6\t0\t6\t6\t6\t60\n"
            )

            run_cli(
                "eval",
                "gapfill",
                "--ordered-fasta",
                str(ordered),
                "--assignments",
                str(assignments),
                "--gfa",
                str(graph),
                "--read-paf",
                str(read_paf),
                "--read-min-anchor-bp",
                "2",
                "--include-fill-sequences",
                "--output-prefix",
                str(eval_prefix),
            )
            review = Path(str(eval_prefix) + ".gapfill_review.tsv")
            rows = read_tsv(review)

            run_cli(
                "gapfill",
                "--ordered-fasta",
                str(ordered),
                "--assignments",
                str(assignments),
                "--gfa",
                str(graph),
                "--reviewed-plan",
                str(review),
                "--output-prefix",
                str(applied_prefix),
                "--apply",
                "--simple-headers",
            )

            records = read_fasta(Path(str(applied_prefix) + ".gapfilled.fa"))

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["task"], "gapfill")
        self.assertEqual(rows[0]["action"], "fill_path")
        self.assertEqual(rows[0]["accept"], "yes")
        self.assertEqual(rows[0]["fill_status"], "fillable")
        self.assertEqual(rows[0]["path_nodes"], "left+,gapper+,right+")
        self.assertEqual(rows[0]["fill_sequence"], "GGGGTT")
        self.assertEqual(rows[0]["longread_bridge_reads"], "1")
        self.assertEqual(records["chr1"], "AAAACCGGGGTTCCCC")


if __name__ == "__main__":
    unittest.main()
