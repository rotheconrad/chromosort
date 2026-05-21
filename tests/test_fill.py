import csv
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GRAPH_DATA = Path(__file__).resolve().parent / "data" / "graph_gotchas"


def write_assignments(path, rows):
    header = [
        "contig",
        "kept",
        "new_name",
        "assigned_ref",
        "order_in_ref",
        "ref_start",
        "ref_end",
        "orientation",
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


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


def run_fill(*args):
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    subprocess.run(
        [sys.executable, "-m", "chromosort.fill", *args],
        cwd=ROOT,
        check=True,
        env=env,
    )


class FillTests(unittest.TestCase):
    def test_fill_applies_unique_sequence_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ordered = tmp_path / "ordered.fa"
            assignments = tmp_path / "assignments.tsv"
            graph = tmp_path / "graph.gfa"
            prefix = tmp_path / "filled"

            ordered.write_text(">chr1_left\nAAAACC\n>chr1_right\nTTCCCC\n")
            write_assignments(
                assignments,
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
                ],
            )
            graph.write_text(
                "H\tVN:Z:1.0\n"
                "S\tleft\tAAAACC\tLN:i:6\n"
                "S\tgapper\tCCGGGGTT\tLN:i:8\n"
                "S\tright\tTTCCCC\tLN:i:6\n"
                "L\tleft\t+\tgapper\t+\t2M\n"
                "L\tgapper\t+\tright\t+\t2M\n"
            )

            run_fill(
                "--ordered-fasta",
                str(ordered),
                "--assignments",
                str(assignments),
                "--gfa",
                str(graph),
                "--output-prefix",
                str(prefix),
                "--apply",
                "--include-fill-sequences",
                "--simple-headers",
            )

            plan = read_tsv(Path(str(prefix) + ".fill_plan.tsv"))[0]
            self.assertEqual(plan["fill_status"], "fillable")
            self.assertEqual(plan["path_nodes"], "left+,gapper+,right+")
            self.assertEqual(plan["fill_sequence"], "GGGGTT")
            self.assertEqual(plan["right_trim_bp"], "2")
            self.assertEqual(plan["applied"], "yes")

            records = read_fasta(Path(str(prefix) + ".filled.fa"))
            self.assertEqual(records["chr1"], "AAAACCGGGGTTCCCC")

    def test_fill_refuses_ambiguous_graph_paths(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ordered = tmp_path / "ordered.fa"
            assignments = tmp_path / "assignments.tsv"
            prefix = tmp_path / "plan"

            ordered.write_text(">chr1_left\nAAAACCCC\n>chr1_right\nGGGGTTTT\n")
            write_assignments(
                assignments,
                [
                    {
                        "contig": "left",
                        "kept": "yes",
                        "new_name": "chr1_left",
                        "assigned_ref": "chr1",
                        "order_in_ref": 1,
                        "ref_start": 1,
                        "ref_end": 8,
                        "orientation": "+",
                    },
                    {
                        "contig": "right",
                        "kept": "yes",
                        "new_name": "chr1_right",
                        "assigned_ref": "chr1",
                        "order_in_ref": 2,
                        "ref_start": 37,
                        "ref_end": 44,
                        "orientation": "+",
                    },
                ],
            )

            run_fill(
                "--ordered-fasta",
                str(ordered),
                "--assignments",
                str(assignments),
                "--gfa",
                str(GRAPH_DATA / "unitigs.gfa"),
                "--output-prefix",
                str(prefix),
            )

            plan = read_tsv(Path(str(prefix) + ".fill_plan.tsv"))[0]
            self.assertEqual(plan["graph_status"], "ambiguous_paths")
            self.assertEqual(plan["fill_status"], "ambiguous_paths")
            self.assertEqual(plan["candidate_paths"], "2")
            self.assertFalse(Path(str(prefix) + ".filled.fa").exists())


if __name__ == "__main__":
    unittest.main()
