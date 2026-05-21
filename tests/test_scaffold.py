import csv
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().parent / "data" / "scaffold"
GRAPH_DATA = Path(__file__).resolve().parent / "data" / "graph_gotchas"


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


def run_custom_scaffold(tmp_path, ordered_fasta, assignments, *extra_args):
    prefix = tmp_path / "custom"
    cmd = [
        sys.executable,
        "-m",
        "chromosort.scaffold",
        "--ordered-fasta",
        str(ordered_fasta),
        "--assignments",
        str(assignments),
        "--output-prefix",
        str(prefix),
        "--simple-headers",
        *extra_args,
    ]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    subprocess.run(cmd, cwd=ROOT, check=True, env=env)
    return prefix


def write_assignment_table(path, rows):
    header = [
        "contig",
        "kept",
        "new_name",
        "assigned_ref",
        "order_in_ref",
        "ref_start",
        "ref_end",
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=header, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_graph_assignment_table(path, rows):
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
        writer.writerows(rows)


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

    def test_scaffold_graph_report_marks_direct_edges(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ordered = tmp_path / "ordered.fa"
            assignments = tmp_path / "assignments.tsv"
            ordered.write_text(
                ">chr1_left\nAAAACCCC\n"
                ">chr1_bridge_good\nCCCCGGGG\n"
                ">chr1_right\nGGGGTTTT\n"
            )
            write_graph_assignment_table(
                assignments,
                [
                    {
                        "contig": "left",
                        "kept": "yes",
                        "new_name": "chr1_left",
                        "assigned_ref": "chr1",
                        "order_in_ref": "1",
                        "ref_start": "1",
                        "ref_end": "8",
                        "orientation": "+",
                    },
                    {
                        "contig": "bridge_good",
                        "kept": "yes",
                        "new_name": "chr1_bridge_good",
                        "assigned_ref": "chr1",
                        "order_in_ref": "2",
                        "ref_start": "19",
                        "ref_end": "26",
                        "orientation": "+",
                    },
                    {
                        "contig": "right",
                        "kept": "yes",
                        "new_name": "chr1_right",
                        "assigned_ref": "chr1",
                        "order_in_ref": "3",
                        "ref_start": "37",
                        "ref_end": "44",
                        "orientation": "+",
                    },
                ],
            )

            prefix = run_custom_scaffold(
                tmp_path,
                ordered,
                assignments,
                "--gfa",
                str(GRAPH_DATA / "unitigs.gfa"),
            )
            rows = read_tsv(Path(str(prefix) + ".graph_gaps.tsv"))

            self.assertEqual([row["graph_status"] for row in rows], ["direct_edge", "direct_edge"])
            self.assertEqual(rows[0]["direct_edge"], "yes")
            self.assertEqual(rows[0]["direct_edge_orientations"], "left+>bridge_good+")
            self.assertEqual(rows[0]["direct_edge_overlap_bp"], "4")
            self.assertEqual(rows[0]["shortest_path_edges"], "1")
            self.assertEqual(rows[0]["shortest_path_nodes"], "left+,bridge_good+")

    def test_scaffold_graph_report_marks_short_paths_and_missing_nodes(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ordered = tmp_path / "ordered.fa"
            assignments = tmp_path / "assignments.tsv"
            ordered.write_text(
                ">chr1_left\nAAAACCCC\n"
                ">chr1_right\nGGGGTTTT\n"
                ">chr1_missing\nAAAA\n"
            )
            write_graph_assignment_table(
                assignments,
                [
                    {
                        "contig": "left",
                        "kept": "yes",
                        "new_name": "chr1_left",
                        "assigned_ref": "chr1",
                        "order_in_ref": "1",
                        "ref_start": "1",
                        "ref_end": "8",
                        "orientation": "+",
                    },
                    {
                        "contig": "right",
                        "kept": "yes",
                        "new_name": "chr1_right",
                        "assigned_ref": "chr1",
                        "order_in_ref": "2",
                        "ref_start": "37",
                        "ref_end": "44",
                        "orientation": "+",
                    },
                    {
                        "contig": "missing",
                        "kept": "yes",
                        "new_name": "chr1_missing",
                        "assigned_ref": "chr1",
                        "order_in_ref": "3",
                        "ref_start": "50",
                        "ref_end": "53",
                        "orientation": "+",
                    },
                ],
            )

            prefix = run_custom_scaffold(
                tmp_path,
                ordered,
                assignments,
                "--gfa",
                str(GRAPH_DATA / "unitigs.gfa"),
            )
            rows = read_tsv(Path(str(prefix) + ".graph_gaps.tsv"))

            self.assertEqual(rows[0]["graph_status"], "short_path")
            self.assertEqual(rows[0]["direct_edge"], "no")
            self.assertEqual(rows[0]["shortest_path_edges"], "2")
            self.assertEqual(rows[0]["shortest_path_nodes"], "left+,bridge_good+,right+")
            self.assertEqual(rows[0]["candidate_intermediate_nodes"], "bridge_good+")
            self.assertEqual(rows[1]["graph_status"], "missing_right_node")
            self.assertEqual(rows[1]["right_graph_node"], ".")

    def test_scaffold_reports_negative_reference_gap_as_overlap(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ordered = tmp_path / "ordered.fa"
            assignments = tmp_path / "assignments.tsv"
            ordered.write_text(">chr1_left\nAAAACCCC\n>chr1_right\nCCCCGGGG\n")
            write_assignment_table(
                assignments,
                [
                    {
                        "contig": "left",
                        "kept": "yes",
                        "new_name": "chr1_left",
                        "assigned_ref": "chr1",
                        "order_in_ref": "1",
                        "ref_start": "1",
                        "ref_end": "8",
                    },
                    {
                        "contig": "right",
                        "kept": "yes",
                        "new_name": "chr1_right",
                        "assigned_ref": "chr1",
                        "order_in_ref": "2",
                        "ref_start": "5",
                        "ref_end": "12",
                    },
                ],
            )

            prefix = run_custom_scaffold(tmp_path, ordered, assignments)

            records = read_fasta(Path(str(prefix) + ".scaffold.fa"))
            self.assertEqual(records["chr1"], "AAAACCCCCCCCGGGG")

            gaps = read_tsv(Path(str(prefix) + ".scaffold_gaps.tsv"))
            self.assertEqual(gaps[0]["raw_inferred_gap_bp"], "-4")
            self.assertEqual(gaps[0]["gap_bp"], "0")
            self.assertEqual(gaps[0]["overlap_bp"], "4")
            self.assertEqual(gaps[0]["overlap_class"], "terminal_overlap")
            self.assertEqual(gaps[0]["overlap_action"], "zero_gap")

            summary = {
                row["scaffold"]: row
                for row in read_tsv(Path(str(prefix) + ".scaffold_summary.tsv"))
            }
            self.assertEqual(summary["chr1"]["overlap_gaps"], "1")
            self.assertEqual(summary["chr1"]["overlap_bp"], "4")
            self.assertEqual(summary["chr1"]["trimmed_bp"], "0")

    def test_scaffold_can_trim_reference_inferred_terminal_overlap(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ordered = tmp_path / "ordered.fa"
            assignments = tmp_path / "assignments.tsv"
            ordered.write_text(">chr1_left\nAAAACCCC\n>chr1_right\nCCCCGGGG\n")
            write_assignment_table(
                assignments,
                [
                    {
                        "contig": "left",
                        "kept": "yes",
                        "new_name": "chr1_left",
                        "assigned_ref": "chr1",
                        "order_in_ref": "1",
                        "ref_start": "1",
                        "ref_end": "8",
                    },
                    {
                        "contig": "right",
                        "kept": "yes",
                        "new_name": "chr1_right",
                        "assigned_ref": "chr1",
                        "order_in_ref": "2",
                        "ref_start": "5",
                        "ref_end": "12",
                    },
                ],
            )

            prefix = run_custom_scaffold(
                tmp_path,
                ordered,
                assignments,
                "--overlap-policy",
                "trim-reference",
            )

            records = read_fasta(Path(str(prefix) + ".scaffold.fa"))
            self.assertEqual(records["chr1"], "AAAACCCCGGGG")

            gaps = read_tsv(Path(str(prefix) + ".scaffold_gaps.tsv"))
            self.assertEqual(gaps[0]["trimmed_bp"], "4")
            self.assertEqual(gaps[0]["overlap_action"], "trimmed_reference")

    def test_scaffold_graph_confirm_can_trim_direct_terminal_overlap(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ordered = tmp_path / "ordered.fa"
            assignments = tmp_path / "assignments.tsv"
            ordered.write_text(">chr1_left\nAAAACCCC\n>chr1_bridge_good\nCCCCGGGG\n")
            write_graph_assignment_table(
                assignments,
                [
                    {
                        "contig": "left",
                        "kept": "yes",
                        "new_name": "chr1_left",
                        "assigned_ref": "chr1",
                        "order_in_ref": "1",
                        "ref_start": "1",
                        "ref_end": "8",
                        "orientation": "+",
                    },
                    {
                        "contig": "bridge_good",
                        "kept": "yes",
                        "new_name": "chr1_bridge_good",
                        "assigned_ref": "chr1",
                        "order_in_ref": "2",
                        "ref_start": "5",
                        "ref_end": "12",
                        "orientation": "+",
                    },
                ],
            )

            prefix = run_custom_scaffold(
                tmp_path,
                ordered,
                assignments,
                "--gfa",
                str(GRAPH_DATA / "unitigs.gfa"),
                "--graph-overlap-policy",
                "confirm",
            )

            records = read_fasta(Path(str(prefix) + ".scaffold.fa"))
            self.assertEqual(records["chr1"], "AAAACCCCGGGG")

            gaps = read_tsv(Path(str(prefix) + ".scaffold_gaps.tsv"))
            self.assertEqual(gaps[0]["graph_overlap_policy"], "confirm")
            self.assertEqual(gaps[0]["overlap_action"], "graph_confirmed_trim_reference")
            self.assertEqual(gaps[0]["graph_overlap_action"], "direct_edge:left+>bridge_good+")
            self.assertEqual(gaps[0]["trimmed_bp"], "4")

    def test_scaffold_can_trim_sequence_confirmed_terminal_overlap(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ordered = tmp_path / "ordered.fa"
            assignments = tmp_path / "assignments.tsv"
            ordered.write_text(">chr1_left\nAAAACCCC\n>chr1_right\nCCCCGGGG\n")
            write_assignment_table(
                assignments,
                [
                    {
                        "contig": "left",
                        "kept": "yes",
                        "new_name": "chr1_left",
                        "assigned_ref": "chr1",
                        "order_in_ref": "1",
                        "ref_start": "1",
                        "ref_end": "8",
                    },
                    {
                        "contig": "right",
                        "kept": "yes",
                        "new_name": "chr1_right",
                        "assigned_ref": "chr1",
                        "order_in_ref": "2",
                        "ref_start": "5",
                        "ref_end": "12",
                    },
                ],
            )

            prefix = run_custom_scaffold(
                tmp_path,
                ordered,
                assignments,
                "--overlap-policy",
                "trim-sequence",
            )

            records = read_fasta(Path(str(prefix) + ".scaffold.fa"))
            self.assertEqual(records["chr1"], "AAAACCCCGGGG")

            gaps = read_tsv(Path(str(prefix) + ".scaffold_gaps.tsv"))
            self.assertEqual(gaps[0]["trimmed_bp"], "4")
            self.assertEqual(gaps[0]["overlap_action"], "trimmed_sequence")
            self.assertEqual(gaps[0]["sequence_overlap_identity"], "1.000")


if __name__ == "__main__":
    unittest.main()
