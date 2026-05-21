import csv
import json
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


def write_tsv(path, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
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


def review_html_data(path):
    text = path.read_text()
    start_marker = '<script id="chromosort-gapfill-review-data" type="application/json">'
    end_marker = "</script>"
    start = text.index(start_marker) + len(start_marker)
    end = text.index(end_marker, start)
    return json.loads(text[start:end])


def run_gapfill(*args):
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    subprocess.run(
        [sys.executable, "-m", "chromosort.gapfill", *args],
        cwd=ROOT,
        check=True,
        env=env,
    )


def write_unique_gapfill_fixture(tmp_path):
    ordered = tmp_path / "ordered.fa"
    assignments = tmp_path / "assignments.tsv"
    graph = tmp_path / "graph.gfa"

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
    return ordered, assignments, graph


class GapfillTests(unittest.TestCase):
    def test_gapfill_applies_unique_sequence_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ordered, assignments, graph = write_unique_gapfill_fixture(tmp_path)
            prefix = tmp_path / "filled"

            run_gapfill(
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

            plan = read_tsv(Path(str(prefix) + ".gapfill_plan.tsv"))[0]
            self.assertEqual(plan["fill_status"], "fillable")
            self.assertEqual(plan["path_nodes"], "left+,gapper+,right+")
            self.assertEqual(plan["fill_sequence"], "GGGGTT")
            self.assertEqual(plan["right_trim_bp"], "2")
            self.assertEqual(plan["accept_fill"], "no")
            self.assertEqual(plan["applied"], "yes")

            records = read_fasta(Path(str(prefix) + ".gapfilled.fa"))
            self.assertEqual(records["chr1"], "AAAACCGGGGTTCCCC")

    def test_reviewed_plan_controls_graph_gapfill_application(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ordered, assignments, graph = write_unique_gapfill_fixture(tmp_path)
            plan_prefix = tmp_path / "planned"
            accepted_plan = tmp_path / "accepted_plan.tsv"
            rejected_plan = tmp_path / "rejected_plan.tsv"

            run_gapfill(
                "--ordered-fasta",
                str(ordered),
                "--assignments",
                str(assignments),
                "--gfa",
                str(graph),
                "--output-prefix",
                str(plan_prefix),
                "--include-fill-sequences",
            )
            rows = read_tsv(Path(str(plan_prefix) + ".gapfill_plan.tsv"))
            self.assertEqual(rows[0]["fill_status"], "fillable")
            self.assertEqual(rows[0]["accept_fill"], "no")

            rows[0]["accept_fill"] = "yes"
            write_tsv(accepted_plan, rows)
            accepted_prefix = tmp_path / "accepted"
            run_gapfill(
                "--ordered-fasta",
                str(ordered),
                "--assignments",
                str(assignments),
                "--gfa",
                str(graph),
                "--reviewed-plan",
                str(accepted_plan),
                "--output-prefix",
                str(accepted_prefix),
                "--apply",
                "--simple-headers",
            )

            accepted_rows = read_tsv(Path(str(accepted_prefix) + ".gapfill_plan.tsv"))
            self.assertEqual(accepted_rows[0]["accept_fill"], "yes")
            self.assertEqual(accepted_rows[0]["applied"], "yes")
            records = read_fasta(Path(str(accepted_prefix) + ".gapfilled.fa"))
            self.assertEqual(records["chr1"], "AAAACCGGGGTTCCCC")

            rows[0]["accept_fill"] = "no"
            write_tsv(rejected_plan, rows)
            rejected_prefix = tmp_path / "rejected"
            run_gapfill(
                "--ordered-fasta",
                str(ordered),
                "--assignments",
                str(assignments),
                "--gfa",
                str(graph),
                "--reviewed-plan",
                str(rejected_plan),
                "--output-prefix",
                str(rejected_prefix),
                "--apply",
                "--simple-headers",
            )

            rejected_rows = read_tsv(Path(str(rejected_prefix) + ".gapfill_plan.tsv"))
            self.assertEqual(rejected_rows[0]["accept_fill"], "no")
            self.assertEqual(rejected_rows[0]["applied"], "no")
            records = read_fasta(Path(str(rejected_prefix) + ".gapfilled.fa"))
            self.assertEqual(records["chr1"], "AAAACC" + "N" * 14 + "TTCCCC")

    def test_review_html_embeds_gapfill_plan_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ordered, assignments, graph = write_unique_gapfill_fixture(tmp_path)
            prefix = tmp_path / "planned"
            review_html = tmp_path / "gapfill_review.html"

            run_gapfill(
                "--ordered-fasta",
                str(ordered),
                "--assignments",
                str(assignments),
                "--gfa",
                str(graph),
                "--output-prefix",
                str(prefix),
                "--include-fill-sequences",
                "--review-html",
                str(review_html),
            )

            self.assertTrue(review_html.exists())
            self.assertIn("Export reviewed TSV", review_html.read_text())
            data = review_html_data(review_html)
            self.assertEqual(data["schema"], "chromosort-gapfill-review-v1")
            self.assertIn("accept_fill", data["columns"])
            self.assertIn("risk_flags", data["columns"])
            self.assertEqual(len(data["rows"]), 1)
            row = data["rows"][0]
            self.assertEqual(row["accept_fill"], "no")
            self.assertEqual(row["fill_status"], "fillable")
            self.assertEqual(row["path_nodes"], "left+,gapper+,right+")
            self.assertEqual(row["fill_sequence"], "GGGGTT")

    def test_reviewed_plan_rejects_stale_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ordered, assignments, graph = write_unique_gapfill_fixture(tmp_path)
            plan_prefix = tmp_path / "planned"
            reviewed_plan = tmp_path / "reviewed_plan.tsv"

            run_gapfill(
                "--ordered-fasta",
                str(ordered),
                "--assignments",
                str(assignments),
                "--gfa",
                str(graph),
                "--output-prefix",
                str(plan_prefix),
            )
            rows = read_tsv(Path(str(plan_prefix) + ".gapfill_plan.tsv"))
            rows[0]["accept_fill"] = "yes"
            rows[0]["path_nodes"] = "left+,stale+,right+"
            write_tsv(reviewed_plan, rows)

            with self.assertRaises(subprocess.CalledProcessError):
                run_gapfill(
                    "--ordered-fasta",
                    str(ordered),
                    "--assignments",
                    str(assignments),
                    "--gfa",
                    str(graph),
                    "--reviewed-plan",
                    str(reviewed_plan),
                    "--output-prefix",
                    str(tmp_path / "stale"),
                    "--apply",
                )

    def test_gapfill_refuses_ambiguous_graph_paths(self):
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

            run_gapfill(
                "--ordered-fasta",
                str(ordered),
                "--assignments",
                str(assignments),
                "--gfa",
                str(GRAPH_DATA / "unitigs.gfa"),
                "--output-prefix",
                str(prefix),
            )

            plan = read_tsv(Path(str(prefix) + ".gapfill_plan.tsv"))[0]
            self.assertEqual(plan["graph_status"], "ambiguous_paths")
            self.assertEqual(plan["fill_status"], "ambiguous_paths")
            self.assertEqual(plan["candidate_paths"], "2")
            self.assertIn("branching", plan["risk_flags"])
            self.assertIn("self_loop", plan["risk_flags"])
            self.assertIn("left", plan["high_degree_nodes"])
            self.assertIn("bridge_alt", plan["self_loop_nodes"])
            self.assertGreater(int(plan["branch_complexity_score"]), 0)
            self.assertFalse(Path(str(prefix) + ".gapfilled.fa").exists())

    def test_gaf_support_resolves_ambiguous_graph_paths(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ordered = tmp_path / "ordered.fa"
            assignments = tmp_path / "assignments.tsv"
            prefix = tmp_path / "supported"

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

            run_gapfill(
                "--ordered-fasta",
                str(ordered),
                "--assignments",
                str(assignments),
                "--gfa",
                str(GRAPH_DATA / "unitigs.gfa"),
                "--gaf",
                str(GRAPH_DATA / "reads_to_graph.gaf"),
                "--min-gaf-path-support",
                "2",
                "--output-prefix",
                str(prefix),
                "--apply",
                "--include-fill-sequences",
                "--simple-headers",
            )

            plan = read_tsv(Path(str(prefix) + ".gapfill_plan.tsv"))[0]
            self.assertEqual(plan["graph_status"], "gaf_resolved_paths")
            self.assertEqual(plan["fill_status"], "fillable")
            self.assertEqual(plan["gaf_path_support"], "2")
            self.assertEqual(plan["gaf_best_alt_support"], "1")
            self.assertEqual(plan["path_nodes"], "left+,bridge_good+,right+")
            self.assertEqual(plan["fill_sequence"], "GGGG")
            self.assertEqual(plan["right_trim_bp"], "4")
            self.assertEqual(plan["applied"], "yes")

            records = read_fasta(Path(str(prefix) + ".gapfilled.fa"))
            self.assertEqual(records["chr1"], "AAAACCCCGGGGTTTT")

    def test_ref_paf_support_resolves_ambiguous_graph_paths(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ordered = GRAPH_DATA / "gapfill_ordered.fa"
            assignments = GRAPH_DATA / "gapfill_assignments.tsv"
            prefix = tmp_path / "ref_supported"
            review_html = tmp_path / "ref_supported.review.html"

            run_gapfill(
                "--ordered-fasta",
                str(ordered),
                "--assignments",
                str(assignments),
                "--gfa",
                str(GRAPH_DATA / "unitigs.gfa"),
                "--ref-paf",
                str(GRAPH_DATA / "unitig_to_ref.paf"),
                "--min-ref-path-support",
                "7",
                "--output-prefix",
                str(prefix),
                "--apply",
                "--include-fill-sequences",
                "--review-html",
                str(review_html),
                "--simple-headers",
            )

            plan = read_tsv(Path(str(prefix) + ".gapfill_plan.tsv"))[0]
            self.assertEqual(plan["graph_status"], "ref_paf_resolved_paths")
            self.assertEqual(plan["fill_status"], "fillable")
            self.assertEqual(plan["ref_path_support"], "8")
            self.assertEqual(plan["ref_best_alt_support"], "6")
            self.assertIn("branching", plan["risk_flags"])
            self.assertIn("high_degree", plan["risk_flags"])
            self.assertEqual(plan["self_loop_nodes"], ".")
            self.assertEqual(plan["path_nodes"], "left+,bridge_good+,right+")
            self.assertEqual(plan["fill_sequence"], "GGGG")
            self.assertEqual(plan["applied"], "yes")

            records = read_fasta(Path(str(prefix) + ".gapfilled.fa"))
            self.assertEqual(records["chr1"], "AAAACCCCGGGGTTTT")
            self.assertIn("Candidate path comparison", review_html.read_text())
            review_data = review_html_data(review_html)
            self.assertIn("candidateColumns", review_data)
            candidates = review_data["rows"][0]["_candidate_paths"]
            self.assertEqual(len(candidates), 2)
            self.assertEqual(candidates[0]["reported"], "yes")
            self.assertEqual(candidates[0]["path_nodes"], "left+,bridge_good+,right+")
            self.assertEqual(candidates[1]["path_nodes"], "left+,bridge_alt+,right+")
            self.assertEqual(candidates[1]["ref_support"], "6")

    def test_hic_support_resolves_ambiguous_graph_paths(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ordered = tmp_path / "ordered.fa"
            assignments = tmp_path / "assignments.tsv"
            prefix = tmp_path / "hic_supported"

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

            run_gapfill(
                "--ordered-fasta",
                str(ordered),
                "--assignments",
                str(assignments),
                "--gfa",
                str(GRAPH_DATA / "unitigs.gfa"),
                "--hic-pairs",
                str(GRAPH_DATA / "hic_pairs.tsv"),
                "--min-hic-path-support",
                "40",
                "--output-prefix",
                str(prefix),
                "--apply",
                "--include-fill-sequences",
                "--simple-headers",
            )

            plan = read_tsv(Path(str(prefix) + ".gapfill_plan.tsv"))[0]
            self.assertEqual(plan["graph_status"], "hic_resolved_paths")
            self.assertEqual(plan["fill_status"], "fillable")
            self.assertEqual(plan["hic_path_support"], "47")
            self.assertEqual(plan["hic_best_alt_support"], "5")
            self.assertEqual(plan["path_nodes"], "left+,bridge_good+,right+")
            self.assertEqual(plan["fill_sequence"], "GGGG")
            self.assertEqual(plan["applied"], "yes")

            records = read_fasta(Path(str(prefix) + ".gapfilled.fa"))
            self.assertEqual(records["chr1"], "AAAACCCCGGGGTTTT")

    def test_conflicting_gaf_and_hic_support_remains_unresolved(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ordered = tmp_path / "ordered.fa"
            assignments = tmp_path / "assignments.tsv"
            hic_pairs = tmp_path / "conflicting_hic.tsv"
            prefix = tmp_path / "conflict"

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
            hic_pairs.write_text(
                "node_a\tnode_b\tcount\n"
                "left\tbridge_good\t1\n"
                "bridge_good\tright\t1\n"
                "left\tbridge_alt\t30\n"
                "bridge_alt\tright\t25\n"
            )

            run_gapfill(
                "--ordered-fasta",
                str(ordered),
                "--assignments",
                str(assignments),
                "--gfa",
                str(GRAPH_DATA / "unitigs.gfa"),
                "--gaf",
                str(GRAPH_DATA / "reads_to_graph.gaf"),
                "--hic-pairs",
                str(hic_pairs),
                "--min-gaf-path-support",
                "2",
                "--min-hic-path-support",
                "30",
                "--output-prefix",
                str(prefix),
            )

            plan = read_tsv(Path(str(prefix) + ".gapfill_plan.tsv"))[0]
            self.assertEqual(plan["graph_status"], "ambiguous_paths")
            self.assertEqual(plan["fill_status"], "ambiguous_paths")
            self.assertEqual(plan["reason"], "conflicting_gaf_hic_support")
            self.assertEqual(plan["gaf_path_support"], "2")
            self.assertEqual(plan["gaf_best_alt_support"], "1")
            self.assertEqual(plan["hic_path_support"], "2")
            self.assertEqual(plan["hic_best_alt_support"], "55")
            self.assertFalse(Path(str(prefix) + ".gapfilled.fa").exists())


if __name__ == "__main__":
    unittest.main()
