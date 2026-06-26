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
MIXED_EVIDENCE_DATA = Path(__file__).resolve().parent / "data" / "mixed_evidence"


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


def write_projected_gapfill_fixture(tmp_path):
    ordered = tmp_path / "projected_ordered.fa"
    assignments = tmp_path / "projected_assignments.tsv"
    graph = tmp_path / "projected_unitigs.noseq.gfa"

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
        "S\tutg_left_a\t*\tLN:i:3\n"
        "S\tutg_left_terminal\t*\tLN:i:3\n"
        "S\tutg_bridge\t*\tLN:i:4\n"
        "S\tutg_right_terminal\t*\tLN:i:3\n"
        "S\tutg_right_b\t*\tLN:i:3\n"
        "L\tutg_left_a\t+\tutg_left_terminal\t+\t0M\n"
        "L\tutg_left_terminal\t+\tutg_bridge\t+\t0M\n"
        "L\tutg_bridge\t+\tutg_right_terminal\t+\t0M\n"
        "L\tutg_right_terminal\t+\tutg_right_b\t+\t0M\n"
        "P\tleft\tutg_left_a+,utg_left_terminal+\t0M\n"
        "P\tright\tutg_right_terminal+,utg_right_b+\t0M\n"
    )
    return ordered, assignments, graph


def write_projected_sequence_gapfill_fixture(tmp_path):
    ordered = tmp_path / "projected_sequence_ordered.fa"
    assignments = tmp_path / "projected_sequence_assignments.tsv"
    graph = tmp_path / "projected_sequence_unitigs.gfa"

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
        "S\tutg_left_body\tAAAA\tLN:i:4\n"
        "S\tutg_left_terminal\tCC\tLN:i:2\n"
        "S\tutg_bridge\tCCGGGGTT\tLN:i:8\n"
        "S\tutg_right_terminal\tTT\tLN:i:2\n"
        "S\tutg_right_body\tCCCC\tLN:i:4\n"
        "L\tutg_left_body\t+\tutg_left_terminal\t+\t0M\n"
        "L\tutg_left_terminal\t+\tutg_bridge\t+\t2M\n"
        "L\tutg_bridge\t+\tutg_right_terminal\t+\t2M\n"
        "L\tutg_right_terminal\t+\tutg_right_body\t+\t0M\n"
        "P\tleft\tutg_left_body+,utg_left_terminal+\t0M\n"
        "P\tright\tutg_right_terminal+,utg_right_body+\t0M\n"
    )
    return ordered, assignments, graph


class EvalTests(unittest.TestCase):
    def test_eval_fix_writes_review_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / "sample"
            gaf = Path(tmp) / "reads.gaf"
            gaf.write_text(
                "read_graph\t40\t0\t40\t+\t>contig_04>graph_neighbor\t40\t0\t40\t40\t40\t60\n"
            )
            run_cli(
                "eval",
                "fix",
                "--assembly-fasta",
                str(DATA / "assembly.fa"),
                "--coords",
                str(DATA / "sample.coords"),
                "--contigs",
                "contig_04",
                "--gaf",
                str(gaf),
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
        self.assertEqual(rows[0]["gaf_node"], "contig_04")
        self.assertEqual(rows[0]["gaf_node_reads"], "1")
        self.assertEqual(rows[0]["gaf_node_orientations"], "+:1")
        self.assertEqual(rows[1]["slice_start"], "21")

    def test_eval_fix_reports_projected_unitig_breakpoint_context(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / "sample"
            graph = Path(tmp) / "unitigs.gfa"
            graph.write_text(
                "S\tutgL\t*\tLN:i:20\n"
                "S\tutgR\t*\tLN:i:20\n"
                "L\tutgL\t+\tutgR\t+\t0M\n"
                "P\tcontig_04\tutgL+,utgR+\t0M\n"
            )
            run_cli(
                "eval",
                "fix",
                "--assembly-fasta",
                str(DATA / "assembly.fa"),
                "--coords",
                str(DATA / "sample.coords"),
                "--contigs",
                "contig_04",
                "--gfa",
                str(graph),
                "--output-prefix",
                str(prefix),
                "--min-segment-bp",
                "5",
                "--mode",
                "sensitive",
            )

            rows = read_tsv(Path(str(prefix) + ".fix_review.tsv"))

        self.assertEqual(rows[0]["graph_node_status"], "missing_node")
        self.assertEqual(rows[0]["graph_unitig"], "utgL")
        self.assertEqual(rows[0]["graph_unitig_orientation"], "+")
        self.assertEqual(rows[0]["graph_unitig_boundary_distance_bp"], "1")
        self.assertEqual(rows[0]["graph_near_path_boundary"], "yes")

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

    def test_eval_scaffold_reports_mixed_gfa_paf_gaf_branch_support(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / "mixed"
            run_cli(
                "eval",
                "scaffold",
                "--ordered-fasta",
                str(SCAFFOLD_DATA / "ordered.fa"),
                "--assignments",
                str(SCAFFOLD_DATA / "assignments.tsv"),
                "--gfa",
                str(MIXED_EVIDENCE_DATA / "branch.gfa"),
                "--gaf",
                str(MIXED_EVIDENCE_DATA / "reads_to_graph.gaf"),
                "--read-paf",
                str(MIXED_EVIDENCE_DATA / "reads_to_assembly.paf"),
                "--read-min-anchor-bp",
                "2",
                "--output-prefix",
                str(prefix),
            )

            rows = read_tsv(Path(str(prefix) + ".scaffold_review.tsv"))

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["task"], "scaffold")
        self.assertEqual(rows[0]["graph_path_nodes"], "contigA+,bridge_good+,contigB+")
        self.assertEqual(rows[0]["gaf_path_nodes"], "contigA+,bridge_good+,contigB+")
        self.assertEqual(rows[0]["gaf_path_support"], "1")
        self.assertEqual(rows[0]["gaf_best_alt_path_nodes"], "contigA+,bridge_alt+,contigB+")
        self.assertEqual(rows[0]["gaf_best_alt_support"], "2")
        self.assertEqual(rows[0]["gaf_support_status"], "supports_alternate")
        self.assertEqual(rows[0]["gaf_selected_reads"], "read_good")
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
        self.assertEqual(rows[0]["gaf_support_status"], ".")
        self.assertEqual(rows[0]["longread_bridge_reads"], "1")
        self.assertEqual(records["chr1"], "AAAACCGGGGTTCCCC")

    def test_eval_gapfill_reports_external_patch_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ordered, assignments, graph = write_gapfill_fixture(tmp_path)
            patch_table = tmp_path / "patches.tsv"
            prefix = tmp_path / "eval_patch"
            patch_table.write_text(
                "scaffold\tleft_contig\tright_contig\tpatch_id\tsource\tpatch_sequence\n"
                "chr1\tchr1_left\tchr1_right\tpatch_eval\tLR_Gapcloser\tGGGGTT\n"
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
                "--patch-table",
                str(patch_table),
                "--include-fill-sequences",
                "--output-prefix",
                str(prefix),
            )
            rows = read_tsv(Path(str(prefix) + ".gapfill_review.tsv"))

        self.assertEqual(rows[0]["patch_candidate_count"], "1")
        self.assertEqual(rows[0]["patch_best_id"], "patch_eval")
        self.assertEqual(rows[0]["patch_best_source"], "LR_Gapcloser")
        self.assertEqual(rows[0]["patch_graph_status"], "exact_graph_match")
        self.assertEqual(rows[0]["patch_best_sequence"], "GGGGTT")

    def test_eval_gapfill_accepts_scaffold_fasta_with_agp(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            _ordered, _assignments, graph = write_gapfill_fixture(tmp_path)
            scaffold = tmp_path / "scaffold.fa"
            agp = tmp_path / "scaffold.agp"
            eval_prefix = tmp_path / "eval_scaffold"
            applied_prefix = tmp_path / "applied_scaffold"

            scaffold.write_text(">chr1\nAAAACC" + "N" * 14 + "TTCCCC\n")
            agp.write_text(
                "##agp-version\t2.1\n"
                "chr1\t1\t6\t1\tW\tleft\t1\t6\t+\n"
                "chr1\t7\t20\t2\tN\t14\tscaffold\tyes\talign_genus\n"
                "chr1\t21\t26\t3\tW\tright\t1\t6\t+\n"
            )

            run_cli(
                "eval",
                "gapfill",
                "--scaffold-fasta",
                str(scaffold),
                "--agp",
                str(agp),
                "--gfa",
                str(graph),
                "--include-fill-sequences",
                "--output-prefix",
                str(eval_prefix),
            )
            review = Path(str(eval_prefix) + ".gapfill_review.tsv")
            rows = read_tsv(review)

            run_cli(
                "gapfill",
                "--scaffold-fasta",
                str(scaffold),
                "--agp",
                str(agp),
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
        self.assertEqual(rows[0]["left_contig"], "left")
        self.assertEqual(rows[0]["right_contig"], "right")
        self.assertEqual(rows[0]["accept"], "yes")
        self.assertEqual(rows[0]["fill_status"], "fillable")
        self.assertEqual(rows[0]["path_nodes"], "left+,gapper+,right+")
        self.assertEqual(records["chr1"], "AAAACCGGGGTTCCCC")

    def test_eval_gapfill_reports_projected_unitig_plan(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ordered, assignments, graph = write_projected_gapfill_fixture(tmp_path)
            prefix = tmp_path / "eval_projected"

            run_cli(
                "eval",
                "gapfill",
                "--ordered-fasta",
                str(ordered),
                "--assignments",
                str(assignments),
                "--gfa",
                str(graph),
                "--project-gfa-paths",
                "--output-prefix",
                str(prefix),
            )
            rows = read_tsv(Path(str(prefix) + ".gapfill_review.tsv"))

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["task"], "gapfill")
        self.assertEqual(rows[0]["accept"], "no")
        self.assertEqual(rows[0]["left_graph_node"], "utg_left_terminal")
        self.assertEqual(rows[0]["right_graph_node"], "utg_right_terminal")
        self.assertEqual(rows[0]["fill_status"], "projected_path_planning_only")
        self.assertEqual(
            rows[0]["path_nodes"],
            "utg_left_terminal+,utg_bridge+,utg_right_terminal+",
        )

    def test_eval_gapfill_projected_sequence_review_can_drive_apply(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ordered, assignments, graph = write_projected_sequence_gapfill_fixture(tmp_path)
            eval_prefix = tmp_path / "eval_projected_sequence"
            applied_prefix = tmp_path / "applied_projected_sequence"

            run_cli(
                "eval",
                "gapfill",
                "--ordered-fasta",
                str(ordered),
                "--assignments",
                str(assignments),
                "--gfa",
                str(graph),
                "--project-gfa-paths",
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
                "--project-gfa-paths",
                "--reviewed-plan",
                str(review),
                "--output-prefix",
                str(applied_prefix),
                "--apply",
                "--simple-headers",
            )
            records = read_fasta(Path(str(applied_prefix) + ".gapfilled.fa"))

        self.assertEqual(rows[0]["accept"], "yes")
        self.assertEqual(rows[0]["fill_status"], "fillable")
        self.assertEqual(rows[0]["fill_sequence"], "GGGGTT")
        self.assertEqual(
            rows[0]["path_nodes"],
            "utg_left_terminal+,utg_bridge+,utg_right_terminal+",
        )
        self.assertEqual(records["chr1"], "AAAACCGGGGTTCCCC")


if __name__ == "__main__":
    unittest.main()
