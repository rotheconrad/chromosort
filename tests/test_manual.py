import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().parent / "data"
GRAPH_DATA = DATA / "graph_gotchas"


def run_manual(
    tmp_path,
    *extra_args,
    alignment_args=None,
    output_name="manual.html",
    ref_fasta=None,
    assembly_fasta=None,
):
    output_html = tmp_path / output_name
    if alignment_args is None:
        alignment_args = ["--paf", str(DATA / "sample.paf")]
    if ref_fasta is None:
        ref_fasta = DATA / "ref.fa"
    if assembly_fasta is None:
        assembly_fasta = DATA / "assembly.fa"
    cmd = [
        sys.executable,
        "-m",
        "chromosort.manual",
        "--ref-fasta",
        str(ref_fasta),
        "--assembly-fasta",
        str(assembly_fasta),
        *alignment_args,
        "--output-html",
        str(output_html),
        *extra_args,
    ]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    subprocess.run(cmd, cwd=ROOT, check=True, env=env)
    return output_html


def run_manual_apply(tmp_path, assembly_fasta, recipe, *extra_args):
    recipe_path = tmp_path / "recipe.json"
    output_fasta = tmp_path / "manual.fa"
    report = tmp_path / "manual.tsv"
    recipe_path.write_text(json.dumps(recipe) + "\n")
    cmd = [
        sys.executable,
        "-m",
        "chromosort.manual",
        "apply",
        "--assembly-fasta",
        str(assembly_fasta),
        "--recipe",
        str(recipe_path),
        "--output-fasta",
        str(output_fasta),
        "--report",
        str(report),
        *extra_args,
    ]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    subprocess.run(cmd, cwd=ROOT, check=True, env=env)
    return output_fasta, report


def run_manual_mode(tmp_path, mode, *extra_args):
    output_html = tmp_path / f"{mode}.html"
    cmd = [
        sys.executable,
        "-m",
        "chromosort.manual",
        mode,
        "--ref-fasta",
        str(DATA / "ref.fa"),
        "--assembly-fasta",
        str(DATA / "assembly.fa"),
        "--paf",
        str(DATA / "sample.paf"),
        "--output-html",
        str(output_html),
        *extra_args,
    ]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    subprocess.run(cmd, cwd=ROOT, check=True, env=env)
    return output_html


def html_data(path):
    text = path.read_text()
    marker = "window.CHROMOSORT_MANUAL_DATA = "
    start = text.index(marker) + len(marker)
    end = text.index(";\n</script>", start)
    return json.loads(text[start:end])


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


def read_agp_rows(path):
    with open(path) as fh:
        return [
            line.rstrip("\n").split("\t")
            for line in fh
            if line.strip() and not line.startswith("#")
        ]


class ManualTests(unittest.TestCase):
    def test_manual_dashboard_orders_all_contigs_without_filtering(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_html = run_manual(Path(tmp), "--suggested-output-fasta", "sample.manual.fa")
            text = output_html.read_text()
            data = html_data(output_html)

            self.assertIn("ChromoSort Manual", text)
            self.assertIn("Remove unaligned", text)
            self.assertIn("Load the assembly FASTA to export edited sequence.", text)
            self.assertFalse(data["settings"]["embedSequences"])
            self.assertEqual(data["stats"]["assemblyContigs"], 6)
            self.assertEqual(data["stats"]["unalignedContigs"], 1)
            self.assertEqual(data["stats"]["embeddedSequenceCount"], 0)
            self.assertEqual(
                [piece["source"] for piece in data["initialPieces"]],
                [
                    "contigA",
                    "contigChimeraDup",
                    "contigDup",
                    "contigB",
                    "contigC",
                    "contigNo",
                ],
            )
            self.assertFalse(data["initialPieces"][-1]["aligned"])
            self.assertFalse(any(piece["removed"] for piece in data["initialPieces"]))

    def test_manual_dashboard_can_embed_sequences_for_single_file_export(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_html = run_manual(Path(tmp), "--embed-sequences")
            data = html_data(output_html)

            self.assertTrue(data["settings"]["embedSequences"])
            self.assertEqual(data["stats"]["embeddedSequenceCount"], 6)
            self.assertIn("contigA", data["sequences"])

    def test_manual_dashboard_embeds_graph_context(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_html = run_manual(
                Path(tmp),
                "--gfa",
                str(GRAPH_DATA / "unitigs.gfa"),
                alignment_args=["--paf", str(GRAPH_DATA / "unitig_to_ref.paf")],
                ref_fasta=GRAPH_DATA / "ref.fa",
                assembly_fasta=GRAPH_DATA / "assembly.fa",
            )
            text = output_html.read_text()
            data = html_data(output_html)

            self.assertIn("Graph present", text)
            self.assertIn("Graph filter", text)
            self.assertIn("graphPanel", text)
            self.assertIn("same ref", text)
            self.assertEqual(data["inputs"]["gfa"], str(GRAPH_DATA / "unitigs.gfa"))
            self.assertEqual(data["stats"]["graphContigsPresent"], 8)
            self.assertEqual(data["stats"]["graphBranchingContigs"], 5)

            graph_by_name = {
                item["name"]: item["graph"]
                for item in data["queryRecords"]
            }
            self.assertEqual(graph_by_name["left"]["graphNode"], "left")
            self.assertEqual(graph_by_name["left"]["graphComplexity"], "branching")
            self.assertEqual(graph_by_name["left"]["graphCoverage"], 20)
            self.assertEqual(graph_by_name["left"]["graphOutDegree"], 3)
            self.assertEqual(
                sorted(edge["otherNode"] for edge in graph_by_name["left"]["outgoing"]),
                ["bridge_alt", "bridge_good", "reverse_only"],
            )
            left_outgoing = {
                edge["otherNode"]: edge
                for edge in graph_by_name["left"]["outgoing"]
            }
            self.assertTrue(left_outgoing["bridge_good"]["sameBestRef"])
            self.assertEqual(left_outgoing["bridge_good"]["otherBestRef"], "chr1")
            self.assertTrue(left_outgoing["bridge_alt"]["sameBestRef"])
            self.assertTrue(left_outgoing["reverse_only"]["otherAligned"])
            self.assertEqual(graph_by_name["isolated"]["graphComplexity"], "simple")
            self.assertTrue(graph_by_name["bridge_alt"]["graphSelfLoop"])

    def test_manual_fix_mode_embeds_review_event_queue(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            review = tmp_path / "fix_review.tsv"
            review.write_text(
                "schema\tevent_id\ttask\taction\ttarget\taccept\tstatus\tconfidence\treason\tnotes\tsource_contig\tnew_contig\n"
                "chromosort-review-event-v1\tfix:contigA:1\tfix\tsplit_piece\tcontigA\tyes\tcandidate\t.\ttest\t\tcontigA\tcontigA_left\n"
            )
            output_html = run_manual_mode(tmp_path, "fix", "--review-table", str(review))
            text = output_html.read_text()
            data = html_data(output_html)

            self.assertIn("eventPanel", text)
            self.assertEqual(data["mode"], "manual-fix")
            self.assertEqual(data["reviewTask"], "fix")
            self.assertEqual(data["stats"]["reviewEvents"], 1)
            self.assertEqual(data["reviewEvents"][0]["task"], "fix")
            self.assertEqual(data["reviewEvents"][0]["sources"], ["contigA", "contigA_left"])

    def test_manual_task_dashboard_embeds_modular_evidence_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            review = tmp_path / "fix_review.tsv"
            read_paf = tmp_path / "reads_to_assembly.paf"
            gaf = tmp_path / "reads_to_graph.gaf"
            read_paf.write_text("")
            gaf.write_text("")
            review.write_text(
                "schema\tevent_id\ttask\taction\ttarget\taccept\tstatus\tconfidence\treason\tnotes\t"
                "source_contig\tgraph_node\tgraph_node_status\tlongread_spanning_reads\t"
                "longread_split_reads\tlongread_left_edge_reads\tlongread_right_edge_reads\t"
                "gaf_node\tgaf_node_reads\tgaf_node_traversals\tgaf_node_orientations\t"
                "gaf_path_examples\n"
                "chromosort-review-event-v1\tfix:contigA:1\tfix\tsplit_piece\tcontigA\tno\t"
                "candidate\t.\ttest\t\tcontigA\tcontigA\tpresent\t3\t1\t2\t2\tcontigA\t4\t5\t+:4\t"
                "contigA+>bridge+\n"
            )
            output_html = run_manual_mode(
                tmp_path,
                "fix",
                "--review-table",
                str(review),
                "--read-paf",
                str(read_paf),
                "--min-read-mapq",
                "15",
                "--gaf",
                str(gaf),
                "--min-gaf-mapq",
                "30",
            )
            text = output_html.read_text()
            data = html_data(output_html)

            self.assertIn("evidencePanels", text)
            self.assertIn("Long-read PAF", text)
            self.assertIn("Long-read GAF", text)
            self.assertEqual(data["inputs"]["readPaf"], str(read_paf))
            self.assertEqual(data["inputs"]["gaf"], str(gaf))
            self.assertEqual(data["settings"]["minReadMapq"], 15)
            self.assertEqual(data["settings"]["minGafMapq"], 30)
            fields = data["reviewEvents"][0]["fields"]
            self.assertEqual(fields["graph_node"], "contigA")
            self.assertEqual(fields["longread_spanning_reads"], "3")
            self.assertEqual(fields["gaf_node_reads"], "4")

    def test_manual_apply_piece_recipe_cuts_inverts_and_removes(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            assembly = tmp_path / "assembly.fa"
            assembly.write_text(
                ">ctg1\nAACCGTTA\n"
                ">ctg2\nGGGCCC\n"
                ">drop\nTTTT\n"
            )
            recipe = {
                "schema": "chromosort-manual-v1",
                "scaffoldEnabled": False,
                "gapBp": 2,
                "pieces": [
                    {
                        "id": "p1",
                        "source": "ctg1",
                        "name": "ctg1_left",
                        "start": 1,
                        "end": 4,
                        "strand": "+",
                        "scaffold": "chr1",
                        "removed": False,
                    },
                    {
                        "id": "p2",
                        "source": "ctg1",
                        "name": "ctg1_right_rc",
                        "start": 5,
                        "end": 8,
                        "strand": "-",
                        "scaffold": "chr1",
                        "removed": False,
                    },
                    {
                        "id": "p3",
                        "source": "ctg2",
                        "name": "ctg2",
                        "start": 1,
                        "end": 6,
                        "strand": "+",
                        "scaffold": "chr2",
                        "removed": False,
                    },
                    {
                        "id": "p4",
                        "source": "drop",
                        "name": "drop",
                        "start": 1,
                        "end": 4,
                        "strand": "+",
                        "scaffold": "unplaced",
                        "removed": True,
                    },
                ],
            }

            output_fasta, report = run_manual_apply(tmp_path, assembly, recipe)

            self.assertEqual(
                read_fasta(output_fasta),
                {
                    "ctg1_left": "AACC",
                    "ctg1_right_rc": "TAAC",
                    "ctg2": "GGGCCC",
                },
            )
            rows = read_report(report)
            self.assertEqual(rows[-1]["removed"], "yes")
            self.assertEqual(rows[-1]["output_name"], ".")
            agp_rows = read_agp_rows(Path(str(output_fasta) + ".agp"))
            self.assertEqual(agp_rows[1][0], "ctg1_right_rc")
            self.assertEqual(agp_rows[1][5], "ctg1")
            self.assertEqual(agp_rows[1][6], "5")
            self.assertEqual(agp_rows[1][7], "8")
            self.assertEqual(agp_rows[1][8], "-")
            self.assertTrue(Path(str(output_fasta) + ".components.tsv").exists())
            checklist = Path(str(output_fasta) + ".submission_checklist.tsv").read_text()
            self.assertIn("fasta_agp_object_match\tok\t0", checklist)

    def test_manual_apply_can_scaffold_active_pieces(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            assembly = tmp_path / "assembly.fa"
            assembly.write_text(">ctg1\nAACCGTTA\n>ctg2\nGGGCCC\n")
            recipe = {
                "schema": "chromosort-manual-v1",
                "scaffoldEnabled": True,
                "gapBp": 3,
                "pieces": [
                    {
                        "id": "p1",
                        "source": "ctg1",
                        "name": "left",
                        "start": 1,
                        "end": 4,
                        "strand": "+",
                        "scaffold": "chr1",
                        "removed": False,
                    },
                    {
                        "id": "p2",
                        "source": "ctg1",
                        "name": "right",
                        "start": 5,
                        "end": 8,
                        "strand": "-",
                        "scaffold": "chr1",
                        "removed": False,
                    },
                    {
                        "id": "p3",
                        "source": "ctg2",
                        "name": "ctg2",
                        "start": 1,
                        "end": 6,
                        "strand": "+",
                        "scaffold": "chr2",
                        "removed": False,
                    },
                ],
            }

            output_fasta, _ = run_manual_apply(tmp_path, assembly, recipe)

            self.assertEqual(
                read_fasta(output_fasta),
                {
                    "chr1": "AACCNNNTAAC",
                    "chr2": "GGGCCC",
                },
            )
            agp_rows = read_agp_rows(Path(str(output_fasta) + ".agp"))
            gap_row = next(row for row in agp_rows if row[4] == "N")
            self.assertEqual(gap_row[0], "chr1")
            self.assertEqual(gap_row[5], "3")
            self.assertEqual(gap_row[6], "contig")
            self.assertEqual(gap_row[7], "no")
            self.assertEqual(gap_row[8], "na")


if __name__ == "__main__":
    unittest.main()
