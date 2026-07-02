import csv
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data"
CHIMERIC_DATA = DATA / "chimeric"
GRAPH_DATA = DATA / "graph_gotchas"
MIXED_EVIDENCE_DATA = DATA / "mixed_evidence"
SCAFFOLD_DATA = DATA / "scaffold"


def run_chromo(*args):
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    subprocess.run(
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


class AgentSmokeWorkflowTests(unittest.TestCase):
    def test_agent_smoke_workflow_writes_expected_artifacts(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)

            sort_prefix = tmp_path / "sort"
            run_chromo(
                "sort",
                "--ref-fasta",
                str(DATA / "ref.fa"),
                "--assembly-fasta",
                str(DATA / "assembly.fa"),
                "--paf",
                str(DATA / "sample.paf"),
                "--output-prefix",
                str(sort_prefix),
                "--min-aligned-bp",
                "10",
                "--min-query-cov",
                "0.10",
                "--min-novel-ref-bp",
                "1",
                "--orient-to-reference",
            )
            self.assertTrue(Path(str(sort_prefix) + ".ordered.fa").exists())
            self.assertTrue(Path(str(sort_prefix) + ".ordered.agp").exists())
            self.assertGreater(
                len(read_tsv(Path(str(sort_prefix) + ".contig_assignments.tsv"))),
                0,
            )

            plot_prefix = tmp_path / "plot"
            run_chromo(
                "plot",
                "--ref-fasta",
                str(DATA / "ref.fa"),
                "--assembly-fasta",
                str(DATA / "assembly.fa"),
                "--paf",
                str(DATA / "sample.paf"),
                "--assignments",
                str(Path(str(sort_prefix) + ".contig_assignments.tsv")),
                "--output-prefix",
                str(plot_prefix),
                "--formats",
                "svg",
            )
            self.assertTrue(Path(str(plot_prefix) + ".svg").exists())

            fixed_fasta = tmp_path / "fixed.fa"
            fix_report = tmp_path / "fix_report.tsv"
            run_chromo(
                "fix",
                "--assembly-fasta",
                str(CHIMERIC_DATA / "assembly.fa"),
                "--coords",
                str(CHIMERIC_DATA / "sample.coords"),
                "--contigs",
                "contig_04",
                "--mode",
                "sensitive",
                "--min-segment-bp",
                "5",
                "--min-piece-bp",
                "1",
                "--output-fasta",
                str(fixed_fasta),
                "--report",
                str(fix_report),
            )
            self.assertTrue(fixed_fasta.exists())
            self.assertTrue(fix_report.exists())

            eval_prefix = tmp_path / "eval_all"
            run_chromo(
                "eval",
                "all",
                "--assembly-fasta",
                str(CHIMERIC_DATA / "assembly.fa"),
                "--coords",
                str(CHIMERIC_DATA / "sample.coords"),
                "--contigs",
                "contig_04",
                "--ordered-fasta",
                str(SCAFFOLD_DATA / "ordered.fa"),
                "--assignments",
                str(SCAFFOLD_DATA / "assignments.tsv"),
                "--gfa",
                str(MIXED_EVIDENCE_DATA / "branch.gfa"),
                "--output-prefix",
                str(eval_prefix),
                "--min-segment-bp",
                "5",
                "--mode",
                "sensitive",
            )
            manifest = read_tsv(Path(str(eval_prefix) + ".eval_all_outputs.tsv"))
            self.assertEqual(
                [row["label"] for row in manifest],
                ["fix_review", "scaffold_review", "gapfill_review"],
            )

            scaffold_prefix = tmp_path / "scaffold"
            run_chromo(
                "scaffold",
                "--ordered-fasta",
                str(SCAFFOLD_DATA / "ordered.fa"),
                "--assignments",
                str(SCAFFOLD_DATA / "assignments.tsv"),
                "--output-prefix",
                str(scaffold_prefix),
                "--simple-headers",
            )
            self.assertTrue(Path(str(scaffold_prefix) + ".scaffold.fa").exists())
            self.assertTrue(Path(str(scaffold_prefix) + ".scaffold.agp").exists())

            gapfill_prefix = tmp_path / "gapfill"
            run_chromo(
                "gapfill",
                "--ordered-fasta",
                str(GRAPH_DATA / "gapfill_ordered.fa"),
                "--assignments",
                str(GRAPH_DATA / "gapfill_assignments.tsv"),
                "--gfa",
                str(GRAPH_DATA / "unitigs.gfa"),
                "--ref-paf",
                str(GRAPH_DATA / "unitig_to_ref.paf"),
                "--min-ref-path-support",
                "7",
                "--output-prefix",
                str(gapfill_prefix),
                "--apply",
                "--apply-all-fillable",
                "--simple-headers",
            )
            self.assertTrue(Path(str(gapfill_prefix) + ".gapfill_plan.tsv").exists())
            self.assertTrue(Path(str(gapfill_prefix) + ".gapfilled.fa").exists())
