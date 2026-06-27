import os
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


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


def run_cli_raw(*args):
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    return subprocess.run(
        [sys.executable, "-m", "chromosort.cli", *args],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
    )


class CliTests(unittest.TestCase):
    def test_top_level_help_lists_subcommands(self):
        result = run_cli("--help")
        self.assertIn("sort", result.stdout)
        self.assertIn("clean", result.stdout)
        self.assertIn("eval", result.stdout)
        self.assertIn("fix", result.stdout)
        self.assertIn("cut", result.stdout)
        self.assertIn("manual", result.stdout)
        self.assertIn("gafprep", result.stdout)
        self.assertIn("scaffold", result.stdout)
        self.assertIn("gapfill", result.stdout)
        self.assertIn("plot", result.stdout)
        self.assertIn("graph-map", result.stdout)

    def test_subcommand_help_dispatches_to_full_commands(self):
        sort_help = run_cli("sort", "--help").stdout
        clean_help = run_cli("clean", "--help").stdout
        eval_help = run_cli("eval", "--help").stdout
        eval_fix_help = run_cli("eval", "fix", "--help").stdout
        eval_scaffold_help = run_cli("eval", "scaffold", "--help").stdout
        eval_gapfill_help = run_cli("eval", "gapfill", "--help").stdout
        fix_help = run_cli("fix", "--help").stdout
        cut_help = run_cli("cut", "--help").stdout
        manual_help = run_cli("manual", "--help").stdout
        manual_apply_help = run_cli("manual", "apply", "--help").stdout
        manual_fix_help = run_cli("manual", "fix", "--help").stdout
        manual_scaffold_help = run_cli("manual", "scaffold", "--help").stdout
        manual_gapfill_help = run_cli("manual", "gapfill", "--help").stdout
        gafprep_help = run_cli("gafprep", "--help").stdout
        scaffold_help = run_cli("scaffold", "--help").stdout
        gapfill_help = run_cli("gapfill", "--help").stdout
        plot_help = run_cli("plot", "--help").stdout
        graph_map_help = run_cli("graph-map", "--help").stdout
        self.assertIn("--output-prefix", sort_help)
        self.assertIn("--paf", sort_help)
        self.assertIn("--gfa", sort_help)
        self.assertIn("--graph-guard", sort_help)
        self.assertIn("--fix-scope", clean_help)
        self.assertIn("--fix-mode", clean_help)
        self.assertIn("--discarded-fasta", clean_help)
        self.assertIn("fix", eval_help)
        self.assertIn("scaffold", eval_help)
        self.assertIn("gapfill", eval_help)
        self.assertIn("--read-paf", eval_fix_help)
        self.assertIn("--gaf", eval_fix_help)
        self.assertIn("--output-prefix", eval_fix_help)
        self.assertIn("--ordered-fasta", eval_scaffold_help)
        self.assertIn("--assignments", eval_scaffold_help)
        self.assertIn("--gaf", eval_scaffold_help)
        self.assertIn("--gfa", eval_gapfill_help)
        self.assertIn("--read-paf", eval_gapfill_help)
        self.assertIn("--all", fix_help)
        self.assertIn("--mode", fix_help)
        self.assertIn("--reviewed-plan", fix_help)
        self.assertNotIn("--auto", fix_help)
        self.assertIn("--paf", fix_help)
        self.assertIn("--gfa", fix_help)
        self.assertIn("--graph-report", fix_help)
        self.assertIn("--graph-guard", fix_help)
        self.assertIn("--cut", cut_help)
        self.assertIn("--cuts-file", cut_help)
        self.assertIn("--output-fasta", cut_help)
        self.assertIn("--output-html", manual_help)
        self.assertIn("--embed-sequences", manual_help)
        self.assertIn("--gfa", manual_help)
        self.assertIn("--read-paf", manual_help)
        self.assertIn("--gaf", manual_help)
        self.assertIn("--recipe", manual_apply_help)
        self.assertIn("--review-table", manual_fix_help)
        self.assertIn("--read-paf", manual_fix_help)
        self.assertIn("--gaf", manual_fix_help)
        self.assertIn("--review-table", manual_scaffold_help)
        self.assertIn("--read-paf", manual_scaffold_help)
        self.assertIn("--gaf", manual_scaffold_help)
        self.assertIn("--review-table", manual_gapfill_help)
        self.assertIn("--read-paf", manual_gapfill_help)
        self.assertIn("--gaf", manual_gapfill_help)
        self.assertIn("--assembly-gfa", gafprep_help)
        self.assertIn("--eval-review-table", gafprep_help)
        self.assertIn("--graphaligner-threads", gafprep_help)
        self.assertIn("--fixed-gap-bp", scaffold_help)
        self.assertIn("--reviewed-plan", scaffold_help)
        self.assertIn("--gfa", scaffold_help)
        self.assertIn("--graph-overlap-policy", scaffold_help)
        self.assertIn("--graph-max-path-edges", scaffold_help)
        self.assertIn("--apply", gapfill_help)
        self.assertIn("--gaf", gapfill_help)
        self.assertIn("--hic-pairs", gapfill_help)
        self.assertIn("--ref-paf", gapfill_help)
        self.assertIn("--reviewed-plan", gapfill_help)
        self.assertIn("--review-html", gapfill_help)
        self.assertIn("--max-path-edges", gapfill_help)
        self.assertIn("--min-gaf-path-support", gapfill_help)
        self.assertIn("--min-hic-path-support", gapfill_help)
        self.assertIn("--min-ref-path-support", gapfill_help)
        self.assertIn("--include-fill-sequences", gapfill_help)
        self.assertIn("--per-ref", plot_help)
        self.assertIn("--sel-ref", plot_help)
        self.assertIn("--paf", plot_help)
        self.assertIn("--formats", plot_help)
        self.assertIn("--gfa-overlay", plot_help)
        self.assertIn("--ctg-fasta", graph_map_help)
        self.assertIn("--utg-gfa", graph_map_help)
        self.assertIn("--output-prefix", graph_map_help)

    def test_alignment_inputs_are_mutually_exclusive(self):
        result = run_cli_raw(
            "sort",
            "--ref-fasta",
            "ref.fa",
            "--assembly-fasta",
            "assembly.fa",
            "--coords",
            "sample.coords",
            "--paf",
            "sample.paf",
            "--output-prefix",
            "sample",
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("not allowed with argument", result.stderr)

    def test_fix_scope_flags_are_mutually_exclusive(self):
        result = run_cli_raw(
            "fix",
            "--assembly-fasta",
            "assembly.fa",
            "--coords",
            "sample.coords",
            "--contigs",
            "contig_01",
            "--all",
            "--output-fasta",
            "fixed.fa",
            "--report",
            "fixed.tsv",
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("use either --all or --contigs/--contigs-file", result.stderr)

    def test_gapfill_replaces_old_fill_subcommand(self):
        result = run_cli_raw("fill", "--help")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("unknown command: fill", result.stderr)

    def test_graph_safety_flags_require_graph_input(self):
        sort_result = run_cli_raw(
            "sort",
            "--ref-fasta",
            "ref.fa",
            "--assembly-fasta",
            "assembly.fa",
            "--paf",
            "sample.paf",
            "--output-prefix",
            "sample",
            "--graph-guard",
        )
        self.assertNotEqual(sort_result.returncode, 0)
        self.assertIn("--graph-guard requires --gfa", sort_result.stderr)

        fix_result = run_cli_raw(
            "fix",
            "--assembly-fasta",
            "assembly.fa",
            "--paf",
            "sample.paf",
            "--contigs",
            "contig_01",
            "--output-fasta",
            "fixed.fa",
            "--report",
            "fixed.tsv",
            "--graph-guard",
        )
        self.assertNotEqual(fix_result.returncode, 0)
        self.assertIn("--graph-guard requires --gfa", fix_result.stderr)

        scaffold_result = run_cli_raw(
            "scaffold",
            "--ordered-fasta",
            "ordered.fa",
            "--assignments",
            "assignments.tsv",
            "--output-prefix",
            "sample",
            "--graph-overlap-policy",
            "warn",
        )
        self.assertNotEqual(scaffold_result.returncode, 0)
        self.assertIn("--graph-overlap-policy warn/confirm requires --gfa", scaffold_result.stderr)


if __name__ == "__main__":
    unittest.main()
