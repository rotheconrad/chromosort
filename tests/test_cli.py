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
        self.assertIn("fix", result.stdout)
        self.assertIn("scaffold", result.stdout)
        self.assertIn("plot", result.stdout)

    def test_subcommand_help_dispatches_to_full_commands(self):
        sort_help = run_cli("sort", "--help").stdout
        fix_help = run_cli("fix", "--help").stdout
        scaffold_help = run_cli("scaffold", "--help").stdout
        plot_help = run_cli("plot", "--help").stdout
        self.assertIn("--output-prefix", sort_help)
        self.assertIn("--paf", sort_help)
        self.assertIn("--auto", fix_help)
        self.assertIn("--paf", fix_help)
        self.assertIn("--fixed-gap-bp", scaffold_help)
        self.assertIn("--per-ref", plot_help)
        self.assertIn("--paf", plot_help)
        self.assertIn("--formats", plot_help)

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


if __name__ == "__main__":
    unittest.main()
