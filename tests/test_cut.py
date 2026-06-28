import csv
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run_cut(tmp_path, assembly_fasta, *extra_args):
    output_fasta = tmp_path / "cut.fa"
    report = tmp_path / "cut.tsv"
    cmd = [
        sys.executable,
        "-m",
        "chromosort.cut",
        "--assembly-fasta",
        str(assembly_fasta),
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


def run_cut_raw(tmp_path, assembly_fasta, *extra_args):
    output_fasta = tmp_path / "cut.fa"
    report = tmp_path / "cut.tsv"
    cmd = [
        sys.executable,
        "-m",
        "chromosort.cut",
        "--assembly-fasta",
        str(assembly_fasta),
        "--output-fasta",
        str(output_fasta),
        "--report",
        str(report),
        *extra_args,
    ]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    return subprocess.run(cmd, cwd=ROOT, env=env, text=True, capture_output=True)


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


class CutTests(unittest.TestCase):
    def write_assembly(self, tmp_path):
        assembly = tmp_path / "assembly.fa"
        assembly.write_text(
            ">contig_1 description\n"
            "ACGTACGTAA\n"
            ">contig_2\n"
            "CCCCGG\n"
            ">contig_3\n"
            "TTTT\n"
        )
        return assembly

    def test_shorthand_cuts_one_contig_at_multiple_positions(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            assembly = self.write_assembly(tmp_path)
            output_fasta, report = run_cut(
                tmp_path,
                assembly,
                "--contig",
                "contig_1",
                "--pos",
                "4",
                "7",
                "--simple-headers",
            )

            self.assertEqual(
                read_fasta(output_fasta),
                {
                    "contig_1_cut001": "ACGT",
                    "contig_1_cut002": "ACG",
                    "contig_1_cut003": "TAA",
                    "contig_2": "CCCCGG",
                    "contig_3": "TTTT",
                },
            )
            rows = read_report(report)
            self.assertEqual(
                [
                    (row["new_contig"], row["slice_start"], row["slice_end"], row["piece_bp"])
                    for row in rows
                ],
                [
                    ("contig_1_cut001", "1", "4", "4"),
                    ("contig_1_cut002", "5", "7", "3"),
                    ("contig_1_cut003", "8", "10", "3"),
                ],
            )
            self.assertTrue(all(row["cut_after_positions"] == "4,7" for row in rows))
            agp_rows = read_agp_rows(Path(str(output_fasta) + ".agp"))
            self.assertEqual(agp_rows[0][0], "contig_1_cut001")
            self.assertEqual(agp_rows[0][5], "contig_1")
            self.assertEqual(agp_rows[0][6], "1")
            self.assertEqual(agp_rows[0][7], "4")
            self.assertEqual(agp_rows[-1][0], "contig_3")
            self.assertTrue(Path(str(output_fasta) + ".components.tsv").exists())
            checklist = Path(str(output_fasta) + ".submission_checklist.tsv").read_text()
            self.assertIn("fasta_agp_length_match\tok\t0", checklist)

    def test_repeated_cut_specs_cut_multiple_contigs(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            assembly = self.write_assembly(tmp_path)
            output_fasta, report = run_cut(
                tmp_path,
                assembly,
                "--cut",
                "contig_1:3,5",
                "--cut",
                "contig_2:4",
                "--name-separator",
                "_part",
                "--simple-headers",
            )

            self.assertEqual(
                read_fasta(output_fasta),
                {
                    "contig_1_part001": "ACG",
                    "contig_1_part002": "TA",
                    "contig_1_part003": "CGTAA",
                    "contig_2_part001": "CCCC",
                    "contig_2_part002": "GG",
                    "contig_3": "TTTT",
                },
            )
            self.assertEqual(
                [row["original_contig"] for row in read_report(report)],
                ["contig_1", "contig_1", "contig_1", "contig_2", "contig_2"],
            )

    def test_cuts_file_accepts_header_and_multiple_positions(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            assembly = self.write_assembly(tmp_path)
            cuts_file = tmp_path / "cuts.tsv"
            cuts_file.write_text("contig\tposition\ncontig_2\t2\t5\n")
            output_fasta, _ = run_cut(
                tmp_path,
                assembly,
                "--cuts-file",
                str(cuts_file),
                "--simple-headers",
            )

            records = read_fasta(output_fasta)
            self.assertEqual(records["contig_2_cut001"], "CC")
            self.assertEqual(records["contig_2_cut002"], "CCG")
            self.assertEqual(records["contig_2_cut003"], "G")

    def test_rejects_terminal_cut_position(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            assembly = self.write_assembly(tmp_path)
            result = run_cut_raw(tmp_path, assembly, "--cut", "contig_1:10")

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("between 1 and 9", result.stderr)

    def test_rejects_duplicate_cut_position(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            assembly = self.write_assembly(tmp_path)
            result = run_cut_raw(tmp_path, assembly, "--cut", "contig_1:4,4")

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Duplicate cut position", result.stderr)


if __name__ == "__main__":
    unittest.main()
