import csv
import gzip
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from chromosort.gafprep import (
    TARGET_COLUMNS,
    build_targets,
    extract_fastq,
    output_paths,
    parse_read_paf,
    sanitize_gfa,
    select_reads,
    write_graphaligner_script,
)


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


def default_args(**overrides):
    values = {
        "target_padding": 10,
        "contig_end_window": 30,
        "target_reads_per_interval": 2,
        "min_mapq": 20,
        "min_aligned_bp": 10,
        "include_secondary_paf": False,
        "background_bin_size": 100,
        "background_reads_per_bin": 1,
        "max_reads": None,
        "max_reads_per_contig": None,
        "seed": 1,
        "graphaligner_bin": "GraphAligner",
        "graphaligner_threads": 16,
        "graphaligner_preset": "vg",
        "precise_clipping": "0.9",
    }
    values.update(overrides)
    return SimpleNamespace(**values)


class GafPrepTests(unittest.TestCase):
    def write_path(self, tmp_path, name, text):
        path = tmp_path / name
        path.write_text(text)
        return path

    def test_parse_read_paf_filters_secondary_mapq_and_aligned_bp(self):
        with tempfile.TemporaryDirectory() as tmp:
            paf = self.write_path(
                Path(tmp),
                "reads.paf",
                "\n".join(
                    [
                        "read1\t100\t0\t80\t+\tctgA\t200\t40\t120\t78\t80\t60\ttp:A:P",
                        "lowmapq\t100\t0\t80\t+\tctgA\t200\t40\t120\t78\t80\t5\ttp:A:P",
                        "short\t100\t0\t8\t+\tctgA\t200\t40\t48\t8\t8\t60\ttp:A:P",
                        "secondary\t100\t0\t80\t+\tctgA\t200\t40\t120\t78\t80\t60\ttp:A:S",
                    ]
                )
                + "\n",
            )

            rows = list(parse_read_paf(paf, min_mapq=20, min_aligned_bp=10, include_secondary=False))

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].read_id, "read1")
        self.assertEqual(rows[0].contig, "ctgA")
        self.assertEqual(rows[0].aln_start, 41)
        self.assertEqual(rows[0].matches, 78)
        self.assertAlmostEqual(rows[0].identity_estimate, 97.5)

    def test_build_targets_for_fix_scaffold_gapfill_and_fallback(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            fix = self.write_path(
                tmp_path,
                "sample.fix_review.tsv",
                "event_id\ttask\taction\ttarget\treason\tsource_contig\tlongread_breakpoint_position\n"
                "fix:ctgA:1\tfix\tsplit_piece\tctgA\tcandidate_breakpoint\tctgA\t100\n"
                "fix:ctgC:no_split\tfix\tno_split\tctgC\tneeds_context\tctgC\t.\n",
            )
            scaffold = self.write_path(
                tmp_path,
                "sample.scaffold_review.tsv",
                "event_id\ttask\taction\tleft_contig\tright_contig\tgraph_status\n"
                "scaffold:1\tscaffold\tscaffold_gap\tctgA\tctgB\tshort_path\n",
            )
            gapfill = self.write_path(
                tmp_path,
                "sample.gapfill_review.tsv",
                "event_id\ttask\taction\tleft_contig\tright_contig\tfill_status\tpath_nodes\n"
                "gapfill:1\tgapfill\tfill_path\tctgB\tctgC\tambiguous\tctgB+,alt+,ctgC+\n",
            )
            lengths = {"ctgA": 200, "ctgB": 120, "ctgC": 80}

            targets = build_targets([fix, scaffold, gapfill], "auto", lengths, default_args())

        rows = [target.__dict__ for target in targets]
        types = [row["target_type"] for row in rows]
        self.assertIn("fix_breakpoint", types)
        self.assertIn("fix_fallback_left_end", types)
        self.assertIn("fix_fallback_right_end", types)
        self.assertIn("scaffold_left_flank", types)
        self.assertIn("scaffold_right_flank", types)
        self.assertIn("gapfill_left_flank", types)
        self.assertIn("gapfill_right_flank", types)
        fix_target = next(row for row in rows if row["target_type"] == "fix_breakpoint")
        self.assertEqual((fix_target["start"], fix_target["end"]), (90, 110))
        left_flank = next(row for row in rows if row["target_type"] == "scaffold_left_flank")
        self.assertEqual((left_flank["contig"], left_flank["start"], left_flank["end"]), ("ctgA", 171, 200))

    def test_select_reads_preserves_links_deduplicates_and_background_is_seeded(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            review = self.write_path(
                tmp_path,
                "sample.scaffold_review.tsv",
                "event_id\ttask\taction\tleft_contig\tright_contig\n"
                "scaffold:1\tscaffold\tscaffold_gap\tctgA\tctgB\n",
            )
            paf = self.write_path(
                tmp_path,
                "reads.paf",
                "\n".join(
                    [
                        "bridge\t200\t0\t70\t+\tctgA\t200\t170\t200\t30\t30\t60\ttp:A:P",
                        "bridge\t200\t80\t140\t+\tctgB\t120\t0\t35\t35\t35\t60\ttp:A:P",
                        "nearA\t100\t0\t40\t+\tctgA\t200\t160\t198\t38\t38\t50\ttp:A:P",
                        "background\t100\t0\t40\t+\tctgC\t500\t250\t290\t40\t40\t60\ttp:A:P",
                    ]
                )
                + "\n",
            )
            targets = build_targets([review], "auto", {"ctgA": 200, "ctgB": 120, "ctgC": 500}, default_args())
            args = default_args(target_reads_per_interval=2, background_reads_per_bin=1, seed=9)

            first = select_reads(paf, targets, args)
            second = select_reads(paf, targets, args)

        self.assertEqual(sorted(first.selected), sorted(second.selected))
        self.assertIn("bridge", first.selected)
        bridge = first.selected["bridge"]
        self.assertIn("target_interval", bridge.selection_reasons)
        self.assertIn("bridge_two_targets", bridge.selection_reasons)
        self.assertEqual(len(bridge.links), 2)
        self.assertEqual(sorted(bridge.target_ids), [targets[0].target_id, targets[1].target_id])
        self.assertTrue(any(row.is_background for row in first.selected.values()))

    def test_extract_fastq_streams_plain_and_gzipped_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            plain = self.write_path(
                tmp_path,
                "reads.fastq",
                "@read1 comment\nAAAA\n+\n!!!!\n@read2\nCCCC\n+\n####\n",
            )
            gz_in = tmp_path / "reads.fastq.gz"
            with gzip.open(gz_in, "wt") as out:
                out.write("@read1 comment\nAAAA\n+\n!!!!\n@read2\nCCCC\n+\n####\n")
            plain_out = tmp_path / "plain.selected.fastq.gz"
            gz_out = tmp_path / "gz.selected.fastq.gz"

            plain_count, plain_missing = extract_fastq(plain, plain_out, {"read2"})
            gz_count, gz_missing = extract_fastq(gz_in, gz_out, {"read1", "absent"})

            with gzip.open(plain_out, "rt") as fh:
                plain_text = fh.read()
            with gzip.open(gz_out, "rt") as fh:
                gz_text = fh.read()

        self.assertEqual((plain_count, plain_missing), (1, []))
        self.assertIn("@read2\nCCCC", plain_text)
        self.assertEqual(gz_count, 1)
        self.assertEqual(gz_missing, ["absent"])
        self.assertIn("@read1 comment\nAAAA", gz_text)

    def test_sanitize_gfa_removes_a_records_and_only_full_consuming_links(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            gfa = self.write_path(
                tmp_path,
                "graph.gfa",
                "H\tVN:Z:1.0\n"
                "S\tctgA\tAAAAA\tLN:i:5\n"
                "S\tctgB\tCCCCCC\tLN:i:6\n"
                "S\tother\tTTTTT\tLN:i:5\n"
                "A\tctgA\t0\tread1\t+\t0\t5\t5M\n"
                "L\tctgA\t+\tctgB\t+\t3M\n"
                "L\tctgA\t+\tother\t+\t5M\n"
                "L\tctgB\t+\tother\t+\t6M\n"
                "P\tpath1\tctgA+,ctgB+\t3M\n"
                "W\tsample\t0\tpath2\t0\t10\t>ctgA>ctgB\n",
            )
            sanitized = tmp_path / "sanitized.gfa"
            dropped = tmp_path / "dropped.tsv"

            summary = sanitize_gfa(gfa, sanitized, dropped, {"ctgA"})
            out_text = sanitized.read_text()
            dropped_rows = read_tsv(dropped)

        self.assertNotIn("\nA\t", out_text)
        self.assertIn("L\tctgA\t+\tctgB\t+\t3M", out_text)
        self.assertNotIn("L\tctgA\t+\tother\t+\t5M", out_text)
        self.assertNotIn("L\tctgB\t+\tother\t+\t6M", out_text)
        self.assertIn("P\tpath1", out_text)
        self.assertIn("W\tsample", out_text)
        self.assertEqual(summary["removed_A_records"], 1)
        self.assertEqual(summary["dropped_links"], 2)
        self.assertEqual([row["near_selected_target"] for row in dropped_rows], ["yes", "no"])

    def test_graphaligner_script_is_executable(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            script = tmp_path / "run_graphaligner.sh"
            args = default_args(graphaligner_bin="GraphAligner", graphaligner_threads=4)

            write_graphaligner_script(
                script,
                args,
                tmp_path / "sample.graphaligner.gfa",
                tmp_path / "sample.selected.fastq.gz",
                tmp_path / "sample.gaf",
            )

            self.assertTrue(os.access(script, os.X_OK))
            text = script.read_text()
            self.assertIn("GraphAligner", text)
            self.assertIn("--precise-clipping 0.9", text)
            self.assertIn("-t 4", text)

    def test_cli_smoke_writes_expected_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            assembly = self.write_path(tmp_path, "assembly.fa", ">ctgA\n" + "A" * 120 + "\n")
            gfa = self.write_path(
                tmp_path,
                "graph.gfa",
                "H\tVN:Z:1.0\n"
                "S\tctgA\t" + "A" * 120 + "\tLN:i:120\n"
                "S\tother\tCCCC\tLN:i:4\n"
                "A\tctgA\t0\tread1\t+\t0\t20\t20M\n"
                "L\tctgA\t+\tother\t+\t120M\n",
            )
            paf = self.write_path(
                tmp_path,
                "reads.paf",
                "read1\t80\t0\t30\t+\tctgA\t120\t45\t75\t30\t30\t60\ttp:A:P\n",
            )
            fastq = self.write_path(tmp_path, "reads.fastq", "@read1\nACGT\n+\n!!!!\n")
            review = self.write_path(
                tmp_path,
                "sample.fix_review.tsv",
                "event_id\ttask\taction\ttarget\tsource_contig\tlongread_breakpoint_position\n"
                "fix:ctgA:1\tfix\tsplit_piece\tctgA\tctgA\t60\n",
            )
            prefix = tmp_path / "out" / "sample.gafprep"

            run_cli(
                "gafprep",
                "--assembly-fasta",
                str(assembly),
                "--assembly-gfa",
                str(gfa),
                "--read-paf",
                str(paf),
                "--reads",
                str(fastq),
                "--eval-review-table",
                str(review),
                "--output-prefix",
                str(prefix),
                "--target-padding",
                "10",
                "--min-aligned-bp",
                "1",
                "--background-reads-per-bin",
                "0",
            )
            paths = output_paths(prefix)

            selected = read_tsv(paths["selected_reads"])
            links = read_tsv(paths["selected_read_review_links"])
            with gzip.open(paths["selected_fastq"], "rt") as fh:
                selected_fastq = fh.read()

            self.assertEqual(selected[0]["read_id"], "read1")
            self.assertEqual(links[0]["review_row_id"], "fix:ctgA:1")
            self.assertIn("@read1", selected_fastq)
            self.assertTrue(paths["graphaligner_script"].exists())
            self.assertTrue(os.access(paths["graphaligner_script"], os.X_OK))
            self.assertTrue(paths["graphaligner_gfa"].exists())
            self.assertEqual(len(read_tsv(paths["dropped_gfa_links"])), 1)


if __name__ == "__main__":
    unittest.main()
