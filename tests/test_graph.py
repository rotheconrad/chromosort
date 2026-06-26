import csv
import os
import subprocess
import tempfile
import unittest
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from chromosort.graph import build_path_projection, parse_gfa_overlap_bp, read_gfa, summarize_projections


DATA = Path(__file__).resolve().parent / "data" / "graph_gotchas"


class GraphParserTests(unittest.TestCase):
    def test_reads_nodes_edges_and_orientations_from_fixture(self):
        graph = read_gfa(DATA / "unitigs.gfa")

        self.assertEqual(len(graph.nodes), 9)
        self.assertEqual(len(graph.edges), 10)
        self.assertEqual(graph.nodes["left"].length, 8)
        self.assertEqual(graph.nodes["unsequenced_gap_probe"].length, 12)
        self.assertFalse(graph.nodes["unsequenced_gap_probe"].has_sequence)

        confident_edges = graph.direct_edges(
            "left",
            "bridge_good",
            source_orientation="+",
            target_orientation="+",
        )
        self.assertEqual(len(confident_edges), 1)
        self.assertEqual(confident_edges[0].overlap, "4M")
        self.assertEqual(confident_edges[0].overlap_bp, 4)

        self.assertFalse(
            graph.has_direct_edge(
                "left",
                "bridge_good",
                source_orientation="-",
                target_orientation="+",
            )
        )
        self.assertTrue(
            graph.has_direct_edge(
                "right",
                "reverse_only",
                source_orientation="-",
                target_orientation="-",
            )
        )

    def test_overlap_parser_keeps_complex_cigars_from_becoming_trim_lengths(self):
        self.assertEqual(parse_gfa_overlap_bp("0M"), 0)
        self.assertEqual(parse_gfa_overlap_bp("4M"), 4)
        self.assertEqual(parse_gfa_overlap_bp("2=1X"), 3)
        self.assertIsNone(parse_gfa_overlap_bp("*"))
        self.assertIsNone(parse_gfa_overlap_bp("5M1I4M"))
        with self.assertRaises(ValueError):
            parse_gfa_overlap_bp("4")

    def test_strict_parser_rejects_bad_links_and_length_mismatches(self):
        with tempfile.TemporaryDirectory() as tmp:
            bad_link = Path(tmp) / "bad_link.gfa"
            bad_link.write_text("S\tutg1\tAAAA\nL\tutg1\t+\tmissing\t+\t0M\n")
            with self.assertRaisesRegex(ValueError, "missing"):
                read_gfa(bad_link)

            bad_length = Path(tmp) / "bad_length.gfa"
            bad_length.write_text("S\tutg1\tAAAA\tLN:i:5\n")
            with self.assertRaisesRegex(ValueError, "LN:i:5"):
                read_gfa(bad_length)

    def test_noseq_p_path_projection_uses_ln_lengths(self):
        with tempfile.TemporaryDirectory() as tmp:
            gfa = Path(tmp) / "noseq.gfa"
            gfa.write_text(
                "S\tutg1\t*\tLN:i:1000\n"
                "S\tutg2\t*\tLN:i:500\n"
                "L\tutg1\t+\tutg2\t+\t0M\n"
                "P\tptg1\tutg1+,utg2+\t0M\n"
            )
            graph = read_gfa(gfa)
            warnings = []
            projections = build_path_projection(graph, path_names=["ptg1"], warnings=warnings)

        self.assertFalse(warnings)
        self.assertEqual(len(graph.paths), 1)
        self.assertFalse(graph.nodes["utg1"].has_sequence)
        self.assertEqual(graph.nodes["utg1"].length, 1000)
        self.assertEqual(
            [(row.segment, row.contig_start_0, row.contig_end_0) for row in projections],
            [("utg1", 0, 1000), ("utg2", 1000, 1500)],
        )

    def test_full_sequence_p_path_projection_uses_sequence_lengths(self):
        with tempfile.TemporaryDirectory() as tmp:
            gfa = Path(tmp) / "full.gfa"
            gfa.write_text(
                "S\tutgA\tAAAAA\tLN:i:5\n"
                "S\tutgB\tCCCC\tLN:i:4\n"
                "P\tptgA\tutgA+,utgB+\t0M\n"
            )
            graph = read_gfa(gfa)
            projections = build_path_projection(graph, path_names=["ptgA"])

        self.assertTrue(graph.nodes["utgA"].has_sequence)
        self.assertEqual([row.segment_length for row in projections], [5, 4])
        self.assertEqual(projections[-1].contig_end_0, 9)

    def test_w_walk_projection_and_reverse_step(self):
        with tempfile.TemporaryDirectory() as tmp:
            gfa = Path(tmp) / "walk.gfa"
            gfa.write_text(
                "S\tutg1\t*\tLN:i:10\n"
                "S\tutg2\t*\tLN:i:8\n"
                "W\tsample\t0\tptgW\t0\t18\t>utg1<utg2\n"
            )
            graph = read_gfa(gfa)
            projections = build_path_projection(graph, path_names=["ptgW"])

        self.assertEqual(graph.paths[0].record_type, "W")
        self.assertEqual(
            [(row.segment, row.segment_orientation, row.contig_start_0, row.contig_end_0) for row in projections],
            [("utg1", "+", 0, 10), ("utg2", "-", 10, 18)],
        )

    def test_reused_unitig_and_length_warnings_are_reported(self):
        with tempfile.TemporaryDirectory() as tmp:
            gfa = Path(tmp) / "reused.gfa"
            gfa.write_text(
                "S\tutgR\t*\tLN:i:7\n"
                "S\tutgX\t*\n"
                "P\tptg1\tutgR+,utgX+\t0M\n"
                "P\tptg2\tutgR-\t*\n"
            )
            graph = read_gfa(gfa)
            warnings = []
            projections = build_path_projection(graph, warnings=warnings)
            summaries = summarize_projections(
                projections,
                fasta_lengths={"ptg1": 20, "ptg2": 7},
                tolerance_bp=0,
                tolerance_frac=0.0,
                warnings=warnings,
            )

        reused = [row for row in projections if row.segment == "utgR"]
        self.assertEqual([row.duplicated_segment_count for row in reused], [2, 2])
        self.assertTrue(all(row.is_reused_segment for row in reused))
        self.assertTrue(any("no sequence and no LN:i" in warning.message for warning in warnings))
        self.assertTrue(any(warning.code == "path_length_mismatch" for warning in warnings))
        self.assertEqual({row.path_name: row.length_status for row in summaries}["ptg1"], "length_mismatch")

    def test_graph_map_command_writes_projection_reports(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            fasta = tmp_path / "ctg.fa"
            gfa = tmp_path / "utg.gfa"
            prefix = tmp_path / "sample.graphmap"
            fasta.write_text(">ptg1\n" + "A" * 1500 + "\n")
            gfa.write_text(
                "S\tutg1\t*\tLN:i:1000\n"
                "S\tutg2\t*\tLN:i:500\n"
                "P\tptg1\tutg1+,utg2+\t0M\n"
            )
            env = os.environ.copy()
            env["PYTHONPATH"] = str(ROOT / "src")
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "chromosort.cli",
                    "graph-map",
                    "--ctg-fasta",
                    str(fasta),
                    "--utg-gfa",
                    str(gfa),
                    "--output-prefix",
                    str(prefix),
                ],
                cwd=ROOT,
                env=env,
                check=True,
            )
            with open(Path(f"{prefix}.utg_to_ctg.tsv"), newline="") as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            with open(Path(f"{prefix}.path_summary.tsv"), newline="") as fh:
                summaries = list(csv.DictReader(fh, delimiter="\t"))

        self.assertEqual([row["unitig"] for row in rows], ["utg1", "utg2"])
        self.assertEqual(rows[1]["contig_start_0"], "1000")
        self.assertEqual(summaries[0]["length_status"], "ok")


if __name__ == "__main__":
    unittest.main()
