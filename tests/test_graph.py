import tempfile
import unittest
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from chromosort.graph import parse_gfa_overlap_bp, read_gfa


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


if __name__ == "__main__":
    unittest.main()
