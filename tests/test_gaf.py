import tempfile
import unittest
from dataclasses import dataclass
from pathlib import Path

from chromosort.gaf import (
    choose_gaf_supported_path,
    contains_oriented_subpath,
    gaf_path_supports,
    parse_gaf_path,
    read_gaf,
    reverse_oriented_nodes,
    summarize_gaf_traversal,
)


@dataclass(frozen=True)
class Edge:
    source_key: tuple
    target_key: tuple


class GafTests(unittest.TestCase):
    def test_parse_gaf_path_and_reverse(self):
        nodes = parse_gaf_path(">left>bridge<right")

        self.assertEqual(nodes, [("left", "+"), ("bridge", "+"), ("right", "-")])
        self.assertEqual(
            reverse_oriented_nodes(nodes),
            [("right", "+"), ("bridge", "-"), ("left", "-")],
        )

    def test_read_gaf_filters_mapq(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "reads.gaf"
            path.write_text(
                "read_keep\t10\t0\t10\t+\t>left>right\t10\t0\t10\t10\t10\t30\n"
                "read_drop\t10\t0\t10\t+\t>left>alt\t10\t0\t10\t10\t10\t5\n"
            )

            records = read_gaf(path, min_mapq=20)

        self.assertEqual([record.query for record in records], ["read_keep"])
        self.assertEqual(records[0].path, [("left", "+"), ("right", "+")])

    def test_path_support_counts_forward_and_reverse_reads(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "reads.gaf"
            path.write_text(
                "read_forward\t10\t0\t10\t+\t>left>bridge>right\t10\t0\t10\t10\t10\t60\n"
                "read_reverse\t10\t0\t10\t+\t<right<bridge<left\t10\t0\t10\t10\t10\t60\n"
                "read_alt\t10\t0\t10\t+\t>left>alt>right\t10\t0\t10\t10\t10\t60\n"
            )
            records = read_gaf(path)
            paths = [
                [
                    Edge(("left", "+"), ("bridge", "+")),
                    Edge(("bridge", "+"), ("right", "+")),
                ],
                [
                    Edge(("left", "+"), ("alt", "+")),
                    Edge(("alt", "+"), ("right", "+")),
                ],
            ]

            supports = gaf_path_supports(records, paths)

        self.assertTrue(
            contains_oriented_subpath(
                [("left", "+"), ("bridge", "+"), ("right", "+")],
                [("bridge", "+"), ("right", "+")],
            )
        )
        self.assertEqual(supports, [2, 1])
        self.assertEqual(choose_gaf_supported_path(paths, supports, min_support=2), 0)

    def test_summarize_gaf_traversal_reports_alternate_support(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "reads.gaf"
            path.write_text(
                "read_selected\t10\t0\t10\t+\t>left>bridge>right\t10\t0\t10\t10\t10\t60\n"
                "read_alt_1\t10\t0\t10\t+\t>left>alt>right\t10\t0\t10\t10\t10\t60\n"
                "read_alt_2\t10\t0\t10\t+\t>left>alt>right\t10\t0\t10\t10\t10\t60\n"
            )
            records = read_gaf(path)
            paths = [
                [
                    Edge(("left", "+"), ("bridge", "+")),
                    Edge(("bridge", "+"), ("right", "+")),
                ],
                [
                    Edge(("left", "+"), ("alt", "+")),
                    Edge(("alt", "+"), ("right", "+")),
                ],
            ]

            summary = summarize_gaf_traversal(records, paths, selected_index=0, min_support=1)

        self.assertEqual(summary.selected_path_nodes, "left+,bridge+,right+")
        self.assertEqual(summary.selected_support, 1)
        self.assertEqual(summary.best_alt_path_nodes, "left+,alt+,right+")
        self.assertEqual(summary.best_alt_support, 2)
        self.assertEqual(summary.best_path_nodes, "left+,alt+,right+")
        self.assertEqual(summary.support_status, "supports_alternate")
        self.assertEqual(summary.selected_reads, ("read_selected",))
        self.assertEqual([item.support for item in summary.path_supports], [1, 2])


if __name__ == "__main__":
    unittest.main()
