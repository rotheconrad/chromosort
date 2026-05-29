import tempfile
import unittest
from pathlib import Path

from chromosort.longreads import (
    read_depth_at,
    read_depth_window,
    read_long_read_paf,
    summarize_breakpoint,
    summarize_contig_bridge,
)


class LongReadEvidenceTests(unittest.TestCase):
    def write_paf(self, rows):
        tmp = tempfile.TemporaryDirectory()
        path = Path(tmp.name) / "reads_to_assembly.paf"
        with open(path, "w") as out:
            out.write("\n".join(rows) + "\n")
        self.addCleanup(tmp.cleanup)
        return path

    def test_reads_long_read_paf_and_filters_rows(self):
        path = self.write_paf(
            [
                "read1\t10000\t0\t5000\t+\tctg1\t9000\t0\t5000\t5000\t5000\t60\ttp:A:P",
                "read2\t10000\t0\t5000\t+\tctg1\t9000\t0\t5000\t4500\t5000\t10\ttp:A:P",
                "read3\t10000\t0\t5000\t+\tctg1\t9000\t0\t5000\t5000\t5000\t60\ttp:A:S",
            ]
        )

        evidence = read_long_read_paf(path, min_mapq=20)

        self.assertEqual(len(evidence.alignments), 1)
        self.assertEqual(evidence.alignments[0].read, "read1")
        self.assertEqual(evidence.alignments[0].contig, "ctg1")
        self.assertEqual(evidence.alignments[0].read_start, 1)
        self.assertEqual(evidence.alignments[0].contig_end, 5000)

    def test_summarizes_breakpoint_spanning_and_split_reads(self):
        path = self.write_paf(
            [
                "span\t12000\t1000\t9000\t+\tctg1\t10000\t999\t8000\t7000\t7001\t60\ttp:A:P",
                "split\t12000\t1000\t4500\t+\tctg1\t10000\t1499\t4900\t3400\t3401\t60\ttp:A:P",
                "split\t12000\t4600\t8500\t+\tctg1\t10000\t5099\t8500\t3400\t3401\t60\ttp:A:P",
                "leftclip\t9000\t0\t4000\t+\tctg1\t10000\t1999\t4950\t2950\t2951\t60\ttp:A:P",
                "rightclip\t9000\t0\t4000\t+\tctg1\t10000\t5049\t8000\t2950\t2951\t60\ttp:A:P",
            ]
        )
        evidence = read_long_read_paf(path)

        support = summarize_breakpoint(
            evidence,
            "ctg1",
            5000,
            window_bp=200,
            min_anchor_bp=1000,
        )

        self.assertEqual(support.spanning_reads, ("span",))
        self.assertEqual(support.split_reads, ("split",))
        self.assertEqual(support.left_edge_reads, ("leftclip", "split"))
        self.assertEqual(support.right_edge_reads, ("rightclip", "split"))
        self.assertEqual(support.nearby_count, 4)
        self.assertEqual(support.as_dict()["split_reads"], "1")

    def test_summarizes_contig_end_bridges(self):
        path = self.write_paf(
            [
                "bridge1\t15000\t1000\t5500\t+\tleft\t10000\t6999\t10000\t3000\t3001\t60\ttp:A:P",
                "bridge1\t15000\t5800\t10000\t+\tright\t8000\t0\t4200\t4200\t4200\t60\ttp:A:P",
                "bridge2\t15000\t8000\t12000\t-\tleft\t10000\t6500\t10000\t3500\t3500\t60\ttp:A:P",
                "bridge2\t15000\t3000\t7900\t+\tright\t8000\t0\t4900\t4900\t4900\t60\ttp:A:P",
                "offtarget\t15000\t1000\t5000\t+\tleft\t10000\t1000\t5000\t4000\t4000\t60\ttp:A:P",
            ]
        )
        evidence = read_long_read_paf(path)

        support = summarize_contig_bridge(
            evidence,
            "left",
            "right",
            terminal_window_bp=4000,
            min_anchor_bp=1000,
        )

        self.assertEqual(support.bridge_count, 2)
        self.assertEqual(support.bridge_reads, ("bridge1", "bridge2"))
        self.assertEqual(support.orientation_summary, "+/+:1,-/+:1")
        self.assertEqual(support.read_order_summary, "left_before_right:1,right_before_left:1")
        self.assertEqual(support.median_read_gap_bp, 200)
        self.assertEqual(support.as_dict()["bridge_reads"], "2")

    def test_counts_read_depth(self):
        path = self.write_paf(
            [
                "read1\t10000\t0\t5000\t+\tctg1\t9000\t999\t5000\t4000\t4001\t60\ttp:A:P",
                "read2\t10000\t0\t5000\t+\tctg1\t9000\t2999\t8000\t5000\t5001\t60\ttp:A:P",
                "read3\t10000\t0\t5000\t+\tctg2\t9000\t2999\t8000\t5000\t5001\t60\ttp:A:P",
            ]
        )
        evidence = read_long_read_paf(path)

        self.assertEqual(read_depth_at(evidence, "ctg1", 3500), 2)
        self.assertEqual(read_depth_window(evidence, "ctg1", 7500, 8500), 1)


if __name__ == "__main__":
    unittest.main()
