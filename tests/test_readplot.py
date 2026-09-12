import json
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from chromosort import readplot
from chromosort.longreads import summarize_breakpoint, read_depth_at
from chromosort.provenance import sha256_file
from chromosort.readalign import load_reads
from test_manifest import make_bundle


def add_read_fixture(manifest):
    root = manifest.parent
    sam = root / "reads.sam"
    sam.write_text("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:copyA\tLN:100\n" + "\n".join([
        "span\t0\tcopyA\t11\t60\t10S80M10S\t*\t0\t0\t*\t*\tNM:i:0",
        "split\t0\tcopyA\t1\t60\t45M55S\t*\t0\t0\t*\t*\tSA:Z:copyA,55,+,54S46M,60,0;",
        "split\t2048\tcopyA\t55\t60\t54S46M\t*\t0\t0\t*\t*\tNM:i:0",
        "deletion\t16\tcopyA\t1\t60\t40M20D40M\t*\t0\t0\t*\t*\tNM:i:20",
        "secondary\t256\tcopyA\t1\t60\t100M\t*\t0\t0\t*\t*",
        "lowmapq\t0\tcopyA\t1\t5\t100M\t*\t0\t0\t*\t*",
        "tiny\t0\tcopyA\t45\t60\t10M\t*\t0\t0\t*\t*",
    ]) + "\n")
    data = json.loads(manifest.read_text())
    data["evidence"].append({"evidence_id": "reads", "kind": "sam", "path": sam.name,
                             "sha256": sha256_file(sam), "assembly_sha256": data["assembly"]["sha256"],
                             "read_sample": "target", "readset_id": "reads", "evidence_class": "controlled_geometry"})
    features = root / "features.tsv"
    features.write_text("contig\tstart\tend\tlabel\ncopyA\t35\t60\tSynthetic feature\n")
    data["evidence"].append({"evidence_id": "features", "kind": "features", "path": features.name,
                             "sha256": sha256_file(features), "assembly_sha256": data["assembly"]["sha256"]})
    variants = root / "variants.vcf"
    variants.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                        "copyA\t30\t.\tA\tT\t.\tPASS\t.\n"
                        "copyA\t50\t.\tA\tAT\t.\tPASS\t.\n"
                        "copyA\t65\t.\tA\t<DEL>\t.\tPASS\tEND=70\n")
    data["evidence"].append({"evidence_id": "variants", "kind": "variants", "path": variants.name,
                             "sha256": sha256_file(variants), "assembly_sha256": data["assembly"]["sha256"]})
    manifest.write_text(json.dumps(data))
    return sam


class ReadPlotTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.manifest = make_bundle(self.root)
        self.sam = add_read_fixture(self.manifest)

    def tearDown(self):
        self.temp.cleanup()

    def test_supplementary_cigar_and_no_artificial_split(self):
        reads = load_reads(self.sam, "sam", min_mapq=20)
        self.assertEqual(reads.counts["secondary_excluded"], 1)
        self.assertEqual(reads.counts["mapq_excluded"], 1)
        self.assertEqual(reads.counts["supplementary_retained"], 1)
        support = summarize_breakpoint(reads.evidence, "copyA", 50, 10, 10)
        self.assertEqual(support.spanning_reads, ("span",))
        self.assertEqual(support.split_reads, ("split",))
        self.assertNotIn("deletion", support.spanning_reads)
        self.assertNotIn("tiny", support.split_reads)
        self.assertEqual(read_depth_at(reads.evidence, "copyA", 50), 2)

    def test_native_reverse_packing_and_coverage(self):
        alignments = load_reads(self.sam, "sam", min_mapq=20).evidence.alignments
        forward = readplot.Region("copyA", 20, 80)
        reverse = readplot.Region("copyA", 20, 80, True)
        self.assertEqual(forward.axis(20), 0)
        self.assertEqual(reverse.axis(20), 60)
        self.assertEqual(readplot.coverage_bins(alignments, forward, 10), readplot.coverage_bins(alignments, reverse, 10))
        for region in (forward, reverse):
            packed = readplot.pack_reads(alignments, region, 100, 100)
            self.assertEqual({name for _, name, _ in packed}, {a.read for a in alignments})
            occupied = {}
            for row, name, segments in packed:
                starts = [region.span(a.contig_start - 1, a.contig_end)[0] for a in segments]
                end = max(region.span(a.contig_start - 1, a.contig_end)[1] for a in segments)
                self.assertGreater(min(starts), occupied.get(row, -1))
                occupied[row] = end

    def test_manifest_rejects_incompatible_sam_header_outside_display_window(self):
        original = self.sam.read_text()
        for header, message in [("@SQ\tSN:foreign\tLN:100\n", "target absent"),
                                ("@SQ\tSN:chimera\tLN:201\n", "target length differs")]:
            with self.subTest(header=header):
                self.sam.write_text(header + original)
                data = json.loads(self.manifest.read_text())
                next(e for e in data["evidence"] if e["evidence_id"] == "reads")["sha256"] = sha256_file(self.sam)
                self.manifest.write_text(json.dumps(data))
                with self.assertRaisesRegex(ValueError, message):
                    readplot.main(["--manifest", str(self.manifest), "--evidence-id", "reads",
                                   "--contig", "copyA", "-o", str(self.root / "invalid")])

    def test_static_exports_replay_and_counts(self):
        options = ["--manifest", str(self.manifest), "--evidence-id", "reads", "--contig", "copyA",
                   "--start", "20", "--end", "80", "--breakpoint", "50", "--min-anchor-bp", "10",
                   "--read-window-bp", "10", "--bin-bp", "5", "--event-id", "synthetic-cut", "--reverse"]
        summary = readplot.main([*options, "-o", str(self.root / "panel")])
        self.assertEqual(summary["spanning_reads"], 1)
        self.assertEqual(summary["split_reads"], 1)
        self.assertEqual(summary["molecules_displayed"], 4)
        self.assertEqual(summary["alignments_displayed"], 5)
        svg = (self.root / "panel.svg").read_text()
        for label in ("SNP", "INDEL", "SV", "Synthetic feature"):
            self.assertIn(label, svg)
        ET.parse(self.root / "panel.svg")
        self.assertTrue((self.root / "panel.pdf").read_bytes().startswith(b"%PDF"))
        self.assertTrue((self.root / "panel.png").read_bytes().startswith(b"\x89PNG"))
        readplot.main([*options, "-o", str(self.root / "replay")])
        self.assertEqual((self.root / "panel.svg").read_bytes(), (self.root / "replay.svg").read_bytes())
        readplot.main(["replay", "--figure-json", str(self.root / "panel.figure.json"), "-o", str(self.root / "saved-replay")])
        self.assertEqual((self.root / "panel.pdf").read_bytes(), (self.root / "saved-replay.pdf").read_bytes())

    def test_paf_missing_detail_is_unavailable_and_cigar_gaps_exclude_coverage(self):
        paf = self.root / "reads.paf"
        paf.write_text("plain\t100\t0\t100\t+\tcopyA\t100\t0\t100\t100\t100\t60\n"
                       "deletion\t80\t0\t80\t-\tcopyA\t100\t0\t100\t80\t100\t60\tcg:Z:40M20D40M\n")
        reads = load_reads(paf)
        self.assertEqual(reads.counts["cigar_unavailable"], 1)
        self.assertIsNone(reads.evidence.alignments[0].is_supplementary)
        self.assertIsNone(reads.evidence.alignments[0].clip_left)
        self.assertEqual(read_depth_at(reads.evidence, "copyA", 50), 1)
        self.assertEqual(read_depth_at(reads.evidence, "copyA", 30), 2)


if __name__ == "__main__":
    unittest.main()
