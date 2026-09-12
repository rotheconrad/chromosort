import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from chromosort import reference_order, scaffold, evaluate, manual, plot
from chromosort.manifest import InputBundle
from chromosort.provenance import sha256_file


def make_bundle(root):
    """Controlled geometry only, not a biological performance benchmark."""
    root = Path(root)
    assembly = root / "assembly.fa"
    assembly.write_text("".join(f">{name}\n{'ACGT' * (length // 4)}\n" for name, length in
                                [("copyA", 100), ("copyB", 100), ("ambiguous", 100), ("novel", 80), ("chimera", 200)]))
    references, labels, evidence = [], [], []
    for rid, group in [("A", "A"), ("B", "B"), ("Aalt", "A")]:
        fasta = root / (rid + ".fa")
        fasta.write_text(">chr1\n" + "ACGT" * 250 + "\n>chr2\n" + "ACGT" * 250 + "\n")
        references.append({"reference_id": rid, "fasta": fasta.name, "sha256": sha256_file(fasta),
                           "role": "alternative" if rid == "Aalt" else "backbone"})
        for chrom in ("chr1", "chr2"):
            labels.append({"reference_id": rid, "sequence_id": chrom, "chromosome": chrom,
                           "output_group": "target", "subgenome": group, "copy_label": "1",
                           "homoeolog_group": chrom, "backbone": rid != "Aalt"})
        lines = []
        for name, identity in [("copyA", 100 if group == "A" else 70),
                               ("copyB", 100 if group == "B" else 70), ("ambiguous", 100)]:
            lines.append(f"{name}\t100\t0\t100\t+\tchr1\t1000\t100\t200\t{identity}\t100\t60\tcg:Z:100M\ttp:A:P\n")
        if group == "A":
            for start, chrom in [(0, "chr1"), (100, "chr2")]:
                lines.append(f"chimera\t200\t{start}\t{start+100}\t+\t{chrom}\t1000\t500\t600\t100\t100\t60\tcg:Z:100M\ttp:A:P\n")
        paf = root / (rid + ".paf")
        paf.write_text("".join(lines))
        evidence.append({"evidence_id": rid, "kind": "assembly_paf", "path": paf.name,
                         "sha256": sha256_file(paf), "assembly_sha256": sha256_file(assembly),
                         "reference_id": rid, "reference_sha256": sha256_file(fasta),
                         "coordinate_system": "PAF", "command": ["synthetic geometry fixture"],
                         "software": {"fixture": "1"}, "evidence_class": "controlled_geometry"})
    manifest = root / "manifest.json"
    manifest.write_text(json.dumps({"schema": "chromosort-input-v1", "assembly": {
        "path": assembly.name, "sha256": sha256_file(assembly), "sample_id": "target", "stage": "input"},
        "references": references, "reference_sequences": labels, "evidence": evidence}, indent=2) + "\n")
    return manifest


class ManifestTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.manifest = make_bundle(self.root)

    def tearDown(self):
        self.tmp.cleanup()

    def sort(self, *extra):
        reference_order.main(["--manifest", str(self.manifest), "-o", str(self.root / "sorted"),
                              "--min-aligned-bp", "1", "--min-query-cov", "0.4", "--retain-all", *extra])
        with open(self.root / "sorted.contig_assignments.tsv") as handle:
            return {r["contig"]: r for r in csv.DictReader(handle, delimiter="\t")}

    def test_copy_retention_ambiguity_and_exact_once(self):
        rows = self.sort()
        self.assertEqual(rows["copyA"]["partition"], "placed")
        self.assertEqual(rows["copyB"]["partition"], "placed")
        self.assertEqual(rows["copyA"]["subgenome"], "A")
        self.assertEqual(rows["copyB"]["subgenome"], "B")
        self.assertEqual(rows["ambiguous"]["partition"], "ambiguous")
        self.assertEqual(rows["novel"]["partition"], "unplaced")
        self.assertAlmostEqual(float(rows["copyA"]["decision_margin"]), 0.3)
        records = list(reference_order.iter_fasta_records(self.root / "sorted.ordered.fa"))
        self.assertEqual(len(records), 5)
        self.assertEqual(sum(len(seq) for _, _, seq in records), 580)
        self.assertEqual(len({name for name, _, _ in records}), 5)
        scaffold.main(["--ordered-fasta", str(self.root / "sorted.ordered.fa"), "--assignments",
                       str(self.root / "sorted.contig_assignments.tsv"), "-o", str(self.root / "scaffold")])
        outputs = dict((name, seq) for name, _, seq in reference_order.iter_fasta_records(self.root / "scaffold.scaffold.fa"))
        self.assertIn("target::chr1::A::1", outputs)
        self.assertIn("target::chr1::B::1", outputs)
        self.assertTrue(any(name.startswith("ambiguous::") for name in outputs))

    def test_stale_evidence_and_nonunique_backbone(self):
        data = json.loads(self.manifest.read_text())
        data["evidence"][0]["assembly_sha256"] = "0" * 64
        self.manifest.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, "different assembly stage"):
            InputBundle(self.manifest)
        data["evidence"][0]["assembly_sha256"] = data["assembly"]["sha256"]
        data["reference_sequences"][-1]["backbone"] = True
        self.manifest.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, "Multiple coordinate backbones"):
            InputBundle(self.manifest)

    def test_alternative_reference_multiplicity_does_not_change_margin(self):
        before = self.sort()
        data = json.loads(self.manifest.read_text())
        data["references"] = [r for r in data["references"] if r["reference_id"] != "Aalt"]
        data["reference_sequences"] = [r for r in data["reference_sequences"] if r["reference_id"] != "Aalt"]
        data["evidence"] = [r for r in data["evidence"] if r["reference_id"] != "Aalt"]
        self.manifest.write_text(json.dumps(data))
        after = self.sort()
        for name in before:
            self.assertEqual(before[name]["decision_margin"], after[name]["decision_margin"])
            self.assertEqual(before[name]["partition"], after[name]["partition"])

    def test_overlapping_copies_are_retained_and_withheld_from_scaffolding(self):
        for rid in ("A", "Aalt", "B"):
            paf = self.root / (rid + ".paf")
            text = paf.read_text()
            if rid == "B":
                text = "\n".join(line for line in text.splitlines() if not line.startswith("ambiguous\t")) + "\n"
            paf.write_text(text)
        data = json.loads(self.manifest.read_text())
        for entry in data["evidence"]:
            entry["sha256"] = sha256_file(self.root / entry["path"])
        self.manifest.write_text(json.dumps(data))
        rows = self.sort()
        for name in ("copyA", "ambiguous"):
            self.assertEqual(rows[name]["partition"], "placed")
            self.assertEqual(rows[name]["kept"], "yes")
            self.assertEqual(rows[name]["scaffold_eligible"], "no")
        self.assertEqual(len(list(reference_order.iter_fasta_records(self.root / "sorted.ordered.fa"))), 5)

    def test_scaffold_rejects_mixed_backbones_and_supplied_copy_groups(self):
        records = [scaffold.FastaRecord(name, ">" + name, "ACGT") for name in ("a", "b")]
        rows = {name: scaffold.AssignmentRow(name, name, "chr1", i, i * 10, i * 10 + 3,
                                            backbone=backbone, copy_label="1")
                for i, (name, backbone) in enumerate((("a", "refA"), ("b", "refB")))}
        with self.assertRaisesRegex(ValueError, "mix reference coordinate backbones"):
            scaffold.group_scaffold_members(records, rows)
        rows["b"].backbone = "refA"
        rows["b"].copy_label = "2"
        with self.assertRaisesRegex(ValueError, "distinct supplied subgenome/copy groups"):
            scaffold.group_scaffold_members(records, rows)

    def test_alternative_reference_never_creates_fix_transitions(self):
        args = evaluate.parse_fix_args(["--manifest", str(self.manifest), "--contigs", "copyA", "chimera",
            "--min-segment-bp", "1", "--min-piece-bp", "1", "--mode", "sensitive", "-o", str(self.root / "eval")])
        events = evaluate.build_fix_events(args)
        copy = [e for e in events if e.target == "copyA"]
        self.assertEqual([e.action for e in copy], ["no_split"])
        split = [e for e in events if e.target == "chimera"]
        self.assertEqual(len(split), 2)
        self.assertTrue(all(e.fields["reference_id"] == "A" for e in split))

    def test_manifest_dashboard_and_dotplot(self):
        args = manual.parse_generate_args(["--manifest", str(self.manifest), "-o", str(self.root / "review.html")])
        data = manual.build_dashboard_data(args)
        self.assertEqual(len(data["refRecords"]), 6)
        self.assertTrue(data["inputs"]["assemblySha256"])
        manual.write_dashboard(self.root / "review.html", data)
        plot.main(["--manifest", str(self.manifest), "-o", str(self.root / "dotplot"), "--formats", "svg"])
        self.assertTrue((self.root / "dotplot.svg").is_file())

    def test_missing_backbone_abstains_and_retains(self):
        data = json.loads(self.manifest.read_text())
        for row in data["reference_sequences"]:
            row["backbone"] = False
        self.manifest.write_text(json.dumps(data))
        rows = self.sort()
        self.assertEqual(rows["copyA"]["status"], "insufficient_backbone_alignment")
        self.assertEqual(rows["copyA"]["assigned_ref"], ".")
        self.assertEqual(rows["copyA"]["kept"], "yes")


if __name__ == "__main__":
    unittest.main()
