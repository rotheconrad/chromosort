"""Independent mixed-plan geometry and exact sequence tests, not biological scoring."""

import json
import random
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from chromosort import fix_contigs, manual, repairs, workflow
from chromosort.multireference import write_tsv
from chromosort.provenance import sha256_file
from chromosort.reference_order import iter_fasta_records, reverse_complement
from test_manifest import make_bundle


class RepairTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.manifest = make_bundle(self.root)
        rng = random.Random(193)
        self.original = {n: "".join(rng.choice("ACGT") for _ in seq)
                         for n, _, seq in iter_fasta_records(self.root / "assembly.fa")}
        (self.root / "assembly.fa").write_text("".join(f">{n}\n{s}\n" for n, s in self.original.items()))
        data = json.loads(self.manifest.read_text())
        data["assembly"]["sha256"] = sha256_file(self.root / "assembly.fa")
        for entry in data["evidence"]:
            entry["assembly_sha256"] = data["assembly"]["sha256"]
            path = self.root / entry["path"]
            if entry["reference_id"] in {"A", "Aalt"}:
                lines = [r for r in path.read_text().splitlines() if not r.startswith("chimera\t")]
                for lo, hi, chrom in [(0, 100, "chr1"), (100, 160, "chr2"), (160, 200, "chr1")]:
                    size = hi - lo
                    lines.append(f"chimera\t200\t{lo}\t{hi}\t+\t{chrom}\t1000\t{lo}\t{hi}\t{size}\t{size}\t60\tcg:Z:{size}M")
                path.write_text("\n".join(lines) + "\n")
                entry["sha256"] = sha256_file(path)
        reads = self.root / "reads.paf"
        reads.write_text("".join(f"m{i}\t100\t0\t100\t+\tchimera\t200\t40\t140\t100\t100\t60\tcg:Z:100M\n" for i in [1, 2]))
        data["evidence"].append({"evidence_id": "reads", "kind": "read_paf", "path": reads.name,
            "sha256": sha256_file(reads), "assembly_sha256": data["assembly"]["sha256"]})
        self.manifest.write_text(json.dumps(data))
        self.scan = self.root / "scan"
        self.options = ["--manifest", str(self.manifest), "--mode", "sensitive", "--min-segment-bp", "1",
                        "--min-piece-bp", "1", "--continuity-anchor-bp", "10", "--read-evidence-id", "reads"]
        workflow.main(["scan", *self.options, "--read-window-bp", "30", "--read-min-anchor-bp", "10",
                       "--output-dir", str(self.scan)])
        self.decisions = workflow.read_decisions(self.scan / "decisions.tsv")
        self.assertEqual([int(r["position"]) for r in self.decisions], [100, 160])
        for row in self.decisions:
            row.update(decision="defer" if int(row["position"]) == 100 else "accept", reviewer="scripted geometry test")
        write_tsv(self.root / "decisions.tsv", self.decisions, workflow.DECISION_COLUMNS)

    def tearDown(self):
        self.temp.cleanup()

    def edit(self, start=20, end=70, action="reverse", name="edit"):
        table = workflow.main(["edit", "--scan-dir", str(self.scan), "--event-id", self.decisions[0]["event_id"],
            "--action", action, "--start", str(start), "--end", str(end), "--notes", "Explicit scripted geometry",
            "--output-dir", str(self.root / name)])
        row = repairs.read_edits(table)[0]
        self.assertEqual(row["decision"], "pending")
        return row

    def apply(self, rows, name="applied", plan_only=False):
        write_tsv(self.root / "reviewed_edits.tsv", rows, repairs.EDIT_COLUMNS)
        return workflow.main(["apply", "--scan-dir", str(self.scan), "--decisions", str(self.root / "decisions.tsv"),
            "--reviewed-edits", str(self.root / "reviewed_edits.tsv"), "--output-dir", str(self.root / name),
            *(["--plan-only"] if plan_only else [])])

    def test_mixed_plan_explicit_reverse_and_selective_cut_preserve_native_remainder(self):
        # Reads defer the whole automatic plan; selective execution is explicitly reviewed.
        fasta = self.root / "automatic.fa"
        fix_contigs.main([*self.options, "--all", "--output-fasta", str(fasta), "--report", str(self.root / "auto.tsv")])
        self.assertEqual({n:s for n,_,s in iter_fasta_records(fasta)}, self.original)
        row = self.edit()
        row.update(decision="accept", reviewer="scripted geometry test")
        preview = self.apply([row], "preview", plan_only=True)
        self.assertFalse((preview.parent / "reviewed.fa").exists())
        result = self.apply([row])
        self.assertEqual(preview.read_bytes(), (result.parent / "recipe.json").read_bytes())
        source = self.original["chimera"]
        records = {n:s for n,_,s in iter_fasta_records(result)}
        repaired = [s for n,s in records.items() if n not in self.original]
        self.assertEqual(repaired, [source[:20] + reverse_complement(source[20:70]) + source[70:160], source[160:]])
        for name, sequence in self.original.items():
            if name != "chimera":
                self.assertEqual(records[name], sequence)
        recipe = json.loads((result.parent / "recipe.json").read_text())
        self.assertEqual([(p["source"], p["start"], p["end"]) for p in recipe["pieces"] if p["strand"] == "-"], [("chimera", 21, 70)])
        workflow.verify_outputs(self.root / "assembly.fa", result, str(result) + ".agp", True)
        replay = self.root / "replay.fa"
        manual.main(["apply", "--assembly-fasta", str(self.root / "assembly.fa"), "--recipe", str(result.parent / "recipe.json"), "-o", str(replay)])
        self.assertEqual(result.read_bytes(), replay.read_bytes())
        self.assertEqual(Path(str(result) + ".agp").read_bytes(), Path(str(replay) + ".agp").read_bytes())
        figures = [json.loads(p.read_text()) for p in (preview.parent / "figures").glob("*.figure.json")]
        self.assertEqual({f["settings"]["breakpoint"] for f in figures}, {20, 70})
        for figure in (preview.parent / "figures").glob("*.svg"):
            self.assertNotIn("cut after base", figure.read_text())
            self.assertIn("boundary", figure.read_text())

    def test_pending_stale_geometry_and_missing_review_fail_before_sequence_output(self):
        row = self.edit()
        variants = [(row, "Explicit edit decision"),
            ({**row, "decision": "accept", "reviewer": "reviewer", "assembly_sha256": "stale"}, "different assembly"),
            ({**row, "decision": "accept", "reviewer": "reviewer", "scan_sha256": "stale"}, "different assembly"),
            ({**row, "decision": "accept", "reviewer": "reviewer", "end": 71}, "geometry changed"),
            ({**row, "decision": "accept"}, "Reviewer label"),
            ({**row, "decision": "accept", "reviewer": "reviewer", "notes": ""}, "rationale"),
            ({**row, "decision": "accept", "reviewer": "reviewer", "target": "copyA"}, "same target")]
        for index, (invalid, message) in enumerate(variants):
            with self.subTest(message=message), self.assertRaisesRegex(ValueError, message):
                self.apply([invalid], f"invalid-{index}")
            self.assertFalse((self.root / f"invalid-{index}").exists())

    def test_conflicting_cuts_overlapping_reversals_and_duplicate_rows_fail(self):
        first = self.edit(20, 70)
        crossing = self.edit(150, 170, name="crossing")
        overlap = self.edit(60, 80, name="overlap")
        duplicate_cut = self.edit(160, 160, action="cut", name="duplicate-cut")
        contradictory_cut = self.edit(100, 100, action="cut", name="contradiction")
        for row in [first, crossing, overlap, duplicate_cut, contradictory_cut]:
            row.update(decision="accept", reviewer="reviewer")
        for rows, message in [([crossing], "inside a reversed"), ([first, overlap], "overlap"),
                              ([first, first], "Duplicate edit_id"), ([duplicate_cut], "same cut"),
                              ([contradictory_cut], "contradicts")]:
            with self.subTest(message=message), self.assertRaisesRegex(ValueError, message):
                self.apply(rows)
        self.assertFalse((self.root / "applied").exists())

    def test_rejected_and_deferred_edits_preserve_sequence_and_are_reported(self):
        row = self.edit()
        for decision in ["reject", "defer"]:
            edited = {**row, "decision": decision, "reviewer": "reviewer"}
            result = self.apply([edited], decision)
            records = [s for n,_,s in iter_fasta_records(result) if n not in self.original]
            self.assertEqual(records, [self.original["chimera"][:160], self.original["chimera"][160:]])
            self.assertIn("retained_" + decision, (result.parent / "decision_outcomes.tsv").read_text())

    def test_whole_contig_reverse_and_terminal_panels(self):
        for row in self.decisions:
            row["decision"] = "reject"
        write_tsv(self.root / "decisions.tsv", self.decisions, workflow.DECISION_COLUMNS)
        row = self.edit(0, 200)
        row.update(decision="accept", reviewer="reviewer")
        result = self.apply([row])
        records = {n:s for n,_,s in iter_fasta_records(result)}
        self.assertEqual(records["chimera"], reverse_complement(self.original["chimera"]))
        figures = [json.loads(p.read_text()) for p in (result.parent / "figures").glob("*.figure.json")]
        self.assertEqual(len(figures), 2)
        self.assertTrue(all(f["settings"]["breakpoint"] is None for f in figures))

    def test_adjacent_edits_and_boundary_cut_have_unambiguous_composition(self):
        rows = [self.edit(20, 40), self.edit(40, 70, name="adjacent"), self.edit(70, 70, action="cut", name="cut")]
        for row in rows:
            row.update(decision="accept", reviewer="reviewer")
        result = self.apply(rows)
        source = self.original["chimera"]
        records = [s for n,_,s in iter_fasta_records(result) if n not in self.original]
        self.assertEqual(records, [source[:20] + reverse_complement(source[20:40]) + reverse_complement(source[40:70]),
                                  source[70:160], source[160:]])


if __name__ == "__main__":
    unittest.main()
