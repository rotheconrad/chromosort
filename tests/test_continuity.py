"""Independent geometry fixtures for review deferral, never biological truth labels."""

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from dataclasses import replace

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from chromosort import continuity, evaluate, fix_contigs, workflow
from chromosort.manifest import InputBundle
from chromosort.multireference import write_tsv
from chromosort.provenance import sha256_file
from chromosort.reference_order import iter_fasta_records
from chromosort.review import write_review_events
from test_manifest import make_bundle


class ContinuityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.manifest = make_bundle(self.root)
        self.args = SimpleNamespace(continuity_min_mapq=20, continuity_anchor_bp=10,
            continuity_min_spanning=2, read_continuity_policy="review", read_evidence_id=None)

    def tearDown(self):
        self.temp.cleanup()

    def reads(self, lines):
        paf = self.root / "reads.paf"
        paf.write_text("".join(lines))
        data = json.loads(self.manifest.read_text())
        data["evidence"] = [e for e in data["evidence"] if e["kind"] != "read_paf"]
        data["evidence"].append({"evidence_id": "reads", "kind": "read_paf", "path": paf.name,
            "sha256": sha256_file(paf), "assembly_sha256": data["assembly"]["sha256"], "readset_id": "one-source"})
        self.manifest.write_text(json.dumps(data))
        return InputBundle(self.manifest)

    @staticmethod
    def line(name, cigar="120M", mapq=60):
        query = 100 if cigar == "50M20D50M" else 120
        tag = "\tcg:Z:" + cigar if cigar else ""
        return f"{name}\t{query}\t0\t{query}\t+\tchimera\t200\t40\t160\t{query}\t120\t{mapq}{tag}\n"

    def plans(self):
        args = fix_contigs.parse_args(["--manifest", str(self.manifest), "--all", "--mode", "sensitive",
            "--min-segment-bp", "1", "--min-piece-bp", "1", "--output-fasta", str(self.root / "fixed.fa"),
            "--report", str(self.root / "fix.tsv")])
        blocks = fix_contigs.collect_blocks(self.manifest, "manifest", None, 1, 0, 50000)
        return fix_contigs.build_plans(self.root / "assembly.fa", ["chimera"], blocks, args)

    def test_continuity_defers_whole_plan_and_preserves_review_proposal(self):
        bundle = self.reads([self.line("molecule1"), self.line("molecule2")])
        plans = self.plans()
        original_pieces = list(plans["chimera"].pieces)
        rows = continuity.assess_plans(plans, bundle, self.args)
        self.assertEqual(plans["chimera"].pieces, original_pieces)
        self.assertEqual(rows[0]["spanning_molecules"], 2)
        continuity.assess_plans(plans, bundle, self.args, apply=True)
        self.assertEqual(plans["chimera"].status, "not_split_read_continuity")
        args = evaluate.parse_fix_args(["--manifest", str(self.manifest), "--all", "--mode", "sensitive",
            "--min-segment-bp", "1", "--min-piece-bp", "1", "--continuity-anchor-bp", "10", "-o", str(self.root / "review")])
        events = evaluate.build_fix_events(args)
        cut = next(e for e in events if e.action == "split_piece")
        self.assertEqual(cut.fields["read_continuity_spanning_molecules"], 2)
        self.assertIn("requires review", cut.reason)
        self.assertFalse(cut.accept)
        self.assertTrue(all(not e.accept for e in events if e.target == "chimera"))

    def test_duplicate_segments_of_one_molecule_do_not_trigger_threshold(self):
        bundle = self.reads([self.line("same"), self.line("same")])
        plans = self.plans()
        rows = continuity.assess_plans(plans, bundle, self.args, apply=True)
        self.assertEqual(rows[0]["spanning_molecules"], 1)
        self.assertEqual(plans["chimera"].status, "split")

    def test_one_conflicted_boundary_defers_entire_multi_boundary_plan(self):
        bundle = self.reads([self.line("first"), self.line("second")])
        plans = self.plans()
        first, last = plans["chimera"].pieces
        plans["chimera"].pieces = [replace(first, slice_start=0, slice_end=100),
                                   replace(last, slice_start=100, slice_end=180),
                                   replace(last, slice_start=180, slice_end=200, orientation="-")]
        rows = continuity.assess_plans(plans, bundle, self.args, apply=True)
        self.assertEqual([r["spanning_molecules"] for r in rows], [2, 0])
        self.assertTrue(all(r["automatic_action"] == "defer_contig_plan_for_review" for r in rows))
        self.assertEqual(plans["chimera"].pieces, [])

    def test_explicit_reviewed_plan_can_override_advisory_deferral(self):
        self.reads([self.line("first"), self.line("second")])
        args = evaluate.parse_fix_args(["--manifest", str(self.manifest), "--all", "--mode", "sensitive",
            "--min-segment-bp", "1", "--min-piece-bp", "1", "--continuity-anchor-bp", "10", "-o", str(self.root / "eval")])
        events = [replace(event, accept=True) for event in evaluate.build_fix_events(args)]
        reviewed = self.root / "reviewed.tsv"
        write_review_events(reviewed, events)
        fasta = self.root / "explicit.fa"
        fix_contigs.main(["--assembly-fasta", str(self.root / "assembly.fa"), "--reviewed-plan", str(reviewed),
                         "--output-fasta", str(fasta), "--report", str(self.root / "explicit.tsv")])
        self.assertEqual(len(list(iter_fasta_records(fasta))), 6)
        workflow.verify_outputs(self.root / "assembly.fa", fasta, Path(str(fasta) + ".agp"), True)

    def test_cigar_deletions_and_missing_cigar_cannot_claim_continuity(self):
        bundle = self.reads([self.line("deletion1", "50M20D50M"), self.line("deletion2", "50M20D50M"), self.line("span_only", "")])
        rows = continuity.assess_plans(self.plans(), bundle, self.args, apply=True)
        self.assertEqual(rows[0]["spanning_molecules"], 0)
        self.assertEqual(rows[0]["alignments_without_cigar"], 1)
        bundle = self.reads([self.line("span_only", "")])
        self.assertEqual(continuity.assess_plans(self.plans(), bundle, self.args)[0]["continuity_status"], "unavailable_cigar")

    def test_mapq_exclusion_and_missing_bridges_remain_neutral(self):
        bundle = self.reads([self.line("low", mapq=5), self.line("unknown", mapq=255)])
        plans = self.plans()
        rows = continuity.assess_plans(plans, bundle, self.args, apply=True)
        self.assertEqual(rows[0]["continuity_status"], "no_spanning_conflict")
        self.assertEqual(plans["chimera"].status, "split")

    def test_manifest_rejects_incompatible_paf_targets_before_mapq_filtering(self):
        original = self.line("low", mapq=5)
        for row, message in [(original.replace("\tchimera\t", "\tforeign\t"), "target absent"),
                             (original.replace("\tchimera\t200\t", "\tchimera\t201\t"), "target length differs")]:
            with self.subTest(message=message):
                bundle = self.reads([row])
                with self.assertRaisesRegex(ValueError, message):
                    continuity.assess_plans(self.plans(), bundle, self.args)

    def test_correlated_formats_not_pooled_and_independent_sources_require_selection(self):
        self.reads([self.line("one")])
        data = json.loads(self.manifest.read_text())
        entry = dict(data["evidence"][-1], evidence_id="representation2")
        data["evidence"].append(entry)
        self.manifest.write_text(json.dumps(data))
        rows = continuity.assess_plans(self.plans(), InputBundle(self.manifest), self.args)
        self.assertEqual(rows[0]["spanning_molecules"], 1)
        data["evidence"][-1]["readset_id"] = "independent-source"
        self.manifest.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, "Select one read source"):
            continuity.assess_plans(self.plans(), InputBundle(self.manifest), self.args)
        self.args.read_evidence_id = "reads"
        self.assertEqual(continuity.assess_plans(self.plans(), InputBundle(self.manifest), self.args)[0]["evidence_id"], "reads")

    def test_default_and_report_only_cli_preserve_sequence_with_different_cut_decisions(self):
        self.reads([self.line("molecule1"), self.line("molecule2")])
        for policy, expected in [("review", 5), ("report", 6)]:
            fasta = self.root / (policy + ".fa")
            report = self.root / (policy + ".tsv")
            fix_contigs.main(["--manifest", str(self.manifest), "--all", "--mode", "sensitive",
                "--min-segment-bp", "1", "--min-piece-bp", "1", "--continuity-anchor-bp", "10",
                "--read-continuity-policy", policy, "--output-fasta", str(fasta), "--report", str(report)])
            self.assertEqual(len(list(iter_fasta_records(fasta))), expected)
            workflow.verify_outputs(self.root / "assembly.fa", fasta, Path(str(fasta) + ".agp"), True)
            with open(str(report) + ".breakpoint_evidence.tsv") as handle:
                self.assertEqual(next(csv.DictReader(handle, delimiter="\t"))["spanning_molecules"], "2")

    def test_small_gap_inspection_has_two_edge_panels_and_no_cut(self):
        data = json.loads(self.manifest.read_text())
        for entry in data["evidence"]:
            paf = self.root / entry["path"]
            if entry["reference_id"] in {"A", "Aalt"}:
                lines = [line for line in paf.read_text().splitlines() if not line.startswith("chimera\t")]
                lines += ["chimera\t200\t0\t60\t+\tchr1\t1000\t500\t560\t60\t60\t60\tcg:Z:60M",
                          "chimera\t200\t100\t200\t+\tchr1\t1000\t600\t700\t100\t100\t60\tcg:Z:100M"]
                paf.write_text("\n".join(lines) + "\n")
                entry["sha256"] = sha256_file(paf)
        self.manifest.write_text(json.dumps(data))
        out = self.root / "inspection"
        workflow.main(["scan", "--manifest", str(self.manifest), "--output-dir", str(out),
            "--min-segment-bp", "50", "--min-piece-bp", "200", "--inspect-min-gap-bp", "20", "--read-window-bp", "50"])
        rows = workflow.read_decisions(out / "decisions.tsv")
        self.assertEqual([r["action"] for r in rows], ["inspect"])
        figures = [json.loads(p.read_text()) for p in (out / "figures").glob("*.figure.json")]
        self.assertEqual({f["settings"]["breakpoint"] for f in figures}, {60, 100})
        for row in rows:
            row.update(decision="accept", reviewer="scripted inspection check")
        write_tsv(out / "accepted.tsv", rows, workflow.DECISION_COLUMNS)
        fasta = workflow.main(["apply", "--scan-dir", str(out), "--decisions", str(out / "accepted.tsv"), "--output-dir", str(self.root / "applied")])
        self.assertEqual({n:s for n,_,s in iter_fasta_records(fasta)}, {n:s for n,_,s in iter_fasta_records(self.root / "assembly.fa")})


if __name__ == "__main__":
    unittest.main()
