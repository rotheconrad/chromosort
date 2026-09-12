import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from chromosort import workflow, manual, scaffold, evaluate, fix_contigs
from chromosort.review import ReviewEvent, write_review_events
from dataclasses import replace
from chromosort.multireference import write_tsv
from chromosort.reference_order import iter_fasta_records
from test_manifest import make_bundle


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.manifest = make_bundle(self.root)
        self.scan = self.root / "scan"
        workflow.main(["scan", "--manifest", str(self.manifest), "--output-dir", str(self.scan),
                       "--mode", "sensitive", "--min-segment-bp", "1", "--min-piece-bp", "1",
                       "--read-window-bp", "50", "--read-min-anchor-bp", "10"])

    def tearDown(self):
        self.temp.cleanup()

    def decisions(self, decision):
        rows = workflow.read_decisions(self.scan / "decisions.tsv")
        self.assertEqual(len(rows), 1)
        for row in rows:
            row["decision"] = decision
            row["reviewer"] = "scripted regression"
        write_tsv(self.scan / "reviewed.tsv", rows, workflow.DECISION_COLUMNS)
        return rows

    def apply(self, suffix):
        return workflow.main(["apply", "--scan-dir", str(self.scan), "--decisions", str(self.scan / "reviewed.tsv"),
                              "--output-dir", str(self.root / suffix)])

    def test_pending_is_not_an_accepted_decision(self):
        self.decisions("pending")
        with self.assertRaisesRegex(ValueError, "Explicit decision required"):
            self.apply("pending")

    def test_rejected_scaffold_join_retains_both_components(self):
        source = Path(__file__).parent / "data" / "scaffold"
        event = ReviewEvent("join", "scaffold", "scaffold_gap", "chr1", accept=False,
                            fields={"scaffold": "chr1", "left_contig": "chr1_contigA",
                                    "right_contig": "chr1_contigB", "gap_bp": 5})
        review = self.root / "joins.tsv"
        write_review_events(review, [event])
        prefix = self.root / "scaffold"
        scaffold.main(["--ordered-fasta", str(source / "ordered.fa"), "--assignments", str(source / "assignments.tsv"),
                       "--reviewed-plan", str(review), "--require-reviewed-joins", "-o", str(prefix)])
        fasta = Path(str(prefix) + ".scaffold.fa")
        self.assertEqual(len(list(iter_fasta_records(fasta))), 3)
        workflow.verify_outputs(source / "ordered.fa", fasta, Path(str(prefix) + ".scaffold.agp"), True)

    def test_partial_legacy_piece_acceptance_cannot_drop_rejected_sequence(self):
        args = evaluate.parse_fix_args(["--manifest", str(self.manifest), "--contigs", "chimera",
            "--min-segment-bp", "1", "--min-piece-bp", "1", "--mode", "sensitive", "-o", str(self.root / "eval")])
        events = evaluate.build_fix_events(args)
        events = [replace(event, accept=(i == 0)) for i, event in enumerate(events)]
        review = self.root / "partial.tsv"
        write_review_events(review, events)
        with self.assertRaisesRegex(ValueError, "exactly once"):
            fix_contigs.build_reviewed_plans(self.root / "assembly.fa", review, args)

    def test_selective_cut_and_replay_preserve_every_base(self):
        self.decisions("accept")
        fasta = self.apply("accepted")
        records = list(iter_fasta_records(fasta))
        self.assertEqual(len(records), 6)
        self.assertEqual(sum(len(seq) for _, _, seq in records), 580)
        self.assertTrue(all(row["omitted_bp"] == row["duplicated_bp"] == 0 for row in
                            json.loads((fasta.parent / "apply.json").read_text())["accounting"]))
        manual.main(["apply", "--assembly-fasta", str(self.root / "assembly.fa"), "--recipe", str(fasta.parent / "recipe.json"),
                     "-o", str(self.root / "replay.fa")])
        self.assertEqual(fasta.read_bytes(), (self.root / "replay.fa").read_bytes())

    def test_rejected_cut_keeps_original_contig(self):
        self.decisions("reject")
        fasta = self.apply("rejected")
        records = {name: seq for name, _, seq in iter_fasta_records(fasta)}
        original = {name: seq for name, _, seq in iter_fasta_records(self.root / "assembly.fa")}
        self.assertEqual(records, original)

    def test_refinement_requires_explicit_decision_and_stale_inputs_fail(self):
        rows = self.decisions("accept")
        rows[0]["position"] = 90
        write_tsv(self.scan / "reviewed.tsv", rows, workflow.DECISION_COLUMNS)
        with self.assertRaisesRegex(ValueError, "requires decision=refine"):
            self.apply("invalid")
        rows[0]["decision"] = "refine"
        write_tsv(self.scan / "reviewed.tsv", rows, workflow.DECISION_COLUMNS)
        self.apply("refined")
        fasta = self.root / "assembly.fa"
        fasta.write_text(fasta.read_text().replace("ACGT", "TCGT", 1))
        with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
            self.apply("stale")


if __name__ == "__main__":
    unittest.main()
