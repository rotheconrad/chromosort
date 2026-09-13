"""Independent coordinate fixtures, not measured development outcomes."""

import itertools
import json
import random
import sys
import tempfile
import unittest
from collections import defaultdict
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from chromosort import evaluate, fix_contigs as fix, manual, repairs, review_proposals as review, workflow
from chromosort.manifest import InputBundle
from chromosort.multireference import write_tsv
from chromosort.provenance import sha256_file
from chromosort.reference_order import iter_fasta_records, reverse_complement


def make_inspection_bundle(root, redundant=False, second_frame=False, graph=False, reads=False):
    rng = random.Random(173)
    sequence = "".join(rng.choice("ACGT") for _ in range(1000))
    assembly = root / "assembly.fa"
    assembly.write_text("".join(f">{name}\n{sequence}\n" for name in ["mixed", "clean", "reverse", "ambiguous"]))
    references, labels, evidence = [], [], []
    geometry = [(0, 300, "+"), (300, 380, "-"), (380, 700, "+"), (700, 850, "-"), (850, 1000, "+")]
    for rid in (["A", "B"] if second_frame else ["A"]):
        ref = root / (rid + ".fa")
        ref.write_text(">chr\n" + sequence * 2 + "\n")
        references.append({"reference_id": rid, "fasta": ref.name, "sha256": sha256_file(ref), "role": "backbone"})
        labels.append({"reference_id": rid, "sequence_id": "chr", "chromosome": "chr", "output_group": rid, "backbone": True})
        rows = []
        for name, spans in [("mixed", geometry), ("clean", [(0, 1000, "+")]),
                            ("reverse", [(0, 1000, "-")]), ("ambiguous", [(0, 500, "+"), (500, 1000, "-")])]:
            for start, end, orientation in spans:
                if rid == "B":
                    orientation = "+" if orientation == "-" else "-"
                n = end - start
                rows.append(f"{name}\t1000\t{start}\t{end}\t{orientation}\tchr\t2000\t{start}\t{end}\t{n}\t{n}\t60\tcg:Z:{n}M\ttp:A:P\n")
        if redundant:
            rows += [rows[1], "mixed\t1000\t310\t370\t-\tchr\t2000\t310\t370\t60\t60\t60\tcg:Z:60M\ttp:A:P\n",
                     "mixed\t1000\t280\t320\t+\tchr\t2000\t280\t320\t40\t40\t60\tcg:Z:40M\ttp:A:P\n"]
        paf = root / (rid + ".paf")
        paf.write_text("".join(rows))
        evidence.append({"evidence_id": rid, "kind": "assembly_paf", "path": paf.name,
            "sha256": sha256_file(paf), "assembly_sha256": sha256_file(assembly),
            "reference_id": rid, "reference_sha256": sha256_file(ref), "dependency_group": "shared_geometry"})
    if reads:
        path = root / "reads.paf"
        path.write_text("".join(f"{name}\t100\t0\t100\t+\tmixed\t1000\t{start}\t{start+100}\t100\t100\t60\tcg:Z:100M\n"
            for name, start in [("edge", 250), ("edge", 250), ("interior", 310)]))
        evidence.append({"evidence_id": "reads", "kind": "read_paf", "path": path.name,
            "sha256": sha256_file(path), "assembly_sha256": sha256_file(assembly), "dependency_group": "readset"})
    if graph:
        path = root / "graph.gfa"
        path.write_text("S\tmixed\t" + sequence + "\n")
        evidence.append({"evidence_id": "graph", "kind": "gfa", "path": path.name,
            "sha256": sha256_file(path), "assembly_sha256": sha256_file(assembly)})
    manifest = root / "manifest.json"
    manifest.write_text(json.dumps({"schema": "chromosort-input-v1", "assembly": {"path": assembly.name,
        "sha256": sha256_file(assembly), "stage": "fixture"}, "references": references,
        "reference_sequences": labels, "evidence": evidence}))
    return manifest


def settings(manifest):
    planner = evaluate.parse_fix_args(["--manifest", str(manifest), "--all", "--mode", "comprehensive",
        "--min-segment-bp", "100", "--min-piece-bp", "1", "-o", str(manifest.parent / "unused")])
    args = SimpleNamespace(**{**vars(planner), "inspect_orientation": True, "inspect_min_orientation_bp": 50,
        "alternative_plans": 5, "review_max_blocks": 128, "read_min_anchor_bp": 10,
        "continuity_anchor_bp": 10, "min_mapq": 20, "contigs": None})
    return args, planner


class InspectionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def build(self, **options):
        manifest = make_inspection_bundle(self.root, **options)
        args, planner = settings(manifest)
        return review.build(InputBundle(manifest), args, planner, [])

    def test_short_and_smoothed_intervals_and_orientation_controls(self):
        record = self.build()
        events = [e for e in record["events"] if e["kind"] == "orientation_discrepancy"]
        self.assertEqual([(e["target"], e["start0"], e["end0"]) for e in events], [("mixed", 300, 380), ("mixed", 700, 850)])
        self.assertIn("below_correction_segment_floor", events[0]["loss_stages"])
        self.assertIn("smoothed_by_objective", events[1]["loss_stages"])
        contexts = {c["target"]: c["orientation"]["frames"][0]["status"] for c in record["contigs"]}
        self.assertEqual(contexts["reverse"], "whole_record_orientation")
        self.assertEqual(contexts["clean"], "uniform_forward_orientation")
        self.assertEqual(contexts["ambiguous"], "ambiguous_terminal_orientation")
        self.assertTrue(all(e["decision"] == "pending" and not e["implicit_edits"] for e in record["events"]))

    def test_redundant_nested_and_opposing_overlap_keep_one_interval(self):
        record = self.build(redundant=True)
        events = [e for e in record["events"] if e["kind"] == "orientation_discrepancy"]
        self.assertEqual(len(events), 2)
        first = events[0]
        self.assertEqual(first["positions"], [300, 380])
        self.assertEqual(first["frame_contexts"][0]["supporting_row_count"], 3)
        self.assertEqual(first["frame_contexts"][0]["supporting_union_bp"], 80)
        self.assertEqual(first["frame_contexts"][0]["supporting_additive_bp"], 220)
        self.assertEqual(first["frame_contexts"][0]["opposing_overlap_bp"], 20)
        self.assertIn("overlapping_opposite_orientation_evidence", first["loss_stages"])
        # Row order and identical duplicate count cannot change event geometry/ID.
        paf = self.root / "A.paf"
        paf.write_text("".join(reversed(paf.read_text().splitlines(keepends=True))))
        manifest = self.root / "manifest.json"
        data = json.loads(manifest.read_text())
        data["evidence"][0]["sha256"] = sha256_file(paf)
        manifest.write_text(json.dumps(data))
        args, planner = settings(manifest)
        repeated = review.build(InputBundle(manifest), args, planner, [])
        self.assertEqual([e["event_id"] for e in events], [e["event_id"] for e in repeated["events"] if e["kind"] == "orientation_discrepancy"])

    def test_opposite_reference_frames_are_context_not_independent_events(self):
        record = self.build(second_frame=True)
        events = [e for e in record["events"] if e["kind"] == "orientation_discrepancy"]
        self.assertEqual(len(events), 2)
        self.assertTrue(all(len(e["reference_frames"]) == 2 for e in events))
        self.assertEqual({c["anchor_orientation"] for c in events[0]["frame_contexts"]}, {"+", "-"})
        self.assertIn("withheld_cross_reference_overlap", events[0]["loss_stages"])
        self.assertTrue(all(not b["correction_eligible"] for b in record["blocks"]))

    def test_read_and_graph_evidence_are_on_actual_edges_and_deduplicate_names(self):
        record = self.build(reads=True, graph=True)
        event = next(e for e in record["events"] if e["positions"] == [300, 380])
        edges = event["coordinate_evidence"]
        self.assertEqual([e["position"] for e in edges], [300, 380])
        self.assertEqual(edges[0]["reads"]["spanning_molecules"], 1)
        self.assertEqual(edges[1]["reads"]["spanning_molecules"], 1)
        self.assertEqual(edges[0]["reads"]["spanning_read_ids"], "edge")
        self.assertEqual(edges[1]["reads"]["spanning_read_ids"], "interior")
        self.assertEqual(edges[0]["graph"]["sides"]["left_base"]["unitig_offset_0"], 299)
        self.assertEqual(edges[0]["graph"]["sides"]["right_base"]["unitig_offset_0"], 300)
        self.assertEqual(edges[1]["graph"]["sides"]["right_base"]["unitig_offset_0"], 380)

    def test_limits_are_explicit_and_never_truncate_correction_inputs(self):
        manifest = make_inspection_bundle(self.root)
        args, planner = settings(manifest)
        args.review_max_blocks = 2
        record = review.build(InputBundle(manifest), args, planner, [])
        mixed = next(c for c in record["contigs"] if c["target"] == "mixed")
        self.assertEqual(mixed["orientation"]["status"], "skipped_block_cap")
        self.assertEqual(mixed["alternatives"]["status"], "skipped_block_cap")
        self.assertEqual(mixed["canonical_alignment_blocks"], 5)
        self.assertEqual(len([b for b in record["blocks"] if b["target"] == "mixed"]), 5)
        args.alternative_plans = 21
        with self.assertRaisesRegex(ValueError, "0–20"):
            review.validate_arguments(args)

    def test_kbest_matches_exhaustive_additive_objective_and_boundary_rule(self):
        blocks = [fix.QueryBlock("c", ref, 0, end-start, start, end, strand, weight, 100*weight, 1)
            for start, end, ref, strand, weight in [(0, 130, "A", "+", 130), (80, 100, "B", "+", 20),
                (110, 240, "C", "+", 130), (240, 280, "C", "-", 40), (280, 400, "B", "+", 120)]]
        args = SimpleNamespace(breakpoint_penalty_bp=31)
        exhaustive = []
        for flags in itertools.product([False, True], repeat=4):
            cuts = tuple(i + 1 for i, yes in enumerate(flags) if yes)
            edges = [0, *cuts, 5]
            discordant = 0
            for start, end in zip(edges, edges[1:]):
                weights = defaultdict(float)
                for block in blocks[start:end]:
                    weights[block.ref, block.orientation] += block.weighted_bp
                discordant += sum(weights.values()) - max(weights.values())
            exhaustive.append((discordant + 31 * len(cuts), cuts))
        expected = sorted(exhaustive, key=lambda p: (p[0], len(p[1]), p[1]))
        self.assertEqual(review.kbest_partitions(blocks, args, True, 16), expected)
        self.assertEqual(review.kbest_partitions(blocks, args, True, 3), expected[:3])
        groups = [fix.summarize_group(blocks[:2], True), fix.summarize_group(blocks[2:], True)]
        self.assertEqual(fix.boundary_between(groups[0].blocks[-1], groups[1].blocks[0]), 105)
        self.assertNotEqual(fix.boundary_between(groups[0].summary, groups[1].summary), 105)

    def test_alternative_plans_are_unique_ranked_and_exclude_existing_plan(self):
        manifest = make_inspection_bundle(self.root)
        bundle = InputBundle(manifest)
        args, planner = settings(manifest)
        planner.breakpoint_penalty_bp = 25
        planner.min_piece_aligned_bp = 1
        planner.min_piece_query_frac = 0
        planner.alternative_plans = 5
        blocks = [fix.QueryBlock("mixed", frame, start, end, start, end, "+", end-start, 100*(end-start), 1)
                  for start, end, frame in [(0, 300, "A::chr"), (300, 400, "B::chr"),
                                             (400, 600, "A::chr"), (600, 1000, "C::chr")]]
        trace = review.partition_trace(blocks, 1000, planner, [])
        events, status = review.alternative_events(bundle, "mixed", blocks, [], trace, [300, 400, 600], planner)
        self.assertTrue(events)
        self.assertLessEqual(len(events), 5)
        self.assertNotIn([300, 400, 600], [e["positions"] for e in events])
        self.assertEqual(len(events), len({tuple(e["positions"]) for e in events}))
        scores = [e["objective"]["total"] for e in events]
        self.assertEqual(scores, sorted(scores))
        for i, event in enumerate(events, 1):
            objective = event["objective"]
            self.assertEqual(objective["rank"], i)
            self.assertEqual(objective["total"], objective["discordant_bp"] + objective["penalty_bp"])
            self.assertEqual(objective["delta_to_best_partition"], objective["total"] - trace["best_partition"]["total"])
        review.coordinate_evidence(events, bundle, args)
        self.assertEqual(status["paths_per_prefix"], 24)
        for event in events:
            self.assertEqual(event["positions"], [r["position"] for r in event["coordinate_evidence"]])
            self.assertTrue(all(r["reads"]["continuity_status"] == "unavailable_no_read_evidence"
                                and r["graph"]["status"] == "unavailable_no_graph" for r in event["coordinate_evidence"]))

    def test_orientation_event_cap_and_unavailable_graph_projection(self):
        manifest = make_inspection_bundle(self.root, graph=True)
        bundle = InputBundle(manifest)
        args, _ = settings(manifest)
        args.inspect_min_orientation_bp = 1
        raw = []
        for index in range(61):
            raw.append({"block_id": str(index), "backbone": True, "frame": "A::chr", "start0": index * 10,
                "end0": (index + 1) * 10, "orientation": "-" if index % 2 else "+", "loss_stages": [],
                "correction_eligible": False, "row_count": 1})
        events, ledger = review.orientation_events(bundle, "mixed", raw,
            {"include_orientation": False, "best_partition": None, "after_terminal_merge": []}, args)
        self.assertEqual(len(events), 25)
        self.assertEqual(ledger["omitted_by_event_cap"], 5)
        # A named graph segment with a different length has no reliable native projection.
        graph = self.root / "graph.gfa"
        graph.write_text("S\tmixed\tACGT\n")
        data = json.loads(manifest.read_text())
        data["evidence"][-1]["sha256"] = sha256_file(graph)
        manifest.write_text(json.dumps(data))
        review.coordinate_evidence(events[:1], InputBundle(manifest), args)
        self.assertTrue(all(c["graph"]["status"] == "unavailable_exact_projection" for c in events[0]["coordinate_evidence"]))

    def test_sensitive_mode_does_not_attribute_automatic_losses_to_diagnostic_dp(self):
        manifest = make_inspection_bundle(self.root)
        args, planner = settings(manifest)
        args.mode = planner.mode = "sensitive"
        automatic = [{"target": e.target, "action": "cut", "position": int(e.fields["slice_end"]), "reason": e.reason}
                     for e in evaluate.build_fix_events(planner) if e.action == "split_piece" and int(e.fields["slice_end"]) < 1000]
        record = review.build(InputBundle(manifest), args, planner, automatic)
        event = next(e for e in record["events"] if e["kind"] == "orientation_discrepancy" and e["positions"] == [700, 850])
        self.assertIn("represented_at_both_automatic_boundaries", event["loss_stages"])
        self.assertFalse(set(event["loss_stages"]) & {"smoothed_by_objective", "merged_with_unsupported_terminal_group", "fails_group_support"})
        for candidate in record["events"]:
            if candidate["objective"] is not None:
                self.assertFalse(candidate["objective"]["objective_comparable_to_automatic"])
                self.assertIn("Sensitive automatic correction bypasses", " ".join(candidate["reasons"]))

    def test_workflow_inspection_acknowledgment_explicit_edit_and_stage_binding(self):
        manifest = make_inspection_bundle(self.root)
        base = ["scan", "--manifest", str(manifest), "--mode", "comprehensive", "--min-segment-bp", "100",
                "--read-window-bp", "50", "--read-min-anchor-bp", "10", "--continuity-anchor-bp", "10"]
        workflow.main([*base, "--output-dir", str(self.root / "default")])
        scan = self.root / "review"
        workflow.main([*base, "--output-dir", str(scan), "--inspect-orientation", "--inspect-min-orientation-bp", "50", "--alternative-plans", "5"])
        baseline = workflow.read_decisions(self.root / "default" / "decisions.tsv")
        decisions = workflow.read_decisions(scan / "decisions.tsv")
        self.assertEqual(baseline, [d for d in decisions if not d["event_id"].startswith("inspection:")])
        record = json.loads((scan / "review_proposals.json").read_text())
        event = next(e for e in record["events"] if e["kind"] == "orientation_discrepancy" and e["start0"] == 300)
        for edge in event["figures"]:
            self.assertTrue((scan / (edge["path"] + ".svg")).is_file())
        for decision in decisions:
            decision.update(decision="accept" if decision["action"] == "inspect" else "reject", reviewer="unit fixture")
        reviewed = self.root / "decisions.tsv"
        write_tsv(reviewed, decisions, workflow.DECISION_COLUMNS)
        workflow.main(["apply", "--scan-dir", str(scan), "--decisions", str(reviewed), "--output-dir", str(self.root / "ack")])
        original = {n: s for n, _, s in iter_fasta_records(self.root / "assembly.fa")}
        self.assertEqual(original, {n: s for n, _, s in iter_fasta_records(self.root / "ack" / "reviewed.fa")})
        edit = self.root / "edit"
        workflow.main(["edit", "--scan-dir", str(scan), "--event-id", event["event_id"], "--action", "reverse",
                       "--start", "300", "--end", "380", "--notes", "Independent coordinate fixture", "--output-dir", str(edit)])
        edits = repairs.read_edits(edit / "reviewed_edits.tsv")
        edits[0].update(decision="accept", reviewer="unit fixture")
        write_tsv(edit / "accepted.tsv", edits, repairs.EDIT_COLUMNS)
        workflow.main(["apply", "--scan-dir", str(scan), "--decisions", str(reviewed), "--reviewed-edits", str(edit / "accepted.tsv"),
                       "--output-dir", str(self.root / "applied")])
        changed = {n: s for n, _, s in iter_fasta_records(self.root / "applied" / "reviewed.fa")}
        expected = original["mixed"][:300] + reverse_complement(original["mixed"][300:380]) + original["mixed"][380:]
        self.assertEqual(changed["mixed"], expected)
        manual.main(["apply", "--assembly-fasta", str(self.root / "assembly.fa"), "--recipe", str(self.root / "applied" / "recipe.json"),
                     "-o", str(self.root / "replay.fa")])
        self.assertEqual((self.root / "replay.fa").read_bytes(), (self.root / "applied" / "reviewed.fa").read_bytes())
        (scan / "review_proposals.json").write_text("{}")
        with self.assertRaisesRegex(ValueError, "Stale review proposals"):
            workflow.load_scan(scan)


if __name__ == "__main__":
    unittest.main()
