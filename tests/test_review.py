import csv
import tempfile
import unittest
from pathlib import Path

from chromosort.review import (
    ReviewEvent,
    accepted_events,
    read_review_events,
    review_columns,
    write_review_events,
)


class ReviewEventTests(unittest.TestCase):
    def test_review_columns_include_core_explicit_and_discovered_fields(self):
        event = ReviewEvent(
            event_id="fix:ctg1:1",
            task="fix",
            action="split_piece",
            target="ctg1",
            fields={"slice_start": 1, "slice_end": 100},
        )

        columns = review_columns([event], extra_columns=["dominant_ref", "slice_start"])

        self.assertEqual(columns[:10], [
            "schema",
            "event_id",
            "task",
            "action",
            "target",
            "accept",
            "status",
            "confidence",
            "reason",
            "notes",
        ])
        self.assertIn("dominant_ref", columns)
        self.assertIn("slice_start", columns)
        self.assertIn("slice_end", columns)

    def test_writes_and_reads_review_events(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "events.tsv"
            events = [
                ReviewEvent(
                    event_id="fix:ctg1:1",
                    task="fix",
                    action="split_piece",
                    target="ctg1",
                    accept=True,
                    confidence="high",
                    reason="planner_and_reads",
                    notes="keep",
                    fields={"slice_start": 1, "slice_end": 100, "dominant_ref": "chr1"},
                ),
                ReviewEvent(
                    event_id="fix:ctg1:2",
                    task="fix",
                    action="split_piece",
                    target="ctg1",
                    fields={"slice_start": 101, "slice_end": 200, "dominant_ref": "chr2"},
                ),
            ]

            write_review_events(
                path,
                events,
                extra_columns=["dominant_ref", "slice_start", "slice_end"],
            )
            loaded = read_review_events(path, expected_task="fix")

        self.assertEqual(len(loaded), 2)
        self.assertTrue(loaded[0].accept)
        self.assertFalse(loaded[1].accept)
        self.assertEqual(loaded[0].fields["dominant_ref"], "chr1")
        self.assertEqual([event.event_id for event in accepted_events(loaded)], ["fix:ctg1:1"])

    def test_rejects_invalid_accept_value(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "events.tsv"
            with open(path, "w", newline="") as out:
                writer = csv.writer(out, delimiter="\t")
                writer.writerow([
                    "schema",
                    "event_id",
                    "task",
                    "action",
                    "target",
                    "accept",
                    "status",
                    "confidence",
                    "reason",
                    "notes",
                ])
                writer.writerow([
                    "chromosort-review-event-v1",
                    "fix:ctg1",
                    "fix",
                    "split_piece",
                    "ctg1",
                    "maybe",
                    "candidate",
                    ".",
                    ".",
                    "",
                ])

            with self.assertRaisesRegex(ValueError, "Invalid accept value"):
                read_review_events(path)

    def test_rejects_duplicate_event_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "events.tsv"
            events = [
                ReviewEvent("fix:ctg1", "fix", "split_piece", "ctg1"),
                ReviewEvent("fix:ctg1", "fix", "split_piece", "ctg1"),
            ]

            with self.assertRaisesRegex(ValueError, "Duplicate review event_id"):
                write_review_events(path, events)

    def test_rejects_wrong_expected_task(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "events.tsv"
            write_review_events(
                path,
                [ReviewEvent("scaffold:j1", "scaffold", "join", "ctgA|ctgB")],
            )

            with self.assertRaisesRegex(ValueError, "expected 'fix'"):
                read_review_events(path, expected_task="fix")


if __name__ == "__main__":
    unittest.main()
