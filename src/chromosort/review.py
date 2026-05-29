"""Shared review-event table helpers."""

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

from .paths import ensure_parent_dir


SCHEMA = "chromosort-review-event-v1"
TASKS = {"fix", "scaffold", "gapfill"}
CORE_COLUMNS = [
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
]
TRUE_VALUES = {"1", "true", "t", "yes", "y", "accept", "accepted"}
FALSE_VALUES = {"0", "false", "f", "no", "n", "reject", "rejected", "", "."}


@dataclass(frozen=True)
class ReviewEvent:
    """One editable candidate decision shared by eval, manual, and executors."""

    event_id: str
    task: str
    action: str
    target: str
    accept: bool = False
    status: str = "candidate"
    confidence: str = "."
    reason: str = "."
    notes: str = ""
    fields: Dict[str, object] = field(default_factory=dict)
    schema: str = SCHEMA

    def __post_init__(self):
        if not self.event_id:
            raise ValueError("review event_id is required")
        if self.task not in TASKS:
            raise ValueError(f"unsupported review task: {self.task!r}")
        if not self.action:
            raise ValueError("review action is required")
        if not self.target:
            raise ValueError("review target is required")
        if self.schema != SCHEMA:
            raise ValueError(f"unsupported review schema: {self.schema!r}")

    def to_row(self, columns: Sequence[str]):
        row = {}
        for column in columns:
            if column == "schema":
                value = self.schema
            elif column == "event_id":
                value = self.event_id
            elif column == "task":
                value = self.task
            elif column == "action":
                value = self.action
            elif column == "target":
                value = self.target
            elif column == "accept":
                value = "yes" if self.accept else "no"
            elif column == "status":
                value = self.status
            elif column == "confidence":
                value = self.confidence
            elif column == "reason":
                value = self.reason
            elif column == "notes":
                value = self.notes
            else:
                value = self.fields.get(column, ".")
            row[column] = format_field(value)
        return row

    @classmethod
    def from_row(cls, row, line_number=None):
        location = f" on line {line_number}" if line_number is not None else ""
        schema = row.get("schema", "")
        if schema != SCHEMA:
            raise ValueError(f"Unsupported review schema{location}: {schema!r}.")
        missing = [column for column in CORE_COLUMNS if column not in row]
        if missing:
            raise ValueError(
                f"Review table is missing required column(s){location}: "
                f"{', '.join(missing)}."
            )
        fields = {
            key: value
            for key, value in row.items()
            if key not in CORE_COLUMNS and value not in {None, ""}
        }
        return cls(
            event_id=row["event_id"],
            task=row["task"],
            action=row["action"],
            target=row["target"],
            accept=parse_accept(row.get("accept", "no"), line_number=line_number),
            status=row.get("status") or "candidate",
            confidence=row.get("confidence") or ".",
            reason=row.get("reason") or ".",
            notes=row.get("notes") or "",
            fields=fields,
            schema=schema,
        )


def parse_accept(value, line_number=None):
    text = str(value).strip().lower()
    if text in TRUE_VALUES:
        return True
    if text in FALSE_VALUES:
        return False
    location = f" on line {line_number}" if line_number is not None else ""
    raise ValueError(f"Invalid accept value{location}: {value!r}.")


def format_field(value):
    if value is None:
        return "."
    if isinstance(value, bool):
        return "yes" if value else "no"
    return str(value)


def review_columns(events: Iterable[ReviewEvent], extra_columns: Optional[Sequence[str]] = None):
    """Return core columns plus task-specific columns used by events."""

    columns = list(CORE_COLUMNS)
    seen = set(columns)
    for column in extra_columns or []:
        if column not in seen:
            columns.append(column)
            seen.add(column)
    discovered = []
    for event in events:
        for column in event.fields:
            if column not in seen:
                discovered.append(column)
                seen.add(column)
    columns.extend(sorted(discovered))
    return columns


def ensure_unique_event_ids(events: Sequence[ReviewEvent]):
    seen = {}
    for index, event in enumerate(events, start=1):
        if event.event_id in seen:
            raise ValueError(
                f"Duplicate review event_id {event.event_id!r} "
                f"at rows {seen[event.event_id]} and {index}."
            )
        seen[event.event_id] = index


def write_review_events(path, events, extra_columns: Optional[Sequence[str]] = None):
    """Write review events as a spreadsheet-friendly TSV."""

    events = list(events)
    ensure_unique_event_ids(events)
    columns = review_columns(events, extra_columns)
    ensure_parent_dir(path)
    with open(path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for event in events:
            writer.writerow(event.to_row(columns))


def read_review_events(path, expected_task=None):
    """Read a review-event TSV and return validated events."""

    events: List[ReviewEvent] = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Review table is empty: {path}")
        missing = [column for column in CORE_COLUMNS if column not in reader.fieldnames]
        if missing:
            raise ValueError(
                f"Review table {path} is missing required column(s): "
                f"{', '.join(missing)}."
            )
        for line_number, row in enumerate(reader, start=2):
            event = ReviewEvent.from_row(row, line_number=line_number)
            if expected_task is not None and event.task != expected_task:
                raise ValueError(
                    f"Review event {event.event_id!r} on line {line_number} "
                    f"has task {event.task!r}; expected {expected_task!r}."
                )
            events.append(event)
    ensure_unique_event_ids(events)
    return events


def accepted_events(events: Iterable[ReviewEvent]):
    return [event for event in events if event.accept]


def events_by_id(events: Iterable[ReviewEvent]):
    grouped = {}
    for event in events:
        if event.event_id in grouped:
            raise ValueError(f"Duplicate review event_id {event.event_id!r}.")
        grouped[event.event_id] = event
    return grouped


def event_table_path(prefix, suffix):
    return Path(str(prefix) + suffix)
