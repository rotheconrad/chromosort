---
title: Production Upgrade Roadmap
description: Planned eval and manual review modes for task-specific human-guided assembly curation.
---

# Production Upgrade Roadmap

This roadmap captures the next major production upgrade direction for
ChromoSort. It is a planning document, not a description of currently available
commands.

## Goal

The next upgrade will add a paired review layer around the sequence-changing
commands that most often need human judgment:

- `fix`
- `scaffold`
- `gapfill`

The existing `fix`, `scaffold`, and `gapfill` behavior should remain available
as the default algorithmic path. Reviewed paths will add optional inputs that
let users accept, reject, or add specific decisions when biological outliers
need more control than a fully automatic run can provide.

## Review Layer Design

`chromo eval` and `chromo manual` should become counterpart interfaces over the
same underlying review-event model:

| Layer | Audience | Primary output |
| --- | --- | --- |
| `chromo eval` | command-line and spreadsheet-first users | task-specific TSV tables |
| `chromo manual` | visual review and GUI-first users | task-specific browser dashboards |

The intent is for both layers to expose the same evidence, candidate decisions,
and row identifiers. `eval` should stay table-only. `manual` should remain
exploratory, but each mode should also zoom directly to the candidate events
that need review.

## Planned Eval Modes

### `chromo eval fix`

Prepare an editable table for candidate contig-fix decisions. Rows should
describe proposed split pieces or candidate breakpoint groups, with supporting
evidence from:

- whole-genome reference-to-assembly alignments,
- optional GFA node context,
- optional long-read backmapping to the assembly.

The reviewed table should become an optional input to `chromo fix`, allowing
the command to apply explicit accepted slices while preserving its current
planner-driven behavior when no reviewed table is supplied.

### `chromo eval scaffold`

Prepare an editable table for adjacent-contig scaffold decisions. Rows should
describe candidate junctions, gaps, overlaps, order/orientation concerns, and
graph or read support. Evidence should include:

- `chromo sort` assignment intervals,
- inferred reference-space gaps and overlaps,
- optional GFA direct edges or short paths,
- optional long-read bridges between contig ends.

The reviewed table should become an optional input to `chromo scaffold`, so
users can pin down outlier junctions while leaving ordinary scaffold decisions
algorithmic.

### `chromo eval gapfill`

Prepare an editable table for graph-supported fill decisions. This should be a
table-only counterpart to the existing gapfill planning behavior, expanded as
needed to include long-read backmapping. Evidence should include:

- candidate GFA paths,
- GAF read-path support when available,
- Hi-C-like graph-node contacts when available,
- reference-placement PAF support,
- long-read evidence mapped back to the assembly.

The reviewed table should become an optional input to `chromo gapfill`, with
the command revalidating accepted rows before applying sequence.

## Planned Manual Modes

`chromo manual` should evolve from one general review dashboard into
task-specific review modes:

- `chromo manual fix`
- `chromo manual scaffold`
- `chromo manual gapfill`

Each mode should keep the useful "browse around and explore" feel, but it
should open around a focused queue of candidate events. Selecting a row should
zoom to the relevant contig, breakpoint, junction, or graph path while still
allowing nearby alignments, graph neighborhoods, and contigs to be inspected.

## Shared Review Events

A shared internal review-event model should keep `eval` and `manual` from
drifting apart. A review event should carry:

- stable identifiers for the source FASTA records and candidate action,
- task type (`fix`, `scaffold`, or `gapfill`),
- proposed action and default acceptance state,
- evidence summaries from alignments, graph inputs, and long-read support,
- validation fields needed by the execution command,
- human-editable accept/reject and note fields.

The execution commands should validate reviewed tables against the current
FASTA, assignments, graph paths, and settings before changing sequence.
Rejected, deleted, stale, or invalid rows should fail safely or fall back to the
current conservative behavior, depending on the command and mode.

## Long-Read Backmapping Evidence

Long reads mapped back to the assembly should become a shared evidence source
for the production upgrade. The first implementation should prefer a lightweight
PAF-based reader to avoid adding heavy required dependencies. BAM or CRAM
support can follow later as an optional capability.

The shared evidence layer should summarize:

- reads spanning candidate breakpoints,
- split or clipped read clusters near breakpoints,
- contig-end read bridges with orientation and estimated gap or overlap,
- coverage anomalies around candidate cuts or joins,
- concordance or conflict with GFA links and graph paths.

## Implementation Phases

1. Add shared long-read evidence parsing and summary helpers.
2. Add the shared review-event data model and TSV serialization.
3. Implement `chromo eval fix` and `chromo fix` reviewed-table application.
4. Implement `chromo eval scaffold` and `chromo scaffold` reviewed-table
   application.
5. Align `chromo eval gapfill` with current gapfill plan semantics and add any
   missing reviewed-table compatibility.
6. Refactor `chromo manual` into task-specific modes over the same review-event
   model.
7. Expand docs and tests around mixed algorithmic-plus-reviewed workflows.

## Guardrails

- Existing default behavior remains the default.
- Reviewed execution paths must be explicit.
- Table schemas should be stable and spreadsheet-friendly.
- Sequence-changing commands must revalidate reviewed rows before applying
  them.
- Ambiguous or stale evidence should remain reviewable rather than guessed.
