---
title: Production Upgrade Roadmap
description: Eval and manual review modes for task-specific human-guided assembly curation.
---

# Production Upgrade Roadmap

This roadmap tracks the production review upgrade for ChromoSort. Completed
rows describe behavior now available in the current development branch; planned
rows describe remaining follow-up work if new review needs appear.

## Goal

The upgrade adds a paired review layer around the sequence-changing
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

## Eval Modes

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

## Manual Modes

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

| Phase | Status | Notes |
| --- | --- | --- |
| 1. Add shared long-read evidence parsing and summary helpers. | Done | Added PAF-backed long-read evidence helpers for breakpoint support, contig-end bridges, and read-depth summaries. |
| 2. Add the shared review-event data model and TSV serialization. | Done | Added a shared review-event schema and TSV reader/writer for `eval`, `manual`, and reviewed execution paths. |
| 3. Implement `chromo eval fix` and `chromo fix` reviewed-table application. | Done | Added `chromo eval fix` table generation and `chromo fix --reviewed-plan` application for accepted `split_piece` rows. |
| 4. Implement `chromo eval scaffold` and `chromo scaffold` reviewed-table application. | Done | Added table-only scaffold junction evaluation with GFA and long-read bridge context, plus accepted gap overrides through `chromo scaffold --reviewed-plan`. |
| 5. Align `chromo eval gapfill` with current gapfill plan semantics and add any missing reviewed-table compatibility. | Done | Added table-only gapfill evaluation from the existing graph path planner, long-read bridge context, and `chromo gapfill --reviewed-plan` support for the shared review-event table. |
| 6. Refactor `chromo manual` into task-specific modes over the same review-event model. | Done | Added `chromo manual fix`, `chromo manual scaffold`, and `chromo manual gapfill` modes with shared review-event table embedding and focused event queues inside the existing exploratory dashboard. |
| 7. Expand docs and tests around mixed algorithmic-plus-reviewed workflows. | Done | Updated command, input, status, changelog, README, and workflow docs and verified the full test suite across the review upgrade. |

## Guardrails

- Existing default behavior remains the default.
- Reviewed execution paths must be explicit.
- Table schemas should be stable and spreadsheet-friendly.
- Sequence-changing commands must revalidate reviewed rows before applying
  them.
- Ambiguous or stale evidence should remain reviewable rather than guessed.

## Next Chapter: GAF Evidence And Modular Manual Panels

The next review upgrade should fully surface long-read GAF graph alignments
across the same `eval` and `manual` task modes. GAF should be treated as a
complementary evidence stream rather than a replacement for GFA or long-read
PAF:

- GFA answers whether the graph topology permits a relationship.
- Long-read PAF answers whether reads support breakpoints or junctions in
  assembly-contig coordinate space.
- Long-read GAF answers whether reads traverse graph paths that support,
  disambiguate, or conflict with a proposed graph relationship.

The main product goal is communication. Users should be able to inspect four
independent evidence panels when the corresponding inputs are provided, while
the command-line `eval` tables expose compact summary columns for the same
evidence.

### Evidence Inputs

| Evidence | Role | Optional input |
| --- | --- | --- |
| Whole-genome alignment | Reference placement, dot plots, fix candidates, sort/scaffold coordinates | `--coords` or whole-genome `--paf` |
| Assembly graph | Graph nodes, direct links, short paths, graph complexity | `--gfa` |
| Long-read PAF to assembly | Breakpoint support, contig-end read bridges, read-space gap/overlap estimates | `--read-paf` |
| Long-read GAF to graph | Read traversal support for graph paths, branch support, path conflicts | `--gaf` |

Each evidence stream should remain optional where possible. The dashboard
should render only the panels backed by provided evidence, and `eval` should
write `.` values for unavailable evidence rather than forcing unnecessary
inputs.

### Eval Mode Expansion

`chromo eval fix`, `chromo eval scaffold`, and `chromo eval gapfill` should all
accept long-read GAF evidence.

For `eval fix`, GAF evidence should remain advisory. It should summarize graph
traversal context near candidate split nodes or breakpoint-associated graph
nodes, but PAF-to-assembly evidence should remain the primary coordinate-level
support for breakpoint decisions.

For `eval scaffold`, GAF evidence should be first-class support. Scaffold
junction rows should report whether reads traverse the direct GFA edge or the
short graph path connecting adjacent contigs, whether support favors an
alternate path, and whether graph traversal conflicts with the reference-space
order.

For `eval gapfill`, GAF is already central to branch resolution. The next work
should move GAF parsing and path-support helpers into a shared evidence module,
then reuse the same summaries in the review-event table and manual dashboard.

### Manual Dashboard Panels

Task-specific manual dashboards should expose a modular evidence layout for the
selected event, contig, junction, or candidate path.

| Panel | Shows | Present when |
| --- | --- | --- |
| Alignment panel | Dot plot, reference placement, local alignment rows, task target context | `--coords` or whole-genome `--paf` |
| GFA panel | Node status, direct links, short paths, graph complexity, neighbor context | `--gfa` |
| Long-read PAF panel | Spanning reads, split/edge reads, contig-end bridges, median read-space gaps | `--read-paf` |
| Long-read GAF panel | Read-path support, candidate-path support, branch/conflict summaries | `--gaf` |

The current manual dashboard already has the alignment view, optional GFA
context, and task-specific event queue. The next dashboard step is to split the
selected-event evidence into these explicit panels instead of relying only on
review-table columns and badges.

### Executor Policy

The sequence-changing commands should remain conservative:

- `chromo fix` may consume GAF-derived review-table rows only when the reviewed
  table explicitly accepts the change; automatic fixing should still rely on
  existing coordinate-based logic.
- `chromo scaffold` may use GAF as report-only support or as reviewed-table
  context for accepted gap decisions, but it should not reorder, orient, or
  trim based on GAF alone.
- `chromo gapfill` can continue using GAF as a branch-support signal, with
  reviewed rows revalidated against the current GFA path before sequence is
  applied.

### GAF Implementation Phases

| Phase | Status | Notes |
| --- | --- | --- |
| G1. Move GAF parsing and path-support helpers into a shared evidence module. | Done | Added `chromosort.gaf` with GAF path parsing, MAPQ-filtered reading, oriented subpath support counting, and gapfill-compatible path support helpers. |
| G2. Add shared GAF summary objects for graph traversal evidence. | Done | Added reusable traversal summary objects with selected-path support, best alternate support, support status, selected reads, and per-path read support details. |
| G3. Add `--gaf` support to `chromo eval scaffold`. | Done | Scaffold review rows now report GAF support for direct or short GFA paths between adjacent contigs, including selected-path support, best alternate support, status, and supporting reads. |
| G4. Add advisory `--gaf` support to `chromo eval fix`. | Done | Fix review rows now summarize GAF reads traversing the candidate contig or resolved graph node without changing automatic fix decisions. |
| G5. Refactor `chromo eval gapfill` to use the shared GAF evidence layer. | Done | Gapfill planning now uses the shared GAF traversal summary layer and exposes support status plus selected supporting reads while preserving existing branch-resolution behavior. |
| G6. Add modular evidence panels to `chromo manual fix/scaffold/gapfill`. | Done | Manual task dashboards now accept optional `--read-paf` and `--gaf` paths and render alignment, GFA, long-read PAF, and long-read GAF evidence panels from provided inputs and selected review-event fields. |
| G7. Add docs, fixtures, and regression tests for mixed GFA/PAF/GAF review. | Done | Added a mixed-evidence scaffold fixture and regression where long-read PAF supports the junction while GAF more strongly supports an alternate GFA branch, plus workflow/eval documentation for the resulting review fields. |

## Next Chapter: Architecture And Documentation Consistency

The architecture documentation now needs a publication-style refresh that
describes not only what each algorithm does, but exactly where each algorithm,
decision rule, data model, and evidence stream is used. The guiding question
for reviewers should be:

> Which subcommand, mode, or parameter activates this algorithm or data model?

The update should keep the existing methods-oriented architecture style:
precise claims, traceable command boundaries, conservative descriptions of
failure modes, and explicit links between code modules, CLI parameters, and
user-facing artifacts. The README and user docs should then be checked against
that architecture map so the project tells one coherent story.

### Architecture Coverage Targets

| Coverage target | Documentation outcome |
| --- | --- |
| Algorithm activation map | A table or section in `docs/architecture.md` that maps algorithms and decision rules to commands, modes, parameters, outputs, and validation tests. |
| Data-model usage map | Clear descriptions of where `Segment`, match metrics, graph records, review events, manual recipes, scaffold gaps, and fill plans enter and leave the workflow. |
| Evidence-source map | Consistent descriptions of whole-genome alignments, GFA, long-read PAF, GAF, Hi-C-like pairs, and reference-placement PAF across architecture, input, eval, manual, scaffold, and gapfill docs. |
| Executor policy map | Explicit notes on which evidence streams are report-only, which can alter reviewed rows, and which can directly change sequence under guarded parameters. |
| Cross-doc consistency | README, status, workflows, command docs, input/output docs, and changelog should agree on current commands, optional inputs, reviewed-table paths, and version metadata. |

### Documentation Refresh Phases

| Phase | Status | Notes |
| --- | --- | --- |
| A1. Add this architecture-docs roadmap chapter. | Done | Track the documentation refresh separately from the completed review and GAF implementation chapters. |
| A2. Refresh `docs/architecture.md` with algorithm and data-model usage maps. | Done | Added command/parameter activation maps for sorting, fixing, eval modes, manual modes, scaffolding, gapfill branch resolution, review events, long-read evidence, GAF summaries, and evidence authority. |
| A3. Synchronize README and command docs with the architecture map. | Done | Updated README and command pages so eval/manual/scaffold/gapfill/fix descriptions align with algorithm activation, evidence authority, reviewed plans, and modular evidence panels. |
| A4. Synchronize workflow, input, output, status, and troubleshooting docs. | Done | Updated evidence contracts, eval outputs, status/architecture links, workflow guidance, index steps, and troubleshooting text to match the architecture map. |
| A5. Run documentation and test consistency checks. | Planned | Use text searches, metadata checks, and the test suite to catch stale version strings, missing command references, and broken assumptions. |
| A6. Bump the patch version and release metadata by 0.0.1. | Planned | Update package metadata, citation, Pixi, conda recipe, README/status/changelog, then tag and push the release. |
