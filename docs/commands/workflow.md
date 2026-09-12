---
title: chromo workflow
---

# `chromo workflow`

This is the analysis handoff contract for researchers and their assistants. It
is separate from the repository's coding-agent instructions. An assistant can
collect evidence, generate proposals and explain them; a researcher records the
biological decisions, or a separately labeled scripted experiment supplies its
own decisions. Saved recipes replay without an agent or browser session.

```bash
chromo workflow scan --manifest inputs.json --output-dir review \
  --mode conservative --read-evidence-id reads_hifi
# Open review/index.html; inspect read panels and the linked alignment/graph
# dashboard; choose accept, reject, or refine; enter a reviewer label and
# download decisions.tsv. Preserve the scan.json and its original inputs.
chromo workflow apply --scan-dir review --decisions reviewed-decisions.tsv \
  --output-dir applied
chromo manual apply --assembly-fasta assembly.fa --recipe applied/recipe.json \
  -o replayed.fa
chromo workflow align --manifest inputs.json --assembly-fasta applied/reviewed.fa \
  --output-dir realigned --stage reviewed --threads 4
chromo workflow validate --manifest realigned/manifest.json \
  --apply-audit applied/apply.json --output validation.json
```

`scan` writes an editable decision TSV, immutable proposal record, assignment
and candidate reports, dot plots, read figures, and two linked dashboards. All
decisions begin `pending`. Each event ID depends on the assembly digest and its
native event geometry. Uncovered internal backbone-alignment intervals are
`inspect` events; acknowledging one makes no sequence change. Reference-only
evidence cannot establish biological error truth.

`apply` requires every proposal exactly once, an explicit decision and a
reviewer label. `defer` records that an event remains unresolved and preserves
that boundary. `accept` keeps the proposed cut; `reject` leaves that boundary
intact; `refine` uses the reviewed native position. Other proposal identity
fields cannot be changed. Separate interval edits use their own explicit
records, described below. Unreviewed contigs remain intact. Accepted cuts
partition each source contig exactly once; no fragment disappears when another
cut is rejected. The recipe, TSV, input closure, software version, FASTA, AGP,
component ledger and accounting report form the audit trail. A new output
directory is required for each scan/apply revision.

Use the explicit interval-repair path below for additional cuts and reversals
bound to the scan. Use `manual` recipes for arbitrary reorderings, renaming and
intentional removals. New dashboard recipes bind the original
assembly digest. Supervised workflow recipes also enforce complete source
coverage. Legacy recipes without digests remain readable but do not provide
that identity guarantee. Direct `fix --reviewed-plan` now rejects partially
accepted piece sets that would omit or duplicate source sequence; reject all
pieces to keep the contig, or use workflow decisions to reject individual cuts.

`align` invokes minimap2 with recorded version and argument arrays, `-c` and
`--secondary=no`, separately for each reference. It writes a fresh manifest
for the exact changed FASTA. Optional `--reads reads.fastq.gz --read-sample ID`
also remaps reads (`--read-preset map-hifi`, `map-ont`, or `map-pb`). Prior graph,
read and variant evidence is not carried into a changed coordinate system.
Prepare fresh graph mappings with `gafprep`/GraphAligner where needed. Mapping
cost is distinct from curation cost.

`validate` verifies changed-FASTA evidence identity, complete source accounting,
and exact AGP reconstruction. It reports biological validation as unresolved;
identity and replay checks are not an independent assessment of correctness.
Continue with fresh-manifest `sort --orient-to-reference`, reviewed scaffolding,
and another realignment after each FASTA change. Whole-contig orientation is a
separate explicit operation from cutting an internal event.

## Selective scaffold review

```bash
chromo eval scaffold --ordered-fasta sorted.ordered.fa \
  --assignments sorted.contig_assignments.tsv -o scaffold-review
# Edit accept, gap_bp, overlap_policy and notes in scaffold-review.scaffold_review.tsv
chromo scaffold --ordered-fasta sorted.ordered.fa \
  --assignments sorted.contig_assignments.tsv \
  --reviewed-plan scaffold-review.scaffold_review.tsv \
  --require-reviewed-joins -o scaffolded
```

With `--require-reviewed-joins`, rejected or absent joins separate scaffold
records. Without it, legacy gap-override behavior is preserved. Each accepted
junction can choose `zero-gap`, `warn`, `trim-reference`, or `trim-sequence`.
`trim-reference` is an explicit deletion policy and is not sequence validation.
Graph-confirmed trimming additionally requires terminal sequence identity;
reference coordinate overlap and a graph edge alone do not justify deletion.
Backbone, subgenome and copy groups cannot be mixed in one scaffold.

Record assistant model/harness, prompts, budgets and human interventions outside
the input-only curation bundle if evaluating assisted workflows. Do not expose
held-out truth to the assistant or reviewer. A scripted diagnostic establishes
replay feasibility; it does not measure human effort, accuracy or biological
superiority.


## Reproducible development diagnostics

From a repository checkout, use the compact [public development fixtures]({{ "/benchmarks/" | relative_url }})
and a fresh installed candidate Python (without `PYTHONPATH=src`):

```bash
python scripts/check_public_development.py \
  --development-root benchmarks/assembly-errors/fixtures/development \
  --output-dir results/development-check \
  --software-archive dist/chromosort-0.4.0rc2.tar.gz
python scripts/check_public_workflow.py \
  --manifest benchmarks/assembly-errors/fixtures/development/D001/manifest.json \
  --read-evidence-id reads_sim01 --mode sensitive \
  --output-dir results/workflow-check --realign
```

The first checks source accounting and AGP reconstruction for every development
manifest. The second explicitly labels its accept-all/reject-all decisions as
scripted diagnostics, checks byte-identical recipe replay, and optionally realigns
and validates the changed FASTA. Neither reads benchmark truth or measures human
judgment. The compact fixtures contain read PAF; larger generated bundles may contain
PAF and BAM representations of the same mapping. Select one explicit evidence
ID for read panels; correlated representations are not independent support.
Use the exact software archive being evaluated, including its build directory
when checking a newer repository snapshot.


## Continuity conflicts and gap inspection

Manifest scans retain the original cut proposals but explain when read continuity
would defer their automatic application. Decisions remain pending. The continuity
source, filters and counts are in the event fields; `--read-evidence-id` chooses
the source for both this assessment and the read panels.

`--inspect-min-gap-bp 1000` independently controls review-only uncovered intervals
between flanks at least `--min-segment-bp` long. It does not change correction
smoothing, `--min-piece-bp`, or the cut planner. Each inspect item links two read
panels, one at each interval edge. Its table position remains the interval center,
which is a locator, not a cut. Accepting an inspect item changes no sequence.
Use `--contigs` to limit a large review queue; raising the inspection floor reduces
queue size without changing the automatic correction algorithm.

See [the policy and its limitations]({{ '/review-refinement/' | relative_url }}).

## Explicit interval repairs

`workflow edit` prepares a **pending** request; it does not detect or correct an
inversion. Link the requested edit to an existing event on the same contig.
For example, after reviewing a gap or a contig-level proposal:

```bash
# EVENT_ID is copied from review/decisions.tsv. Coordinates here are examples.
chromo workflow edit --scan-dir review --event-id "$EVENT_ID" \
  --action reverse --start 20000 --end 70000 \
  --notes 'Explicitly proposed interval; review junction evidence before applying' \
  --output-dir interval-review
# Inspect interval-review/edit-review.html and both native-boundary panels.
# In interval-review/reviewed_edits.tsv, set decision and reviewer explicitly.
# Resolve every scan event in reviewed-decisions.tsv, using defer when unresolved.
chromo workflow apply --scan-dir review --decisions reviewed-decisions.tsv \
  --reviewed-edits interval-review/reviewed_edits.tsv \
  --plan-only --output-dir repair-preview
# Inspect repair-preview/recipe.json, plan.json, decision_outcomes.tsv and figures.
chromo workflow apply --scan-dir review --decisions reviewed-decisions.tsv \
  --reviewed-edits interval-review/reviewed_edits.tsv --output-dir repaired
chromo manual apply --assembly-fasta assembly.fa --recipe repaired/recipe.json \
  -o replayed.fa
chromo workflow align --manifest inputs.json --assembly-fasta repaired/reviewed.fa \
  --output-dir repaired-alignment --stage explicitly-reviewed
chromo workflow validate --manifest repaired-alignment/manifest.json \
  --apply-audit repaired/apply.json --output repaired-validation.json
```

Native intervals are **0-based, half-open**: `reverse --start s --end e` reverse
complements exactly `assembly[s:e]` and retains its prefix and suffix. An extra
cut uses `--action cut --start p --end p`, where `p` is an internal boundary.
These operations never infer a missing alignment block, copy number or error
origin. Reference disagreement and absent spanning reads do not establish that
a reversal is justified.

The edit TSV records `schema`, `edit_id`, `event_id`, `target`, `action`, `start`,
`end`, `decision`, `reviewer`, `notes`, `assembly_sha256` and `scan_sha256`. Its
identity binds the exact geometry and original scan. Change only the decision,
reviewer and rationale; generate a new record to change geometry or action.
`accept` applies an explicit edit; `reject` and `defer` record distinct reasons
to leave that intervention unapplied. `pending` blocks apply, including a
plan-only compilation. The initial `workflow edit` packet is available for
evidence inspection while a decision remains pending. Multiple generated edit
rows can be combined under one header, with each edit ID appearing once.

The preview writes a complete recipe and both-edge read panels without writing
a FASTA. Actual apply uses a new output directory. `plan.json` records requested
operations; `apply.json` and `decision_outcomes.tsv` record completed operations
after sequence/AGP checks pass. The original source coordinates remain in the
recipe, component ledger and audit even when output coordinates change. A
terminal reverse interval is supported; its terminal panel has no internal
breakpoint support count.

### Composition rules

- Accepted reverse intervals must be disjoint. Adjacent intervals are separate
  reversals, not one combined reversal.
- An accepted cut may coincide with a reverse endpoint, but cannot lie strictly
  inside the reversed interval. Resolve that composition explicitly before apply.
- Duplicate accepted cuts and an extra cut contradicting a rejected/deferred
  scan cut at the same boundary fail. Resolve an existing proposed cut in its
  scan decision row.
- Within each cut-delimited output record, source slices remain in native order.
  Only the explicitly selected interval is reverse-complemented; all other
  slices keep positive orientation. Joining these slices inserts no bases.
- Every input base occurs exactly once. Unselected copies and contigs remain
  present. No edit can implicitly remove sequence or normalize another variant.

An automatic continuity conflict still defers the whole contig plan. Selective
execution is available through explicit review, where each boundary decision
and interval reversal is recorded. An accepted `inspect` row alone remains an
acknowledgement and makes no sequence change. Scripted repairs using known
development coordinates demonstrate repair/replay capability under supplied
decisions; they do not establish autonomous discovery, biological classification
or human/agent review performance.
