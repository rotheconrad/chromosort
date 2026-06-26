---
title: Spreadsheet Review Tables
description: Use chromo eval review-event TSVs to accept, reject, audit, and apply reviewed fix, scaffold, and gapfill decisions.
---

# Spreadsheet Review Tables

Use this guide when ChromoSort has proposed candidate decisions and you want a
spreadsheet-first review step before sequence changes.

The key question is:

> Which rows are accepted, which executor can apply them, and how does
> ChromoSort prevent stale rows from changing a new FASTA?

## The Core Idea

`chromo eval` writes editable TSV review tables. These tables are decision
queues, not FASTA outputs. A reviewer can accept, reject, annotate, delete, or
add rows, then pass the reviewed table to the matching sequence-changing
executor.

All shared eval tables use schema `chromosort-review-event-v1`. The required
core columns are:

| Column | Meaning |
| --- | --- |
| `schema` | Must be `chromosort-review-event-v1`. |
| `event_id` | Unique row identifier used for validation and audit. |
| `task` | One of `fix`, `scaffold`, or `gapfill`. |
| `action` | Candidate action, such as `split_piece`, `scaffold_gap`, or `fill_path`. |
| `target` | Human-readable row target. |
| `accept` | The review decision used by the executor. |
| `status` | Candidate or planner status. |
| `confidence` | Planner or reviewer confidence field when available. |
| `reason` | Why the row exists or why it was accepted or rejected. |
| `notes` | Free-form reviewer notes. |

Task-specific columns hold source contigs, slices, gap lengths, graph paths,
support metrics, risk flags, and provenance. Keep those columns intact unless
the command reference says they are user-editable.

## Accept Values

The `accept` field is parsed deliberately but flexibly.

| Accepted as yes | Accepted as no |
| --- | --- |
| `yes`, `y`, `true`, `t`, `1`, `accept`, `accepted` | `no`, `n`, `false`, `f`, `0`, `reject`, `rejected`, empty, `.` |

Any other value is an error. That is useful: typos should stop an executor
instead of silently applying a row.

## Eval Modes And Executors

| Review table | Proposed action | Matching executor | What accepted rows can do |
| --- | --- | --- | --- |
| `<prefix>.fix_review.tsv` | `split_piece` | `chromo fix --reviewed-plan` | Emit exact reviewed source slices. |
| `<prefix>.scaffold_review.tsv` | `scaffold_gap` | `chromo scaffold --reviewed-plan` | Override matching junction gap lengths. |
| `<prefix>.gapfill_review.tsv` | `fill_path` | `chromo gapfill --reviewed-plan --apply` | Restrict graph fills to accepted, revalidated paths. |

Task-specific manual dashboards can load the same tables with `--review-table`.
The browser queue is a review view. The executor still makes the sequence
change.

## Fix Review Workflow

Create the table:

```bash
chromo eval fix \
  --assembly-fasta assembly.fa \
  --paf paf/sample.paf \
  --contigs suspect_contig_1 suspect_contig_2 \
  --mode conservative \
  --gfa assembly_graph.gfa \
  --read-paf reads_to_assembly.paf \
  --gaf reads_to_graph.gaf \
  --output-prefix review/sample.fix
```

Review `review/sample.fix.fix_review.tsv`.

Apply accepted rows:

```bash
chromo fix \
  --assembly-fasta assembly.fa \
  --reviewed-plan review/sample.fix.fix_review.tsv \
  --output-fasta results/sample.fixed.fa \
  --report results/sample.fixed.tsv
```

With a reviewed plan, `chromo fix` does not need coords, PAF, contig targets,
or `--all`. The accepted rows supply explicit source slices.

## Scaffold Review Workflow

Create the table:

```bash
chromo eval scaffold \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --gfa assembly_graph.gfa \
  --gaf reads_to_graph.gaf \
  --read-paf reads_to_assembly.paf \
  --output-prefix review/sample.scaffold
```

Apply accepted gap overrides:

```bash
chromo scaffold \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --reviewed-plan review/sample.scaffold.scaffold_review.tsv \
  --output-prefix results/sample.reviewed_scaffold
```

Accepted scaffold rows pin `gap_bp` for matching
`scaffold`/`left_contig`/`right_contig` junctions. Evidence columns remain
provenance. They do not reorder contigs, choose graph branches, or fill gaps.

## Gapfill Review Workflow

Create the table:

```bash
chromo eval gapfill \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --gfa assembly_graph.gfa \
  --gaf reads_to_graph.gaf \
  --ref-paf paf/sample.ref_vs_asm.paf \
  --read-paf reads_to_assembly.paf \
  --output-prefix review/sample.gapfill
```

For scaffold-first workflows, replace `--ordered-fasta --assignments` with
`--scaffold-fasta --agp` so the review table is built from AGP component
flanks around scaffold N gaps.

Apply accepted fill paths:

```bash
chromo gapfill \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --gfa assembly_graph.gfa \
  --reviewed-plan review/sample.gapfill.gapfill_review.tsv \
  --output-prefix results/sample.reviewed_gapfill \
  --apply
```

Gapfill rows are rechecked against the current scaffold, contig pair, path
nodes, and fillability status before sequence is inserted. Stale accepted paths
fail instead of being silently applied.

## Spreadsheet Review Habits

Use a spreadsheet to make review faster, but keep the table machine-readable:

1. Freeze the identifier columns and `accept`.
2. Filter by `action`, `status`, `confidence`, and evidence columns.
3. Edit `accept` and `notes` first.
4. Avoid changing provenance fields such as source contigs, path nodes, slices,
   or gap keys unless you are intentionally authoring a new reviewed row.
5. Keep rejected rows in the file when useful; unaccepted rows are ignored by
   executors but remain part of the audit trail.
6. Save as tab-delimited text, not a spreadsheet workbook.
7. Re-run the executor and read its report.

## Evidence Columns Are Context

Optional evidence fields can be very helpful, but their job is different from
the `accept` field.

| Evidence family | Examples | How to use it |
| --- | --- | --- |
| `graph_*` | GFA node status, graph paths, direct edges, unitig boundary context. | Decide whether graph topology supports or complicates the row. |
| `gaf_*` | Graph read traversal status, selected path support, alternate path support. | Check whether graph-aligned reads support the chosen path or a competitor. |
| `longread_*` | Spanning, split, bridge, edge, or nearby read counts. | Distinguish support at breakpoints and contig-end junctions. |
| Gap and overlap fields | `gap_bp`, `raw_inferred_gap_bp`, `overlap_class`, `risk_flags`. | Decide whether the row should stay algorithmic, be overridden, or be left unresolved. |

Evidence columns do not apply sequence changes by themselves.

## Stale-Row Protection

Reviewed executors validate more than `accept=yes`.

| Executor | Stale protection |
| --- | --- |
| `chromo fix --reviewed-plan` | Reads accepted `fix` task rows and applies explicit source slices from the table. Task and schema validation prevent wrong-table use. |
| `chromo scaffold --reviewed-plan` | Applies only accepted `scaffold_gap` rows matching the current scaffold, left contig, and right contig. |
| `chromo gapfill --reviewed-plan --apply` | Recomputes the current graph path and fillability for accepted rows; mismatched `path_nodes` or unfillable paths are rejected. |

Regenerate review tables after changing FASTA files, assignments, GFA, GAF,
path-search settings, or the primary coords/PAF evidence.

## Cheat Sheet

| If you want... | Do this |
| --- | --- |
| A row to apply | Set `accept=yes` and keep its task-specific key fields valid. |
| A row to stay visible but not apply | Set `accept=no` and add `notes`. |
| Manual visual review of the same rows | Open the table with task-specific `chromo manual`. |
| Gapfill from legacy plan TSV | Use `accept_fill=yes` in the gapfill plan, or use eval's shared `accept=yes` table. |
| Strict reproducibility | Keep the eval table, edited reviewed table, executor report, and output FASTA together. |

## Common Traps

Do not edit `schema`, `task`, or `event_id` casually. They are validation and
audit fields.

Do not use a `fix` table with scaffold or gapfill executors. Task validation
will reject mismatches.

Do not assume `accept=yes` overrides physics or graph validation. Gapfill still
refuses ambiguous, stale, unsequenced, oversized, or flank-mismatched paths.

Do not delete rejected rows if you need an audit trail. `accept=no` is clearer
than disappearing evidence.

Do not use spreadsheet software that silently changes large IDs, contig names,
or tab delimiters.

## What To Look At Next In ChromoSort

- Use [Manual Dashboard Review]({{ '/guides/manual-dashboard-review/' | relative_url }})
  when table rows need visual context.
- Use [Scaffolding, Gaps, And Overlaps]({{ '/guides/scaffolding-gaps-overlaps/' | relative_url }})
  for scaffold-junction review.
- Use [Reading ChromoSort Audit Tables]({{ '/guides/audit-tables/' | relative_url }})
  for row-level report interpretation.
- Use [chromo eval]({{ '/commands/eval/' | relative_url }}) for the complete
  review-table reference.
