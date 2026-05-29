---
title: chromo eval
description: Usage, outputs, parameters, and reasoning for chromo eval review tables.
---

# chromo eval

`chromo eval` prepares task-specific TSV review tables. These tables are meant
for spreadsheet-first review workflows where users want the algorithm to
propose candidate decisions, then keep, reject, delete, or add rows before a
sequence-changing command applies the reviewed decisions.

At present, `chromo eval fix`, `chromo eval scaffold`, and
`chromo eval gapfill` are available.

## Run `chromo eval fix`

```bash
chromo eval fix \
  --assembly-fasta assembly.fa \
  --coords mummer/raw.coords \
  --contigs suspect_contig_1 suspect_contig_2 \
  --output-prefix results/sample.eval_fix
```

The command writes:

```text
results/sample.eval_fix.fix_review.tsv
```

Each accepted `split_piece` row can be applied with:

```bash
chromo fix \
  --assembly-fasta assembly.fa \
  --reviewed-plan results/sample.eval_fix.fix_review.tsv \
  --output-fasta results/sample.fixed.fa \
  --report results/sample.fixed.tsv
```

When `--reviewed-plan` is used, `chromo fix` does not require `--coords`,
`--paf`, `--contigs`, or `--all`. The reviewed table supplies the exact source
contigs and slices to apply.

## Optional Evidence

`chromo eval fix` can include report-only context from:

- `--gfa`: assembly graph node status and local complexity fields.
- `--read-paf`: long-read-to-assembly PAF evidence around candidate
  breakpoints.

Long-read evidence is summarized as spanning-read, split-read, edge-read, and
nearby-read counts for candidate breakpoint rows.

## Run `chromo eval scaffold`

```bash
chromo eval scaffold \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --output-prefix results/sample.eval_scaffold \
  --gfa assembly_graph.gfa \
  --read-paf reads_to_assembly.paf
```

The command writes:

```text
results/sample.eval_scaffold.scaffold_review.tsv
```

Each accepted `scaffold_gap` row can override the gap length for one matching
adjacent scaffold junction:

```bash
chromo scaffold \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --reviewed-plan results/sample.eval_scaffold.scaffold_review.tsv \
  --output-prefix results/sample.reviewed_scaffold
```

`chromo scaffold` still builds every ordinary junction algorithmically. The
reviewed table only pins accepted rows whose `scaffold`, `left_contig`, and
`right_contig` match the current inputs. Stale accepted rows are rejected so an
old spreadsheet cannot silently change a new scaffold.

`chromo eval scaffold` can include report-only context from:

- `--gfa`: direct-link and short-path graph context for adjacent contigs.
- `--gaf`: long-read graph traversal support for the direct or short GFA path.
- `--read-paf`: long-read-to-assembly bridge evidence between contig ends.

Long-read scaffold evidence is summarized as bridge-read counts, orientation
summaries, read-order summaries, and median read-space gap estimates.

## Run `chromo eval gapfill`

```bash
chromo eval gapfill \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --gfa assembly_graph.gfa \
  --gaf reads_to_graph.gaf \
  --ref-paf paf/sample.ref_vs_asm.paf \
  --read-paf reads_to_assembly.paf \
  --output-prefix results/sample.eval_gapfill
```

The command writes:

```text
results/sample.eval_gapfill.gapfill_review.tsv
```

Each accepted `fill_path` row can be applied with:

```bash
chromo gapfill \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --gfa assembly_graph.gfa \
  --reviewed-plan results/sample.eval_gapfill.gapfill_review.tsv \
  --output-prefix results/sample.reviewed_gapfill \
  --apply
```

`chromo eval gapfill` uses the same graph path planner as `chromo gapfill`.
Fillable rows are accepted by default in the review-event table; unresolved,
ambiguous, stale, or unfillable rows remain visible but are not accepted by
default. `chromo gapfill --reviewed-plan` validates accepted rows against the
current scaffold, contig pair, and `path_nodes` before applying sequence.

`chromo eval gapfill` can include report-only context from:

- `--gaf`: read-path support for candidate GFA paths.
- `--hic-pairs`: graph-node contact support.
- `--ref-paf`: reference-placement support for intermediate graph nodes.
- `--read-paf`: long-read-to-assembly bridge evidence between contig ends.

## Parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--assembly-fasta` | required | Assembly FASTA containing contigs to evaluate. |
| `--coords` / `--paf` | required | Whole-genome reference-to-assembly alignment used by the fix planner. |
| `--contigs`, `--contigs-file`, `--all` | none | Evaluate selected contigs or all split-signal contigs. |
| `--output-prefix` | required | Prefix for `<prefix>.fix_review.tsv`. |
| `--gfa` | none | Optional assembly graph GFA for node context fields. |
| `--read-paf` | none | Optional long-read-to-assembly PAF for breakpoint support fields. |
| `--mode` | `conservative` | Fix planner mode used to prepare candidate rows. |
| `--ordered-fasta` | required for scaffold mode | Final ordered FASTA from `chromo sort`. |
| `--assignments` | required for scaffold mode | Matching `<prefix>.contig_assignments.tsv` from `chromo sort`. |
| `--fixed-gap-bp` | none | In scaffold mode, prepare rows using fixed gaps instead of inferred gaps. |
| `--gfa` | required for gapfill mode | Assembly graph GFA for candidate fill paths. |
| `--gaf` | none | Optional long-read graph traversal evidence for scaffold and gapfill modes. Scaffold mode requires `--gfa` when `--gaf` is provided. |
| `--hic-pairs`, `--ref-paf` | none | Optional graph-path support evidence for gapfill mode. |
| `--include-fill-sequences` | off | In gapfill mode, include candidate fill sequence in the review table. |

The planner threshold options mirror the corresponding executor command.

## Reasoning

`eval` is evidence-first and table-only. It does not change FASTA records. The
sequence-changing step remains explicit: review the table, then pass it to the
corresponding executor command.
