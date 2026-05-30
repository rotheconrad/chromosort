---
title: chromo eval
description: Usage, outputs, parameters, and reasoning for chromo eval review tables.
---

# chromo eval

`chromo eval` prepares task-specific TSV review tables. These tables are meant
for spreadsheet-first review workflows where users want the algorithm to
propose candidate decisions, then keep, reject, delete, or add rows before a
sequence-changing command applies the reviewed decisions.

The available modes are `chromo eval fix`, `chromo eval scaffold`, and
`chromo eval gapfill`. Each mode reuses the algorithm or planner from the
matching executor, then writes a task-scoped `chromosort-review-event-v1` table.

| Eval mode | Reuses algorithm from | Reviewed executor path | Optional evidence columns |
| --- | --- | --- | --- |
| `fix` | `chromo fix` breakpoint planner | `chromo fix --reviewed-plan` | `graph_*`, `gaf_*`, `longread_*` |
| `scaffold` | `chromo scaffold` gap/overlap and graph-junction reporting | `chromo scaffold --reviewed-plan` | `graph_*`, `gaf_*`, `longread_*` |
| `gapfill` | `chromo gapfill` path planner and fillability checks | `chromo gapfill --reviewed-plan --apply` | GAF, Hi-C, reference-placement PAF, and long-read bridge fields |

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
- `--gaf`: advisory graph traversal context for candidate contig graph nodes.
- `--read-paf`: long-read-to-assembly PAF evidence around candidate
  breakpoints.

Long-read evidence is summarized as spanning-read, split-read, edge-read, and
nearby-read counts for candidate breakpoint rows.

For same-reference inversion review, use `--mode comprehensive` to ask the fix
planner to expose orientation-change rows in the review table without applying
them. Pair the table with `--read-paf`, `--gfa`, and `--gaf` when you need to
decide whether the inversion is a real assembly feature or an assembly error.
The [Agent and Review Playbook]({{ '/review-playbook/' | relative_url }}#review-same-reference-inversions)
describes this distinction and why real inversions should usually remain native
for pangenome graph inputs.

## Run `chromo eval scaffold`

```bash
chromo eval scaffold \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --output-prefix results/sample.eval_scaffold \
  --gfa assembly_graph.gfa \
  --gaf reads_to_graph.gaf \
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
GAF scaffold evidence reports the first selected GFA path, its supporting read
count, the best alternate path and count when present, and a compact status such
as `supports_selected`, `supports_alternate`, `tied_support`, or `no_support`.
These fields are advisory for `chromo scaffold`; accepted reviewed rows pin gap
lengths, not graph branches.

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

| Mode | Parameter | Default | Meaning |
| --- | --- | --- | --- |
| `fix` | `--assembly-fasta` | required | Assembly FASTA containing contigs to evaluate. |
| `fix` | `--coords` / `--paf` | required | Whole-genome reference-to-assembly alignment used by the fix planner. |
| `fix` | `--contigs`, `--contigs-file`, `--all` | none | Evaluate selected contigs or all split-signal contigs. |
| `fix` | `--mode` | `conservative` | Fix planner mode used to prepare candidate rows. |
| `scaffold`, `gapfill` | `--ordered-fasta` | required | Final ordered FASTA from `chromo sort`. |
| `scaffold`, `gapfill` | `--assignments` | required | Matching `<prefix>.contig_assignments.tsv` from `chromo sort`. |
| `scaffold` | `--fixed-gap-bp` | none | Prepare rows using fixed gaps instead of inferred gaps. |
| `scaffold` | `--gfa` | none | Optional assembly graph GFA for direct-link and short-path junction context; required when `--gaf` is provided. |
| `gapfill` | `--gfa` | required | Assembly graph GFA for candidate fill paths. |
| all modes | `--output-prefix` | required | Prefix for the mode-specific review TSV. |
| all modes | `--read-paf` | none | Optional long-read-to-assembly PAF for breakpoint or contig-end bridge support fields. |
| all modes | `--gaf` | none | Optional long-read graph traversal evidence. In `fix`, GAF is advisory node context; in `scaffold` and `gapfill`, it reports candidate path support. |
| `gapfill` | `--hic-pairs`, `--ref-paf` | none | Optional graph-path support evidence for gapfill branch resolution. |
| `gapfill` | `--include-fill-sequences` | off | Include candidate fill sequence in the review table. |

The planner threshold options mirror the corresponding executor command.

## Reasoning

`eval` is evidence-first and table-only. It does not change FASTA records. The
sequence-changing step remains explicit: review the table, then pass it to the
corresponding executor command.
