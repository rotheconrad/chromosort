---
title: chromo eval
description: Usage, outputs, parameters, and reasoning for chromo eval review tables.
---

# chromo eval

`chromo eval` prepares task-specific TSV review tables. These tables are meant
for spreadsheet-first review workflows where users want the algorithm to
propose candidate decisions, then keep, reject, delete, or add rows before a
sequence-changing command applies the reviewed decisions.

At present, `chromo eval fix` is available. Scaffold and gapfill review-table
modes are planned in the production roadmap.

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

The planner threshold options mirror `chromo fix`.

## Reasoning

`eval` is evidence-first and table-only. It does not change FASTA records. The
sequence-changing step remains explicit: review the table, then pass it to the
corresponding executor command.
