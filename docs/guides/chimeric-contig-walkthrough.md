---
title: Suspected Chimeric Contig Review
description: A worked review path from dot-plot signal to eval fix, manual review, reviewed split, re-alignment, sort, and validation plot.
---

# Suspected Chimeric Contig Review

Use this walkthrough when a raw dot plot or sort report suggests that one
contig may join sequence from different references or incompatible regions.

The goal is:

> Review the candidate as raw-assembly evidence, apply only accepted edits, then
> re-align the edited FASTA before sorting again.

## Step 1: Start With Raw Evidence

Draw raw plots from the exact raw FASTA and raw alignment:

```bash
chromo plot \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.raw.fa \
  --paf paf/raw.ref_vs_asm.paf \
  --output-prefix plots/raw \
  --per-ref
```

If the suspicious signal is limited to a few references, use `--sel-ref` to
keep the review focused.

Look for:

| Pattern | First interpretation |
| --- | --- |
| One contig has strong blocks on two references | Possible misjoin, translocation, repeat, or reference difference. |
| One clean reverse-strand contig | Usually orientation, not a split. |
| Internal reverse block on one reference | Possible inversion; review carefully before editing. |
| Many tiny off-target hits | Often repeats or aligner noise. |

## Step 2: Use Sort Reports As A Candidate List

Run sort as a review pass if you need table evidence:

```bash
chromo sort \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.raw.fa \
  --paf paf/raw.ref_vs_asm.paf \
  --output-prefix review/raw_sort
```

Inspect:

- `review/raw_sort.contig_assignments.tsv`,
- `review/raw_sort.contig_ref_matches.tsv`,
- any rows with `kept_split_candidate`, `ambiguous_ref_match`, or strong
  secondary reference support.

This sort run helps select contigs. It does not split them.

## Step 3: Create A Fix Review Table

For one or a few suspicious contigs:

```bash
chromo eval fix \
  --assembly-fasta assembly.raw.fa \
  --paf paf/raw.ref_vs_asm.paf \
  --contigs suspect_contig_1 suspect_contig_2 \
  --mode conservative \
  --gfa assembly_graph.gfa \
  --read-paf reads/raw.reads_to_assembly.paf \
  --gaf graph_reads/raw.reads_to_graph.gaf \
  --output-prefix review/sample.fix
```

The output is:

```text
review/sample.fix.fix_review.tsv
```

Review `split_piece` rows beside the dot plot. Keep `accept=no` for rows that
look weak, repetitive, inverted but biological, or unsupported by independent
evidence.

## Step 4: Use Manual Review For Hard Cases

Open the same review table in a browser dashboard when context matters:

```bash
chromo manual fix \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.raw.fa \
  --paf paf/raw.ref_vs_asm.paf \
  --review-table review/sample.fix.fix_review.tsv \
  --gfa assembly_graph.gfa \
  --read-paf reads/raw.reads_to_assembly.paf \
  --gaf graph_reads/raw.reads_to_graph.gaf \
  --output-html review/sample.manual_fix.html
```

Use the dashboard to compare dot-plot blocks, graph neighborhoods, read
evidence fields, and candidate pieces. Export notes or a manual recipe if the
right edit is not captured by the table.

## Step 5: Apply Accepted Split Rows

When the table is reviewed:

```bash
chromo fix \
  --assembly-fasta assembly.raw.fa \
  --reviewed-plan review/sample.fix.fix_review.tsv \
  --output-fasta results/sample.fixed.fa \
  --report results/sample.fixed.tsv
```

With `--reviewed-plan`, `chromo fix` uses accepted source slices from the table.
It does not need the alignment file, `--contigs`, or `--all`.

If exact cut positions are known from another review source, use `chromo cut`
instead of asking the planner to rediscover them.

## Step 6: Re-Align The Fixed FASTA

The fixed FASTA has new records and new query coordinates:

```bash
minimap2 -x asm5 -c -t 16 --secondary=no \
  reference.fa results/sample.fixed.fa \
  > paf/sample.fixed.ref_vs_asm.paf
```

Validate with a fresh plot:

```bash
chromo plot \
  --ref-fasta reference.fa \
  --assembly-fasta results/sample.fixed.fa \
  --paf paf/sample.fixed.ref_vs_asm.paf \
  --output-prefix plots/sample.fixed \
  --per-ref
```

The fixed plot should show that accepted pieces now place coherently. If not,
return to review rather than stacking more edits on stale evidence.

## Step 7: Sort The Fixed Assembly

```bash
chromo sort \
  --ref-fasta reference.fa \
  --assembly-fasta results/sample.fixed.fa \
  --paf paf/sample.fixed.ref_vs_asm.paf \
  --output-prefix results/sample.fixed \
  --orient-to-reference
```

Now scaffold or gapfill only from this fixed-stage sort output.

## Decision Log

Keep a short decision table beside the output:

| Contig | Evidence | Decision | Output |
| --- | --- | --- | --- |
| `suspect_contig_1` | Two strong reference blocks, accepted split rows, read support near breakpoints. | Apply reviewed fix. | `sample.fixed.fa` |
| `suspect_contig_2` | Internal inversion with read/graph support. | Leave native. | No split. |
| `suspect_contig_3` | Weak off-target repeat hits only. | Reject split rows. | No split. |

This is small, but it prevents future reviewers from re-litigating the same
rows without context.

## Common Traps

Do not run `chromo fix --all` and skip review on a new dataset.

Do not split whole-contig reverse orientation. Use sorting orientation when
that is the goal.

Do not treat same-reference inversions as errors by default.

Do not use raw alignment files to validate `sample.fixed.fa`.

Do not lose the reviewed TSV. It is the reason the sequence changed.

## What To Look At Next In ChromoSort

- Use [Chimeric Contig And Breakpoint Review]({{ '/guides/breakpoint-review/' | relative_url }})
  for planner modes and breakpoint status labels.
- Use [Manual Dashboard Review]({{ '/guides/manual-dashboard-review/' | relative_url }})
  for browser review and recipes.
- Use [Long-Read PAF And GAF Support]({{ '/guides/long-read-paf-gaf/' | relative_url }})
  when read evidence is present.
- Use [chromo eval]({{ '/commands/eval/' | relative_url }}) and
  [chromo fix]({{ '/commands/fix/' | relative_url }}) for command details.
