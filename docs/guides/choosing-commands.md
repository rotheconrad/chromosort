---
title: Sort, Clean, Fix, Cut, Or Manual?
description: Choose the least invasive ChromoSort command that fits the evidence and review state.
---

# Sort, Clean, Fix, Cut, Or Manual?

Use this guide when you know something looks wrong or unfinished, but you are
not sure which ChromoSort command should act on it.

The safest habit is:

> Choose the least invasive command that matches the evidence you already have.

## The Core Idea

ChromoSort commands sit on a spectrum. Some only report evidence. Some filter
or reorder complete contigs. Some split, cut, invert, scaffold, or fill
sequence. The more a command changes sequence, the more explicit review you
should require before running it at scale.

Start with evidence and review state, not with the command name.

## Command Roles

| Command | Use it when... | It changes sequence by... | Main caution |
| --- | --- | --- | --- |
| `chromo plot` | You need to see alignment structure. | It does not write FASTA. | Plots describe the exact FASTA pair used to make the alignment. |
| `chromo sort` | Contigs should be assigned, filtered, ordered, and optionally whole-contig oriented. | It writes retained complete contigs in reference order. | It does not cut chimeric contigs. |
| `chromo clean` | The assembly is mostly correct and you want one conservative cleanup pass. | It combines initial sorting, targeted conservative fixing, and final ordering. | Re-align the cleaned FASTA before validation or further alignment-dependent steps. |
| `chromo eval` | You want a spreadsheet queue before applying fixes, scaffolds, or fills. | It does not write sequence outputs. | Accepted rows still need a matching executor. |
| `chromo fix` | Reviewed contigs need alignment-planned splitting. | It emits copied or split contig pieces. | Use targeted contigs or reviewed plans for production work. |
| `chromo cut` | You already know exact cut positions. | It cuts after exact 1-based positions. | It does not discover breakpoints from alignments. |
| `chromo manual` | A human should browse, remove, restore, split, invert, reorder, or label pieces. | Browser export or `manual apply` writes reviewed FASTA. | Re-align manual FASTA before treating it as a new assembly input. |
| `chromo gafprep` | Targeted GAF evidence is needed without aligning every read to the full graph. | It does not write FASTA; it writes a selected FASTQ, sanitized GFA, GraphAligner script, and audit TSVs. | It prepares GraphAligner inputs; GraphAligner still writes the actual GAF. |
| `chromo scaffold` | Final sorted contigs should become one record per reference sequence. | It joins ordered contigs with N gaps and optional reviewed overlap handling. | It needs the ordered FASTA and matching assignment TSV from the same sort run. |
| `chromo gapfill` | Reviewed graph paths should replace scaffold N gaps. | It can insert graph path sequence and trim the right flank when applied. | Apply only when graph path support and stale-row checks pass. |

## Decision Tree

### I Only Need To Understand The Pattern

Use `chromo plot` and the guide pages.

```bash
chromo plot \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.fa \
  --paf assembly.paf \
  --output-prefix plots/sample
```

No FASTA is changed. This is the right first move when you do not yet know
whether a signal is a misjoin, orientation difference, repeat, or reference
difference.

### I Need A Reference-Ordered Contig FASTA

Use `chromo sort`.

```bash
chromo sort \
  --assembly-fasta assembly.fa \
  --paf assembly.paf \
  --output-prefix results/sample \
  --orient-to-reference
```

This is appropriate when contigs are mostly trustworthy units and the main job
is to assign, filter, order, and optionally orient whole contigs.

### The Assembly Is Mostly Correct, But Needs Conservative Cleanup

Use `chromo clean`.

```bash
chromo clean \
  --assembly-fasta raw.fa \
  --paf raw.paf \
  --output-prefix results/sample.clean \
  --orient-to-reference
```

`clean` is useful when you want a standard cleanup pass: discard weak or
redundant contigs, run conservative fix planning on retained raw contigs, then
write a cleaned ordered FASTA. The cleaned FASTA needs fresh alignment evidence
for validation.

### A Contig Might Be Chimeric, But I Need Review First

Use `chromo eval fix` or `chromo manual fix`.

```bash
chromo eval fix \
  --assembly-fasta raw.fa \
  --paf raw.paf \
  --contigs suspect_contig \
  --mode conservative \
  --output-prefix review/sample.fix
```

Then edit or inspect `review/sample.fix.fix_review.tsv` before applying it with
`chromo fix --reviewed-plan`.

### I Have Reviewed Split Candidates

Use targeted `chromo fix` or a reviewed plan.

```bash
chromo fix \
  --assembly-fasta raw.fa \
  --reviewed-plan review/sample.fix.fix_review.tsv \
  --output-fasta results/sample.fixed.fa \
  --report results/sample.fixed.tsv
```

When a reviewed plan is supplied, the table provides the source contigs and
slices. Re-align `sample.fixed.fa` before sorting or plotting it as the current
assembly.

### I Know The Exact Cut Positions

Use `chromo cut`.

```bash
chromo cut \
  --assembly-fasta assembly.fa \
  --cut contigA:450000,910000 \
  --output-fasta results/sample.cut.fa \
  --report results/sample.cut.tsv
```

Cut positions are 1-based and mean "cut after this base." Terminal cuts are
rejected. This is the right tool for externally reviewed coordinates, not for
discovering breakpoints.

### The Case Needs Human Curation

Use `chromo manual`.

```bash
chromo manual \
  --ref-fasta reference.fa \
  --assembly-fasta raw.fa \
  --paf raw.paf \
  --gfa assembly_graph.gfa \
  --output-html review/sample.manual.html
```

Manual review is the best fit when there are several interacting decisions:
remove or restore contigs, split one region, invert a piece, label scaffold
groups, or compare graph and long-read evidence before exporting a recipe.

## Stage-Aware Workflow Patterns

### Conservative Cleanup Path

```text
raw.fa + raw.paf
  -> chromo clean
  -> clean.fa
  -> re-align clean.fa
  -> validation plots and final sort/scaffold decisions
```

Use this path for mostly correct assemblies where you want a compact first pass.

### Reviewed Chimeric-Contig Path

```text
raw.fa + raw.paf
  -> chromo sort reports and dot plots
  -> chromo eval fix or chromo manual fix
  -> chromo fix --reviewed-plan
  -> fixed.fa
  -> re-align fixed.fa
  -> chromo sort
```

Use this path when split candidates need traceable review.

### Final Scaffolding Path

```text
final_contigs.fa + final_contigs.paf
  -> chromo sort
  -> ordered.fa + contig_assignments.tsv
  -> chromo scaffold
  -> scaffold.fa
  -> optional chromo gapfill with reviewed graph paths
```

Use this path after contig-level decisions are settled.

## Cheat Sheet

| Situation | First command to consider |
| --- | --- |
| Need to inspect patterns | `chromo plot` |
| Need ordered retained contigs | `chromo sort` |
| Need a conservative one-pass cleanup | `chromo clean` |
| Need a review spreadsheet | `chromo eval` |
| Need to split suspected chimeras from alignment evidence | `chromo fix` |
| Need to cut exact known coordinates | `chromo cut` |
| Need a browser review dashboard | `chromo manual` |
| Need chromosome-scale records from final ordered contigs | `chromo scaffold` |
| Need graph-supported sequence through N gaps | `chromo gapfill` |

## Common Traps

Do not use `chromo fix` just because a contig is reversed. Whole-contig reverse
orientation is usually a sort or manual orientation decision.

Do not use `chromo cut` when you only know that a contig looks suspicious. Cut
needs exact reviewed positions.

Do not run sequence-changing commands and then reuse old alignment evidence as
if it describes the new FASTA.

Do not treat `chromo clean` as a substitute for validation. It is a conservative
cleanup command, but its output still needs fresh plots or alignments.

Do not scaffold before the ordered FASTA and assignment TSV come from the same
final sort run.

## What To Look At Next In ChromoSort

- Use [Alignment Evidence And The Exact FASTA Rule]({{ '/guides/alignment-evidence/' | relative_url }})
  before chaining commands across FASTA stages.
- Use [Sorting Decisions And Duplicate-Overlap Filtering]({{ '/guides/sorting-decisions/' | relative_url }})
  when deciding whether a sort report is acceptable.
- Use [Chimeric Contig And Breakpoint Review]({{ '/guides/breakpoint-review/' | relative_url }})
  before applying `chromo fix`.
- Use the [command reference]({{ '/commands/' | relative_url }}) for exact
  parameters and outputs.
