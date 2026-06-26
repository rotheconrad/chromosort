---
title: Choosing PAF Or MUMmer Coords
description: Practical guide for choosing minimap2 PAF or MUMmer coords as ChromoSort alignment evidence.
---

# Choosing PAF Or MUMmer Coords

Use this guide when you are deciding which whole-genome alignment format should
be the primary evidence for a ChromoSort run.

For most new ChromoSort runs, minimap2 PAF is the recommended starting point.
It is fast, compact, supports MAPQ filtering, and works across `chromo sort`,
`chromo fix`, `chromo manual`, and `chromo plot`. MUMmer coords remains a good
choice when a project already has a tuned nucmer workflow, when you want a
second aligner view, or when minimap2 gives a surprising result.

## The Core Idea

ChromoSort normalizes MUMmer coords and minimap2 PAF rows into the same internal
alignment model before sorting, plotting, and fixing. The remaining differences
usually come from the aligners, not from separate ChromoSort decision logic.

Expect differences in:

- row fragmentation,
- chaining behavior,
- primary versus secondary alignment handling,
- MAPQ availability,
- identity fields,
- repeat and paralog placement.

Those differences are useful evidence. They should not be treated as automatic
proof that one file is wrong.

## What ChromoSort Reads

| Format | ChromoSort uses | Most useful when |
| --- | --- | --- |
| MUMmer coords | Reference/query names, reference/query coordinates, lengths, row lengths, percent identity, and orientation. | You already have filtered nucmer outputs, want a second aligner perspective, or are comparing against older MUMmer-based workflows. |
| minimap2 PAF | Query name/length/start/end, target name/length/start/end, strand, matching bases, block length, and MAPQ. | You want a fast default input with MAPQ filtering and compact whole-genome alignments. |

Provide exactly one of `--coords` or `--paf` to commands that accept alignment
evidence.

## Starting Recommendation

For same-species or close-reference production runs, start with:

```bash
minimap2 \
  -x asm5 \
  -c \
  -t 16 \
  --secondary=no \
  reference.fa \
  assembly.fa \
  > paf/sample.ref_vs_asm.paf
```

Then use that PAF anywhere you would otherwise use MUMmer coords:

```bash
chromo sort \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.fa \
  --paf paf/sample.ref_vs_asm.paf \
  --output-prefix results/sample
```

The `-c` flag matters. It asks minimap2 to perform base-level alignment and
write CIGAR-bearing PAF rows. ChromoSort does not parse the CIGAR string
directly, but the base-level alignment changes the PAF match and block-length
columns used for identity summaries and identity filtering.

## Preset Gallery

Choose the strictest minimap2 assembly preset that recovers the expected
chromosome-scale alignments.

| Preset | Consider it when | Watch for |
| --- | --- | --- |
| `asm5` | Same species, same breeding pool, cultivar or line comparisons, or a new assembly against a close reference. | If expected contigs are missing or fragmented, the preset may be too stringent. |
| `asm10` | More divergent same-species material, wild or exotic accessions, pangenome references, or close relatives where `asm5` misses syntenic blocks. | More repeat, paralog, or ambiguous signal may appear; inspect plots and best-reference shares. |
| `asm20` | Related species or highly divergent reference choices where `asm10` still misses expected chromosome-scale alignments. | This is the noisiest choice for ChromoSort and may require stricter MAPQ or segment filters. |

The presets are starting points, not hard biological cutoffs. If a preset
changes the conclusion for a marginal contig, review the plot and reports
rather than choosing by habit.

## When MUMmer Coords Is The Better First Choice

Use MUMmer coords first when:

- the project already has a trusted nucmer plus `delta-filter` workflow,
- the reference and assembly are close enough that your MUMmer settings recover
  the expected chromosome-scale alignments cleanly,
- you need continuity with previous reports or published workflows,
- PAF results look surprising and you want a second aligner perspective.

A typical coords export is:

```bash
nucmer \
  -t 16 \
  -c 500 \
  -p mummer/sample \
  reference.fa \
  assembly.fa

delta-filter \
  -i 95 \
  -l 10000 \
  -1 \
  mummer/sample.delta \
  > mummer/sample.filter

show-coords \
  -r \
  -c \
  -l \
  mummer/sample.filter \
  > mummer/sample.coords
```

## Interpreting Disagreements

Small PAF-vs-coords disagreements are expected. In soybean `chromo fix`
testing, split counts differed by about 5-10%, while marginal split-contig sets
differed by about 20-30% depending on mode. Use those numbers as practical
expectations, not universal guarantees.

When the two formats disagree, ask:

- Did both files align the same exact FASTA pair?
- Was minimap2 run with `-c --secondary=no`?
- Were secondary PAF rows included or skipped?
- Are low-MAPQ PAF rows driving the call?
- Did MUMmer filtering remove short or lower-identity blocks that PAF retained?
- Are the discordant blocks long, coherent, and visible on a per-reference plot?
- Does graph, GAF, or long-read PAF evidence support either interpretation?

Use `chromo eval` and `chromo manual` when the decision affects sequence.

## Filter Cheat Sheet

| Signal | Useful filter or setting |
| --- | --- |
| Too many small repeat-like PAF rows | Increase `--min-segment-bp` or keep `--secondary=no`. |
| Low-confidence PAF placements | Increase `--min-mapq`. |
| PAF identity looks unexpectedly low | Confirm minimap2 was run with `-c`. |
| Expected syntenic blocks are missing | Try `asm10` after `asm5`, then inspect plots before using `asm20`. |
| MUMmer output is very fragmented | Review `delta-filter -i`, `-l`, and `-1` settings. |
| Coords and PAF disagree on a marginal split | Use dot plots, `chromo eval`, graph evidence, and long-read evidence before editing sequence. |

## Common Traps

Do not run both formats and then silently mix decisions. Choose one primary
alignment for a workflow stage, then use the other as review evidence when
needed.

Do not assume `asm20` is better because it is more permissive. More alignments
can mean more repeat noise and more ambiguous placements.

Do not include secondary PAF rows by default. ChromoSort skips rows marked
`tp:A:S` unless the relevant `--include-secondary-*` option is set, and
`--secondary=no` keeps the file cleaner from the start.

Do not compare old coords against new PAF from a different FASTA stage. Format
choice and FASTA-stage choice are separate questions.

Do not set identity filters before looking at the identity distribution from
your chosen aligner and preset.

## What To Look At Next In ChromoSort

- Use [Alignment Evidence And The Exact FASTA Rule]({{ '/guides/alignment-evidence/' | relative_url }})
  before reusing any alignment file.
- Use [Input Files]({{ '/input-files/' | relative_url }}#creating-input-files-with-minimap2)
  for the current minimap2 recipe.
- Use [Troubleshooting]({{ '/troubleshooting/' | relative_url }}#coords-and-paf-disagree)
  when the two formats disagree.
- Use [How to Interpret Dot Plots]({{ '/dot-plots/' | relative_url }}) to
  decide whether a disagreement is biologically meaningful.
