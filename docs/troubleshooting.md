---
title: Troubleshooting
description: Common ChromoSort setup, input, plotting, graph, and gapfill checks.
---

# Troubleshooting

This page collects first-pass checks for common setup, alignment, and graph-review problems.

## `chromo` Command Not Found

Activate the environment that installed ChromoSort:

```bash
mamba activate chromosort
chromo --help
```

If you are using Pixi, run commands through Pixi:

```bash
pixi run help
pixi run test
```

## MUMmer Coordinates Do Not Parse

Use `show-coords -r -c -l` on a filtered delta file. ChromoSort expects reference names, query names, reference/query coordinates, lengths, coverage, percent identity, and sequence names.

```bash
show-coords -r -c -l mummer/sample.filter > mummer/sample.coords
```

## PAF Rows Are Missing

By default ChromoSort skips secondary PAF rows marked `tp:A:S`. Add `--include-secondary-paf` only when secondary rows are part of the review plan. For noisy alignments, also check `--min-mapq`.

## Edited FASTA Does Not Match Old Alignments

After `chromo sort` writes `ordered.fa`, `chromo clean` writes `clean.fa`,
`chromo fix` writes `fixed.fa`, `chromo cut` writes a cut FASTA,
`chromo manual apply` writes a manual FASTA, or `chromo scaffold` writes
scaffold records, re-align the edited FASTA before running another
alignment-dependent command on that edited FASTA. Old raw assembly alignments
usually do not describe the edited sequence names, membership, orientation, or
coordinates.

It is fine to reuse raw coords or PAF for another decision about the same raw
assembly, such as running `chromo fix --assembly-fasta raw.fa --coords
raw.coords` after inspecting `chromo sort --assembly-fasta raw.fa --coords
raw.coords`. It is not safe to run that same `raw.coords` against
`sample.ordered.fa` or `sample.fixed.fa`.

## Plot Shows The Old Assembly

`chromo plot --assignments` draws the original alignment rows and uses a
`chromo sort` assignment table to order the query axis. This is useful for
reviewing sort decisions. It does not validate `ordered.fa`, `fixed.fa`, or a
manual FASTA unless the coords or PAF were generated from that exact FASTA. The
same rule applies to `clean.fa` from `chromo clean`.

## GFA Nodes Are Reported Missing

Graph-aware commands expect GFA segment names to match the assembly FASTA or the original contig names in the assignment report. If another tool renamed, polished, split, or scaffolded the FASTA after graph export, keep a name map or regenerate graph evidence for the renamed sequences.

## Gapfill Plan Is Ambiguous

`chromo gapfill` refuses to guess through ambiguous graph branches. Add GAF,
Hi-C-like pair evidence, or reference-placement PAF when those evidence layers
are trustworthy, then review the candidate-path table, `chromo eval gapfill`
table, or HTML reviewer. Ties, weak support, conflicting evidence, missing
sequence, invalid overlaps, and flank mismatches intentionally remain
unresolved.

## Plots Are Empty Or Sparse

Check that the reference and assembly FASTA IDs match the alignment file exactly. If you applied strict filters upstream, try plotting the unfiltered or less strictly filtered coords/PAF file first, then add `--min-segment-bp`, `--min-segment-idy`, or `--min-mapq` only after confirming the expected rows are present.

If the plot is not empty but the pattern is hard to classify, use
[How to Interpret Dot Plots]({{ '/dot-plots/' | relative_url }}) to compare it
against common clean, reversed, chimeric, inversion, duplicate, gap, and
repeat-like examples.
