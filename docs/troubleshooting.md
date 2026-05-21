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

## Fixed FASTA Does Not Match Old Alignments

After `chromo fix`, `chromo cut`, or `chromo manual apply`, re-align the edited FASTA before running `chromo sort`, `chromo plot`, `chromo scaffold`, or `chromo gapfill` with coordinate-dependent reports. Old raw-assembly alignments usually do not describe the edited sequence names or coordinates.

## GFA Nodes Are Reported Missing

Graph-aware commands expect GFA segment names to match the assembly FASTA or the original contig names in the assignment report. If another tool renamed, polished, split, or scaffolded the FASTA after graph export, keep a name map or regenerate graph evidence for the renamed sequences.

## Gapfill Plan Is Ambiguous

`chromo gapfill` refuses to guess through ambiguous graph branches. Add GAF, Hi-C-like pair evidence, or reference-placement PAF when those evidence layers are trustworthy, then review the candidate-path table or HTML reviewer. Ties, weak support, conflicting evidence, missing sequence, invalid overlaps, and flank mismatches intentionally remain unresolved.

## Plots Are Empty Or Sparse

Check that the reference and assembly FASTA IDs match the alignment file exactly. If you applied strict filters upstream, try plotting the unfiltered or less strictly filtered coords/PAF file first, then add `--min-segment-bp`, `--min-segment-idy`, or `--min-mapq` only after confirming the expected rows are present.
