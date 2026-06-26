---
title: Troubleshooting
description: Common ChromoSort setup, input, plotting, graph, and gapfill checks.
---

# Troubleshooting

This page collects first-pass checks for common setup, alignment, and graph-review problems.

## Troubleshooting Entry Points

Use these entry points when a command runs but the evidence, plot, or reviewed
plan does not look right.

| Symptom | First check | Deeper guide |
| --- | --- | --- |
| `chromo` is not available | [Activate the installed environment](#chromo-command-not-found) | [Installation]({{ '/installation/' | relative_url }}) |
| MUMmer coords fail to parse or PAF rows disappear | [Check coords and PAF format assumptions](#mummer-coordinates-do-not-parse) | [Choosing PAF or MUMmer coords]({{ '/guides/paf-vs-coords/' | relative_url }}) |
| Coords and PAF disagree | [Compare aligner settings and filters](#coords-and-paf-disagree) | [Choosing PAF or MUMmer coords]({{ '/guides/paf-vs-coords/' | relative_url }}) |
| Edited FASTA was paired with old coords or PAF | [Re-align the exact edited FASTA](#edited-fasta-does-not-match-old-alignments) | [Alignment evidence and the exact FASTA rule]({{ '/guides/alignment-evidence/' | relative_url }}) |
| Plot shows an old assembly | [Check plot provenance](#plot-shows-the-old-assembly) | [Alignment evidence and the exact FASTA rule]({{ '/guides/alignment-evidence/' | relative_url }}) |
| Plot is empty or sparse | [Check plot inputs and sequence names](#plots-are-empty-or-sparse) | [FASTA and evidence name matching]({{ '/guides/name-matching/' | relative_url }}) |
| Audit or review status is confusing | Compare status labels before changing sequence | [Reading ChromoSort audit tables]({{ '/guides/audit-tables/' | relative_url }}) |
| GFA nodes, overlays, or projections are missing | [Check graph and FASTA coordinate systems](#gfa-overlay-or-projection-is-empty) | [hifiasm unitig-to-contig projection]({{ '/guides/hifiasm-graph-projection/' | relative_url }}) |
| Gapfill candidates remain unresolved | [Review ambiguity and support evidence](#gapfill-plan-is-ambiguous) | [Graph-supported gap filling]({{ '/guides/graph-gap-filling/' | relative_url }}) |

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

For broader alignment-evidence choices, see
[Choosing PAF or MUMmer coords]({{ '/guides/paf-vs-coords/' | relative_url }}).

## Coords And PAF Disagree

Small disagreements are expected. In soybean `chromo fix` testing, coords and
PAF split counts differed by about 5-10%, while marginal split-contig sets
differed by about 20-30%. ChromoSort normalizes both formats before decision
logic, so first check aligner settings, minimap2 preset, `-c --secondary=no`,
MAPQ filters, MUMmer `delta-filter` settings, row counts, and dot plots. Use
`chromo eval` with long-read PAF, GFA, or GAF for stronger event evidence.
The guide-level comparison is
[Choosing PAF or MUMmer coords]({{ '/guides/paf-vs-coords/' | relative_url }}).

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

For the full rule and stage-by-stage examples, see
[Alignment evidence and the exact FASTA rule]({{ '/guides/alignment-evidence/' | relative_url }}).

## Plot Shows The Old Assembly

`chromo plot --assignments` draws the original alignment rows and uses a
`chromo sort` assignment table to order the query axis. This is useful for
reviewing sort decisions. It does not validate `ordered.fa`, `fixed.fa`, or a
manual FASTA unless the coords or PAF were generated from that exact FASTA. The
same rule applies to `clean.fa` from `chromo clean`.

Use [How to interpret dot plots]({{ '/dot-plots/' | relative_url }}) for plot
patterns, and
[Alignment evidence and the exact FASTA rule]({{ '/guides/alignment-evidence/' | relative_url }})
when the plotted FASTA pair is uncertain.

## GFA Nodes Are Reported Missing

Graph-aware commands expect GFA segment names to match the assembly FASTA or the original contig names in the assignment report. If another tool renamed, polished, split, or scaffolded the FASTA after graph export, keep a name map or regenerate graph evidence for the renamed sequences.

See [FASTA and evidence name matching]({{ '/guides/name-matching/' | relative_url }})
and [Assembly graph evidence]({{ '/guides/assembly-graph-evidence/' | relative_url }})
for the guide-level view of these checks.

## GFA Overlay Or Projection Is Empty

For hifiasm outputs, check whether your FASTA and GFA use the same coordinate
system. Dot plots made from `p_ctg.fa` or `hap*.p_ctg.fa` are in contig
coordinates, while `p_utg.gfa`, `r_utg.gfa`, and their `.noseq.gfa`
equivalents are usually in unitig coordinates. ChromoSort can project unitigs
onto contigs only when the GFA includes `P` path records or `W` walk records
whose names match the contig FASTA records.

Run `chromo graph-map` first when you are unsure:

```bash
chromo graph-map \
  --ctg-fasta assembly.p_ctg.fa \
  --utg-gfa assembly.p_utg.noseq.gfa \
  --output-prefix review/sample.graphmap
```

If the warning table says paths are missing, the GFA may contain only `S`, `L`,
and read-alignment `A` records. Those files are still useful for topology and
junction context, but they cannot define contig-axis unitig intervals for
`chromo plot --gfa-overlay`.

See
[hifiasm unitig-to-contig projection]({{ '/guides/hifiasm-graph-projection/' | relative_url }})
for the coordinate-system model behind this warning.

## Gapfill Plan Is Ambiguous

`chromo gapfill` refuses to guess through ambiguous graph branches. Add GAF,
Hi-C-like pair evidence, or reference-placement PAF when those evidence layers
are trustworthy, then review the candidate-path table, `chromo eval gapfill`
table, or HTML reviewer. Ties, weak support, conflicting evidence, missing
sequence, invalid overlaps, and flank mismatches intentionally remain
unresolved.

See [Graph-supported gap filling]({{ '/guides/graph-gap-filling/' | relative_url }})
for the candidate-path and risk-flag interpretation workflow.

## Plots Are Empty Or Sparse

Check that the reference and assembly FASTA IDs match the alignment file exactly. If you applied strict filters upstream, try plotting the unfiltered or less strictly filtered coords/PAF file first, then add `--min-segment-bp`, `--min-segment-idy`, or `--min-mapq` only after confirming the expected rows are present.

If the plot is not empty but the pattern is hard to classify, use
[How to Interpret Dot Plots]({{ '/dot-plots/' | relative_url }}) to compare it
against common clean, reversed, chimeric, inversion, duplicate, gap, and
repeat-like examples.

If the plot is empty because rows do not connect to the FASTA records, use
[FASTA and evidence name matching]({{ '/guides/name-matching/' | relative_url }}).
