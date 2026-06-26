---
title: hifiasm Unitig-To-Contig Projection
description: Project hifiasm unitig GFA coordinates onto contig FASTA coordinates for plot overlays and breakpoint review.
---

# hifiasm Unitig-To-Contig Projection

Use this guide when your dot plots are built from hifiasm contig FASTA records,
but your graph evidence is in hifiasm unitig GFA coordinates.

The main question is:

> Can this unitig graph evidence be compared to this contig-coordinate plot?

## The Core Idea

hifiasm often gives you both contig FASTA and unitig graph files. A coordinate
on `assembly.p_ctg.fa` is not the same coordinate system as a coordinate on
`assembly.p_utg.gfa`. ChromoSort can bridge them only when the GFA tells how
unitigs are arranged along contig paths or walks.

That bridge comes from GFA `P` path records or `W` walk records whose names
match the contig FASTA records.

```text
contig FASTA record
  -> matching GFA path or walk
  -> ordered unitig steps
  -> projected contig-coordinate intervals
```

## When Projection Is Needed

Projection is needed when:

- the plotted FASTA is `p_ctg.fa`, `hap*.p_ctg.fa`, or another contig FASTA,
- the graph is `p_utg.gfa`, `r_utg.gfa`, or `.noseq.gfa`,
- graph features are unitig-local,
- you want to compare graph boundaries to dot-plot breakpoints or contig
  coordinates.

Projection is not needed when GFA `S` segment names directly match the FASTA
records and graph evidence is being used only as node-level context.

## Run Graph Projection

```bash
chromo graph-map \
  --ctg-fasta assembly.p_ctg.fa \
  --utg-gfa assembly.p_utg.noseq.gfa \
  --output-prefix results/sample.graphmap
```

This writes:

| Output | What to read |
| --- | --- |
| `<prefix>.utg_to_ctg.tsv` | Projected unitig steps with contig coordinates, orientation, segment length, overlap fields, and reuse flags. |
| `<prefix>.path_summary.tsv` | One row per requested contig/path with projected bp, FASTA length, length-difference status, and warning counts. |
| `<prefix>.warnings.tsv` | Missing paths, missing segments, no-sequence/no-length records, and path-vs-FASTA length mismatches. |

By default, graph-map reports overlaps but does not subtract them from
projected spans. Add `--trim-overlaps` only when that is the projection model
you want to inspect.

## Draw Unitig Boundaries On Dot Plots

For visual review, use the same projection idea through `chromo plot`:

```bash
chromo plot \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.p_ctg.fa \
  --paf paf/sample.ref_vs_p_ctg.paf \
  --gfa-overlay assembly.p_utg.noseq.gfa \
  --gfa-overlay-mode unitig-boundaries \
  --output-prefix plots/sample.with_graph \
  --formats svg pdf
```

The overlay is drawn on the query axis. If the GFA lacks matching path/walk
records and segment names do not directly match the query FASTA, ChromoSort
writes a warning and an empty overlay report instead of guessing.

## Read The Projection Tables

### Projection Rows

`utg_to_ctg.tsv` answers:

| Field family | Why it matters |
| --- | --- |
| Contig/path name | Confirms which FASTA record was projected. |
| Contig start/end | Places the unitig interval on the plot coordinate system. |
| Unitig name and orientation | Identifies the GFA segment and direction. |
| Segment length | Confirms whether length came from sequence or `LN:i`. |
| Overlap fields | Shows GFA overlaps between adjacent path steps. |
| Reuse flag | Highlights repeated unitig usage in paths or walks. |

### Summary Rows

`path_summary.tsv` is the first place to look for projection health:

| Pattern | Interpretation |
| --- | --- |
| Projected length close to FASTA length | Projection is probably usable for review. |
| Large length difference | Path/walk and FASTA may not describe the same sequence. |
| Missing path | GFA does not contain a path/walk named for that FASTA contig. |
| Zero-length steps | Segment has no sequence and no `LN:i` length. |
| Reused segments | The path/walk repeats unitigs; inspect before using boundaries as simple evidence. |

### Warning Rows

Treat warnings as review prompts, not noise. A missing segment, missing length,
or path/FASTA mismatch can make a projected breakpoint misleading.

## What `.noseq.gfa` Can Do

`.noseq.gfa` is often ideal for projection and review because it keeps:

- segment names,
- segment lengths through `LN:i`,
- links,
- paths and walks,
- topology tags.

It does not keep nucleotide sequence. That is fine for boundary overlays and
unitig-neighborhood review. It is not enough for graph gapfill application,
where the graph must provide sequence to insert.

## Breakpoint Review Pattern

A useful graph-aware breakpoint review sequence is:

```text
assembly.p_ctg.fa + reference PAF
  -> chromo plot with GFA overlay
  -> inspect whether dot-plot breakpoint falls near projected unitig boundary
  -> chromo eval fix --gfa for graph_unitig fields
  -> manual or spreadsheet review
```

An alignment transition at a real unitig boundary is different from a transition
through the middle of a simple, well-supported unitig. It is still not proof by
itself. Use read evidence, alignment evidence, and graph context together.

## Cheat Sheet

| If you see... | Think... |
| --- | --- |
| Matching `P` or `W` path names | Projection can map unitigs to contig coordinates. |
| Only `S` and `L` records | Topology exists, but contig-coordinate projection may not. |
| `S node * LN:i:1000` | noseq segment has usable length. |
| `S node *` without length | Coordinates cannot be projected reliably for that segment. |
| Large path-vs-FASTA mismatch | Do not trust boundary locations without further review. |
| Reused unitig step | Inspect whether the graph path is repeat-like or circular. |

## Common Traps

Do not compare unitig-local GFA positions directly to contig FASTA positions.

Do not assume every hifiasm GFA has `P` or `W` records. Some graph files are
topology-only for this purpose.

Do not use `.noseq.gfa` for sequence insertion. It can support projection and
topology review, not gapfill sequence.

Do not ignore projection warnings just because a plot overlay was generated.
The warnings explain where the coordinate bridge is weak.

Do not treat a unitig boundary as proof of a misjoin. It is context for review.

## What To Look At Next In ChromoSort

- Use [Assembly Graph Evidence]({{ '/guides/assembly-graph-evidence/' | relative_url }})
  for GFA record types and graph guardrails.
- Use [chromo graph-map]({{ '/commands/graph-map/' | relative_url }}) for the
  projection command reference.
- Use [chromo plot]({{ '/commands/plot/' | relative_url }}) for graph overlay
  options.
- Use [Chimeric Contig And Breakpoint Review]({{ '/guides/breakpoint-review/' | relative_url }})
  when projected boundaries line up with suspicious alignment transitions.
