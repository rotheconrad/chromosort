---
title: chromo graph-map
description: Usage, outputs, parameters, and reasoning for chromo graph-map.
---

# chromo graph-map

Use `chromo graph-map` when graph evidence is in unitig coordinates but your
FASTA, alignments, and dot plots are in contig coordinates. This is common with
hifiasm outputs where dot plots are made from `p_ctg.fa` or `hap*.p_ctg.fa`,
while graph evidence comes from `p_utg.gfa`, `r_utg.gfa`, or their
`.noseq.gfa` equivalents.

Unitig and contig coordinates are not directly comparable. A `p_ctg.fa`
coordinate is local to the contig sequence. A `p_utg.gfa` coordinate is local to
one unitig segment. ChromoSort projects unitig intervals onto contig coordinates
only through GFA `P` path records or `W` walk records whose path names match
FASTA contig names.

## Run `chromo graph-map`

```bash
chromo graph-map \
  --ctg-fasta assembly.p_ctg.fa \
  --utg-gfa assembly.p_utg.noseq.gfa \
  --output-prefix results/sample.graphmap
```

This writes:

| Output | Description |
| --- | --- |
| `<prefix>.utg_to_ctg.tsv` | One projected interval per path step, including contig coordinates, unitig name, orientation, segment length, overlap fields, source GFA, and reuse flags. |
| `<prefix>.path_summary.tsv` | One row per requested contig/path with projected length, FASTA length, length-difference status, zero-length step count, and reuse count. |
| `<prefix>.warnings.tsv` | Clear warnings for missing paths, missing segments, no-sequence/no-`LN:i` segments, and path-vs-FASTA length mismatches. |

By default, ChromoSort reports GFA overlaps but does not subtract them from
projected spans. Add `--trim-overlaps` only when you intentionally want each
left overlap removed from the following path step.

## hifiasm noseq GFA

For graph topology, unitig boundaries, and junction context, `.noseq.gfa` is
usually preferred because it preserves segment names, lengths, links,
paths/walks, and tags without embedding large sequence strings. Use full `.gfa`
or FASTA only when nucleotide sequence is needed.

If a noseq segment has `S	utgX	*	LN:i:1000`, ChromoSort can project it as a
1000 bp segment. If the segment has `S	utgX	*` and no `LN:i`, ChromoSort keeps
the length as 0 and writes:

```text
Segment utgX has no sequence and no LN:i length; cannot project coordinates reliably.
```

Some hifiasm GFA files contain only `S`, `L`, and read-alignment `A` records.
Those files are still useful as graph topology inputs, but they do not contain
the path/walk records needed to map unitig coordinates onto contig FASTA
coordinates. In that case `chromo graph-map` writes a warning instead of
pretending the coordinate systems are interchangeable.

## Parameters

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `--ctg-fasta` | required | Contig FASTA whose record names should match GFA path/walk names. |
| `--ctg-fai` | none | Optional FASTA index; defaults to `<ctg-fasta>.fai` when present. |
| `--utg-gfa` | required | Unitig or contig GFA. `.noseq.gfa` is preferred when sequence is not needed. |
| `--path-name` | all FASTA contigs | Specific GFA path/walk name to project; may be repeated. |
| `--output-prefix` | required | Prefix for projection, path-summary, and warning TSVs. |
| `--trim-overlaps` | off | Subtract each left overlap from the following path-step span. |
| `--length-tolerance-bp` | `1000` | Absolute path-vs-FASTA length tolerance. |
| `--length-tolerance-frac` | `0.01` | Fractional path-vs-FASTA length tolerance. |
