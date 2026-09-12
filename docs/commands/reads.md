---
title: chromo reads
---

# `chromo reads`

Render a native-coordinate read evidence window, optionally reversed for
comparison, from a digest-checked [input manifest]({{ '/labeled-references/' | relative_url }}).
SVG, PDF and PNG share the same drawing model. No remapping or variant calling
is performed by the rendering command.

```bash
chromo reads --manifest inputs.json --evidence-id reads_hifi \
  --contig contig_1 --start 90000 --end 110000 --breakpoint 100000 \
  --event-id decision_example --action review --min-mapq 20 \
  --min-anchor-bp 1000 --read-window-bp 2000 --bin-bp 100 \
  -o figures/breakpoint
chromo reads replay --figure-json figures/breakpoint.figure.json \
  -o figures/replayed
```

Window intervals are 0-based half-open. `--breakpoint N` marks the boundary
after N native bases; it must lie strictly inside the window. `--reverse`
reverses display coordinates and arrows while preserving native labels and
native-strand colors. `--action inspect` labels an inspection boundary rather
than a proposed cut. Omit `--breakpoint` for an overview.

PAF, SAM and indexed BAM are supported. BAM uses `samtools view` and requires a
matching index; use the standard adjacent `.bai` name. CRAM is deliberately
unsupported until a reference-aware adapter is validated. No new Python plotting
dependency is required: SVG/PDF use ChromoSort's existing writers, PNG uses
Pillow. Missing read evidence produces an explicitly empty panel. Select
`--evidence-id` if several read sources are present.

## Evidence semantics

- Secondary alignments are excluded by default; supplementary SAM/BAM segments
  are retained and drawn as dashed bars. Segments of one read share a packed
  row and are connected. SA tags are annotations, never invented extra rows.
- Coverage is mean unique-read-name coverage per bin using aligned CIGAR blocks
  when available; D/N gaps contribute no aligned coverage. Without CIGAR it is
  explicitly alignment-span coverage, not a base pileup.
- Spanning counts require one alignment to cover both anchors. When CIGAR is
  available, an aligned block must bridge the anchors; a deletion crossing the
  boundary is not counted as continuity. Split-across-boundary counts require
  distinct, nonoverlapping segments of the same read on the two sides, with
  nonoverlapping read coordinates. Inversion partners may map elsewhere and
  need supplementary/clipping inspection even when this local split count is 0.
- Terminal clips use the original SAM CIGAR endpoints, not viewport edges.
  `--min-clip-bp` (default 100) controls the near-boundary clipping count.
  PAF does not identify SAM supplementary flags or terminal clipping; both are
  shown as unavailable. A missing SAM NM tag means identity unavailable.
- MAPQ 255 is unavailable and fails a positive MAPQ threshold. MAPQ from
  different reference alignments is not a comparable posterior probability.
- `*.figure.json` records filter counts, unique molecules versus alignments,
  displayed counts, exact inputs, stage/sample/evidence role and settings.
  `*.reads.tsv` lists every displayed segment. Spanning/split reads receive
  selection priority, followed by deterministic read-name hashes; selected
  molecules are packed in display-coordinate order. Limits affect the read
  drawing, never the full filtered coverage or evidence summaries.

Orientation balance, clipped alignments and coverage differences do not by
themselves establish an assembly error. Same-sample and cross-sample read
evidence must be labeled separately. BAM region access retrieves overlapping
alignments only; partners outside the region are not silently reconstructed.

## Features and variants

Feature TSV columns: `contig`, `start`, `end`, `label` (0-based half-open).
Variant TSV adds `kind` (`SNP`, `INDEL`, `SV`). A `variants` entry ending in
`.vcf` or `.vcf.gz` is also accepted; each supplied ALT allele is drawn in its
class lane, retaining the supplied FILTER value in the adapter. Filtering or
calling variants is an upstream responsibility; no PASS-only selection is
silently added. Empty/missing overlays are labeled unavailable.

The [workflow dashboard]({{ '/commands/workflow/' | relative_url }}) links each event to its read SVG/PDF.
Replay requires the saved ChromoSort version and the same input digests;
`--manifest` on `reads replay` supports a relocated, unchanged manifest.
