---
title: chromo plot
description: Usage, outputs, parameters, and reasoning for chromo plot.
---

# chromo plot

Use `chromo plot` when you already have a MUMmer coords or minimap2 PAF file
and want a visual check without running `mummerplot` or re-aligning only for
the plot. It writes PDF by default and can also write SVG or PNG. Forward-strand
alignments are blue and reverse-strand alignments are red.

<figure>
  <img src="https://rotheconrad.github.io/chromosort/assets/chromo_plot_example.png" alt="ChromoSort whole-genome dot plot with forward and reverse alignments." width="900">
  <figcaption><strong>Figure 2. chromo plot whole-genome view.</strong> Example PNG output from the synthetic coords fixture, showing all reference sequences on the x-axis and assembly contigs on the y-axis.</figcaption>
</figure>

<figure>
  <img src="https://rotheconrad.github.io/chromosort/assets/chromo_plot_example.chr1.png" alt="ChromoSort per-reference dot plot for chr1." width="900">
  <figcaption><strong>Figure 3. chromo plot per-reference view.</strong> Example `--per-ref` output for `chr1`, useful for inspecting a chromosome-level slice without re-running an aligner.</figcaption>
</figure>

## Run `chromo plot`

Whole-genome plot from MUMmer coords:

```bash
chromo plot \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.fa \
  --coords mummer/sample.coords \
  --output-prefix plots/sample
```

This writes `plots/sample.pdf`.

Whole-genome and per-reference plots from PAF, ordered by a `chromo sort`
assignment report:

```bash
chromo plot \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.fa \
  --paf paf/sample.paf \
  --assignments results/sample.contig_assignments.tsv \
  --output-prefix plots/sample \
  --per-ref \
  --formats pdf svg
```

When `--assignments` is provided, the query axis is ordered by the kept contigs
in the assignment report. This is useful for reviewing a sorted ChromoSort
order from the same alignment file that powered the sort.

## `chromo plot` Outputs

| Output | Description |
| --- | --- |
| `<prefix>.pdf` | Whole-genome PDF dot plot by default. |
| `<prefix>.svg` | Whole-genome SVG dot plot when `--formats svg` is set. |
| `<prefix>.png` | Whole-genome PNG dot plot when `--formats png` is set. |
| `<prefix>.<ref>.<format>` | Per-reference plots when `--per-ref` is set. |

## `chromo plot` Parameters

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `--coords` | required unless `--paf` | MUMmer `show-coords` alignment file. |
| `--paf` | required unless `--coords` | minimap2 PAF alignment file. |
| `--formats` | `pdf` | One or more output formats: `pdf`, `svg`, `png`. |
| `--assignments` | none | Optional `chromo sort` assignment report for ordering the query axis by kept sorted contigs. |
| `--per-ref` | off | Also write one plot per reference sequence with plotted alignments. |
| `--per-ref-query-order` | `fasta` | Use FASTA order or first reference-hit order for per-reference query axes. |
| `--min-segment-bp` | `0` | Minimum query-aligned bp for a row to be drawn. |
| `--min-segment-idy` | `0.0` | Ignore individual alignment rows below this percent identity. |
| `--min-mapq` | `0` | Ignore PAF rows below this MAPQ. Ignored for coords. |
| `--include-secondary-paf` | off | Include PAF rows marked `tp:A:S`; skipped by default. |
| `--max-segments` | `0` | Maximum drawn alignment rows after filtering; 0 means no limit. |
