---
title: chromo gafprep
description: Usage, outputs, parameters, and reasoning for targeted GraphAligner GAF preparation.
---

# chromo gafprep

`chromo gafprep` prepares targeted GraphAligner inputs from existing
long-read-to-assembly PAF alignments and ChromoSort review tables. It does not
run GraphAligner and it does not create the GAF itself. Instead, it selects
review-relevant reads, writes a smaller FASTQ, sanitizes a GFA for
GraphAligner, and writes a runnable shell script.

Use it when full-depth read-to-graph alignment would be too expensive, but you
only need GAF evidence near candidate breakpoints, scaffold junctions, gapfill
paths, graph branches, or contig ends.

For broad review, first run `chromo eval all` on the same FASTA naming stage
used by the read-to-assembly PAF. Then pass the emitted fix, scaffold, and
gapfill review tables to `chromo gafprep` so one selected FASTQ can support all
three downstream review paths.

## Run `chromo gafprep`

```bash
chromo eval all \
  --assembly-fasta results/sample.ordered.fa \
  --coords mummer/ordered.coords \
  --all \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --gfa assembly_graph.gfa \
  --read-paf reads_to_ordered.paf \
  --output-prefix review/sample.eval

chromo gafprep \
  --assembly-fasta results/sample.ordered.fa \
  --assembly-gfa assembly_graph.gfa \
  --read-paf reads_to_ordered.paf \
  --reads reads.fastq.gz \
  --eval-review-table review/sample.eval.fix_review.tsv \
  --eval-review-table review/sample.eval.scaffold_review.tsv \
  --eval-review-table review/sample.eval.gapfill_review.tsv \
  --output-prefix results/sample.gafprep \
  --seed 1
```

Then run the generated script on an interactive node or scheduler job:

```bash
bash results/sample.gafprep.graphaligner.sh
```

The script writes:

```text
results/sample.gafprep.gaf
```

That GAF can be passed back into `chromo eval fix`, `chromo eval scaffold`,
`chromo eval gapfill`, `chromo manual`, or `chromo gapfill --gaf`.

## Inputs

`--read-paf` must be a read-to-assembly PAF, with reads as queries and assembly
contigs as targets:

```bash
minimap2 -x map-hifi -c -t 16 --secondary=no assembly.fa reads.fastq.gz \
  > reads_to_assembly.paf
```

Do not pass a reference-to-assembly PAF here. Reference PAF supports sorting,
plotting, and fixing; read PAF supports read selection around assembly
coordinates.

`--eval-review-table` accepts one or more ChromoSort review TSVs from
`chromo eval fix`, `chromo eval scaffold`, `chromo eval gapfill`, or the three
tables emitted by `chromo eval all`. The command infers the review type from
the `task` column, filename, or table columns. Use `--review-type fix`,
`--review-type scaffold`, or `--review-type gapfill` only when a custom table
cannot be inferred.

The review tables, `--assembly-fasta`, and `--read-paf` should refer to the
same assembly naming stage. For example, if scaffold and gapfill rows name
records from `results/sample.ordered.fa`, align reads to that ordered FASTA and
pass that same FASTA to `chromo gafprep`.

## Outputs

For prefix `results/sample.gafprep`, the command writes:

| Output | Meaning |
| --- | --- |
| `results/sample.gafprep.targets.tsv` | Target intervals derived from review rows. |
| `results/sample.gafprep.selected_reads.tsv` | One row per selected read with best assembly alignment, reasons, target IDs, and review row IDs. |
| `results/sample.gafprep.selected_read_review_links.tsv` | Read-to-review-row links preserving target, event, distance, and overlap provenance. |
| `results/sample.gafprep.selected_read_ids.txt` | Selected read IDs, one per line. |
| `results/sample.gafprep.selected.fastq.gz` | FASTQ subset streamed from the original reads. |
| `results/sample.gafprep.graphaligner.gfa` | GFA copy sanitized for GraphAligner. |
| `results/sample.gafprep.graphaligner.sh` | Executable GraphAligner shell script. |
| `results/sample.gafprep.summary.tsv` | Run counters, selected-read counts, PAF filter counts, and warnings. |
| `results/sample.gafprep.dropped_gfa_links.tsv` | Audit rows for pathological full-consuming GFA links removed from the GraphAligner GFA. |
| `results/sample.gafprep.gfa_sanitize_summary.tsv` | GFA sanitization counters and warnings. |

The link table is the key audit file when interpreting targeted GAF evidence:
it explains which selected reads came from which review rows.

## Target Selection

Fix review rows generate windows around candidate breakpoint positions, planner
breakpoints, split-piece boundaries, or alignment-query boundaries when those
columns are present. If exact coordinates are unavailable, `chromo gafprep`
falls back to both contig ends.

Scaffold and gapfill review rows generate windows at the right end of the left
flank and the left end of the right flank. This captures reads that may bridge
or distinguish adjacent contigs, graph paths, and alternate branches.

Default windows are intentionally broad:

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `--target-padding` | `50000` | Padding around exact breakpoint-like coordinates. |
| `--contig-end-window` | `100000` | Contig-end fallback and junction flank window. |
| `--target-reads-per-interval` | `20` | Best reads retained per target interval. |
| `--min-mapq` | `20` | Minimum read-to-assembly PAF MAPQ. |
| `--min-aligned-bp` | `5000` | Minimum PAF alignment block length. |
| `--background-bin-size` | `1000000` | Assembly bin size for background sampling. |
| `--background-reads-per-bin` | `1` | Deterministic background reads retained per bin. |
| `--seed` | `1` | Seed for deterministic background sampling. |

Optional `--max-reads` and `--max-reads-per-contig` cap the final selected read
set after target and background reasons have been assigned.

## GFA Sanitization

The GraphAligner GFA preparation is conservative:

- `S` segment records and sequence strings are kept unchanged.
- Node IDs are kept unchanged.
- Ordinary `L` links are kept unchanged.
- `H`, `P`, and `W` records are kept.
- hifiasm-style `A` records are removed by default.
- A simple `L` overlap is dropped only when its overlap length fully consumes
  either endpoint node.
- Sequences are never trimmed and replacement links are never invented.

Every dropped link is written to `.dropped_gfa_links.tsv` with node names,
orientations, node lengths, CIGAR, overlap length, and reason. If a dropped link
touches a selected target contig, summaries include:

```text
GAF evidence limited by graph sanitization near target
```

Treat targeted GAF from that region as limited evidence.

## GraphAligner Script Options

`chromo gafprep` writes a script like:

```bash
GraphAligner \
  -g results/sample.gafprep.graphaligner.gfa \
  -f results/sample.gafprep.selected.fastq.gz \
  -a results/sample.gafprep.gaf \
  -t 16 \
  -x vg \
  --precise-clipping 0.9
```

Use these options to customize the generated script:

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--graphaligner-bin` | `GraphAligner` | Executable name or path. |
| `--graphaligner-threads` | `16` | Thread count written to `-t`. |
| `--graphaligner-preset` | `vg` | Preset written to `-x`. |
| `--precise-clipping` | `0.9` | GraphAligner clipping setting. |

GraphAligner is intentionally not a Python dependency of ChromoSort.

## ctg.gfa Or utg.gfa?

Use the graph that matches your review goal.

For contig-level node context and direct contig junctions, a `ctg.gfa` whose
segment names match the assembly FASTA is usually easiest. It lets selected
reads align to graph nodes that review tables can name directly.

For unitig-level branch evidence, a `utg.gfa` or `utg.noseq.gfa` can be useful,
especially when hifiasm paths or walks connect unitigs into contigs. Remember
that unitig coordinates are not contig coordinates. If you need to compare
unitig graph context to contig breakpoints, use `chromo graph-map` or a GFA
with matching `P`/`W` paths to understand the projection first.

Full-sequence GFA is needed for GraphAligner. `.noseq.gfa` files can be useful
for topology review, but they cannot produce meaningful read-to-sequence graph
alignments unless the graph aligner has access to sequence through another
supported route.

## Reasoning

Targeted GAF is review evidence, not automatic sequence validation. It helps
answer whether selected reads traverse a candidate graph path or alternate
branch. It does not prove that a fill sequence is valid, and it does not change
FASTA records by itself.
