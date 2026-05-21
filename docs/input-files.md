---
title: Input Files
description: Prepare MUMmer, minimap2, GFA, GAF, and Hi-C-like inputs for ChromoSort.
---

# Input Files

## Creating Input Files With MUMmer

ChromoSort uses `show-coords` output, ideally generated from a filtered nucmer
delta. The commands below are general defaults; tune them for genome size,
repeat content, assembly quality, and evolutionary distance.

```bash
mkdir -p mummer

ref=reference.fa
asm=assembly.fa
name=sample

# 1. Whole-assembly alignment.
nucmer \
  -t 16 \
  -c 500 \
  -p "mummer/${name}" \
  "$ref" \
  "$asm"

# 2. Keep one best alignment chain per query/reference region.
delta-filter \
  -i 95 \
  -l 10000 \
  -1 \
  "mummer/${name}.delta" \
  > "mummer/${name}.filter"

# 3. Export coordinates used by ChromoSort.
show-coords \
  -r \
  -c \
  -l \
  "mummer/${name}.filter" \
  > "mummer/${name}.coords"

# Optional visual inspection from the existing coords file.
chromo plot \
  --ref-fasta "$ref" \
  --assembly-fasta "$asm" \
  --coords "mummer/${name}.coords" \
  --output-prefix "mummer/plot_${name}" \
  --per-ref
```

### Why These MUMmer Choices?

`nucmer` aligns the reference and assembly at whole-genome scale. The `-c`
minimum cluster length removes very small seeds that are often unhelpful for
chromosome-scale contig ordering and chimeric-contig splitting.

`delta-filter -1` is used because these workflows usually want a primary
placement for each contig segment rather than every local repeat hit. It reduces
redundant alignments before ChromoSort applies interval merging and transition
detection.

`delta-filter -i` and `-l` enforce minimum identity and alignment length before
reporting. Use stricter values for very similar assemblies, and looser values
for distant species or more fragmented assemblies.

`show-coords -r -c -l` reports reference coordinates, query coordinates, lengths,
coverage, percent identity, and sequence names. ChromoSort reads those fields
and recomputes merged coverage itself.

## Creating Input Files With minimap2

ChromoSort can use minimap2 PAF directly. Choose the minimap2 preset that
matches the expected divergence between the reference and assembly.

```bash
mkdir -p paf

ref=reference.fa
asm=assembly.fa
name=sample

minimap2 \
  -x asm5 \
  -t 16 \
  --secondary=no \
  "$ref" \
  "$asm" \
  > "paf/${name}.paf"

chromo sort \
  --ref-fasta "$ref" \
  --assembly-fasta "$asm" \
  --paf "paf/${name}.paf" \
  --output-prefix "results/${name}" \
  --orient-to-reference

chromo plot \
  --ref-fasta "$ref" \
  --assembly-fasta "$asm" \
  --paf "paf/${name}.paf" \
  --assignments "results/${name}.contig_assignments.tsv" \
  --output-prefix "plots/${name}" \
  --per-ref
```

`--coords` and `--paf` are mutually exclusive for `chromo sort`,
`chromo fix`, and `chromo plot`. For PAF input, ChromoSort computes percent
identity from the PAF match and block-length columns, uses the PAF strand for
orientation, and skips rows marked `tp:A:S` unless `--include-secondary-paf` is
set. Use `--min-mapq` to ignore low-MAPQ PAF rows.

## Graph Input Files

Graph-aware ChromoSort commands use these graph-related evidence files:

- GFA: the assembly graph, used by `chromo sort --gfa`, `chromo manual --gfa`,
  `chromo fix --gfa`, `chromo scaffold --gfa`, and `chromo gapfill --gfa`.
- reference-to-assembly PAF: the same minimap2 alignment format used by
  `chromo sort`, `chromo fix`, `chromo manual`, `chromo gapfill --ref-paf`, and
  `chromo plot`.
- GAF: optional read-to-graph alignments used by `chromo gapfill --gaf` to resolve
  otherwise ambiguous graph paths.
- Hi-C pairs: optional graph-node contact counts used by
  `chromo gapfill --hic-pairs` as an additional conservative branch-support
  signal.

### Where to Find the GFA

The GFA usually comes from the assembler, not from ChromoSort. Look in the
original assembly output directory before any post-processing or renaming step.
Common examples are hifiasm primary/haplotype graph files, Verkko graph files,
or graph outputs produced while converting unitig/contig graphs to FASTA. The
most important practical rule is that GFA segment names must still match the
sequence names ChromoSort sees in the assembly FASTA or in the `chromo sort`
assignment report. If the FASTA was exported from the same graph, this usually
works naturally. If the FASTA was renamed, polished, split, or scaffolded by
another tool, keep a name map or regenerate graph evidence for the renamed
sequences.

For graph review, use the graph closest to the FASTA being sorted or filled:

```text
assembly.fa              # FASTA passed to chromo sort/fix/gapfill
assembly_graph.gfa       # GFA whose S records match assembly.fa sequence IDs
```

ChromoSort currently reads GFA `S` segment records and `L` link records. Segment
sequences are required only when `chromo gapfill --apply` may insert graph
sequence. Report-only graph evidence can still use segments with `*` sequence
fields when lengths are provided with `LN:i`.

### Which PAF Files to Keep

The main PAF file for ChromoSort is a reference-to-assembly whole-genome
alignment:

```bash
minimap2 \
  -x asm5 \
  -t 16 \
  --secondary=no \
  reference.fa \
  assembly.fa \
  > paf/sample.ref_vs_asm.paf
```

Use this PAF anywhere you would otherwise use MUMmer coords:

```bash
chromo sort --ref-fasta reference.fa --assembly-fasta assembly.fa \
  --paf paf/sample.ref_vs_asm.paf --output-prefix results/sample

chromo fix --assembly-fasta assembly.fa --paf paf/sample.ref_vs_asm.paf \
  --contigs suspect_1 suspect_2 --output-fasta results/sample.fixed.fa \
  --report results/sample.fixed.tsv

chromo manual --ref-fasta reference.fa --assembly-fasta assembly.fa \
  --paf paf/sample.ref_vs_asm.paf --output-html results/sample.manual.html

chromo plot --ref-fasta reference.fa --assembly-fasta assembly.fa \
  --paf paf/sample.ref_vs_asm.paf --output-prefix plots/sample

chromo gapfill --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --gfa assembly_graph.gfa --ref-paf paf/sample.ref_vs_asm.paf \
  --output-prefix results/sample.gapfill
```

If you run `chromo fix` and create a new fixed FASTA, re-align the fixed FASTA
and keep that second PAF beside the fixed results:

```bash
minimap2 -x asm5 -t 16 --secondary=no reference.fa results/sample.fixed.fa \
  > paf/sample.fixed.ref_vs_asm.paf
```

The original PAF explains the raw assembly; the fixed PAF explains the edited
assembly. Do not mix a fixed FASTA with an old raw-assembly PAF unless you are
only doing a very specific manual comparison. For `chromo gapfill --ref-paf`, use
the PAF whose query names still match the GFA segment names being evaluated.

### Creating GAF Read-to-Graph Alignments

GAF is a graph-alignment format. ChromoSort uses it in `chromo gapfill --gaf` as
optional read-path evidence for candidate graph fills. A typical source is a
long-read-to-GFA alignment from a graph aligner:

```bash
GraphAligner \
  -g assembly_graph.gfa \
  -f reads.fastq.gz \
  -a graph_alignments/sample.reads_to_graph.gaf
```

`GraphAligner` is an optional external tool; it is not needed for the core
sorting/fixing/scaffolding workflow. In the current `chromo gapfill` implementation,
GAF is used as path support only. ChromoSort reads the query name, path string,
and MAPQ columns, filters with `--min-gaf-mapq`, and counts how many read paths
contain each candidate graph path. If one candidate path has unique support
above `--min-gaf-path-support`, that path can resolve an otherwise ambiguous
GFA branch. Tied, weak, or absent support keeps the fill unresolved for review.

### Optional Hi-C Pair Evidence

`chromo gapfill --hic-pairs` accepts a simple tab-delimited graph-node contact
table:

```text
node_a  node_b  count
left  bridge_good  25
bridge_good  right  22
left  bridge_alt  3
bridge_alt  right  2
```

The first row may be a header. Node names must match GFA segment names.
ChromoSort treats contacts as undirected and scores a candidate fill path by
summing contacts between adjacent graph nodes along that path. An ambiguous
branch can be resolved by Hi-C support only when one candidate has unique summed
support at or above `--min-hic-path-support`. If GAF and Hi-C both uniquely
support different paths, ChromoSort keeps the junction unresolved for manual
review instead of choosing between conflicting evidence.
