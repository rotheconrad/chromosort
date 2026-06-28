---
title: Input Files
description: Prepare MUMmer, minimap2, GFA, GAF, Hi-C-like, and patch-candidate inputs for ChromoSort.
---

# Input Files

ChromoSort does not run aligners or graph tools for you. It consumes files from
MUMmer, minimap2, assembly graph outputs, graph aligners, optional contact
tables, and optional patch-candidate workflows. The most important rule is that
each evidence file must describe the same FASTA records used by the command you
are running.

## Input Sets By Task

| Goal | Required inputs | Optional inputs |
| --- | --- | --- |
| Sort contigs with `chromo sort` | Reference FASTA, assembly FASTA, and either MUMmer coords or minimap2 PAF | GFA for report-only graph context, FASTA indexes |
| Clean mostly-correct assemblies with `chromo clean` | Reference FASTA, raw assembly FASTA, and either MUMmer coords or minimap2 PAF from that raw assembly | FASTA indexes, optional fix-target list |
| Prepare fix review tables with `chromo eval fix` | Assembly FASTA and either MUMmer coords or minimap2 PAF | GFA for graph context, long-read-to-assembly PAF, long-read-to-graph GAF |
| Fix chimeric contigs with `chromo fix` | Assembly FASTA and either MUMmer coords or minimap2 PAF | GFA for report-only graph context |
| Cut reviewed coordinates with `chromo cut` | Assembly FASTA and explicit cut positions | Assembly FAI |
| Review manually with `chromo manual` | Reference FASTA, assembly FASTA, and either MUMmer coords or minimap2 PAF | GFA, long-read-to-assembly PAF, long-read-to-graph GAF, FASTA indexes, embedded sequences, `chromo eval` review table for task modes |
| Plot alignments with `chromo plot` | Reference FASTA, assembly FASTA, and either MUMmer coords or minimap2 PAF | Assignment TSV for ChromoSort query ordering, GFA overlay with matching path/walk records, FASTA indexes. See the [dot-plot guide]({{ '/dot-plots/' | relative_url }}) for interpretation examples. |
| Map graph coordinates with `chromo graph-map` | Contig FASTA and GFA with path/walk records | FASTA index, selected path names |
| Prepare scaffold review tables with `chromo eval scaffold` | Ordered FASTA and matching `chromo sort` assignment TSV | GFA for graph junction context, long-read-to-assembly PAF, long-read-to-graph GAF |
| Scaffold sorted contigs with `chromo scaffold` | Ordered FASTA and matching `chromo sort` assignment TSV | GFA for report-only graph junction evidence |
| Prepare gapfill review tables with `chromo eval gapfill` | Ordered FASTA plus matching assignment TSV, or scaffold FASTA plus AGP; GFA | `--project-gfa-paths` for unitig-level GFA `P`/`W` projection, GAF read paths, Hi-C-like graph-node contacts, reference-placement PAF, long-read-to-assembly PAF, external patch table/FASTA |
| Fill graph-supported gaps with `chromo gapfill` | Ordered FASTA plus matching assignment TSV, or scaffold FASTA plus AGP; GFA | `--project-gfa-paths` for unitig projection, GAF read paths, Hi-C-like graph-node contacts, reference-placement PAF, external patch table/FASTA, reviewed plan TSV |
| Prepare targeted GraphAligner inputs with `chromo gafprep` | Assembly FASTA, assembly GFA, long-read-to-assembly PAF, read FASTQ, and one or more ChromoSort review TSVs | Optional target/background sampling limits and GraphAligner script settings |

When a command accepts `--coords` or `--paf`, provide exactly one of them. For
most new runs, use minimap2 PAF as the primary alignment input because it is
fast, compact, and carries MAPQ; keep MUMmer coords as a good alternative when
you want a different aligner perspective. When a command takes an assignment
report, use the report written by the same `chromo sort` run as the ordered
FASTA.

## Name Matching Rules

Most confusing input-file problems come from sequence names drifting apart
between FASTA, alignments, reports, and graph files.

- Reference names in coords or PAF must match records in `--ref-fasta`.
- Query names in coords or PAF must match records in `--assembly-fasta`.
- If you run `chromo clean`, `chromo fix`, `chromo cut`, or `chromo manual apply`,
  re-align the new FASTA before running downstream alignment-dependent commands.
- GFA `S` record names should match the assembly FASTA records used for direct
  graph-node review. For scaffold and gapfill, ChromoSort also tries ChromoSort
  `new_name` values from the assignment report when resolving graph nodes.
- Unitig-level GFA coordinates are not contig FASTA coordinates. For hifiasm
  `p_utg.gfa`, `r_utg.gfa`, or `.noseq.gfa` evidence on a `p_ctg.fa` dot plot,
  use GFA `P` path or `W` walk records to project unitigs onto matching contig
  names with `chromo graph-map` or `chromo plot --gfa-overlay`.
- `chromo gafprep --read-paf` expects reads as PAF queries and assembly contigs
  as PAF targets. This is different from the reference-to-assembly PAF used by
  `chromo sort`, `chromo fix`, `chromo manual`, and `chromo plot`.
- `chromo gapfill --ref-paf` scores intermediate graph nodes, so the PAF query
  names must match GFA segment names for those nodes.
- Reviewed gapfill plans are tied to current `scaffold`, `left_contig`,
  `right_contig`, and `path_nodes` values. If the graph, ordered FASTA, or
  assignments change, regenerate the plan.

For commands that expose `--ref-fai` or `--assembly-fai`, FASTA indexes are
optional. ChromoSort uses `<fasta>.fai` when present and otherwise scans the
FASTA length directly. The index must describe the exact FASTA file being used.

## FASTA And Alignment Compatibility

A MUMmer coords file or minimap2 PAF file is tied to the exact reference FASTA
and assembly FASTA that produced it. The sequence IDs and coordinates in that
alignment do not automatically follow later FASTA edits.

You can reuse an alignment file for multiple decisions about the same assembly
FASTA. For example, if `raw.coords` was generated from `raw.fa`, it can support
both a first-pass `chromo sort --assembly-fasta raw.fa` and a later
`chromo fix --assembly-fasta raw.fa` on reviewed or automatically detected
contigs from that same raw assembly.

Re-run MUMmer or minimap2 before using a changed FASTA as the assembly input to
another alignment-dependent command. This applies after any step that removes,
splits, cuts, reverse-complements, renames, scaffolds, or manually exports
records, including:

- `chromo sort` when you want to operate on `<prefix>.ordered.fa`.
- `chromo clean` when you want to validate or continue from `<prefix>.clean.fa`.
- `chromo fix` when you want to operate on `fixed.fa`.
- `chromo cut` when you want to operate on the cut FASTA.
- `chromo manual apply` or browser FASTA export when you want to validate or
  continue with the manually edited FASTA.
- `chromo scaffold` when later evidence should describe scaffold records rather
  than the pre-scaffold contigs.

`chromo plot --assignments` is the main exception to "changed FASTA means
changed alignment" expectations. It still plots the original coords or PAF
rows, but uses a `chromo sort` assignment report to order the query axis by kept
sorted contigs. This is useful for reviewing sort decisions without re-running
an aligner. It is not a fresh alignment of `<prefix>.ordered.fa`. The
[dot-plot guide]({{ '/dot-plots/' | relative_url }}) shows how to interpret
the visual patterns that this review plot can reveal.

## File Format Contracts

### FASTA and FAI

FASTA record IDs are read from the first whitespace-delimited token after `>`.
Keep IDs unique. Optional `.fai` files speed length lookup and are used by
commands that expose `--ref-fai` or `--assembly-fai`. If an index is stale,
delete it or regenerate it before running ChromoSort.

### MUMmer coords

Use `show-coords` output from a reference-vs-assembly nucmer alignment. The
recommended export is:

```bash
show-coords -r -c -l sample.filter > sample.coords
```

ChromoSort reads reference coordinates, query coordinates, row lengths,
percent identity, reference length, query length, and sequence names. Coordinates
are normalized internally before interval merging.

MUMmer coords is a good primary input when a project already has a stable
nucmer workflow or when minimap2 PAF gives a surprising result that deserves a
second aligner view. It may produce a more fragmented row set than PAF on large
plant genomes, so `chromo fix` can take longer even when the final biological
interpretation is similar.

### minimap2 PAF

ChromoSort expects standard PAF rows with at least the first 12 columns. It uses
query name, query length, query start/end, strand, target name, target length,
target start/end, matching bases, block length, and MAPQ. Percent identity is
computed as matching bases divided by block length. Rows marked `tp:A:S` are
skipped unless `--include-secondary-paf` or `--include-secondary-ref-paf` is
set for the command reading that PAF.

For normal reference ordering, minimap2 should be run with the reference FASTA
as target and the assembly FASTA as query. This is the recommended starting
point for most ChromoSort production runs:

```bash
minimap2 -x asm5 -c -t 16 --secondary=no reference.fa assembly.fa \
  > paf/sample.ref_vs_asm.paf
```

`-c` asks minimap2 to perform base-level alignment and write CIGAR-bearing PAF
rows. ChromoSort does not parse the CIGAR string directly, but the base-level
alignment changes the PAF match and block-length columns that ChromoSort uses
for identity summaries and optional identity filters. Without `-c`, long
assembly alignments can have misleadingly low PAF column identity even when the
minimap2 divergence tags indicate a close match.

### GFA

ChromoSort reads GFA segment (`S`), link (`L`), path (`P`), and walk (`W`)
records. Unknown record types are ignored. Segment sequences may be `*` for
report-only graph context when `LN:i` length tags are present. Segment sequences
are required for `chromo gapfill --apply`, because the command must validate
flank sequences and construct inserted graph sequence.

Contig FASTA/GFA coordinates and unitig GFA coordinates are different coordinate
systems. A `p_ctg.fa` dot plot uses contig coordinates. A `p_utg.gfa` feature is
local to one unitig `S` segment. ChromoSort only overlays unitig features on a
contig dot plot after projecting them through matching GFA `P` path or `W` walk
records. If those path/walk records are absent, ChromoSort warns instead of
treating the coordinates as interchangeable.

For hifiasm graph topology, boundaries, and junction evidence, `.noseq.gfa` is
usually preferred when nucleotide sequence is not needed. It preserves segment
names, `LN:i` lengths, links, paths/walks, and tags while avoiding huge embedded
sequence strings. Use a full `.gfa` or FASTA only when sequence is required.

Only simple link overlap CIGARs made from `M`, `=`, and `X` operations are
treated as exact overlap lengths. Complex overlaps are preserved as unknown so
sequence-changing commands cannot use them as trim lengths.

### GAF

Commands with `--gaf` read graph-alignment rows with at least 12 columns. They
use the query name, path string, and MAPQ. The path string must encode oriented
graph nodes, for example `>left>bridge_good>right`. In `chromo eval fix`, GAF is
advisory node context. In `chromo eval scaffold`, `chromo eval gapfill`, and
`chromo gapfill`, it reports candidate graph traversal support. GAF support
does not insert sequence by itself; `chromo gapfill --apply` still requires a
validated GFA path with usable segment sequences and overlaps.

`chromo gafprep` is different: it prepares a smaller read FASTQ, a sanitized
GFA, and a GraphAligner shell script. GraphAligner still produces the actual
GAF file.

### Hi-C Pair Table

`chromo gapfill --hic-pairs` expects a tab-delimited table with graph node names
and non-negative integer contact counts:

```text
node_a  node_b  count
left  bridge_good  25
bridge_good  right  22
```

The first data row may be a header. Contacts are treated as undirected and are
summed across adjacent node pairs on each candidate graph path.

### External Patch Candidate Table

`chromo eval gapfill --patch-table` and `chromo gapfill --patch-table` import
candidate patches from external gap closers or patching workflows as review
evidence. Typical sources include TGS-GapCloser, LR_Gapcloser, Sealer, RagTag
patch, or a custom consensus workflow. ChromoSort compares each candidate patch
to the graph-derived fill sequence for the same gap, but it does not insert
external patch sequence in this mode.

The table is tab-delimited and keyed by scaffold object plus the two flanking
components:

```text
scaffold  left_contig  right_contig  patch_id  source  patch_sequence
chr1      ctgA         ctgB          patch_01  TGS-GapCloser  ACGTACGT
```

Required fields are `scaffold`, `left_contig`, and `right_contig`. Provide
either an inline `patch_sequence` column or a `patch_id` column plus
`--patch-fasta` containing matching FASTA records. Optional `source` and
`notes` columns are copied into the plan. Common aliases such as `object`,
`left_component`, `right_component`, `tool`, `method`, `sequence`, and `seq`
are also accepted.

Use the component names ChromoSort reports for the current gapfill mode. In
ordered-FASTA mode, these are the adjacent ordered contig IDs. In
scaffold/AGP mode, they are AGP component IDs around each N gap. A patch table
keyed only by scaffold record names cannot identify which N block it describes.

Patch evidence adds columns such as `patch_candidate_count`, `patch_best_id`,
`patch_best_source`, `patch_best_bp`, `patch_graph_status`, and
`patch_best_notes`. With `--include-fill-sequences`, the best external sequence
is also written as `patch_best_sequence`. Graph comparison statuses are
`exact_graph_match`, `same_length_graph_mismatch`, `graph_mismatch`, or
`patch_only_no_graph_sequence`.

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

# Optional: add --sel-ref Gm6 Gm12 Gm15 for focused replotting.
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

ChromoSort can use minimap2 PAF directly. Choose the strictest minimap2 preset
that still recovers the expected chromosome-scale alignments.

```bash
mkdir -p paf

ref=reference.fa
asm=assembly.fa
name=sample

minimap2 \
  -x asm5 \
  -c \
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

# Optional: add --sel-ref Gm6 Gm12 Gm15 for focused replotting.
```

`--coords` and `--paf` are mutually exclusive for `chromo sort`,
`chromo fix`, and `chromo plot`. For PAF input, ChromoSort computes percent
identity from the PAF match and block-length columns, uses the PAF strand for
orientation, and skips rows marked `tp:A:S` unless `--include-secondary-paf` is
set. `--secondary=no` keeps those skipped secondary rows out of the PAF file in
the first place, reducing file size and making downstream review less noisy. Use
`--min-mapq` to ignore low-MAPQ PAF rows.

PAF and coords are not expected to be byte-for-byte interchangeable. ChromoSort
normalizes both formats into the same internal alignment records, but minimap2
and MUMmer can differ in chaining, row fragmentation, secondary/primary
classification, MAPQ availability, and identity reporting. In a soybean
coords-vs-PAF `chromo fix` benchmark, split counts were within about 5-10%,
while the exact set of marginal split contigs differed by about 20-30% depending
on mode. Larger disagreements are a prompt to inspect plots, row counts, MAPQ,
secondary rows, and preset/filter choices before treating the event as biology.

### Choosing `asm5`, `asm10`, or `asm20`

minimap2's assembly presets tune seeding, chaining, and alignment scoring for
different assembly-to-reference distances. The
[minimap2 manual](https://lh3.github.io/minimap2/minimap2.html) describes
`asm5` as appropriate for very close average divergence, `asm10` for roughly
1% average divergence, and `asm20` for several-percent average divergence. The
[minimap2 cookbook](https://github.com/lh3/minimap2/blob/master/cookbook.md#full-genome-alignment)
uses broader wording for full-genome comparisons. Treat these descriptions as
starting points, not hard ANI cutoffs.

For ChromoSort, choose the most specific preset that still gives the expected
chromosome-scale alignments:

| Preset | Consider it when | ChromoSort-specific caution |
| --- | --- | --- |
| `-x asm5` | Same species, same breeding pool, cultivar/line comparisons, or a new assembly against a close reference. This is the safest default for reference-guided sorting when high identity is expected. | Too stringent for distant references: real contigs may appear fragmented or unaligned. If many expected contigs are `no_alignment` or break into short blocks, try `asm10`. |
| `-x asm10` | More divergent same-species material, wild or exotic accessions, pangenome references, or close relatives where `asm5` misses obvious syntenic blocks. | More permissive than `asm5`; ambiguous repeat or paralog matches may increase. Inspect dot plots, assignment reports, and best-reference shares. |
| `-x asm20` | Related species, highly divergent reference choices, or difficult graph-node/reference placement where `asm10` still misses expected chromosome-scale alignments. | Use with care. This is the noisiest of the three for ChromoSort because it can increase low-MAPQ rows, ambiguous placements, and misleading short matches. |

A practical progression is to start with `asm5` for same-species work, move to
`asm10` only if expected contigs are missing or highly fragmented, and reserve
`asm20` for clearly more divergent references. When using `asm10` or `asm20`,
keep `-c --secondary=no` in the minimap2 command, review `chromo plot` output,
and consider stricter ChromoSort filters such as `--min-mapq` and larger minimum
aligned-bp thresholds. Identity filtering can also be useful, but check the PAF
identity distribution from your chosen preset before setting `--min-segment-idy`
because PAF column identity is not always directly comparable to MUMmer coords
percent identity.

ChromoSort does not need minimap2's `--cs` tag or SAM output. Recommended PAF
comes from:

```bash
minimap2 -x asm5 -c --secondary=no reference.fa assembly.fa
```

Use `asm10` or `asm20` in place of `asm5` when the stricter preset misses
expected chromosome-scale alignments. ChromoSort reads the PAF coordinates,
lengths, strand, match count, block length, MAPQ, and names.

## Graph Input Files

Graph-aware ChromoSort commands use these graph-related evidence files:

- GFA: the assembly graph, used by `chromo sort --gfa`, `chromo manual --gfa`,
  `chromo eval fix/scaffold/gapfill/all`, `chromo fix --gfa`,
  `chromo scaffold --gfa`, `chromo gapfill --gfa`, `chromo graph-map`,
  `chromo gafprep --assembly-gfa`, and `chromo plot --gfa-overlay`.
- reference-to-assembly PAF: the minimap2 alignment format used by
  `chromo sort`, `chromo fix`, `chromo manual`, and `chromo plot`; for
  `chromo gapfill --ref-paf`, the PAF query names must match the GFA graph
  nodes being scored.
- long-read-to-assembly PAF: optional read alignments used by
  `chromo eval fix`, `chromo eval scaffold`, `chromo eval gapfill`,
  `chromo eval all`, and
  `chromo manual --read-paf` task dashboards to report split, bridge, and
  contig-end support. `chromo gapfill --read-paf` also reports bridge evidence
  across adjacent contig ends, but does not use reads as the inserted sequence.
- GAF: optional read-to-graph alignments used by `chromo eval fix`,
  `chromo eval scaffold`, `chromo eval gapfill`, `chromo manual --gaf`, and
  `chromo gapfill --gaf` to report or resolve graph traversal support.
- Hi-C pairs: optional graph-node contact counts used by
  `chromo gapfill --hic-pairs` as an additional conservative branch-support
  signal.
- external patch candidates: optional table/FASTA pairs used by
  `chromo eval gapfill --patch-table` and `chromo gapfill --patch-table` to
  compare outside patch sequences with graph-derived fills as review evidence.

The [Architecture]({{ '/architecture/' | relative_url }}) page maps these
evidence files to the subcommands, modes, and parameters that activate them,
including which uses are report-only and which can affect sequence output.

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

For graph review, use the graph closest to the FASTA being sorted or filled.
The simplest case is a contig-level GFA whose `S` segment names match the FASTA:

```text
assembly.fa              # FASTA passed to chromo sort/fix/gapfill
assembly_graph.gfa       # GFA whose S records match assembly.fa sequence IDs
```

ChromoSort reads GFA `S` segment records, `L` link records, and `P`/`W`
path/walk records. Segment sequences are required only when
`chromo gapfill --apply` may insert graph sequence. Report-only graph evidence
can still use segments with `*` sequence fields when lengths are provided with
`LN:i`.

When alignments and dot plots use `p_ctg.fa` but graph evidence comes from
`p_utg.gfa` or `p_utg.noseq.gfa`, first check whether the GFA has path/walk
records whose names match the contig FASTA:

```bash
chromo graph-map \
  --ctg-fasta assembly.p_ctg.fa \
  --utg-gfa assembly.p_utg.noseq.gfa \
  --output-prefix results/sample.graphmap
```

If the GFA contains only `S`, `L`, and assembler read-alignment records, it can
still be useful as topology evidence, but it cannot by itself project unitig
coordinates onto contigs.

For graph-aware gap filling, the same rule applies more strictly:

| Graph input | Planning use | Applied sequence use |
| --- | --- | --- |
| Contig-level GFA with matching segment names and sequences | Direct graph-path planning. | Yes. |
| Contig-level `.noseq.gfa` with matching segment names and `LN:i` lengths | Topology and report-only context. | No. |
| Unitig-level GFA with matching `P`/`W` contig paths and sequences | Projection-backed planning with `chromo gapfill --project-gfa-paths`, optionally confirmed first with `chromo graph-map`. | Yes, when terminal unitig sequence validates against component ends. |
| Unitig-level `.noseq.gfa` with matching `P`/`W` paths | Projection and topology review. | No. |
| Unitig-level GFA without matching `P`/`W` paths | Limited node-level context only. | No. |

If you already have a scaffold FASTA, keep or create an AGP/component map that
connects scaffold gaps back to the original contigs. The scaffold record itself
usually will not match the assembler graph, but `chromo gapfill` with
`--scaffold-fasta --agp` uses AGP to recover the component identities needed
for direct contig-level graph filling or, with `--project-gfa-paths`, the path
names needed to project component ends onto terminal unitigs. Scaffold FASTA
without AGP is rejected for graph-aware filling because the N gaps can no longer
be tied safely to graph segments.

### Which PAF Files to Keep

The main PAF file for ChromoSort is a reference-to-assembly whole-genome
alignment:

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
```

If you run `chromo fix` and create a new fixed FASTA, re-align the fixed FASTA
and keep that second PAF beside the fixed results:

```bash
minimap2 -x asm5 -c -t 16 --secondary=no reference.fa results/sample.fixed.fa \
  > paf/sample.fixed.ref_vs_asm.paf
```

The original PAF explains the raw assembly; the fixed PAF explains the edited
assembly. Do not mix a fixed FASTA with an old raw-assembly PAF unless you are
only doing a very specific manual comparison.

For `chromo gapfill --ref-paf`, the PAF is used to place intermediate graph
nodes into the expected reference-space gap. If your GFA segment names are the
same as your assembly FASTA record names, the ordinary reference-vs-assembly
PAF above can be used. If your graph has unitig names that are not the same as
the FASTA records passed to `chromo sort`, export a graph-node FASTA from the
GFA or assembler output and align those graph nodes to the reference:

```bash
minimap2 \
  -x asm5 \
  -c \
  -t 16 \
  --secondary=no \
  reference.fa \
  graph_nodes.fa \
  > paf/graph_nodes_to_ref.paf

chromo gapfill \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --gfa assembly_graph.gfa \
  --ref-paf paf/graph_nodes_to_ref.paf \
  --output-prefix results/sample.gapfill
```

Use the PAF whose query names match the GFA segment names being evaluated by
gapfill.

For `chromo gapfill --read-paf`, use long reads mapped back to the pre-scaffold
contig FASTA or ordered FASTA records. The target names should match the
assignment `new_name`, original `contig`, or ordered FASTA names. In
`--scaffold-fasta --agp` mode, target names should match AGP component IDs.
This PAF is used to count reads that anchor on both sides of an adjacent
component junction. It is bridge evidence only; inserted bases still come from
a sequence-bearing GFA path. PAF targets that are only scaffold record names are
not currently interpreted as reads spanning scaffold N blocks.

### Creating GAF Read-to-Graph Alignments

GAF is a graph-alignment format. ChromoSort uses it as optional read-path
evidence for fix, scaffold, manual, and gapfill review. A typical source is a
long-read-to-GFA alignment from a graph aligner:

```bash
GraphAligner \
  -g assembly_graph.gfa \
  -f reads.fastq.gz \
  -a graph_alignments/sample.reads_to_graph.gaf
```

`GraphAligner` is an optional external tool; it is not needed for the core
sorting/fixing/scaffolding workflow. ChromoSort reads the query name, path
string, and MAPQ columns, filters with `--min-gaf-mapq`, and counts how many
read paths contain a graph node or candidate graph traversal. In `chromo eval
fix`, GAF is advisory node/traversal evidence. In `chromo eval scaffold`, it
reports support for selected and alternate graph paths between adjacent
contigs. In `chromo eval gapfill` and `chromo gapfill --gaf`, one candidate
path with unique support above `--min-gaf-path-support` can resolve an
otherwise ambiguous GFA branch. Tied, weak, absent, or conflicting support
keeps the event reviewable instead of forcing a hidden choice.

For large full-depth HiFi datasets, running GraphAligner against the whole GFA
can be expensive. Use `chromo eval all` plus `chromo gafprep` when you want GAF
evidence near fix, scaffold, and gapfill review targets from one selected read
subset:

```bash
minimap2 -x map-hifi -c -t 16 --secondary=no results/sample.ordered.fa reads.fastq.gz \
  > reads_to_ordered.paf

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
  --output-prefix results/sample.gafprep

bash results/sample.gafprep.graphaligner.sh
```

This workflow writes selected reads and a sanitized GFA first; GraphAligner
still writes `results/sample.gafprep.gaf`. Targeted GAF is review evidence, not
automatic sequence validation. The read PAF, `gafprep --assembly-fasta`, and
review tables should name the same assembly FASTA records.

Choose the graph for the question being asked. A contig-level `ctg.gfa` whose
`S` names match the assembly FASTA is best for direct contig-node and junction
review. A unitig-level `utg.gfa` can be more informative for branch structure,
but unitig coordinates are not contig coordinates; use matching `P`/`W` paths,
`chromo graph-map`, or `chromo plot --gfa-overlay` when you need to relate
unitig context back to contig breakpoints. A `.noseq.gfa` is useful for
topology reports, but a sequence-bearing graph is normally needed for
GraphAligner read-to-graph alignment.

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

## Preflight Checks

Before starting a long workflow, inspect a few lines from each file and confirm
the names describe the same records:

```bash
# FASTA IDs.
grep '^>' reference.fa | head
grep '^>' assembly.fa | head

# MUMmer coords should end each data row with reference and query names.
head -40 mummer/sample.coords

# PAF query names are column 1; reference/target names are column 6.
awk 'BEGIN{FS="\t"} !/^#/ && NF>=12 {print $1, $6; exit}' paf/sample.paf

# GFA graph nodes are S records. These names should match graph-aware inputs.
awk 'BEGIN{FS="\t"} $1=="S" {print $2; count++} count==10 {exit}' assembly_graph.gfa

# GAF path strings are column 6; MAPQ is column 12.
awk 'BEGIN{FS="\t"} !/^#/ && NF>=12 {print $1, $6, $12; exit}' reads_to_graph.gaf

# Hi-C pair files should have node names and integer counts.
head graph_contacts.tsv
```

Common symptoms and likely input causes:

| Symptom | Likely cause | What to check |
| --- | --- | --- |
| Many `no_alignment` rows in `chromo sort` | Alignment query names do not match assembly FASTA IDs, or the wrong FASTA was aligned | Compare `grep '^>' assembly.fa` with PAF column 1 or coords query names. |
| Reference names are missing or all contigs map unexpectedly | Alignment target/reference names do not match `--ref-fasta`, or reference/query order was swapped when aligning | Compare reference FASTA IDs with PAF column 6 or coords reference names. |
| `chromo scaffold` says the ordered FASTA is missing kept assignments | The ordered FASTA and assignment TSV come from different `chromo sort` runs | Use `<prefix>.ordered.fa` with the matching `<prefix>.contig_assignments.tsv`. |
| Graph reports show many missing nodes | GFA segment names do not match assembly names, ordered FASTA names, or assignment names | Compare GFA `S` names against original contig IDs and ChromoSort `new_name` values. |
| `chromo gapfill` reports flank sequence mismatch | The GFA sequence does not match the ordered FASTA flank records | Use a GFA from the same assembly stage or regenerate graph evidence after editing. |
| `--ref-paf` support is always zero | PAF query names do not match GFA intermediate node names, or placements fall outside the expected reference gap | Align the graph-node/unitig FASTA to the reference and inspect PAF column 1. |
| Reviewed gapfill plan is rejected as stale | The current graph path no longer matches the exported plan | Regenerate the gapfill plan after changing FASTA, assignments, GFA, or path-search settings. |

When in doubt, keep input files grouped by assembly stage:

```text
raw/
  assembly.fa
  sample.coords
  sample.ref_vs_asm.paf

fixed/
  sample.fixed.fa
  sample.fixed.coords
  sample.fixed.ref_vs_asm.paf
  sample.fixed.contig_assignments.tsv
```

This makes it harder to accidentally use raw-assembly evidence with a fixed or
manually edited FASTA.
