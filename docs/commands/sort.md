---
title: chromo sort
description: Usage, outputs, parameters, and reasoning for chromo sort.
---

# chromo sort

Use `chromo sort` when the goal is to find the best matched assembly contigs for
each reference chromosome and write them in reference order.

The coords or PAF file must describe the same assembly FASTA passed to
`--assembly-fasta`. If you later want to use `<prefix>.ordered.fa` as the input
FASTA for `chromo fix`, `chromo plot`, or another alignment-dependent command,
re-run MUMmer or minimap2 against `<prefix>.ordered.fa` first. The original
alignment remains valid for reviewing the original assembly and for
`chromo plot --assignments`, but it is not a new alignment of the ordered FASTA.
For help reading those review plots, see the
[dot-plot guide]({{ '/dot-plots/' | relative_url }}).

## What `chromo sort` Does

Given a reference FASTA, assembly FASTA, and MUMmer coords or minimap2 PAF file,
`chromo sort`:

1. Parses alignment rows for reference chromosome, assembly contig,
   coordinates, alignment length, percent identity, and orientation.
2. Merges overlapping alignment intervals so coverage is not inflated by
   repeated rows.
3. Assigns each contig to the reference chromosome with the strongest merged
   query coverage.
4. Applies thresholds for aligned bp, query coverage, and best-reference share.
5. Flags strong multi-reference contigs as possible `chromo fix` candidates.
6. Classifies overlap shape against already-kept reference intervals.
7. Removes contained/internal duplicate overlaps that add little or no new
   reference coverage.
8. Keeps or rescues terminal overlaps when they contribute enough one-sided
   extension.
9. Protects flagged split candidates from silent duplicate-overlap removal.
10. Sorts retained contigs by reference FASTA order and reference start.
11. Writes an ordered FASTA with names like `chromosome_contig`.
12. Optionally writes GFA graph evidence for assignments and overlap calls.
13. Writes TSV reports that explain every keep/reject decision.

## Run `chromo sort`

```bash
chromo sort \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.fa \
  --coords mummer/sample.coords \
  --output-prefix results/sample \
  --orient-to-reference
```

Optional discarded FASTA:

```bash
chromo sort \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.fa \
  --coords mummer/sample.coords \
  --output-prefix results/sample \
  --discarded-fasta results/sample.discarded.fa
```

Optional graph evidence for assignment review:

```bash
chromo sort \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.fa \
  --paf minimap2/sample.paf \
  --output-prefix results/sample \
  --gfa assembly_graph.gfa
```

## `chromo sort` Outputs

| Output | Description |
| --- | --- |
| `<prefix>.ordered.fa` | Retained contigs, ordered by reference chromosome and position. |
| `<prefix>.contig_assignments.tsv` | One row per assembly contig with final status and assignment metrics. |
| `<prefix>.contig_ref_matches.tsv` | One row per contig-reference match before final assignment. |
| `<prefix>.chromosome_summary.tsv` | One row per reference sequence with ordered contig lists and covered reference bp. |
| `<prefix>.graph_assignments.tsv` | Optional report-only graph evidence for assignment and duplicate-overlap decisions when `--gfa` is provided. |
| `<prefix>.run_summary.txt` | Inputs, thresholds, output paths, status counts, and PAF diagnostics when `--paf` is used. |

### Example `chromo sort` Output

**Table 1. Example `contig_assignments.tsv` rows.** Selected columns from the
synthetic fixture show the main status classes: kept contigs, a duplicate
overlap, and an unaligned contig.

| contig | status | kept | new_name | assigned_ref | query_cov | avg_identity | overlap_class |
| --- | --- | --- | --- | --- | ---: | ---: | --- |
| `contigA` | `kept` | `yes` | `chr1_contigA` | `chr1` | `1.0000` | `99.000` | `.` |
| `contigB` | `kept` | `yes` | `chr1_contigB` | `chr1` | `1.0000` | `98.500` | `.` |
| `contigDup` | `duplicate_overlap` | `no` | `.` | `chr1` | `1.0000` | `99.000` | `contained_overlap` |
| `contigNo` | `no_alignment` | `no` | `.` | `.` | `0.0000` | `0.000` | `.` |

**Table 2. Example `chromosome_summary.tsv` rows.** The summary table groups
retained contigs by reference and records the covered reference span.

| ref | ref_length | kept_contigs | covered_ref_bp | ref_cov | ordered_new_names |
| --- | ---: | ---: | ---: | ---: | --- |
| `chr1` | `120` | `2` | `120` | `1.0000` | `chr1_contigA,chr1_contigB` |
| `chr2` | `100` | `1` | `60` | `0.6000` | `chr2_contigC` |

**Listing 1. Example ordered FASTA headers.** Headers retain the original contig
name, assigned reference interval, orientation, query coverage, and mean
identity used for audit.

```text
>chr1_contigA original=contigA ref=chr1 ref_start=1 ref_end=80 orientation=+ reverse_complemented=no query_cov=1.0000 avg_identity=99.000
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>chr2_contigC original=contigC ref=chr2 ref_start=10 ref_end=69 orientation=- reverse_complemented=yes query_cov=1.0000 avg_identity=97.000
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
```

## `chromo sort` Parameters

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `--coords` | required unless `--paf` | MUMmer `show-coords` alignment file. |
| `--paf` | required unless `--coords` | minimap2 PAF alignment file. |
| `--gfa` | none | Optional assembly graph GFA for report-only evidence about resolved graph nodes, node degree, self loops, and direct graph links to overlap-best contigs. |
| `--graph-guard` | off | Requires `--gfa`; emits conservative warnings when graph links contradict or complicate duplicate/terminal-overlap decisions without changing sorting output. |
| `--min-aligned-bp` | `100000` | Minimum merged query-aligned bp required before a contig can be kept. |
| `--min-query-cov` | `0.50` | Minimum fraction of the contig covered by its best reference match. |
| `--min-best-ref-share` | `0.50` | Minimum fraction of all matched bp that must belong to the best reference chromosome. |
| `--large-alignment-min-bp` | `10000000` | Rescue a near-threshold contig when its best reference match has at least this many merged query-aligned bp. Set 0 to disable. |
| `--large-alignment-min-query-cov` | `0.45` | Minimum query coverage required by large-alignment rescue. |
| `--min-segment-idy` | `0.0` | Ignore individual alignment rows below this percent identity. |
| `--min-mapq` | `0` | Ignore PAF rows below this MAPQ. Ignored for coords. |
| `--include-secondary-paf` | off | Include PAF rows marked `tp:A:S`; skipped by default. |
| `--min-novel-ref-bp` | `50000` | Keep an otherwise-good contig if it adds at least this many new reference bp. |
| `--min-novel-ref-frac` | `0.20` | Keep an otherwise-good contig if this fraction of its reference span is novel. |
| `--overlap-mode` | `span` | Use broad first-to-last reference spans for duplicate filtering. Set `alignment` to use exact merged alignment intervals. |
| `--novel-ref-criteria` | `both` | Require both novel-bp and novel-fraction thresholds during duplicate filtering. Set `either` for the older permissive behavior. |
| `--min-terminal-extension-bp` | `100000` | Rescue a terminally overlapping contig if its one-sided novel extension has at least this many reference bp. |
| `--min-terminal-extension-frac` | `0.02` | Rescue a terminally overlapping contig if its one-sided novel extension covers at least this fraction of the overlap-filter interval. |
| `--no-terminal-overlap-rescue` | off | Report terminal overlaps without rescuing contigs that fail the standard novel-reference thresholds. |
| `--split-candidate-min-aligned-bp` | `100000` | Minimum merged query-aligned bp on at least two references before split-candidate protection applies. |
| `--split-candidate-min-query-frac` | `0.05` | Minimum query-length fraction on at least two references before split-candidate protection applies. |
| `--split-candidate-max-best-share` | `0.95` | Do not flag contigs whose best reference accounts for more than this share of total matched bp. |
| `--orient-to-reference` | off | Reverse-complement retained contigs whose dominant alignment is reverse strand. |
| `--no-overlap-filter` | off | Keep all contigs passing basic match thresholds, even if they overlap better contigs. |
| `--no-protect-split-candidates` | off | Let strong multi-reference split candidates be filtered like ordinary contigs. |

For small microbial genomes or tiny test fixtures, lower `--min-aligned-bp`
and `--min-novel-ref-bp`.
For large plant or animal assemblies, the defaults are intentionally
conservative.

For assemblies with many short alternate/contaminant fragments around strong
chromosome-scale contigs, a stricter cleanup pass is often appropriate:

```bash
chromo sort \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.fa \
  --coords mummer/sample.coords \
  --output-prefix results/sample \
  --min-aligned-bp 1000000 \
  --min-novel-ref-frac 0.5
```

## `chromo sort` Status Values

`kept`: written to the ordered FASTA.

`kept_split_candidate`: written to the ordered FASTA and flagged as a strong
multi-reference candidate for `chromo fix` review. These contigs are protected
from duplicate-overlap removal by default.

`kept_large_alignment`: written to the ordered FASTA because it has a very large
best-reference alignment even though its query coverage is slightly below
`--min-query-cov`.

`kept_terminal_overlap`: written to the ordered FASTA because it overlaps an
already-kept contig at one end but still contributes enough novel terminal
reference span. This status is also used when the terminal-extension rescue keeps
a contig that would otherwise fail the standard novel-reference fraction.

`no_alignment`: no usable alignment rows were found for this contig.

`below_min_aligned_bp`: best match did not meet `--min-aligned-bp`.

`below_min_query_cov`: best match did not meet `--min-query-cov`.

`ambiguous_ref_match`: the best chromosome did not dominate the contig's total
matched bp enough to pass `--min-best-ref-share`.

`terminal_overlap`: the contig passed the basic match thresholds, but a stronger
contig already covered most of its reference span and the one-sided terminal
extension did not pass the keep or rescue thresholds.

`duplicate_overlap`: the contig passed the basic match thresholds, but better
contigs on the same reference chromosome already covered nearly all of its
reference span in a contained or internal-overlap pattern. These contigs are not
written to the ordered FASTA.

## Reasoning Behind `chromo sort`

### Use Segment Coordinates, Not Tiling Summaries

`show-tiling` can be useful for MUMmer workflows, but it is another derived
representation and is not always produced. `show-coords` from a filtered delta
and minimap2 PAF both contain the required information: reference name, query
name, reference coordinates, query coordinates, alignment length, sequence
length, percent identity, and strand. ChromoSort normalizes either format into
the same internal segment representation before sorting, fixing, or plotting.

### Merge Intervals Before Calculating Coverage

Whole-genome aligners can report overlapping rows for the same contig-reference
pair. Summing raw row lengths can produce apparent coverage greater than 100
percent. By merging query intervals and reference intervals first, ChromoSort
estimates coverage from the unique aligned span instead of from repeated rows.

### Assign Contigs by Best Reference Share

Many genomes contain repeats, paralogous regions, translocations, or retained
haplotigs. A contig may have alignments to more than one chromosome. ChromoSort
chooses the chromosome with the largest merged query-aligned bp, then requires
that match to represent a configurable share of all matched bp. This keeps clear
placements and flags ambiguous ones.

### Filter Duplicate Overlaps After Assignment

The first assignment pass asks, "Does this contig have a good placement?" The
overlap pass asks, "Does this contig add useful new reference coverage beyond
better contigs already kept?" This second question is important for assemblies
that include short duplicate fragments, alternate haplotigs, or local redundant
contigs. Rejected contigs are marked `duplicate_overlap`, with novel coverage
and the strongest overlapping kept contig reported.

By default, duplicate filtering uses each contig's broad first-to-last reference
span and requires both the novel-bp and novel-fraction thresholds before
retaining a lower-ranked contig. This is intentionally stricter than exact
alignment-block overlap: whole-genome alignments can be fragmented by repeats,
local variation, or filtering, and short duplicate fragments often land in the
internal gaps of a stronger chromosome-scale contig. Use
`--overlap-mode alignment --novel-ref-criteria either` for the older, more
permissive behavior.

Terminal overlaps are separated from contained/internal duplicates. If the novel
reference interval sits at one end of the lower-ranked contig, the assignment
report records `overlap_class=terminal_overlap`, the extension side, and the
extension bp/fraction. Terminal overlaps that pass the standard novel-reference
thresholds are kept as `kept_terminal_overlap`. Mostly overlapping terminal
contigs can still be rescued when their one-sided extension passes
`--min-terminal-extension-bp` and `--min-terminal-extension-frac`.

Strong multi-reference contigs are treated differently. If at least two
references carry substantial query support and no single reference explains
nearly all matched bp, the contig is marked `kept_split_candidate` and retained
even if its best-reference interval overlaps a better contig or its best single
reference is below `--min-query-cov`. Its secondary supported reference spans
also block lower-ranked duplicate fragments during overlap filtering. This keeps
likely `chromo fix` targets available for review instead of removing them during
sorting while still reducing clutter around those loci.

Very large single-reference alignments are also rescued by default when they
land just below `--min-query-cov`. This prevents chromosome-scale contigs with
fragmented alignments from being discarded in favor of smaller redundant pieces.

### Treat Graph Evidence as Review Context

When `--gfa` is provided, `chromo sort` writes
`<prefix>.graph_assignments.tsv` beside the normal assignment report. The graph
report resolves each assembly contig to a GFA segment when possible, records
node length, sequence availability, in/out degree, self loops, and direct graph
links to the contig that drove a duplicate-overlap decision. It does not change
which contigs are kept, renamed, or written to FASTA.

With `--graph-guard`, `chromo sort` also writes stderr warnings when a contig
classified as `duplicate_overlap` or `terminal_overlap` has a direct GFA link to
the overlap-best contig. This is a review signal only; it does not rescue or
discard contigs.

### Preserve Full Contigs

`chromo sort` does not trim contigs to alignment spans. It writes the full input
contig sequence because unaligned ends and terminal extensions may be real
assembly sequence. The output is an ordered contig FASTA, not a hard
reference-guided reconstruction. Optional overlap trimming is handled only by
`chromo scaffold`, where the sequence-changing policy is explicit.

### Make Orientation Optional

`--orient-to-reference` reverse-complements retained contigs whose dominant
alignment is on the reverse strand. This is helpful for downstream inspection
and plotting. It is optional because some workflows prefer to preserve original
assembly orientation exactly.

## Batch Sorting Example

```bash
mkdir -p results

for asm in assemblies/*.fa; do
  sample=$(basename "$asm" .fa)
  chromo sort \
    --ref-fasta reference.fa \
    --assembly-fasta "$asm" \
    --coords "mummer/${sample}.coords" \
    --output-prefix "results/${sample}" \
    --orient-to-reference
done
```
