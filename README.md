# ChromoSort

Reference-order genome assembly contigs from MUMmer alignments.

ChromoSort takes an assembly FASTA plus MUMmer `show-coords` output against a
reference genome, assigns each contig to its best reference chromosome, removes
low-value duplicate overlaps, and writes a new FASTA ordered by the reference.
It is meant for reusable, cross-species genome assembly cleanup and reporting,
not for one specific organism or naming convention.

## What ChromoSort Does

Given:

- a reference genome FASTA,
- an assembly FASTA,
- a MUMmer `show-coords` file from the reference-vs-assembly alignment,

ChromoSort:

1. Parses `show-coords` rows for reference chromosome, assembly contig, position,
   length, percent identity, and orientation.
2. Merges overlapping alignment intervals so coverage is not inflated by
   repeated MUMmer rows.
3. Assigns each contig to the reference chromosome with the strongest merged
   query coverage.
4. Applies basic quality thresholds for aligned bp, query coverage, and
   best-reference share.
5. Removes duplicate-overlap contigs that add little or no new reference
   coverage after better contigs on the same chromosome have claimed intervals.
6. Sorts retained contigs by reference FASTA order and reference start.
7. Writes an ordered FASTA with names like `chromosome_contig`.
8. Writes TSV reports that explain every keep/reject decision.

ChromoSort does not build a pseudomolecule, insert gaps, call variants, polish
sequence, or infer unaligned sequence placement. It keeps full contigs, in a
reference-informed order, with transparent reporting.

## Installation With Mamba

Install Mambaforge, Miniforge, or another conda-compatible distribution with
`mamba`, then create the environment:

```bash
git clone https://github.com/rotheconrad/chromosort.git
cd chromosort

mamba env create -f environment.yml
mamba activate chromosort
```

The environment installs:

- Python
- MUMmer 4 (`nucmer`, `delta-filter`, `show-coords`)
- pytest for the test suite
- ChromoSort in editable mode

Check the command:

```bash
chromosort --help
chromosort-fix-contigs --help
```

Run tests:

```bash
pytest
```

The tests are also compatible with the Python standard-library runner:

```bash
python -m unittest discover -s tests -v
```

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

# Optional visual inspection.
mummerplot \
  -p "mummer/plot_${name}" \
  -R "$ref" \
  --postscript \
  --large \
  --layout \
  --fat \
  "mummer/${name}.filter"
```

### Why These MUMmer Choices?

`nucmer` aligns the reference and assembly at whole-genome scale. The `-c`
minimum cluster length removes very small seeds that are often unhelpful for
chromosome-scale contig ordering.

`delta-filter -1` is used because this workflow wants a primary placement for
each contig rather than every local repeat hit. It reduces redundant alignments
before ChromoSort applies its own interval merging and overlap logic.

`delta-filter -i` and `-l` enforce minimum identity and alignment length before
reporting. Use stricter values for very similar assemblies, and looser values
for distant species or more fragmented assemblies.

`show-coords -r -c -l` reports alignments sorted by reference and includes
coverage and sequence lengths. ChromoSort reads the coordinates and lengths,
then recomputes merged coverage itself.

## Running ChromoSort

```bash
chromosort \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.fa \
  --coords mummer/sample.coords \
  --output-prefix results/sample \
  --orient-to-reference
```

This writes:

| Output | Description |
| --- | --- |
| `results/sample.ordered.fa` | Retained contigs, ordered by reference chromosome and position. |
| `results/sample.contig_assignments.tsv` | One row per assembly contig with final status and assignment metrics. |
| `results/sample.contig_ref_matches.tsv` | One row per contig-reference match before final assignment. |
| `results/sample.chromosome_summary.tsv` | One row per reference sequence with ordered contig lists and covered reference bp. |
| `results/sample.run_summary.txt` | Inputs, thresholds, output paths, and status counts. |

Optional discarded FASTA:

```bash
chromosort \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.fa \
  --coords mummer/sample.coords \
  --output-prefix results/sample \
  --discarded-fasta results/sample.discarded.fa
```

## Key Parameters

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `--min-aligned-bp` | `100000` | Minimum merged query-aligned bp required before a contig can be kept. |
| `--min-query-cov` | `0.50` | Minimum fraction of the contig covered by its best reference match. |
| `--min-best-ref-share` | `0.50` | Minimum fraction of all matched bp that must belong to the best reference chromosome. |
| `--min-segment-idy` | `0.0` | Ignore individual `show-coords` rows below this percent identity. |
| `--min-novel-ref-bp` | `50000` | Keep an otherwise-good contig if it adds at least this many new reference bp. |
| `--min-novel-ref-frac` | `0.20` | Keep an otherwise-good contig if this fraction of its reference span is novel. |
| `--orient-to-reference` | off | Reverse-complement retained contigs whose dominant alignment is reverse strand. |
| `--no-overlap-filter` | off | Keep all contigs passing basic match thresholds, even if they overlap better contigs. |

For small microbial genomes or tiny test fixtures, lower `--min-aligned-bp`.
For large plant or animal assemblies, the defaults are intentionally
conservative.

## Assignment Status Values

`kept`: written to the ordered FASTA.

`no_alignment`: no usable `show-coords` rows were found for this contig.

`below_min_aligned_bp`: best match did not meet `--min-aligned-bp`.

`below_min_query_cov`: best match did not meet `--min-query-cov`.

`ambiguous_ref_match`: the best chromosome did not dominate the contig's total
matched bp enough to pass `--min-best-ref-share`.

`duplicate_overlap`: the contig passed the basic match thresholds, but better
contigs on the same reference chromosome already covered nearly all of its
reference span. These contigs are not written to the ordered FASTA.

## Reasoning Behind the Decisions

### Use `show-coords`, Not `show-tiling`

`show-tiling` can be useful, but it is another derived representation and is not
always produced. `show-coords` from a filtered delta contains the required
information: reference name, query name, reference coordinates, query
coordinates, alignment length, sequence length, percent identity, and strand.
Using `show-coords` keeps the workflow closer to standard MUMmer output and
avoids an extra preprocessing step.

### Merge Intervals Before Calculating Coverage

MUMmer can report overlapping rows for the same contig-reference pair. Summing
raw `LEN 2` values can produce apparent coverage greater than 100 percent. By
merging query intervals and reference intervals first, ChromoSort estimates
coverage from the unique aligned span instead of from repeated rows.

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

### Preserve Full Contigs

ChromoSort does not trim contigs to alignment spans. It writes the full input
contig sequence because unaligned ends may be real assembly sequence. The output
is an ordered contig FASTA, not a hard reference-guided reconstruction.

### Make Orientation Optional

`--orient-to-reference` reverse-complements retained contigs whose dominant
alignment is on the reverse strand. This is helpful for downstream inspection
and plotting. It is optional because some workflows prefer to preserve original
assembly orientation exactly.

## Batch Example

```bash
mkdir -p results

for asm in assemblies/*.fa; do
  sample=$(basename "$asm" .fa)
  chromosort \
    --ref-fasta reference.fa \
    --assembly-fasta "$asm" \
    --coords "mummer/${sample}.coords" \
    --output-prefix "results/${sample}" \
    --orient-to-reference
done
```

## Fixing User-Nominated Chimeric Contigs

ChromoSort also includes a focused contig-splitting command:

```bash
chromosort-fix-contigs \
  --assembly-fasta assembly.fa \
  --coords mummer/sample.coords \
  --contigs contig_04 contig_12 \
  --output-fasta results/sample.fixed.fa \
  --report results/sample.fixed_contigs.tsv
```

This command is intentionally user-directed. It does not scan the whole
assembly and decide to cut contigs on its own. Instead, you provide the contig
names you already suspect are chimeric, usually after inspecting dot plots,
assignment reports, or other QC evidence.

For each requested contig, `chromosort-fix-contigs`:

1. Reads passing `show-coords` alignment segments for that contig.
2. Sorts those segments by query-coordinate order along the assembly contig.
3. Merges nearby neighboring rows that map to the same reference sequence and
   orientation.
4. Looks for reference transitions along the contig.
5. Places breakpoints halfway between neighboring reference blocks.
6. Replaces the original contig with two or more pieces in the output FASTA.
7. Writes a TSV report with slice coordinates, reference labels, orientation,
   identity, and split status.

By default, unrequested contigs are copied unchanged, producing a full fixed
assembly FASTA. Use `--pieces-only` to write only the split pieces.

### Chimeric Contig Naming

Split pieces are named:

```text
REFERENCE-CONTIG-PART
```

For example, a contig named `contig_04` with its first half matching `chrom02`
and second half matching `chrom07` becomes:

```text
chrom02-contig_04-a
chrom07-contig_04-b
```

If a contig has more than one breakpoint, ChromoSort emits as many pieces as
the ordered reference blocks require. For example, if `contig_12` has a small
`chrom04` block, a larger `chrom05` block, and another `chrom04` block, it
becomes:

```text
chrom04-contig_12-a
chrom05-contig_12-b
chrom04-contig_12-c
```

The reference names and contig names are not hard-coded. Whatever identifiers
appear in your FASTA and MUMmer output are used. Change the separator with
`--name-separator`.

### Chimeric Contig Parameters

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `--contigs` | none | Space-separated names of contigs to inspect and split. |
| `--contigs-file` | none | Optional file with one contig name per line. |
| `--min-segment-bp` | `10000` | Minimum alignment segment length used to infer split blocks. |
| `--min-segment-idy` | `0.0` | Minimum percent identity for split-informing alignment rows. |
| `--max-merge-gap` | `1000` | Merge nearby same-reference rows separated by this many query bp or less. |
| `--min-piece-bp` | `1` | Do not emit split pieces shorter than this length. |
| `--orient-to-reference` | off | Reverse-complement split pieces from reverse-strand blocks. |
| `--pieces-only` | off | Write only split pieces instead of a full fixed assembly FASTA. |

### Why User-Nominated Splitting?

Cutting contigs is a stronger intervention than ordering contigs. A reference
transition can reflect a real assembly chimera, but it can also reflect
structural variation, assembly graph complexity, misassembly in the reference,
or poor alignment around repeats. Requiring an explicit contig list keeps this
step auditable: ChromoSort proposes breakpoints from MUMmer coordinates, but the
user decides which contigs are appropriate to cut.

The synthetic test data under `tests/data/chimeric` includes two chimeric
contigs: one split roughly half-and-half between two reference chromosomes, and
one with 25 percent of its sequence matching one chromosome and 75 percent
matching another.

## Development

```bash
mamba env create -f environment.yml
mamba activate chromosort
pytest
```

The tests use tiny synthetic FASTA and coords files under `tests/data`. They do
not require running MUMmer; MUMmer is included in the environment so users can
generate real input files.

## Citation

If you use ChromoSort, cite this repository and cite MUMmer for the alignment
files used by the workflow.
