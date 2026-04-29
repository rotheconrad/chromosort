# ChromoSort

Reference-guided genome assembly utilities for sorting contigs and splitting
user-flagged or auto-detected chimeric contigs from MUMmer alignments.

ChromoSort provides one command, `chromo`, with two subcommands:

- `chromo sort` assigns assembly contigs to reference chromosomes, removes
  low-value duplicate overlaps, and writes a reference-ordered FASTA.
- `chromo fix` splits contigs with chromosome transitions or inversion blocks
  into reference-labeled pieces, either from a user-supplied contig list or with
  `--auto` detection.

Both workflows use standard MUMmer `show-coords` output and are designed for
reuse across species and genome assembly projects. ChromoSort does not build
pseudomolecules, insert gaps, polish sequence, call variants, or force contigs
to match a reference. It keeps full sequence pieces and writes TSV reports so
each keep, reject, or split decision is auditable.

## Table of Contents

- [Quick Start](#quick-start)
- [Installation With Mamba](#installation-with-mamba)
- [Creating Input Files With MUMmer](#creating-input-files-with-mummer)
  - [Why These MUMmer Choices?](#why-these-mummer-choices)
- [chromo sort](#chromo-sort)
  - [What `chromo sort` Does](#what-chromo-sort-does)
  - [Run `chromo sort`](#run-chromo-sort)
  - [`chromo sort` Outputs](#chromo-sort-outputs)
  - [`chromo sort` Parameters](#chromo-sort-parameters)
  - [`chromo sort` Status Values](#chromo-sort-status-values)
  - [Reasoning Behind `chromo sort`](#reasoning-behind-chromo-sort)
  - [Batch Sorting Example](#batch-sorting-example)
- [chromo fix](#chromo-fix)
  - [What `chromo fix` Does](#what-chromo-fix-does)
  - [Run `chromo fix` With User-Nominated Contigs](#run-chromo-fix-with-user-nominated-contigs)
  - [Run `chromo fix` With Auto Detection](#run-chromo-fix-with-auto-detection)
  - [`chromo fix` Outputs](#chromo-fix-outputs)
  - [`chromo fix` Naming](#chromo-fix-naming)
  - [`chromo fix` Parameters](#chromo-fix-parameters)
  - [Reasoning Behind `chromo fix`](#reasoning-behind-chromo-fix)
- [Development](#development)
- [Citation](#citation)
- [Contact](#contact)
- [Funding](#funding)
- [Acknowledgements](#acknowledgements)
- [Version History](#version-history)

## Quick Start

```bash
git clone https://github.com/rotheconrad/chromosort.git
cd chromosort

mamba env create -f environment.yml
mamba activate chromosort

chromo --help
chromo sort --help
chromo fix --help
```

Typical reference-ordering run:

```bash
chromo sort \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.fa \
  --coords mummer/sample.coords \
  --output-prefix results/sample \
  --orient-to-reference
```

Typical chimeric-contig fixing run:

```bash
chromo fix \
  --assembly-fasta assembly.fa \
  --coords mummer/sample.coords \
  --contigs contig_04 contig_12 \
  --output-fasta results/sample.fixed.fa \
  --report results/sample.fixed_contigs.tsv
```

Auto-detect contigs with chromosome or orientation transitions:

```bash
chromo fix \
  --assembly-fasta assembly.fa \
  --coords mummer/sample.coords \
  --auto \
  --output-fasta results/sample.fixed.fa \
  --report results/sample.fixed_contigs.tsv
```

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

Legacy command aliases are retained for compatibility:

- `chromosort` is equivalent to `chromo sort`
- `chromosort-fix-contigs` is equivalent to `chromo fix`

New workflows should use `chromo sort` and `chromo fix`.

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

## chromo sort

Use `chromo sort` when the goal is to find the best matched assembly contigs for
each reference chromosome and write them in reference order.

### What `chromo sort` Does

Given a reference FASTA, assembly FASTA, and MUMmer `show-coords` file,
`chromo sort`:

1. Parses `show-coords` rows for reference chromosome, assembly contig,
   coordinates, alignment length, percent identity, and orientation.
2. Merges overlapping alignment intervals so coverage is not inflated by
   repeated MUMmer rows.
3. Assigns each contig to the reference chromosome with the strongest merged
   query coverage.
4. Applies thresholds for aligned bp, query coverage, and best-reference share.
5. Removes contigs that duplicate already-kept reference intervals and add
   little or no new reference coverage.
6. Sorts retained contigs by reference FASTA order and reference start.
7. Writes an ordered FASTA with names like `chromosome_contig`.
8. Writes TSV reports that explain every keep/reject decision.

### Run `chromo sort`

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

### `chromo sort` Outputs

| Output | Description |
| --- | --- |
| `<prefix>.ordered.fa` | Retained contigs, ordered by reference chromosome and position. |
| `<prefix>.contig_assignments.tsv` | One row per assembly contig with final status and assignment metrics. |
| `<prefix>.contig_ref_matches.tsv` | One row per contig-reference match before final assignment. |
| `<prefix>.chromosome_summary.tsv` | One row per reference sequence with ordered contig lists and covered reference bp. |
| `<prefix>.run_summary.txt` | Inputs, thresholds, output paths, and status counts. |

### `chromo sort` Parameters

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

### `chromo sort` Status Values

`kept`: written to the ordered FASTA.

`no_alignment`: no usable `show-coords` rows were found for this contig.

`below_min_aligned_bp`: best match did not meet `--min-aligned-bp`.

`below_min_query_cov`: best match did not meet `--min-query-cov`.

`ambiguous_ref_match`: the best chromosome did not dominate the contig's total
matched bp enough to pass `--min-best-ref-share`.

`duplicate_overlap`: the contig passed the basic match thresholds, but better
contigs on the same reference chromosome already covered nearly all of its
reference span. These contigs are not written to the ordered FASTA.

### Reasoning Behind `chromo sort`

#### Use `show-coords`, Not `show-tiling`

`show-tiling` can be useful, but it is another derived representation and is not
always produced. `show-coords` from a filtered delta contains the required
information: reference name, query name, reference coordinates, query
coordinates, alignment length, sequence length, percent identity, and strand.
Using `show-coords` keeps the workflow closer to standard MUMmer output and
avoids an extra preprocessing step.

#### Merge Intervals Before Calculating Coverage

MUMmer can report overlapping rows for the same contig-reference pair. Summing
raw `LEN 2` values can produce apparent coverage greater than 100 percent. By
merging query intervals and reference intervals first, ChromoSort estimates
coverage from the unique aligned span instead of from repeated rows.

#### Assign Contigs by Best Reference Share

Many genomes contain repeats, paralogous regions, translocations, or retained
haplotigs. A contig may have alignments to more than one chromosome. ChromoSort
chooses the chromosome with the largest merged query-aligned bp, then requires
that match to represent a configurable share of all matched bp. This keeps clear
placements and flags ambiguous ones.

#### Filter Duplicate Overlaps After Assignment

The first assignment pass asks, "Does this contig have a good placement?" The
overlap pass asks, "Does this contig add useful new reference coverage beyond
better contigs already kept?" This second question is important for assemblies
that include short duplicate fragments, alternate haplotigs, or local redundant
contigs. Rejected contigs are marked `duplicate_overlap`, with novel coverage
and the strongest overlapping kept contig reported.

#### Preserve Full Contigs

`chromo sort` does not trim contigs to alignment spans. It writes the full input
contig sequence because unaligned ends may be real assembly sequence. The output
is an ordered contig FASTA, not a hard reference-guided reconstruction.

#### Make Orientation Optional

`--orient-to-reference` reverse-complements retained contigs whose dominant
alignment is on the reverse strand. This is helpful for downstream inspection
and plotting. It is optional because some workflows prefer to preserve original
assembly orientation exactly.

### Batch Sorting Example

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

## chromo fix

Use `chromo fix` when the goal is to split chimeric or structurally inconsistent
assembly contigs into reference-labeled pieces.

### What `chromo fix` Does

For each requested or auto-detected contig, `chromo fix`:

1. Reads passing `show-coords` alignment segments for that contig.
2. Sorts those segments by query-coordinate order along the assembly contig.
3. Merges nearby neighboring rows that map to the same reference sequence and
   orientation.
4. Looks for reference transitions along the contig.
5. Looks for orientation transitions along the contig, including same-reference
   inversion blocks.
6. Places breakpoints halfway between neighboring alignment blocks.
7. Replaces the original contig with two or more pieces in the output FASTA.
8. Writes a TSV report with slice coordinates, reference labels, orientation,
   identity, and split status.

By default, unrequested contigs are copied unchanged, producing a full fixed
assembly FASTA. Use `--pieces-only` to write only the split pieces.

### Run `chromo fix` With User-Nominated Contigs

```bash
chromo fix \
  --assembly-fasta assembly.fa \
  --coords mummer/sample.coords \
  --contigs contig_04 contig_12 \
  --output-fasta results/sample.fixed.fa \
  --report results/sample.fixed_contigs.tsv
```

This mode is intentionally user-directed. Provide contigs that you already
suspect are chimeric, usually after inspecting dot plots, assignment reports, or
other QC evidence.

### Run `chromo fix` With Auto Detection

```bash
chromo fix \
  --assembly-fasta assembly.fa \
  --coords mummer/sample.coords \
  --auto \
  --output-fasta results/sample.fixed.fa \
  --report results/sample.fixed_contigs.tsv
```

`--auto` scans all contigs for passing alignment blocks that change reference
sequence or orientation. It can detect interchromosomal chimeras and
same-reference inversion blocks after `--min-segment-bp`, `--min-segment-idy`,
and `--max-merge-gap` have been tuned for the assembly.

### `chromo fix` Outputs

| Output | Description |
| --- | --- |
| `--output-fasta` | Full fixed assembly FASTA by default, with split pieces replacing fixed contigs. |
| `--report` | TSV report describing split pieces and unsplit requested contigs. |

The report includes original contig name, split status, new contig name,
dominant reference, slice coordinates, alignment coordinates, orientation,
reverse-complement status, identity, segment count, and the reason for the
decision.

### `chromo fix` Naming

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
the ordered query blocks require. For example, if `contig_12` has a small
`chrom04` block, a larger `chrom05` block, and another `chrom04` block, it
becomes:

```text
chrom04-contig_12-a
chrom05-contig_12-b
chrom04-contig_12-c
```

The same naming pattern is used for inversions. A contig with a large inverted
block in the middle of a `chrom06` match might become:

```text
chrom06-contig_21-a
chrom06-contig_21-b
chrom06-contig_21-c
```

The report records each piece's orientation so the inverted block is explicit.
The reference names and contig names are not hard-coded. Whatever identifiers
appear in your FASTA and MUMmer output are used. Change the separator with
`--name-separator`.

### `chromo fix` Parameters

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `--contigs` | none | Space-separated names of contigs to inspect and split. |
| `--contigs-file` | none | Optional file with one contig name per line. |
| `--auto` | off | Scan all contigs and split those with reference or orientation transitions. |
| `--min-segment-bp` | `10000` | Minimum alignment segment length used to infer split blocks. |
| `--min-segment-idy` | `0.0` | Minimum percent identity for split-informing alignment rows. |
| `--max-merge-gap` | `1000` | Merge nearby same-reference rows separated by this many query bp or less. |
| `--min-piece-bp` | `1` | Do not emit split pieces shorter than this length. |
| `--orient-to-reference` | off | Reverse-complement split pieces from reverse-strand blocks. |
| `--pieces-only` | off | Write only split pieces instead of a full fixed assembly FASTA. |

### Reasoning Behind `chromo fix`

#### User-Nominated Splitting By Default

Cutting contigs is a stronger intervention than ordering contigs. A reference
transition can reflect a real assembly chimera, but it can also reflect
structural variation, assembly graph complexity, misassembly in the reference,
or poor alignment around repeats. Requiring an explicit contig list keeps this
step auditable: ChromoSort proposes breakpoints from MUMmer coordinates, but the
user decides which contigs are appropriate to cut.

#### Auto Detection As An Opt-In Workflow

`--auto` is useful after the alignment filters have been tuned. It scans for
both chromosome transitions and orientation transitions, so it can detect
classical interchromosomal chimeras as well as same-reference inversion blocks.
It is opt-in because automatic contig cutting should be reviewed carefully.

#### Breakpoints Between Alignment Blocks

`chromo fix` places breakpoints halfway between neighboring query-ordered
alignment blocks. When blocks are adjacent, the breakpoint lands at the
alignment boundary. When there is an unaligned gap, the gap is divided between
the neighboring pieces instead of being discarded.

#### Synthetic Test Cases

The synthetic test data under `tests/data/chimeric` includes four fix cases:
one contig split roughly half-and-half between two reference chromosomes, one
with 25 percent of its sequence matching one chromosome and 75 percent matching
another, one with a large inverted block in the middle, and one with an inverted
block at the end.

## Development

```bash
mamba env create -f environment.yml
mamba activate chromosort
pytest
```

The tests are also compatible with the Python standard-library runner:

```bash
python -m unittest discover -s tests -v
```

The tests use tiny synthetic FASTA and coords files under `tests/data`. They do
not require running MUMmer; MUMmer is included in the environment so users can
generate real input files.

## Citation

If you use ChromoSort, cite this repository and cite MUMmer for the alignment
files used by the workflow.

## Contact

Please use the GitHub issue tracker for bug reports, feature requests, and
questions:

<https://github.com/rotheconrad/chromosort/issues>

## Funding

Funding information has not yet been specified. Add grant numbers, institutional
support, or project-specific funding statements here when they are available.

## Acknowledgements

ChromoSort depends on the MUMmer ecosystem for whole-genome alignment and
coordinate export. Thanks to the genome assembly and comparative genomics
communities whose workflows motivated transparent reference-guided contig
sorting and splitting tools.

## Version History

| Version | Notes |
| --- | --- |
| `0.1.0` | Initial public package with `chromo sort`, `chromo fix`, duplicate-overlap filtering, user-nominated contig splitting, auto detection, and synthetic tests. |
