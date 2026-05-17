# ChromoSort

Reference-guided genome assembly utilities for sorting contigs, splitting
user-flagged or auto-detected chimeric contigs, and scaffolding final ordered
contigs with N gaps.

ChromoSort provides one command, `chromo`, with four subcommands:

- `chromo sort` assigns assembly contigs to reference chromosomes, removes
  low-value duplicate overlaps, protects strong multi-reference split
  candidates, and writes a reference-ordered FASTA.
- `chromo fix` splits contigs with chromosome transitions or reviewed inversion
  blocks into reference-labeled pieces, either from a user-supplied contig list
  or with conservative `--auto` detection that smooths over small
  SV/INDEL-like noise and rejects implausibly break-heavy auto plans.
- `chromo scaffold` joins final sorted contigs into per-reference scaffold
  records with inferred reference-space N gaps by default, or with a fixed
  user-specified number of Ns.
- `chromo plot` draws PDF/SVG/PNG dot plots from existing MUMmer `show-coords` or
  minimap2 PAF alignments, so each fix/sort/scaffold step can be visually
  reviewed without re-running an aligner just to make a plot.

The sorting, fixing, and plotting workflows use standard MUMmer `show-coords`
output or minimap2 PAF. They are designed for reuse across species and genome
assembly projects. ChromoSort does not polish sequence, call variants, or force
contigs to match a reference. It keeps full sequence pieces and writes TSV/plot
reports so each keep, reject, split, plot, or scaffold-gap decision is
auditable.

## Table of Contents

- [Quick Start](#quick-start)
- [Installation With Mamba](#installation-with-mamba)
- [Creating Input Files With MUMmer](#creating-input-files-with-mummer)
  - [Why These MUMmer Choices?](#why-these-mummer-choices)
- [Creating Input Files With minimap2](#creating-input-files-with-minimap2)
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
- [chromo plot](#chromo-plot)
  - [Run `chromo plot`](#run-chromo-plot)
  - [`chromo plot` Outputs](#chromo-plot-outputs)
  - [`chromo plot` Parameters](#chromo-plot-parameters)
- [chromo scaffold](#chromo-scaffold)
  - [What `chromo scaffold` Does](#what-chromo-scaffold-does)
  - [Run `chromo scaffold` With Inferred Gaps](#run-chromo-scaffold-with-inferred-gaps)
  - [Run `chromo scaffold` With Fixed Gaps](#run-chromo-scaffold-with-fixed-gaps)
  - [`chromo scaffold` Outputs](#chromo-scaffold-outputs)
  - [`chromo scaffold` Parameters](#chromo-scaffold-parameters)
  - [Reasoning Behind `chromo scaffold`](#reasoning-behind-chromo-scaffold)
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
chromo plot --help
chromo scaffold --help
```

Recommended order when raw dot plots show possible misjoined contigs:

```bash
# 1. Fix only reviewed/suspect raw contigs, or run conservative auto mode.
chromo fix \
  --assembly-fasta assembly.fa \
  --coords mummer/raw.coords \
  --contigs suspect_contig_1 suspect_contig_2 \
  --output-fasta results/sample.fixed.fa \
  --report results/sample.fixed_contigs.tsv

# 2. Re-align the fixed FASTA with MUMmer, then sort the fixed assembly.
chromo sort \
  --ref-fasta reference.fa \
  --assembly-fasta results/sample.fixed.fa \
  --coords mummer/fixed.coords \
  --output-prefix results/sample.fixed \
  --orient-to-reference

# 3. Plot from the existing coordinates; no new MUMmer run is needed.
chromo plot \
  --ref-fasta reference.fa \
  --assembly-fasta results/sample.fixed.fa \
  --coords mummer/fixed.coords \
  --assignments results/sample.fixed.contig_assignments.tsv \
  --output-prefix plots/sample.fixed \
  --per-ref
```

Running `chromo fix` before `chromo sort` is safest for suspected misjoins
because sorting is a placement/filtering step, not a splitter. `chromo sort`
now protects strong multi-reference split candidates by default, but reviewed
raw-assembly fixing is still the cleaner workflow when the dot plot already
shows which contigs need attention.

Typical reference-ordering run:

```bash
chromo sort \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.fa \
  --coords mummer/sample.coords \
  --output-prefix results/sample \
  --orient-to-reference
```

The same command accepts minimap2 PAF instead of MUMmer coords:

```bash
chromo sort \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.fa \
  --paf paf/sample.paf \
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

Auto-detect contigs with chromosome transitions:

```bash
chromo fix \
  --assembly-fasta assembly.fa \
  --coords mummer/sample.coords \
  --auto \
  --output-fasta results/sample.fixed.fa \
  --report results/sample.fixed_contigs.tsv
```

Typical scaffolding run after final sorting:

```bash
chromo scaffold \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --output-prefix results/sample
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
- minimap2
- MUMmer 4 (`nucmer`, `delta-filter`, `show-coords`)
- Pillow for optional PNG plot output
- pytest for the test suite
- ChromoSort in editable mode

Legacy command aliases are retained for compatibility:

- `chromosort` is equivalent to `chromo sort`
- `chromosort-fix-contigs` is equivalent to `chromo fix`
- `chromosort-scaffold` is equivalent to `chromo scaffold`

New workflows should use `chromo sort`, `chromo fix`, `chromo plot`, and
`chromo scaffold`.

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

## chromo sort

Use `chromo sort` when the goal is to find the best matched assembly contigs for
each reference chromosome and write them in reference order.

### What `chromo sort` Does

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
6. Removes contigs that duplicate already-kept reference intervals and add
   little or no new reference coverage.
7. Protects flagged split candidates from silent duplicate-overlap removal.
8. Sorts retained contigs by reference FASTA order and reference start.
9. Writes an ordered FASTA with names like `chromosome_contig`.
10. Writes TSV reports that explain every keep/reject decision.

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
| `--coords` | required unless `--paf` | MUMmer `show-coords` alignment file. |
| `--paf` | required unless `--coords` | minimap2 PAF alignment file. |
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

### `chromo sort` Status Values

`kept`: written to the ordered FASTA.

`kept_split_candidate`: written to the ordered FASTA and flagged as a strong
multi-reference candidate for `chromo fix` review. These contigs are protected
from duplicate-overlap removal by default.

`kept_large_alignment`: written to the ordered FASTA because it has a very large
best-reference alignment even though its query coverage is slightly below
`--min-query-cov`.

`no_alignment`: no usable alignment rows were found for this contig.

`below_min_aligned_bp`: best match did not meet `--min-aligned-bp`.

`below_min_query_cov`: best match did not meet `--min-query-cov`.

`ambiguous_ref_match`: the best chromosome did not dominate the contig's total
matched bp enough to pass `--min-best-ref-share`.

`duplicate_overlap`: the contig passed the basic match thresholds, but better
contigs on the same reference chromosome already covered nearly all of its
reference span. These contigs are not written to the ordered FASTA.

### Reasoning Behind `chromo sort`

#### Use Segment Coordinates, Not Tiling Summaries

`show-tiling` can be useful for MUMmer workflows, but it is another derived
representation and is not always produced. `show-coords` from a filtered delta
and minimap2 PAF both contain the required information: reference name, query
name, reference coordinates, query coordinates, alignment length, sequence
length, percent identity, and strand. ChromoSort normalizes either format into
the same internal segment representation before sorting, fixing, or plotting.

#### Merge Intervals Before Calculating Coverage

Whole-genome aligners can report overlapping rows for the same contig-reference
pair. Summing raw row lengths can produce apparent coverage greater than 100
percent. By merging query intervals and reference intervals first, ChromoSort
estimates coverage from the unique aligned span instead of from repeated rows.

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

By default, duplicate filtering uses each contig's broad first-to-last reference
span and requires both the novel-bp and novel-fraction thresholds before
retaining a lower-ranked contig. This is intentionally stricter than exact
alignment-block overlap: whole-genome alignments can be fragmented by repeats,
local variation, or filtering, and short duplicate fragments often land in the
internal gaps of a stronger chromosome-scale contig. Use
`--overlap-mode alignment --novel-ref-criteria either` for the older, more
permissive behavior.

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
4. Collapses adjacent same-reference/orientation runs so ordinary alignment
   gaps do not become breakpoints.
5. Looks for reference transitions along the contig.
6. In user-nominated mode, also looks for orientation transitions along the
   contig, including same-reference inversion blocks.
7. In `--auto` mode, scores candidate breakpoints with a breakpoint-penalty
   segmentation model that smooths over weak local discordance.
8. Ranks auto split plans by evidence and accepts them until the run-level
   `--auto-max-breakpoints` budget is used.
9. Places accepted breakpoints halfway between neighboring alignment blocks.
10. Replaces the original contig with two or more pieces in the output FASTA.
11. Writes a TSV report with slice coordinates, reference labels, orientation,
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
sequence. Candidate contigs are then segmented with a conservative
breakpoint-penalty model: keeping a small discordant block inside a larger
piece has a cost, but adding a breakpoint has a larger fixed penalty. This lets
`chromo fix` smooth over real-looking small SVs, local repeat hits, fragmented
alignments, and INDEL-sized gaps while still splitting strong large-scale
chimeras. Auto mode also reviews same-reference orientation events: simple
contiguous inversions are left unchanged, while complex/nested orientation
events with large overlapping reference spans can be split.

By default, at most four auto breakpoints are accepted across the whole run.
Candidate plans are ranked by evidence, with reference/chromosome transitions
ahead of same-reference complex orientation events. Lower-priority plans that
do not fit the remaining budget are left unchanged and reported for manual
review.

Use `--auto-split-inversions` only when all same-reference orientation changes
should be considered for splitting. Without it, conservative auto mode limits
same-reference auto cuts to complex/nested events because large cultivar
inversions may be biological or reference-relative structure rather than
misjoined contigs.

Use `--auto-sensitive` to recover the earlier behavior that splits every
passing reference or orientation transition. That mode is useful for debugging
and exploratory scans, but it is more likely to introduce extra breakpoints.
Even in sensitive mode, adjacent same-reference/orientation runs are collapsed
before cutting.

### `chromo fix` Outputs

| Output | Description |
| --- | --- |
| `--output-fasta` | Full fixed assembly FASTA by default, with split pieces replacing fixed contigs. |
| `--report` | TSV report describing split pieces and unsplit requested contigs. |

The report includes original contig name, split status, new contig name,
dominant reference, slice coordinates, alignment coordinates, orientation,
reverse-complement status, identity, segment count, and the reason for the
decision. Auto candidates that contain discordant blocks but are not cut are
reported as `not_split_auto_smooth`. Candidates rejected by the auto breakpoint
cap are reported as `not_split_auto_too_many_breakpoints`.

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
appear in your FASTA and alignment output are used. Change the separator with
`--name-separator`.

### `chromo fix` Parameters

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `--coords` | required unless `--paf` | MUMmer `show-coords` alignment file. |
| `--paf` | required unless `--coords` | minimap2 PAF alignment file. |
| `--contigs` | none | Space-separated names of contigs to inspect and split. |
| `--contigs-file` | none | Optional file with one contig name per line. |
| `--auto` | off | Scan all contigs and conservatively split those with strong reference transitions. |
| `--auto-sensitive` | off | Disable auto smoothing and split every passing reference/orientation transition. |
| `--auto-split-inversions` | off | In conservative auto mode, also split same-reference orientation transitions. |
| `--auto-complex-inversions` | on | In conservative auto mode, split complex/nested same-reference orientation events while leaving simple contiguous inversions unchanged. |
| `--min-segment-bp` | `10000` | Minimum alignment segment length used to infer split blocks. |
| `--min-segment-idy` | `0.0` | Minimum percent identity for split-informing alignment rows. |
| `--max-merge-gap` | `1000` | Merge nearby same-reference rows separated by this many query bp or less. |
| `--min-mapq` | `0` | Ignore PAF rows below this MAPQ. Ignored for coords. |
| `--include-secondary-paf` | off | Include PAF rows marked `tp:A:S`; skipped by default. |
| `--min-piece-bp` | `1` | Do not emit split pieces shorter than this length. |
| `--auto-breakpoint-penalty-bp` | `50000` | Identity-weighted aligned bp cost charged for each auto breakpoint. |
| `--auto-min-piece-aligned-bp` | `50000` | Minimum dominant aligned bp required in each auto-split piece. |
| `--auto-min-piece-query-frac` | `0.05` | Minimum query-span fraction required in each auto-split piece. |
| `--auto-complex-inversion-min-piece-aligned-bp` | `1000000` | Minimum dominant aligned bp for pieces used to classify a same-reference orientation event as complex. |
| `--auto-complex-inversion-min-overlap-frac` | `0.50` | Minimum reference-span overlap fraction for classifying a same-reference orientation event as complex. |
| `--auto-max-breakpoints` | `4` | Maximum accepted auto breakpoints across the whole run. Set negative to disable. |
| `--orient-to-reference` | off | Reverse-complement split pieces from reverse-strand blocks. |
| `--pieces-only` | off | Write only split pieces instead of a full fixed assembly FASTA. |

### Reasoning Behind `chromo fix`

#### User-Nominated Splitting By Default

Cutting contigs is a stronger intervention than ordering contigs. A reference
transition can reflect a real assembly chimera, but it can also reflect
structural variation, assembly graph complexity, misassembly in the reference,
or poor alignment around repeats. Requiring an explicit contig list keeps this
step auditable: ChromoSort proposes breakpoints from alignment coordinates, but the
user decides which contigs are appropriate to cut.

#### Auto Detection As An Opt-In Workflow

`--auto` is useful after the alignment filters have been tuned. It prioritizes
chromosome/reference transitions because these are the strongest signal for
misjoined contigs. Same-reference orientation events are handled more carefully:
simple contiguous inversions are ignored by default, while complex/nested
events with large overlapping reference spans can be split. Use
`--auto-split-inversions` to consider all same-reference orientation changes.
Auto fixing is opt-in because automatic contig cutting should be reviewed
carefully.

#### Collapse Same-Target Runs Before Cutting

Whole-genome alignments often contain many neighboring rows for the same
reference and orientation, separated by local gaps, repeats, or small assembly
differences. `chromo fix` collapses adjacent same-target runs before placing
breakpoints. A contig is cut at accepted reference transitions, or at
same-reference orientation transitions only when they are complex enough or
explicitly enabled, not at every ordinary alignment-row boundary.

#### Breakpoint-Penalty Segmentation

Auto mode uses a small dynamic-programming segmentation model. For each contig,
`chromo fix` asks whether the query-ordered alignment blocks are better
explained as one smoothed piece or as multiple pieces separated by breakpoints.

Keeping a discordant block inside a larger piece costs that block's
identity-weighted aligned bp. Adding a breakpoint costs
`--auto-breakpoint-penalty-bp`. A breakpoint is accepted only when the reduction
in discordant support is worth paying the penalty and every resulting piece has
enough dominant aligned support and spans at least 5% of the contig by default.
This makes the default behavior breakpoint-averse: small terminal off-target
blocks, small inversions, short transposed/repeat-like hits, fragmented
same-chromosome alignments, and INDEL-sized gaps are smoothed over instead of
cut.

Auto mode also caps the number of accepted breakpoints across the run with
`--auto-max-breakpoints`. The default of four is meant as a practical guardrail
for soybean-scale samples: if many candidate plans compete for that budget,
ChromoFix keeps the highest-priority/highest-score plans and reports the rest
as `not_split_auto_too_many_breakpoints` for manual dot plot review.

#### Breakpoints Between Alignment Blocks

After segmentation, `chromo fix` places accepted breakpoints halfway between
neighboring query-ordered alignment blocks. When blocks are adjacent, the
breakpoint lands at the alignment boundary. When there is an unaligned gap, the
gap is divided between the neighboring pieces instead of being discarded.

#### Synthetic Test Cases

The synthetic test data under `tests/data/chimeric` includes direct fix cases:
one contig split roughly half-and-half between two reference chromosomes, one
with 25 percent of its sequence matching one chromosome and 75 percent matching
another, one with a large inverted block in the middle, and one with an inverted
block at the end.

The noisier benchmark under `tests/data/noisy_fix` adds INDEL-like gaps, a small
local inversion, a short repeat-like hit to another chromosome, true large-scale
chimeras, a complex chimera with internal gaps, and terminal/internal inversion
cases. The expected behavior is conservative: split the large-scale chromosome
transition patterns, split only complex same-reference orientation events by
default, split all strong inversions only when inversion splitting is enabled,
and report weaker discordance as `not_split_auto_smooth`.

## chromo plot

Use `chromo plot` when you already have a MUMmer coords or minimap2 PAF file
and want a visual check without running `mummerplot` or re-aligning only for
the plot. It writes PDF by default and can also write SVG or PNG. Forward-strand
alignments are blue and reverse-strand alignments are red.

### Run `chromo plot`

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

### `chromo plot` Outputs

| Output | Description |
| --- | --- |
| `<prefix>.pdf` | Whole-genome PDF dot plot by default. |
| `<prefix>.svg` | Whole-genome SVG dot plot when `--formats svg` is set. |
| `<prefix>.png` | Whole-genome PNG dot plot when `--formats png` is set. |
| `<prefix>.<ref>.<format>` | Per-reference plots when `--per-ref` is set. |

### `chromo plot` Parameters

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

## chromo scaffold

Use `chromo scaffold` when the final sorted and filtered contigs look good and
you want one FASTA record per reference chromosome or linkage group.

### What `chromo scaffold` Does

Given a final `chromo sort` ordered FASTA and its matching
`<prefix>.contig_assignments.tsv`, `chromo scaffold`:

1. Reads kept contigs from the assignment report.
2. Reads the ordered FASTA records by their `new_name` values.
3. Groups contigs by assigned reference sequence in ordered FASTA order.
4. Joins neighboring contigs with N gaps.
5. Infers gap length from adjacent reference coordinates by default.
6. Optionally uses a fixed user-provided number of Ns between every neighboring
   contig.
7. Writes scaffold FASTA, gap report, scaffold summary, and run summary files.

The intended input is the final ordered FASTA from the same `chromo sort` run as
the assignment report. If you run `chromo fix`, re-run `chromo sort` on the
fixed assembly before scaffolding so the coordinates and FASTA names match the
final contigs.

### Run `chromo scaffold` With Inferred Gaps

```bash
chromo scaffold \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --output-prefix results/sample
```

For adjacent contigs on the same reference, inferred gaps are calculated as:

```text
next_ref_start - previous_ref_end - 1
```

Negative values, which indicate overlapping reference spans, are written as
zero-length gaps in the FASTA and reported in the gap TSV.

### Run `chromo scaffold` With Fixed Gaps

```bash
chromo scaffold \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --output-prefix results/sample.fixed100 \
  --fixed-gap-bp 100
```

Fixed-gap mode ignores inferred gap length for FASTA construction and inserts
the requested number of Ns between every neighboring contig on the same
scaffold. The report still records the raw inferred gap for comparison.

### `chromo scaffold` Outputs

| Output | Description |
| --- | --- |
| `<prefix>.scaffold.fa` | One FASTA record per assigned reference sequence, with ordered contigs joined by Ns. |
| `<prefix>.scaffold_gaps.tsv` | One row per inserted gap with flanking contigs, inferred gap, written gap, and gap mode. |
| `<prefix>.scaffold_summary.tsv` | One row per scaffold with contig count, scaffold length, sequence bp, gap bp, and ordered contig list. |
| `<prefix>.run_summary.txt` | Inputs, gap model, output paths, and total scaffold counts. |

### `chromo scaffold` Parameters

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `--ordered-fasta` | required | Final ordered FASTA from `chromo sort`. |
| `--assignments` | required | Matching `<prefix>.contig_assignments.tsv` report from `chromo sort`. |
| `--output-prefix` | required | Prefix for scaffold FASTA and reports. |
| `--fixed-gap-bp` | none | Insert this many Ns between neighboring contigs instead of inferred gaps. |
| `--simple-headers` | off | Write scaffold FASTA headers containing only the scaffold ID. |

### Reasoning Behind `chromo scaffold`

#### Require the Assignment Report

The ordered FASTA contains sequence, but the assignment report contains the
reference coordinates needed to infer gap lengths. Requiring both files keeps
scaffolding explicit and prevents guessing chromosome names from FASTA IDs,
which can be fragile when contig names contain separators.

#### Infer Gaps Conservatively

Inferred gaps are based only on adjacent retained contigs on the same reference
sequence. ChromoSort does not add leading or trailing Ns for chromosome ends,
and overlapping reference spans become zero-length FASTA gaps rather than
negative gaps. The raw inferred value is still reported so users can inspect
overlap or coordinate oddities.

#### Keep Fixed Gaps Available

Some downstream tools and submission workflows prefer a constant gap size. The
`--fixed-gap-bp` option supports that convention while preserving the inferred
gap estimate in the report for transparency.

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

The tests use tiny synthetic FASTA, coords, and PAF files under `tests/data`.
They do not require running MUMmer or minimap2; both are included in the
environment so users can generate real input files.

## Citation

If you use ChromoSort, cite this repository and cite MUMmer or minimap2 for the
alignment files used by the workflow.

## Contact

Please use the GitHub issue tracker for bug reports, feature requests, and
questions:

<https://github.com/rotheconrad/chromosort/issues>

## Funding

Funding information has not yet been specified. Add grant numbers, institutional
support, or project-specific funding statements here when they are available.

## Acknowledgements

ChromoSort can consume MUMmer and minimap2 whole-genome alignments. Thanks to
the genome assembly and comparative genomics communities whose workflows
motivated transparent reference-guided contig sorting, splitting, plotting, and
scaffolding tools.

## Version History

| Version | Notes |
| --- | --- |
| `0.2.1` | Tightened `chromo sort` duplicate filtering for contaminated/alternate-fragment assemblies by using span-based overlap by default, requiring both novel coverage thresholds, rescuing very large near-threshold alignments, and letting split candidates protect their secondary reference spans. |
| `0.2.0` | Added minimap2 PAF input for `chromo sort` and `chromo fix`, plus `chromo plot` PDF/SVG/PNG dot plots for coords/PAF with optional assignment-report query ordering. |
| `0.1.2` | Raised the default auto-split query-span support threshold to 5% so small terminal off-target blocks are reported for review instead of being cut automatically. |
| `0.1.1` | Tightened `chromo fix` breakpoint placement by collapsing adjacent same-reference/orientation runs, added complex same-reference orientation detection, added a run-level auto breakpoint budget, protected strong multi-reference split candidates during `chromo sort`, and documented the fix-before-sort workflow for suspected misjoins. |
| `0.1.0` | Initial public package with `chromo sort`, `chromo fix`, `chromo scaffold`, duplicate-overlap filtering, user-nominated contig splitting, conservative auto smoothing, inferred/fixed-gap scaffolding, and synthetic tests. |
