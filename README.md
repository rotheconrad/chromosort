# ChromoSort

Reference-guided genome assembly utilities for sorting contigs, splitting
selected contigs or all detected chimeric contigs, and scaffolding final
ordered contigs with N gaps.

ChromoSort provides one command, `chromo`, with five subcommands:

- `chromo sort` assigns assembly contigs to reference chromosomes, removes
  low-value duplicate overlaps, labels and can rescue terminal overlaps,
  protects strong multi-reference split candidates, and writes a
  reference-ordered FASTA.
- `chromo fix` splits contigs with chromosome transitions or reviewed inversion
  blocks into reference-labeled pieces. Use `--contigs` to inspect a reviewed
  subset or `--all` to scan every contig; both scopes use the same `--mode`
  planner and per-contig breakpoint guardrails.
- `chromo cut` cuts reviewed contigs at exact user-provided positions and writes
  a new FASTA plus a TSV report. It is useful when the breakpoint is already
  known from manual dot-plot review.
- `chromo scaffold` joins final sorted contigs into per-reference scaffold
  records with inferred reference-space N gaps by default, reports negative-gap
  overlaps, and can optionally trim reviewed terminal overlaps.
- `chromo plot` draws PDF/SVG/PNG dot plots from existing MUMmer `show-coords` or
  minimap2 PAF alignments, so each fix/sort/scaffold step can be visually
  reviewed without re-running an aligner just to make a plot.

The sorting, fixing, and plotting workflows use standard MUMmer `show-coords`
output or minimap2 PAF. Manual cuts do not require alignments. These commands
are designed for reuse across species and genome assembly projects. ChromoSort
does not polish sequence, call variants, or force contigs to match a reference.
It keeps full sequence pieces and writes TSV/plot reports so each keep, reject,
split, manual cut, plot, or scaffold-gap decision is auditable.

## Table of Contents

- [Quick Start](#quick-start)
- [Handling Overlapping Contigs](#handling-overlapping-contigs)
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
- [chromo cut](#chromo-cut)
  - [Run `chromo cut`](#run-chromo-cut)
  - [`chromo cut` Outputs](#chromo-cut-outputs)
  - [`chromo cut` Parameters](#chromo-cut-parameters)
- [chromo fix](#chromo-fix)
  - [What `chromo fix` Does](#what-chromo-fix-does)
  - [Run `chromo fix` With Selected Contigs](#run-chromo-fix-with-selected-contigs)
  - [Run `chromo fix` Across All Contigs](#run-chromo-fix-across-all-contigs)
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
  - [Run `chromo scaffold` With Overlap Trimming](#run-chromo-scaffold-with-overlap-trimming)
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
chromo cut --help
chromo plot --help
chromo scaffold --help
```

Recommended order when raw dot plots show possible misjoined contigs:

```bash
# 1. Fix only reviewed/suspect raw contigs, or use --all to scan every contig.
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

Manual cut at reviewed positions:

```bash
chromo cut \
  --assembly-fasta assembly.fa \
  --cut contig_1:234567,450000 \
  --cut contig_2:10000 \
  --output-fasta results/sample.cut.fa \
  --report results/sample.cut_contigs.tsv
```

Scan all contigs with the same conservative planner:

```bash
chromo fix \
  --assembly-fasta assembly.fa \
  --coords mummer/sample.coords \
  --all \
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

## Handling Overlapping Contigs

Large contigs sometimes have strong reference support but overlap at their ends.
This can happen with alternate graph paths, small assembly redundancies, or real
terminal dovetails that another scaffolding workflow might trim or merge. In
ChromoSort, overlap handling is deliberately split across commands so the
sequence-changing step stays explicit.

`chromo sort` assigns and filters contigs in reference space. It lets stronger
contigs claim reference intervals first, then asks whether lower-ranked contigs
add useful novel reference coverage. Fully contained or internal redundant
fragments are still reported as `duplicate_overlap`. Dovetail-like one-sided
overlaps are now classified separately as `terminal_overlap`; if the contig is
kept, its status is `kept_terminal_overlap`. A mostly overlapping contig can also
be rescued when its one-sided extension passes
`--min-terminal-extension-bp` and `--min-terminal-extension-frac`.

`chromo fix` does not resolve overlap between two separate contigs. It only
splits within a selected contig when query-ordered alignment blocks support a
reference or eligible orientation transition. If two already-separate contigs
overlap at their ends on the same reference, `fix` will not trim, merge, or
deduplicate them.

`chromo scaffold` joins the final sorted contigs. In the default
`--overlap-policy zero-gap` mode, adjacent negative reference gaps are written as
zero-length gaps: no Ns are inserted and neither contig is trimmed. The raw
negative inferred gap, overlap bp, overlap class, overlap fractions, policy, and
action are reported in `<prefix>.scaffold_gaps.tsv`, and scaffold-level overlap
and trimming totals are reported in `<prefix>.scaffold_summary.tsv`.

When you want sequence surgery at scaffolding time, choose an explicit overlap
policy:

```bash
# Current conservative behavior, with stderr warnings for negative gaps.
chromo scaffold \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --output-prefix results/sample.warn \
  --overlap-policy warn

# Trim the right contig by the reference-inferred terminal overlap.
chromo scaffold \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --output-prefix results/sample.trim_ref \
  --overlap-policy trim-reference

# Trim only when the left suffix and right prefix confirm the overlap sequence.
chromo scaffold \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --output-prefix results/sample.trim_seq \
  --overlap-policy trim-sequence \
  --trim-sequence-min-identity 0.98
```

The safest default is still not to trim: overlapping reference spans can reflect
assembly redundancy, true variation, reference differences, or alignment
artifacts. The new reports make the case easy to find, and the trimming policies
make the intervention deliberate when a reviewed dataset needs it.

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

New workflows should use `chromo sort`, `chromo fix`, `chromo cut`,
`chromo plot`, and `chromo scaffold`.

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
6. Classifies overlap shape against already-kept reference intervals.
7. Removes contained/internal duplicate overlaps that add little or no new
   reference coverage.
8. Keeps or rescues terminal overlaps when they contribute enough one-sided
   extension.
9. Protects flagged split candidates from silent duplicate-overlap removal.
10. Sorts retained contigs by reference FASTA order and reference start.
11. Writes an ordered FASTA with names like `chromosome_contig`.
12. Writes TSV reports that explain every keep/reject decision.

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

### `chromo sort` Status Values

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

#### Preserve Full Contigs

`chromo sort` does not trim contigs to alignment spans. It writes the full input
contig sequence because unaligned ends and terminal extensions may be real
assembly sequence. The output is an ordered contig FASTA, not a hard
reference-guided reconstruction. Optional overlap trimming is handled only by
`chromo scaffold`, where the sequence-changing policy is explicit.

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

## chromo cut

Use `chromo cut` when you already know the exact breakpoint position for one or
more contigs and want a transparent FASTA edit without running the alignment
planner in `chromo fix`.

Positions are 1-based and mean "cut after this base." For example, cutting a
10 bp contig at position 4 emits bases 1-4 and 5-10. Terminal cuts are rejected
because they would create empty pieces.

### Run `chromo cut`

Repeat `--cut` for multiple contigs, and use commas for multiple positions on
the same contig:

```bash
chromo cut \
  --assembly-fasta assembly.fa \
  --cut contig_1:234567,450000 \
  --cut contig_2:10000 \
  --output-fasta results/sample.cut.fa \
  --report results/sample.cut_contigs.tsv
```

For a single contig, the original explicit form is also supported:

```bash
chromo cut \
  --assembly-fasta assembly.fa \
  --contig contig_1 \
  --pos 234567 450000 \
  --output-fasta results/sample.cut.fa \
  --report results/sample.cut_contigs.tsv
```

For batch review, provide a simple cuts file:

```text
contig	position
contig_1	234567	450000
contig_2	10000
```

Then run:

```bash
chromo cut \
  --assembly-fasta assembly.fa \
  --cuts-file reviewed_cuts.tsv \
  --output-fasta results/sample.cut.fa \
  --report results/sample.cut_contigs.tsv
```

### `chromo cut` Outputs

| Output | Description |
| --- | --- |
| `--output-fasta` | Full assembly FASTA with requested contigs replaced by cut pieces. Uncut contigs are copied unchanged. |
| `--report` | TSV report describing every emitted cut piece, including original contig, new contig, slice coordinates, piece length, and cut positions. |

Cut pieces are named `CONTIG_cut001`, `CONTIG_cut002`, and so on by default.
Change the inserted text with `--name-separator`, or use `--simple-headers` to
write only the new FASTA IDs in cut-piece headers.

### `chromo cut` Parameters

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `--assembly-fasta` | required | Assembly FASTA containing the contigs to cut. |
| `--assembly-fai` | auto | Optional FASTA index used for length validation. Defaults to `<assembly-fasta>.fai` when present. |
| `--cut` | none | Reviewed cut as `CONTIG:POS[,POS...]`; may be repeated. |
| `--cuts-file` | none | Optional file with contig and one or more cut positions per line. |
| `--contig` / `--pos` | none | Convenience form for one contig with one or more positions. |
| `--output-fasta` | required | New FASTA path. The input FASTA is never modified. |
| `--report` | required | TSV report path for the cut pieces. |
| `--min-piece-bp` | `1` | Reject cut plans that create pieces shorter than this length. |
| `--name-separator` | `_cut` | Text inserted before the numeric suffix. |
| `--simple-headers` | off | Write cut-piece FASTA headers containing only the new sequence ID. |

## chromo fix

Use `chromo fix` when the goal is to split chimeric or structurally inconsistent
assembly contigs into reference-labeled pieces.

### What `chromo fix` Does

For each selected contig, `chromo fix`:

1. Reads passing `show-coords` or PAF alignment segments for that contig.
2. Sorts those segments by query-coordinate order along the assembly contig.
3. Merges nearby neighboring rows that map to the same reference sequence and
   orientation.
4. Collapses adjacent same-reference/orientation runs so ordinary alignment
   gaps do not become breakpoints.
5. Applies the selected `--mode` to decide which reference/orientation
   transitions are eligible.
6. In smoothed modes, scores candidate breakpoints with a breakpoint-penalty
   segmentation model that filters weak local discordance.
7. Rejects any plan exceeding `--max-breakpoints-per-contig`.
8. Places accepted breakpoints halfway between neighboring alignment blocks.
9. Replaces the original contig with two or more pieces in the output FASTA.
10. Writes a TSV report with slice coordinates, reference labels, orientation,
   identity, and split status.

By default, unrequested contigs are copied unchanged, producing a full fixed
assembly FASTA. Use `--pieces-only` to write only the split pieces.

`chromo fix` is not a cross-contig overlap resolver. It does not merge two
separate contigs, trim a terminal overlap between neighboring contigs, or choose
one contig over another. Use `chromo sort` reports to identify duplicate or
terminal overlap relationships, then use `chromo scaffold --overlap-policy` only
when you want an explicit scaffolding-time trim.

### Run `chromo fix` With Selected Contigs

```bash
chromo fix \
  --assembly-fasta assembly.fa \
  --coords mummer/sample.coords \
  --contigs contig_04 contig_12 \
  --output-fasta results/sample.fixed.fa \
  --report results/sample.fixed_contigs.tsv
```

`--contigs` and `--contigs-file` only choose which contigs to inspect. They do
not switch to a different splitting algorithm. By default, selected contigs use
the same conservative smoothing and breakpoint penalties as `--all`, which is
useful when you want Benning-style targeted fixes without allowing off-target
contigs to receive a break.

### Run `chromo fix` Across All Contigs

```bash
chromo fix \
  --assembly-fasta assembly.fa \
  --coords mummer/sample.coords \
  --all \
  --output-fasta results/sample.fixed.fa \
  --report results/sample.fixed_contigs.tsv
```

`--all` scans every contig and queues those with passing split signals for the
same planner used by `--contigs`. The default `--mode conservative` smooths over
real-looking small SVs, local repeat hits, fragmented alignments, and INDEL-sized
gaps while still splitting strong large-scale chimeras.

Use `--mode chromosome` when only reference/chromosome transitions should be
eligible. Use `--mode comprehensive` when all same-reference orientation changes
should also be considered. Use `--mode sensitive` for the earlier direct behavior
that cuts every passing reference/orientation transition after collapsing
adjacent same-target runs.

### `chromo fix` Outputs

| Output | Description |
| --- | --- |
| `--output-fasta` | Full fixed assembly FASTA by default, with split pieces replacing fixed contigs. |
| `--report` | TSV report describing split pieces and unsplit requested contigs. |

The report includes original contig name, split status, new contig name,
dominant reference, slice coordinates, alignment coordinates, orientation,
reverse-complement status, identity, segment count, and the reason for the
decision. Candidates that contain discordant blocks but are not cut by
breakpoint smoothing are reported as `not_split_smooth`. Candidates rejected by
the per-contig breakpoint cap are reported as `not_split_too_many_breakpoints`.

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
| `--all` | off | Inspect all contigs with passing split signals. |
| `--mode` | `conservative` | Planner used for `--contigs`, `--contigs-file`, or `--all`: `conservative`, `chromosome`, `comprehensive`, or `sensitive`. |
| `--min-segment-bp` | `10000` | Minimum alignment segment length used to infer split blocks. |
| `--min-segment-idy` | `0.0` | Minimum percent identity for split-informing alignment rows. |
| `--max-merge-gap` | `1000` | Merge nearby same-reference rows separated by this many query bp or less. |
| `--min-mapq` | `0` | Ignore PAF rows below this MAPQ. Ignored for coords. |
| `--include-secondary-paf` | off | Include PAF rows marked `tp:A:S`; skipped by default. |
| `--min-piece-bp` | `1` | Do not emit split pieces shorter than this length. |
| `--breakpoint-penalty-bp` | `50000` | Identity-weighted aligned bp cost charged for each smoothed breakpoint. |
| `--min-piece-aligned-bp` | `50000` | Minimum dominant aligned bp required in each smoothed split piece. |
| `--min-piece-query-frac` | `0.05` | Minimum query-span fraction required in each smoothed split piece. |
| `--complex-inversion-min-piece-aligned-bp` | `1000000` | Minimum dominant aligned bp for pieces used to classify a same-reference orientation event as complex. |
| `--complex-inversion-min-overlap-frac` | `0.50` | Minimum reference-span overlap fraction for classifying a same-reference orientation event as complex. |
| `--max-breakpoints-per-contig` | `4` | Maximum accepted breakpoints per contig. Set negative to disable. |
| `--orient-to-reference` | off | Reverse-complement split pieces from reverse-strand blocks. |
| `--pieces-only` | off | Write only split pieces instead of a full fixed assembly FASTA. |

### Reasoning Behind `chromo fix`

#### Scope Is Separate From Mode

Cutting contigs is a stronger intervention than ordering contigs. A reference
transition can reflect a real assembly chimera, but it can also reflect
structural variation, assembly graph complexity, misassembly in the reference,
or poor alignment around repeats. `--contigs` and `--contigs-file` keep the
operation auditable by limiting which contigs can receive a break. `--all` uses
the same planner across the whole assembly after the alignment filters have
been tuned.

`--mode conservative` prioritizes chromosome/reference transitions because
these are the strongest signal for misjoined contigs. Same-reference orientation
events are handled more carefully: simple contiguous inversions are ignored by
default, while complex/nested events with large overlapping reference spans can
be split. `--mode comprehensive` considers all same-reference orientation
changes. `--mode chromosome` ignores same-reference orientation changes.
`--mode sensitive` disables breakpoint-penalty smoothing and is mainly useful
for debugging or intentionally aggressive exploratory scans.

#### Collapse Same-Target Runs Before Cutting

Whole-genome alignments often contain many neighboring rows for the same
reference and orientation, separated by local gaps, repeats, or small assembly
differences. `chromo fix` collapses adjacent same-target runs before placing
breakpoints. A contig is cut at accepted reference transitions, or at
same-reference orientation transitions only when they are complex enough or
explicitly enabled, not at every ordinary alignment-row boundary.

#### Breakpoint-Penalty Segmentation

Smoothed modes use a small dynamic-programming segmentation model. For each contig,
`chromo fix` asks whether the query-ordered alignment blocks are better
explained as one smoothed piece or as multiple pieces separated by breakpoints.

Keeping a discordant block inside a larger piece costs that block's
identity-weighted aligned bp. Adding a breakpoint costs
`--breakpoint-penalty-bp`. A breakpoint is accepted only when the reduction
in discordant support is worth paying the penalty and every resulting piece has
enough dominant aligned support and spans at least 5% of the contig by default.
This makes the default behavior breakpoint-averse: small terminal off-target
blocks, small inversions, short transposed/repeat-like hits, fragmented
same-chromosome alignments, and INDEL-sized gaps are smoothed over instead of
cut.

`--max-breakpoints-per-contig` caps accepted breakpoints independently for each
contig. The default of four is meant as a practical guardrail for soybean-scale
samples: a contig that appears to need many breaks is more likely to need manual
dot plot review than automatic sequence surgery. Those plans are reported as
`not_split_too_many_breakpoints`.

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
default, split all strong inversions only in `--mode comprehensive`, and report
weaker discordance as `not_split_smooth`.

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
6. Reports negative inferred gaps as adjacent reference overlaps.
7. Optionally trims reviewed terminal overlaps according to `--overlap-policy`.
8. Optionally uses a fixed user-provided number of Ns between every neighboring
   contig.
9. Writes scaffold FASTA, gap report, scaffold summary, and run summary files.

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
zero-length gaps in the FASTA and reported in the gap TSV by default. Use
`--overlap-policy warn` to keep the same FASTA behavior while emitting stderr
warnings, `--overlap-policy trim-reference` to trim the right contig by the
reference-inferred terminal overlap, or `--overlap-policy trim-sequence` to trim
only when the left suffix and right prefix confirm the overlap sequence at
`--trim-sequence-min-identity`.

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

### Run `chromo scaffold` With Overlap Trimming

```bash
chromo scaffold \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --output-prefix results/sample.trim_seq \
  --overlap-policy trim-sequence
```

`trim-reference` removes the reference-inferred overlap from the left side of
the right contig when the adjacent reference spans form a terminal overlap.
`trim-sequence` is more conservative: it trims only when the left contig suffix
and right contig prefix match across the inferred overlap with at least
`--trim-sequence-min-identity` identity. Non-terminal overlaps are reported but
not trimmed by either trimming policy.

### `chromo scaffold` Outputs

| Output | Description |
| --- | --- |
| `<prefix>.scaffold.fa` | One FASTA record per assigned reference sequence, with ordered contigs joined by Ns. |
| `<prefix>.scaffold_gaps.tsv` | One row per inserted gap with flanking contigs, inferred gap, written gap, overlap bp/class/fractions, overlap policy/action, trimmed bp, and sequence-overlap identity when checked. |
| `<prefix>.scaffold_summary.tsv` | One row per scaffold with contig count, scaffold length, sequence bp, gap bp, overlap totals, trimming totals, and ordered contig list. |
| `<prefix>.run_summary.txt` | Inputs, gap model, output paths, and total scaffold counts. |

### `chromo scaffold` Parameters

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `--ordered-fasta` | required | Final ordered FASTA from `chromo sort`. |
| `--assignments` | required | Matching `<prefix>.contig_assignments.tsv` report from `chromo sort`. |
| `--output-prefix` | required | Prefix for scaffold FASTA and reports. |
| `--fixed-gap-bp` | none | Insert this many Ns between neighboring contigs instead of inferred gaps. |
| `--overlap-policy` | `zero-gap` | Handling for negative inferred gaps: `zero-gap`, `warn`, `trim-reference`, or `trim-sequence`. |
| `--trim-sequence-min-identity` | `0.98` | Minimum suffix/prefix identity required by `--overlap-policy trim-sequence`. |
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
negative gaps by default. The raw inferred value, overlap bp, overlap class,
overlap fractions, policy, action, and trimming amount are reported so users can
inspect overlap or coordinate oddities.

#### Trim Only When Asked

Overlaps are not trimmed by default because a negative reference gap can mean
several different things: a real dovetail, retained alternate sequence,
reference/assembly structural difference, or an alignment artifact. The
`trim-reference` policy trims only terminal overlaps and trusts the reference
coordinate estimate. The `trim-sequence` policy trims only terminal overlaps
whose left suffix and right prefix agree at high identity. Contained/internal
overlaps are reported for review rather than trimmed automatically.

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
| `0.2.4` | Added `chromo cut` for exact reviewed breakpoint cuts, with repeatable `--cut CONTIG:POS[,POS...]`, single-contig `--contig/--pos`, batch `--cuts-file`, cut-piece FASTA output, and an audit TSV report. |
| `0.2.3` | Added explicit terminal-overlap classification/rescue in `chromo sort`, richer scaffold overlap reporting, and `chromo scaffold --overlap-policy` modes for warn-only, reference-coordinate trimming, and sequence-confirmed trimming. |
| `0.2.2` | Reworked `chromo fix` so `--contigs`/`--contigs-file` only select the inspection subset, `--all` scans every candidate contig, `--mode` controls planner behavior for both scopes, and breakpoint limits apply per contig. |
| `0.2.1` | Tightened `chromo sort` duplicate filtering for contaminated/alternate-fragment assemblies by using span-based overlap by default, requiring both novel coverage thresholds, rescuing very large near-threshold alignments, and letting split candidates protect their secondary reference spans. |
| `0.2.0` | Added minimap2 PAF input for `chromo sort` and `chromo fix`, plus `chromo plot` PDF/SVG/PNG dot plots for coords/PAF with optional assignment-report query ordering. |
| `0.1.2` | Raised the default auto-split query-span support threshold to 5% so small terminal off-target blocks are reported for review instead of being cut automatically. |
| `0.1.1` | Tightened `chromo fix` breakpoint placement by collapsing adjacent same-reference/orientation runs, added complex same-reference orientation detection, added a run-level auto breakpoint budget, protected strong multi-reference split candidates during `chromo sort`, and documented the fix-before-sort workflow for suspected misjoins. |
| `0.1.0` | Initial public package with `chromo sort`, `chromo fix`, `chromo scaffold`, duplicate-overlap filtering, user-nominated contig splitting, conservative auto smoothing, inferred/fixed-gap scaffolding, and synthetic tests. |
