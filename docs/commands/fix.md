---
title: chromo fix
description: Usage, outputs, parameters, and reasoning for chromo fix.
---

# chromo fix

Use `chromo fix` when the goal is to split chimeric or structurally inconsistent
assembly contigs into reference-labeled pieces.

The coords or PAF file must describe the same assembly FASTA passed to
`--assembly-fasta`. It is valid to inspect `chromo sort` reports from `raw.fa`
and then run `chromo fix --assembly-fasta raw.fa --coords raw.coords` on selected
original contigs. It is not valid to run `chromo fix` on `sample.ordered.fa`
with coords or PAF that were generated from `raw.fa`. After `chromo fix` writes
`fixed.fa`, re-run MUMmer or minimap2 before sorting or plotting that fixed
FASTA. Use the [dot-plot guide]({{ '/dot-plots/' | relative_url }}) when
checking whether the repaired contigs now place cleanly.

## What `chromo fix` Does

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

## Run `chromo fix` With Selected Contigs

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

## Run `chromo fix` Across All Contigs

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

Optional graph context for the reviewed contigs:

```bash
chromo fix \
  --assembly-fasta assembly.fa \
  --paf minimap2/sample.paf \
  --contigs contig_04 contig_12 \
  --output-fasta results/sample.fixed.fa \
  --report results/sample.fixed_contigs.tsv \
  --gfa assembly_graph.gfa
```

## Run `chromo fix` From an Eval Table

Use `chromo eval fix` when you want a spreadsheet-editable table first:

```bash
chromo eval fix \
  --assembly-fasta assembly.fa \
  --paf minimap2/sample.paf \
  --contigs contig_04 contig_12 \
  --gfa assembly_graph.gfa \
  --read-paf reads_to_assembly.paf \
  --gaf reads_to_graph.gaf \
  --output-prefix results/sample.eval_fix
```

The optional `--gfa`, `--read-paf`, and `--gaf` inputs add review columns to
the table. They do not make `chromo fix` apply graph or long-read evidence
directly; the accepted table rows still supply explicit source slices.

After review, apply accepted `split_piece` rows:

```bash
chromo fix \
  --assembly-fasta assembly.fa \
  --reviewed-plan results/sample.eval_fix.fix_review.tsv \
  --output-fasta results/sample.fixed.fa \
  --report results/sample.fixed_contigs.tsv
```

With `--reviewed-plan`, the table supplies target contigs and exact slices, so
`--coords`/`--paf` and `--contigs`/`--all` are omitted.

## `chromo fix` Outputs

| Output | Description |
| --- | --- |
| `--output-fasta` | Full fixed assembly FASTA by default, with split pieces replacing fixed contigs. |
| `--report` | TSV report describing split pieces and unsplit requested contigs. |
| `--graph-report` | Optional graph context TSV when `--gfa` is provided. Defaults to the `--report` path with a `.graph.tsv` suffix. |

The report includes original contig name, split status, new contig name,
dominant reference, slice coordinates, alignment coordinates, orientation,
reverse-complement status, identity, segment count, and the reason for the
decision. Candidates that contain discordant blocks but are not cut by
breakpoint smoothing are reported as `not_split_smooth`. Candidates rejected by
the per-contig breakpoint cap are reported as `not_split_too_many_breakpoints`.

### Example `chromo fix` Output

**Table 1. Example split report rows.** Selected columns from a sensitive-mode
fixture show how one source contig can become multiple reference-labeled pieces.

| original_contig | status | new_contig | part_index | dominant_ref | slice_start | slice_end | piece_bp | orientation |
| --- | --- | --- | ---: | --- | ---: | ---: | ---: | --- |
| `contig_04` | `split` | `chrom02-contig_04-a` | `1` | `chrom02` | `1` | `20` | `20` | `+` |
| `contig_04` | `split` | `chrom07-contig_04-b` | `2` | `chrom07` | `21` | `40` | `20` | `+` |
| `contig_12` | `split` | `chrom04-contig_12-a` | `1` | `chrom04` | `1` | `5` | `5` | `+` |
| `contig_12` | `split` | `chrom05-contig_12-b` | `2` | `chrom05` | `6` | `35` | `30` | `+` |

**Listing 1. Example fixed FASTA records.** Split-piece FASTA headers carry the
original contig, slice interval, alignment interval, and orientation so the edit
can be audited after the sequence file leaves ChromoSort.

```text
>chrom02-contig_04-a original=contig_04 ref=chrom02 slice=1-20 alignment=1-20 orientation=+ reverse_complemented=no avg_identity=100.000
AAAAAAAAAAAAAAAAAAAA
>chrom07-contig_04-b original=contig_04 ref=chrom07 slice=21-40 alignment=21-40 orientation=+ reverse_complemented=no avg_identity=100.000
CCCCCCCCCCCCCCCCCCCC
```

## `chromo fix` Naming

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

## `chromo fix` Parameters

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `--coords` | required unless `--paf` | MUMmer `show-coords` alignment file. |
| `--paf` | required unless `--coords` | minimap2 PAF alignment file. |
| `--reviewed-plan` | none | Reviewed `chromo eval fix` table. Accepted `split_piece` rows are applied directly and replace the alignment-driven planner path; graph and long-read columns remain provenance for those reviewed slices. |
| `--gfa` | none | Optional assembly graph GFA for report-only context about reviewed source contigs. |
| `--graph-report` | report path with `.graph.tsv` suffix | Optional path for the `--gfa` graph context report. |
| `--graph-guard` | off | Requires `--gfa`; emits conservative warnings for graph-simple planned splits and graph-complex unsplit contigs without changing the fixed FASTA. |
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

## Reasoning Behind `chromo fix`

### Scope Is Separate From Mode

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

### Collapse Same-Target Runs Before Cutting

Whole-genome alignments often contain many neighboring rows for the same
reference and orientation, separated by local gaps, repeats, or small assembly
differences. `chromo fix` collapses adjacent same-target runs before placing
breakpoints. A contig is cut at accepted reference transitions, or at
same-reference orientation transitions only when they are complex enough or
explicitly enabled, not at every ordinary alignment-row boundary.

### Breakpoint-Penalty Segmentation

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

### Keep Graph Context Beside Split Decisions

When `--gfa` is provided, `chromo fix` writes a graph context report for the
requested contigs. This is useful after manual review because the split report
shows the alignment-supported edit, while the graph report shows whether the
source contig is present in the assembly graph and whether it sits in a simple
or tangled local graph neighborhood. The graph report does not alter breakpoint
planning or FASTA output.

With `--graph-guard`, `chromo fix` also writes stderr warnings when a planned
split sits in a simple graph neighborhood or an unsplit target sits in a
high-degree/self-loop graph neighborhood. The warning is meant to send those
cases back through manual review, not to rewrite breakpoints automatically.

### Breakpoints Between Alignment Blocks

After segmentation, `chromo fix` places accepted breakpoints halfway between
neighboring query-ordered alignment blocks. When blocks are adjacent, the
breakpoint lands at the alignment boundary. When there is an unaligned gap, the
gap is divided between the neighboring pieces instead of being discarded.

### Synthetic Test Cases

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
