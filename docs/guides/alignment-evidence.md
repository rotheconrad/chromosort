---
title: Alignment Evidence And The Exact FASTA Rule
description: Learn when ChromoSort coords and PAF evidence can be reused and when edited FASTA outputs need fresh alignments.
---

# Alignment Evidence And The Exact FASTA Rule

Use this guide when you are asking:

> Can I reuse this coords or PAF file, or do I need to re-align?

The short answer is: a MUMmer coords file or minimap2 PAF file describes one
exact reference FASTA and one exact assembly FASTA. You can use that evidence
for more than one decision about the same assembly. You should not treat it as
evidence for a FASTA that has been renamed, filtered, split, cut, manually
exported, scaffolded, reverse-complemented, or gapfilled.

## The Core Idea

Alignment evidence is not a floating set of genome annotations. It is a set of
rows in two coordinate systems:

- reference coordinates from the reference FASTA used by the aligner,
- query coordinates from the assembly FASTA used by the aligner.

If either coordinate system changes, the rows no longer describe the new file.
Sometimes the biological sequence is still mostly the same, but the names,
membership, order, orientation, or coordinate offsets have changed enough that
alignment-dependent commands need new evidence.

This rule protects three things:

- **Coordinates:** a breakpoint at base `400000` in `raw_contig_7` is not
  automatically base `400000` in a split, trimmed, or scaffolded record.
- **Names:** PAF query names and coords query names must exist in the assembly
  FASTA passed to the command.
- **Interpretation:** a plot drawn from old rows can be useful review context,
  but it is not validation of an edited FASTA.

## What ChromoSort Reads Or Writes

| Stage | Can reuse the old alignment? | Why |
| --- | --- | --- |
| `raw.fa` with `raw.coords` or `raw.paf` | Yes, for more decisions about `raw.fa`. | The query FASTA still matches the alignment rows. |
| `chromo sort` reports from `raw.fa` | Yes, for reviewing `raw.fa` decisions. | Assignment reports summarize the same raw alignment evidence. |
| `<prefix>.ordered.fa` | No, if it becomes the assembly input to another alignment-dependent command. | Records may be renamed, filtered, reordered, and optionally reverse-complemented. |
| `<prefix>.clean.fa` | No. | `chromo clean` can discard records, emit split pieces, orient records, and rename outputs. |
| `fixed.fa` from `chromo fix` | No. | Split pieces have new names and new query coordinate systems. |
| Cut FASTA from `chromo cut` | No. | Exact cuts create new records and slice coordinates. |
| Browser or `manual apply` FASTA | No. | Manual edits can remove, invert, split, reorder, and scaffold records. |
| `<prefix>.scaffold.fa` | No, if later evidence should describe scaffold records. | Scaffold records are joined composites, not the original contigs. |
| `<prefix>.gapfilled.fa` | No. | Accepted graph fills and right-flank trims change scaffold sequence. |

The main exception is `chromo plot --assignments`. It still plots the original
alignment rows, but orders the query axis using a `chromo sort` assignment
report. That is a review view of the original alignment, not a fresh alignment
of `<prefix>.ordered.fa`.

## Safe Workflow Patterns

### Review And Fix The Same Raw Assembly

This is safe because both commands use the same `raw.fa` and the same raw
alignment evidence:

```text
raw.fa + raw.coords
  -> chromo sort for assignment and split-candidate review
  -> chromo fix on raw.fa using raw.coords
```

After `chromo fix` writes `fixed.fa`, the safe path changes:

```text
fixed.fa
  -> re-align with MUMmer or minimap2
  -> chromo sort, chromo plot, or chromo manual using fixed.coords/fixed.paf
```

### Clean A Mostly Correct Assembly

`chromo clean` uses raw alignment evidence to make cleanup decisions. Its output
needs fresh validation evidence:

```text
raw.fa + raw.paf
  -> chromo clean
  -> clean.fa
  -> re-align clean.fa
  -> validation plots from clean.paf
```

### Review A Sort Without Re-aligning

This is useful but easy to overread:

```bash
chromo plot \
  --ref-fasta reference.fa \
  --assembly-fasta raw.fa \
  --paf raw.paf \
  --assignments results/sample.contig_assignments.tsv \
  --output-prefix plots/sample.sorted_review \
  --per-ref
```

This plot reviews how the raw alignment rows support the sort decisions. It
does not prove that `results/sample.ordered.fa` aligns cleanly. To validate the
ordered FASTA itself, align `results/sample.ordered.fa` and plot the new PAF or
coords.

## Practical Review Workflow

1. Name the exact reference FASTA and assembly FASTA that produced the alignment.
2. Ask whether the command you want to run will read that same assembly FASTA.
3. If the command will read a changed FASTA, re-run MUMmer or minimap2 first.
4. Keep output folders grouped by assembly stage, such as `raw/`, `fixed/`, and
   `clean/`.
5. Use plots from the same FASTA stage as the decision you are making.
6. Treat old plots as historical review context, not final validation.

## Cheat Sheet

| Question | Answer |
| --- | --- |
| Can I run `chromo fix` on `raw.fa` after inspecting `chromo sort` reports from `raw.fa`? | Yes, if the same raw coords or PAF is used. |
| Can I run `chromo fix` on `sample.ordered.fa` using `raw.coords`? | No. Re-align `sample.ordered.fa` first. |
| Can I draw `chromo plot --assignments` with raw alignment rows? | Yes, as a sort-review plot of the raw assembly. |
| Does `chromo plot --assignments` validate `ordered.fa`? | No. It does not re-align the ordered FASTA. |
| Can I scaffold immediately after sorting? | Yes, with the matching `ordered.fa` and assignment TSV from the same sort run. |
| Can I use old contig-level PAF to validate scaffold records? | No. Align the scaffold FASTA if scaffold-level validation is needed. |

## Common Traps

Do not mix `raw.coords` with `fixed.fa`. The split pieces in `fixed.fa` did not
exist when `raw.coords` was made.

Do not assume unchanged-looking sequence means unchanged coordinates. Renaming,
removing records, and reverse-complementing retained records are enough to break
alignment-dependent assumptions.

Do not use `chromo plot --assignments` as final validation of an edited FASTA.
It is a review convenience that redraws old alignment rows.

Do not let FASTA indexes drift. If a `.fai` file was generated for an earlier
FASTA version, regenerate it or remove it before running tools that rely on
length lookup.

Do not keep every file in one flat directory with the same sample prefix. Stage
folders such as `raw/`, `fixed/`, and `clean/` make mistakes easier to see.

## What To Look At Next In ChromoSort

- Use [Input Files]({{ '/input-files/' | relative_url }}#fasta-and-alignment-compatibility)
  for the formal file compatibility rules.
- Use [Choosing PAF Or MUMmer Coords]({{ '/guides/paf-vs-coords/' | relative_url }})
  when deciding how to create the next alignment.
- Use [FASTA And Evidence Name Matching]({{ '/guides/name-matching/' | relative_url }})
  before starting a long run.
- Use [How to Interpret Dot Plots]({{ '/dot-plots/' | relative_url }}) after
  re-aligning an edited FASTA.
