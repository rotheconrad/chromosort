---
title: Proposal: chromo clean
description: Draft design for a conservative one-command cleanup workflow.
---

# Proposal: `chromo clean`

Status: design proposal, not implemented.

`chromo clean` would provide a conservative one-command cleanup path for
assemblies that are mostly correct, but contain short unaligned or redundant
fragments plus one or a few reference-scale misjoins or complex inversions.

The command should not hide the FASTA/alignment compatibility rule. It can make
multiple decisions from one raw alignment because all decisions are about the
same raw assembly before the final FASTA is written. The cleaned FASTA should
still be re-aligned for final validation.

## Target Use Case

Use `chromo clean` when:

- The assembly is already close to chromosome scale.
- The reference is high quality and close enough for confident whole-genome
  alignment.
- Most contigs are expected to place cleanly.
- The main cleanup needs are removing unaligned or redundant fragments,
  orienting/order placed contigs, and conservatively splitting rare
  chimeric contigs or large complex inversions.

Do not use `chromo clean` as the first choice when:

- The raw dot plot shows many widespread rearrangements.
- Many contigs appear to need multiple breaks.
- The reference may be structurally different enough that automatic splitting
  would erase real biology.
- Manual curation or graph evidence is required before sequence edits.

## Proposed CLI

```bash
chromo clean \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.fa \
  --coords mummer/raw.coords \
  --output-prefix results/sample \
  --orient-to-reference
```

Equivalent PAF input:

```bash
chromo clean \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.fa \
  --paf paf/raw.paf \
  --output-prefix results/sample \
  --orient-to-reference
```

Important options:

```text
--fix-mode conservative|chromosome|comprehensive|sensitive
--fix-scope kept|split-candidates|kept-and-split-candidates|file
--fix-targets-file FILE
--discarded-fasta FILE
--reports-only
--projected-plot
--projected-alignment-prefix PREFIX
```

Defaults:

```text
--fix-mode conservative
--fix-scope kept
--orient-to-reference off, matching chromo sort and chromo fix conventions
```

`--fix-scope kept` means the command first uses sort logic to decide which raw
contigs survive placement and duplicate-overlap filtering, then runs the
conservative fix planner only on those retained original contigs. This matches
the cleanup goal: discard obvious spurious records before allowing automatic
sequence surgery.

## Proposed Algorithm

1. Parse the reference FASTA, raw assembly FASTA, and raw coords or PAF.
2. Run the same assignment and overlap-filter logic as `chromo sort` on the raw
   assembly, but do not write `ordered.fa` yet.
3. Mark raw contigs as retained or discarded.
4. Select fix targets from retained raw contigs according to `--fix-scope`.
5. Run the same planner as `chromo fix --mode conservative` on selected original
   contigs, using the same raw alignment rows.
6. Keep unselected retained contigs unchanged.
7. Replace split retained contigs with accepted split pieces.
8. Orient each emitted contig or piece to the reference if
   `--orient-to-reference` is set.
9. Order emitted records by reference FASTA order and reference placement.
10. Write the cleaned FASTA and audit reports.

This should be implemented as one coordinated command rather than literally
writing `ordered.fa` and then running `chromo fix` on it. The internal sequence
of decisions is:

```text
raw alignment -> sort decisions on raw contigs
              -> fix decisions on retained raw contigs
              -> final orientation/order of emitted records
              -> clean.fa
```

The command must not feed raw coords or PAF into a sub-step that expects the
already-renamed or already-oriented `ordered.fa` records.

## Why Discard Before Fix

Discarding first keeps the automatic fixer focused on contigs that already have
enough reference support to belong in the cleaned assembly. It avoids spending
breakpoint logic on contaminant fragments, alternate haplotigs, or short
unplaced records that the final FASTA would not retain anyway.

The retained-contig restriction should be visible in the report. A discarded
contig should never silently disappear without a sort status explaining why.

## Why Fix Before Final Re-Order

Accepted splits create new pieces whose natural order should be based on their
dominant reference placement, not on the original source-contig order. Therefore
the final order should be assigned after fix planning. This also lets large
reference-supported inversions become separate oriented pieces when the selected
fix mode allows them.

## Proposed Outputs

```text
<prefix>.clean.fa
<prefix>.discarded.fa                         optional
<prefix>.initial_sort.contig_assignments.tsv
<prefix>.initial_sort.contig_ref_matches.tsv
<prefix>.initial_sort.chromosome_summary.tsv
<prefix>.fix_targets.txt
<prefix>.fix_report.tsv
<prefix>.clean_contigs.tsv
<prefix>.clean_chromosome_summary.tsv
<prefix>.run_summary.txt
```

Optional projected visual-review outputs:

```text
<prefix>.projected.paf
<prefix>.projected.coords
<prefix>.projected_plot.pdf
<prefix>.projected_plot.<ref>.pdf
```

Projected alignments would be derived from the raw alignment rows after applying
discard, split, orientation, and naming transforms. They would be useful for
quick review of the clean command's own decisions, but must be labeled clearly
as projected evidence rather than a fresh MUMmer or minimap2 alignment.

## Report Semantics

`<prefix>.clean_contigs.tsv` should be the main user-facing audit table. Useful
columns:

```text
source_contig
source_length
clean_status
clean_name
kept_by_sort
sort_status
sort_assigned_ref
sort_query_cov
sort_best_ref_share
sort_overlap_class
fix_selected
fix_status
part_index
slice_start
slice_end
piece_bp
dominant_ref
ref_start
ref_end
orientation
reverse_complemented
order_in_ref
reason
```

Suggested `clean_status` values:

```text
kept_unsplit
kept_split_piece
discarded_no_alignment
discarded_below_threshold
discarded_ambiguous
discarded_duplicate_overlap
discarded_terminal_overlap
not_split_smooth
not_split_too_many_breakpoints
```

## Safety Defaults

`chromo clean` should be conservative by default:

- Use current `chromo sort` defaults for assignment and duplicate-overlap
  filtering.
- Use `--fix-mode conservative`.
- Use current `chromo fix` defaults for breakpoint smoothing and
  `--max-breakpoints-per-contig 4`.
- Run fix only on retained contigs by default.
- Preserve complete initial sort and fix reports.
- Warn when any retained contig has more than the breakpoint cap, many
  secondary references, or a complex same-reference orientation pattern that is
  not cut in conservative mode.

## Validation Guidance

The command should end its run summary with explicit next steps:

```text
The cleaned FASTA was generated from raw alignment evidence.
For final validation, re-run MUMmer or minimap2 against:
  <prefix>.clean.fa
Then inspect chromo plot or mummerplot from the clean-FASTA alignment.
```

If `--projected-plot` is used, the summary should add:

```text
Projected plots show transformed raw alignment evidence.
They are useful for audit, but they are not a substitute for aligning clean.fa.
```

## Implementation Notes

The clean command should reuse internal sort and fix functions rather than
shelling out to `chromo sort` and `chromo fix`. The key new code is a reconciler
that combines:

- retained/discarded decisions from sort assignments,
- split-piece decisions from fix planning,
- final reference order and orientation for emitted records,
- unified clean reports.

Projected alignment output can be a second phase. The first useful
implementation can write `clean.fa` plus reports and require users to re-align
for plots.
