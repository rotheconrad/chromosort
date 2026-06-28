---
title: chromo clean
description: Usage, outputs, parameters, and reasoning for chromo clean.
---

# chromo clean

Use `chromo clean` when an assembly is already mostly correct but needs
reference-guided cleanup: removing unaligned or redundant fragments, splitting a
small number of strong misjoins, optionally orienting records to the reference,
and writing one cleaned FASTA in reference order.

`chromo clean` makes all decisions from the raw assembly and the raw alignment
file you provide. It does not run MUMmer or minimap2. After it writes
`<prefix>.clean.fa`, re-run MUMmer or minimap2 on that cleaned FASTA for final
validation plots.

## What `chromo clean` Does

Given a reference FASTA, raw assembly FASTA, and MUMmer coords or minimap2 PAF,
`chromo clean`:

1. Runs the same assignment and duplicate-overlap filtering logic as
   `chromo sort` on the raw assembly.
2. Marks raw contigs as retained or discarded.
3. Selects retained raw contigs for fix planning. The default is all retained
   contigs.
4. Runs the same conservative breakpoint planner as `chromo fix` on selected
   retained raw contigs.
5. Replaces accepted split contigs with reference-labeled pieces.
6. Keeps retained unsplit contigs unchanged.
7. Optionally reverse-complements emitted records with `--orient-to-reference`.
8. Orders emitted records by reference FASTA order and reference placement.
9. Writes a cleaned FASTA plus initial sort, fix, clean, and run-summary
   reports.

The command is not implemented as "write `ordered.fa`, then run `chromo fix` on
that file." Raw coords or PAF are never fed into a sub-step expecting renamed or
oriented `ordered.fa` records.

## Run `chromo clean`

```bash
chromo clean \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.fa \
  --coords mummer/raw.coords \
  --output-prefix results/sample \
  --orient-to-reference \
  --discarded-fasta results/sample.discarded.fa
```

The same workflow can use minimap2 PAF:

```bash
chromo clean \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.fa \
  --paf paf/raw.paf \
  --output-prefix results/sample \
  --orient-to-reference
```

If you only want to run fix planning on split candidates surfaced by the sort
step:

```bash
chromo clean \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.fa \
  --coords mummer/raw.coords \
  --output-prefix results/sample \
  --fix-scope split-candidates
```

## `chromo clean` Outputs

| Output | Description |
| --- | --- |
| `<prefix>.clean.fa` | Retained raw contigs and accepted split pieces, oriented if requested and ordered by reference placement. |
| `<prefix>.clean.agp` | AGP 2.1 map for the cleaned FASTA, mapping each emitted record back to its raw source-contig slice and orientation. |
| `<prefix>.clean_components.tsv` | Human-readable provenance table mirroring the clean AGP rows, with clean status for each component. |
| `<prefix>.initial_sort.contig_assignments.tsv` | Full first-pass sort assignment report for raw contigs. |
| `<prefix>.initial_sort.contig_ref_matches.tsv` | Per-contig/reference match report from the raw alignment. |
| `<prefix>.initial_sort.chromosome_summary.tsv` | Chromosome summary from the initial sort decision. |
| `<prefix>.fix_targets.txt` | Original raw contig IDs selected for fix planning. |
| `<prefix>.fix_report.tsv` | `chromo fix`-style report for selected retained contigs. |
| `<prefix>.clean_contigs.tsv` | Main unified audit table for discarded contigs, retained unsplit contigs, and retained split pieces. |
| `<prefix>.clean_chromosome_summary.tsv` | Final cleaned-record summary grouped by reference sequence. |
| `<prefix>.submission_checklist.tsv` | FASTA/AGP consistency checks and file manifest for downstream submission packaging. |
| `<prefix>.run_summary.txt` | Inputs, outputs, sort/fix/clean status counts, and validation guidance. |
| `--discarded-fasta` path | Optional FASTA of raw contigs discarded by the initial sort step. |

### Example `clean_contigs.tsv` Rows

| source_contig | clean_status | clean_name | kept_by_sort | sort_status | fix_selected | fix_status | dominant_ref | slice_start | slice_end |
| --- | --- | --- | --- | --- | --- | --- | --- | ---: | ---: |
| `contig_01` | `discarded_no_alignment` | `.` | `no` | `no_alignment` | `no` | `.` | `.` | `.` | `.` |
| `contig_04` | `kept_split_piece` | `chrom02_contig_04_a` | `yes` | `kept_split_candidate` | `yes` | `split` | `chrom02` | `1` | `20` |
| `contig_04` | `kept_split_piece` | `chrom07_contig_04_b` | `yes` | `kept_split_candidate` | `yes` | `split` | `chrom07` | `21` | `40` |
| `contig_inv_mid` | `not_split_single_target` | `chrom06_contig_inv_mid` | `yes` | `kept` | `yes` | `not_split_single_target` | `chrom06` | `1` | `55` |

## `chromo clean` Parameters

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `--coords` | required unless `--paf` | MUMmer `show-coords` file generated from the raw assembly FASTA. |
| `--paf` | required unless `--coords` | minimap2 PAF generated from the raw assembly FASTA. |
| `--output-prefix` | required | Prefix for clean FASTA and audit reports. |
| `--discarded-fasta` | none | Optional FASTA of raw contigs discarded by sort filtering. |
| `--orient-to-reference` | off | Reverse-complement emitted records whose dominant alignment is reverse-strand. |
| `--fix-scope` | `kept` | Which raw contigs to inspect with fix planning: `kept`, `split-candidates`, `kept-and-split-candidates`, or `file`. |
| `--fix-targets-file` | none | Original raw contig IDs to inspect when `--fix-scope file` is used. |
| `--fix-mode` | `conservative` | Breakpoint planner mode passed to the internal fix step. |
| `--min-segment-bp` | `10000` | Minimum query-aligned bp for an alignment segment to inform splitting. |
| `--breakpoint-penalty-bp` | `50000` | Smoothed breakpoint penalty for conservative/chromosome/comprehensive modes. |
| `--min-piece-aligned-bp` | `50000` | Minimum dominant aligned bp required in each smoothed split piece. |
| `--min-piece-query-frac` | `0.05` | Minimum query-span fraction required for each smoothed split piece. |
| `--max-breakpoints-per-contig` | `4` | Maximum accepted breakpoints per source contig. |

All standard `chromo sort` assignment and duplicate-overlap thresholds are also
available, including `--min-aligned-bp`, `--min-query-cov`,
`--min-best-ref-share`, `--min-novel-ref-bp`, `--min-novel-ref-frac`,
`--overlap-mode`, `--novel-ref-criteria`, terminal-overlap rescue thresholds,
and split-candidate thresholds.

## Reasoning Behind `chromo clean`

### Discard Before Fix

The cleanup target is a mostly good assembly. In that situation, obvious
unaligned fragments, redundant overlaps, and alternate fragments should not
receive automatic breakpoint surgery. `chromo clean` first asks which raw
contigs are worth retaining, then runs fix planning only on the selected
retained raw contigs by default.

### Fix Before Final Order

Accepted splits create new pieces whose reference placement may differ from the
source contig's best reference. `chromo clean` therefore plans splits before
the final cleaned FASTA is ordered. Split pieces and unsplit contigs are then
sorted together by final dominant reference placement.

### Still Re-Align The Cleaned FASTA

The cleaned FASTA is derived from raw alignment evidence. That is enough for an
auditable cleanup decision, but it is not a substitute for final validation.
After `chromo clean`, align `<prefix>.clean.fa` to the reference and inspect
fresh `chromo plot` or `mummerplot` output from that clean-FASTA alignment. The
[dot-plot guide]({{ '/dot-plots/' | relative_url }}) gives examples of the
patterns to look for during that validation step.
