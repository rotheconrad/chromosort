---
title: Mostly-Correct Assembly Cleanup
description: A worked chromo clean path from raw FASTA and raw alignments to cleaned FASTA validation plots.
---

# Mostly-Correct Assembly Cleanup

Use this walkthrough when an assembly is already broadly correct, but still has
small unaligned fragments, redundant contigs, or a few strong split candidates
that should be handled in one conservative pass.

The goal is:

> Use raw evidence to make an auditable cleanup FASTA, then re-align that
> cleaned FASTA before trusting it downstream.

## When This Path Fits

`chromo clean` is a good first pass when:

- whole-genome dot plots are mostly syntenic,
- most contigs have one dominant reference placement,
- questionable sequence is limited to small fragments or a few obvious
  candidates,
- you want discarded raw contigs, split pieces, retained contigs, and final
  records in one report family.

Use a more explicit reviewed path instead when many contigs are chimeric, when
same-reference inversions need biological review, or when exact manual cuts are
already known.

## Inputs

Start with one exact FASTA pair and one primary alignment:

```text
reference.fa
assembly.raw.fa
paf/raw.ref_vs_asm.paf
```

PAF is shown below, but the same workflow can use MUMmer coords.

## Step 1: Make A Raw Review Plot

```bash
chromo plot \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.raw.fa \
  --paf paf/raw.ref_vs_asm.paf \
  --output-prefix plots/raw \
  --per-ref
```

The raw plot helps confirm that this is really a cleanup case. Look for mostly
clean diagonals, a small number of split-like contigs, and ordinary fragments
rather than widespread rearrangement.

## Step 2: Run Conservative Cleanup

```bash
chromo clean \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.raw.fa \
  --paf paf/raw.ref_vs_asm.paf \
  --output-prefix results/sample.clean \
  --orient-to-reference \
  --discarded-fasta results/sample.clean.discarded.fa
```

By default, `clean`:

1. runs sort assignment and duplicate-overlap filtering on raw contigs,
2. selects retained raw contigs for conservative fix planning,
3. emits accepted split pieces and retained unsplit contigs,
4. optionally orients emitted records to the reference,
5. orders the cleaned records by reference placement.

The command does not re-align. It uses the raw assembly and raw alignment to
write an edited FASTA.

## Step 3: Read The Cleanup Reports

Open these files before using `clean.fa`:

| File | Question |
| --- | --- |
| `sample.clean.initial_sort.contig_assignments.tsv` | Which raw contigs were kept, discarded, or flagged as split candidates? |
| `sample.clean.fix_targets.txt` | Which original raw contigs went through fix planning? |
| `sample.clean.fix_report.tsv` | Which selected contigs split, copied unchanged, or were smoothed over? |
| `sample.clean.clean_contigs.tsv` | Which final records were emitted, discarded, or retained unsplit? |
| `sample.clean.run_summary.txt` | Which inputs, thresholds, and status counts define the run? |
| `sample.clean.discarded.fa` | Which raw sequences were excluded from the cleaned FASTA? |

Review these status patterns:

| Pattern | What to check |
| --- | --- |
| Many `discarded_no_alignment` rows | Name matching, aligner sensitivity, or true unplaced sequence. |
| Many `duplicate_overlap` rows | Possible haplotigs, repeats, or redundant fragments. |
| `kept_split_piece` rows | Whether the split was expected and each piece has enough support. |
| `not_split_smooth` rows | Whether the planner smoothed over weak discordance appropriately. |
| `not_split_single_target` rows | Whether the contig stayed as one reference-consistent record. |

## Step 4: Narrow Fix Scope When Needed

If the default run inspects too many retained contigs, rerun with split
candidates only:

```bash
chromo clean \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.raw.fa \
  --paf paf/raw.ref_vs_asm.paf \
  --output-prefix results/sample.clean_candidates \
  --fix-scope split-candidates \
  --orient-to-reference
```

Use this when the assembly is very clean and you want breakpoint planning only
where the initial sort report already raised a split-candidate signal.

## Step 5: Re-Align The Cleaned FASTA

The cleaned FASTA has new names, membership, orientation, and sometimes split
piece coordinates. Make fresh alignment evidence:

```bash
minimap2 -x asm5 -c -t 16 --secondary=no \
  reference.fa results/sample.clean.clean.fa \
  > paf/sample.clean.ref_vs_asm.paf
```

Then validate the cleaned assembly:

```bash
chromo plot \
  --ref-fasta reference.fa \
  --assembly-fasta results/sample.clean.clean.fa \
  --paf paf/sample.clean.ref_vs_asm.paf \
  --output-prefix plots/sample.clean \
  --per-ref
```

The validation plot is the evidence for the cleaned FASTA. The raw plot is
historical context.

## Step 6: Continue From The Cleaned Stage

If validation looks good, downstream commands should use the cleaned FASTA and
fresh cleaned alignment:

```bash
chromo sort \
  --ref-fasta reference.fa \
  --assembly-fasta results/sample.clean.clean.fa \
  --paf paf/sample.clean.ref_vs_asm.paf \
  --output-prefix results/sample.final \
  --orient-to-reference
```

Then scaffold only when the final `ordered.fa` and `contig_assignments.tsv`
come from the same sort run.

## Common Traps

Do not treat `chromo clean` as validation. It writes a cleaned FASTA from raw
evidence; the cleaned FASTA needs fresh evidence.

Do not mix `raw.paf` with `sample.clean.clean.fa`. Split and oriented records
have different query coordinates.

Do not assume every discarded contig is biologically irrelevant. Read the
discarded FASTA and assignment report when haplotypes, novel sequence, or
repeats matter.

Do not use this path for complex manual curation. Switch to eval/manual review
when many decisions need human judgment.

## What To Look At Next In ChromoSort

- Use [chromo clean]({{ '/commands/clean/' | relative_url }}) for exact
  parameters and output names.
- Use [Alignment Evidence And The Exact FASTA Rule]({{ '/guides/alignment-evidence/' | relative_url }})
  before chaining outputs.
- Use [Reading ChromoSort Audit Tables]({{ '/guides/audit-tables/' | relative_url }})
  while inspecting cleanup reports.
- Use [Suspected Chimeric Contig Review]({{ '/guides/chimeric-contig-walkthrough/' | relative_url }})
  when `clean` surfaces a contig that deserves focused review.
