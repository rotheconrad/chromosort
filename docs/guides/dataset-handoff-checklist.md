---
title: Dataset Handoff Checklist
description: Package ChromoSort inputs, evidence files, review tables, plots, reports, decisions, and notes for reproducible handoff.
---

# Dataset Handoff Checklist

Use this guide when a ChromoSort review is moving to another person, another
project folder, or another chat.

The goal is:

> A reviewer should know which FASTA each evidence file describes, what was
> changed, what was left native, and what still needs review.

## The Core Package

Create a small manifest before handing off files:

| Item | Include |
| --- | --- |
| Reference FASTA | Path, version, source, and `.fai` status. |
| Assembly FASTA stages | Raw, clean, fixed, manual, sorted, scaffolded, or gapfilled FASTA paths. |
| Primary alignment | Exact coords or PAF file for each FASTA stage being reviewed. |
| ChromoSort commands | Commands used to produce reports, plots, FASTA outputs, and dashboards. |
| Run summaries | Every `<prefix>.run_summary.txt`. |
| Main audit tables | Assignment, fix, clean, manual, scaffold, and gapfill reports. |
| Review tables | `fix_review.tsv`, `scaffold_review.tsv`, `gapfill_review.tsv`, and edited reviewed plans. |
| Visual outputs | Dot plots, manual dashboards, gapfill review HTML, and focused per-reference panels. |
| Graph/read evidence | GFA, GAF, read-to-assembly PAF, graph-node reference PAF, Hi-C pair tables. |
| Decision notes | What changed, what was rejected, what remains unresolved. |

## Suggested Folder Layout

```text
handoff/
  README.md
  manifest.tsv
  commands.sh
  raw/
    reference.fa
    assembly.raw.fa
    raw.ref_vs_asm.paf
    plots/
  clean/
    sample.clean.fa
    sample.clean_contigs.tsv
    sample.run_summary.txt
    sample.clean.ref_vs_asm.paf
  fixed/
    sample.fixed.fa
    sample.fixed.tsv
    sample.fixed.ref_vs_asm.paf
  sort/
    sample.ordered.fa
    sample.contig_assignments.tsv
    sample.contig_ref_matches.tsv
  review/
    sample.fix_review.tsv
    sample.manual_fix.html
    sample.gapfill_review.tsv
  graph/
    assembly_graph.gfa
    reads_to_graph.gaf
    graph_nodes_to_ref.paf
  scaffold/
    sample.scaffold.fa
    sample.scaffold_gaps.tsv
    sample.graph_gaps.tsv
    sample.submission_checklist.tsv
  gapfill/
    sample.gapfill_plan.tsv
    sample.gapfill.review.html
    sample.reviewed_gapfill.gapfilled.fa
    sample.reviewed_gapfill.submission_checklist.tsv
```

The exact folder names do not matter. Stage separation matters.

## Manifest Fields

A simple TSV is enough:

| Column | Example |
| --- | --- |
| `stage` | `raw`, `clean`, `fixed`, `sort`, `scaffold`, `gapfill` |
| `file` | `clean/sample.clean.fa` |
| `file_type` | `fasta`, `paf`, `gfa`, `review_tsv`, `plot`, `html`, `summary` |
| `describes_fasta` | `clean/sample.clean.fa` |
| `created_by` | `chromo clean`, `minimap2`, `chromo plot` |
| `created_from` | `raw/assembly.raw.fa + raw/raw.ref_vs_asm.paf` |
| `review_status` | `accepted`, `rejected`, `needs_review`, `context_only` |
| `notes` | Short human note. |

The `describes_fasta` column is the important one. It prevents raw alignments
from being mistaken for validation of cleaned, fixed, or scaffolded FASTA.

## Decision Table

Include a small decision table for contested events:

| Target | Evidence files | Decision | Output affected | Notes |
| --- | --- | --- | --- | --- |
| `contig_04` | raw plot, `fix_review.tsv`, reads PAF | Accepted split | `fixed/sample.fixed.fa` | Two strong reference blocks. |
| `contig_inv_mid` | dot plot, GFA, reads PAF | Leave native | No edit | Likely real inversion or unresolved reference difference. |
| `chr1:left|right` | graph plan, GAF, Hi-C | Accepted fill | `gapfill/sample.reviewed_gapfill.gapfilled.fa` | Path `left+,bridge_good+,right+`. |
| `contigDup` | assignment table | Excluded from ordered FASTA | `sort/sample.ordered.fa` | `duplicate_overlap`; keep source FASTA for audit. |

This table should be boring. Boring is good: future reviewers can see what
happened without reconstructing it from filenames.

## Commands File

Keep the runnable commands in order:

```bash
# 1. Primary raw alignment.
minimap2 -x asm5 -c -t 16 --secondary=no reference.fa assembly.raw.fa \
  > raw/raw.ref_vs_asm.paf

# 2. Raw review plot.
chromo plot \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.raw.fa \
  --paf raw/raw.ref_vs_asm.paf \
  --output-prefix raw/plots/sample.raw \
  --per-ref

# 3. Cleanup or reviewed fix.
chromo clean \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.raw.fa \
  --paf raw/raw.ref_vs_asm.paf \
  --output-prefix clean/sample \
  --orient-to-reference
```

Do not rely on shell history. Put the commands in a file that travels with the
data.

## Handoff Notes For Common Cases

### Cleaned Assembly

Include:

- raw FASTA and raw PAF or coords,
- `clean.fa`,
- `clean_contigs.tsv`,
- `initial_sort.*` reports,
- `fix_report.tsv`,
- cleaned FASTA re-alignment and validation plots.

State whether `clean.fa` has already been re-aligned.

### Reviewed Chimeric Contig

Include:

- raw plot around the candidate,
- `fix_review.tsv`,
- any manual dashboard,
- `fix --reviewed-plan` report,
- fixed FASTA,
- fixed FASTA re-alignment and validation plot.

State which candidates were rejected and why.

### Graph-Aware Gapfill

Include:

- ordered FASTA and matching assignment TSV,
- GFA used for planning,
- GAF, graph-node reference PAF, and Hi-C pairs when used,
- gapfill plan,
- review HTML or edited reviewed plan,
- applied gapfill report and FASTA.

State whether the GFA had sequence and whether any accepted rows were rejected
as stale.

## Preflight Checks For The Receiver

Before continuing a handoff, the next reviewer should ask:

1. Which FASTA does this PAF or coords file describe?
2. Do the FASTA IDs match the alignment query names?
3. Do GFA segment names match the assembly or assignment records being used?
4. Were any FASTA records renamed, split, oriented, scaffolded, or gapfilled?
5. Are reviewed tables from the same FASTA, assignment, graph, and path-search
   settings as the executor run?
6. Which rows changed sequence, and where is the report?

If any answer is unclear, regenerate evidence before changing sequence.

## Common Traps

Do not hand off only the final FASTA. ChromoSort is built around audit tables.

Do not bundle every stage in one flat folder with reused sample prefixes.

Do not omit rejected reviewed rows when those rows explain why sequence was
left native.

Do not assume a graph or read file is self-explanatory. Record which command
used it and whether it was report-only or sequence-changing.

Do not pass a browser manual FASTA without the recipe JSON or apply report when
the edit needs to be reproduced.

## What To Look At Next In ChromoSort

- Use [Alignment Evidence And The Exact FASTA Rule]({{ '/guides/alignment-evidence/' | relative_url }})
  to audit stage compatibility.
- Use [Reading ChromoSort Audit Tables]({{ '/guides/audit-tables/' | relative_url }})
  for output report interpretation.
- Use [Manual Dashboard Review]({{ '/guides/manual-dashboard-review/' | relative_url }})
  when a handoff includes browser recipes.
- Use the [Agent And Review Playbook]({{ '/review-playbook/' | relative_url }})
  for broader review practices and delivery labels.
