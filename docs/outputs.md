---
title: Output Files
description: Guide to ChromoSort FASTA, TSV, HTML, plot, and run-summary outputs.
---

# Output Files

Most ChromoSort commands write both sequence files and audit tables. The FASTA outputs are intended for downstream assembly workflows; the TSV and text reports are intended to make every keep, reject, split, cut, scaffold, or graph-fill decision inspectable.

## Command Outputs

| Command | Primary outputs |
| --- | --- |
| [`chromo sort`]({{ '/commands/sort/' | relative_url }}) | `<prefix>.ordered.fa`, `<prefix>.contig_assignments.tsv`, `<prefix>.contig_ref_matches.tsv`, `<prefix>.chromosome_summary.tsv`, optional `<prefix>.graph_assignments.tsv`, and `<prefix>.run_summary.txt`. |
| [`chromo fix`]({{ '/commands/fix/' | relative_url }}) | Reviewed fixed FASTA at `--output-fasta`, split report at `--report`, and optional graph report. |
| [`chromo cut`]({{ '/commands/cut/' | relative_url }}) | Cut FASTA at `--output-fasta` and cut-piece report at `--report`. |
| [`chromo manual`]({{ '/commands/manual/' | relative_url }}) | Self-contained HTML dashboard, browser FASTA download, recipe JSON download, and reproducible `manual apply` FASTA/report outputs. |
| [`chromo plot`]({{ '/commands/plot/' | relative_url }}) | Whole-genome and optional per-reference dot plots in PDF, SVG, or PNG. |
| [`chromo scaffold`]({{ '/commands/scaffold/' | relative_url }}) | `<prefix>.scaffold.fa`, `<prefix>.scaffold_gaps.tsv`, optional `<prefix>.graph_gaps.tsv`, `<prefix>.scaffold_summary.tsv`, and `<prefix>.run_summary.txt`. |
| [`chromo gapfill`]({{ '/commands/gapfill/' | relative_url }}) | `<prefix>.gapfill_plan.tsv`, optional review HTML, optional `<prefix>.gapfilled.fa`, and `<prefix>.run_summary.txt`. |

## Audit Tables

The report tables are deliberately redundant. They keep original contig IDs, new contig IDs, reference assignments, coordinates, status labels, and decision metrics close to the sequence output they explain. Prefer using these tables instead of parsing FASTA headers when downstream scripts need assignment, split, or gap metadata.

## Sequence-Changing Outputs

Only a subset of commands change sequence:

- `chromo fix` replaces selected contigs with split pieces.
- `chromo cut` replaces selected contigs with pieces cut at exact reviewed positions.
- `chromo manual apply` reproduces a reviewed browser recipe.
- `chromo scaffold` joins sorted contigs with inferred or fixed N gaps and only trims overlaps when an explicit overlap policy asks it to.
- `chromo gapfill --apply` inserts graph sequence only for fillable and accepted paths.

`chromo sort`, `chromo plot`, and report-only graph evidence modes do not trim, polish, or rewrite sequence beyond optional orientation of retained contigs in `chromo sort --orient-to-reference`.
