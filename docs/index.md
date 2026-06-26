---
title: Documentation
description: User guide for ChromoSort.
---

# ChromoSort Documentation

ChromoSort is a reference-guided assembly curation toolkit for sorting contigs, conservatively cleaning mostly-correct assemblies, splitting reviewed chimeric contigs, making manual cuts, reviewing dot plots, scaffolding ordered contigs, applying reviewed graph-supported gap fills, and plotting existing whole-genome alignments.

This documentation is organized as a user guide for installing ChromoSort, preparing inputs, running common workflows, and reading command outputs.

The most important workflow rule is that alignment evidence matches one exact
FASTA pair. MUMmer coords or minimap2 PAF generated from `raw.fa` can support
multiple decisions about `raw.fa`, but edited outputs such as `ordered.fa`,
`fixed.fa`, manual FASTA exports, and scaffold FASTAs need fresh alignments
before they are used as inputs to another alignment-dependent command.

## Start Here

- [Installation](https://rotheconrad.github.io/chromosort/installation/)
- [Input files](https://rotheconrad.github.io/chromosort/input-files/)
- [Workflows](https://rotheconrad.github.io/chromosort/workflows/)
- [Agent and review playbook]({{ '/review-playbook/' | relative_url }})
- [Guides]({{ '/guides/' | relative_url }}) for the FASTA/alignment rule,
  PAF-vs-coords choices, name matching, audit tables, sorting decisions,
  command choice, breakpoint review, inversions, and dot plots.
- [How to interpret dot plots]({{ '/dot-plots/' | relative_url }})
- [Command reference](https://rotheconrad.github.io/chromosort/commands/)
- [Production upgrade roadmap]({{ '/roadmap/' | relative_url }})

## Commands

- [chromo sort](https://rotheconrad.github.io/chromosort/commands/sort/)
- [chromo clean](https://rotheconrad.github.io/chromosort/commands/clean/)
- [chromo eval](https://rotheconrad.github.io/chromosort/commands/eval/)
- [chromo fix](https://rotheconrad.github.io/chromosort/commands/fix/)
- [chromo cut](https://rotheconrad.github.io/chromosort/commands/cut/)
- [chromo manual](https://rotheconrad.github.io/chromosort/commands/manual/)
- [chromo plot](https://rotheconrad.github.io/chromosort/commands/plot/) and
  [dot-plot interpretation]({{ '/dot-plots/' | relative_url }})
- [chromo graph-map](https://rotheconrad.github.io/chromosort/commands/graph-map/)
- [chromo scaffold](https://rotheconrad.github.io/chromosort/commands/scaffold/)
- [chromo gapfill](https://rotheconrad.github.io/chromosort/commands/gapfill/)

## Running Analyses

- [Output files](https://rotheconrad.github.io/chromosort/outputs/)
- [Guides]({{ '/guides/' | relative_url }}) for educational explanations of
  evidence files, audit reports, command decisions, and visual review.
- [Agent and review playbook]({{ '/review-playbook/' | relative_url }}) for
  choosing one primary coords or PAF alignment, inversion review,
  long-read/GFA/GAF evidence, and handoffs between chats or projects.
- [Architecture](https://rotheconrad.github.io/chromosort/architecture/) for
  algorithm, data-model, evidence, and command-activation details.
- [Production upgrade roadmap]({{ '/roadmap/' | relative_url }})
- [Troubleshooting](https://rotheconrad.github.io/chromosort/troubleshooting/)

## Project Information

- [Current status and version history](https://rotheconrad.github.io/chromosort/status/)
- [Contributing](https://rotheconrad.github.io/chromosort/contributing/)

## Main Workflow

1. Generate MUMmer `show-coords` or minimap2 PAF alignments between the reference and assembly.
2. For mostly-correct assemblies, use `chromo clean` for conservative
   sort-filter-fix cleanup, then re-align the cleaned FASTA.
3. For more complex cases, use `chromo plot` and the
   [dot-plot guide]({{ '/dot-plots/' | relative_url }}), then use
   `chromo manual` when helpful to inspect suspicious contigs.
4. Use `chromo eval` and task-specific `chromo manual` dashboards when outlier
   `fix`, `scaffold`, or `gapfill` decisions need table or GUI review.
5. Use `chromo fix` or `chromo cut` for reviewed sequence edits.
6. Re-align the fixed assembly.
7. Use `chromo sort` to assign, filter, orient, and order contigs.
8. Re-align again if a later command should operate on `ordered.fa` rather than
   on the original assembly and its assignment report.
9. Use `chromo graph-map` or `chromo plot --gfa-overlay` when unitig-level GFA
   features need to be projected onto contig FASTA coordinates.
10. Use `chromo scaffold` to build one scaffold per reference sequence.
11. Use `chromo gapfill` only when reviewed GFA paths should replace N gaps.
