---
title: Guides
description: ChromoSort-specific guides for commands, reports, and reviewed plans.
---

# Guides

These guides explain how ChromoSort turns alignment, graph, read, and review
evidence into command decisions, reports, reviewed plans, and sequence-changing
outputs.

For general sequence-analysis ideas, start with
[Concepts]({{ '/concepts/' | relative_url }}). For full review paths, use
[Examples]({{ '/examples/' | relative_url }}). If a command output looks wrong,
start with [Troubleshooting]({{ '/troubleshooting/' | relative_url }}).

| Guide | Use it when |
| --- | --- |
| [Alignment evidence and the exact FASTA rule]({{ '/guides/alignment-evidence/' | relative_url }}) | You need to know when coords/PAF can be reused and when edited FASTA outputs require fresh alignments. |
| [FASTA and evidence name matching]({{ '/guides/name-matching/' | relative_url }}) | Plots, graph overlays, or command reports look empty because FASTA, PAF, coords, GFA, GAF, or Hi-C names may not match. |
| [Reading ChromoSort audit tables]({{ '/guides/audit-tables/' | relative_url }}) | You need to translate TSV status labels into review decisions. |
| [Sorting decisions and duplicate-overlap filtering]({{ '/guides/sorting-decisions/' | relative_url }}) | You need to understand kept, duplicate-overlap, terminal-overlap, split-candidate, and graph-guarded sort outcomes. |
| [Sort, clean, fix, cut, or manual?]({{ '/guides/choosing-commands/' | relative_url }}) | You need the least invasive ChromoSort action for the evidence in front of you. |
| [Chimeric contig and breakpoint review]({{ '/guides/breakpoint-review/' | relative_url }}) | You need to review planner modes, smoothing, breakpoint penalties, selected contigs, and manual split boundaries. |
| [Manual dashboard review]({{ '/guides/manual-dashboard-review/' | relative_url }}) | You need to move from browser review to a reproducible manual recipe. |
| [Spreadsheet review tables]({{ '/guides/spreadsheet-review-tables/' | relative_url }}) | You need to accept, reject, or annotate evaluated fix, scaffold, or gapfill rows safely. |
| [Scaffolding, gaps, and overlaps]({{ '/guides/scaffolding-gaps-overlaps/' | relative_url }}) | You need to review inferred gaps, fixed gaps, negative gaps, overlaps, trim policies, or gap overrides. |
| [Graph-supported gap filling]({{ '/guides/graph-gap-filling/' | relative_url }}) | You need to decide whether graph paths are strong enough to replace scaffold N gaps. |
