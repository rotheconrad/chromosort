---
title: Guides
description: ChromoSort-specific guides, evidence guides, and worked review examples.
---

# Guides

The [Concepts]({{ '/concepts/' | relative_url }}) explain why a pattern might
justify action. These guides explain how to carry that reasoning into
ChromoSort inputs, commands, reports, reviewed tables, and reproducible
outputs.

If a command output looks wrong, start with
[Troubleshooting]({{ '/troubleshooting/' | relative_url }}). If you need the
visual grammar first, start with
[How to interpret dot plots]({{ '/dot-plots/' | relative_url }}).

## Evidence And Interpretation

| Guide | Use it when |
| --- | --- |
| [Alignment evidence and the exact FASTA rule]({{ '/guides/alignment-evidence/' | relative_url }}) | You need to know when coords/PAF can be reused and when edited FASTA outputs require fresh alignments. |
| [Choosing PAF or MUMmer coords]({{ '/guides/paf-vs-coords/' | relative_url }}) | You need to choose, compare, or sanity-check whole-genome alignment evidence for a ChromoSort stage. |
| [FASTA and evidence name matching]({{ '/guides/name-matching/' | relative_url }}) | Plots, graph overlays, or command reports look empty because FASTA, PAF, coords, GFA, GAF, or Hi-C names may not match. |
| [Reading ChromoSort audit tables]({{ '/guides/audit-tables/' | relative_url }}) | You need to translate TSV status labels into review decisions. |
| [Inversions and orientation changes]({{ '/guides/inversions-orientation/' | relative_url }}) | You need to connect orientation patterns to ChromoSort review choices. |
| [Assembly graph evidence]({{ '/guides/assembly-graph-evidence/' | relative_url }}) | You need to read GFA segments, links, paths, walks, graph complexity labels, and graph guardrails. |
| [hifiasm unitig-to-contig projection]({{ '/guides/hifiasm-graph-projection/' | relative_url }}) | You need to project unitig graph coordinates before they can annotate contig-axis plots. |
| [Long-read PAF and GAF support]({{ '/guides/long-read-paf-gaf/' | relative_url }}) | You need to interpret read bridges, breakpoint support, graph traversals, MAPQ filters, and advisory evidence. |

## Command Decisions And Review

| Guide | Use it when |
| --- | --- |
| [Sorting decisions and duplicate-overlap filtering]({{ '/guides/sorting-decisions/' | relative_url }}) | You need to understand kept, duplicate-overlap, terminal-overlap, split-candidate, and graph-guarded sort outcomes. |
| [Sort, clean, fix, cut, or manual?]({{ '/guides/choosing-commands/' | relative_url }}) | You need the least invasive ChromoSort action for the evidence in front of you. |
| [Chimeric contig and breakpoint review]({{ '/guides/breakpoint-review/' | relative_url }}) | You need to review planner modes, smoothing, breakpoint penalties, selected contigs, and manual split boundaries. |
| [Manual dashboard review]({{ '/guides/manual-dashboard-review/' | relative_url }}) | You need to move from browser review to a reproducible manual recipe. |
| [Spreadsheet review tables]({{ '/guides/spreadsheet-review-tables/' | relative_url }}) | You need to accept, reject, or annotate evaluated fix, scaffold, or gapfill rows safely. |
| [Scaffolding, gaps, and overlaps]({{ '/guides/scaffolding-gaps-overlaps/' | relative_url }}) | You need to review inferred gaps, fixed gaps, negative gaps, overlaps, trim policies, or gap overrides. |
| [Graph-supported gap filling]({{ '/guides/graph-gap-filling/' | relative_url }}) | You need to decide whether graph paths are strong enough to replace scaffold N gaps. |

## Worked Review Guides

| Guide | Use it when |
| --- | --- |
| [Mostly-correct assembly cleanup]({{ '/guides/mostly-correct-cleanup/' | relative_url }}) | You want a conservative `chromo clean` path from raw FASTA through validation plots. |
| [Suspected chimeric contig review]({{ '/guides/chimeric-contig-walkthrough/' | relative_url }}) | You want to go from a dot-plot signal to `eval fix`, manual review, reviewed splitting, re-alignment, and sorting. |
| [Graph-aware scaffold and gapfill review]({{ '/guides/graph-scaffold-gapfill-walkthrough/' | relative_url }}) | You want a combined GFA, GAF, Hi-C, reference-placement, scaffold, and gapfill review path. |
| [Dataset handoff checklist]({{ '/guides/dataset-handoff-checklist/' | relative_url }}) | You need to package FASTA, alignments, graph evidence, review tables, plots, summaries, and notes for a new review session. |
