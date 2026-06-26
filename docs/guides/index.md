---
title: Guides
description: Concept and ChromoSort guides for reference-guided assembly review.
---

# Guides

The guides are split by how portable the lesson is.

**Concept guides** teach genome-assembly and sequence-analysis ideas that are
useful even outside ChromoSort. **ChromoSort guides** explain how those ideas
map onto ChromoSort commands, reports, reviewed plans, and sequence-changing
outputs. **Worked examples** show complete review paths.

For command failures, stale inputs, empty plots, graph projection problems, or
ambiguous gap fills, start with
[Troubleshooting]({{ '/troubleshooting/' | relative_url }}).

## Concept Guides

| Guide | Use it when |
| --- | --- |
| [How to interpret dot plots]({{ '/dot-plots/' | relative_url }}) | You need the visual grammar for reference/query alignment plots before deciding what action is justified. |
| [Choosing PAF or MUMmer coords]({{ '/guides/paf-vs-coords/' | relative_url }}) | You need to choose, compare, or sanity-check whole-genome alignment evidence. |
| [Inversions and orientation changes]({{ '/guides/inversions-orientation/' | relative_url }}) | You need to distinguish reverse orientation, internal inversions, complex same-reference patterns, and reference differences. |
| [Assembly graph evidence]({{ '/guides/assembly-graph-evidence/' | relative_url }}) | You need to read GFA segments, links, paths, walks, graph complexity labels, and graph guardrails. |
| [hifiasm unitig-to-contig projection]({{ '/guides/hifiasm-graph-projection/' | relative_url }}) | You need to understand why unitig graph coordinates must be projected before they can annotate contig-axis plots. |
| [Long-read PAF and GAF support]({{ '/guides/long-read-paf-gaf/' | relative_url }}) | You need to interpret read bridges, breakpoint support, graph traversals, MAPQ filters, and advisory evidence. |

## ChromoSort Guides

These are program-specific guides for running ChromoSort, reading its reports,
and deciding which command should or should not change sequence.

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

## Worked Examples

| Example | Use it when |
| --- | --- |
| [Mostly-correct assembly cleanup]({{ '/guides/mostly-correct-cleanup/' | relative_url }}) | You want a conservative `chromo clean` path from raw FASTA through validation plots. |
| [Suspected chimeric contig review]({{ '/guides/chimeric-contig-walkthrough/' | relative_url }}) | You want to go from a dot-plot signal to `eval fix`, manual review, reviewed splitting, re-alignment, and sorting. |
| [Graph-aware scaffold and gapfill review]({{ '/guides/graph-scaffold-gapfill-walkthrough/' | relative_url }}) | You want a combined GFA, GAF, Hi-C, reference-placement, scaffold, and gapfill review path. |
| [Dataset handoff checklist]({{ '/guides/dataset-handoff-checklist/' | relative_url }}) | You need to package FASTA, alignments, graph evidence, review tables, plots, summaries, and notes for a new review session. |

## Troubleshooting Paths

Some guides are most useful when a run looks wrong. The troubleshooting page is
now the shortest entry point for those cases.

| Symptom | Start here |
| --- | --- |
| Empty or sparse plots | [Plots are empty or sparse]({{ '/troubleshooting/#plots-are-empty-or-sparse' | relative_url }}) |
| Edited FASTA with old coords/PAF | [Edited FASTA does not match old alignments]({{ '/troubleshooting/#edited-fasta-does-not-match-old-alignments' | relative_url }}) |
| Coords and PAF disagree | [Coords and PAF disagree]({{ '/troubleshooting/#coords-and-paf-disagree' | relative_url }}) |
| GFA nodes or overlays are missing | [GFA overlay or projection is empty]({{ '/troubleshooting/#gfa-overlay-or-projection-is-empty' | relative_url }}) |
| Gapfill candidates are unresolved | [Gapfill plan is ambiguous]({{ '/troubleshooting/#gapfill-plan-is-ambiguous' | relative_url }}) |
