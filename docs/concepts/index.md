---
title: Concepts
description: Portable genome-assembly and sequence-analysis concepts used by ChromoSort.
---

# Concepts

These guides teach ideas that are useful beyond ChromoSort. Start here when
you are learning how to read alignment plots, compare evidence files, or reason
about graph and long-read support before choosing any particular command.

| Concept | Use it when |
| --- | --- |
| [How to interpret dot plots]({{ '/dot-plots/' | relative_url }}) | You need the visual grammar for reference/query alignment plots before deciding what action is justified. |
| [Choosing PAF or MUMmer coords]({{ '/guides/paf-vs-coords/' | relative_url }}) | You need to choose, compare, or sanity-check whole-genome alignment evidence. |
| [Inversions and orientation changes]({{ '/guides/inversions-orientation/' | relative_url }}) | You need to distinguish reverse orientation, internal inversions, complex same-reference patterns, and reference differences. |
| [Assembly graph evidence]({{ '/guides/assembly-graph-evidence/' | relative_url }}) | You need to read GFA segments, links, paths, walks, graph complexity labels, and graph guardrails. |
| [hifiasm unitig-to-contig projection]({{ '/guides/hifiasm-graph-projection/' | relative_url }}) | You need to understand why unitig graph coordinates must be projected before they can annotate contig-axis plots. |
| [Long-read PAF and GAF support]({{ '/guides/long-read-paf-gaf/' | relative_url }}) | You need to interpret read bridges, breakpoint support, graph traversals, MAPQ filters, and advisory evidence. |

When the question becomes "which ChromoSort command should act on this
evidence?", move to the [Guides]({{ '/guides/' | relative_url }}).
