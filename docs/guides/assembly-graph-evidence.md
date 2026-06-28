---
title: Assembly Graph Evidence
description: Understand how ChromoSort uses GFA graph records, graph context, noseq graphs, and report-only graph guardrails.
---

# Assembly Graph Evidence

Use this guide when a ChromoSort command accepts `--gfa` and you need to know
what that graph evidence can and cannot decide.

The main question is:

> Is the graph describing context, confirming a reviewed junction, or supplying
> sequence that can actually be applied?

## The Core Idea

Assembly graphs are powerful because they keep adjacency, orientation, paths,
walks, coverage tags, and sometimes sequence. They are also easy to overread.
ChromoSort treats graph evidence conservatively: most graph-aware modes report
context, while only explicit graph-fill application can insert graph sequence.

The safe mental model is:

```text
GFA topology -> review context
GFA sequence + validated path + reviewed apply -> possible FASTA change
```

## GFA Records ChromoSort Reads

| GFA record | What it represents | How ChromoSort uses it |
| --- | --- | --- |
| `S` | Segment or node. May contain sequence or `*` plus `LN:i:length`. | Node names, lengths, optional sequence, coverage tags, and node-level context. |
| `L` | Oriented link between two segments, with an overlap. | Direct adjacency, graph-neighborhood reports, path search, and overlap context. |
| `P` | Named path through oriented segments. | Unitig-to-contig projection when path names match contig FASTA names. |
| `W` | Named walk through oriented segments. | Unitig-to-contig projection for assemblies that encode contigs as walks. |

Segments with sequence are required when `chromo gapfill --apply` may insert
graph sequence. Report-only commands can use no-sequence graphs when segment
lengths are available through `LN:i`.

## Full GFA Versus `.noseq.gfa`

| Graph type | Good for | Limitation |
| --- | --- | --- |
| Full sequence GFA | Gapfill application, flank validation, graph sequence insertion. | Larger files; more expensive to move around. |
| `.noseq.gfa` with `LN:i` lengths | Topology review, graph-map projection, boundary overlays, graph neighborhood context. | Cannot supply inserted sequence for gapfill. |
| `.noseq.gfa` without lengths | Link topology only. | Coordinates cannot be projected reliably. |
| GFA without `P` or `W` records | Local node/link context. | Cannot project unitig coordinates onto contig FASTA coordinates. |

Use the graph closest to the FASTA stage being reviewed. If the FASTA has been
renamed, split, manually exported, or scaffolded, an old graph may no longer
match its record names or flank sequence.

## Where Graph Evidence Appears

| Command | Graph input | What graph evidence can do |
| --- | --- | --- |
| `chromo sort --gfa` | Optional GFA. | Adds report-only graph assignment context. |
| `chromo manual --gfa` | Optional GFA. | Shows node badges, graph complexity, neighbors, links, overlaps, and same-reference context. |
| `chromo eval fix --gfa` | Optional GFA. | Adds graph node, unitig boundary, and local topology fields to review rows. |
| `chromo fix --gfa` | Optional GFA. | Writes graph context for fix review; split decisions remain alignment/review driven. |
| `chromo scaffold --gfa` | Optional GFA. | Writes report-only junction context by default. |
| `chromo eval scaffold --gfa` | Optional GFA. | Adds direct-link and short-path context to scaffold review rows. |
| `chromo graph-map` | Required GFA. | Projects unitig/path coordinates to contig coordinates when possible. |
| `chromo plot --gfa-overlay` | Optional GFA. | Draws projected graph intervals on the query axis. |
| `chromo gafprep --assembly-gfa` | Required GFA. | Writes a GraphAligner-oriented GFA for targeted read-to-graph alignment, preserving sequence while dropping audited hifiasm `A` records and pathological full-consuming links. |
| `chromo gapfill --gfa` | Required GFA. | Plans and, with reviewed `--apply`, can insert validated graph path sequence. |

## Report-Only Versus Sequence-Changing

Graph evidence is report-only unless a command explicitly says it changes
sequence.

| Evidence use | Sequence-changing? |
| --- | --- |
| Manual dashboard graph neighborhood panel | No |
| Eval `graph_*` columns | No |
| Sort graph assignment context | No |
| Scaffold `graph_gaps.tsv` with default graph policy | No |
| `gafprep` sanitized GraphAligner GFA and audit tables | No |
| Scaffold graph-confirmed terminal overlap trim | Yes, only with explicit `--graph-overlap-policy confirm` and compatible overlap policy |
| Gapfill plan TSV | No |
| `chromo gapfill --apply` with fillable accepted paths | Yes |

This separation is intentional. Graph context can make a row suspicious or more
believable, but the sequence-changing step stays visible.

## Graph Complexity Labels

Graph reports and dashboards highlight complexity that should slow review:

| Signal | Meaning |
| --- | --- |
| Missing graph node | FASTA or assignment names do not resolve to GFA `S` records. |
| Direct edge | Adjacent records have an orientation-aware GFA link. |
| Short path | A small explicit GFA path connects the records within the configured search depth. |
| Branching | Multiple candidate graph paths exist. |
| High-degree node | A node participates in many graph links and may be repeat-like. |
| Self-loop node | The graph has a node connected to itself. |
| Unsequenced node | A path uses a segment without sequence, preventing sequence insertion. |
| Cycle guard | Path enumeration avoided cycles while searching. |

These labels do not mean "bad." They mean the graph is asking for review
instead of blind application.

## Practical Review Workflow

1. Confirm GFA `S` names match the FASTA stage or assignment rows being used.
2. Check whether sequence is needed. Use full sequence GFA for gapfill
   application; `.noseq.gfa` is often fine for topology review.
3. If contig dot plots use `p_ctg.fa` but graph evidence is unitig-level, run
   `chromo graph-map` or use `chromo plot --gfa-overlay`.
4. Read graph report fields beside alignment, read, and audit-table evidence.
5. Treat graph conflicts as review prompts.
6. Use `chromo gapfill` only when graph sequence should actually replace Ns.

Preflight the graph:

```bash
awk 'BEGIN{FS="\t"} $1=="S" {print $2; count++} count==10 {exit}' assembly_graph.gfa
awk 'BEGIN{FS="\t"} $1=="L" {print $2, $3, $4, $5, $6; count++} count==5 {exit}' assembly_graph.gfa
awk 'BEGIN{FS="\t"} $1=="P" || $1=="W" {print; count++} count==5 {exit}' assembly_graph.gfa
```

## Cheat Sheet

| If you need... | Use... |
| --- | --- |
| Local graph neighborhood context | `chromo manual --gfa` or `chromo scaffold --gfa`. |
| Unitig boundaries on contig dot plots | `chromo plot --gfa-overlay` or `chromo graph-map`. |
| Graph context in spreadsheet review | `chromo eval fix/scaffold/gapfill --gfa`. |
| Sequence through an N gap | `chromo gapfill --gfa --apply` after reviewed planning. |
| Topology only from hifiasm | `.noseq.gfa` with `LN:i` lengths and path/walk records. |

## Common Traps

Do not expect graph node names to survive arbitrary renaming. Compare GFA `S`
names with FASTA IDs and ChromoSort `new_name` values.

Do not use a no-sequence graph for gapfill application. It can report topology
but cannot supply inserted bases.

Do not treat direct graph links as automatic permission to trim, reorder, or
fill. Those are separate sequence-changing choices.

Do not project unitig coordinates onto contig plots unless the GFA has matching
`P` or `W` records or direct segment-to-contig names.

Do not resolve conflicting evidence silently. Graph, GAF, Hi-C, and reference
placement can disagree; ChromoSort keeps those cases reviewable.

## What To Look At Next In ChromoSort

- Use [hifiasm Unitig-To-Contig Projection]({{ '/guides/hifiasm-graph-projection/' | relative_url }})
  when unitig GFA coordinates need to line up with contig FASTA plots.
- Use [Long-Read PAF And GAF Support]({{ '/guides/long-read-paf-gaf/' | relative_url }})
  to understand read evidence around graph paths and breakpoints.
- Use [Graph-Supported Gap Filling]({{ '/guides/graph-gap-filling/' | relative_url }})
  before applying graph path sequence.
- Use [Input Files]({{ '/input-files/' | relative_url }}#graph-input-files)
  for graph input contracts.
