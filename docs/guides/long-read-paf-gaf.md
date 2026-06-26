---
title: Long-Read PAF And GAF Support
description: Interpret long-read-to-assembly PAF and read-to-graph GAF support fields in ChromoSort review workflows.
---

# Long-Read PAF And GAF Support

Use this guide when review tables include `longread_*` or `gaf_*` fields and
you need to understand what those counts can decide.

The main question is:

> Are reads supporting a breakpoint, a contig-end bridge, or a graph traversal,
> and is that support advisory or sequence-changing?

## The Core Idea

ChromoSort uses long-read evidence in two different coordinate systems:

- long-read-to-assembly PAF explains how reads align to contig sequences,
- read-to-graph GAF explains how reads traverse GFA graph nodes.

Both evidence streams are most useful in review tables and manual dashboards.
They can make a candidate more believable, point to an alternate graph path, or
warn that evidence is weak or tied. They do not generally change sequence by
themselves. In `chromo gapfill`, unique non-conflicting GAF support can help
choose among candidate graph paths before the chosen path still has to pass
sequence validation.

## Long-Read-To-Assembly PAF

Create a read-to-assembly PAF with long reads as queries and assembly contigs as
targets:

```bash
minimap2 -x map-hifi -c -t 16 --secondary=no assembly.fa reads.fastq.gz \
  > reads/sample.reads_to_assembly.paf
```

ChromoSort indexes those alignments by read and contig, filters by MAPQ and
identity settings, and summarizes support around candidate events.

### Breakpoint Support

`chromo eval fix --read-paf` can add:

| Field | Meaning |
| --- | --- |
| `longread_breakpoint_position` | Candidate cut coordinate on the source contig. |
| `longread_spanning_reads` | Reads with alignment anchors on both sides of the candidate breakpoint. |
| `longread_split_reads` | Reads represented by split alignments around the breakpoint. |
| `longread_left_edge_reads` | Reads ending or clipping near the left side of the event. |
| `longread_right_edge_reads` | Reads ending or clipping near the right side of the event. |
| `longread_nearby_reads` | Reads near the review window. |

Spanning reads can support continuity through a candidate breakpoint. Split and
edge reads can support a break, but they can also reflect repeats, clipping, or
alignment artifacts. Use them with the dot plot and graph context.

### Contig-End Bridge Support

`chromo eval scaffold --read-paf` and `chromo eval gapfill --read-paf` can add:

| Field | Meaning |
| --- | --- |
| `longread_bridge_reads` | Reads that anchor near both contig ends in the junction. |
| `longread_orientation_summary` | Orientation pattern counts for the bridging alignments. |
| `longread_read_order_summary` | Read-space order counts, such as left-before-right. |
| `longread_median_read_gap_bp` | Median read-space gap estimate when bridge reads support one. |

Bridge support is assembly-coordinate evidence. It can support adjacency, but
it does not insert graph sequence or override a stale reviewed plan.

## Read-To-Graph GAF

Create GAF with a graph aligner:

```bash
GraphAligner \
  -g assembly_graph.gfa \
  -f reads.fastq.gz \
  -a graph_alignments/sample.reads_to_graph.gaf
```

GAF path strings encode oriented graph traversal, such as:

```text
>left>bridge_good>right
<right<bridge_good<left
```

ChromoSort parses those paths, filters by `--min-gaf-mapq`, and counts reads
that contain selected or alternate oriented graph paths.

## GAF Support Status

Review tables can include:

| Field | Meaning |
| --- | --- |
| `gaf_path_nodes` | Selected path in oriented node notation. |
| `gaf_path_support` | Reads supporting the selected path. |
| `gaf_best_alt_path_nodes` | Best alternate path, when available. |
| `gaf_best_alt_support` | Reads supporting the best alternate. |
| `gaf_support_status` | Compact comparison of selected and alternate support. |
| `gaf_selected_reads` | Supporting read names, when written. |

Status values include:

| Status | Interpretation |
| --- | --- |
| `supports_selected` | One path has unique support and it is the selected path. |
| `supports_alternate` | Another candidate path has stronger unique support. |
| `tied_support` | More than one path has the same best support. |
| `weak_support` | Support exists but is below the configured threshold. |
| `no_support` | Candidate paths have no supporting GAF reads. |
| `no_paths` or `no_graph_path` | No graph path was available for comparison. |
| `missing_gfa` or missing-node statuses | Graph inputs did not resolve the needed nodes. |

## Where Read Evidence Can Affect Actions

| Workflow | Long-read PAF role | GAF role |
| --- | --- | --- |
| `eval fix` | Breakpoint support fields for review. | Advisory node/traversal context. |
| `manual fix` | Displays existing review-table evidence panels. | Displays existing review-table evidence panels. |
| `eval scaffold` | Contig-end bridge support fields. | Selected versus alternate graph path support. |
| `manual scaffold` | Displays bridge and graph support fields. | Displays path support fields. |
| `eval gapfill` | Contig-end bridge support fields. | Candidate graph path support. |
| `gapfill --gaf` | Not used to insert sequence directly. | Can resolve an otherwise ambiguous graph branch only with unique support above threshold. |

Even when GAF helps select a candidate path in gapfill, the path still must
have sequence, valid overlaps, matching FASTA flanks, and acceptable fill
length.

## Practical Review Workflow

1. Confirm which PAF is which. Reference-to-assembly PAF is not the same as
   long-read-to-assembly PAF.
2. Check MAPQ thresholds and whether secondary alignments were included.
3. For breakpoints, compare spanning and split-read counts to the dot plot.
4. For scaffold junctions, compare bridge reads to inferred gap or overlap.
5. For graph paths, compare selected and alternate GAF support.
6. If GAF, Hi-C, and reference-placement support point to different paths,
   leave the junction unresolved for review.

## Cheat Sheet

| Evidence | Best use |
| --- | --- |
| Read-to-assembly PAF near breakpoint | Check whether reads span, split, or end near a candidate cut. |
| Read-to-assembly PAF at contig ends | Check whether reads bridge adjacent sorted contigs. |
| GAF selected path support | Check whether reads traverse a candidate graph path. |
| GAF alternate support | Detect when another graph branch has stronger support. |
| GAF tied or weak support | Keep the event reviewable; do not force a path. |
| Conflicting evidence sources | Refuse automatic resolution and review manually. |

## Common Traps

Do not pass the reference-to-assembly PAF as `--read-paf`. Read PAF uses reads
as queries and assembly contigs as targets.

Do not assume many nearby reads equal breakpoint support. Look for the type of
support: spanning, split, edge, or bridge.

Do not treat GAF support as sequence validation. Gapfill still checks graph
sequence and flank matches.

Do not compare GAF node names to FASTA names unless the graph and FASTA are at
the same naming stage.

Do not let one evidence source override a conflict. ChromoSort keeps
conflicting graph-path support unresolved.

## What To Look At Next In ChromoSort

- Use [Assembly Graph Evidence]({{ '/guides/assembly-graph-evidence/' | relative_url }})
  for GFA and graph-context basics.
- Use [Graph-Supported Gap Filling]({{ '/guides/graph-gap-filling/' | relative_url }})
  for how GAF support can affect graph path selection.
- Use [Spreadsheet Review Tables]({{ '/guides/spreadsheet-review-tables/' | relative_url }})
  for accepting or rejecting rows that contain read evidence.
- Use [Input Files]({{ '/input-files/' | relative_url }}#creating-gaf-read-to-graph-alignments)
  for GAF creation notes.
