---
title: Graph-Aware Scaffold And Gapfill Review
description: A worked synthetic graph fixture walkthrough for scaffold graph reports, gapfill planning, review, and accepted graph path application.
---

# Graph-Aware Scaffold And Gapfill Review

Use this walkthrough when sorted contigs are ready for scaffolding and a graph
may explain one or more gaps between adjacent contigs.

The goal is:

> Compare GFA, GAF, Hi-C-like contacts, and reference-placement support before
> applying any graph sequence.

This walkthrough uses the repository fixture under `tests/data/graph_gotchas`.
It is tiny by design, but it exercises the same file types used in real runs.

## Fixture Files

| File | Role |
| --- | --- |
| `ref.fa` | Toy two-chromosome reference. |
| `assembly.fa` | Toy contig/unitig sequences for sort and manual review. |
| `unitigs.gfa` | Assembly graph with sequence, branches, a cycle, one length-only segment, and overlaps. |
| `unitig_to_ref.paf` | Unitig-to-reference placements. |
| `reads_to_graph.gaf` | Read-to-graph paths. |
| `hic_pairs.tsv` | Simple graph-node contact counts. |
| `gapfill_ordered.fa` | Tiny two-flank ordered FASTA for the gapfill example. |
| `gapfill_assignments.tsv` | Matching assignment table for `gapfill_ordered.fa`. |

The key graph choice is between `bridge_good` and `bridge_alt`. The expected
best path for the confident gap is:

```text
left+,bridge_good+,right+
```

## Step 1: Sort With Graph Context

```bash
mkdir -p results/graph_gotchas
DATA=tests/data/graph_gotchas

chromo sort \
  --ref-fasta "$DATA/ref.fa" \
  --assembly-fasta "$DATA/assembly.fa" \
  --paf "$DATA/unitig_to_ref.paf" \
  --gfa "$DATA/unitigs.gfa" \
  --output-prefix results/graph_gotchas/sort
```

Read:

- `sort.contig_assignments.tsv` for retained contigs and statuses,
- graph assignment context if written for the run,
- `sort.run_summary.txt` for inputs and thresholds.

## Step 2: Open A Manual Graph Dashboard

```bash
chromo manual \
  --ref-fasta "$DATA/ref.fa" \
  --assembly-fasta "$DATA/assembly.fa" \
  --paf "$DATA/unitig_to_ref.paf" \
  --gfa "$DATA/unitigs.gfa" \
  --embed-sequences \
  --output-html results/graph_gotchas/manual.html
```

Open `results/graph_gotchas/manual.html` and inspect:

- contig dot plots,
- graph node badges,
- simple, branching, self-loop, and missing graph filters,
- selected-node neighbors and oriented links.

This dashboard is review context. It does not fill gaps.

## Step 3: Scaffold With Graph Junction Reports

```bash
chromo scaffold \
  --ordered-fasta results/graph_gotchas/sort.ordered.fa \
  --assignments results/graph_gotchas/sort.contig_assignments.tsv \
  --gfa "$DATA/unitigs.gfa" \
  --output-prefix results/graph_gotchas/scaffold
```

Review:

| Output | What it tells you |
| --- | --- |
| `scaffold.scaffold_gaps.tsv` | Inferred gaps, overlaps, gap modes, and overlap actions. |
| `scaffold.graph_gaps.tsv` | Direct edges, short graph paths, missing nodes, and orientations. |
| `scaffold.scaffold_summary.tsv` | Scaffold lengths, gap totals, and ordered contig lists. |
| `scaffold.submission_checklist.tsv` | FASTA/AGP consistency, gap counts, and handoff checks. |

Graph scaffold reports are still report-only by default.

## Step 4: Plan Graph Fills

Use the focused two-flank gapfill fixture:

```bash
chromo gapfill \
  --ordered-fasta "$DATA/gapfill_ordered.fa" \
  --assignments "$DATA/gapfill_assignments.tsv" \
  --gfa "$DATA/unitigs.gfa" \
  --ref-paf "$DATA/unitig_to_ref.paf" \
  --gaf "$DATA/reads_to_graph.gaf" \
  --hic-pairs "$DATA/hic_pairs.tsv" \
  --output-prefix results/graph_gotchas/gapfill \
  --include-fill-sequences \
  --review-html results/graph_gotchas/gapfill.review.html
```

Open:

```text
results/graph_gotchas/gapfill.gapfill_plan.tsv
results/graph_gotchas/gapfill.review.html
```

In this toy case, the plan should show `bridge_good` as the fillable selected
path and `bridge_alt` as the weaker competing branch.

## Step 5: Review The Candidate Row

Read these columns together:

| Column | Expected lesson |
| --- | --- |
| `path_nodes` | The selected graph path, such as `left+,bridge_good+,right+`. |
| `candidate_paths` | Whether the graph had more than one possible bridge. |
| `gaf_path_support` and `gaf_best_alt_support` | Whether graph-aligned reads support the selected path or an alternate. |
| Hi-C support columns | Whether contacts support the same graph branch. |
| `ref_path_support` and `ref_best_alt_support` | Whether graph nodes place in the expected reference gap. |
| `risk_flags` | Branching, high-degree, self-loop, unsequenced, or conflict warnings. |
| `fill_status` | Whether the selected path is actually fillable. |
| `fill_bp` and `right_trim_bp` | How much graph sequence is inserted and how much right flank overlap is trimmed. |

Only accept rows where the graph path is fillable and the evidence fits the
biological question.

## Step 6: Apply Accepted Fills

After reviewing the plan or HTML, export a reviewed plan with
`accept_fill=yes` only for the accepted fillable row. Then apply:

```bash
chromo gapfill \
  --ordered-fasta "$DATA/gapfill_ordered.fa" \
  --assignments "$DATA/gapfill_assignments.tsv" \
  --gfa "$DATA/unitigs.gfa" \
  --ref-paf "$DATA/unitig_to_ref.paf" \
  --gaf "$DATA/reads_to_graph.gaf" \
  --hic-pairs "$DATA/hic_pairs.tsv" \
  --reviewed-plan results/graph_gotchas/gapfill.reviewed_plan.tsv \
  --output-prefix results/graph_gotchas/gapfill.reviewed \
  --apply \
  --simple-headers
```

Review:

- `gapfill.reviewed.gapfill_plan.tsv` for `applied=yes`,
- `gapfill.reviewed.gapfilled.fa` for the filled scaffold record,
- `gapfill.reviewed.submission_checklist.tsv` for FASTA/AGP handoff checks,
- `gapfill.reviewed.run_summary.txt` for fill-status counts.

## What This Example Teaches

| Scenario | Lesson |
| --- | --- |
| `confident_gap_path` | A fillable path can be applied after review when support agrees. |
| `ambiguous_branch` | Alternate graph paths should stay reviewable until support separates them. |
| `cycle_guard` | Path search avoids graph cycles rather than following them indefinitely. |
| `orientation_specific` | Oriented GFA links matter. |
| `disconnected_mapped_node` | Reference placement alone does not create a graph path. |
| `repeat_or_duplicate_warning` | Branches and repeat-like nodes deserve caution even when alignments exist. |

## Common Traps

Do not confuse scaffold graph reports with gap filling. `chromo scaffold --gfa`
does not insert graph sequence.

Do not use a graph-node PAF whose query names do not match GFA segment names.

Do not apply a reviewed plan after changing ordered FASTA, assignments, GFA,
GAF, Hi-C pairs, reference-placement PAF, or path-search settings.

Do not ignore alternate path support in the review HTML. The whole point of the
toy branch is to make the alternative visible.

## What To Look At Next In ChromoSort

- Use [Graph-Supported Gap Filling]({{ '/guides/graph-gap-filling/' | relative_url }})
  for the general decision rules.
- Use [Assembly Graph Evidence]({{ '/guides/assembly-graph-evidence/' | relative_url }})
  for GFA basics.
- Use [Long-Read PAF And GAF Support]({{ '/guides/long-read-paf-gaf/' | relative_url }})
  for read-path support fields.
- Use [Dataset Handoff Checklist]({{ '/guides/dataset-handoff-checklist/' | relative_url }})
  before passing a graph-aware run to another reviewer.
