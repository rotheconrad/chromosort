---
title: Graph-Supported Gap Filling
description: Review chromo gapfill path enumeration, support evidence, risk flags, flank validation, reviewed application, and unresolved fallbacks.
---

# Graph-Supported Gap Filling

Use this guide when a final sorted scaffold has N gaps and an assembly graph may
contain sequence that should replace some of them.

The main question is:

> Is this graph path unique, supported, sequenced, validated against the FASTA
> flanks, and explicitly accepted?

## The Core Idea

`chromo gapfill` is the graph sequence application command. It is separate from
`chromo scaffold --gfa`, which reports graph junction evidence but does not
insert graph sequence.

Gap filling is intentionally gated:

```text
ordered FASTA + assignment TSV + GFA
  -> enumerate candidate graph paths
  -> score support and risk
  -> review plan
  -> apply only accepted fillable paths
```

Without `--apply`, no FASTA is written. With `--apply`, unresolved junctions
fall back to inferred or fixed N gaps.

## What Gapfill Reads Or Writes

`chromo gapfill` reads:

- final `chromo sort` ordered FASTA,
- matching `contig_assignments.tsv`,
- GFA with segment sequence for applied fills,
- optional GAF read-to-graph alignments,
- optional Hi-C graph-node contact table,
- optional reference-placement PAF for graph nodes,
- optional reviewed plan.

It writes:

- `<prefix>.gapfill_plan.tsv`,
- optional review HTML,
- optional `<prefix>.gapfilled.fa` with `--apply`,
- `<prefix>.run_summary.txt`.

## Path Status Gallery

| Status or reason | Meaning | Usual action |
| --- | --- | --- |
| `fillable` | Selected graph path passed sequence and flank validation. | Review, then accept if the biological context fits. |
| `missing_node` | One or both scaffold flanks did not resolve to GFA nodes. | Check names and graph stage. |
| `no_graph_path` | No path was found within the configured search depth. | Leave Ns or inspect graph inputs. |
| `ambiguous_paths` | More than one candidate path remains. | Use GAF, Hi-C, reference-placement evidence, or manual review. |
| `left_flank_sequence_mismatch` or `right_flank_sequence_mismatch` | GFA flank sequence does not match ordered FASTA. | Stop; graph and FASTA likely come from different stages. |
| `unsequenced_flank_node` or `unsequenced_intermediate_node` | A needed GFA segment lacks sequence. | Use a full sequence GFA or leave unresolved. |
| `conflicting_gaf_hic_support` or related conflict | Evidence sources support different paths. | Keep unresolved and review manually. |
| oversized or invalid overlap reasons | Fill cannot be reconstructed safely. | Review graph and thresholds. |

## Candidate Path Evidence

When multiple graph paths exist, ChromoSort can resolve a branch only when one
candidate has unique non-conflicting support.

| Evidence source | What it scores | Resolution rule |
| --- | --- | --- |
| GAF | Reads traversing candidate graph paths. | One path must have unique support at or above `--min-gaf-path-support`. |
| Hi-C pairs | Summed graph-node contacts along candidate paths. | One path must have unique support at or above `--min-hic-path-support`. |
| Reference-placement PAF | Intermediate graph nodes placed inside the expected reference-space gap. | One path must have uniquely stronger support at or above `--min-ref-path-support`. |

If different evidence sources uniquely support different paths, ChromoSort
keeps the junction unresolved.

## Risk Flags

Read `risk_flags` before accepting a fill:

| Risk flag | Why it matters |
| --- | --- |
| `branching` | More than one graph path exists. |
| `high_degree` | Candidate path includes graph nodes with many links. |
| `self_loop` | Candidate path includes self-loop nodes. |
| `unsequenced` | Candidate path includes nodes without sequence. |
| `cycles_avoided` | Path search encountered cycles and avoided them. |
| `weak_gaf_support`, `weak_hic_support`, `weak_ref_paf_support` | Evidence exists but does not meet threshold. |
| `tied_gaf_support`, `tied_hic_support`, `tied_ref_paf_support` | Candidate paths tie under that evidence source. |
| `conflicting_support` | Evidence sources disagree. |
| `sequence_validation_failed` | Path selection happened, but sequence reconstruction failed. |

Risk flags do not always mean "reject." They mean "review with care."

## Plan First

```bash
chromo gapfill \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --gfa assembly_graph.gfa \
  --gaf reads_to_graph.gaf \
  --hic-pairs graph_contacts.tsv \
  --ref-paf paf/graph_nodes_to_ref.paf \
  --output-prefix results/sample.gapfill \
  --review-html results/sample.gapfill.review.html \
  --include-fill-sequences
```

Planning writes the TSV and optional review HTML. It does not write a filled
FASTA. The HTML can filter rows, compare candidate paths, toggle accepted
fillable paths, and export a reviewed-plan TSV.

## Apply Reviewed Fills

After review, apply only accepted rows:

```bash
chromo gapfill \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --gfa assembly_graph.gfa \
  --gaf reads_to_graph.gaf \
  --hic-pairs graph_contacts.tsv \
  --ref-paf paf/graph_nodes_to_ref.paf \
  --reviewed-plan results/sample.gapfill.reviewed_plan.tsv \
  --output-prefix results/sample.reviewed_gapfill \
  --apply
```

With a reviewed plan, ChromoSort rechecks the current scaffold, contig pair,
path nodes, and fillability before applying the row. If the graph path has
changed, a stale accepted row fails instead of being silently applied.

For fillable paths, ChromoSort inserts graph sequence after the left flank and
trims the right flank by the final GFA overlap so the graph path is not
duplicated at the join.

## Read The Gapfill Plan

| Column family | What to inspect |
| --- | --- |
| Junction IDs | `scaffold`, `left_contig`, `right_contig`. |
| Graph resolution | `graph_status`, `path_nodes`, `candidate_paths`, `intermediate_nodes`. |
| Support | `gaf_*`, Hi-C support, and `ref_*` support columns. |
| Risk | `risk_flags`, `branch_complexity_score`, high-degree, self-loop, and unsequenced node lists. |
| Sequence action | `fill_status`, `fill_bp`, `right_trim_bp`, `fallback_gap_bp`. |
| Review state | `accept_fill` or shared review-table `accept`, plus `applied`. |

Keep the plan TSV, reviewed plan, run summary, and gapfilled FASTA together.

## Cheat Sheet

| If you see... | Think... | Action |
| --- | --- | --- |
| `fillable` with low risk | Candidate can be accepted if the biology fits. | Review and apply with `--reviewed-plan --apply`. |
| `ambiguous_paths` | Graph has more than one possible bridge. | Add support evidence or leave Ns. |
| Strong alternate GAF support | Selected path may not be best. | Review candidate details. |
| `conflicting_support` | Evidence sources disagree. | Do not auto-fill. |
| `unsequenced` | Needed graph node lacks sequence. | Use full GFA or leave unresolved. |
| Flank mismatch | Graph and ordered FASTA disagree. | Regenerate evidence from the correct stage. |
| Reviewed plan rejected as stale | Inputs changed since review. | Regenerate the plan. |

## Common Traps

Do not assume `chromo scaffold --gfa` filled gaps. It only reports graph
junction context.

Do not apply all fillable rows blindly on a new dataset. Review the plan first,
especially rows with risk flags.

Do not use a no-sequence GFA for applied gap filling.

Do not let GAF, Hi-C, or reference-placement evidence resolve a branch when
they conflict. ChromoSort deliberately refuses those cases.

Do not use an old reviewed plan after changing ordered FASTA, assignments, GFA,
GAF, PAF, or path-search settings.

## What To Look At Next In ChromoSort

- Use [Assembly Graph Evidence]({{ '/guides/assembly-graph-evidence/' | relative_url }})
  for graph context and GFA record basics.
- Use [Long-Read PAF And GAF Support]({{ '/guides/long-read-paf-gaf/' | relative_url }})
  for graph traversal and bridge evidence.
- Use [Spreadsheet Review Tables]({{ '/guides/spreadsheet-review-tables/' | relative_url }})
  for reviewed `fill_path` rows.
- Use [chromo gapfill]({{ '/commands/gapfill/' | relative_url }}) for the full
  parameter and output reference.
