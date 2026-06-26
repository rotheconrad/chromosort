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

or:

scaffold FASTA + AGP + GFA
  -> validate AGP/component/N-gap layout
  -> enumerate candidate graph paths for AGP gaps
  -> score support and risk
  -> review plan
  -> apply only accepted fillable paths
```

Without `--apply`, no FASTA is written. With `--apply`, ChromoSort now requires
one of two explicit application modes: `--reviewed-plan` for accepted reviewed
rows, or `--apply-all-fillable` for deliberate exploratory or benchmark runs
that apply every currently fillable path. Unresolved or unaccepted junctions
fall back to inferred, AGP, or fixed N gaps.

## What Gapfill Reads Or Writes

`chromo gapfill` reads:

- final `chromo sort` ordered FASTA plus matching `contig_assignments.tsv`,
  or existing scaffold FASTA plus AGP,
- GFA with segment sequence for applied fills,
- optional GAF read-to-graph alignments,
- optional Hi-C graph-node contact table,
- optional reference-placement PAF for graph nodes,
- optional long-read-to-contig PAF for bridge evidence,
- optional external patch candidate TSV/FASTA for graph-fill comparison,
- optional reviewed plan.

It writes:

- `<prefix>.gapfill_plan.tsv`,
- `<prefix>.gapfilled.agp`,
- `<prefix>.gapfilled_components.tsv`,
- optional review HTML,
- optional `<prefix>.gapfilled.fa` with `--apply`,
- `<prefix>.submission_checklist.tsv`,
- `<prefix>.run_summary.txt`.

The submission checklist is a review aid, not a replacement for NCBI validation.
In planning mode it marks FASTA-dependent checks as not checked because
`gapfilled.fa` has not been written yet. After `--apply`, it checks the final
FASTA against the AGP map and reports unresolved gap and graph-fill totals.

`chromo gapfill` can be run instead of `chromo scaffold` when you still have the
ordered FASTA, matching assignments, and a graph that resolves the ordered
contigs. It writes scaffold-like records itself, filling accepted graph paths
and using N gaps elsewhere. If you already have a scaffold FASTA, run with
`--scaffold-fasta` and `--agp`. The AGP is required because scaffold FASTA
alone loses the component identities around each N gap. In scaffold/AGP mode,
AGP rows must cover each scaffold object exactly, and gap rows must match
N-only spans in the scaffold FASTA.

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
| Long-read PAF | Reads with terminal anchors on both adjacent contigs. | Reported as bridge evidence only; does not insert sequence or choose a path by itself. |
| External patch table | Candidate patch sequences from outside tools keyed to the same scaffold/component gap. | Reported as concordance or conflict with the selected graph fill; does not choose or insert a path. |

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
| `multiple_external_patches` | More than one imported patch candidate maps to this gap. |
| `external_patch_graph_mismatch` | The best imported patch does not match the graph-derived fill sequence. |
| `external_patch_only` | A patch candidate exists but no graph fill sequence is available for comparison. |

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
  --read-paf paf/reads_to_ordered_contigs.paf \
  --patch-table results/patch_candidates.tsv \
  --patch-fasta results/patch_candidates.fa \
  --output-prefix results/sample.gapfill \
  --review-html results/sample.gapfill.review.html \
  --include-fill-sequences
```

Planning writes the TSV and optional review HTML. It does not write a filled
FASTA. The HTML can filter rows, compare candidate paths, toggle accepted
fillable paths, and export a reviewed-plan TSV.

For scaffold-first workflows, use the scaffold FASTA and AGP instead of the
ordered FASTA and assignment table:

```bash
chromo gapfill \
  --scaffold-fasta results/sample.scaffold.fa \
  --agp results/sample.scaffold.agp \
  --gfa assembly_graph.gfa \
  --read-paf paf/reads_to_agp_components.paf \
  --output-prefix results/sample.scaffold_gapfill \
  --review-html results/sample.scaffold_gapfill.review.html \
  --include-fill-sequences
```

`--read-paf` remains component-level bridge evidence in this mode: PAF target
names should match AGP component IDs. A PAF whose targets are only scaffold
record names is not currently interpreted as reads spanning N blocks.

For unitig-level assembler graphs whose `P` path or `W` walk records are named
for the ordered contigs or AGP components, add `--project-gfa-paths`:

```bash
chromo gapfill \
  --scaffold-fasta results/sample.scaffold.fa \
  --agp results/sample.scaffold.agp \
  --gfa hifiasm.p_utg.noseq.gfa \
  --project-gfa-paths \
  --output-prefix results/sample.projected_gapfill \
  --review-html results/sample.projected_gapfill.review.html
```

Projected rows use terminal unitig names in `left_graph_node`,
`right_graph_node`, and `path_nodes`. With sequence-bearing GFA segments, rows
can become `fillable` after terminal unitig sequence validation and can be
applied through the normal reviewed apply path. With `.noseq.gfa`, rows remain
`projected_path_planning_only` topology evidence and are not applied.

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
  --read-paf paf/reads_to_ordered_contigs.paf \
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

For non-reviewed exploratory runs, add `--apply-all-fillable` to `--apply`.
That mode still refuses ambiguous or unverifiable paths, but it does not require
human acceptance.

## Read The Gapfill Plan

| Column family | What to inspect |
| --- | --- |
| Junction IDs | `scaffold`, `left_contig`, `right_contig`. |
| Graph resolution | `graph_status`, `path_nodes`, `candidate_paths`, `intermediate_nodes`. |
| Support | `gaf_*`, Hi-C support, and `ref_*` support columns. |
| Long-read bridge evidence | `longread_*` bridge count, orientation, read order, and median read-gap columns. |
| External patch comparison | `patch_candidate_count`, `patch_best_*`, and `patch_graph_status`. |
| Risk | `risk_flags`, `branch_complexity_score`, high-degree, self-loop, and unsequenced node lists. |
| Sequence action | `fill_status`, `fill_bp`, `right_trim_bp`, `fallback_gap_bp`. |
| Review state | `accept_fill` or shared review-table `accept`, plus `applied`. |
| Submission provenance | AGP/component rows describing ordered contigs, fallback gaps, graph-fill slices, and trimmed coordinates. |

Keep the plan TSV, reviewed plan, run summary, and gapfilled FASTA together.

## Cheat Sheet

| If you see... | Think... | Action |
| --- | --- | --- |
| `fillable` with low risk | Candidate can be accepted if the biology fits. | Review and apply with `--reviewed-plan --apply`. |
| `ambiguous_paths` | Graph has more than one possible bridge. | Add support evidence or leave Ns. |
| Strong alternate GAF support | Selected path may not be best. | Review candidate details. |
| `conflicting_support` | Evidence sources disagree. | Do not auto-fill. |
| `patch_graph_status=exact_graph_match` | An external patch agrees with the graph fill. | Useful supporting evidence, still review before apply. |
| `external_patch_graph_mismatch` | External patch sequence and graph-fill sequence disagree. | Inspect both sources and leave unresolved if uncertain. |
| `unsequenced` | Needed graph node lacks sequence. | Use full GFA or leave unresolved. |
| Flank mismatch | Graph and ordered FASTA disagree. | Regenerate evidence from the correct stage. |
| Reviewed plan rejected as stale | Inputs changed since review. | Regenerate the plan. |

## Common Traps

Do not assume `chromo scaffold --gfa` filled gaps. It only reports graph
junction context.

Do not apply all fillable rows blindly on a new dataset. `--apply` alone is
rejected; use reviewed application for production, and reserve
`--apply-all-fillable` for deliberate exploratory or benchmark runs.

Do not use a no-sequence GFA for applied gap filling.

Do not expect `.noseq.gfa` projection to apply sequence. `--project-gfa-paths`
can apply only when the projected path uses sequence-bearing GFA segments and
the terminal unitig sequence validates against the component ends.

Do not use scaffold FASTA alone for graph-aware filling. Keep AGP with the
scaffold FASTA so ChromoSort can recover component identities around each N gap.

Do not expect long-read PAF to fill a gap by itself. It is bridge evidence. The
inserted sequence must come from a validated graph path in the current
implementation.

Do not treat `--patch-table` as an external sequence application mode. Imported
patches are comparison evidence; ChromoSort still applies only accepted,
validated graph-path sequence.

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
