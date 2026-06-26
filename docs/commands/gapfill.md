---
title: chromo gapfill
description: Usage, outputs, parameters, and reasoning for chromo gapfill.
---

# chromo gapfill

Use `chromo gapfill` after final sorting and manual review when a GFA graph
gives a validated sequence path between adjacent sorted contigs, or when an
otherwise ambiguous graph branch has unique non-conflicting support evidence.

## What `chromo gapfill` Does

Given either a final `chromo sort` ordered FASTA plus its matching assignment
report, or an existing scaffold FASTA plus AGP, together with a GFA assembly
graph and optional GAF graph alignments, Hi-C pair evidence,
reference-placement PAF, or long-read bridge evidence, `chromo gapfill`:

1. Builds adjacent flanks from retained sorted contigs, or from AGP component
   rows around AGP N-gap rows.
2. In scaffold/AGP mode, validates that AGP rows cover scaffold FASTA records
   contiguously, component lengths match their object spans, and AGP gap spans
   are N-only in the scaffold FASTA.
3. Resolves each flank to a GFA segment using original and renamed contig IDs,
   or, with `--project-gfa-paths`, projects component IDs through matching GFA
   `P`/`W` paths and uses the terminal unitigs for path search.
4. Enumerates graph paths up to `--max-path-edges`.
5. Uses GAF read-path support, Hi-C contact support, and reference-placement
   PAF support to resolve an otherwise ambiguous graph branch only when one
   candidate path has unique support above threshold and evidence sources do
   not conflict.
6. Reports optional long-read PAF bridge support across the adjacent contig
   ends. Long-read PAF is evidence; it does not provide inserted sequence.
7. Annotates candidate-path risk, including high-degree graph nodes, self-loop
   nodes, unsequenced nodes, cycle guards, weak/tied/conflicting support, and a
   branch-complexity score.
8. Rejects missing nodes, disconnected flanks, unresolved ambiguous paths,
   unknown or invalid overlaps, oversized fills, and flank sequence mismatches
   for applied fills. Projected unitig-GFA paths are currently reported for
   planning/review only.
9. Writes `<prefix>.gapfill_plan.tsv` for review with `accept_fill=no` by default,
   can write a self-contained HTML reviewer with `--review-html`, and can accept
   the shared review-event table from `chromo eval gapfill`.
10. Writes AGP and component-provenance tables describing ordered contig
    components, fallback N gaps, and any applied graph-fill components.
11. With `--apply`, writes `<prefix>.gapfilled.fa`. Application must be explicit:
   either provide `--reviewed-plan` to apply only accepted rows, or add
   `--apply-all-fillable` when you deliberately want every currently fillable
   path applied.

## Plan Graph Fills

```bash
chromo gapfill \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --gfa assembly_graph.gfa \
  --ref-paf paf/sample.ref_vs_asm.paf \
  --read-paf paf/reads_to_ordered_contigs.paf \
  --output-prefix results/sample.gapfill \
  --review-html results/sample.gapfill.review.html
```

Planning mode writes the gapfill plan but does not create a FASTA. Add
`--include-fill-sequences` when you want short candidate sequences embedded in
the TSV for manual review. To make application explicitly reviewed, edit the
`accept_fill` column from `no` to `yes` only for rows you want to apply, then
pass that edited table back with `--reviewed-plan`. When `--review-html` is
provided, the HTML table can filter rows, toggle accepted fillable paths, show
side-by-side candidate path comparisons for ambiguous branches, and export a
reviewed-plan TSV with the same columns.

For spreadsheet-first review using the shared review-event schema, run
`chromo eval gapfill` instead. It writes `<prefix>.gapfill_review.tsv` with
accepted `fill_path` rows that can also be passed to `--reviewed-plan`.

## Plan From An Existing Scaffold And AGP

```bash
chromo gapfill \
  --scaffold-fasta results/sample.scaffold.fa \
  --agp results/sample.scaffold.agp \
  --gfa assembly_graph.gfa \
  --read-paf paf/reads_to_agp_components.paf \
  --output-prefix results/sample.scaffold_gapfill \
  --review-html results/sample.scaffold_gapfill.review.html
```

This mode is for the standard scaffold-first workflow. The scaffold FASTA gives
the emitted sequence, and AGP restores the component identities around each
N-gap so those flanks can be matched back to the graph. AGP rows must cover each
scaffold object exactly; gap rows must correspond to N-only scaffold spans.
Adjacent AGP component rows without an intervening gap are carried forward as
zero-gap joins and are not treated as graph-fill candidates.

## Plan With Unitig-Level GFA Paths

```bash
chromo gapfill \
  --scaffold-fasta results/sample.scaffold.fa \
  --agp results/sample.scaffold.agp \
  --gfa hifiasm.p_utg.noseq.gfa \
  --project-gfa-paths \
  --output-prefix results/sample.projected_gapfill \
  --review-html results/sample.projected_gapfill.review.html
```

Use `--project-gfa-paths` when the GFA `S` records are unitigs but the ordered
or AGP component IDs appear as GFA `P` path or `W` walk names. ChromoSort
projects each component onto its path/walk, finds the right terminal unitig of
the left component and the left terminal unitig of the right component, and
then searches the graph between those unitigs.

With sequence-bearing GFA segments, projected paths can become normal `fillable`
rows and can be applied through the same reviewed or all-fillable apply gates as
direct contig-level graph paths. ChromoSort validates that the projected left
terminal unitig matches the left component suffix and the projected right
terminal unitig matches the right component prefix before marking the row
fillable. A `.noseq.gfa` with `LN:i` lengths remains useful for topology review;
those rows stay `projected_path_planning_only` and are not applied.

## How This Differs From `chromo scaffold`

`chromo scaffold` creates per-reference FASTA records with N gaps and can report
GFA context, but it does not insert graph sequence. `chromo gapfill` also creates
per-reference records from the ordered FASTA and assignment table, or it can
start from an existing scaffold FASTA and AGP. In both cases it can replace
selected junction gaps with validated GFA path sequence and trim the right flank
by the terminal graph overlap.

That means `chromo gapfill` can skip a prior `chromo scaffold` run when you
still have the final `ordered.fa`, matching `contig_assignments.tsv`, and a GFA
that resolves the ordered contigs. It is effectively "scaffold plus reviewed
graph fills." If you already have a scaffold FASTA, use `--scaffold-fasta` with
`--agp`; scaffold FASTA alone is too lossy for graph-aware filling because the
scaffold record no longer names the original components around each N gap.

## Apply Reviewed Graph Fills

```bash
chromo gapfill \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --gfa assembly_graph.gfa \
  --ref-paf paf/sample.ref_vs_asm.paf \
  --read-paf paf/reads_to_ordered_contigs.paf \
  --gaf reads_to_graph.gaf \
  --hic-pairs graph_contacts.tsv \
  --reviewed-plan results/sample.eval_gapfill.gapfill_review.tsv \
  --output-prefix results/sample.reviewed_gapfill \
  --apply
```

Applied mode requires either `--reviewed-plan` or `--apply-all-fillable`.
Reviewed application is the production path: only accepted, revalidated rows are
filled, and all other junctions fall back to Ns. `--apply-all-fillable` is an
explicit exploratory or benchmark mode that applies every currently fillable
path without human acceptance. In both modes, ambiguous or unverifiable paths
are refused. If GAF is provided,
an ambiguous GFA branch can be filled only when one candidate path has unique
support of at least `--min-gaf-path-support` reads after `--min-gaf-mapq`
filtering. If Hi-C pair evidence is provided, one candidate path must have
unique summed contact support of at least `--min-hic-path-support`. If
`--ref-paf` is provided, one candidate path can be chosen when its intermediate
graph nodes have uniquely stronger same-reference placement support inside the
expected reference-space gap. When evidence sources uniquely support different
paths, the branch remains unresolved. When `--reviewed-plan` is used,
ChromoSort accepts either the legacy gapfill-plan TSV with `accept_fill=yes` or
the shared `chromo eval gapfill` review-event TSV with `accept=yes`. It rechecks
the current scaffold, contig pair, and `path_nodes` before applying an accepted
row, so stale reviewed paths fail instead of being applied. For a fillable path,
ChromoSort inserts the graph sequence after the left flank and trims the right
flank prefix by the final GFA overlap so the joined sequence follows the graph
path without duplicating the overlap. Unfilled junctions receive the inferred
reference-space N gap in ordered-FASTA mode, the AGP gap length in scaffold/AGP
mode, or `--fixed-gap-bp` when provided.

`--read-paf` should be long reads mapped back to the pre-scaffold contig FASTA
used by the ordered FASTA, or to records whose names match the `contig`,
`new_name`, or FASTA record names in the assignment table. In scaffold/AGP mode,
the PAF target names should match AGP component IDs. ChromoSort currently treats
long-read PAF as component-to-component bridge evidence; it does not interpret a
scaffold-record-only PAF as reads spanning N blocks. ChromoSort counts reads
with terminal anchors on both sides of each candidate junction and reports
bridge counts, orientation summaries, read-order summaries, and median read-gap
estimates. These fields can support or question a fill decision, but they do
not choose a graph branch by themselves and do not provide inserted bases.

## Compare External Patch Candidates

```bash
chromo gapfill \
  --scaffold-fasta results/sample.scaffold.fa \
  --agp results/sample.scaffold.agp \
  --gfa assembly_graph.gfa \
  --patch-table results/patch_candidates.tsv \
  --patch-fasta results/patch_candidates.fa \
  --output-prefix results/sample.patch_review \
  --review-html results/sample.patch_review.html
```

Use `--patch-table` to import candidate gap patches from external tools such as
TGS-GapCloser, LR_Gapcloser, Sealer, or RagTag patch as review evidence. The
table is keyed by `scaffold`, `left_contig`, and `right_contig`; in
scaffold/AGP mode these are AGP component IDs around the gap. Rows can include
`patch_sequence` directly, or `patch_id` plus `--patch-fasta`.

Minimum columns:

```text
scaffold  left_contig  right_contig  patch_id  source  patch_sequence
```

Accepted aliases include `object` for `scaffold`, `left_component` and
`right_component`, `tool` for `source`, and `sequence` or `seq` for
`patch_sequence`. The gapfill plan reports candidate count, best patch ID,
source, length, notes, and whether the patch exactly matches, length-matches
but disagrees with, or differs from the selected graph fill sequence.

External patches are not inserted by `chromo gapfill` in this mode. They are
evidence for review and comparison; the inserted sequence source remains an
accepted, validated graph path.

## `chromo gapfill` Outputs

| Output | Description |
| --- | --- |
| `<prefix>.gapfill_plan.tsv` | One row per adjacent sorted contig pair, or adjacent AGP component pair in scaffold/AGP mode, with graph status, path nodes, GAF support counts/status/supporting reads, Hi-C and reference-placement support counts, long-read bridge fields, external patch comparison fields, risk flags, branch-complexity score, high-degree/self-loop/unsequenced node lists, fill status, inserted bp, right-trim bp, fallback gap bp, editable `accept_fill`, and whether the fill was applied. |
| `<prefix>.gapfilled.agp` | AGP 2.1 map for the current gapfill output model, including ordered or AGP-derived components, fallback N gaps, and graph-fill slices. In apply mode it describes the emitted `gapfilled.fa`; in planning mode it describes the all-fallback scaffold map. |
| `<prefix>.gapfilled_components.tsv` | Human-readable provenance table mirroring the AGP rows, including ordered or AGP-derived components, fallback gaps, applied graph-fill segment slices, trimmed right-flank coordinates, and status fields. |
| `--review-html` path | Optional self-contained HTML table for reviewing gapfill-plan rows, comparing candidate paths, and exporting a reviewed-plan TSV. |
| `<prefix>.gapfilled.fa` | Optional FASTA written only with `--apply`, containing one record per assigned reference plus unassigned records. |
| `<prefix>.submission_checklist.tsv` | Submission-oriented checklist with FASTA/AGP consistency checks, unresolved gap counts, graph-filled span totals, non-ACGTN counts, and files to keep together. In planning mode, FASTA-dependent checks are marked as not checked because no final FASTA is written. |
| `<prefix>.run_summary.txt` | Inputs, parameters, output paths, and fill-status counts. |

### Example `chromo gapfill` Output

**Table 1. Example `gapfill_plan.tsv` row.** Selected columns from a graph fixture
show a junction resolved by reference-placement PAF support. The default
`accept_fill=no` makes planning review explicit before strict reviewed
application.

| scaffold | left_contig | right_contig | graph_status | path_nodes | candidate_paths | ref_path_support | ref_best_alt_support | risk_flags | fill_status | fill_bp | right_trim_bp | accept_fill | applied |
| --- | --- | --- | --- | --- | ---: | ---: | ---: | --- | --- | ---: | ---: | --- | --- |
| `chr1` | `chr1_left` | `chr1_right` | `ref_paf_resolved_paths` | `left+,bridge_good+,right+` | `2` | `8` | `6` | `branching,high_degree` | `fillable` | `4` | `4` | `no` | `no` |

**Listing 1. Example applied gapfilled FASTA output.** With
`--apply --reviewed-plan`, accepted fillable paths insert graph sequence and
trim the right flank by the terminal GFA overlap; unresolved or unaccepted
junctions use fallback N gaps.

```text
>chr1 contigs=2 filled_gaps=1 fallback_gaps=0 fill_bp=4 fallback_gap_bp=0 trimmed_bp=4
AAAACCCCGGGGTTTT
```

## `chromo gapfill` Parameters

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `--ordered-fasta` | input mode | Final ordered FASTA from `chromo sort`; use with `--assignments`. |
| `--assignments` | input mode | Matching `<prefix>.contig_assignments.tsv` report from `chromo sort`; use with `--ordered-fasta`. |
| `--scaffold-fasta` | input mode | Existing scaffold FASTA; use with `--agp` for post-scaffold gap filling. |
| `--agp` | input mode | AGP map for `--scaffold-fasta`, used to recover component identities and N-gap spans. |
| `--gfa` | required | Assembly graph GFA containing segment sequences and links. |
| `--project-gfa-paths` | off | Use matching GFA `P`/`W` path records to project components to terminal unitigs before path search. Sequence-bearing projected paths can be fillable; `.noseq.gfa` paths remain topology-only. |
| `--projection-trim-overlaps` | off | With `--project-gfa-paths`, subtract path-step overlaps while building the path projection. |
| `--gaf` | none | Optional GAF graph alignments used to resolve otherwise ambiguous candidate graph paths. |
| `--hic-pairs` | none | Optional TSV of graph-node contact counts with `node_a`, `node_b`, and `count` columns. |
| `--ref-paf` | none | Optional reference-to-assembly PAF used to score intermediate graph nodes against the expected reference-space gap. |
| `--read-paf` | none | Optional long-read-to-contig, ordered-FASTA, or AGP-component PAF used as bridge evidence across adjacent component ends; it does not provide inserted sequence. |
| `--patch-table` | none | Optional TSV of external patch candidates keyed by `scaffold`, `left_contig`, and `right_contig`; used as review evidence. |
| `--patch-fasta` | none | Optional FASTA containing `patch_id` records referenced by `--patch-table`. |
| `--output-prefix` | required | Prefix for gapfill plan, run summary, and optional gapfilled FASTA. |
| `--apply` | off | Write `<prefix>.gapfilled.fa`; requires `--reviewed-plan` or `--apply-all-fillable`. |
| `--apply-all-fillable` | off | With `--apply`, apply every currently fillable path without a reviewed plan. Use deliberately for exploratory or benchmark runs. |
| `--reviewed-plan` | none | Optional edited gapfill plan TSV or `chromo eval gapfill` review-event TSV. With `--apply`, only accepted rows are applied after the current path is rechecked. |
| `--review-html` | none | Optional self-contained HTML review dashboard for the gapfill plan. |
| `--fixed-gap-bp` | none | Use this many Ns for unresolved gaps in `--apply` output instead of inferred reference-space gaps. |
| `--max-path-edges` | `4` | Maximum GFA link depth searched between adjacent sorted contigs. |
| `--max-candidate-paths` | `2` | Stop path enumeration after this many candidates. The default distinguishes unique from ambiguous paths. |
| `--min-gaf-mapq` | `20` | Minimum GAF MAPQ for a read path to support a candidate fill. |
| `--min-gaf-path-support` | `1` | Minimum supporting GAF read paths required to resolve an ambiguous branch. |
| `--min-hic-path-support` | `1` | Minimum summed Hi-C contact support required to resolve an ambiguous branch. |
| `--min-ref-path-support` | `1` | Minimum expected-gap reference-placement support required to resolve an ambiguous branch. |
| `--min-ref-paf-mapq` | `0` | Minimum MAPQ for PAF rows used by `--ref-paf`. |
| `--min-ref-paf-idy` | `0.0` | Minimum percent identity for PAF rows used by `--ref-paf`. |
| `--include-secondary-ref-paf` | off | Include secondary PAF rows marked `tp:A:S` when reading `--ref-paf`. |
| `--min-read-mapq` | `0` | Minimum MAPQ for rows from `--read-paf`. |
| `--read-terminal-window-bp` | `5000` | Terminal window used when detecting long-read bridges between adjacent contigs. |
| `--read-min-anchor-bp` | `500` | Minimum aligned bases required at each contig end for a read to count as a bridge. |
| `--max-fill-bp` | `1000000` | Maximum inserted graph sequence allowed for one fill. Set negative to disable. |
| `--include-fill-sequences` | off | Include candidate fill sequences in the TSV plan. |
| `--simple-headers` | off | Write gapfilled FASTA headers containing only the scaffold ID. |

## Reasoning Behind `chromo gapfill`

### Filling Is Explicit

`chromo scaffold --gfa` remains report-only. `chromo gapfill` is the explicit
sequence-changing command, and it only changes sequence when `--apply` is paired
with either `--reviewed-plan` or `--apply-all-fillable`.
This keeps evidence review separate from FASTA construction.

### Reviewed Plan Gate

For strict reviewed application, run either `chromo eval gapfill` or
`chromo gapfill` once in planning mode, mark only approved rows as accepted,
then rerun with `--apply` and `--reviewed-plan`. ChromoSort recomputes the graph
path and validates the accepted row before applying it. Accepted rows whose
current `path_nodes` or fillable status no longer match are rejected with an
error. `--review-html` writes a browser-based table for the legacy plan-review
step; `chromo eval gapfill` writes the table-only counterpart.

For deliberate exploratory or benchmark runs, `--apply --apply-all-fillable`
applies all currently fillable paths. That mode still refuses ambiguous,
unsequenced, oversized, stale, or flank-mismatched paths, but it does not require
human acceptance.

### Unique Paths Or Unique Evidence

Assembly graphs often contain repeats, bubbles, and alternate paths. If more
than one candidate path is found within the search limit, `chromo gapfill`
usually marks the junction `ambiguous_paths` and falls back to Ns in applied
output. GAF read paths, Hi-C contacts, and reference-placement PAF are supported
tie-breakers: an ambiguous branch can be resolved only when one candidate has
unique support above the configured threshold and no other evidence source
uniquely supports a different path. Ties, weak support, or conflicting evidence
remain unresolved for manual review.

### Verify the Flanks

Before applying a graph fill, ChromoSort checks that the ordered FASTA flank
sequences match the oriented GFA segments used by the path. This protects
against applying a graph path to a FASTA that has been renamed, trimmed,
reverse-complemented, or otherwise edited without matching graph coordinates.

### Use The Right Graph Level

When GFA segment names directly match the ordered FASTA contigs, use `--gfa` as
direct contig-graph mode. This is the current applied gapfill mode.

For hifiasm-style unitig graphs, the FASTA may be contig-level while the GFA is
unitig-level. In that case, projection through GFA `P` path or `W` walk records
is required to map contig ends back to terminal unitigs. Use `--project-gfa-paths`
for that mode, and use `chromo graph-map` when you want a separate projection
QC table first. A sequence-bearing unitig GFA can support applied filling when
terminal unitig sequence validates against the component ends. A `.noseq.gfa`
with `LN:i` lengths can support projection and topology review, but it cannot
support applied sequence filling.
If no GFA is available, `chromo gapfill` can report long-read bridge evidence
only after future non-graph patch modes are added; current graph-aware filling
cannot insert sequence without a graph path or accepted sequence source.
External patch candidates imported with `--patch-table` are comparison evidence,
not an accepted sequence source.
