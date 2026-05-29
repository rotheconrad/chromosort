---
title: chromo gapfill
description: Usage, outputs, parameters, and reasoning for chromo gapfill.
---

# chromo gapfill

Use `chromo gapfill` after final sorting and manual review when a GFA graph gives a
simple sequence path between adjacent sorted contigs.

## What `chromo gapfill` Does

Given a final `chromo sort` ordered FASTA, the matching assignment report, a
GFA assembly graph, and optional GAF graph alignments, Hi-C pair evidence, or
reference-placement PAF evidence, `chromo gapfill`:

1. Groups retained contigs by assigned reference sequence.
2. Looks at adjacent contig pairs in sorted order.
3. Resolves each flank to a GFA segment using original and renamed contig IDs.
4. Enumerates graph paths up to `--max-path-edges`.
5. Uses GAF read-path support, Hi-C contact support, and reference-placement
   PAF support to resolve an otherwise ambiguous graph branch only when one
   candidate path has unique support above threshold and evidence sources do
   not conflict.
6. Annotates candidate-path risk, including high-degree graph nodes, self-loop
   nodes, unsequenced nodes, cycle guards, weak/tied/conflicting support, and a
   branch-complexity score.
7. Rejects missing nodes, disconnected flanks, unresolved ambiguous paths,
   unsequenced nodes, unknown or invalid overlaps, oversized fills, and flank
   sequence mismatches.
8. Writes `<prefix>.gapfill_plan.tsv` for review with `accept_fill=no` by default,
   can write a self-contained HTML reviewer with `--review-html`, and can accept
   the shared review-event table from `chromo eval gapfill`.
9. With `--apply`, writes `<prefix>.gapfilled.fa`. Without `--reviewed-plan`, all
   currently fillable paths are applied; with `--reviewed-plan`, only rows with
   accepted fill decisions are applied and other junctions fall back to inferred
   or fixed N gaps.

## Plan Graph Fills

```bash
chromo gapfill \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --gfa assembly_graph.gfa \
  --ref-paf paf/sample.ref_vs_asm.paf \
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

## Apply Reviewed Graph Fills

```bash
chromo gapfill \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --gfa assembly_graph.gfa \
  --ref-paf paf/sample.ref_vs_asm.paf \
  --gaf reads_to_graph.gaf \
  --hic-pairs graph_contacts.tsv \
  --reviewed-plan results/sample.eval_gapfill.gapfill_review.tsv \
  --output-prefix results/sample.reviewed_gapfill \
  --apply
```

Applied mode still refuses ambiguous or unverifiable paths. If GAF is provided,
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
reference-space N gap, or `--fixed-gap-bp` when provided.

## `chromo gapfill` Outputs

| Output | Description |
| --- | --- |
| `<prefix>.gapfill_plan.tsv` | One row per adjacent sorted contig pair with graph status, path nodes, GAF support counts/status/supporting reads, Hi-C and reference-placement support counts, risk flags, branch-complexity score, high-degree/self-loop/unsequenced node lists, fill status, inserted bp, right-trim bp, fallback gap bp, editable `accept_fill`, and whether the fill was applied. |
| `--review-html` path | Optional self-contained HTML table for reviewing gapfill-plan rows, comparing candidate paths, and exporting a reviewed-plan TSV. |
| `<prefix>.gapfilled.fa` | Optional FASTA written only with `--apply`, containing one record per assigned reference plus unassigned records. |
| `<prefix>.run_summary.txt` | Inputs, parameters, output paths, and fill-status counts. |

### Example `chromo gapfill` Output

**Table 1. Example `gapfill_plan.tsv` row.** Selected columns from a graph fixture
show a junction resolved by reference-placement PAF support. The default
`accept_fill=no` makes planning review explicit before strict reviewed
application.

| scaffold | left_contig | right_contig | graph_status | path_nodes | candidate_paths | ref_path_support | ref_best_alt_support | risk_flags | fill_status | fill_bp | right_trim_bp | accept_fill | applied |
| --- | --- | --- | --- | --- | ---: | ---: | ---: | --- | --- | ---: | ---: | --- | --- |
| `chr1` | `chr1_left` | `chr1_right` | `ref_paf_resolved_paths` | `left+,bridge_good+,right+` | `2` | `8` | `6` | `branching,high_degree` | `fillable` | `4` | `4` | `no` | `no` |

**Listing 1. Example applied gapfilled FASTA output.** With `--apply`, fillable
paths insert graph sequence and trim the right flank by the terminal GFA
overlap; unresolved junctions use fallback N gaps.

```text
>chr1 contigs=2 filled_gaps=1 fallback_gaps=0 fill_bp=4 fallback_gap_bp=0 trimmed_bp=4
AAAACCCCGGGGTTTT
```

## `chromo gapfill` Parameters

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `--ordered-fasta` | required | Final ordered FASTA from `chromo sort`. |
| `--assignments` | required | Matching `<prefix>.contig_assignments.tsv` report from `chromo sort`. |
| `--gfa` | required | Assembly graph GFA containing segment sequences and links. |
| `--gaf` | none | Optional GAF graph alignments used to resolve otherwise ambiguous candidate graph paths. |
| `--hic-pairs` | none | Optional TSV of graph-node contact counts with `node_a`, `node_b`, and `count` columns. |
| `--ref-paf` | none | Optional reference-to-assembly PAF used to score intermediate graph nodes against the expected reference-space gap. |
| `--output-prefix` | required | Prefix for gapfill plan, run summary, and optional gapfilled FASTA. |
| `--apply` | off | Write `<prefix>.gapfilled.fa` using only accepted graph paths. |
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
| `--max-fill-bp` | `1000000` | Maximum inserted graph sequence allowed for one fill. Set negative to disable. |
| `--include-fill-sequences` | off | Include candidate fill sequences in the TSV plan. |
| `--simple-headers` | off | Write gapfilled FASTA headers containing only the scaffold ID. |

## Reasoning Behind `chromo gapfill`

### Filling Is Explicit

`chromo scaffold --gfa` remains report-only. `chromo gapfill` is the explicit
sequence-changing command, and it only changes sequence when `--apply` is set.
This keeps evidence review separate from FASTA construction.

### Reviewed Plan Gate

For strict reviewed application, run either `chromo eval gapfill` or
`chromo gapfill` once in planning mode, mark only approved rows as accepted,
then rerun with `--apply` and `--reviewed-plan`. ChromoSort recomputes the graph
path and validates the accepted row before applying it. Accepted rows whose
current `path_nodes` or fillable status no longer match are rejected with an
error. `--review-html` writes a browser-based table for the legacy plan-review
step; `chromo eval gapfill` writes the table-only counterpart.

### Unique Paths Or Unique Evidence

Assembly graphs often contain repeats, bubbles, and alternate paths. If more
than one candidate path is found within the search limit, `chromo gapfill` usually
marks the junction `ambiguous_paths` and falls back to Ns in applied output. GAF
read paths and Hi-C contacts are supported tie-breakers: an ambiguous branch can
be resolved only when one candidate has unique support above the configured
threshold. Ties, weak support, or conflicting GAF/Hi-C support still remain
unresolved for manual review.

### Verify the Flanks

Before applying a graph fill, ChromoSort checks that the ordered FASTA flank
sequences match the oriented GFA segments used by the path. This protects
against applying a graph path to a FASTA that has been renamed, trimmed,
reverse-complemented, or otherwise edited without matching graph coordinates.
