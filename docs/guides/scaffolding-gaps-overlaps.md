---
title: Scaffolding, Gaps, And Overlaps
description: Review chromo scaffold gap inference, fixed gaps, negative gaps, overlap classes, trimming policies, graph evidence, and reviewed gap overrides.
---

# Scaffolding, Gaps, And Overlaps

Use this guide when final sorted contigs look ready to join into chromosome or
linkage-group records, but the gap and overlap report needs interpretation.

The key question is:

> Should this junction be joined with Ns, trimmed, reviewed, or left as a
> reported overlap?

## The Core Idea

`chromo scaffold` joins retained, sorted contigs into one FASTA record per
assigned reference sequence. It does not discover new contigs, reorder them, or
fill graph sequence. Its main job is to convert adjacent sorted contigs into
scaffold records while recording every inferred gap, overlap, trim, and
reviewed override.

The command requires the final ordered FASTA and the matching
`contig_assignments.tsv` from the same `chromo sort` run. The FASTA provides
sequence. The assignment report provides reference coordinates.

## What ChromoSort Reads Or Writes

`chromo scaffold` reads:

- `<prefix>.ordered.fa`,
- matching `<prefix>.contig_assignments.tsv`,
- optional reviewed scaffold table,
- optional GFA graph context for report-only junction evidence.

It writes:

- `<prefix>.scaffold.fa`,
- `<prefix>.scaffold_gaps.tsv`,
- optional `<prefix>.graph_gaps.tsv`,
- `<prefix>.scaffold_summary.tsv`,
- `<prefix>.run_summary.txt`.

After scaffolding, align the scaffold FASTA if scaffold-level validation is
needed.

## Gap Inference

For adjacent contigs on the same assigned reference, inferred gap length is:

```text
next_ref_start - previous_ref_end - 1
```

The result appears as `raw_inferred_gap_bp` in `scaffold_gaps.tsv`.

| Raw inferred value | Meaning | Default FASTA behavior |
| --- | --- | --- |
| Positive | Reference coordinates leave a gap between adjacent contigs. | Insert that many Ns unless fixed or reviewed gap mode overrides it. |
| Zero | Adjacent reference spans touch. | Insert zero Ns. |
| Negative | Reference spans overlap. | Insert zero Ns and report the overlap unless an explicit overlap policy trims. |

`--fixed-gap-bp` changes the FASTA gap length to a constant value, but the
report still keeps `raw_inferred_gap_bp` for comparison.

## Gap Modes

| `gap_mode` | How `gap_bp` was chosen |
| --- | --- |
| `inferred` | From adjacent reference coordinates, with negative values converted to zero. |
| `fixed` | From `--fixed-gap-bp`. |
| `reviewed` | From an accepted `scaffold_gap` row in a reviewed table. |

Reviewed gap rows override only matching junction gap lengths. They do not
change order, orientation, trimming, graph branches, or gapfill sequence.

## Overlap Classes

Negative inferred gaps are classified in the report:

| `overlap_class` | Pattern | Review stance |
| --- | --- | --- |
| `no_overlap` | No negative gap. | Ordinary gap or adjacent contigs. |
| `terminal_overlap` | The right contig starts inside the left span but extends beyond it. | Possible dovetail or terminal redundancy; eligible for trimming policies. |
| `contained_overlap` | The right contig lies inside the left reference span. | Usually review carefully; not automatically trimmed. |
| `spanning_overlap` | The right contig spans over the left reference span. | Review as a coordinate or structural oddity. |
| `internal_overlap` | Partial nonterminal overlap. | Report for review; not automatically trimmed. |

Only terminal overlaps are trim candidates. Contained, spanning, and internal
overlaps are recorded but trimming is skipped.

## Overlap Policies

| Policy | FASTA behavior | Report clues |
| --- | --- | --- |
| `zero-gap` | Convert negative inferred gaps to zero Ns. | `overlap_action=zero_gap`. |
| `warn` | Same FASTA behavior as zero-gap, with warnings. | `overlap_action=zero_gap` plus stderr warnings. |
| `trim-reference` | Trim the right contig by the reference-inferred terminal overlap. | `overlap_action=trimmed_reference`. |
| `trim-sequence` | Trim only if left suffix and right prefix match at the identity threshold. | `overlap_action=trimmed_sequence` or `trim_skipped_sequence_identity`. |

The `trim-sequence` policy is more conservative than `trim-reference` because
it asks the sequence to confirm the overlap before trimming. The report records
`trimmed_bp` and `sequence_overlap_identity` when sequence identity was checked.

## Graph Junction Evidence

With `--gfa`, `chromo scaffold` writes `<prefix>.graph_gaps.tsv`. By default,
this is report-only evidence:

- direct oriented links between adjacent contigs,
- short graph paths up to the configured depth,
- orientation mismatches,
- missing graph nodes,
- graph overlap details.

Graph evidence does not fill gaps. It can only affect trimming when
`--graph-overlap-policy confirm` is used with `zero-gap` or `warn`, a direct
orientation-matching graph edge exists, and the overlap is terminal. In that
case the scaffold gap report records
`overlap_action=graph_confirmed_trim_reference`.

Use `chromo gapfill` when the goal is to insert graph path sequence through an
N gap.

## Practical Scaffolding Workflow

Start from final sorted contigs:

```bash
chromo scaffold \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --output-prefix results/sample
```

Review:

1. `results/sample.scaffold_gaps.tsv` for every junction.
2. `raw_inferred_gap_bp`, `gap_bp`, and `gap_mode`.
3. `overlap_class`, `overlap_action`, `trimmed_bp`, and sequence identity.
4. `results/sample.scaffold_summary.tsv` for scaffold length and gap totals.
5. Optional `results/sample.graph_gaps.tsv` for report-only graph context.

If a few junctions need human gap-length overrides, create a review table:

```bash
chromo eval scaffold \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --gfa assembly_graph.gfa \
  --gaf reads_to_graph.gaf \
  --read-paf reads_to_assembly.paf \
  --output-prefix review/sample.scaffold
```

Then apply accepted rows:

```bash
chromo scaffold \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --reviewed-plan review/sample.scaffold.scaffold_review.tsv \
  --output-prefix results/sample.reviewed_scaffold
```

## Which Command Handles Which Junction Problem?

| Problem | Best command |
| --- | --- |
| Need one FASTA record per assigned reference | `chromo scaffold`. |
| Need constant N gaps for downstream convention | `chromo scaffold --fixed-gap-bp`. |
| Need to override a few gap lengths after review | `chromo eval scaffold`, then `chromo scaffold --reviewed-plan`. |
| Need to inspect graph links at adjacent contigs | `chromo scaffold --gfa` or `chromo eval scaffold --gfa`. |
| Need to insert graph sequence through a gap | `chromo gapfill --apply` after reviewed graph-fill planning. |
| Need to reorder or remove contigs | Go back to sort/manual review before scaffolding. |

## Cheat Sheet

| If you see... | Think... | Action |
| --- | --- | --- |
| `gap_mode=inferred` | Normal reference-coordinate gap model. | Check whether raw gap values look plausible. |
| `gap_mode=fixed` | Constant N gap mode. | Confirm this matches downstream requirements. |
| `gap_mode=reviewed` | Human-accepted override. | Keep the reviewed table with the output. |
| `raw_inferred_gap_bp < 0` | Adjacent reference spans overlap. | Read `overlap_class` and `overlap_action`. |
| `terminal_overlap` | Possible trim candidate. | Use explicit trimming policy only after review. |
| `trim_skipped_nonterminal` | Overlap was not terminal. | Review; do not expect automatic trimming. |
| `graph_gaps.tsv` direct edge | Graph supports adjacency. | Useful context; not a gap fill by itself. |

## Common Traps

Do not scaffold with an ordered FASTA and assignment table from different sort
runs.

Do not treat negative inferred gaps as automatic sequence errors. They can
reflect dovetails, alternate sequence, true structural differences, or
alignment artifacts.

Do not turn on trimming just to make the report quieter. Trimming is a
sequence-changing choice.

Do not expect `chromo scaffold --gfa` to insert graph sequence. It reports graph
context; `chromo gapfill` applies reviewed graph fills.

Do not use old contig-level PAF to validate scaffold records. Align the
scaffold FASTA if scaffold-level validation matters.

## What To Look At Next In ChromoSort

- Use [Spreadsheet Review Tables]({{ '/guides/spreadsheet-review-tables/' | relative_url }})
  for reviewed gap overrides.
- Use [Manual Dashboard Review]({{ '/guides/manual-dashboard-review/' | relative_url }})
  when junction rows need visual evidence.
- Use [Sorting Decisions And Duplicate-Overlap Filtering]({{ '/guides/sorting-decisions/' | relative_url }})
  before scaffolding questionable sort outputs.
- Use [chromo scaffold]({{ '/commands/scaffold/' | relative_url }}) and
  [chromo gapfill]({{ '/commands/gapfill/' | relative_url }}) for full command
  references.
