---
title: chromo scaffold
description: Usage, outputs, parameters, and reasoning for chromo scaffold.
---

# chromo scaffold

Use `chromo scaffold` when the final sorted and filtered contigs look good and
you want one FASTA record per reference chromosome or linkage group.

## What `chromo scaffold` Does

Given a final `chromo sort` ordered FASTA and its matching
`<prefix>.contig_assignments.tsv`, `chromo scaffold`:

1. Reads kept contigs from the assignment report.
2. Reads the ordered FASTA records by their `new_name` values.
3. Groups contigs by assigned reference sequence in ordered FASTA order.
4. Joins neighboring contigs with N gaps.
5. Infers gap length from adjacent reference coordinates by default.
6. Reports negative inferred gaps as adjacent reference overlaps.
7. Optionally trims reviewed terminal overlaps according to `--overlap-policy`.
8. Optionally uses a fixed user-provided number of Ns between every neighboring
   contig.
9. Optionally writes a report-only GFA graph-evidence table for adjacent
   scaffold junctions.
10. Optionally applies accepted gap overrides from `chromo eval scaffold`.
11. Writes scaffold FASTA, AGP, component provenance, gap report, scaffold
    summary, and run summary files.

The intended input is the final ordered FASTA from the same `chromo sort` run as
the assignment report. If you run `chromo fix`, re-run `chromo sort` on the
fixed assembly before scaffolding so the coordinates and FASTA names match the
final contigs.

## Run `chromo scaffold` With Inferred Gaps

```bash
chromo scaffold \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --output-prefix results/sample
```

For adjacent contigs on the same reference, inferred gaps are calculated as:

```text
next_ref_start - previous_ref_end - 1
```

Negative values, which indicate overlapping reference spans, are written as
zero-length gaps in the FASTA and reported in the gap TSV by default. Use
`--overlap-policy warn` to keep the same FASTA behavior while emitting stderr
warnings, `--overlap-policy trim-reference` to trim the right contig by the
reference-inferred terminal overlap, or `--overlap-policy trim-sequence` to trim
only when the left suffix and right prefix confirm the overlap sequence at
`--trim-sequence-min-identity`.

## Run `chromo scaffold` With Fixed Gaps

```bash
chromo scaffold \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --output-prefix results/sample.fixed100 \
  --fixed-gap-bp 100
```

Fixed-gap mode ignores inferred gap length for FASTA construction and inserts
the requested number of Ns between every neighboring contig on the same
scaffold. The report still records the raw inferred gap for comparison.

## Run `chromo scaffold` With A Reviewed Plan

```bash
chromo eval scaffold \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --output-prefix results/sample.eval_scaffold \
  --gfa assembly_graph.gfa \
  --gaf reads_to_graph.gaf \
  --read-paf reads_to_assembly.paf

chromo scaffold \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --reviewed-plan results/sample.eval_scaffold.scaffold_review.tsv \
  --output-prefix results/sample.reviewed_scaffold
```

The reviewed table is optional. When supplied, accepted `scaffold_gap` rows
override `gap_bp` for matching `scaffold`/`left_contig`/`right_contig`
junctions, and the gap report marks those junctions with `gap_mode=reviewed`.
Junctions without accepted rows keep the normal inferred or fixed-gap behavior.
Accepted rows that no longer match the current ordered FASTA and assignment TSV
are rejected as stale.

GFA, GAF, and long-read PAF evidence in the eval table is review context for
the junction. `chromo scaffold --reviewed-plan` applies accepted gap lengths; it
does not reorder, orient, trim, or gapfill contigs based on GAF or long-read
support alone.

## Run `chromo scaffold` With Overlap Trimming

```bash
chromo scaffold \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --output-prefix results/sample.trim_seq \
  --overlap-policy trim-sequence
```

`trim-reference` removes the reference-inferred overlap from the left side of
the right contig when the adjacent reference spans form a terminal overlap.
`trim-sequence` is more conservative: it trims only when the left contig suffix
and right contig prefix match across the inferred overlap with at least
`--trim-sequence-min-identity` identity. Non-terminal overlaps are reported but
not trimmed by either trimming policy.

## Run `chromo scaffold` With GFA Graph Evidence

```bash
chromo scaffold \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --output-prefix results/sample.graph \
  --gfa assembly_graph.gfa
```

When `--gfa` is provided, `chromo scaffold` writes
`<prefix>.graph_gaps.tsv`. With the default `--graph-overlap-policy report`,
this is report-only: scaffold FASTA construction, gap lengths, and overlap
policies are unchanged. The graph report resolves each adjacent scaffold contig
to a GFA segment when possible, records direct orientation-aware GFA links, and
searches for a short explicit GFA path up to `--graph-max-path-edges`.

Set `--graph-overlap-policy warn` to emit warnings for graph-confirmed terminal
overlaps. Set `--graph-overlap-policy confirm` only after review; it allows a
direct orientation-matching GFA edge to trim a terminal reference-space overlap
when the normal overlap policy is `zero-gap` or `warn`. Nonterminal overlaps,
missing nodes, orientation mismatches, and indirect paths remain report-only.

## `chromo scaffold` Outputs

| Output | Description |
| --- | --- |
| `<prefix>.scaffold.fa` | One FASTA record per assigned reference sequence, with ordered contigs joined by Ns. |
| `<prefix>.scaffold.agp` | AGP 2.1 scaffold map with ordered contig components and written N gaps. |
| `<prefix>.scaffold_components.tsv` | Human-readable provenance table mirroring the AGP rows, including source, status, component coordinates, gap lengths, linkage evidence, and notes. |
| `<prefix>.scaffold_gaps.tsv` | One row per inserted gap with flanking contigs, inferred gap, written gap, overlap bp/class/fractions, overlap policy/action, trimmed bp, and sequence-overlap identity when checked. |
| `<prefix>.graph_gaps.tsv` | Optional report-only GFA evidence for adjacent scaffold junctions when `--gfa` is provided, including direct links, short paths, orientations, overlap bp, and missing-node statuses. |
| `<prefix>.scaffold_summary.tsv` | One row per scaffold with contig count, scaffold length, sequence bp, gap bp, overlap totals, trimming totals, and ordered contig list. |
| `<prefix>.submission_checklist.tsv` | Submission-oriented checklist with FASTA character checks, AGP continuity checks, gap/component counts, unresolved gap counts, and files to keep together. |
| `<prefix>.run_summary.txt` | Inputs, gap model, output paths, and total scaffold counts. |

### Example `chromo scaffold` Output

**Table 1. Example `scaffold_gaps.tsv` row.** The gap report records the flanking
contigs, inferred reference-space gap, FASTA gap actually written, overlap
classification, and overlap policy action.

| scaffold | left_contig | right_contig | raw_inferred_gap_bp | gap_bp | gap_mode | overlap_class | overlap_action |
| --- | --- | --- | ---: | ---: | --- | --- | --- |
| `chr1` | `chr1_contigA` | `chr1_contigB` | `5` | `5` | `inferred` | `no_overlap` | `none` |

**Table 2. Example `scaffold_summary.tsv` rows.** The summary table gives one
row per emitted scaffold record.

| scaffold | contigs | scaffold_bp | sequence_bp | gap_bp | gaps | ordered_contigs |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| `chr1` | `2` | `12` | `7` | `5` | `1` | `chr1_contigA,chr1_contigB` |
| `chr2` | `1` | `2` | `2` | `0` | `0` | `chr2_contigC` |

**Listing 1. Example scaffold FASTA output.** Scaffold headers summarize the
number of source contigs, sequence bases, gap bases, and gap mode.

```text
>chr1 contigs=2 sequence_bp=7 gap_bp=5 gap_mode=inferred
AAAANNNNNTTT
>chr2 contigs=1 sequence_bp=2 gap_bp=0 gap_mode=inferred
GG
```

## `chromo scaffold` Parameters

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `--ordered-fasta` | required | Final ordered FASTA from `chromo sort`. |
| `--assignments` | required | Matching `<prefix>.contig_assignments.tsv` report from `chromo sort`. |
| `--output-prefix` | required | Prefix for scaffold FASTA and reports. |
| `--fixed-gap-bp` | none | Insert this many Ns between neighboring contigs instead of inferred gaps. |
| `--reviewed-plan` | none | Optional edited table from `chromo eval scaffold`; accepted `scaffold_gap` rows override matching junction gap lengths while evidence columns remain provenance. |
| `--overlap-policy` | `zero-gap` | Handling for negative inferred gaps: `zero-gap`, `warn`, `trim-reference`, or `trim-sequence`. |
| `--trim-sequence-min-identity` | `0.98` | Minimum suffix/prefix identity required by `--overlap-policy trim-sequence`. |
| `--simple-headers` | off | Write scaffold FASTA headers containing only the scaffold ID. |
| `--gfa` | none | Optional assembly graph GFA for report-only graph evidence at scaffold junctions. |
| `--graph-overlap-policy` | `report` | Graph safety mode for negative-gap overlaps: `report`, `warn`, or `confirm`. `confirm` only trims direct oriented GFA-confirmed terminal overlaps under `zero-gap`/`warn`. |
| `--graph-max-path-edges` | `4` | Maximum explicit GFA link depth searched for short paths in the graph gap report. |

## Reasoning Behind `chromo scaffold`

### Require the Assignment Report

The ordered FASTA contains sequence, but the assignment report contains the
reference coordinates needed to infer gap lengths. Requiring both files keeps
scaffolding explicit and prevents guessing chromosome names from FASTA IDs,
which can be fragile when contig names contain separators.

### Infer Gaps Conservatively

Inferred gaps are based only on adjacent retained contigs on the same reference
sequence. ChromoSort does not add leading or trailing Ns for chromosome ends,
and overlapping reference spans become zero-length FASTA gaps rather than
negative gaps by default. The raw inferred value, overlap bp, overlap class,
overlap fractions, policy, action, and trimming amount are reported so users can
inspect overlap or coordinate oddities.

### Trim Only When Asked

Overlaps are not trimmed by default because a negative reference gap can mean
several different things: a real dovetail, retained alternate sequence,
reference/assembly structural difference, or an alignment artifact. The
`trim-reference` policy trims only terminal overlaps and trusts the reference
coordinate estimate. The `trim-sequence` policy trims only terminal overlaps
whose left suffix and right prefix agree at high identity. Contained/internal
overlaps are reported for review rather than trimmed automatically.

### Keep Fixed Gaps Available

Some downstream tools and submission workflows prefer a constant gap size. The
`--fixed-gap-bp` option supports that convention while preserving the inferred
gap estimate in the report for transparency.

### Keep Submission Provenance

`chromo scaffold` writes AGP and component-provenance tables beside the FASTA.
The AGP maps each emitted scaffold record back to ordered contig components and
written N gaps. When overlap trimming is enabled, the component coordinates
reflect the trimmed component slice. Keep the AGP, component TSV, gap report,
and input ordered FASTA together when preparing NCBI-style submission material
or handing the scaffold to another reviewer.

### Keep Graph Evidence Report-Only

The optional GFA report is evidence, not an assembler. It can show that two
adjacent scaffold contigs are directly linked in the assembly graph, connected
through a short path, connected only in a different orientation, absent from the
graph, or disconnected within the configured search depth. It does not fill
gaps, trim sequence, or reorder contigs. Those operations remain explicit
review steps.
