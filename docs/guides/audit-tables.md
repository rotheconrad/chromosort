---
title: Reading ChromoSort Audit Tables
description: Learn how to interpret ChromoSort TSV reports, status labels, review tables, and sequence-changing output provenance.
---

# Reading ChromoSort Audit Tables

Use this guide when you have a ChromoSort TSV output and need to decide what it
means before trusting a FASTA, plot, manual recipe, scaffold, or graph fill.

The important habit is simple: read the table before reading the FASTA header.
ChromoSort FASTA headers carry useful provenance, but the TSV reports are the
decision log.

## The Core Idea

Most ChromoSort commands write both sequence outputs and audit tables. The
sequence file is for downstream tools. The table explains why each record was
kept, rejected, split, cut, joined, reviewed, or filled.

Audit tables are deliberately redundant. They keep original IDs, new IDs,
reference labels, coordinates, status fields, and evidence summaries together
so that users can review decisions without reverse-engineering them from FASTA
headers.

## How To Read Any Audit Table

Start with these columns or column families:

| What to find | Why it matters |
| --- | --- |
| Original ID | Tells you which input contig, piece, scaffold junction, or graph path the row describes. |
| New ID | Tells you what name was written to the output FASTA, if any. |
| Status or action | The command's decision: kept, discarded, split, not split, gap mode, fill status, and so on. |
| Accepted or applied field | In reviewed workflows, tells you whether a row can change sequence. |
| Reference and coordinate fields | Let you compare the table to dot plots and alignment evidence. |
| Evidence summaries | Coverage, identity, MAPQ-derived filtering, graph support, GAF support, long-read support, and risk flags. |
| Reason fields | Explain why a candidate was rejected, smoothed over, left unresolved, or marked stale. |

Then ask two questions:

1. Did this row change sequence or only report evidence?
2. If it changed sequence, which input FASTA and evidence files made that change
   valid?

## Output Families

| Command | Tables to read first | Main question |
| --- | --- | --- |
| `chromo sort` | `<prefix>.contig_assignments.tsv`, `<prefix>.contig_ref_matches.tsv`, `<prefix>.chromosome_summary.tsv` | Which contigs were kept, assigned, filtered, or flagged for split review? |
| `chromo clean` | `<prefix>.initial_sort.contig_assignments.tsv`, `<prefix>.fix_report.tsv`, `<prefix>.clean_contigs.tsv` | Which raw contigs were discarded, inspected, split, or emitted into the cleaned FASTA? |
| `chromo eval` | `<prefix>.fix_review.tsv`, `<prefix>.scaffold_review.tsv`, `<prefix>.gapfill_review.tsv` | Which candidate decisions need human accept/reject review? |
| `chromo fix` | The `--report` TSV | Which requested or detected contigs were split, copied, or left unsplit? |
| `chromo cut` | The `--report` TSV | Which exact requested cut positions produced which pieces? |
| `chromo manual apply` | The `--report` TSV | Which browser-reviewed pieces were emitted or removed? |
| `chromo gafprep` | `<prefix>.targets.tsv`, `<prefix>.selected_reads.tsv`, `<prefix>.selected_read_review_links.tsv`, `<prefix>.dropped_gfa_links.tsv` | Which reads were selected for targeted GraphAligner, which review rows selected them, and did GFA sanitization limit evidence? |
| `chromo graph-map` | `<prefix>.utg_to_ctg.tsv`, `<prefix>.path_summary.tsv`, `<prefix>.warnings.tsv` | Did unitig graph coordinates project cleanly onto contig FASTA coordinates? |
| `chromo scaffold` | `<prefix>.scaffold_gaps.tsv`, `<prefix>.scaffold_summary.tsv`, `<prefix>.submission_checklist.tsv`, optional `<prefix>.graph_gaps.tsv` | What gaps, overlaps, trims, graph context, and FASTA/AGP handoff checks were recorded? |
| `chromo gapfill` | `<prefix>.gapfill_plan.tsv`, `<prefix>.submission_checklist.tsv` | Which graph paths are fillable, ambiguous, risky, accepted, or applied, and what final handoff checks remain? |

## Status Gallery

### Sort Assignment Rows

The `status` and `kept` fields in `contig_assignments.tsv` tell you whether a
contig entered the ordered FASTA.

| Status | Meaning | Usual next question |
| --- | --- | --- |
| `kept` | Passed placement and overlap filters. | Does the dot plot support this order and orientation? |
| `kept_split_candidate` | Retained and flagged as a strong multi-reference candidate. | Should this contig go through `chromo eval fix`, `chromo manual`, or `chromo fix` review? |
| `kept_large_alignment` | Rescued because the best reference match was very large despite slightly low query coverage. | Is the missing coverage due to fragmentation, repeats, or a real issue? |
| `kept_terminal_overlap` | Retained because it contributes enough one-sided terminal reference span. | Should scaffolding later report or trim the overlap? |
| `no_alignment` | No usable alignment rows were found. | Do names match, and was the aligner too strict? |
| `below_min_aligned_bp` or `below_min_query_cov` | Alignment support did not pass thresholds. | Should thresholds change, or is the contig truly weakly supported? |
| `ambiguous_ref_match` | No reference dominated enough to assign confidently. | Is this repeat signal, a real translocation, or a split candidate? |
| `duplicate_overlap` | A better contig already covers nearly all of the reference span. | Is this an alternate fragment, haplotig, repeat, or real duplicated sequence? |
| `terminal_overlap` | A one-sided overlap did not pass keep or rescue thresholds. | Does the extension matter biologically or for scaffolding? |

### Fix Report Rows

`chromo fix` reports both split and not-split outcomes. A `split` row describes
one emitted piece. A `not_split_*` row records why a reviewed or detected contig
was copied unchanged.

Pay attention to:

- `original_contig`,
- `status`,
- `new_contig`,
- `slice_start` and `slice_end`,
- `dominant_ref`,
- `orientation`,
- reason fields such as smoothing or breakpoint-budget rejection.

If a candidate is `not_split_smooth`, the planner saw discordance but decided it
was not worth a breakpoint under the current mode and thresholds. That is often
a prompt for manual review rather than a bug.

### Eval Review Tables

`chromo eval` writes shared review-event tables. These are not sequence outputs.
They are editable decision queues.

Look for:

- `event_type`, such as `split_piece`, `scaffold_gap`, or `fill_path`,
- `accept` or task-specific accepted fields,
- source contigs and coordinates,
- evidence columns such as `graph_*`, `gaf_*`, and `longread_*`,
- notes or review fields added by the user.

Rejected, deleted, stale, or unaccepted rows should not change sequence. The
matching executor revalidates accepted rows before applying them.

### Scaffold Gap Reports

`scaffold_gaps.tsv` explains every join between adjacent sorted contigs.

Key fields include:

- `raw_inferred_gap_bp`,
- `gap_bp`,
- `gap_mode`,
- `overlap_class`,
- `overlap_action`,
- trimming and sequence-identity fields when overlap policies check sequence.

A negative inferred gap means adjacent reference spans overlap. By default,
ChromoSort writes a zero-length FASTA gap and reports the overlap. Trimming
happens only when an explicit overlap policy asks for it.

### Gapfill Plans

`gapfill_plan.tsv` is a review table first and a sequence application log when
`--apply` is used.

Read these fields together:

- `graph_status`,
- `path_nodes`,
- `candidate_paths`,
- GAF, Hi-C, and reference-placement support columns,
- `risk_flags`,
- `fill_status`,
- `accept_fill`,
- `applied`.

Without a reviewed plan, `chromo gapfill --apply` applies fillable paths. With a
reviewed plan, only accepted rows are applied and all accepted rows are
rechecked against the current graph path and fillability status.

## Practical Review Workflow

1. Open the run summary to confirm inputs and thresholds.
2. Open the main row-level table for the command.
3. Sort or filter by status fields.
4. Inspect rows that changed sequence.
5. Inspect rows that were rejected or left unresolved.
6. Compare suspicious rows to dot plots or manual dashboard evidence.
7. Keep the table beside the FASTA in downstream folders.

For spreadsheet review, freeze the identifier and status columns before editing
accept/reject fields. Avoid changing provenance columns unless the command docs
explicitly say a field is user-editable.

## Cheat Sheet

| If you want to know... | Read... |
| --- | --- |
| Why a contig was kept or discarded | `contig_assignments.tsv` |
| Which reference each contig matched before final assignment | `contig_ref_matches.tsv` |
| Which raw contigs were emitted by `chromo clean` | `clean_contigs.tsv` |
| Which fix pieces replaced a source contig | `fix` report or `fix_report.tsv` |
| Whether a reviewed eval row can change sequence | `accept` plus `event_type` |
| How many Ns were inserted between scaffold contigs | `scaffold_gaps.tsv` |
| Whether an overlap was trimmed | `overlap_action` in `scaffold_gaps.tsv` |
| Whether graph context changed a scaffold | It does not by default; read `graph_gaps.tsv` as report-only evidence. |
| Why a graph fill did not apply | `fill_status`, `risk_flags`, and `applied` in `gapfill_plan.tsv` |

## Common Traps

Do not parse FASTA headers when a TSV report exists. Headers are helpful, but
tables are more complete and easier to audit.

Do not treat report-only graph evidence as a sequence change. `graph_*` columns
often explain context without changing the FASTA.

Do not assume `accept_fill=no` means the candidate is biologically false. It
means the planning table has not accepted that sequence-changing action.

Do not edit provenance fields in reviewed tables unless you are intentionally
creating a new reviewed decision row and understand the executor validation
rules.

Do not apply an old reviewed table after changing FASTA, assignments, graph
inputs, or path-search settings. Regenerate the table from current inputs.

## What To Look At Next In ChromoSort

- Use [Output Files]({{ '/outputs/' | relative_url }}) for command-by-command
  output definitions.
- Use [FASTA And Evidence Name Matching]({{ '/guides/name-matching/' | relative_url }})
  when a report contains missing or stale identifiers.
- Use [How to Interpret Dot Plots]({{ '/dot-plots/' | relative_url }}) when a
  table status needs visual context.
- Use [chromo eval]({{ '/commands/eval/' | relative_url }}) and
  [chromo manual]({{ '/commands/manual/' | relative_url }}) when a table row
  needs human review before sequence changes.
