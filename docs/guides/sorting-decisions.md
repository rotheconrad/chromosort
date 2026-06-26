---
title: Sorting Decisions And Duplicate-Overlap Filtering
description: Learn how chromo sort assigns, keeps, filters, rescues, and flags contigs before writing an ordered FASTA.
---

# Sorting Decisions And Duplicate-Overlap Filtering

Use this guide when `chromo sort` kept some contigs, filtered others, or marked
split candidates, and you need to decide whether the result is trustworthy.

The main question is:

> Did ChromoSort remove noise, or did it hide something I should review?

## The Core Idea

`chromo sort` turns raw whole-genome alignment evidence into a conservative
ordered FASTA. It does this by assigning each contig to the reference sequence
that explains the most merged query span, then filtering contigs that are weak,
ambiguous, redundant, or mostly contained in a better-supported contig.

The command is not a chimeric-contig fixer. It does not cut a contig at an
alignment transition. When a contig has strong evidence for more than one
reference, `chromo sort` can keep it as `kept_split_candidate` so a human or a
review workflow can inspect it before `chromo fix`, `chromo cut`, or
`chromo manual` changes sequence.

## What ChromoSort Reads Or Writes

`chromo sort` reads:

- one assembly FASTA,
- one primary whole-genome alignment for that exact assembly FASTA, either
  MUMmer coords or minimap2 PAF,
- optional GFA graph context for report-only assignment context.

It writes:

- an ordered FASTA containing retained contigs,
- `<prefix>.contig_assignments.tsv`,
- `<prefix>.contig_ref_matches.tsv`,
- `<prefix>.chromosome_summary.tsv`,
- optional graph-aware report fields.

The ordered FASTA may be filtered, renamed, reordered, and optionally
reverse-complemented with `--orient-to-reference`. Re-align it before using it
as the assembly FASTA for another alignment-dependent decision.

## How Sort Decides

Start with these report fields:

| Field or concept | What it tells you |
| --- | --- |
| Total merged query coverage | How much of the contig is explained by usable alignment rows after overlapping query intervals are merged. |
| Best reference | The reference sequence with the strongest merged support for that contig. |
| `best_ref_share` | Whether one reference dominates the contig enough for confident assignment. |
| Per-reference matches | Which other reference sequences also explain substantial parts of the contig. |
| Novel reference span | How much new reference-space coverage the contig contributes after stronger contigs are considered. |
| Status | The final keep, filter, rescue, or review decision. |

Merged query coverage protects against double-counting overlapping alignments.
Best-reference share protects against assigning a contig when two references
are nearly tied. Duplicate-overlap filtering protects against writing a shorter
or weaker contig that mostly covers reference span already covered by a better
contig.

## Status Gallery

### `kept`

The contig passed the support, assignment, and overlap filters.

Next check: confirm that the dot plot and assignment report agree about order
and orientation. A `kept` status means the row is usable for sorting, not that
every local structural difference has been biologically interpreted.

### `kept_split_candidate`

The contig was retained, but it has substantial support for more than one
reference or an otherwise suspicious split-like signal.

Next check: inspect the dot plot and per-reference rows. Use
`chromo eval fix`, `chromo manual`, or targeted `chromo fix` only after the
candidate is reviewed.

### `kept_large_alignment`

The contig had a very large best-reference match and was rescued even though it
missed a standard coverage threshold.

Next check: ask why query coverage was low. Common causes include repeat-rich
sequence, fragmentation, aligner filtering, or a real assembly/reference
difference.

### `kept_terminal_overlap`

The contig overlapped a stronger contig but still contributed enough one-sided
terminal reference span to keep.

Next check: review it later during scaffolding. A terminal overlap can be a
useful extension, a dovetail, an alternate fragment, or an alignment artifact.

### `duplicate_overlap`

A stronger retained contig already covers nearly all of the reference span this
contig would contribute.

Next check: decide whether this is redundant assembly, an alternate haplotype,
a repeat-driven match, or real duplicated biology. `chromo sort` excludes the
contig from the ordered FASTA, but it does not delete the source FASTA.

### `terminal_overlap`

The contig had a one-sided overlap, but not enough new terminal span to trigger
the terminal-overlap rescue.

Next check: compare this row with nearby retained contigs and scaffold overlap
reports. Do not treat it as a within-contig breakpoint.

### `ambiguous_ref_match`

No reference sequence dominated enough to assign the contig confidently.

Next check: inspect per-reference evidence and dot plots. Ambiguity can come
from repeats, homeologous regions, contamination, a true translocation, or a
chimeric contig.

### `below_min_aligned_bp`, `below_min_query_cov`, Or `no_alignment`

The contig did not have enough usable alignment evidence.

Next check: verify FASTA names, aligner settings, PAF MAPQ filtering, coords
format, and whether the contig is real unplaced sequence.

## Decision Patterns

| Pattern | Usually means | Best next step |
| --- | --- | --- |
| One long clean hit to one reference | Straightforward sort candidate | Keep, order, and optionally orient with `--orient-to-reference`. |
| One whole-contig reverse-strand hit | Orientation difference | Use sort orientation if desired; do not split. |
| Two strong hits on different references | Split candidate or structural difference | Review with dot plots and `chromo eval fix`; avoid automatic `--all` fixes unless you have already validated the run. |
| Short off-target hits plus one dominant hit | Repeat, paralog, or aligner noise | Check thresholds and leave unsplit unless other evidence supports a real event. |
| Short contig contained in a longer retained contig | Redundant assembly or haplotig-like signal | Read `duplicate_overlap` rows before deciding whether the excluded contig is acceptable. |
| Adjacent retained contigs overlap in reference space | Scaffold or overlap issue | Review with `chromo scaffold` gap and overlap reports. |

## Practical Review Workflow

1. Open `<prefix>.contig_assignments.tsv`.
2. Filter by `status`.
3. Review all `kept_split_candidate`, `ambiguous_ref_match`,
   `duplicate_overlap`, and `terminal_overlap` rows.
4. Open `<prefix>.contig_ref_matches.tsv` for suspicious contigs.
5. Draw dot plots from the same raw FASTA and same alignment evidence.
6. Use `chromo eval fix` or `chromo manual` for candidates that might change
   sequence.
7. After any fix, cut, or manual FASTA is written, re-align that FASTA before
   sorting again.

Example review plot:

```bash
chromo plot \
  --ref-fasta reference.fa \
  --assembly-fasta raw.fa \
  --paf raw.paf \
  --assignments results/sample.contig_assignments.tsv \
  --output-prefix plots/sample.sort_review \
  --per-ref
```

This reviews the raw alignment rows in sorted order. It does not validate the
ordered FASTA as a new assembly.

## Cheat Sheet

| If you see... | Think... | Action |
| --- | --- | --- |
| `kept` | Confident enough for ordering | Inspect plots, then scaffold or re-align as needed. |
| `kept_split_candidate` | Keep now, review before cutting | Use eval/manual/fix on selected contigs. |
| `duplicate_overlap` | Redundant reference span | Check whether this is acceptable loss from ordered output. |
| `terminal_overlap` | One-sided overlap not rescued | Review as scaffold/overlap context. |
| `ambiguous_ref_match` | No dominant reference | Review per-reference matches and possible repeats. |
| `no_alignment` | No usable rows | Check names, aligner settings, or true unplaced sequence. |

## Common Traps

Do not treat `duplicate_overlap` as proof that a sequence is biologically
unimportant. It is a reference-space redundancy decision in this sort run.

Do not split every `kept_split_candidate`. The status is a review flag, not a
verdict.

Do not mix sort reports and FASTA files from different stages. The assignment
report must match the ordered FASTA used by `chromo scaffold`.

Do not assume graph context changed the sort decision unless the command docs
state that a setting is sequence-changing. Most graph fields around sorting are
guardrails or report-only context.

Do not validate `ordered.fa` with the raw PAF or coords. Re-align the ordered
FASTA if the edited file itself needs validation.

## What To Look At Next In ChromoSort

- Use [Reading ChromoSort Audit Tables]({{ '/guides/audit-tables/' | relative_url }})
  when a status label needs table context.
- Use [How To Interpret Dot Plots]({{ '/dot-plots/' | relative_url }}) for
  visual patterns behind split candidates and duplicate overlaps.
- Use [Chimeric Contig And Breakpoint Review]({{ '/guides/breakpoint-review/' | relative_url }})
  before changing a suspicious contig.
- Use [chromo sort]({{ '/commands/sort/' | relative_url }}) for the full
  parameter reference.
