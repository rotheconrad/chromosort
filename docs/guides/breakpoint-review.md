---
title: Chimeric Contig And Breakpoint Review
description: Review split candidates, chromo fix modes, breakpoint planning, and reviewed-plan workflows before changing contigs.
---

# Chimeric Contig And Breakpoint Review

Use this guide when a contig may be chimeric or structurally inconsistent, and
you need to decide whether `chromo fix` should split it.

The key question is:

> Is the discordance strong enough, local enough, and well supported enough to
> cut this contig?

## The Core Idea

`chromo fix` is a breakpoint planner and executor for selected contigs. It
reads alignment evidence for one exact assembly FASTA, finds substantial
changes in reference assignment or orientation, and emits copied or split
pieces.

That does not mean every alignment transition should become a cut. Repeats,
low-MAPQ alignments, reference differences, true structural variants, and noisy
alignment fragments can all look like breakpoints. Good review separates the
observation from the action.

## What ChromoSort Reads Or Writes

For planning, `chromo fix` reads:

- `--assembly-fasta`,
- either `--coords` or `--paf` from that same assembly FASTA,
- selected contigs from `--contigs`, `--contigs-file`, or `--all`,
- optional GFA, GAF, or long-read PAF evidence for report context.

It writes:

- an output FASTA with unchanged copied contigs and split pieces,
- a fix report TSV describing source contigs, emitted pieces, slices, status,
  dominant reference, orientation, and reasons.

`chromo eval fix` uses the same planner but writes a review table instead of a
FASTA. Accepted `split_piece` rows can later be applied with
`chromo fix --reviewed-plan`.

## What Counts As A Split Candidate

| Evidence pattern | What it might mean | Review stance |
| --- | --- | --- |
| One contig has large blocks on different references | Misjoin, translocation, repeat, contamination, or reference difference | High-priority review candidate. |
| Same reference, distant blocks with inconsistent order | Assembly error, structural variant, repeat, or stale alignment | Review plot boundaries and read support. |
| Same reference, internal orientation switch | Inversion, local orientation error, or reference difference | Do not automatically split; review as an inversion case. |
| Many tiny off-target hits around a strong primary hit | Repeat or aligner noise | Raise thresholds or leave unchanged unless support improves. |
| Transition near contig ends only | Terminal overlap or partial alignment issue | Usually a sort/scaffold review issue, not a fix breakpoint. |

## Fix Modes

`chromo fix` modes define which transitions can become candidate breakpoints.

| Mode | What it exposes | Smoothing | Typical use |
| --- | --- | --- | --- |
| `chromosome` | Reference or chromosome changes only. | Yes | Very conservative multi-reference review. |
| `conservative` | Reference/chromosome changes plus complex same-reference orientation events. | Yes | Default production-oriented mode. |
| `comprehensive` | Reference/chromosome transitions and same-reference orientation changes after smoothing. | Yes | Review broader structural signals, including possible inversions. |
| `sensitive` | Every passing reference or orientation transition after same-target collapse. | No | Diagnostic pass when you want to see noisy candidates. |

Comprehensive mode is orientation-aware after smoothing, but it is not a
guaranteed superset of conservative mode. Compare reports when a case is
borderline.

## How Breakpoints Are Planned

The planner:

1. reads usable alignment segments,
2. orders them along the query contig,
3. merges nearby compatible segments,
4. collapses adjacent same-target evidence,
5. applies the selected mode,
6. smooths noisy transitions when the mode uses smoothing,
7. caps the number of breakpoints,
8. places each accepted breakpoint halfway between adjacent alignment blocks.

Important thresholds include minimum segment support, minimum piece support,
breakpoint penalty, and maximum breakpoints. If a contig reports
`not_split_smooth`, the planner saw discordance but smoothed it away under the
current settings. If it reports `not_split_too_many_breakpoints`, the event is
too fragmented for the current breakpoint budget.

## Review Workflow

### 1. Start From Selected Contigs

Prefer targeted contig lists for production review:

```bash
chromo eval fix \
  --assembly-fasta raw.fa \
  --paf raw.paf \
  --contigs suspect_contig_1 suspect_contig_2 \
  --mode conservative \
  --output-prefix review/sample.fix
```

Use `--all` for discovery or controlled test runs, then narrow the list before
applying sequence changes.

### 2. Inspect The Review Table

Open `review/sample.fix.fix_review.tsv`.

Read these fields together:

| Field family | Why it matters |
| --- | --- |
| Source contig and slice coordinates | Confirm which interval would be emitted. |
| Event type and accepted fields | Confirm whether a row is eligible for application. |
| Dominant reference and orientation | Compare each piece to the dot plot. |
| Planner reason or status | Understand smoothing, threshold, or breakpoint-budget decisions. |
| `graph_*`, `gaf_*`, and `longread_*` fields | Use optional evidence as context, not automatic permission. |

### 3. Compare Against The Plot

The table should match the visual pattern:

- blocks should be long enough to trust,
- boundaries should be reasonably sharp,
- pieces should have enough aligned support,
- off-target repeat hits should not drive the cut,
- same-reference inversions should get extra caution.

### 4. Apply Reviewed Rows

```bash
chromo fix \
  --assembly-fasta raw.fa \
  --reviewed-plan review/sample.fix.fix_review.tsv \
  --output-fasta results/sample.fixed.fa \
  --report results/sample.fixed.tsv
```

When `--reviewed-plan` is used, the reviewed table supplies exact source slices.
The executor revalidates rows so stale or unaccepted decisions do not silently
change sequence.

### 5. Re-Align The Fixed FASTA

```text
results/sample.fixed.fa
  -> minimap2 or MUMmer against the reference
  -> validation plot
  -> chromo sort if the fixed assembly looks right
```

Old raw alignment rows are useful provenance, but they do not validate the
fixed FASTA.

## When To Use Manual Review Instead

Use `chromo manual` when:

- the candidate has several possible breakpoints,
- an internal inversion might be real biology,
- graph evidence and read evidence disagree,
- pieces need to be removed, restored, reordered, or inverted by judgment,
- the case is easier to inspect in a browser than in a spreadsheet.

`chromo manual fix --review-table` can load the eval table as a focused queue
while still letting you browse nearby contigs and graph context.

## Cheat Sheet

| If you see... | First response |
| --- | --- |
| Strong two-reference contig | `chromo eval fix --mode conservative` on that contig. |
| Whole-contig reverse alignment | Use sort orientation; do not fix. |
| Same-reference internal inversion | Review with comprehensive eval, reads, graph context, and the inversion guide. |
| Many tiny candidate pieces | Raise thresholds, inspect repeats, or use manual review. |
| `not_split_smooth` | Planner decided the transition was not worth a breakpoint under current settings. |
| `not_split_too_many_breakpoints` | Candidate is too fragmented for the breakpoint budget. |
| Exact externally known breakpoint | Use `chromo cut`, not planner discovery. |

## Common Traps

Do not run `chromo fix --all` on a new dataset and treat the output as final.
Use it first as discovery, then review.

Do not split true biological inversions just to make the dot plot look like the
reference.

Do not ignore `not_split_*` rows. They often explain why a suspicious plot did
not become a FASTA edit.

Do not apply a reviewed table after changing the assembly FASTA or regenerating
evidence with different source names.

Do not forget that optional graph and long-read evidence is usually advisory in
fix review. It should inform the decision, not bypass it.

## What To Look At Next In ChromoSort

- Use [Sort, Clean, Fix, Cut, Or Manual?]({{ '/guides/choosing-commands/' | relative_url }})
  if you are unsure whether a split is the right action.
- Use [Inversions And Orientation Changes]({{ '/guides/inversions-orientation/' | relative_url }})
  before cutting same-reference orientation events.
- Use [chromo eval]({{ '/commands/eval/' | relative_url }}) for reviewed table
  workflows.
- Use [chromo fix]({{ '/commands/fix/' | relative_url }}) for the full planner
  and parameter reference.
