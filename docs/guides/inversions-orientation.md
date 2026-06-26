---
title: Inversions And Orientation Changes
description: Distinguish whole-contig reverse orientation, internal inversions, reference differences, and fixable orientation errors.
---

# Inversions And Orientation Changes

Use this guide when a dot plot has red reverse-strand alignments, an internal
orientation switch, or a contig that appears backwards relative to the
reference.

The decision is:

> Is this only orientation, a real inversion, a reference difference, or a
> sequence error that should be edited?

## The Core Idea

Orientation is not automatically an assembly error. Assemblers can emit either
DNA strand. A whole contig may be reverse-complemented relative to the
reference and still be perfectly valid. Internal orientation changes need more
care: they can represent a real inversion allele, a reference/sample
difference, an assembly error, or repeat-mediated alignment noise.

ChromoSort can orient whole retained contigs to the reference during sorting.
It can also expose orientation-change candidates through `chromo eval fix` or
`chromo fix` modes. But true inversions should usually be reviewed, labeled,
or preserved rather than automatically split or reference-normalized.

## Pattern Gallery

### Whole-Contig Reverse Orientation

What it looks like: one long red alignment segment spans most of a contig and
one reference region.

Most likely interpretation: the contig is reverse-complemented relative to the
reference.

Best action:

- keep the contig as one unit,
- use `chromo sort --orient-to-reference` if you want reference-oriented output,
- do not use `chromo fix` only because the slope is negative.

### Mixed Orientation Inside One Contig

What it looks like: forward flanks with a reverse internal block, or several
orientation switches inside the same reference assignment.

Most likely interpretation: possible internal inversion, local assembly error,
repeat-driven alignment, or reference difference.

Best action:

- inspect segment length, identity, MAPQ, and breakpoint sharpness,
- compare coords and PAF if both exist,
- add read-to-assembly PAF, GFA, and GAF evidence when available,
- use `chromo eval fix --mode comprehensive` to expose orientation-change rows
  for review without applying them.

### Reverse Block On A Different Reference

What it looks like: one part of a contig aligns to one reference and another
part aligns in either orientation to a different reference.

Most likely interpretation: split candidate, translocation relative to the
reference, unresolved repeat, contamination, or imperfect reference.

Best action:

- treat it as a chimeric-contig review candidate,
- inspect per-reference match rows and dot plots,
- use conservative fix review before applying sequence changes.

### Many Small Red And Blue Fragments

What it looks like: scattered short alignments in both orientations.

Most likely interpretation: repeats, paralogs, low-specificity alignment, or
overly permissive aligner settings.

Best action:

- check PAF MAPQ and secondary alignment settings,
- raise minimum segment length or identity filters for review plots,
- avoid making breakpoint decisions from tiny fragments alone.

## Command Choices

| Goal | Command | Why |
| --- | --- | --- |
| Orient whole retained contigs to match the reference | `chromo sort --orient-to-reference` | Changes strand of whole retained records after assignment. |
| See an orientation pattern without changing FASTA | `chromo plot` | Draws the evidence in the original coordinate systems. |
| Review internal orientation transitions as candidates | `chromo eval fix --mode comprehensive` | Writes review rows without applying sequence edits. |
| Apply reviewed split pieces | `chromo fix --reviewed-plan` | Uses accepted source slices from a reviewed table. |
| Manually invert a piece for a specific curated output | `chromo manual` or `manual apply` | Keeps the decision explicit and reproducible. |
| Cut exact known coordinates from outside evidence | `chromo cut` | Applies known positions without using an alignment planner. |

## Review Questions

Ask these in order:

1. Is the entire contig reversed, or only an internal block?
2. Does one reference assignment dominate the whole contig?
3. Are the alignment blocks long, high identity, and high MAPQ?
4. Are the putative breakpoints sharp and consistent across evidence sources?
5. Do long reads span the breakpoints in the assembly coordinate system?
6. Does graph context support the assembly path or suggest a junction problem?
7. Is the goal to preserve the sample haplotype or make a
   reference-normalized experimental FASTA?

The last question matters. For pangenome graph inputs, preserving real
haplotype structure is usually the right default. A reference-normalized FASTA
can be useful for a specific comparison, but it should be labeled as such.

## Example Review Workflow

Generate a review table instead of immediately editing sequence:

```bash
chromo eval fix \
  --assembly-fasta assembly.fa \
  --paf assembly.paf \
  --contigs contig_with_orientation_switch \
  --mode comprehensive \
  --read-paf reads_to_assembly.paf \
  --gfa assembly_graph.gfa \
  --gaf reads_to_graph.gaf \
  --output-prefix review/sample.inversion
```

Then review:

- the dot plot for the contig,
- `review/sample.inversion.fix_review.tsv`,
- long-read spanning and split-read support near candidate breakpoints,
- graph node and path context,
- whether another assembly or reference shows the same orientation.

Only apply a sequence edit if the reviewed goal requires it.

## Cheat Sheet

| Pattern | Usually do |
| --- | --- |
| One long reverse-strand contig | Orient with `chromo sort --orient-to-reference` if desired. |
| Internal reverse block with strong support | Review as possible inversion; preserve unless the goal says otherwise. |
| Internal reverse block with weak or tiny rows | Treat as noise until stronger evidence appears. |
| Orientation switch plus different reference | Review as a split candidate. |
| Confirmed biological inversion | Report or preserve in native assembly; avoid silent reference normalization. |
| Curated reference-normalized output needed | Use manual or reviewed fix workflow and label the output clearly. |

## Common Traps

Do not call every red segment an error. Red means reverse-strand alignment, not
bad sequence.

Do not split a whole-contig reverse alignment. Whole-record orientation is a
sort or manual-orientation decision.

Do not use reference order as the only truth. A sample can carry real structure
that differs from the reference.

Do not let tiny repeat hits create an internal inversion story. Require strong
blocks and independent support.

Do not silently overwrite a native haplotype assembly with a
reference-normalized edit. Keep that output separate and documented.

## What To Look At Next In ChromoSort

- Use [How To Interpret Dot Plots]({{ '/dot-plots/' | relative_url }}#internal-inversion)
  for visual examples.
- Use [Chimeric Contig And Breakpoint Review]({{ '/guides/breakpoint-review/' | relative_url }})
  when orientation changes might become fix candidates.
- Use the [Agent And Review Playbook]({{ '/review-playbook/' | relative_url }}#review-same-reference-inversions)
  for pangenome-aware inversion review.
- Use [chromo manual]({{ '/commands/manual/' | relative_url }}) when the case
  needs browser review and an explicit recipe.
