---
title: How To Decide When To Fix A Contig
description: Conceptual guide to deciding when a contig should be split, cut, or left unchanged.
---

# How To Decide When To Fix A Contig

Fixing a contig means changing sequence. It is the most consequential kind of
assembly curation, so the question is not simply "does the plot look odd?" The
better question is:

> Is there strong evidence that one assembled sequence joins pieces that should
> not be joined in this assembly?

A good fix decision separates observation from action. A dot plot may show a
multi-reference pattern, an internal orientation change, or a sudden jump in
reference position. Those are observations. The action depends on whether the
pattern is coherent, local, supported, and unlikely to be explained by repeats,
reference differences, true structural variation, or stale alignments.

## The Core Idea

A contig is a candidate for fixing when its own coordinate order says one
continuous assembled molecule, but the evidence says two or more incompatible
genomic neighborhoods.

<figure class="concept-figure">
<svg viewBox="0 0 920 420" role="img" aria-labelledby="fix-title fix-desc">
  <title id="fix-title">Conceptual diagram of a chimeric contig decision</title>
  <desc id="fix-desc">A contig has one segment aligned to chromosome 3 and another segment aligned to chromosome 11, with a breakpoint between them. The reviewed fix emits two pieces.</desc>
  <defs>
    <marker id="fix-arrow" markerWidth="9" markerHeight="9" refX="7" refY="3" orient="auto" markerUnits="strokeWidth">
      <path d="M0,0 L0,6 L8,3 z" fill="#11443b" />
    </marker>
  </defs>
  <rect x="20" y="20" width="880" height="380" rx="12" fill="#fbfbf8" stroke="#deded8" />
  <text x="52" y="62" font-size="22" font-weight="700" fill="#202124">One contig, two incompatible neighborhoods</text>
  <text x="52" y="94" font-size="15" fill="#5f6368">A fix is justified only when the split model explains the evidence better than leaving one joined contig.</text>

  <text x="72" y="145" font-size="16" font-weight="700" fill="#202124">Reference evidence</text>
  <rect x="72" y="165" width="300" height="22" rx="11" fill="#dceee7" />
  <text x="86" y="181" font-size="13" fill="#11443b">Chr03 region</text>
  <rect x="500" y="165" width="300" height="22" rx="11" fill="#eee2d4" />
  <text x="514" y="181" font-size="13" fill="#6c4622">Chr11 region</text>

  <text x="72" y="237" font-size="16" font-weight="700" fill="#202124">Assembly contig</text>
  <rect x="72" y="257" width="728" height="36" rx="18" fill="#f4f6f5" stroke="#deded8" />
  <rect x="88" y="265" width="300" height="20" rx="10" fill="#206a5d" />
  <rect x="486" y="265" width="300" height="20" rx="10" fill="#b66b33" />
  <line x1="437" y1="246" x2="437" y2="304" stroke="#202124" stroke-width="2" stroke-dasharray="5 5" />
  <text x="392" y="327" font-size="13" fill="#202124">candidate breakpoint</text>

  <path d="M238 191 C238 220, 238 232, 238 257" fill="none" stroke="#206a5d" stroke-width="3" marker-end="url(#fix-arrow)" />
  <path d="M650 191 C650 220, 650 232, 650 257" fill="none" stroke="#b66b33" stroke-width="3" marker-end="url(#fix-arrow)" />

  <text x="72" y="365" font-size="16" font-weight="700" fill="#202124">Reviewed fix output</text>
  <rect x="270" y="348" width="180" height="28" rx="14" fill="#206a5d" />
  <text x="292" y="367" font-size="13" fill="#ffffff">piece on Chr03</text>
  <rect x="500" y="348" width="180" height="28" rx="14" fill="#b66b33" />
  <text x="522" y="367" font-size="13" fill="#ffffff">piece on Chr11</text>
  <line x1="450" y1="362" x2="500" y2="362" stroke="#11443b" stroke-width="2" marker-end="url(#fix-arrow)" />
</svg>
<figcaption><strong>Figure 1. Fixing is a model choice.</strong> The joined contig is one model. The split pieces are another. A fix is appropriate when the split model is better supported by long coherent alignments, local breakpoints, and independent evidence.</figcaption>
</figure>

## Why We Think A Fix Might Work

Reference-guided fixing works because chromosome-scale alignments are often
long enough to reveal incompatible neighborhoods. If the left side of a contig
aligns cleanly to one reference region and the right side aligns cleanly to a
different region, a single breakpoint can explain the pattern with fewer
assumptions than keeping the contig joined.

That reasoning is conservative only when the evidence is conservative:

- the alignment file describes the exact FASTA being reviewed,
- the suspicious blocks are long and coherent,
- the breakpoint is local rather than smeared across many tiny hits,
- repeats and secondary hits are not driving the call,
- optional graph or read evidence does not contradict the split,
- the output pieces would each have enough support to stand alone.

## A Decision Ladder

| Question | Why it matters |
| --- | --- |
| Does the evidence match this exact FASTA? | Old alignments can make already-edited contigs look suspicious. |
| Is the pattern large enough to matter? | Tiny off-target hits are often repeats or aligner noise. |
| Are there two or more coherent neighborhoods? | A fix should explain real blocks, not scattered speckles. |
| Is there a local breakpoint? | Sequence edits need boundaries that can be reviewed and reproduced. |
| Could this be real biology or a reference difference? | True inversions, translocations, or presence/absence differences should not be erased casually. |
| Would the pieces be useful after the fix? | A split that creates unsupported fragments is rarely an improvement. |

## When To Fix

Consider a reviewed fix when a contig has strong blocks on different
references, distant incompatible locations on one reference, or a sharp
same-contig transition that is also supported by graph or read evidence.

Do not fix just because a contig is reversed relative to the reference. A
whole-contig reverse alignment is usually an orientation problem. Do not split
a clean internal inversion simply to make the assembly look like the reference.
An inversion may be real biology, a reference difference, or a case for
special reporting rather than sequence cutting.

## Example

A soybean contig has 11 Mb of high-identity alignment to chromosome 3, then a
sharp transition, then 8 Mb of high-identity alignment to chromosome 11. Both
blocks use most of their local contig spans, and the transition is not made of
many small repeat hits.

The concept-level interpretation is:

1. One contig is carrying evidence from two incompatible neighborhoods.
2. A single reviewed breakpoint explains the pattern.
3. The split pieces should be checked as independent contigs.
4. The fixed FASTA must be re-aligned before sorting or validating it.

The ChromoSort-specific next steps are
[Chimeric contig and breakpoint review]({{ '/guides/breakpoint-review/' | relative_url }}),
[Sort, clean, fix, cut, or manual?]({{ '/guides/choosing-commands/' | relative_url }}),
and [Spreadsheet review tables]({{ '/guides/spreadsheet-review-tables/' | relative_url }}).
