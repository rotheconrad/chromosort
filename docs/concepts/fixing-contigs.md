---
title: How To Decide When To Fix A Contig
description: Conceptual lesson for deciding when a contig should be split, cut, or left unchanged.
---

# How To Decide When To Fix A Contig

Fixing a contig means changing sequence. That makes it one of the most useful
and most dangerous parts of assembly curation. A suspicious plot is not enough.
The better question is:

> Is there strong evidence that one assembled sequence joins pieces that should
> not be joined in this assembly?

This lesson explains the reasoning before the ChromoSort command details. The
goal is to teach when a split model is more believable than the current contig,
why the evidence can support that conclusion, and why many odd-looking patterns
should still be left alone.

## The Core Idea

A contig is a candidate for fixing when its own coordinate order says "one
continuous assembled molecule," but the evidence says "two or more incompatible
genomic neighborhoods."

<figure class="concept-figure">
<svg viewBox="0 0 920 420" role="img" aria-labelledby="fix-model-title fix-model-desc">
  <title id="fix-model-title">Joined contig model compared with split contig model</title>
  <desc id="fix-model-desc">One contig contains a left interval aligned to chromosome 3 and a right interval aligned to chromosome 11. A reviewed split emits two pieces.</desc>
  <rect x="20" y="20" width="880" height="380" rx="12" fill="#fbfbf8" stroke="#deded8" />
  <text x="52" y="62" font-size="22" font-weight="700" fill="#202124">Fixing is a model choice</text>
  <text x="52" y="94" font-size="15" fill="#5f6368">The joined contig is one model. The split pieces are another.</text>

  <text x="72" y="142" font-size="16" font-weight="700" fill="#202124">Reference evidence</text>
  <rect x="72" y="163" width="300" height="24" rx="12" fill="#dceee7" />
  <text x="86" y="181" font-size="13" fill="#11443b">Chr03 neighborhood</text>
  <rect x="500" y="163" width="300" height="24" rx="12" fill="#eee2d4" />
  <text x="514" y="181" font-size="13" fill="#6c4622">Chr11 neighborhood</text>

  <text x="72" y="235" font-size="16" font-weight="700" fill="#202124">Current assembly model</text>
  <rect x="72" y="257" width="728" height="36" rx="18" fill="#f4f6f5" stroke="#deded8" />
  <rect x="88" y="265" width="300" height="20" rx="10" fill="#206a5d" />
  <rect x="486" y="265" width="300" height="20" rx="10" fill="#b66b33" />
  <line x1="437" y1="246" x2="437" y2="304" stroke="#202124" stroke-width="2" stroke-dasharray="5 5" />
  <text x="386" y="327" font-size="13" fill="#202124">candidate breakpoint</text>

  <path d="M238 191 C238 220, 238 232, 238 257" fill="none" stroke="#206a5d" stroke-width="3" />
  <polygon points="238,257 232,246 244,246" fill="#206a5d" />
  <path d="M650 191 C650 220, 650 232, 650 257" fill="none" stroke="#b66b33" stroke-width="3" />
  <polygon points="650,257 644,246 656,246" fill="#b66b33" />

  <text x="72" y="365" font-size="16" font-weight="700" fill="#202124">Reviewed fix model</text>
  <rect x="270" y="348" width="180" height="28" rx="14" fill="#206a5d" />
  <text x="292" y="367" font-size="13" fill="#ffffff">piece on Chr03</text>
  <rect x="500" y="348" width="180" height="28" rx="14" fill="#b66b33" />
  <text x="522" y="367" font-size="13" fill="#ffffff">piece on Chr11</text>
</svg>
<figcaption><strong>Figure 1. Fixing is a model choice.</strong> The split model is appropriate only when it explains the evidence better than the original joined contig.</figcaption>
</figure>

The split model becomes convincing when the alignment blocks are long, coherent,
and separated by a local boundary. It becomes weak when the evidence is made of
small repeat-like hits, stale alignments, low-confidence mappings, or patterns
that are better explained as true structural variation.

## Observation, Interpretation, Action

Good review has three steps:

1. Observation: what pattern is visible?
2. Interpretation: what biological or technical situations could create it?
3. Action: should sequence be changed, reviewed manually, or left alone?

The same observation can lead to different actions. A contig with blocks on two
references might be a misjoin, a real translocation relative to the reference,
an unresolved repeat, contamination, or a stale alignment. Fixing is justified
only after the competing explanations have been narrowed.

## The Evidence Ladder

Do not jump from "odd plot" to "split contig." Climb the evidence ladder.

<figure class="concept-figure">
<svg viewBox="0 0 920 420" role="img" aria-labelledby="fix-ladder-title fix-ladder-desc">
  <title id="fix-ladder-title">Evidence ladder for fixing contigs</title>
  <desc id="fix-ladder-desc">A ladder shows increasing confidence from exact FASTA alignment, coherent blocks, local boundary, independent support, reviewed action, and re-alignment after editing.</desc>
  <rect x="20" y="20" width="880" height="380" rx="12" fill="#fbfbf8" stroke="#deded8" />
  <text x="52" y="62" font-size="22" font-weight="700" fill="#202124">Evidence should accumulate before sequence changes</text>
  <text x="52" y="94" font-size="15" fill="#5f6368">Each step removes a common false-positive explanation.</text>

  <g font-size="14">
    <rect x="90" y="322" width="740" height="42" rx="8" fill="#eef6f2" stroke="#cfe0da" />
    <text x="112" y="348" fill="#11443b">6. Re-align the edited FASTA before validation or downstream sorting</text>

    <rect x="120" y="270" width="680" height="42" rx="8" fill="#eef6f2" stroke="#cfe0da" />
    <text x="142" y="296" fill="#11443b">5. Accept an explicit reviewed action: split, manual review, or leave unchanged</text>

    <rect x="150" y="218" width="620" height="42" rx="8" fill="#eef6f2" stroke="#cfe0da" />
    <text x="172" y="244" fill="#11443b">4. Check independent support: graph, long reads, read pairs, or another reference</text>

    <rect x="180" y="166" width="560" height="42" rx="8" fill="#eef6f2" stroke="#cfe0da" />
    <text x="202" y="192" fill="#11443b">3. Locate a plausible local boundary rather than a smear of tiny hits</text>

    <rect x="210" y="114" width="500" height="42" rx="8" fill="#eef6f2" stroke="#cfe0da" />
    <text x="232" y="140" fill="#11443b">2. Confirm long coherent blocks in incompatible neighborhoods</text>

    <rect x="240" y="62" width="440" height="42" rx="8" fill="#eef6f2" stroke="#cfe0da" />
    <text x="262" y="88" fill="#11443b">1. Confirm the alignment matches this exact FASTA</text>
  </g>
</svg>
<figcaption><strong>Figure 2. The evidence ladder.</strong> The higher the consequence of the action, the more provenance and support the decision needs.</figcaption>
</figure>

The first rung matters most. If an alignment was generated from `raw.fa`, it
describes `raw.fa`. It does not validate `fixed.fa`, `ordered.fa`, a manual
FASTA export, or a scaffold FASTA.

## Pattern Gallery

### Strong Fix Candidate

A strong candidate has a small number of large blocks that disagree with the
single-contig model.

<figure class="concept-figure">
<svg viewBox="0 0 920 420" role="img" aria-labelledby="fix-good-title fix-good-desc">
  <title id="fix-good-title">Strong fix candidate patterns</title>
  <desc id="fix-good-desc">Two cartoon dot plots show a multi-reference contig and a same-reference jump with a sharp boundary.</desc>
  <rect x="20" y="20" width="880" height="380" rx="12" fill="#fbfbf8" stroke="#deded8" />
  <text x="52" y="62" font-size="22" font-weight="700" fill="#202124">Patterns that deserve split review</text>

  <text x="92" y="112" font-size="16" font-weight="700" fill="#202124">A. Multi-reference blocks</text>
  <rect x="92" y="132" width="320" height="220" fill="#ffffff" stroke="#deded8" />
  <line x1="122" y1="322" x2="382" y2="322" stroke="#5f6368" />
  <line x1="122" y1="322" x2="122" y2="162" stroke="#5f6368" />
  <text x="170" y="345" font-size="12" fill="#5f6368">reference position</text>
  <text x="64" y="244" font-size="12" fill="#5f6368" transform="rotate(-90 64 244)">same contig</text>
  <line x1="140" y1="292" x2="236" y2="230" stroke="#206a5d" stroke-width="8" stroke-linecap="round" />
  <line x1="270" y1="210" x2="366" y2="152" stroke="#b66b33" stroke-width="8" stroke-linecap="round" />
  <rect x="134" y="174" width="110" height="22" rx="11" fill="#dceee7" />
  <text x="151" y="190" font-size="12" fill="#11443b">Chr03</text>
  <rect x="262" y="250" width="110" height="22" rx="11" fill="#eee2d4" />
  <text x="279" y="266" font-size="12" fill="#6c4622">Chr11</text>

  <text x="500" y="112" font-size="16" font-weight="700" fill="#202124">B. Same-reference jump</text>
  <rect x="500" y="132" width="320" height="220" fill="#ffffff" stroke="#deded8" />
  <line x1="530" y1="322" x2="790" y2="322" stroke="#5f6368" />
  <line x1="530" y1="322" x2="530" y2="162" stroke="#5f6368" />
  <line x1="548" y1="292" x2="645" y2="230" stroke="#206a5d" stroke-width="8" stroke-linecap="round" />
  <line x1="676" y1="206" x2="775" y2="154" stroke="#206a5d" stroke-width="8" stroke-linecap="round" />
  <line x1="658" y1="170" x2="658" y2="330" stroke="#202124" stroke-width="2" stroke-dasharray="5 5" />
  <text x="610" y="353" font-size="12" fill="#202124">sharp transition</text>
  <path d="M648 230 L676 206" stroke="#202124" stroke-width="2" stroke-dasharray="5 5" />
</svg>
<figcaption><strong>Figure 3. Strong candidates have coherent blocks.</strong> A multi-reference split or a sharp same-reference jump deserves review when the blocks are long and the boundary is local.</figcaption>
</figure>

Review stance: evaluate a split plan, inspect the boundary, and ask whether
each emitted piece would still have enough support.

### Usually Not A Fix

Some patterns are real and important, but still should not be cut.

<figure class="concept-figure">
<svg viewBox="0 0 920 420" role="img" aria-labelledby="fix-not-title fix-not-desc">
  <title id="fix-not-title">Patterns that usually should not be fixed by splitting</title>
  <desc id="fix-not-desc">Two cartoon dot plots show a whole-contig reverse alignment and an internal inversion pattern.</desc>
  <rect x="20" y="20" width="880" height="380" rx="12" fill="#fbfbf8" stroke="#deded8" />
  <text x="52" y="62" font-size="22" font-weight="700" fill="#202124">Odd is not the same as broken</text>

  <text x="92" y="112" font-size="16" font-weight="700" fill="#202124">A. Whole-contig reverse alignment</text>
  <rect x="92" y="132" width="320" height="220" fill="#ffffff" stroke="#deded8" />
  <line x1="122" y1="322" x2="382" y2="322" stroke="#5f6368" />
  <line x1="122" y1="322" x2="122" y2="162" stroke="#5f6368" />
  <line x1="146" y1="170" x2="360" y2="304" stroke="#a44747" stroke-width="9" stroke-linecap="round" />
  <text x="126" y="352" font-size="13" fill="#5f6368">Usually orient, do not split.</text>

  <text x="500" y="112" font-size="16" font-weight="700" fill="#202124">B. Internal inversion pattern</text>
  <rect x="500" y="132" width="320" height="220" fill="#ffffff" stroke="#deded8" />
  <line x1="530" y1="322" x2="790" y2="322" stroke="#5f6368" />
  <line x1="530" y1="322" x2="530" y2="162" stroke="#5f6368" />
  <line x1="548" y1="292" x2="615" y2="250" stroke="#206a5d" stroke-width="8" stroke-linecap="round" />
  <line x1="628" y1="220" x2="704" y2="268" stroke="#a44747" stroke-width="8" stroke-linecap="round" />
  <line x1="716" y1="234" x2="774" y2="190" stroke="#206a5d" stroke-width="8" stroke-linecap="round" />
  <text x="520" y="352" font-size="13" fill="#5f6368">Review biology before changing sequence.</text>
</svg>
<figcaption><strong>Figure 4. Some discordance is not a split request.</strong> Whole-contig reverse orientation is usually an orientation choice. Internal inversions may be real biology or reference difference.</figcaption>
</figure>

Review stance: use orientation, manual review, additional references, graph
context, and read evidence before deciding whether any sequence edit is needed.

### Repeat Noise Or Stale Evidence

Repeat-rich genomes can produce small off-target hits that look dramatic when
compressed into a whole-genome plot. Stale alignments can make already-edited
FASTA records look broken because names or coordinates no longer match.

<figure class="concept-figure">
<svg viewBox="0 0 920 360" role="img" aria-labelledby="fix-noise-title fix-noise-desc">
  <title id="fix-noise-title">Repeat noise and stale evidence can mimic fix candidates</title>
  <desc id="fix-noise-desc">A strong main alignment is surrounded by small speckled hits and a warning that stale evidence compares the wrong FASTA stage.</desc>
  <rect x="20" y="20" width="880" height="320" rx="12" fill="#fbfbf8" stroke="#deded8" />
  <text x="52" y="62" font-size="22" font-weight="700" fill="#202124">Common false positives</text>

  <rect x="82" y="92" width="350" height="210" fill="#ffffff" stroke="#deded8" />
  <line x1="112" y1="272" x2="392" y2="272" stroke="#5f6368" />
  <line x1="112" y1="272" x2="112" y2="122" stroke="#5f6368" />
  <line x1="134" y1="252" x2="356" y2="132" stroke="#206a5d" stroke-width="9" stroke-linecap="round" />
  <g fill="#b66b33">
    <circle cx="166" cy="152" r="4" />
    <circle cx="228" cy="204" r="4" />
    <circle cx="314" cy="242" r="4" />
    <circle cx="370" cy="186" r="4" />
    <circle cx="194" cy="242" r="4" />
    <circle cx="284" cy="158" r="4" />
  </g>
  <text x="104" y="324" font-size="13" fill="#5f6368">Repeat speckles should not drive a cut.</text>

  <rect x="512" y="112" width="300" height="64" rx="8" fill="#f4f6f5" stroke="#deded8" />
  <text x="534" y="139" font-size="14" fill="#202124">raw.fa + raw.paf</text>
  <text x="534" y="160" font-size="13" fill="#5f6368">valid evidence for raw contigs</text>
  <path d="M662 181 L662 219" stroke="#a44747" stroke-width="3" />
  <polygon points="662,224 655,212 669,212" fill="#a44747" />
  <rect x="512" y="226" width="300" height="64" rx="8" fill="#fff3ef" stroke="#e1b9aa" />
  <text x="534" y="253" font-size="14" fill="#202124">fixed.fa + raw.paf</text>
  <text x="534" y="274" font-size="13" fill="#a44747">stale pairing: do not interpret</text>
</svg>
<figcaption><strong>Figure 5. False positives are common.</strong> Repeats, secondary hits, and stale FASTA/alignment pairings can make a harmless contig look suspicious.</figcaption>
</figure>

Review stance: raise filters, inspect per-reference plots, confirm the exact
FASTA pair, and avoid edits until the pattern remains under better evidence.

## A Practical Decision Table

| If you see... | First interpretation | Conservative action |
| --- | --- | --- |
| Strong blocks on different references | Possible misjoin, translocation, repeat, or contamination | Evaluate a reviewed split; inspect graph/read support. |
| Same reference, distant jump | Possible local misassembly or structural difference | Review boundary and compare with independent evidence. |
| Whole-contig reverse alignment | Orientation difference | Orient during sorting; do not split. |
| Blue-red-blue internal pattern | Possible inversion or reference difference | Review as inversion; do not automatically cut. |
| Many short off-target hits | Repeats, paralogs, or secondary alignments | Filter and inspect; do not cut from speckles. |
| Suspicious pattern after editing FASTA | Possible stale evidence | Re-align the exact edited FASTA before interpreting. |

## Example Walkthrough

Imagine a soybean contig with 11 Mb of high-identity alignment to chromosome 3,
then a sharp transition, then 8 Mb of high-identity alignment to chromosome 11.
Both blocks use most of their local contig spans. The transition is not made of
many tiny repeat hits.

The concept-level review is:

1. The current contig asserts a single joined molecule.
2. The evidence places the left and right intervals in incompatible reference
   neighborhoods.
3. A single split explains the pattern better than the joined model.
4. Each output piece would still have substantial support.
5. The fixed FASTA must be re-aligned before sorting or final validation.

That is a reasonable fix candidate. It is not yet an automatic edit. The
reviewer still needs an explicit accepted plan and provenance for the evidence
used.

## Common Traps

Do not split a contig just to make a plot prettier. A prettier reference-normal
plot can be a worse assembly if the sample truly differs from the reference.

Do not treat all same-reference orientation changes as errors. Inversions and
complex haplotype differences require more careful interpretation.

Do not apply a fix table after the FASTA changed. Breakpoint coordinates belong
to a particular source sequence.

Do not let a single evidence stream overrule obvious contradiction from graph,
read, or manual review evidence.

## Brief History And Further Reading

Early genome assembly quality work made an important point that still matters:
contiguity is not the same as correctness. Large contigs can contain structural
errors, and breaking or editing them should be justified by evidence rather
than by N50-style metrics alone.

Reference-based tools such as QUAST helped standardize language around
misassemblies, relocations, inversions, and translocations. Read-backed tools
such as REAPR emphasized that mapped reads can reveal local assembly problems
that reference comparison alone may miss. Modern assembly review usually
combines both ideas: compare to references, but keep read, graph, and k-mer
evidence nearby.

Further reading:

- Gurevich et al. 2013. [QUAST: quality assessment tool for genome assemblies](https://doi.org/10.1093/bioinformatics/btt086).
- Hunt et al. 2013. [REAPR: a universal tool for genome assembly evaluation](https://doi.org/10.1186/gb-2013-14-5-r47).
- Li 2018. [Minimap2: pairwise alignment for nucleotide sequences](https://doi.org/10.1093/bioinformatics/bty191).
- Rhie et al. 2020. [Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies](https://doi.org/10.1038/s41587-020-00702-1).

## What To Read Next

- [How to interpret dot plots]({{ '/dot-plots/' | relative_url }})
- [Chimeric contig and breakpoint review]({{ '/guides/breakpoint-review/' | relative_url }})
- [Sort, clean, fix, cut, or manual?]({{ '/guides/choosing-commands/' | relative_url }})
- [Spreadsheet review tables]({{ '/guides/spreadsheet-review-tables/' | relative_url }})
