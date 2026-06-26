---
title: How Scaffolding Works
description: Conceptual guide to ordering contigs, orienting them, and representing unknown gaps with scaffold Ns.
---

# How Scaffolding Works

Scaffolding turns ordered contigs into larger records. In a chromosome-scale
assembly workflow, the scaffold is often the file a researcher wants to inspect
or submit: one record per chromosome, linkage group, or reference sequence.

The key point is that scaffolding joins existing sequence. It does not discover
the missing bases between contigs. When the sequence between adjacent contigs
is unknown, the scaffold uses Ns as a placeholder.

## The Core Idea

Scaffolding uses placement evidence to decide three things for adjacent contigs:

1. order,
2. orientation,
3. gap or overlap representation.

<figure class="concept-figure">
<svg viewBox="0 0 920 420" role="img" aria-labelledby="scaffold-title scaffold-desc">
  <title id="scaffold-title">Conceptual diagram of scaffolding ordered contigs</title>
  <desc id="scaffold-desc">Three contigs are placed along a reference, oriented, joined into a scaffold, and separated by N gaps where sequence is unknown.</desc>
  <defs>
    <marker id="scaf-arrow" markerWidth="9" markerHeight="9" refX="7" refY="3" orient="auto" markerUnits="strokeWidth">
      <path d="M0,0 L0,6 L8,3 z" fill="#11443b" />
    </marker>
  </defs>
  <rect x="20" y="20" width="880" height="380" rx="12" fill="#fbfbf8" stroke="#deded8" />
  <text x="52" y="62" font-size="22" font-weight="700" fill="#202124">Scaffolding joins placed contigs with explicit uncertainty</text>
  <text x="52" y="94" font-size="15" fill="#5f6368">The output is a larger coordinate frame, not proof that every base between contigs is known.</text>

  <text x="72" y="145" font-size="16" font-weight="700" fill="#202124">Reference placement</text>
  <rect x="72" y="165" width="720" height="20" rx="10" fill="#e8f1ee" stroke="#cfe0da" />
  <rect x="106" y="158" width="170" height="34" rx="8" fill="#206a5d" />
  <text x="150" y="180" font-size="13" fill="#ffffff">contig A</text>
  <rect x="346" y="158" width="145" height="34" rx="8" fill="#4f6d8a" />
  <text x="382" y="180" font-size="13" fill="#ffffff">contig B</text>
  <rect x="598" y="158" width="156" height="34" rx="8" fill="#b66b33" />
  <text x="636" y="180" font-size="13" fill="#ffffff">contig C</text>
  <line x1="276" y1="175" x2="346" y2="175" stroke="#5f6368" stroke-width="2" stroke-dasharray="5 5" />
  <line x1="491" y1="175" x2="598" y2="175" stroke="#5f6368" stroke-width="2" stroke-dasharray="5 5" />
  <text x="288" y="211" font-size="13" fill="#5f6368">unknown gap</text>
  <text x="515" y="211" font-size="13" fill="#5f6368">unknown gap</text>

  <path d="M430 225 L430 260" stroke="#11443b" stroke-width="3" marker-end="url(#scaf-arrow)" />

  <text x="72" y="295" font-size="16" font-weight="700" fill="#202124">Scaffold output</text>
  <rect x="192" y="318" width="170" height="34" rx="8" fill="#206a5d" />
  <text x="236" y="340" font-size="13" fill="#ffffff">contig A</text>
  <rect x="362" y="322" width="66" height="26" rx="4" fill="#f3f5f4" stroke="#b8c3bf" stroke-dasharray="4 4" />
  <text x="383" y="340" font-size="13" fill="#5f6368">Ns</text>
  <rect x="428" y="318" width="145" height="34" rx="8" fill="#4f6d8a" />
  <text x="464" y="340" font-size="13" fill="#ffffff">contig B</text>
  <rect x="573" y="322" width="66" height="26" rx="4" fill="#f3f5f4" stroke="#b8c3bf" stroke-dasharray="4 4" />
  <text x="594" y="340" font-size="13" fill="#5f6368">Ns</text>
  <rect x="639" y="318" width="156" height="34" rx="8" fill="#b66b33" />
  <text x="677" y="340" font-size="13" fill="#ffffff">contig C</text>
</svg>
<figcaption><strong>Figure 1. Scaffolding is ordering plus uncertainty.</strong> Contigs become one larger scaffold record. Ns mark sequence that is not known from the contigs themselves.</figcaption>
</figure>

## Why We Think Scaffolding Works

Scaffolding works when independent evidence agrees on the relative placement of
contigs. A close reference can provide approximate order and orientation. Hi-C,
linkage maps, long reads, or assembly graph links can add support when they are
available. The scaffold is strongest when those evidence layers point to the
same adjacency.

The result is still a model. It is useful because it gives downstream analyses
a chromosome-scale coordinate system, but it should keep uncertainty visible
with gap records, Ns, AGP rows, or equivalent provenance.

## What Scaffolding Changes

Scaffolding can:

- order contigs,
- orient contigs,
- join contigs into larger FASTA records,
- insert N gaps for unknown sequence,
- report inferred gaps and overlaps.

Scaffolding should not:

- pretend Ns are known bases,
- silently remove real sequence,
- resolve an ambiguous graph branch,
- fix a chimeric contig that should have been reviewed earlier.

## Gaps And Overlaps

Adjacent contigs may leave a positive gap, touch exactly, or overlap in
reference coordinates. A positive gap usually becomes Ns. A negative gap is an
overlap and needs a policy: leave a zero-length gap, warn, trim by reference
coordinates, or trim only when the contig sequences confirm the overlap.

This is why scaffolding reports matter. The FASTA alone cannot tell you whether
a junction came from an inferred reference gap, a fixed convention, a reviewed
override, or a detected overlap.

## Example

Three contigs align to the same chromosome in order: A, B, C. The right end of
A maps before the left end of B, leaving an inferred 12 kb reference gap. B and
C overlap by 400 bp in reference coordinates, but their sequences do not match
well at the overlap.

The concept-level scaffold would join A and B with Ns, then join B and C with
a zero-length or reviewed overlap decision rather than trimming sequence
automatically.

The ChromoSort-specific next steps are
[Scaffolding, gaps, and overlaps]({{ '/guides/scaffolding-gaps-overlaps/' | relative_url }}),
[Graph-supported gap filling]({{ '/guides/graph-gap-filling/' | relative_url }}),
and [chromo scaffold]({{ '/commands/scaffold/' | relative_url }}).
