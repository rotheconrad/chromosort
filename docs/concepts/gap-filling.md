---
title: How Gap Filling Works
description: Conceptual guide to replacing scaffold Ns with supported sequence from graph or read evidence.
---

# How Gap Filling Works

Gap filling tries to replace unknown scaffold sequence with actual bases. That
makes it more powerful than scaffolding, and more risky. A gap should be filled
only when the candidate sequence is connected to the correct flanks, supported
by evidence, and not confused with an equally plausible alternative.

The question is:

> Do we know the sequence that belongs between these two scaffold flanks?

## The Core Idea

Scaffolding says "these contigs are adjacent, but the sequence between them is
unknown." Gap filling says "this candidate path supplies the missing sequence."

<figure class="concept-figure">
<svg viewBox="0 0 920 430" role="img" aria-labelledby="gap-title gap-desc">
  <title id="gap-title">Conceptual diagram of graph-supported gap filling</title>
  <desc id="gap-desc">A scaffold has an N gap between left and right contigs. A graph has two possible paths, but read support selects one path that replaces the Ns.</desc>
  <defs>
    <marker id="gap-arrow" markerWidth="9" markerHeight="9" refX="7" refY="3" orient="auto" markerUnits="strokeWidth">
      <path d="M0,0 L0,6 L8,3 z" fill="#11443b" />
    </marker>
  </defs>
  <rect x="20" y="20" width="880" height="390" rx="12" fill="#fbfbf8" stroke="#deded8" />
  <text x="52" y="62" font-size="22" font-weight="700" fill="#202124">Gap filling replaces uncertainty with reviewed sequence</text>
  <text x="52" y="94" font-size="15" fill="#5f6368">A fill is safe only when one candidate path connects the correct flanks and alternatives are ruled out.</text>

  <text x="72" y="143" font-size="16" font-weight="700" fill="#202124">Scaffold gap</text>
  <rect x="112" y="164" width="190" height="34" rx="8" fill="#206a5d" />
  <text x="158" y="186" font-size="13" fill="#ffffff">left contig</text>
  <rect x="302" y="168" width="110" height="26" rx="4" fill="#f3f5f4" stroke="#b8c3bf" stroke-dasharray="4 4" />
  <text x="343" y="186" font-size="13" fill="#5f6368">Ns</text>
  <rect x="412" y="164" width="190" height="34" rx="8" fill="#b66b33" />
  <text x="456" y="186" font-size="13" fill="#ffffff">right contig</text>

  <path d="M357 205 L357 232" stroke="#11443b" stroke-width="3" marker-end="url(#gap-arrow)" />

  <text x="72" y="263" font-size="16" font-weight="700" fill="#202124">Candidate graph paths</text>
  <circle cx="145" cy="316" r="22" fill="#206a5d" />
  <text x="132" y="321" font-size="13" fill="#ffffff">L</text>
  <circle cx="300" cy="286" r="22" fill="#6b8f71" />
  <text x="287" y="291" font-size="13" fill="#ffffff">A</text>
  <circle cx="455" cy="286" r="22" fill="#6b8f71" />
  <text x="442" y="291" font-size="13" fill="#ffffff">B</text>
  <circle cx="610" cy="316" r="22" fill="#b66b33" />
  <text x="597" y="321" font-size="13" fill="#ffffff">R</text>
  <line x1="167" y1="312" x2="278" y2="290" stroke="#6b8f71" stroke-width="5" />
  <line x1="322" y1="286" x2="433" y2="286" stroke="#6b8f71" stroke-width="5" />
  <line x1="477" y1="290" x2="588" y2="312" stroke="#6b8f71" stroke-width="5" />
  <text x="288" y="256" font-size="13" fill="#11443b">supported path</text>

  <circle cx="300" cy="356" r="18" fill="#d9d9d4" />
  <text x="290" y="361" font-size="12" fill="#5f6368">X</text>
  <circle cx="455" cy="356" r="18" fill="#d9d9d4" />
  <text x="445" y="361" font-size="12" fill="#5f6368">Y</text>
  <line x1="164" y1="322" x2="282" y2="350" stroke="#b8c3bf" stroke-width="3" stroke-dasharray="6 5" />
  <line x1="318" y1="356" x2="437" y2="356" stroke="#b8c3bf" stroke-width="3" stroke-dasharray="6 5" />
  <line x1="473" y1="350" x2="591" y2="322" stroke="#b8c3bf" stroke-width="3" stroke-dasharray="6 5" />
  <text x="286" y="389" font-size="13" fill="#5f6368">alternative path left unresolved</text>

  <path d="M650 316 L720 316" stroke="#11443b" stroke-width="3" marker-end="url(#gap-arrow)" />
  <rect x="736" y="282" width="118" height="68" rx="8" fill="#e8f1ee" stroke="#cfe0da" />
  <text x="756" y="307" font-size="13" fill="#11443b">replace Ns</text>
  <text x="756" y="329" font-size="13" fill="#11443b">after review</text>
</svg>
<figcaption><strong>Figure 1. Gap filling is path selection plus sequence validation.</strong> The graph may contain more than one path. Filling is appropriate only when one path is selected and its sequence fits the scaffold flanks.</figcaption>
</figure>

## Why Gap Filling Is Harder Than Scaffolding

Scaffolding can be useful even when the gap sequence is unknown. Gap filling
must choose actual bases. That means the evidence has to answer more questions:

- Do the scaffold flanks map to graph nodes or path ends?
- Is there a path between those flanks?
- Is there exactly one best path, or several plausible paths?
- Does the graph path have sequence?
- Do the path ends match the scaffold FASTA flanks?
- Do reads, graph topology, Hi-C contacts, or reference placement support the
  same path?

If those questions do not converge, the honest result is usually to leave Ns.

## Why We Think It Works

Assembly graphs preserve sequence and adjacency possibilities that may not be
present in a contig-only FASTA. If the graph has a unique sequenced path
between two scaffold flanks, and independent evidence supports that path, the
path can explain the missing interval better than an arbitrary N block.

The method works best when the graph and scaffold come from compatible assembly
stages. If contig names, coordinates, or sequences have drifted, a fill can
appear plausible but belong to the wrong version of the assembly.

## When To Fill

Fill a gap only when the candidate path is:

- connected to both flanks,
- unique or uniquely supported,
- sequenced,
- validated against the flank sequences,
- accepted as part of a reviewed plan.

Leave the gap unresolved when candidate paths tie, evidence sources conflict,
nodes lack sequence, flanks do not validate, or the graph has no path within a
reasonable search depth.

## Example

A scaffold has a 20 kb N gap between contig A and contig B. In the assembly
graph, A and B connect through two possible paths. Long read alignments to the
graph traverse only path 1, and the path sequence matches the suffix of A and
the prefix of B. Path 2 has no read support and includes an unsequenced node.

The concept-level interpretation is:

1. There is a candidate sequence between the correct flanks.
2. One path is better supported than the alternative.
3. The sequence validates against the scaffold boundaries.
4. The fill can be accepted if the review goal is to replace Ns with graph
   sequence.

The ChromoSort-specific next steps are
[Graph-supported gap filling]({{ '/guides/graph-gap-filling/' | relative_url }}),
[Assembly graph evidence]({{ '/guides/assembly-graph-evidence/' | relative_url }}),
and [chromo gapfill]({{ '/commands/gapfill/' | relative_url }}).
