---
title: How Gap Filling Works
description: Conceptual lesson for replacing scaffold Ns with supported sequence from graph or read evidence.
---

# How Gap Filling Works

Gap filling tries to replace unknown scaffold sequence with actual bases. That
makes it more powerful than scaffolding, and more risky. A gap should be filled
only when the candidate sequence is connected to the correct flanks, supported
by evidence, and not confused with an equally plausible alternative.

The key question is:

> Do we know the sequence that belongs between these two scaffold flanks?

## The Core Idea

Scaffolding says "these contigs are adjacent, but the sequence between them is
unknown." Gap filling says "this candidate path supplies the missing sequence."

<figure class="concept-figure">
<svg viewBox="0 0 920 430" role="img" aria-labelledby="gap-core-title gap-core-desc">
  <title id="gap-core-title">Graph-supported gap filling replaces scaffold Ns</title>
  <desc id="gap-core-desc">A scaffold has an N gap between left and right contigs. A graph has two possible paths, but read support selects one path that replaces the Ns.</desc>
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

  <path d="M357 205 L357 232" stroke="#11443b" stroke-width="3" />
  <polygon points="357,237 350,225 364,225" fill="#11443b" />

  <text x="72" y="263" font-size="16" font-weight="700" fill="#202124">Candidate graph paths</text>
  <circle cx="145" cy="316" r="22" fill="#206a5d" />
  <text x="139" y="321" font-size="13" fill="#ffffff">L</text>
  <circle cx="300" cy="286" r="22" fill="#6b8f71" />
  <text x="294" y="291" font-size="13" fill="#ffffff">A</text>
  <circle cx="455" cy="286" r="22" fill="#6b8f71" />
  <text x="449" y="291" font-size="13" fill="#ffffff">B</text>
  <circle cx="610" cy="316" r="22" fill="#b66b33" />
  <text x="604" y="321" font-size="13" fill="#ffffff">R</text>
  <line x1="167" y1="312" x2="278" y2="290" stroke="#6b8f71" stroke-width="5" />
  <line x1="322" y1="286" x2="433" y2="286" stroke="#6b8f71" stroke-width="5" />
  <line x1="477" y1="290" x2="588" y2="312" stroke="#6b8f71" stroke-width="5" />
  <text x="288" y="256" font-size="13" fill="#11443b">supported path</text>

  <circle cx="300" cy="356" r="18" fill="#d9d9d4" />
  <text x="296" y="361" font-size="12" fill="#5f6368">X</text>
  <circle cx="455" cy="356" r="18" fill="#d9d9d4" />
  <text x="451" y="361" font-size="12" fill="#5f6368">Y</text>
  <line x1="164" y1="322" x2="282" y2="350" stroke="#b8c3bf" stroke-width="3" stroke-dasharray="6 5" />
  <line x1="318" y1="356" x2="437" y2="356" stroke="#b8c3bf" stroke-width="3" stroke-dasharray="6 5" />
  <line x1="473" y1="350" x2="591" y2="322" stroke="#b8c3bf" stroke-width="3" stroke-dasharray="6 5" />
  <text x="286" y="389" font-size="13" fill="#5f6368">alternative path left unresolved</text>

  <path d="M650 316 L720 316" stroke="#11443b" stroke-width="3" />
  <polygon points="725,316 713,309 713,323" fill="#11443b" />
  <rect x="736" y="282" width="118" height="68" rx="8" fill="#e8f1ee" stroke="#cfe0da" />
  <text x="756" y="307" font-size="13" fill="#11443b">replace Ns</text>
  <text x="756" y="329" font-size="13" fill="#11443b">after review</text>
</svg>
<figcaption><strong>Figure 1. Gap filling is path selection plus sequence validation.</strong> The graph may contain more than one path. Filling is appropriate only when one path is selected and its sequence fits the scaffold flanks.</figcaption>
</figure>

The difference from scaffolding is simple but important: scaffolding can be
honest with Ns. Gap filling must choose bases.

## Why Gap Filling Is Harder Than Scaffolding

Scaffolding can be useful even when the gap sequence is unknown. Gap filling
must answer more questions:

- Do the scaffold flanks map to graph nodes or path ends?
- Is there a path between those flanks?
- Is there exactly one best path, or several plausible paths?
- Does the graph path have sequence?
- Do the path ends match the scaffold FASTA flanks?
- Do reads, graph topology, Hi-C contacts, or reference placement support the
  same path?

If those questions do not converge, the honest result is usually to leave Ns.

## Path Statuses

Most gap-filling decisions reduce to a small set of patterns.

<figure class="concept-figure">
<svg viewBox="0 0 920 450" role="img" aria-labelledby="gap-status-title gap-status-desc">
  <title id="gap-status-title">Common gap filling path statuses</title>
  <desc id="gap-status-desc">Three rows show no path, one unique path, and ambiguous branching paths between scaffold flanks.</desc>
  <rect x="20" y="20" width="880" height="410" rx="12" fill="#fbfbf8" stroke="#deded8" />
  <text x="52" y="62" font-size="22" font-weight="700" fill="#202124">Different graph situations need different actions</text>

  <g font-size="13">
    <text x="80" y="130" font-size="16" font-weight="700" fill="#202124">No path</text>
    <circle cx="260" cy="124" r="22" fill="#206a5d" />
    <text x="254" y="129" fill="#ffffff">L</text>
    <circle cx="520" cy="124" r="22" fill="#b66b33" />
    <text x="514" y="129" fill="#ffffff">R</text>
    <line x1="290" y1="124" x2="490" y2="124" stroke="#b8c3bf" stroke-width="2" stroke-dasharray="6 5" />
    <text x="610" y="130" fill="#5f6368">Leave Ns or inspect inputs.</text>

    <text x="80" y="235" font-size="16" font-weight="700" fill="#202124">Unique path</text>
    <circle cx="260" cy="229" r="22" fill="#206a5d" />
    <text x="254" y="234" fill="#ffffff">L</text>
    <circle cx="390" cy="229" r="22" fill="#6b8f71" />
    <text x="384" y="234" fill="#ffffff">A</text>
    <circle cx="520" cy="229" r="22" fill="#b66b33" />
    <text x="514" y="234" fill="#ffffff">R</text>
    <line x1="282" y1="229" x2="368" y2="229" stroke="#6b8f71" stroke-width="5" />
    <line x1="412" y1="229" x2="498" y2="229" stroke="#6b8f71" stroke-width="5" />
    <text x="610" y="235" fill="#11443b">Fill after sequence validation and review.</text>

    <text x="80" y="340" font-size="16" font-weight="700" fill="#202124">Ambiguous paths</text>
    <circle cx="260" cy="334" r="22" fill="#206a5d" />
    <text x="254" y="339" fill="#ffffff">L</text>
    <circle cx="390" cy="304" r="19" fill="#6b8f71" />
    <text x="384" y="309" fill="#ffffff">A</text>
    <circle cx="390" cy="364" r="19" fill="#6b8f71" />
    <text x="384" y="369" fill="#ffffff">B</text>
    <circle cx="520" cy="334" r="22" fill="#b66b33" />
    <text x="514" y="339" fill="#ffffff">R</text>
    <line x1="280" y1="326" x2="372" y2="308" stroke="#6b8f71" stroke-width="4" />
    <line x1="408" y1="308" x2="502" y2="326" stroke="#6b8f71" stroke-width="4" />
    <line x1="280" y1="342" x2="372" y2="360" stroke="#6b8f71" stroke-width="4" />
    <line x1="408" y1="360" x2="502" y2="342" stroke="#6b8f71" stroke-width="4" />
    <text x="610" y="340" fill="#a44747">Do not guess. Add evidence or leave unresolved.</text>
  </g>
</svg>
<figcaption><strong>Figure 2. Path status controls the decision.</strong> A unique path is not the same as an ambiguous branch. No path and ambiguous paths should usually remain gaps.</figcaption>
</figure>

Ambiguity is not failure. In many genomes, especially repeat-rich plant
genomes, refusing to fill an ambiguous gap is the correct scientific choice.

## Flank Validation

A graph path is not enough. The path must also fit the scaffold sequence at
both ends.

<figure class="concept-figure">
<svg viewBox="0 0 920 420" role="img" aria-labelledby="gap-flanks-title gap-flanks-desc">
  <title id="gap-flanks-title">Gap fill flank validation</title>
  <desc id="gap-flanks-desc">A good graph path matches the left and right scaffold flanks, while a bad path has mismatched ends and should not be applied.</desc>
  <rect x="20" y="20" width="880" height="380" rx="12" fill="#fbfbf8" stroke="#deded8" />
  <text x="52" y="62" font-size="22" font-weight="700" fill="#202124">The path must fit the scaffold flanks</text>

  <text x="82" y="116" font-size="16" font-weight="700" fill="#202124">Pass</text>
  <rect x="82" y="140" width="190" height="28" rx="8" fill="#206a5d" />
  <text x="126" y="159" font-size="13" fill="#ffffff">left flank</text>
  <rect x="272" y="140" width="230" height="28" rx="8" fill="#6b8f71" />
  <text x="338" y="159" font-size="13" fill="#ffffff">graph fill sequence</text>
  <rect x="502" y="140" width="190" height="28" rx="8" fill="#b66b33" />
  <text x="548" y="159" font-size="13" fill="#ffffff">right flank</text>
  <path d="M245 186 C290 210, 340 210, 385 186" fill="none" stroke="#11443b" stroke-width="3" />
  <path d="M490 186 C535 210, 585 210, 630 186" fill="none" stroke="#11443b" stroke-width="3" />
  <text x="305" y="232" font-size="13" fill="#11443b">path ends match FASTA flanks</text>

  <text x="82" y="285" font-size="16" font-weight="700" fill="#202124">Fail</text>
  <rect x="82" y="309" width="190" height="28" rx="8" fill="#206a5d" />
  <rect x="272" y="309" width="230" height="28" rx="8" fill="#d9d9d4" />
  <rect x="502" y="309" width="190" height="28" rx="8" fill="#b66b33" />
  <line x1="252" y1="290" x2="292" y2="350" stroke="#a44747" stroke-width="4" />
  <line x1="292" y1="290" x2="252" y2="350" stroke="#a44747" stroke-width="4" />
  <line x1="482" y1="290" x2="522" y2="350" stroke="#a44747" stroke-width="4" />
  <line x1="522" y1="290" x2="482" y2="350" stroke="#a44747" stroke-width="4" />
  <text x="305" y="372" font-size="13" fill="#a44747">graph and scaffold likely come from different stages</text>
</svg>
<figcaption><strong>Figure 3. Flank validation protects against stale evidence.</strong> A path from the wrong graph or assembly stage can look plausible until its ends are compared with the actual scaffold sequence.</figcaption>
</figure>

If flank validation fails, stop. It usually means names, coordinates, sequence
versions, or graph stages drifted.

## Evidence Can Support Or Conflict

When several paths exist, evidence should select the same path. Conflicting
evidence is a warning, not a tie-breaker to ignore.

<figure class="concept-figure">
<svg viewBox="0 0 920 430" role="img" aria-labelledby="gap-evidence-title gap-evidence-desc">
  <title id="gap-evidence-title">Evidence support and conflict in gap filling</title>
  <desc id="gap-evidence-desc">GAF reads and reference placement support path A, while Hi-C weakly supports path B, creating a review decision.</desc>
  <rect x="20" y="20" width="880" height="390" rx="12" fill="#fbfbf8" stroke="#deded8" />
  <text x="52" y="62" font-size="22" font-weight="700" fill="#202124">Evidence should converge on one path</text>

  <circle cx="130" cy="214" r="22" fill="#206a5d" />
  <text x="124" y="219" font-size="13" fill="#ffffff">L</text>
  <circle cx="300" cy="154" r="21" fill="#6b8f71" />
  <text x="294" y="159" font-size="13" fill="#ffffff">A</text>
  <circle cx="300" cy="274" r="21" fill="#d9d9d4" />
  <text x="294" y="279" font-size="13" fill="#5f6368">B</text>
  <circle cx="470" cy="214" r="22" fill="#b66b33" />
  <text x="464" y="219" font-size="13" fill="#ffffff">R</text>
  <line x1="151" y1="206" x2="280" y2="162" stroke="#6b8f71" stroke-width="5" />
  <line x1="320" y1="162" x2="449" y2="206" stroke="#6b8f71" stroke-width="5" />
  <line x1="151" y1="222" x2="280" y2="266" stroke="#b8c3bf" stroke-width="4" stroke-dasharray="6 5" />
  <line x1="320" y1="266" x2="449" y2="222" stroke="#b8c3bf" stroke-width="4" stroke-dasharray="6 5" />

  <rect x="570" y="118" width="260" height="58" rx="8" fill="#eef6f2" stroke="#cfe0da" />
  <text x="590" y="142" font-size="13" font-weight="700" fill="#11443b">GAF read traversals</text>
  <text x="590" y="162" font-size="13" fill="#5f6368">support path A</text>
  <rect x="570" y="198" width="260" height="58" rx="8" fill="#eef6f2" stroke="#cfe0da" />
  <text x="590" y="222" font-size="13" font-weight="700" fill="#11443b">Reference placement</text>
  <text x="590" y="242" font-size="13" fill="#5f6368">supports path A</text>
  <rect x="570" y="278" width="260" height="58" rx="8" fill="#fff3ef" stroke="#e1b9aa" />
  <text x="590" y="302" font-size="13" font-weight="700" fill="#a44747">Hi-C contacts</text>
  <text x="590" y="322" font-size="13" fill="#5f6368">weakly favor path B</text>
</svg>
<figcaption><strong>Figure 4. Support can converge or conflict.</strong> A fill is easier to trust when reads, graph topology, and placement evidence point to the same path.</figcaption>
</figure>

The conservative choice is to fill only when support is unique enough for the
review goal. A benchmark run may apply all fillable paths to test behavior. A
production curation run should usually require explicit accepted rows.

## Possible Outcomes

Gap filling is not simply pass or fail.

<figure class="concept-figure">
<svg viewBox="0 0 920 390" role="img" aria-labelledby="gap-outcome-title gap-outcome-desc">
  <title id="gap-outcome-title">Gap filling outcomes</title>
  <desc id="gap-outcome-desc">Three outcomes show accepted fill, reviewed leave gap, and unresolved due to ambiguity.</desc>
  <rect x="20" y="20" width="880" height="350" rx="12" fill="#fbfbf8" stroke="#deded8" />
  <text x="52" y="62" font-size="22" font-weight="700" fill="#202124">A good result can still leave Ns</text>

  <rect x="80" y="118" width="220" height="160" rx="8" fill="#ffffff" stroke="#deded8" />
  <text x="112" y="150" font-size="16" font-weight="700" fill="#11443b">Accepted fill</text>
  <rect x="108" y="178" width="60" height="24" rx="6" fill="#206a5d" />
  <rect x="168" y="178" width="44" height="24" rx="6" fill="#6b8f71" />
  <rect x="212" y="178" width="60" height="24" rx="6" fill="#b66b33" />
  <text x="108" y="238" font-size="13" fill="#5f6368">Unique path, validated flanks, accepted review row.</text>

  <rect x="350" y="118" width="220" height="160" rx="8" fill="#ffffff" stroke="#deded8" />
  <text x="382" y="150" font-size="16" font-weight="700" fill="#202124">Reviewed gap</text>
  <rect x="378" y="178" width="60" height="24" rx="6" fill="#206a5d" />
  <rect x="438" y="180" width="44" height="20" rx="4" fill="#f3f5f4" stroke="#b8c3bf" stroke-dasharray="4 4" />
  <rect x="482" y="178" width="60" height="24" rx="6" fill="#b66b33" />
  <text x="378" y="238" font-size="13" fill="#5f6368">Fill exists, but reviewer keeps Ns for caution.</text>

  <rect x="620" y="118" width="220" height="160" rx="8" fill="#ffffff" stroke="#deded8" />
  <text x="652" y="150" font-size="16" font-weight="700" fill="#a44747">Unresolved</text>
  <line x1="676" y1="190" x2="774" y2="168" stroke="#b8c3bf" stroke-width="4" stroke-dasharray="6 5" />
  <line x1="676" y1="190" x2="774" y2="212" stroke="#b8c3bf" stroke-width="4" stroke-dasharray="6 5" />
  <text x="648" y="238" font-size="13" fill="#5f6368">No path, missing sequence, or ambiguous support.</text>
</svg>
<figcaption><strong>Figure 5. Leaving Ns can be the correct output.</strong> The goal is not to fill every gap. The goal is to fill only gaps whose sequence is defensible.</figcaption>
</figure>

This mindset is especially important in repeat-rich genomes. A confident
unresolved gap is better than a confident wrong sequence.

## Example Walkthrough

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

If either flank failed validation, or if read support tied between paths, the
better result would be to leave the Ns and report why.

## Common Traps

Do not fill gaps from a graph that no longer matches the scaffold FASTA stage.

Do not treat the shortest path as the correct path. Repeats and bubbles can
make the shortest graph path biologically wrong.

Do not override conflicting evidence just because a path can be reconstructed.
Reconstruction asks "can sequence be built?" Review asks "should this sequence
be placed here?"

Do not forget that unsequenced graph nodes cannot contribute bases. Topology
without sequence can support review, but it cannot fill the FASTA interval.

## Brief History And Further Reading

Gap filling began as an attempt to use additional reads or read pairs to close
N gaps left by short-read assemblies. As long reads became common, tools such
as PBJelly used long-read alignments to bridge gaps and improve draft genomes.
With modern graph-based assemblies, the same conceptual question appears in a
new form: which path through the graph, if any, belongs between these scaffold
flanks?

The history is a shift from "find any sequence that spans the gap" toward
"choose a supported sequence path and preserve uncertainty when the evidence is
ambiguous."

Further reading:

- Boetzer and Pirovano 2012. [Toward almost closed genomes with GapFiller](https://doi.org/10.1186/gb-2012-13-6-r56).
- English et al. 2012. [Mind the gap: upgrading genomes with Pacific Biosciences RS long-read sequencing technology](https://doi.org/10.1371/journal.pone.0047768).
- Walker et al. 2014. [Pilon: an integrated tool for comprehensive microbial variant detection and genome assembly improvement](https://doi.org/10.1371/journal.pone.0112963).
- Kolmogorov et al. 2019. [Assembly of long, error-prone reads using repeat graphs](https://doi.org/10.1038/s41587-019-0072-8).
- GFA working group. [Graphical Fragment Assembly format specification](https://github.com/GFA-spec/GFA-spec).

## What To Read Next

- [How scaffolding works]({{ '/concepts/scaffolding/' | relative_url }})
- [Graph-supported gap filling]({{ '/guides/graph-gap-filling/' | relative_url }})
- [Assembly graph evidence]({{ '/guides/assembly-graph-evidence/' | relative_url }})
- [chromo gapfill]({{ '/commands/gapfill/' | relative_url }})
