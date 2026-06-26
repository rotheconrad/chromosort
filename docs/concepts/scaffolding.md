---
title: How Scaffolding Works
description: Conceptual lesson for ordering contigs, orienting them, and representing unknown gaps with scaffold Ns.
---

# How Scaffolding Works

Scaffolding turns ordered contigs into larger records. In a chromosome-scale
assembly workflow, the scaffold is often the file a researcher wants to inspect,
compare, or submit: one record per chromosome, linkage group, or reference
sequence.

The key point is that scaffolding joins existing sequence. It does not discover
the missing bases between contigs. When the sequence between adjacent contigs is
unknown, the scaffold uses Ns as an explicit placeholder.

## The Core Idea

Scaffolding uses placement evidence to decide three things for adjacent
contigs:

1. order,
2. orientation,
3. gap or overlap representation.

<figure class="concept-figure">
<svg viewBox="0 0 920 420" role="img" aria-labelledby="scaffold-core-title scaffold-core-desc">
  <title id="scaffold-core-title">Ordered contigs joined into a scaffold with N gaps</title>
  <desc id="scaffold-core-desc">Three contigs are placed along a reference, oriented, and joined into one scaffold with N gaps where sequence is unknown.</desc>
  <rect x="20" y="20" width="880" height="380" rx="12" fill="#fbfbf8" stroke="#deded8" />
  <text x="52" y="62" font-size="22" font-weight="700" fill="#202124">Scaffolding joins placed contigs with visible uncertainty</text>
  <text x="52" y="94" font-size="15" fill="#5f6368">The output is a larger coordinate frame, not proof that every base between contigs is known.</text>

  <text x="72" y="145" font-size="16" font-weight="700" fill="#202124">Placement evidence</text>
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

  <path d="M430 225 L430 260" stroke="#11443b" stroke-width="3" />
  <polygon points="430,265 423,253 437,253" fill="#11443b" />

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

The scaffold is a useful coordinate system. It should also be honest about what
is known and unknown.

## Why We Think Scaffolding Works

Scaffolding works when independent evidence agrees on the relative placement of
contigs. A close reference can provide approximate order and orientation.
Paired reads, Hi-C contacts, long reads, optical maps, linkage maps, and
assembly graph links can add support when they are available.

<figure class="concept-figure">
<svg viewBox="0 0 920 430" role="img" aria-labelledby="scaffold-evidence-title scaffold-evidence-desc">
  <title id="scaffold-evidence-title">Evidence types used to support scaffolding</title>
  <desc id="scaffold-evidence-desc">Reference placement, Hi-C contacts, long reads, and graph links all support the same contig adjacency.</desc>
  <rect x="20" y="20" width="880" height="390" rx="12" fill="#fbfbf8" stroke="#deded8" />
  <text x="52" y="62" font-size="22" font-weight="700" fill="#202124">Good scaffolds are supported by converging evidence</text>
  <text x="52" y="94" font-size="15" fill="#5f6368">Different data types can support the same adjacency for different reasons.</text>

  <rect x="110" y="158" width="210" height="36" rx="8" fill="#206a5d" />
  <text x="176" y="181" font-size="14" fill="#ffffff">contig A</text>
  <rect x="602" y="158" width="210" height="36" rx="8" fill="#b66b33" />
  <text x="668" y="181" font-size="14" fill="#ffffff">contig B</text>
  <line x1="320" y1="176" x2="602" y2="176" stroke="#5f6368" stroke-width="2" stroke-dasharray="5 5" />
  <text x="424" y="207" font-size="13" fill="#5f6368">candidate junction</text>

  <rect x="95" y="250" width="160" height="54" rx="8" fill="#eef6f2" stroke="#cfe0da" />
  <text x="126" y="272" font-size="13" font-weight="700" fill="#11443b">Reference</text>
  <text x="116" y="291" font-size="12" fill="#5f6368">same chromosome order</text>
  <path d="M255 270 C335 230, 390 220, 448 184" fill="none" stroke="#206a5d" stroke-width="3" />

  <rect x="280" y="298" width="160" height="54" rx="8" fill="#eef6f2" stroke="#cfe0da" />
  <text x="319" y="320" font-size="13" font-weight="700" fill="#11443b">Hi-C</text>
  <text x="306" y="339" font-size="12" fill="#5f6368">contact enrichment</text>
  <path d="M440 316 C478 286, 512 250, 556 184" fill="none" stroke="#4f6d8a" stroke-width="3" />

  <rect x="480" y="298" width="160" height="54" rx="8" fill="#eef6f2" stroke="#cfe0da" />
  <text x="523" y="320" font-size="13" font-weight="700" fill="#11443b">Reads</text>
  <text x="506" y="339" font-size="12" fill="#5f6368">bridges or links</text>
  <path d="M480 316 C442 286, 408 250, 364 184" fill="none" stroke="#b66b33" stroke-width="3" />

  <rect x="665" y="250" width="160" height="54" rx="8" fill="#eef6f2" stroke="#cfe0da" />
  <text x="712" y="272" font-size="13" font-weight="700" fill="#11443b">Graph</text>
  <text x="696" y="291" font-size="12" fill="#5f6368">oriented link/path</text>
  <path d="M665 270 C585 230, 530 220, 472 184" fill="none" stroke="#6b8f71" stroke-width="3" />
</svg>
<figcaption><strong>Figure 2. Evidence can converge on an adjacency.</strong> A scaffold is more credible when multiple independent signals point to the same order and orientation.</figcaption>
</figure>

No single evidence type is perfect. References can be wrong or biologically
different. Hi-C can join chromosome arms but is lower resolution near repeats.
Long reads and graph paths can be missing in repetitive regions. Good
scaffolding treats evidence as support, not magic.

## What Scaffolding Changes

Scaffolding can:

- order contigs,
- orient contigs,
- join contigs into larger FASTA records,
- insert N gaps for unknown sequence,
- report inferred gaps and overlaps,
- create provenance records such as AGP or scaffold-gap tables.

Scaffolding should not:

- pretend Ns are known bases,
- silently remove real sequence,
- resolve an ambiguous graph branch,
- fix a chimeric contig that should have been reviewed earlier.

## Gaps And Overlaps

Adjacent contigs can have a positive gap, no inferred gap, or a negative gap
that indicates overlap in the placement coordinate system.

<figure class="concept-figure">
<svg viewBox="0 0 920 440" role="img" aria-labelledby="scaffold-gaps-title scaffold-gaps-desc">
  <title id="scaffold-gaps-title">Positive gaps, zero gaps, and overlaps in scaffolding</title>
  <desc id="scaffold-gaps-desc">Three rows show a positive gap, touching contigs, and an overlap that should be reviewed before trimming.</desc>
  <rect x="20" y="20" width="880" height="400" rx="12" fill="#fbfbf8" stroke="#deded8" />
  <text x="52" y="62" font-size="22" font-weight="700" fill="#202124">Junctions have different meanings</text>
  <text x="52" y="94" font-size="15" fill="#5f6368">A scaffold-gap report is how the uncertainty stays visible.</text>

  <g font-size="13">
    <text x="82" y="145" font-size="16" font-weight="700" fill="#202124">Positive gap</text>
    <rect x="230" y="125" width="160" height="30" rx="8" fill="#206a5d" />
    <rect x="490" y="125" width="160" height="30" rx="8" fill="#b66b33" />
    <line x1="390" y1="140" x2="490" y2="140" stroke="#5f6368" stroke-width="2" stroke-dasharray="5 5" />
    <text x="414" y="174" fill="#5f6368">insert Ns</text>

    <text x="82" y="235" font-size="16" font-weight="700" fill="#202124">Touching spans</text>
    <rect x="230" y="215" width="210" height="30" rx="8" fill="#206a5d" />
    <rect x="440" y="215" width="210" height="30" rx="8" fill="#b66b33" />
    <text x="401" y="264" fill="#5f6368">zero inferred gap</text>

    <text x="82" y="325" font-size="16" font-weight="700" fill="#202124">Negative gap</text>
    <rect x="230" y="305" width="255" height="30" rx="8" fill="#206a5d" />
    <rect x="430" y="305" width="220" height="30" rx="8" fill="#b66b33" opacity="0.85" />
    <rect x="430" y="300" width="55" height="40" fill="none" stroke="#202124" stroke-width="2" stroke-dasharray="5 5" />
    <text x="398" y="354" fill="#5f6368">overlap: review before trimming</text>
  </g>
</svg>
<figcaption><strong>Figure 3. Gaps and overlaps are junction annotations.</strong> Positive gaps usually become Ns. Overlaps require a policy and should be reported clearly.</figcaption>
</figure>

This is why the scaffold FASTA alone is not enough. The FASTA contains
sequence and Ns. It does not explain why a gap length was chosen, whether an
overlap was trimmed, or whether a gap was manually reviewed.

## Provenance Matters

A scaffold is easiest to trust when it has a map from scaffold coordinates back
to the original components.

<figure class="concept-figure">
<svg viewBox="0 0 920 410" role="img" aria-labelledby="scaffold-provenance-title scaffold-provenance-desc">
  <title id="scaffold-provenance-title">Scaffold provenance with component and gap records</title>
  <desc id="scaffold-provenance-desc">A scaffold FASTA is paired with AGP-like component and gap rows that describe contigs and N intervals.</desc>
  <rect x="20" y="20" width="880" height="370" rx="12" fill="#fbfbf8" stroke="#deded8" />
  <text x="52" y="62" font-size="22" font-weight="700" fill="#202124">A scaffold should keep a component map</text>

  <text x="72" y="118" font-size="16" font-weight="700" fill="#202124">Scaffold FASTA</text>
  <rect x="72" y="140" width="170" height="30" rx="8" fill="#206a5d" />
  <rect x="242" y="144" width="60" height="22" rx="4" fill="#f3f5f4" stroke="#b8c3bf" stroke-dasharray="4 4" />
  <rect x="302" y="140" width="150" height="30" rx="8" fill="#4f6d8a" />
  <rect x="452" y="144" width="80" height="22" rx="4" fill="#f3f5f4" stroke="#b8c3bf" stroke-dasharray="4 4" />
  <rect x="532" y="140" width="165" height="30" rx="8" fill="#b66b33" />
  <text x="252" y="161" font-size="12" fill="#5f6368">Ns</text>
  <text x="480" y="161" font-size="12" fill="#5f6368">Ns</text>

  <text x="72" y="220" font-size="16" font-weight="700" fill="#202124">Map records</text>
  <rect x="72" y="240" width="760" height="112" rx="8" fill="#ffffff" stroke="#deded8" />
  <text x="94" y="267" font-size="13" fill="#202124">1-170: contig A, forward</text>
  <text x="94" y="292" font-size="13" fill="#202124">171-230: gap, 60 Ns, inferred</text>
  <text x="94" y="317" font-size="13" fill="#202124">231-380: contig B, reverse</text>
  <text x="430" y="267" font-size="13" fill="#202124">381-460: gap, 80 Ns, reviewed</text>
  <text x="430" y="292" font-size="13" fill="#202124">461-625: contig C, forward</text>
  <text x="430" y="317" font-size="13" fill="#5f6368">Every scaffold base has a source or a gap reason.</text>
</svg>
<figcaption><strong>Figure 4. Provenance makes scaffolds reviewable.</strong> Component maps and gap reports explain where each scaffold interval came from.</figcaption>
</figure>

AGP files, scaffold-gap reports, and run summaries are not administrative
afterthoughts. They are what let another person reproduce, inspect, and revise
the scaffold later.

## A Bad Scaffold Can Look Convincing

Wrong order or orientation can create a chromosome-scale record that looks
tidy in FASTA but creates problems in alignments, variant calling, and
downstream annotation.

<figure class="concept-figure">
<svg viewBox="0 0 920 390" role="img" aria-labelledby="scaffold-bad-title scaffold-bad-desc">
  <title id="scaffold-bad-title">Correct and incorrect scaffold order</title>
  <desc id="scaffold-bad-desc">One scaffold order follows the reference A B C while another places C between A and B and creates a false long-range rearrangement.</desc>
  <rect x="20" y="20" width="880" height="350" rx="12" fill="#fbfbf8" stroke="#deded8" />
  <text x="52" y="62" font-size="22" font-weight="700" fill="#202124">Scaffolding errors become coordinate errors</text>

  <text x="86" y="122" font-size="16" font-weight="700" fill="#202124">Supported order</text>
  <rect x="86" y="148" width="150" height="32" rx="8" fill="#206a5d" />
  <rect x="250" y="148" width="150" height="32" rx="8" fill="#4f6d8a" />
  <rect x="414" y="148" width="150" height="32" rx="8" fill="#b66b33" />
  <text x="146" y="169" font-size="13" fill="#ffffff">A</text>
  <text x="310" y="169" font-size="13" fill="#ffffff">B</text>
  <text x="474" y="169" font-size="13" fill="#ffffff">C</text>
  <text x="600" y="169" font-size="13" fill="#11443b">consistent downstream coordinates</text>

  <text x="86" y="240" font-size="16" font-weight="700" fill="#202124">Unsupported order</text>
  <rect x="86" y="266" width="150" height="32" rx="8" fill="#206a5d" />
  <rect x="250" y="266" width="150" height="32" rx="8" fill="#b66b33" />
  <rect x="414" y="266" width="150" height="32" rx="8" fill="#4f6d8a" />
  <text x="146" y="287" font-size="13" fill="#ffffff">A</text>
  <text x="310" y="287" font-size="13" fill="#ffffff">C</text>
  <text x="474" y="287" font-size="13" fill="#ffffff">B</text>
  <path d="M610 282 C660 240, 720 240, 770 282" fill="none" stroke="#a44747" stroke-width="3" />
  <text x="600" y="320" font-size="13" fill="#a44747">false rearrangement signal</text>
</svg>
<figcaption><strong>Figure 5. A scaffold is a hypothesis about coordinate order.</strong> Wrong order can create artificial structural signals downstream.</figcaption>
</figure>

This is why high-confidence scaffolding is usually late in the workflow. Fix
bad contigs first, sort and orient contigs carefully, then scaffold the
reviewed set.

## Example Walkthrough

Three contigs align to the same chromosome in order: A, B, C. The right end of
A maps before the left end of B, leaving an inferred 12 kb reference gap. B and
C overlap by 400 bp in reference coordinates, but their sequences do not match
well at the overlap.

The concept-level scaffold is:

1. Join A and B with an N gap because the missing sequence is unknown.
2. Record the inferred 12 kb gap in a gap report.
3. Treat the B/C negative gap as an overlap decision, not as permission to
   delete sequence automatically.
4. Preserve a component map so the scaffold can be audited.
5. Re-align the scaffold if scaffold-level validation is needed.

## Common Traps

Do not confuse scaffolding with gap filling. Scaffolding can place Ns. Gap
filling replaces some Ns with sequence.

Do not trim overlaps just because reference coordinates overlap. Sequence
confirmation and the biological context matter.

Do not scaffold before resolving obvious chimeric contigs. A scaffold cannot
make a bad component safe.

Do not forget that reference-guided scaffolding inherits reference assumptions.
If the sample has real structural differences, the reference may be a guide
rather than a truth source.

## Brief History And Further Reading

Scaffolding has been part of genome assembly since whole-genome shotgun
projects needed to connect contigs using mate pairs and clone-end information.
As sequencing changed, scaffold evidence expanded: long reads can span repeats
and gaps, Hi-C can provide chromosome-scale contact patterns, optical maps and
linked reads can add long-range constraints, and reference-guided methods can
use conserved synteny when an appropriate reference exists.

Modern chromosome-scale assemblies often combine these ideas. The conceptual
lesson has stayed stable: a scaffold is a supported ordering of known sequence
plus explicit uncertainty.

Further reading:

- Batzoglou et al. 2002. [ARACHNE: a whole-genome shotgun assembler](https://doi.org/10.1101/gr.208902).
- Burton et al. 2013. [Chromosome-scale scaffolding of de novo genome assemblies based on chromatin interactions](https://doi.org/10.1038/nbt.2727).
- Dudchenko et al. 2017. [De novo assembly of the Aedes aegypti genome using Hi-C yields chromosome-length scaffolds](https://doi.org/10.1126/science.aal3327).
- Alonge et al. 2019. [RaGOO: fast and accurate reference-guided scaffolding of draft genomes](https://doi.org/10.1186/s13059-019-1829-6).
- NCBI. [AGP specification](https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/).

## What To Read Next

- [How gap filling works]({{ '/concepts/gap-filling/' | relative_url }})
- [Scaffolding, gaps, and overlaps]({{ '/guides/scaffolding-gaps-overlaps/' | relative_url }})
- [Graph-supported gap filling]({{ '/guides/graph-gap-filling/' | relative_url }})
- [chromo scaffold]({{ '/commands/scaffold/' | relative_url }})
