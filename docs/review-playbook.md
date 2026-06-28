---
title: Agent And Review Playbook
description: Reproducible ChromoSort review patterns for new datasets, crops, species, and assistant handoffs.
---

# Agent And Review Playbook

This playbook turns common ChromoSort review work into a portable procedure.
It is written for new project chats, coding agents, and users moving from one
crop, species, or assembly batch to another. The examples use generic file
names; replace reference, assembly, sample, and chromosome names with the names
from your dataset.

The goal is not to make every decision automatic. The goal is to keep the
inputs, evidence classes, and intervention choices explicit enough that another
person or agent can reproduce the same reasoning later.

## Start With A Dataset Manifest

Keep one small manifest per analysis batch. A plain TSV is enough:

```text
sample  reference_fasta  assembly_fasta  alignment_format  alignment_path  gfa  read_paf  gaf  notes
sampleA reference.fa assembly.fa paf paf/sampleA.paf graphs/sampleA.gfa reads/sampleA.reads_to_assembly.paf graph_reads/sampleA.gaf raw assembly
```

At minimum, record:

- the exact reference FASTA and assembly FASTA used by the coords or PAF file,
- the aligner command and preset,
- whether PAF was generated with `-c --secondary=no`,
- the ChromoSort command and output prefix,
- any FASTA-changing step that requires a fresh downstream alignment.

The most important invariant is the
[FASTA/alignment compatibility rule]({{ '/input-files/' | relative_url }}#fasta-and-alignment-compatibility):
do not use an old raw-assembly alignment as evidence for a changed FASTA.

## Choose One Primary Alignment Evidence

MUMmer coords and minimap2 PAF are alternative ways to provide the same class
of evidence: a whole-genome reference-to-assembly alignment. You usually do
not need both for a production run. Pick one primary alignment source, run
`sort`, `plot`, `fix`, and `eval fix` from that source, and use independent
evidence streams such as long-read PAF, GFA, and GAF when a biological decision
needs more support.

PAF is the recommended primary source for most new runs because it is much
faster in large plant-genome tests and supports MAPQ filtering. MUMmer coords
remains a good alternative, especially for projects with existing nucmer
pipelines or for a second aligner perspective on a surprising candidate.

Running both coords and PAF is still useful when benchmarking ChromoSort,
checking parser parity, or tuning minimap2/MUMmer settings for a new genome
group. Treat that comparison as a diagnostic, not as independent biological
validation.

As a loose expectation from the soybean coords-vs-PAF `chromo fix` benchmark,
split counts were close, differing by about 5-10%, while the exact set of
marginal split contigs differed by about 20-30%. Those differences appeared to
come from aligner behavior and output structure, such as row fragmentation,
secondary/primary handling, MAPQ, and identity fields, not from ChromoSort
applying different post-normalization logic to coords and PAF.

If you choose MUMmer, use a filtered reference-vs-assembly coords file:

```bash
nucmer -t 16 -c 500 -p mummer/${sample} reference.fa assembly.fa
delta-filter -i 95 -l 10000 -1 mummer/${sample}.delta > mummer/${sample}.filter
show-coords -r -c -l mummer/${sample}.filter > mummer/${sample}.coords
```

If you choose minimap2 PAF, keep base-level PAF output and suppress secondary
rows:

```bash
minimap2 -x asm5 -c -t 16 --secondary=no reference.fa assembly.fa \
  > paf/${sample}.paf
```

Choose the strictest minimap2 preset that still recovers the expected
chromosome-scale alignments. `asm5` is the safest same-species default. Move to
`asm10` or `asm20` only when expected syntenic blocks are missing or badly
fragmented. When using permissive presets, inspect plots and consider
`--min-mapq`; avoid PAF identity filters until you have checked the PAF
identity distribution.

## Run Sort And Plot From The Primary Alignment

For a PAF-backed run:

```bash
chromo sort \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.fa \
  --paf paf/${sample}.paf \
  --output-prefix results/${sample} \
  --orient-to-reference

chromo plot \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.fa \
  --paf paf/${sample}.paf \
  --assignments results/${sample}.contig_assignments.tsv \
  --output-prefix plots/${sample} \
  --per-ref
```

For a coords-backed run, replace `--paf paf/${sample}.paf` with
`--coords mummer/${sample}.coords`.

For a focused replot, add `--sel-ref` with the reference IDs involved in the
decision.

## Review Primary Alignment Decisions

Start from three reports per sample:

- `<prefix>.run_summary.txt`: settings, input paths, status counts, and PAF
  diagnostics when PAF is the primary alignment.
- `<prefix>.contig_assignments.tsv`: one row per contig with assignment,
  status, second reference, and split-candidate flags.
- `<prefix>.contig_ref_matches.tsv`: per-contig, per-reference evidence.

Use these decision classes:

| Pattern | Interpretation | Next action |
| --- | --- | --- |
| Contig is `kept_split_candidate` with two substantial reference matches | Review target | Inspect the plot and match rows; fix only reviewed contigs with `chromo fix --mode conservative`; re-align the fixed FASTA. |
| One strong same-reference placement plus tiny off-target matches | Usually repeat or paralog noise | Leave unsplit unless other evidence supports a real event. |
| Whole contig is reverse relative to the reference | Orientation issue | Use `chromo sort --orient-to-reference`; do not split. |
| Same-reference internal inversion | Biological inversion, reference difference, or assembly error | Evaluate with read and graph evidence; do not automatically reference-normalize. |
| Terminal overlap between separate contigs | Scaffold/overlap problem, not a within-contig fix | Review sort/scaffold overlap reports and use explicit scaffold overlap policies only when justified. |

If you also ran the other whole-genome alignment format, use it as a diagnostic
cross-check:

| Optional coords-vs-PAF pattern | Interpretation |
| --- | --- |
| Both inputs flag the same contig and references | The alignment-format paths are coherent for that event; still use plot/read/graph evidence for the biological call. |
| Candidate appears only in one input | Aligner-specific evidence or threshold effect; inspect plots and eval evidence before fixing. |
| Dominant assigned reference differs, but the same references are involved | Usually a near-tie assignment difference; decide on the event, not the winner label alone. |

For production fixes, prefer targeted contig lists over `--all`:

```bash
chromo fix \
  --assembly-fasta assembly.fa \
  --paf paf/${sample}.paf \
  --contigs contig_a contig_b \
  --mode conservative \
  --min-mapq 20 \
  --orient-to-reference \
  --output-fasta results/fix/${sample}.fixed.fa \
  --report results/fix/${sample}.fix_report.tsv
```

If you use coords instead of PAF, replace `--paf` with `--coords` and drop the
PAF-specific `--min-mapq` filter.

## Review Same-Reference Inversions

An internal inversion is different from a multi-reference chimeric contig. A
dot plot may show one contig assigned cleanly to one reference, but with a
large internal block in the opposite orientation. That pattern can be a real
sample allele, a difference between the sample and reference, or an assembly
error.

Ask two separate questions:

1. Is the inversion real in this assembly or haplotype?
2. If it is real, should the delivered FASTA preserve it or be deliberately
   reference-normalized?

For pangenome graph inputs, preserve real haplotype structure. Whole-contig
orientation to the reference is often useful, but flipping a validated internal
inversion removes the allele from that assembly path. Reference-normalized
pseudoassemblies can be useful for a specific downstream comparison, but they
should be labeled separately from biological assembly inputs.

Use `chromo eval fix --mode comprehensive` to make a review table without
immediately applying a fix:

```bash
chromo eval fix \
  --assembly-fasta assembly.fa \
  --paf paf/${sample}.paf \
  --contigs contig_with_inversion \
  --mode comprehensive \
  --min-mapq 20 \
  --gfa graphs/${sample}.gfa \
  --gaf graph_reads/${sample}.gaf \
  --read-paf reads/${sample}.reads_to_assembly.paf \
  --read-window-bp 10000 \
  --read-min-anchor-bp 2000 \
  --orient-to-reference \
  --output-prefix review/${sample}.inversion
```

Evidence that makes a same-reference inversion more believable includes:

- the primary alignment showing a sharp orientation switch with strong blocks,
- an optional second whole-genome alignment reproducing the same breakpoints,
- long alignments and high unique support on both sides of the event,
- long reads spanning the candidate breakpoints in the read-to-assembly PAF,
- graph paths or GFA topology consistent with the assembly sequence,
- confirmation in another assembly, map, or related accession when available.

Evidence that argues for caution includes low MAPQ rows, repeat-rich or
homeologous placements, many small off-target hits, breakpoint pileups of
edge-only reads without spanning reads, or disagreement between graph and
read-to-assembly evidence.

## Add Long-Read And Graph Evidence

For read-to-assembly evidence:

```bash
minimap2 -x map-hifi -c -t 16 --secondary=no assembly.fa reads.fastq.gz \
  > reads/${sample}.reads_to_assembly.paf
```

For read-to-graph evidence, use GraphAligner or another graph mapper that writes
GAF:

```bash
GraphAligner \
  -g graphs/${sample}.gfa \
  -f reads.fastq.gz \
  -a graph_reads/${sample}.gaf
```

When full-depth read-to-graph alignment is too expensive, prepare one targeted
GraphAligner input bundle for all review paths:

```bash
chromo eval all \
  --assembly-fasta results/${sample}.ordered.fa \
  --coords mummer/${sample}.ordered.coords \
  --all \
  --ordered-fasta results/${sample}.ordered.fa \
  --assignments results/${sample}.contig_assignments.tsv \
  --gfa graphs/${sample}.gfa \
  --read-paf reads/${sample}.reads_to_ordered.paf \
  --output-prefix review/${sample}.eval

chromo gafprep \
  --assembly-fasta results/${sample}.ordered.fa \
  --assembly-gfa graphs/${sample}.gfa \
  --read-paf reads/${sample}.reads_to_ordered.paf \
  --reads reads.fastq.gz \
  --eval-review-table review/${sample}.eval.fix_review.tsv \
  --eval-review-table review/${sample}.eval.scaffold_review.tsv \
  --eval-review-table review/${sample}.eval.gapfill_review.tsv \
  --output-prefix graph_reads/${sample}.gafprep

bash graph_reads/${sample}.gafprep.graphaligner.sh
```

The review tables and read-to-assembly PAF should name the same FASTA records;
for a one-run GAF bundle, that usually means running `eval all` and minimap2
against the ordered FASTA used for scaffold and gapfill review.

`chromo eval fix` summarizes read support near proposed breakpoints as
`longread_spanning_reads`, `longread_split_reads`, edge-read counts, and nearby
read counts. GFA and GAF fields are advisory context for fix review. They do
not change sequence by themselves.

When using hifiasm unitig graphs with contig FASTA dot plots, first confirm
that the GFA can bridge the coordinate systems. `p_ctg.fa` and
`hap*.p_ctg.fa` dot plots are in contig coordinates; `p_utg.gfa`, `r_utg.gfa`,
and `.noseq.gfa` files are often in unitig coordinates. Use `chromo graph-map`
to verify that the GFA contains `P` path or `W` walk records matching the
contig FASTA names before treating unitig intervals as contig-axis evidence.
When projection is available, `chromo eval fix` reports the containing graph
unitig and nearest projected unitig boundaries around proposed breakpoints, and
`chromo plot --gfa-overlay` can draw those projected unitig spans on the query
axis.

For scaffold and gapfill review, GFA and GAF become more central because the
question is often whether graph paths connect adjacent contigs or explain a
gap. Use `chromo eval scaffold` and `chromo eval gapfill` when the decision is
between contig-end joins or candidate graph paths rather than within-contig
breakpoints.

## Decide What To Deliver

Use explicit delivery labels:

| Delivery class | Use when | Typical command path |
| --- | --- | --- |
| Native assembly | Structure appears real or unresolved | Sort, orient whole contigs if desired, scaffold conservatively. |
| Fixed assembly | Evidence supports an assembly misjoin or wrong breakpoint | Targeted `chromo fix`, `chromo cut`, or reviewed-plan application, then re-align. |
| Review-only table | Evidence is interesting but not ready for sequence editing | `chromo eval fix/scaffold/gapfill` plus focused plots or manual dashboard. |
| Reference-normalized experimental FASTA | Downstream analysis intentionally requires reference orientation at local events | Run explicit comprehensive/sensitive fix or manual edits into a clearly named output folder. |
| Pangenome graph input | Building a graph intended to represent real alleles | Preserve validated inversions and other real structure; avoid unlabelled reference normalization. |

Document the choice in the manifest or a short decision table. Future agents
should be able to answer: what was changed, what was left native, which evidence
supported that choice, and which FASTA each downstream alignment describes.

## Handoff Checklist For A New Chat Or Agent

Before handing off a review, provide:

- the manifest or a list of exact input paths,
- commands used to generate the primary MUMmer coords or minimap2 PAF,
- ChromoSort sort and plot commands,
- the relevant `run_summary`, `contig_assignments`, and `contig_ref_matches`
  files,
- focused plots for any contested contig or reference,
- any GFA, read-to-assembly PAF, and read-to-graph GAF evidence,
- a decision table with rows such as `fix now`, `review before fix`, `leave
  native`, and `reference-normalized experiment`,
- notes about any dataset-specific preset choices such as `asm20` for a
  divergent reference.

Avoid crop-specific assumptions in the handoff. The same evidence logic works
with soybean chromosomes, grass chromosomes, fungal scaffolds, animal
chromosomes, or any other reference/assembly naming convention.
