---
title: Workflows
description: Recommended ChromoSort workflows for sorting, fixing, scaffolding, and graph-aware review.
---

# Workflows

## Quick Start

```bash
git clone https://github.com/rotheconrad/chromosort.git
cd chromosort

mamba env create -f environment.yml
mamba activate chromosort

chromo --help
chromo sort --help
chromo fix --help
chromo cut --help
chromo manual --help
chromo plot --help
chromo scaffold --help
chromo gapfill --help
```

The sections below start after you have prepared a whole-genome alignment file
with
[MUMmer coords]({{ '/input-files/' | relative_url }}#creating-input-files-with-mummer)
or
[minimap2 PAF]({{ '/input-files/' | relative_url }}#creating-input-files-with-minimap2).
See [Input Files]({{ '/input-files/' | relative_url }}) for alignment commands
and graph-related inputs.

## Workflow 1: Reference-Order a Mostly Clean Assembly

Use this workflow when the assembly is already close to chromosome scale and
raw dot plots do not show obvious contig-scale misjoins. The goal is to place
contigs against the reference, write an ordered FASTA, inspect the placement,
and optionally build one scaffold record per reference sequence.

Inputs:

- `reference.fa`
- `assembly.fa`
- one alignment file, such as `mummer/sample.coords` or `paf/sample.paf`.
  If you still need to make it, see the Input Files recipes for
  [MUMmer coords]({{ '/input-files/' | relative_url }}#creating-input-files-with-mummer)
  or
  [minimap2 PAF]({{ '/input-files/' | relative_url }}#creating-input-files-with-minimap2).

Run the placement step:

```bash
chromo sort \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.fa \
  --coords mummer/sample.coords \
  --output-prefix results/sample \
  --orient-to-reference
```

The same workflow can use minimap2 PAF instead of MUMmer coords:

```bash
chromo sort \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.fa \
  --paf paf/sample.paf \
  --output-prefix results/sample \
  --orient-to-reference
```

Then plot the same alignment, using the assignment table to order the query
axis by kept ChromoSort contigs:

```bash
chromo plot \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.fa \
  --coords mummer/sample.coords \
  --assignments results/sample.contig_assignments.tsv \
  --output-prefix plots/sample \
  --per-ref
```

If the ordered contigs look reasonable, make chromosome-scale scaffold records:

```bash
chromo scaffold \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --output-prefix results/sample
```

Review these outputs before treating the result as final:

- `results/sample.contig_assignments.tsv`: placement status, overlap class, and
  kept/discarded decisions for each contig.
- `results/sample.match_report.tsv`: per-reference alignment support for each
  contig.
- `plots/sample.pdf` and per-reference plots: visual placement check.
- `results/sample.scaffold_gaps.tsv`: inferred gaps, overlaps, and scaffold
  overlap actions.

## Workflow 2: Fix Misjoined Contigs Before Sorting

Use this workflow when a raw dot plot shows contigs jumping between references,
orientation blocks, or otherwise looking chimeric. In this case, fix the raw
assembly first, then re-align the fixed FASTA and run `chromo sort` on the
updated assembly. Sorting protects strong split candidates, but it is still a
placement/filtering step, not a splitter.

Start with a raw plot so you can choose suspect contigs:

```bash
chromo plot \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.fa \
  --coords mummer/raw.coords \
  --output-prefix plots/sample.raw \
  --per-ref
```

Fix only reviewed contigs when you know which records are suspect:

```bash
chromo fix \
  --assembly-fasta assembly.fa \
  --coords mummer/raw.coords \
  --contigs suspect_contig_1 suspect_contig_2 \
  --output-fasta results/sample.fixed.fa \
  --report results/sample.fixed_contigs.tsv
```

If you want ChromoSort to scan all contigs with the conservative planner:

```bash
chromo fix \
  --assembly-fasta assembly.fa \
  --coords mummer/raw.coords \
  --all \
  --output-fasta results/sample.fixed.fa \
  --report results/sample.fixed_contigs.tsv
```

After writing `results/sample.fixed.fa`, re-run
[MUMmer]({{ '/input-files/' | relative_url }}#creating-input-files-with-mummer)
or
[minimap2]({{ '/input-files/' | relative_url }}#creating-input-files-with-minimap2)
against that fixed FASTA to create a fresh alignment, such as
`mummer/fixed.coords` or `paf/fixed.paf`. Then sort and plot the fixed assembly:

```bash
chromo sort \
  --ref-fasta reference.fa \
  --assembly-fasta results/sample.fixed.fa \
  --coords mummer/fixed.coords \
  --output-prefix results/sample.fixed \
  --orient-to-reference

chromo plot \
  --ref-fasta reference.fa \
  --assembly-fasta results/sample.fixed.fa \
  --coords mummer/fixed.coords \
  --assignments results/sample.fixed.contig_assignments.tsv \
  --output-prefix plots/sample.fixed \
  --per-ref
```

When the breakpoint is already known from manual review, use an explicit cut
instead of automatic split planning:

```bash
chromo cut \
  --assembly-fasta assembly.fa \
  --cut contig_1:234567,450000 \
  --cut contig_2:10000 \
  --output-fasta results/sample.cut.fa \
  --report results/sample.cut_contigs.tsv
```

For difficult cases, generate a browser review dashboard and apply the exported
recipe:

```bash
chromo manual \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.fa \
  --coords mummer/raw.coords \
  --gfa assembly_graph.gfa \
  --output-html results/sample.manual.html \
  --suggested-output-fasta sample.manual.fa

chromo manual apply \
  --assembly-fasta assembly.fa \
  --recipe chromosort.manual.recipe.json \
  --output-fasta results/sample.manual.fa \
  --report results/sample.manual.tsv
```

Review these outputs before moving downstream:

- `results/sample.fixed_contigs.tsv`: split status, emitted pieces, slice
  coordinates, dominant references, orientations, and reasons.
- `results/sample.fixed.fa`: the FASTA that must be re-aligned before sorting.
- `plots/sample.fixed.*`: confirmation that the repaired contigs now place as
  expected.

## Workflow 3: Scaffold and Fill Graph-Supported Gaps

Use this workflow after final sorting when you want one FASTA record per
reference sequence and you have graph evidence that may explain gaps between
adjacent sorted contigs. The default scaffold step reports graph context but
does not insert graph sequence. The gapfill step plans graph-supported fills
first, then applies only fillable paths, optionally after manual review.

First scaffold the final sorted contigs and write graph-junction evidence:

```bash
chromo scaffold \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --gfa assembly_graph.gfa \
  --output-prefix results/sample
```

Inspect `results/sample.graph_gaps.tsv` and `results/sample.scaffold_gaps.tsv`.
If graph evidence supports candidate gaps, plan fills and write a review page.
For the `--ref-paf` evidence file, see the Input Files notes on
[which PAF files to keep]({{ '/input-files/' | relative_url }}#which-paf-files-to-keep)
so the PAF query names match the graph nodes being scored.

```bash
chromo gapfill \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --gfa assembly_graph.gfa \
  --gaf reads_to_graph.gaf \
  --hic-pairs graph_contacts.tsv \
  --ref-paf paf/graph_nodes_to_ref.paf \
  --output-prefix results/sample.gapfill \
  --review-html results/sample.gapfill.review.html \
  --include-fill-sequences
```

Open `results/sample.gapfill.review.html`, review candidate paths, and export a
reviewed TSV. Then apply only accepted fillable rows:

```bash
chromo gapfill \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --gfa assembly_graph.gfa \
  --gaf reads_to_graph.gaf \
  --hic-pairs graph_contacts.tsv \
  --ref-paf paf/graph_nodes_to_ref.paf \
  --reviewed-plan chromosort.gapfill.reviewed_plan.tsv \
  --output-prefix results/sample.reviewed_gapfill \
  --apply
```

Review these outputs before publishing a filled scaffold:

- `results/sample.graph_gaps.tsv`: graph adjacency, missing graph nodes, direct
  edges, and short paths at scaffold junctions.
- `results/sample.gapfill.gapfill_plan.tsv`: candidate path status, support
  counts, risk flags, fill length, right-trim bp, and apply status.
- `results/sample.gapfill.review.html`: side-by-side review of candidate paths.
- `results/sample.reviewed_gapfill.gapfilled.fa`: final FASTA after accepted
  fills only.

## Handling Overlapping Contigs

Large contigs sometimes have strong reference support but overlap at their ends.
This can happen with alternate graph paths, small assembly redundancies, or real
terminal dovetails that another scaffolding workflow might trim or merge. In
ChromoSort, overlap handling is deliberately split across commands so the
sequence-changing step stays explicit.

`chromo sort` assigns and filters contigs in reference space. It lets stronger
contigs claim reference intervals first, then asks whether lower-ranked contigs
add useful novel reference coverage. Fully contained or internal redundant
fragments are still reported as `duplicate_overlap`. Dovetail-like one-sided
overlaps are now classified separately as `terminal_overlap`; if the contig is
kept, its status is `kept_terminal_overlap`. A mostly overlapping contig can also
be rescued when its one-sided extension passes
`--min-terminal-extension-bp` and `--min-terminal-extension-frac`.

`chromo fix` does not resolve overlap between two separate contigs. It only
splits within a selected contig when query-ordered alignment blocks support a
reference or eligible orientation transition. If two already-separate contigs
overlap at their ends on the same reference, `fix` will not trim, merge, or
deduplicate them.

`chromo scaffold` joins the final sorted contigs. In the default
`--overlap-policy zero-gap` mode, adjacent negative reference gaps are written as
zero-length gaps: no Ns are inserted and neither contig is trimmed. The raw
negative inferred gap, overlap bp, overlap class, overlap fractions, policy, and
action are reported in `<prefix>.scaffold_gaps.tsv`, and scaffold-level overlap
and trimming totals are reported in `<prefix>.scaffold_summary.tsv`.

When you want sequence surgery at scaffolding time, choose an explicit overlap
policy:

```bash
# Current conservative behavior, with stderr warnings for negative gaps.
chromo scaffold \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --output-prefix results/sample.warn \
  --overlap-policy warn

# Trim the right contig by the reference-inferred terminal overlap.
chromo scaffold \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --output-prefix results/sample.trim_ref \
  --overlap-policy trim-reference

# Trim only when the left suffix and right prefix confirm the overlap sequence.
chromo scaffold \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --output-prefix results/sample.trim_seq \
  --overlap-policy trim-sequence \
  --trim-sequence-min-identity 0.98
```

The safest default is still not to trim: overlapping reference spans can reflect
assembly redundancy, true variation, reference differences, or alignment
artifacts. The new reports make the case easy to find, and the trimming policies
make the intervention deliberate when a reviewed dataset needs it.

## Synthetic Graph Workflow

The repo ships a tiny graph-aware fixture under `tests/data/graph_gotchas`.
It is intentionally small enough to inspect by eye, but it exercises the same
file types used in real graph-aware ChromoSort runs: FASTA, PAF, GFA, GAF, and
Hi-C-like contact counts. The same file-type expectations are summarized on the
[Input Files]({{ '/input-files/' | relative_url }}#graph-input-files) page.

```bash
mkdir -p results/graph_gotchas
DATA=tests/data/graph_gotchas

chromo sort \
  --ref-fasta "$DATA/ref.fa" \
  --assembly-fasta "$DATA/assembly.fa" \
  --paf "$DATA/unitig_to_ref.paf" \
  --gfa "$DATA/unitigs.gfa" \
  --output-prefix results/graph_gotchas/sort

chromo manual \
  --ref-fasta "$DATA/ref.fa" \
  --assembly-fasta "$DATA/assembly.fa" \
  --paf "$DATA/unitig_to_ref.paf" \
  --gfa "$DATA/unitigs.gfa" \
  --embed-sequences \
  --output-html results/graph_gotchas/manual.html

chromo scaffold \
  --ordered-fasta results/graph_gotchas/sort.ordered.fa \
  --assignments results/graph_gotchas/sort.contig_assignments.tsv \
  --gfa "$DATA/unitigs.gfa" \
  --output-prefix results/graph_gotchas/scaffold
```

Open `results/graph_gotchas/manual.html` to inspect dot plots beside graph
neighbors. The scaffold graph report at
`results/graph_gotchas/scaffold.graph_gaps.tsv` shows report-only GFA context
for adjacent sorted contigs.

For a focused gap-fill example, the fixture includes `gapfill_ordered.fa` and
`gapfill_assignments.tsv`, a two-flank chr1 case where `bridge_good` and
`bridge_alt` are both possible graph paths:

```bash
chromo gapfill \
  --ordered-fasta "$DATA/gapfill_ordered.fa" \
  --assignments "$DATA/gapfill_assignments.tsv" \
  --gfa "$DATA/unitigs.gfa" \
  --ref-paf "$DATA/unitig_to_ref.paf" \
  --gaf "$DATA/reads_to_graph.gaf" \
  --hic-pairs "$DATA/hic_pairs.tsv" \
  --output-prefix results/graph_gotchas/gapfill \
  --include-fill-sequences \
  --review-html results/graph_gotchas/gapfill.review.html
```

The gapfill plan should mark `left+,bridge_good+,right+` as fillable, while the
review HTML shows both candidate paths side by side with PAF, GAF, Hi-C, and
risk annotations. After reviewing the HTML, export a reviewed TSV, or script the
expected toy approval:

```bash
python - <<'PY'
import csv

src = "results/graph_gotchas/gapfill.gapfill_plan.tsv"
dst = "results/graph_gotchas/gapfill.reviewed_plan.tsv"

with open(src, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    rows = list(reader)
    fieldnames = reader.fieldnames

for row in rows:
    row["accept_fill"] = (
        "yes"
        if row["fill_status"] == "fillable"
        and row["path_nodes"] == "left+,bridge_good+,right+"
        else "no"
    )

with open(dst, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
PY

chromo gapfill \
  --ordered-fasta "$DATA/gapfill_ordered.fa" \
  --assignments "$DATA/gapfill_assignments.tsv" \
  --gfa "$DATA/unitigs.gfa" \
  --ref-paf "$DATA/unitig_to_ref.paf" \
  --gaf "$DATA/reads_to_graph.gaf" \
  --hic-pairs "$DATA/hic_pairs.tsv" \
  --reviewed-plan results/graph_gotchas/gapfill.reviewed_plan.tsv \
  --output-prefix results/graph_gotchas/gapfill.reviewed \
  --apply \
  --simple-headers
```

The reviewed gapfilled FASTA should contain chr1 with the graph-supported bridge
inserted and the right flank overlap trimmed.
