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

Recommended order when raw dot plots show possible misjoined contigs:

```bash
# 1. Fix only reviewed/suspect raw contigs, or use --all to scan every contig.
chromo fix \
  --assembly-fasta assembly.fa \
  --coords mummer/raw.coords \
  --contigs suspect_contig_1 suspect_contig_2 \
  --output-fasta results/sample.fixed.fa \
  --report results/sample.fixed_contigs.tsv

# 2. Re-align the fixed FASTA with MUMmer, then sort the fixed assembly.
chromo sort \
  --ref-fasta reference.fa \
  --assembly-fasta results/sample.fixed.fa \
  --coords mummer/fixed.coords \
  --output-prefix results/sample.fixed \
  --orient-to-reference

# 3. Plot from the existing coordinates; no new MUMmer run is needed.
chromo plot \
  --ref-fasta reference.fa \
  --assembly-fasta results/sample.fixed.fa \
  --coords mummer/fixed.coords \
  --assignments results/sample.fixed.contig_assignments.tsv \
  --output-prefix plots/sample.fixed \
  --per-ref
```

Running `chromo fix` before `chromo sort` is safest for suspected misjoins
because sorting is a placement/filtering step, not a splitter. `chromo sort`
now protects strong multi-reference split candidates by default, but reviewed
raw-assembly fixing is still the cleaner workflow when the dot plot already
shows which contigs need attention.

Typical reference-ordering run:

```bash
chromo sort \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.fa \
  --coords mummer/sample.coords \
  --output-prefix results/sample \
  --orient-to-reference
```

The same command accepts minimap2 PAF instead of MUMmer coords:

```bash
chromo sort \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.fa \
  --paf paf/sample.paf \
  --output-prefix results/sample \
  --orient-to-reference
```

Typical chimeric-contig fixing run:

```bash
chromo fix \
  --assembly-fasta assembly.fa \
  --coords mummer/sample.coords \
  --contigs contig_04 contig_12 \
  --output-fasta results/sample.fixed.fa \
  --report results/sample.fixed_contigs.tsv
```

Manual cut at reviewed positions:

```bash
chromo cut \
  --assembly-fasta assembly.fa \
  --cut contig_1:234567,450000 \
  --cut contig_2:10000 \
  --output-fasta results/sample.cut.fa \
  --report results/sample.cut_contigs.tsv
```

Manual GUI review and editing:

```bash
chromo manual \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.fa \
  --coords mummer/sample.coords \
  --gfa assembly_graph.gfa \
  --output-html results/sample.manual.html \
  --suggested-output-fasta sample.manual.fa
```

Scan all contigs with the same conservative planner:

```bash
chromo fix \
  --assembly-fasta assembly.fa \
  --coords mummer/sample.coords \
  --all \
  --output-fasta results/sample.fixed.fa \
  --report results/sample.fixed_contigs.tsv
```

Typical scaffolding run after final sorting:

```bash
chromo scaffold \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --output-prefix results/sample
```

Plan and apply reviewed graph-supported fills:

```bash
chromo gapfill \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --gfa assembly_graph.gfa \
  --output-prefix results/sample.gapfill \
  --review-html results/sample.gapfill.review.html

chromo gapfill \
  --ordered-fasta results/sample.ordered.fa \
  --assignments results/sample.contig_assignments.tsv \
  --gfa assembly_graph.gfa \
  --gaf reads_to_graph.gaf \
  --hic-pairs graph_contacts.tsv \
  --reviewed-plan chromosort.gapfill.reviewed_plan.tsv \
  --output-prefix results/sample.reviewed_gapfill \
  --apply
```

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
Hi-C-like contact counts.

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
