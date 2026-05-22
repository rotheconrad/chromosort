---
title: chromo cut
description: Usage, outputs, parameters, and reasoning for chromo cut.
---

# chromo cut

Use `chromo cut` when you already know the exact breakpoint position for one or
more contigs and want a transparent FASTA edit without running the alignment
planner in `chromo fix`.

Positions are 1-based and mean "cut after this base." For example, cutting a
10 bp contig at position 4 emits bases 1-4 and 5-10. Terminal cuts are rejected
because they would create empty pieces.

## Run `chromo cut`

Repeat `--cut` for multiple contigs, and use commas for multiple positions on
the same contig:

```bash
chromo cut \
  --assembly-fasta assembly.fa \
  --cut contig_1:234567,450000 \
  --cut contig_2:10000 \
  --output-fasta results/sample.cut.fa \
  --report results/sample.cut_contigs.tsv
```

For a single contig, the original explicit form is also supported:

```bash
chromo cut \
  --assembly-fasta assembly.fa \
  --contig contig_1 \
  --pos 234567 450000 \
  --output-fasta results/sample.cut.fa \
  --report results/sample.cut_contigs.tsv
```

For batch review, provide a simple cuts file:

```text
contig	position
contig_1	234567	450000
contig_2	10000
```

Then run:

```bash
chromo cut \
  --assembly-fasta assembly.fa \
  --cuts-file reviewed_cuts.tsv \
  --output-fasta results/sample.cut.fa \
  --report results/sample.cut_contigs.tsv
```

## `chromo cut` Outputs

| Output | Description |
| --- | --- |
| `--output-fasta` | Full assembly FASTA with requested contigs replaced by cut pieces. Uncut contigs are copied unchanged. |
| `--report` | TSV report describing every emitted cut piece, including original contig, new contig, slice coordinates, piece length, and cut positions. |

Cut pieces are named `CONTIG_cut001`, `CONTIG_cut002`, and so on by default.
Change the inserted text with `--name-separator`, or use `--simple-headers` to
write only the new FASTA IDs in cut-piece headers.

### Example `chromo cut` Output

**Table 1. Example cut report rows.** Cutting `contigA` after bases 20 and 50
emits three adjacent slices and records the complete cut set on each row.

| original_contig | new_contig | part_index | slice_start | slice_end | piece_bp | cut_after_positions |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| `contigA` | `contigA_cut001` | `1` | `1` | `20` | `20` | `20,50` |
| `contigA` | `contigA_cut002` | `2` | `21` | `50` | `30` | `20,50` |
| `contigA` | `contigA_cut003` | `3` | `51` | `80` | `30` | `20,50` |

**Listing 1. Example cut FASTA records.** Requested contigs are replaced by cut
pieces; unrequested contigs are copied through unchanged.

```text
>contigA_cut001
AAAAAAAAAAAAAAAAAAAA
>contigA_cut002
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>contigA_cut003
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>contigB
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
```

## `chromo cut` Parameters

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `--assembly-fasta` | required | Assembly FASTA containing the contigs to cut. |
| `--assembly-fai` | auto | Optional FASTA index used for length validation. Defaults to `<assembly-fasta>.fai` when present. |
| `--cut` | none | Reviewed cut as `CONTIG:POS[,POS...]`; may be repeated. |
| `--cuts-file` | none | Optional file with contig and one or more cut positions per line. |
| `--contig` / `--pos` | none | Convenience form for one contig with one or more positions. |
| `--output-fasta` | required | New FASTA path. The input FASTA is never modified. |
| `--report` | required | TSV report path for the cut pieces. |
| `--min-piece-bp` | `1` | Reject cut plans that create pieces shorter than this length. |
| `--name-separator` | `_cut` | Text inserted before the numeric suffix. |
| `--simple-headers` | off | Write cut-piece FASTA headers containing only the new sequence ID. |
