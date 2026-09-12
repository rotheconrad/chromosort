---
title: Labeled references and copy-preserving curation
---

# Labeled references and copy-preserving curation

`sort`, `eval fix`, `fix`, `plot`, `manual`, `reads`, and `workflow` accept
`--manifest inputs.json`. Legacy `--ref-fasta` and `--paf`/`--coords` commands
remain supported. A manifest binds one exact assembly stage, user labels and
every evidence file to SHA-256 digests. It can describe one or several reference
genomes; a reference genome is distinct from its chromosome records.

```json
{
  "schema": "chromosort-input-v1",
  "assembly": {"path": "assembly.fa", "sha256": "<64 lowercase hex digits>",
               "sample_id": "target", "stage": "input"},
  "references": "references.tsv",
  "reference_sequences": "reference_sequences.tsv",
  "evidence": []
}
```

Paths, including paths inside TSV cells, resolve relative to the JSON manifest.
The two tables can instead be inline JSON row arrays. A complete runnable
controlled fixture is in `tests/data/labeled_references/`.

| Table | Required columns | Optional metadata |
| --- | --- | --- |
| `references` | `reference_id`, `fasta`, `sha256` | `role` (`backbone`, `supporting`, `alternative`), accession |
| `reference_sequences` | `reference_id`, `sequence_id`, `chromosome`, `output_group`, `backbone` | `subgenome`, `homoeolog_group`, `copy_label` |

Every reference FASTA record requires exactly one label row. IDs within a FASTA
must be unique. Empty optional labels mean unknown. `backbone` accepts booleans
or `yes`/`no`. A group is the tuple `(output_group, chromosome, subgenome,
copy_label)`. At most one record is its coordinate backbone. Without a usable
backbone alignment, placement remains unresolved. The reference-level `role`
is descriptive; the sequence-level `backbone` flag is authoritative.

Each evidence entry requires `evidence_id`, `kind`, `path`, `sha256`, and
`assembly_sha256`. Primary reference alignments use `assembly_paf` or
`assembly_coords` and also require `reference_id` and `reference_sha256`.
Choose one primary alignment file per reference genome. Alternative aligners
are separate sensitivity runs. MUMmer input must include reference/query lengths.
PAF names, lengths, bounds and malformed rows are checked. An alignment file's
declared digest establishes its identity; it cannot independently prove the
truthfulness of a user-supplied alignment command.

Read kinds are `read_paf`, `sam`, and `bam`; graph kinds are `gfa` and `gaf`.
Optional `features` and `variants` entries supply figure overlays. BAM indexes
use `index_path` and `index_sha256`. Record command argument arrays, software
versions, read sample/readset identity, `evidence_class`, and `dependency_group`
to distinguish same-sample evidence, cross-sample comparisons and correlated
sources. These metadata do not turn graph fixtures or simulated reads into
independent validation.

## Assignment and ambiguity

Reference keys encode `reference_id::sequence_id`, escaping reserved characters.
Repeated chromosome IDs in different genomes therefore remain distinct. For
query contig q and reference record r, the score is the integral of the maximum
alignment identity over covered query bases: each query base contributes at
most one identity-weighted base. Overlapping or repeated alignment rows cannot
increase that base's vote. For a labeled output group, take the maximum score
among its reference frames. Never sum alternate references as extra support.

Report every candidate in `*.reference_candidates.tsv`. Rank output groups by
their scores, with reference key as a deterministic tie-breaker. The relative
margin is `(best - alternate) / best`. A margin at or below
`--ambiguity-margin` (default 0.05) yields an ambiguous placement. This is an
operational score, not a posterior probability or inferred dosage. A tied best
reference within the same output group does not itself imply group ambiguity.

The selected group must pass `--min-aligned-bp` and `--min-query-cov`, and have
a passing alignment to its explicit backbone. Supporting reference coordinates
are reported but never substituted into backbone gap arithmetic. Orientation
disagreement between reference frames is reported, not declared an error.

```bash
chromo sort --manifest inputs.json --retain-all --orient-to-reference -o results/sort
chromo scaffold --ordered-fasta results/sort.ordered.fa \
  --assignments results/sort.contig_assignments.tsv -o results/scaffold
```

Manifest sorting preserves copies by default. `--filter-within-copy` enables
the legacy overlap filter only where explicit copy labels are supplied. It is
an intentional assembly-representation choice; similarity alone is not evidence
that a genuine copy is redundant. Contigs sharing at least half the shorter
aligned reference footprint are retained but withheld from automatic
scaffolding (`scaffold_eligible=no`), since overlapping copies/repeats require
review. This is a conservative abstention rule, not copy-number inference.

`--retain-all` emits every input contig exactly once, retaining unresolved
records in `ambiguous::` or `unplaced::` partitions. A filtering recommendation
does not delete sequence in this mode. Without `--retain-all`, excluded records
remain explicitly listed in `*.accounting.tsv`; `--discarded-fasta` can save
their sequence. No contig is duplicated into multiple output groups.

`*.contig_assignments.tsv` keeps legacy columns and adds best/alternate scores,
margin, supplied labels, partition, backbone, scaffold eligibility and digests.
Its `assigned_ref` and legacy `ref_start/ref_end` describe the backbone and use
the existing 0-based closed convention. The new candidate table and accounting
ledger use 0-based half-open intervals. AGP uses standard 1-based closed
coordinates. Component notes preserve labels; each stage has its own AGP.

## Correction and review

Manifest correction uses only explicit backbone frames. Overlapping backbone
segments from different groups are withheld when a competitor covers at least
half of the focal segment and has identity within two percentage points or
higher. Supporting genomes cannot manufacture transitions by interleaving their
coordinates. This conservative segment-level rule may reduce sensitivity near
close homoeologs; the candidate table retains their original evidence.

The existing chromosome, conservative, comprehensive and sensitive modes then
operate on the retained segments. Subgenome ancestry switches and inversions
are review candidates, not error truth. `workflow scan` additionally exposes
internal uncovered alignment intervals as inspection items without proposing a
cut. Biological variants, absent parental references, collapsed copies and
autopolyploid ambiguity can remain unresolved.

User labels are supported; de novo subgenome inference, dosage inference and
full phasing are outside this release. Use the
[supervised workflow]({{ '/commands/workflow/' | relative_url }}) to review individual decisions and
the [read panels]({{ '/commands/reads/' | relative_url }}) to explain their evidence.

Manifest correction additionally uses [read-continuity review safeguards]({{ '/review-refinement/' | relative_url }}) in rc2. The reference transition remains a proposal when reads support continuity.
