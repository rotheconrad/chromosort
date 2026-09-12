# Independent scoring protocol v1

[`scripts/score.py`](../scripts/score.py) imports generic file readers only. It never calls ChromoSort,
RagTag, alignment decision code, or the case generator. Truth and input maps
are frozen separately before execution. The scorer validates every input-map
base against target FASTA and every reported output component against the
original input FASTA, including reverse complements. A false provenance ledger
fails validation rather than receiving an accuracy score.

## Input and output

Provide `--truth fixtures/truth/D003` (or generated `work/evaluator/development/D003`), `--assembly` matching its digest,
`--output-fasta`, and `--ledger`. Ledger TSV fields are `output_id`,
`output_start`, `output_end`, `input_id`, `input_start`, `input_end`, `strand`.
Coordinates are 0-based, half-open. `input_id=.` means an N-only output gap.
All output bases must be explained exactly once. Unknown orientations are
rejected. [`scripts/adapt_agp.py`](../scripts/adapt_agp.py) adapts standard AGP; `--ragtag-correct` handles
the reversed object/component direction of RagTag correction AGP. `--parent`
composes a scaffold AGP through an intermediate correction ledger.

The current scorer covers input-derived sequence and N gaps. Filling new
sequence requires an independently verified donor/target map extension before
it can be scored here. It must not receive invented accuracy from an unsupported
ledger. Whole-record renaming, ordering and reverse complements are free.

## Metrics and interpretation

* **Input accounting:** union of retained input coordinates, unaccounted input
  bases, and duplicated input bases. `retain_all_exact` requires every input
  base exactly once. Deliberate overlap trimming may appropriately reduce this
  metric and must be explained in an explicit discard/policy record. Do not
  call a loss unexplained if an external reviewed trim record accounts for it.
* **Target content:** union of covered ACGT target coordinates divided by all
  ACGT target coordinates. Missing and extra-copy bases are also counted.
  Non-ACGT target positions are reported separately. Provenance identity keeps
  near-identical tandem copies distinct, avoiding alignment reassignment of a
  missing copy onto its retained homologue.
* **Cuts:** newly introduced interior input boundaries with retained sequence
  on both sides. Match cuts one-to-one to predeclared error boundaries within
  500 bp. Report precision, recall, localization distances, false-cut counts
  and false cuts per callable Mb. Empty precision/recall denominators are null.
  Terminal deletion is content loss, not automatically a correction cut.
* **Structural relationships:** compose each output path into biological-target
  coordinates. Adjacent pieces with changed chromosome/orientation, overlaps,
  or a target distance inconsistent with output gap length are discrepancies.
  Correct separation into contigs removes adjacency obligations but does not
  restore continuity. Gap error is output N length minus target interval length.
* **Preservation:** per-control retained bases and whether a single exact,
  contiguous target path spans the control. An intact but fragmented control
  can retain all its bases while losing continuity.
* **Assignments (optional):** an output-level TSV with `output_id`, `status`,
  `output_group`. Report assigned bases, correct assigned bases, accuracy among
  assigned bases and inappropriate forced placement of ambiguous controls.
  Missing assignment tables yield `not_supplied`, not zero accuracy. The first
  implementation scores group assignment; chromosome/subgenome-specific
  refinements need explicit output columns and independent validation.

Cut recovery is not called complete error repair. Read cut metrics alongside
retention, copy content, structural discrepancies and preserved continuity.
Collapsed missing sequence is reported separately and is not required to be
reconstructed without an available sequence source. This benchmark does not
reward reference agreement for its own sake.

## Validation and statistics

Analytic tests construct adversarial sequence/provenance examples independently
of the generator: changed orientation, omitted and duplicated identical copies,
unnecessary cuts, wrong gap distances, false provenance and chained AGP maps.
Bundle verification checks checksums and exact PAF query/target lengths and
CIGAR spans. The unedited baseline scores every frozen map.

Report paired per-case results, with scenarios/seeds nested within source loci.
The expanded panel varies read model, target divergence and event settings on
shared loci; its seeds are technical/scenario replicates, not independent
genomes. Do not use events or reads as independent replicates for confidence
intervals, pool correction/scaffolding arms, or claim population-wide
superiority from the pilot. Reviewed/oracle workflows remain separate from
automatic workflows.
