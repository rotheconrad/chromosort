# Public assembly-error benchmark

Reproducible, public-sequence development cases for testing assembly correction,
sequence preservation and review evidence. **Three runnable fixtures** are
included; their public intervals can regenerate **all 14 first-cohort cases
offline**. The [second development panel](cycle2/README.md) adds eight cases
from new wheat loci and a separate paired evaluation. Mapping needs minimap2.
The [third development panel](cycle3/README.md) evaluates review-proposal
coverage and burden with automatic correction held fixed. The independent
scorers and thirty-one analytic tests use the Python standard library.

This is a repository dataset, separate from the installable ChromoSort package.
It contains no private project sequence, held-out inputs/truth, evaluator key,
tool installation, third-party software archive or bulk corrected assemblies.

## Start here

Run from this directory in a clone of the repository, using Python 3.9–3.12:

```sh
python3 scripts/verify.py --fixtures
python3 -m unittest discover -s tests -v
python3 scripts/summarize_history.py
python3 scripts/run_benchmark.py --tool unedited --output-root work/unedited
```

The checked-in examples have exact input/reference FASTA, mapped assembly/read
PAF, checksums, reference labels and development truth:

| Example | Manifest | Purpose |
| --- | --- | --- |
| D001 | [Wheat input](fixtures/development/D001/manifest.json) | Reconstruct a published scaffold-orientation decision; retain original bases and Ns |
| D003 | [Barley input](fixtures/development/D003/manifest.json) | Small documented-pattern error plus a separate simulated biological inversion |
| D010 | [Labeled wheat input](fixtures/development/D010/manifest.json) | Preserve two simulated biological exchanges despite reference discordance |

`fixtures/truth/` is openly available development truth. The compact manifests
omit raw FASTQ/BAM and duplicate combined-reference files; they keep the actual
PAF used by rc2 and the identities needed to regenerate the reads. The
[fixture index](fixtures/index.json) records every included byte and the
original-to-public manifest identities. Full generation restores raw reads.

## Regenerate the first development cohort

The first command reconstructs the five canonical GenBank intervals from the
fixtures and verifies their original accession/checksum registry. It needs no
network and does not invent or replace a sequence source.

```sh
python3 scripts/retrieve.py --output-root work
python3 scripts/generate.py --output-root work
python3 scripts/evidence.py work/development/D* --threads 2
python3 scripts/verify.py --regenerated work
```

Use minimap2 **2.31-r1302** for the recorded PAF identities. `--cases D003 D010`
on the generator selects a subset; `--bam` on the evidence mapper additionally
requires samtools. Generated data and outputs belong in ignored `work/` or an
explicit external directory. Existing cases/results are not overwritten.

The [reproduction guide](docs/reproduction.md) includes isolated setup,
NCBI download audits, paired ChromoSort commands, RagTag commands and expected
checks. Supply an explicitly identified software installation for new runs.
Do not compare an unrecorded moving checkout with historical results.

## Cases and evidence discipline

First-cohort definitions and seeds are in [configs/development.json](configs/development.json).
The [case table](docs/cases.md) distinguishes documented reconstructions,
documented-pattern simulations, actual tandem copies and simulated biology.
Reads are generated from the biological target **before** its assembly copy
is corrupted. The iid substitution model is transparent, not calibrated to
an empirical sequencing platform. No graph assembly or inferred phasing is
claimed by this fixture set.

Assembly-to-reference alignments describe one exact sequence stage. Re-align
changed FASTA before using it in a later alignment-dependent operation. Read
PAF and BAM derived from the same reads are correlated evidence. A conflicting
reference is not by itself proof of an assembly error.

## Recorded results and limits

The repository preserves **294 historical development score records**: 154 pilot, 56 rc1
and 84 rc2. Held-out records and scores are excluded. [Results](docs/results.md)
explain the source and cohort limits, paired read policies and RagTag fairness.
[summary.tsv](reports/summary.tsv) is regenerated from the original checksummed
[score records](reports/historical_scores.jsonl); it does not require private
archives or an installed assembly tool.

The rc2 review policy avoided observed false cuts on the first development cohort
while retaining the previously localized boundaries. Comprehensive/sensitive
still localized only **3 of 9** expected error boundaries. D003 can now be
inspected but is not automatically repaired. The development examples informed
the policy; these results are not independent validation or general superiority.

The [second-panel paired report](cycle2/results.md) adds **60 output records**.
Published baseline and rc3 have identical automatic outputs on all sixteen
paired arms, with no complete automatic repair in eight new cases. Both the
baseline manual API and rc3's explicit reviewed workflow complete eight new
cases plus D003/D012 using supplied edit coordinates. Those scripted repairs
measure replay feasibility, separately from discovery or human review accuracy.

The [third-panel paired report](cycle3/results.md) records one frozen rc4
evaluation. Both versions admitted seven of eight primary cases. Review
proposal coverage increased from 0/4 to 4/4 evaluable errors, with queue size
increasing from 1 to 19 and biological-control inspection burden from 1 to 14
events. All fourteen successful primary automatic FASTA/ledger pairs remain
identical; none of the four evaluable error cases was completely repaired.
The emitted-coordinate D201 replay validates and restores the exact target DNA
sequence, but fails the frozen complete-target provenance test at its two
equivalent edge bases. This sequence comparison is a labeled post-score
diagnostic. Three known diagnostics are reported separately, and the frozen
D208 manifest-admission failure remains explicit.

The [scoring protocol](docs/scoring.md) keeps input accounting, target copies,
cut localization and preserved continuity distinct. The current scorer covers
input-derived sequence and N gaps; donor fills and explicit trim-policy
accounting need additional validated truth formats.

See [source attribution](docs/sources.md), [data/software notice](NOTICE.md),
[software identities](sources/software_identities.json) and the
[next-validation boundary](docs/next-validation.md).
