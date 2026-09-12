# Second development panel: mixed structure and explicit repair

This additive panel preserves the first fourteen development cases and their
scorer. It contains eight scenarios from new public Chinese Spring wheat loci,
plus separate scripted repair diagnostics for D003 and D012. These are
**development experiments**, not unseen validation, independent genomes, or a
study of human review accuracy.

The [prospective protocol](protocol.json), [case definitions](cases.json), and
[source admission](sources/source_admission.json) were hashed in
[the design freeze](prospective-freeze.json) before generator implementation or
software outputs. The baseline is public commit
`dd66e3d987e775d99e612f4cdc19adec039b43e6`. Candidate identity and execution records
are supplied separately; neither baseline nor candidate is tuned after scoring.

| Case | Construction | Read-evidence condition |
|---|---|---|
| D101 | Reciprocal internal A/D exchange is simulated target biology; the two 200 kb targets are then incorrectly joined into one 400 kb input | 12× target-derived 10 kb reads |
| D102 | Same target and input geometry | 1× coverage |
| D103 | Same target and input geometry | 12×, 600 bp reads; shorter than required two 1 kb anchors |
| D104 | Same target, input and mapped reads as D101 | Every read-PAF CIGAR tag removed |
| D105 | Same mixed contig plus one artificial input-derived chimera across the true misjoin | Its mapped row appears 20 times under one name |
| D106 | Same parent chimera as D105 | Two copies of its mapped row have different names |
| D107 | Reported scaffold139 inversion in original wheat chromosome 1A; original gaps retained | 12× reads from the reconstructed target, excluding gap-spanning reads |
| D108 | New 300 kb 1A locus with an injected 8 kb assembly inversion and a separate simulated biological 15 kb inversion | 12× target-derived 10 kb reads |

D101–D106 exchange `[60000,120000)` between the two target sequences. The
assembly error joins complete targets at input position 200000. Four 10 kb
preservation intervals surround the two exchange edges in each target. The six
conditions share sequence and target geometry; they are paired evidence
contrasts, not biological replicates. Reads use the existing transparent
substitution simulator: uniform eligible starts, random strands, 0.001
substitutions, no indels, constant nominal Q30, and a fixed seed. Requested depth
uses ACGT runs long enough for a read. It is not an empirical sequencing model.

D105/D106 intentionally violate independent-molecule support. Their added
sequence is drawn from the corrupted assembly across the known misjoin. The
PAF transformations are explicitly recorded in the input manifests; raw FASTQ
contains the single parent chimera. D106's second name is an alias introduced
in PAF. Neither duplicate rows nor distinct names prove physical independence.
No tool is expected to infer the hidden common origin from sequence identifiers.

## New source admission

Four versioned public intervals are bundled in
[the sequence registry](sources/sequence_registry.json). Two 200 kb A/D intervals
supply the mixed cases. The other two contain the original scaffold139 region
and its corrected guide. D108 uses a gap-free 300 kb interval of that guide.
No original first-cohort locus is substituted into these eight cases.

[Zhu et al. (2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8360199/) report
optical-map-guided scaffold corrections. Table S3, sheet `v2.1 vs v1.0`, row 149,
places original `chr1A_scf_139` at 579062096–579708165 and gives its new
orientation as reverse. Independent source-admission mapping confirms a reverse
central interval with forward flanks in the versioned corrected sequence;
[alignment coordinates](sources/source_alignment_summary.json) are retained.
The inclusive endpoints imply 646070 bp, one more than the table's stated
length. D107 therefore uses an operational orientation reconstruction and the
unchanged 500 bp cut tolerance. It does not reconstruct every change in the
corrected assembly, reanalyse optical-map molecules, or claim base-exact
historical breakpoint truth. Its original gap Ns remain, and simulated reads
cannot validate joins across them.

The source paper and supplemental workbook are cited, not redistributed here.
Sequence and software rights follow the benchmark [notice](../NOTICE.md).

## Outcomes and interpretation

The primary comparison pairs comprehensive-mode `report` and `review` policies
on the same inputs and evidence. It measures cut localization, extra cuts,
retention, biological continuity, and complete target structure separately.
A deferred whole-contig plan may preserve biological variation while leaving a
true misjoin. A pending inspection item is neither an error confirmation nor a
sequence repair.

[The added structure scorer](../scripts/score_structure.py) first invokes the
unchanged exact provenance scorer, then requires one complete, contiguous output
path per target identity. Record order and whole-target reverse complements are
free. Missing or duplicated copies, fragmentation, internal orientation/order
errors, and added gaps fail completeness. This prevents a split assembly with
zero remaining adjacent discrepancies from receiving full-repair credit.
Analytic tests use independent small geometries and an exact ledger example.

The explicit-repair arm supplies the frozen edit coordinates: one true cut for
each mixed case, and the error interval reversal for each orientation case.
The baseline uses its existing manual recipe API. The candidate uses a real
scan event and an explicit reviewed edit or exact existing cut decision, then
previews, applies and replays the recipe. Other proposed changes are rejected.
This tests **scripted, oracle-coordinate repair feasibility**, not discovery,
reviewer effort or autonomous biological judgment. If no usable event exists,
that workflow is recorded as unavailable. D003/D012 are governed by
[their separate prespecified diagnostics](prior-case-diagnostics.json), outside
the eight-case primary denominator.

## Regeneration

From the benchmark directory, with Python 3.9–3.12 and minimap2 2.31-r1302:

```bash
python3 scripts/generate_cycle2.py --output-root work/cycle2-panel
python3 scripts/verify_cycle2.py --generated-root work/cycle2-panel
python3 -m unittest discover -s tests
```

Generation creates `development/` inputs and a separate `evaluator/development/`
truth tree. Give curation sessions only the input tree. FASTQ and full runtime
outputs remain in ignored work directories. The generator records mapper
commands, versions, exact FASTA identities, and evidence derivations. The
[expected file identities](generated_expected.json) cover all eight generated
cases. [Regeneration verification](reproduction-verification.json) checks those
identities and the explicit PAF transformations, without scoring outputs.

Use `scripts/run_cycle2.py --help` for installed, immutable software execution.
It runs commands without opening evaluator truth. Both versions' output hash
ledgers must exist before independent scoring. The public baseline can be
installed from its commit; candidate archive identities do not imply that a
public release archive is available.

The completed [paired results](results.md) include all failures to discover or
fully repair structure, exact small output artifacts and portable score checks.

For the paired run, install the public baseline commit and the candidate whose
identities are in [the execution freeze](execution-freeze.json), using separate
Python environments. `BASELINE_PYTHON` and `CANDIDATE_PYTHON` below name those
installed interpreters, not a moving source checkout. Regenerate the two prior
diagnostics separately:

```bash
python3 scripts/retrieve.py --output-root work/cycle2-prior
python3 scripts/generate.py --output-root work/cycle2-prior --cases D003 D012
python3 scripts/evidence.py work/cycle2-prior/development/D003 work/cycle2-prior/development/D012 --threads 2
python3 scripts/run_cycle2.py --python "$BASELINE_PYTHON" --role baseline \
  --software-id dd66e3d987e775d99e612f4cdc19adec039b43e6 \
  --case-root work/cycle2-panel/development --prior-case-root work/cycle2-prior/development \
  --output-root work/cycle2-baseline
python3 scripts/run_cycle2.py --python "$CANDIDATE_PYTHON" --role candidate \
  --software-id 6e090a5174ff32b301666aba5f16c0c6d3789ffbe6e65984052c7d914c96f173 \
  --case-root work/cycle2-panel/development --prior-case-root work/cycle2-prior/development \
  --output-root work/cycle2-candidate
python3 scripts/score_cycle2.py --baseline-root work/cycle2-baseline \
  --candidate-root work/cycle2-candidate --panel-root work/cycle2-panel \
  --prior-root work/cycle2-prior --output-root work/cycle2-scores
python3 scripts/export_cycle2.py --baseline-root work/cycle2-baseline \
  --candidate-root work/cycle2-candidate --results-root work/cycle2-scores
```

For a software-independent verification of the **published** score artifacts:

```bash
python3 scripts/verify_cycle2_results.py --panel-root work/cycle2-panel --prior-root work/cycle2-prior
python3 scripts/check_clean_checkout.py --with-cycle2
```

The clean-copy check uses only this public benchmark directory, Python and
minimap2. It regenerates all eight cases, verifies the first cohort's compact
fixtures and 294 historical scores, and reproduces all sixty second-cycle
metrics from the exported component ledgers. It does not run an installed
candidate or claim a second independent candidate pass.
