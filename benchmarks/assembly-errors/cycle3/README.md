# Cycle 3: review-proposal coverage and burden

This development cycle measures what a researcher can inspect and select.
Automatic correction defaults remain fixed. Review proposals, applied cuts,
complete target structure and a supplementary scripted replay are separate
outcomes. The eight primary scenarios reuse registered public wheat sequences;
they are not independent genomes, empirical biological variants or blinded
reviewer experiments. D101, D012 and D108 are known diagnostics, reported
separately.

[The paired results](results.md) report the completed single evaluation,
including the additional inspection burden, D208 admission failure and
emitted-coordinate replay that restores exact target DNA but fails the frozen
provenance-path test. The sequence comparison is labeled as a post-score
interpretation diagnostic.

The public baseline is `c7f1a9fd11f1ea79e0b86d695984f70b63c2a8f7`.
[The prospective freeze](prospective-freeze.json) binds case definitions,
matching/scoring rules, commands and metric code before baseline or candidate
curation outputs. [The accepted event contract](review-event-contract.json) and
[its schema clarification](schema-clarifications.json) specify how the two
versions' events are compared. [Execution identities](execution-freeze.json)
record the exact rc4 archive and installed module hashes before candidate
execution. One candidate and one paired evaluation are authorized;
there is no post-score tuning or replacement case.

## Cases and construction

| Case | Fixed construction | Role |
|---|---|---|
| D201 | 8 kb input inversion at `[70000,78000)` plus separate target-biological inversion at `[210000,225000)` | Short aligned error and biological burden |
| D202 | 25 kb input inversion at `[110000,135000)` plus the same biological inversion | Larger aligned error and smoothing |
| D203 | Biological inversion only | Preservation control; inspection is burden, not false correction |
| D204 | Unmodified target and input | Clean control |
| D205 | Complete input record reverse-complemented | Whole-record orientation control, structurally equivalent target |
| D206 | D201 sequence/read geometry, three copies of every assembly PAF row and a nested 7 kb reverse row | Redundant and overlapping representations of the same mapping |
| D207 | Exchange `[50000,100000)` between two 200 kb A/D targets, then join complete targets at 200000 | Competing biological and error transitions |
| D208 | D201 geometry with native and reverse-complement guide frames declared under the same output label | Frozen multiple-frame construction; rejected at manifest admission |

All intervals use native input coordinates, zero-based and half-open. The
orientation cases use `[124441,424441)` of registered `wheat_v2_scf139`, a
300 kb ACGT-only interval. The mixed case uses the two registered 200 kb A/D
sequences. Source accessions, genomic intervals, hashes and attribution remain
in [the unchanged cycle 2 registry](../cycle2/sources/sequence_registry.json).
No new published historical breakpoint reconstruction is claimed.

Target biology precedes assembly corruption. Reads use the existing transparent
IID substitution simulator: 12× requested coverage, 10 kb length, 0.001
substitution rate, no indels, random strands and constant nominal Q30. This is
not a calibrated PacBio or ONT model. Related scenarios share seeds and target
reads where applicable. Reference mappings use minimap2 `asm5`; read mappings
use `map-hifi`, both with CIGAR and `--secondary=no`, two threads and recorded
version. No graph evidence is supplied; unavailable graph support must remain
explicit.

D206 introduces no independent evidence. Its duplicated PAF rows and nested
`[70500,77500)` reverse row have the same source mapping. During construction,
the mapped inversion trimmed one identical base from each edge, producing
`[70001,77999)` and `7998M`. [Amendment A001](construction-amendments.json), made
before protocol freeze or curation output, replaced an overly strict `8000M`
guard with verification of the exact covering row. The sequences, error
coordinates, nested interval and read model did not change. The failed first
construction remains a local runtime record.

D208's two coordinate backbones conflict under the existing manifest rules.
The [baseline status record](baseline-output-status.json) preserves that
admission failure. Its sequence/PAF/truth integrity can pass the independent
input verifier while the software's manifest contract rejects its labels.
The case remains unchanged. Unavailable proposal and burden measurements are
excluded from numeric aggregates and are never represented as zero-queue
success. Completion and successful-pair denominators must expose the failure.

## Frozen measurements

One event ID corresponds to one decision row. New structured metadata enriches
that row; it does not add another copy of the event. Pending `cut` and `inspect`
rows count toward the actionable queue. `keep` acknowledgments are counted
separately. Whole-record contexts and a zero-cut optimum are explanatory
metadata, not alternative queue events.

An orientation error is covered only when **one interval event matches both
native endpoints within 500 bp each**. Two cuts or plan points can contribute
to edge visibility but cannot substitute for that interval. A misjoin is
covered by a point event or explicit plan cut within 500 bp. Matching is greedy
and one-to-one, ordered by endpoint distance, event ID and error ID. A multi-cut
plan remains one queue item and can match only one error under this metric.
Additional points can still create biological or clean-region inspection
burden. These choices are explicit measurement conventions, not claims about
human cognition.

Biological burden counts events that touch declared preservation regions and
the number of regions touched. The primary inversion controls include 2 kb
flanks (for example `[208000,227000)` around the 15 kb biological inversion);
D207 uses a 4 kb region centered on each biological exchange boundary.
Clean-region burden counts events touching
regions outside the error tolerance windows and declared biological regions.
Plan burden uses its actual cut points, not the interval between its first and
last cut. Burden categories can overlap for a plan with several cuts. A flag
on preserved biological variation is a review task, not an automatically false
sequence correction.

Geometric redundancy uses deterministic greedy representatives, sorted by
contig, shape, start, end and ID. The same point/interval geometry, or an
equal-length ordered plan cut set, within 500 bp counts as redundant regardless
of event kind or reference frame. Point, interval and multi-cut plan shapes
remain distinct, and grouping does not chain through successive near matches.
Report the total queue and unique geometries alongside redundancy.

The independent [proposal scorer](../scripts/score_review_proposals.py) also
records stage-explanation availability. Reported stage labels and objective
components are descriptive evidence; their presence is not proof of biological
correctness. The original sequence/provenance scorer and complete-target-path
scorer remain unchanged. Successful automatic FASTA and component-ledger pairs
must be byte-identical, with failed runs and absent outputs reported separately.
Already-correct controls are distinguished from repaired error cases.

## Review and replay contract

Baseline scans use the fixed comprehensive settings in [the protocol](protocol.json).
Candidate scans add:

```text
--inspect-orientation --inspect-min-orientation-bp 1000
--alternative-plans 5 --review-max-blocks 128
```

The producer exposes native intervals or explicit plan cut points, supporting,
opposing and overlapping blocks, reference frames, mapping origin/dependency
identities, loss stages, per-coordinate evidence and uncalibrated objective
components. Events remain pending inspection and make no implicit edit. Input
bounds, skipped contigs and omitted-event counts remain visible. The bounded
alternative search is not an exhaustive inventory of possible segmentations.

The supplementary D201 replay uses known error geometry to select a matching
emitted interval by distance then ID. That is **oracle event selection**. It
uses the emitted endpoints without refining them to truth. If no matching
interval exists, the outcome is unavailable; no `keep` event or invented
interval is used as fallback. Otherwise an explicit reverse edit is reviewed,
previewed, applied, replayed and checked against fresh alignments. Residual
structural errors from imprecise emitted endpoints remain in the result.
This diagnostic does not measure autonomous discovery or blinded reviewer
accuracy.

## Reproduction

From the benchmark directory with Python 3.9–3.12 and minimap2 2.31-r1302:

```bash
python3 scripts/generate_cycle3.py --output-root work/cycle3-panel
python3 scripts/verify_cycle3_inputs.py --generated-root work/cycle3-panel
python3 scripts/prepare_cycle3_known.py --output-root work/cycle3-known
python3 -m unittest discover -s tests
python3 scripts/verify_cycle3_results.py --panel-root work/cycle3-panel --known-root work/cycle3-known
python3 scripts/verify_cycle3_context.py --panel-root work/cycle3-panel --known-root work/cycle3-known
python3 scripts/summarize_cycle3.py
python3 scripts/diagnose_cycle3_sequence_equality.py --panel-root work/cycle3-panel --out work/cycle3-sequence-diagnostic.json
```

[Expected file identities](generated_expected.json) and
[reproduction verification](reproduction-verification.json) cover all eight
primary constructions. [Known diagnostic identities](known-expected.json) are
separate. Only `development/` inputs go to curation sessions; the generated
`evaluator/development/` tree stays evaluator-side.

Use separate immutable installed interpreters for the public baseline and the
validated candidate. `scripts/run_cycle3.py --help` describes fixed automatic
and review execution without opening proposal truth. After both output hash
ledgers exist, `scripts/evaluate_cycle3.py --help` describes paired commitment,
the prespecified event-selection replay and independent scoring. Bulk reads,
figures, logs, software archives and sequence outputs remain in ignored work
directories; portable definitions, event data and provenance ledgers are public.

To repeat software execution, install the exact baseline archive and rc4 wheel
identified in [execution-freeze.json](execution-freeze.json) into separate
environments. The baseline archive is the complete `git archive` of the stated
commit; the rc4 wheel SHA-256 is
`aba814bc0538e701ed7d97f55ab7044f81c64e99d77dd2365c5e677aa9b825fb`.
These software archives are not bundled with the benchmark. Example commands
use installation paths supplied by the person reproducing the run:

```bash
python3 scripts/run_cycle3.py --python /path/to/baseline/bin/python \
  --software-id c7f1a9fd11f1ea79e0b86d695984f70b63c2a8f7 --role baseline \
  --panel-root work/cycle3-panel --known-root work/cycle3-known \
  --output-root work/cycle3-baseline
python3 scripts/run_cycle3.py --python /path/to/rc4/bin/python \
  --software-id 133a54b516f0aa3517f4574528f39fe60f822b31318b6b6c162a97464a6a1d3a --role candidate \
  --panel-root work/cycle3-panel --known-root work/cycle3-known \
  --output-root work/cycle3-candidate
python3 scripts/evaluate_cycle3.py \
  --baseline-root work/cycle3-baseline --candidate-root work/cycle3-candidate \
  --baseline-python /path/to/baseline/bin/python --candidate-python /path/to/rc4/bin/python \
  --panel-root work/cycle3-panel --known-root work/cycle3-known \
  --runtime-root work/cycle3-replays --output-root work/cycle3-results
```

All output directories must be new. Failed D208 arms are expected and retained.
The recorded evaluation binds the runner, collector and replay code hashes;
reproduction does not authorize another development candidate or parameter
search. The bulk output commitments preserve original file identities, while
public normalized JSON changes absolute paths to portable labels and joins
structured metadata to existing event IDs. The artifact index records exact
copies and the one recipe path removal.

`python3 scripts/check_cycle3_clean_copy.py` copies only public benchmark
content to an isolated temporary directory, regenerates inputs, reproduces
metrics and summary tables, and verifies earlier fixtures and historical
summaries. It needs no installed ChromoSort or private files. Original output
FASTA byte hashes remain recorded separately because portable ledger
reconstruction need not reproduce original headers or wrapping.
