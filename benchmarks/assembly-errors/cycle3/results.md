# Cycle 3 paired development results

The frozen rc4 review options exposed proposals covering all four evaluable
primary errors, compared with none in the baseline. The actionable queue grew
from 1 to 19 events and biological-control inspection burden grew from 1 to 14
events. Automatic sequence output did not improve: all fourteen successful
primary FASTA/component-ledger pairs were byte-identical, and neither version
completely repaired any of the four evaluable error cases.

These are constructed development scenarios using previously registered loci,
with related cases sharing sequence and simulated reads. They do not estimate
accuracy across independent genomes, biological error classification, reviewer
performance or improvement over another assembly tool. No held-out evaluation
or post-score tuning was performed.

## Frozen execution and admission

The baseline is public commit `c7f1a9fd11f1ea79e0b86d695984f70b63c2a8f7`
(0.4.0rc3); the candidate is one independently installed 0.4.0rc4 wheel.
[The execution freeze](execution-freeze.json) binds the source archive, wheel,
installed modules, protocol, evaluator and replay code. Both installations use
Python 3.12.14. Case geometry, commands and scoring rules were frozen before
baseline or candidate curation outputs. [Paired output hashes](results/paired-output-commitment.json)
were recorded before event selection and scoring; [replay hashes](results/replay-output-commitment.json)
precede output-structure scoring.

Both versions completed **7/8 primary scans and 14/16 primary automatic arms**.
D208 declares two coordinate backbones under the same output label and both
versions reject it with `Multiple coordinate backbones for wheat_A::1A::A::.`
Its raw sequence/PAF/truth integrity passes independently, but its manifest
construction is inadmissible. This is a failed benchmark construction, not
evidence about review coverage on valid conflicting reference frames. No input
was rescued or replaced after the freeze. The fifth prespecified primary error
and ninth biological control therefore have unavailable review measurements.
All proposal aggregates below exclude D208, while completion retains the
original denominator of eight. No unavailable queue is counted as zero.

[Failure records](results/failures.json) preserve both versions' scan exit code
1 and combined scan logs, plus separate automatic stdout/stderr. The legacy
automatic runner recorded failed status but did not retain the numeric
subcommand exit code; it is explicitly null. There are no candidate event-schema
failures. The three known diagnostics all execute and remain separate.

## Primary proposal coverage and burden

The [frozen rules](protocol.json) require one interval matching both inversion
endpoints within 500 bp; a pair of plan cuts is only edge visibility. A point
error can match a point proposal or an explicit cut in a multi-cut plan.
Assignment is greedy and one-to-one by distance, event ID and error ID. Each
multi-cut plan is one queue item and can match at most one error.

| Metric, completed primary cases only | Baseline | rc4 |
|---|---:|---:|
| Completed primary scans | 7/8 | 7/8 |
| Evaluable errors covered | 0/4 | 4/4 |
| Error edges visible | 0/7 | 7/7 |
| Actionable pending queue | 1 | 19 |
| Events touching biological controls | 1 | 14 |
| Biological regions inspected | 1/8 | 8/8 |
| Events touching clean regions | 0 | 0 |
| Unique geometries / redundant events | 1 / 0 | 19 / 0 |
| Queue items with explanatory text | 1/1 | 19/19 |
| Separate keep acknowledgment rows | 4 | 4 |

| Case | Baseline → rc4 queue | Baseline → rc4 error coverage | rc4 biological-burden events | Main observation |
|---|---:|---:|---:|---|
| D201 | 0 → 3 | 0/1 → 1/1 | 2 | Error interval `[70001,77999)` matches both edges within 1 bp |
| D202 | 0 → 5 | 0/1 → 1/1 | 3 | Error interval `[110000,135000)` matches exactly |
| D203 | 0 → 2 | No error | 2 | Biological inversion creates an interval and a competing plan |
| D204 | 0 → 0 | No error | 0 | Clean control remains an empty queue |
| D205 | 0 → 0 | No error | 0 | Whole-record reverse orientation is context only |
| D206 | 0 → 3 | 0/1 → 1/1 | 2 | Duplicated/nested mapping rows retain D201's event geometries |
| D207 | 1 → 6 | 0/1 → 1/1 | 5 | First alternative proposes the true 200000 join |
| D208 | Unavailable | Unavailable | Unavailable | Same admission failure in both versions |

rc4's nineteen primary queue items comprise seven orientation intervals,
eleven alternative segmentations and the existing automatic cut proposal.
The fourteen biological-burden events are inspection work, not false automatic
corrections. A plan may touch both an error and biology. The zero redundant
count uses the declared shape-aware geometry rule: an interval and a plan
with similar edges remain distinct. It does not mean a researcher sees no
overlapping or repetitive biological questions.

[Per-case scores](results/proposal-scores.jsonl), [queue geometry and stage labels](results/candidate-events.tsv)
and [cohort aggregates](results/cohort-summary.json) make the matching and burden
assignments inspectable. Coverage is an existence metric, not confidence or
ranking accuracy.

## What the new context explains

D201 and D206 retain an aligned 7,998 bp reverse interval below the unchanged
10 kb correction-segment floor. D202's 25 kb reverse interval survives that
floor and is labelled as smoothed by the objective. Their separate 15 kb
biological inversions also appear as orientation discrepancies. Explanatory
labels expose why a review question exists; they do not establish whether
its cause is an assembly error or target biology.

D207's automatic plan remains the biological transition at 300000. Its first
review alternative is the true join at 200000, with objective 150000 versus
149830 for the best automatic partition: a difference of 170 identity-weighted
query-base equivalents. This is the existing additive discordance plus 50000
per cut, not a probability or a new threshold. The other four alternatives
add biological questions; one also includes 200000. A single alternative is
credited by the one-to-one matching rule.

In D206, sixteen input PAF rows collapse to six canonical blocks and two
orientation intervals plus one alternative, the same three geometries as
D201. At the true error, four redundant/nested supporting rows have additive
span 30994 bp but union span 7998 bp. They share the mapping source and do not
provide four independent observations. The unchanged additive correction
objective remains sensitive to duplication: the biological discordance rises
from 15000 to 45000 and the alternative's objective difference falls from
85000 to 55000. Automatic output happens to remain unchanged in this case;
this is not general proof of duplicate-invariant scoring.

All seven successful primary contigs stay below the 128-block bound (one to
six canonical and one to six merged correction blocks). No orientation event
is omitted by the 25-event cap. D207 emits five of twelve distinct alternatives
found by the bounded search; the others are not shown as queue items. The
24-path-per-prefix search is not exhaustive. [Bounds and frame states](results/candidate-bounds.tsv)
retain these limits; D208 supplies no valid test of conflicting coordinate frames.

The new structured events provide read evidence at each actual proposed
coordinate and explicitly unavailable graph evidence. On D201 the error edges
have zero qualifying spanning molecules, while the biological edges have
seven and nine. On D207 the true join has zero, while the automatic biological
cut has twelve. These are observations from the declared idealized simulator,
not validated biological classifiers. [Coordinate records](results/candidate-coordinate-evidence.tsv)
include the exact assembly/read-evidence hashes and dependency group.
[Independent context verification](context-verification.json) checks 74 input
records, 42 mapping blocks, 29 structured events and 60 coordinate records over
the seven primary and three known scans. Identity checks do not certify stage
labels as biological truth.

## Automatic output and the declared replay

For each automatic read policy separately, both versions completely preserve
all three already-correct primary controls but repair **0/4** evaluable error
cases. They localize 0/7 error cuts. The report policy makes one false cut at
D207's biological transition and breaks one preservation control; the review
policy makes no false cuts and breaks no controls. Every successful output
retains each input base exactly once. These outcomes and FASTA/ledger bytes
are identical across versions. [Pairwise equality](results/paired-automatic-equality.json)
covers fourteen primary and six known pairs; two failed D208 policy pairs
have null equality. [Automatic scores](results/automatic-scores.jsonl) preserve
the unchanged provenance and complete-target metrics.

The only prespecified selection replay is D201. Baseline has no matching
interval, so its replay is unavailable without a fallback. rc4's matching
interval is selected using evaluator-known geometry (**oracle event selection**),
then its emitted `[70001,77999)` endpoints are used without refinement.
The linked explicit reverse edit is acknowledged, reviewed, previewed, applied,
replayed byte-for-byte, freshly aligned and validated. All 300000 input bases
remain exactly once and the biological control remains continuous.

The result is **not a complete target path under the frozen provenance scorer**.
Both emitted edges are 1 bp inside the declared `[70000,78000)` error, leaving
four target-relation discrepancies at the two ends. Boundary localization is
2/2 within tolerance with no false cuts, while strict complete structure is
false. Successful workflow validation certifies identity/component
reconstruction, not biological correctness or full repair.
[Replay scores and the selected event](results/scripted-selection.json),
[explicit edits](results/artifacts/candidate/D201/scripted_selection/reviewed-edits.tsv)
and [portable recipe](results/artifacts/candidate/D201/scripted_selection/recipe.json)
preserve that distinction. This is a feasibility diagnostic, not unattended
curation or a blinded researcher study.

**Post-score sequence interpretation diagnostic.** After the primary scores
were known, the coordinator requested a direct nucleotide comparison of the
existing replay with the frozen target. Independently comparing the original
committed FASTA, public-ledger reconstruction and target confirms that all
three **300000 bp sequences are exactly equal, with zero differing positions**.
The common sequence-only SHA-256 is
`fcc9fc97b5bab8025bdb4baa75f30806eaf4cba55e0b67427c99e8e1d49dfa00`.
The declared error interval begins with G and ends with C, a reverse-complement
pair: reversing the emitted interior can restore the exact target DNA while
leaving a different provenance assignment at those bases. Thus the failed
strict path score does **not** imply residual nucleotide error in this replay.
These edge bases alone cannot distinguish the two reversal extents, although
the simulated edit coordinates are known.

[The separately labeled diagnostic](post-score-sequence-diagnostic.json)
records the input, target, original output and ledger identities, both sequence
hashes and the unchanged false complete-structure score. This diagnostic was
not prospectively registered and does not replace the frozen success criterion.
No case, candidate output, proposal score, replay selection or emitted endpoint
changed. The oracle-selection limitation still applies.

## Known diagnostics, reported separately

| Diagnostic | Baseline → rc4 queue | Error coverage | rc4 evidence |
|---|---:|---:|---|
| D101 | 1 → 6 | 0/1 → 1/1 | First alternative at 200000, 95 weighted-bp objective difference |
| D012 | 1 → 4 | 0/1 → 1/1 | Exact `[151000,176000)` orientation interval, with opposite-orientation overlap exposed |
| D108 | 0 → 3 | 0/1 → 1/1 | Exact `[151000,159000)` interval below the correction floor |

Known-case edge visibility rises from 0/5 to 5/5, biological-burden events from
2 to 9, and biological regions inspected from 2/8 to 6/8. Both versions have
zero clean-region burden and zero redundant geometries under the frozen rule.
All six successful automatic pairs are identical, with 0/3 complete repairs
per policy. These previously inspected failures are diagnostic confirmations
and do not enlarge the primary sample.

## Reproduction and remaining limits

[Portable verification](portable-verification.json) reproduces all twenty
successful proposal scores, forty automatic scores and the one replay score
from regenerated inputs and 165 small artifacts. The two unavailable review
records and four unavailable automatic records retain null scores. The
component ledgers reconstruct sequence/provenance metrics; original output
FASTA byte hashes are recorded separately because headers and wrapping may
differ during reconstruction. [The isolated clean-copy check](clean-copy-verification.json)
also reproduces the four summary exports, checks all 31 analytic tests and
retains the earlier fixtures and 294 historical score records. Its additional
post-score diagnostic check reproduces the exact nucleotide comparison
separately from the frozen primary scores.

See [reproduction commands and software identities](README.md#reproduction).
Bulk output commitments retain original logs and figures by hash; those bulk
files and software archives are not redistributed here. Local subprocess
elapsed sums are recorded in [resource use](results/resource-use.json): scan
time was 1.665 s baseline versus 5.455 s candidate over eleven attempted cases,
with extra work enabled in the candidate. This is not a comparative speed
benchmark.

The evidence supports improved proposal availability on these admitted
development scenarios, with a measured increase in biological inspection
burden. It supplies no automatic complete-repair gain, no valid conflicting-
frame benchmark result, no human-review accuracy estimate, no raw-read
platform calibration and no independent-genome generalization claim. Further
evaluation would require a separately authorized, prospectively fixed study.
