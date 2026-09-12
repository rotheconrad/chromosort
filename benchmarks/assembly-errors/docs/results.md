# Recorded development results

These records cover 14 development scenarios derived from a small set of public
wheat/barley loci. The repository includes **294 original development scores**:
154 pilot outputs, 56 rc1 outputs and 84 rc2 outputs. Held-out data, scores,
seeds and evaluator material are excluded. Source records, scenarios and
boundaries are dependent; counts below are descriptive mechanism tests, not
genome-wide accuracy estimates.

The [full score ledger](../reports/historical_scores.jsonl) retains each score,
the selected run identity and the original score/run hashes. The original score
JSON hash can be reconstructed from each entry. `scripts/summarize_history.py`
checks all 294 score identities and generates [summary.tsv](../reports/summary.tsv).
The ledger is compact result evidence, not a substitute for unavailable bulk
historical output/software archives; exact replay limits are in the
[reproduction guide](reproduction.md).

## Bounded rc2 refinement

The [rc2 protocol](../configs/candidate-0.4.0rc2.json) pairs the same cut planner
with two read-evidence policies. `review` defers the whole contig plan when a
proposed boundary has at least two distinct CIGAR-contiguous spanning molecules
with 1 kb anchors and MAPQ ≥20. `report` records the evidence but emits the
algorithm's cuts. Both are automatic arms. Neither represents completed expert
or agent curation; missing spans are neutral.

| Policy | Mode | New boundaries | Localized error boundaries | False boundaries |
| --- | --- | ---: | ---: | ---: |
| report | conservative | 3 | 1 / 9 | 2 |
| review | conservative | 1 | 1 / 9 | 0 |
| report | comprehensive | 5 | 3 / 9 | 2 |
| review | comprehensive | 3 | 3 / 9 | 0 |
| report | sensitive | 8 | 3 / 9 | 5 |
| review | sensitive | 3 | 3 / 9 | 0 |

All 84 rc2 outputs retained every input base exactly once. All 42 report-policy
FASTAs and scores are identical to the corresponding rc1 manifest-mode results.
The paired comparison therefore retains the same manifest/alignment semantics.
Review policy removed the observed false cuts without losing a localized error
boundary in this development cohort. It introduced no new preservation-control
fragmentation relative to the unedited input.

The threshold and inspection changes were motivated by these development cases.
This is explicitly **development-informed refinement**, not independent-cohort
validation. Most error boundaries remain unlocalized: eight of nine under
conservative, six under comprehensive/sensitive. All arms retain ten target-path
discrepancies, 84,420 missing callable target bases and 10,000 extra biological
copy bases across the panel. Those defects include existing input problems;
they differ from input bases lost/duplicated by the tool. Cut localization is
not complete repair, reorientation, dosage recovery or restored continuity.

Specific findings:

* D010's two simulated biological exchange junctions have 19 and 9 spanning
  molecules. Review policy preserves both chromosomes; report policy makes
  two false cuts. Whole-contig deferral can suppress useful cuts elsewhere in
  a mixed plan, a tradeoff not resolved by these examples.
* D011's misjoin remains localized at 250,000 with no read-continuity conflict.
  Absence of spans did not add a cut or independently confirm an error.
* D001 comprehensive/sensitive localize both operational boundaries within
  500 bp. Ns and read sampling explain absent bridges. The published optical
  maps, not simulated lack of bridging, support the original error claim.
* D003 remains unrepaired. The smaller [151,002, 158,998) alignment gap is now
  inspectable at a separate 1 kb floor, with a panel at each edge. The two edges
  show zero spanning and one split molecule each. A distinct simulated
  biological inversion's edges show 10 and 7 spanning molecules. These are
  evidence for inspection, not inferred truth labels or automatic repair.
* D012's 25 kb orientation error remains unresolved. Sensitive's proposed
  155,618/175,252 cuts miss its 151,000/176,000 boundaries. Seven spans at the
  first proposal defer the whole plan, avoiding both false cuts. This is an
  orientation-only injected error on a native-repeat background, not collapse.
* D014 sensitive's false cut is deferred with four spanning molecules; the
  separate 15 kb target-specific fragment remains retained.

The original rc1 scans used BAM; rc2 scans explicitly used PAF. Both automatic
rc2 policies used the same PAF. Inspection, scripted replay, automatic correction
and future human/agent review must remain separate comparison categories.

## Original RagTag and pre-feature comparisons

The [v1](../configs/comparators.json) and [base-alignment amendment](../configs/comparators.v2.json)
preserve the actual comparator parameters. Every row below covers all 14
development cases, including preservation/clean controls.

| Automatic correction arm | New boundaries | Localized error boundaries | False boundaries |
| --- | ---: | ---: | ---: |
| RagTag default | 3 | 1 / 9 | 2 |
| RagTag default with reads | 1 | 0 / 9 | 1 |
| RagTag local-distance/read diagnostic | 1 | 0 / 9 | 1 |
| RagTag with base alignments | 9 | 4 / 9 | 5 |
| RagTag with base alignments and reads | 7 | 2 / 9 | 5 |
| ChromoSort ac4 conservative | 4 | 1 / 9 | 3 |
| ChromoSort ac4 comprehensive | 6 | 3 / 9 | 3 |
| ChromoSort ac4 sensitive | 8 | 3 / 9 | 5 |

The ac4 commit predates manifest-aware/reference-label/review features. Its
flattened per-reference PAF union has different input semantics from rc1/rc2
manifest correction. rc1 conservative/comprehensive avoid D009's one false
cut; this cannot isolate an algorithm effect from the changed evidence rules.

RagTag defaults omit minimap2 base alignment. On D003 a default approximate
chain spans rearranged intervals, while `-c` splits it. Therefore paired
base-alignment arms were added in a documented development amendment; original
defaults remain separate. The pinned RagTag parser rejects `--secondary=no`
through its SAM-flag check, so the supported parameters are `-x asm5 -c -t 2`.
Secondary/MAPQ behavior and joint mapping therefore differ from the independently
mapped reference PAFs. No identical-alignment-set claim is made.

RagTag read-validation arms receive the same underlying simulated target reads
as ChromoSort, but each tool uses its declared interpretation. RagTag's local
distance/window arm is a short-locus diagnostic, not a genome-optimized policy.
The rc2 two-molecule review safeguard is not equivalent to RagTag correction
validation. Read-assisted and reference-only rows must be labeled separately.

Fixed/inferred-gap scaffold arms operate on original inputs and appear separately
in the ledger; do not pool them with correction. D006 retain-all scaffolding
keeps both overlapping fragments and their copy excess. Correct overlap trimming
requires an explicit discard policy beyond raw input-retention scoring.
Multi-reference concatenation is a diagnostic, not RagTag merge or an inferred
subgenome-labeling study.

Earlier rejected-flag and memory-wrapper attempts are excluded from the scored
record set. They were harness compatibility failures, not biological failures
attributed to RagTag. Per-run timings describe tiny local loci and exclude
shared evidence preparation; they do not support genome-scale speed claims.

## Limits that survive this checkpoint

The read model is iid substitutions, not calibrated PacBio/ONT. Actual tandem
copies are distinguished from simulated inversion/exchange alleles. Only the
wheat orientation has the specific reported reconstruction status described in
the [source audit](sources.md); exact historical barley edits are unrecovered.
There is no read-assembled graph, validated donor-fill scoring, completed
balanced blind review, empirical review-time benefit, or demonstrated overall
superiority in this record. Independent new loci and controlled validation are
required before extending the claims.

## Additive second development panel

The [second-panel report](../cycle2/results.md) preserves these historical
records and adds a paired comparison of the published `dd66e3d` baseline with
immutable 0.4.0rc3 on eight new-locus scenarios. It separately records automatic
correction, pending inspection and supplied-coordinate repair, with a stricter
complete-target-path metric. Automatic outputs are unchanged between versions;
the explicit workflow completes all supplied repairs. This is another
development experiment, not an unseen validation.
