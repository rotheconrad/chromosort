---
title: Read continuity and gap inspection
---

# Read continuity and gap inspection

The rc2 refinement adds two conservative review policies motivated by the frozen
development evaluation of rc1. It preserves the rc1 artifacts and parameters for
comparison. Neither policy identifies biological truth or guarantees repair.

## Defer a proposed cut when reads support continuity

For manifest-based automatic `fix`, a proposed boundary triggers review when at
least two distinct read names have one continuous aligned CIGAR block covering
1000 bp on both sides, after MAPQ 20 filtering and exclusion of secondary rows.
Supplementary segments remain eligible, but multiple segments of the same
molecule count once. CIGAR deletions/skips cannot bridge the boundary. PAF without
CIGAR is unavailable for this safeguard; span-only coverage is insufficient.

The default `--read-continuity-policy review` defers the entire contig plan if any
boundary conflicts with continuity. This preserves its original sequence and
orientation in the normal full-assembly output. Automatically merging a subset
of oriented pieces would invent a replacement plan. Selective interventions
belong in explicit review. The original proposed boundaries and their evidence
remain in sidecars and evaluation tables; all accept flags for a deferred plan
are cleared. `workflow scan` retains the proposals with pending decisions.

`--read-continuity-policy report` records the same evidence without altering the
algorithm's cuts. A fully reviewed decision table or workflow recipe can override
the advisory deferral. Legacy non-manifest correction commands remain unchanged.
The configurable thresholds are `--continuity-min-mapq`, `--continuity-anchor-bp`
and `--continuity-min-spanning`; they are operational review settings, not
validated confidence levels.

### Evidence selection

`--read-evidence-id` explicitly selects one manifest PAF, SAM or BAM source.
Without it, a single source is used. Multiple representations are treated as
correlated only when they share a nonempty `dependency_group` or `readset_id`;
then portable PAF is preferred, followed by SAM and BAM, with evidence ID as a
stable tie-breaker. They are never pooled to reach the molecule threshold.
Different source groups require explicit selection. Missing PAF CIGAR does not
silently trigger a fallback to a differently represented alignment.

The selected source's assembly digest, target names and lengths must match the
current stage. The current repository validates target names and lengths before
read filtering, including SAM/BAM sequence headers; incompatible evidence fails
instead of appearing as absent support.
Reports retain the declared read sample and dependence metadata. Unknown sample
identity remains unknown. Continuity from cross-sample reads can motivate review
but does not establish the target genotype. Simulated reads, repeated/chimeric
reads and mapper errors also limit interpretation.

**No spanning reads is neutral.** Low coverage, Ns, sampling rules or an omitted
alignment can explain absence. The safeguard never uses absence to add a cut or
increase confidence that an error exists. Without usable read evidence, existing
reference-only automatic behavior remains; use supervised review when its
interventions cannot be justified independently.

## Inspect smaller gaps without weakening correction smoothing

`workflow scan --inspect-min-gap-bp 1000` exposes uncovered query intervals between
backbone-aligned flanks at least `--min-segment-bp` long. The 1 kb floor is separate
from the existing alignment-block and emitted-piece thresholds. It expands the
review queue; it does not change segmentation scores or propose additional cuts.
Small genuine variants can therefore remain smoothed in the automatic planner
while their unexplained alignment geometry is inspectable. Queue size can be
controlled with the inspection floor or a targeted `--contigs` list.

Each inspect item produces read panels at **both interval boundaries**. A read
can span the interior of an inverted segment without bridging its junctions;
midpoint continuity is not boundary continuity. The table's center position is
only a locator. Accepting an inspect item never edits sequence. Adding a cut or
reorientation requires a separate explicit manual decision and fresh alignment.

## Development observations and residual limits

Before the refinement, D010's proposed biological-exchange cuts had 19 and 9
spanning molecules under the stated filters; D011's misjoin proposal had none.
This motivates deferral, not a general ability to classify exchanges or errors.
D001 also lacks bridges because of retained Ns and simulator sampling rules;
that absence is not independent confirmation of its error.

D003's roughly 8 kb unexplained interval was excluded by rc1's 10 kb inspection
floor. Inspecting both of its edges can expose different read continuity from
the larger simulated biological inversion, but does not determine the correct
repair or prove either interval's origin. The correction algorithm remains
unable to reconstruct its missing reference-alignment blocks automatically.

D012 injects an orientation-only event in a native-repeat background; the
observed mapping does not localize its inversion boundaries. It does not inject
a collapse. A continuity conflict can defer other cuts in the same plan,
including potentially useful ones. D014's real target-specific sequence remains
a preservation case; new review evidence does not infer missing copies, dosage
or phasing. Broader accuracy, review workload and generalization require separate
evaluation. All policy choices here were informed by development data; held-out
inputs/truth were not used for refinement.

## Paired development commands

```bash
chromo fix --manifest CASE/manifest.json --all --mode MODE \
  --read-evidence-id reads_sim01 --read-continuity-policy review \
  --continuity-min-mapq 20 --continuity-anchor-bp 1000 \
  --continuity-min-spanning 2 --output-fasta OUTPUT/output.fa \
  --report OUTPUT/fix.tsv --agp OUTPUT/output.agp
# Repeat in a new directory with --read-continuity-policy report.
chromo workflow scan --manifest CASE/manifest.json --mode MODE \
  --read-evidence-id reads_sim01 --inspect-min-gap-bp 1000 \
  --output-dir OUTPUT/review
```

Use each of the existing conservative/comprehensive/sensitive modes without
changing their smoothing settings or adding orientation normalization. Keep
automatic correction, review-only packets and scripted replay as separate arms.

## Reproducibility record

See the [public development results](https://github.com/rotheconrad/chromosort/blob/main/benchmarks/assembly-errors/docs/results.md)
for paired rc1/rc2 measurements, artifact identities and residual failures.
The frozen rc2 archive predates the public-readiness input-contract checks and
this documentation correction: D012 is an orientation-only injection in a
native-repeat background, not an injected collapse. Rebuilds of the current
repository must be identified separately from that frozen archive.
