---
title: Aligned orientation and competing-plan review
---

# Aligned orientation and competing-plan review

ChromoSort 0.4.0rc4 adds two optional inspection paths to `workflow scan`:
aligned internal orientation discrepancies and competing segmentation plans.
Both produce pending `inspect` events. Acknowledging an inspection makes no
sequence change. An inversion may be genuine variation; the event is evidence
to examine, not an inferred assembly error.

```bash
chromo workflow scan --manifest inputs.json --output-dir review \
  --mode comprehensive --inspect-orientation \
  --inspect-min-orientation-bp 1000 --alternative-plans 5 \
  --review-max-blocks 128
```

The new flags default to off (`--inspect-orientation` absent and
`--alternative-plans 0`). Automatic correction uses its existing defaults,
alignment weights, smoothing, support checks and whole-contig read-continuity
policy. The existing [unaligned-gap inspection]({{ '/review-refinement/' | relative_url }})
remains available. None of these inspection paths writes a corrected FASTA.

Open `review/index.html` for decisions and native read panels, and follow
**Orientation and alternative-plan evidence** to `review-proposals.html`.
That report shows coordinates, the read source and filters, graph availability,
frame context, objective components, alignment blocks and exclusion reasons.
`review_proposals.json` is the complete machine-readable record;
`review_proposals.tsv` is an event table with nested values encoded as JSON.

## Aligned orientation inspection

The input is the manifest's primary alignment representation for each
reference. Existing name, length, coordinate and digest validation applies.
Passing primary mappings use the scan's MAPQ threshold. Orientation candidates
use sequences explicitly labelled `backbone`, before the correction segment
size floor, cross-reference withholding and breakpoint smoothing. Sequences not labelled `backbone` remain context and do not generate events,
including the usual supporting and alternative references. The sequence-level
backbone flag is authoritative, as it is for correction. This review eligibility is deliberately broader than
correction eligibility, which is recorded separately for every block.

Each namespaced reference sequence is a separate coordinate frame. The
inspection procedure is:

1. Coalesce identical alignment rows into a block with its original row count
   and evidence identities. Keep overlapping nonidentical rows as separate
   evidence blocks.
2. Determine the strands at the earliest mapped start and latest mapped end
   within that frame. A frame has an anchor orientation only if both terminal
   strand sets contain the same single strand. Mixed terminal strands remain
   explicitly ambiguous in the contig ledger and do not generate an interval.
3. Union overlapping or touching spans on the contrary strand. Retain intervals
   at least `--inspect-min-orientation-bp` long with anchor-strand evidence
   extending to the native left and right. The flank evidence may overlap the
   interval; such overlap is shown and flagged, not treated as independent
   confirmation. No maximum flank-distance rule is imposed.
4. Combine identical native intervals across frames into one event, preserving
   each frame's anchor and interval strands and supporting/opposing blocks.
   Nonidentical intervals remain separate; there is no fuzzy cross-frame merge.

A uniformly reverse-aligned record has `whole_record_orientation` context and
no internal orientation event. Whole-record reference orientation remains a
separate explicit operation. Terminally ambiguous, terminal discrepancies,
unmapped changes and signals below the review floor can therefore be missed.
Events on genuine inversions add review work and must not be counted as false
automatic corrections.

Every event reports both native edges, frame identities, supporting union
length, additive support length and original row count. Overlapping opposite
strands and competing reference mappings are explicit. Identical or nested
rows do not multiply an interval event, but source labels and row names cannot
prove biological independence. The underlying automatic planner's additive
alignment weights are unchanged and can still count redundant overlap.

Stage explanations include `below_correction_segment_floor`,
`withheld_cross_reference_overlap`, `orientation_not_used_in_mode`,
`smoothed_by_objective`, terminal merging and group support where applicable.
The contig ledger retains the automatic scan proposals and the best partition
before guards. Exact agreement at both emitted automatic boundary positions
is recorded separately; it is not a tolerance-based accuracy score. A mapped
orientation interval can be absent from gap inspection because its bases are
covered by the alignment union.

## Competing segmentation plans

Alternative plans start from the **same merged correction blocks** as the
automatic planner, after its segment and cross-reference eligibility filters.
No extra short blocks are inserted and no overlap weights are changed. A block
signature is its namespaced reference sequence and, when enabled by the
selected mode, its orientation. In conservative and chromosome modes,
same-reference orientation alternatives may therefore be absent. Sensitive
mode does not optimize the smoothing objective; its alternative-objective
comparison is explicitly marked as not comparable to automatic optimization.

For a group `G` of consecutive blocks, let `W(s,G)` be the sum of
`identity_percent × aligned_query_bp / 100` for blocks with signature `s`.
The discordance cost is:

```text
D(G) = sum over signatures s of W(s,G) - max over s of W(s,G)
J(G1,...,Gk) = sum over groups i of D(Gi) + lambda × (k - 1)
```

`lambda` is the existing breakpoint penalty (50,000 by default). `J` is an
additive cost in identity-weighted aligned-query-base equivalents. Overlapping
alignments can contribute more than their native union length. Neither `J`
nor a difference between two values is a probability, confidence interval or
error likelihood. No objective-difference threshold is used to admit a plan.

The implementation extends the interval dynamic program to keep at most
`K = 4 × (N + 1)` paths at each prefix, where `N = --alternative-plans`.
An edge joins a prefix to a consecutive block group and contributes `D(G)`
plus one breakpoint penalty if the group is not first. Paths are ordered by
cost, cut count and the tuple of block-boundary indices. Adjacent groups with
the same dominant signature are then merged by the existing planner helper.
The objective is recalculated for the merged groups.

A native boundary is the integer midpoint of the last block end in the left
group and first block start in the right group, matching the smoothed planner's
boundary rule. Partitions with noninternal, duplicate or nonincreasing native
boundaries are discarded and counted. Zero-cut partitions remain objective
context, not events. The exact existing scan cut set is excluded. Identical
native cut sets are combined, with deterministic cost/signature/index ties.
The best `N` remaining plans are reported by cost, cut count and coordinates.
The bounded prefix search followed by merging and deduplication can return
fewer than `N` events; it is not an exhaustive ranking of all possible native
segmentations.

Each event contains the group signatures, block IDs, discordance, penalty,
total objective and difference to the best unconstrained block partition.
Group support and minimum emitted-piece failures are shown. Alternative
generation precedes terminal merging, breakpoint guards and read-continuity
deferral; an alternative may fail those automatic checks. These plans are
inspection choices, not additional accepted automatic cuts. A plan can have
several boundaries, all of which need review.

## Bounds and evidence availability

`--alternative-plans` accepts 0–20 and `--review-max-blocks` accepts 1–256.
The default cap is 128. Orientation inspection requires both canonical mapping
and merged correction block counts to fit the cap, so its stage trace is also
bounded. Alternative search uses the merged correction block cap. Exceeding
a cap skips that inspection path for the contig, records the counts and
`skipped_block_cap`, and leaves existing automatic inputs intact. Orientation
events are sorted by native start/end and capped at 25 per contig, with the
omitted count recorded. The registry retains the mapping blocks even when
inspection is skipped.

For `n` merged blocks, interval costs use the existing group-support function,
requiring O(n³) block visits in this implementation. The search examines O(K n²) path extensions; each can copy an O(n)
boundary tuple, and each prefix sorts up to O(K n) candidates. Including
worst-case tuple comparisons gives an O(K n³ log(K n)) time bound. Retained
path tuples add up to O(K n²) space. Mapping validation and the existing correction eligibility checks occur
before these caps. Reading evidence and producing figures can dominate scan
time, particularly when many events share coordinates.

Read evidence is assessed at **every proposed native boundary**, using the
existing selected-source, MAPQ, CIGAR-contiguous anchor and unique-read-name
rules. The machine record uses the `continuity-*` thresholds. Read panels use
the existing `read-*` panel settings, which are recorded with each figure;
keep their anchors consistent when comparing counts. Independent read sources
are not silently pooled. Missing reads, missing CIGAR and insufficient spans
are distinguishable. Duplicate names count once, while renamed dependent
alignments may still overstate molecule independence. Lack of spanning reads
does not establish an error. Inspection evidence always uses report policy and
does not invoke automatic deferral or cutting.

With one GFA source, graph paths or direct segments are projected independently
for each contig. A projection must have the exact native contig length.
Boundary `p` shows graph context for both 1-based bases `p` and `p+1`.
Missing, conflicting-source or length-incompatible projections are unavailable;
context from another proposed boundary is never substituted. Graph degree and
proximity are coordinate context, not confirmation of an assembly error.
Sequence identity and path correctness still depend on the supplied graph's
provenance; matching names and length do not establish them.

## Event and edit contract

The root schema is `chromosort-review-proposals-v1`, containing
`assembly_sha256`, `manifest_sha256`, `inputs`, `parameters`, `limits`,
`objective_definition`, `blocks`, `contigs` and `events`.

| Field | Meaning |
| --- | --- |
| `event_id`, `kind`, `target` | Stable stage/geometry ID; `orientation_discrepancy` or `alternative_segmentation`; native contig name. |
| `start0`, `end0`, `positions` | Orientation interval `[start0,end0)` and both edges; for segmentation, all nonempty cut points with start/end equal to the first/last point. A one-cut plan has equal start/end. |
| `action`, `decision`, `implicit_edits` | Always `inspect`, `pending`, `false` in the original proposal. |
| `supporting_block_ids`, `opposing_block_ids`, `overlapping_block_ids` | References into the block registry; inspect frame-local membership and overlaps before interpreting support. |
| `reference_frames`, `frame_contexts` | Namespaced reference sequences and, for orientation intervals, separate strand/flank contexts. |
| `loss_stages`, `reasons`, `objective` | Stage explanations, readable rationale and objective terms (null for orientation events). |
| `coordinate_evidence` | One entry per `position`, with `reads` and `graph` results, source identities and assembly digest. |
| `figures` | Native boundary positions and relative paths to SVG/PDF/PNG figure bundles. |

Block records use `block_id`, native and reference intervals, frame, strand,
identity, MAPQ, reference metadata/digest, correction eligibility, exclusion
stages, `row_count` and source evidence identities. Stable event IDs depend on
assembly bytes, target, kind and native positions, so identical geometry is
stable across row order or redundant mappings. IDs alone do not identify the
evidence revision: `scan.json` binds the detailed JSON by SHA-256 and explicit
edit IDs bind the complete scan digest. Replacing the details invalidates
subsequent edit/apply operations. Keep the source manifest and its inputs.

The same event IDs enter the existing `decisions.tsv` with action `inspect`.
The decision table's point field is an interval midpoint for display; use
`positions` and `start0/end0` from the detailed record for deliberate edits.
Accepting an inspection acknowledges it only. Prepare a chosen interval with
`workflow edit --event-id ID --action reverse --start START0 --end END0`,
or prepare each desired point with `--action cut --start P --end P`.
Supply a rationale, then record explicit accepted/rejected/deferred edits and
a reviewer before `workflow apply`. If a cut already occurs in the primary
scan decision table, resolve it there; duplicate or contradictory edits fail.

See the [explicit edit and replay commands]({{ '/commands/workflow/' | relative_url }})
for preview, apply and fresh-stage validation. Every source base must remain
accounted for exactly once, accepted reverse intervals must be disjoint, and
accepted cuts cannot lie inside them. Changed FASTA always requires new
alignments before further alignment-dependent interpretation.

## Implementation and validation boundaries

`review_proposals.py` owns candidate generation, objective tracing, evidence
collection and detailed reports. `workflow.py` connects those records to the
existing decision and figure flow. `manifest.py` can return provenance alongside
normalized segments without changing its default iterator. `fix_contigs.py`
supplies the unchanged automatic blocks, signatures, objective and boundary
helpers. `continuity.py`, `graph.py` and `readplot.py` supply stage-bound evidence.

Small independent fixtures in `tests/test_review_proposals.py` check short and
smoothed intervals, clean/whole-record controls, ambiguous anchors, overlapping
and duplicated rows, legal multiple frames, exact boundary evidence, unavailable
projections, bounded search against exhaustive small partitions, event caps,
acknowledgment, explicit edits, replay and stale-detail rejection. They establish
implementation behavior, not biological accuracy or reviewer effort.
Measured proposal coverage, added burden, admission failures and automatic
output compatibility belong to the separately frozen development evaluation.
Interactive browser QA remains a release check; static figures, generated HTML
and JavaScript syntax checks do not establish interactive correctness.
