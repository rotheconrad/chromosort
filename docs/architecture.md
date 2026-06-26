---
title: Architecture
description: ChromoSort software architecture, formal data model, algorithms, and validation strategy.
math: true
mermaid: true
---

# Architecture

ChromoSort is a reference-guided genome assembly curation toolkit. It consumes
existing whole-genome alignments, assembly graph records, and optional
cross-evidence tables; it does not assemble reads, run aligners, polish bases,
or infer variants. Its core responsibility is to convert heterogeneous evidence
into deterministic, auditable contig-level decisions: order, split, cut,
scaffold, review, plot, and optionally fill graph-supported gaps.

The implementation is intentionally organized around command boundaries. Each
command writes explicit FASTA, TSV, HTML, or plot artifacts rather than sharing
hidden state. This makes the pipeline inspectable by users and reviewers, and
it makes sequence-changing steps depend on documented inputs and reports.

## Design Principles

1. **Evidence first, sequence edits second.** Alignment metrics, graph context,
   gap decisions, and manual-review choices are written as reports. Sequence
   changes occur only in commands whose purpose is to write a new FASTA:
   `sort`, `clean`, `fix`, `cut`, `manual apply`, `scaffold`, and
   `gapfill --apply`.

2. **Conservative defaults.** The default behavior favors retaining ambiguous
   evidence for review over making irreversible edits. Examples include
   split-candidate protection in `sort`, breakpoint smoothing in `fix`,
   report-only graph evidence in `sort`, `fix`, and `scaffold`, and
   ambiguous-path refusal in `gapfill`.

3. **Simple, explicit data models.** FASTA records, alignment segments,
   merged intervals, assignment rows, graph nodes/edges, split pieces, scaffold
   gaps, and fill plans are represented with small dataclasses or TSV schemas.
   The project avoids an internal database and uses files as the public
   contract between stages.

4. **Format normalization before decision logic.** MUMmer `show-coords` and
   minimap2 PAF records are normalized into a shared `Segment` representation
   before coverage, identity, orientation, and ordering decisions are made.
   GFA graph support is similarly normalized into oriented node and link
   records.

5. **Graph evidence is guarded.** GFA, GAF, Hi-C-like contact tables, and
   reference-placement PAF can refine review and gap-filling decisions, but
   sequence insertion requires a validated graph path. Complex or unknown GFA
   overlaps are preserved as uncertainty rather than converted into trim
   lengths.

6. **Reviewer-readable provenance.** Important decisions are surfaced as status
   codes, reason strings, support counts, interval lengths, risk flags, and
   summary tables. The goal is to make each accepted or rejected decision
   reconstructable from input evidence.

## System Map

<figure class="method-diagram">
<div class="mermaid">
flowchart TD
  A[FASTA] --> B[Alignments]
  C[Graph] --> D[Review]
  B --> D
  B --> E[Edit]
  B --> F[Sort]
  C --> F
  F --> G[OrderedFASTA]
  G --> H[Scaffold]
  C --> H
  H --> I[ScaffoldFASTA]
  G --> J[Gapfill]
  C --> J
  D --> J
  J --> K[GapfilledFASTA]
</div>
</figure>

The commands can be used as a pipeline, but the architecture does not require
a single prescribed run order. Users can stop at any report, inspect it, edit
manual inputs, re-align edited FASTA records, and resume with later commands.

## Command/API Responsibilities

The public command surface is the `chromo` dispatcher in
`src/chromosort/cli.py`. `pyproject.toml` also exposes compatibility scripts:
`chromosort`, `chromosort-clean`, `chromosort-fix-contigs`,
`chromosort-scaffold`, and `chromosort-gapfill`.

| Command/API | Primary module | Main inputs | User-facing responsibility | Sequence-changing behavior |
| --- | --- | --- | --- | --- |
| `chromo sort` | `reference_order.py` | Reference FASTA, assembly FASTA, coords or PAF, optional GFA | Assign contigs to reference sequences, classify low-confidence and duplicate-overlap records, order kept contigs by reference position, and write match/assignment summaries. | Writes an ordered FASTA for kept contigs. Optionally reverse-complements records with `--orient-to-reference`. |
| `chromo clean` | `clean.py` | Reference FASTA, raw assembly FASTA, coords or PAF | Apply sort assignment/filtering to raw contigs, select retained raw contigs for conservative fix planning, reconcile retained unsplit contigs and split pieces, and write unified audit reports. | Writes a cleaned FASTA. Optionally reverse-complements emitted records with `--orient-to-reference`. |
| `chromo eval fix` | `evaluate.py` and `fix_contigs.py` | Assembly FASTA, coords or PAF, selected contigs or `--all`, optional GFA, long-read PAF, and GAF | Prepare an editable review-event TSV for candidate split pieces using the fix planner plus optional graph and read evidence. | Does not change FASTA. Accepted rows can later be applied by `chromo fix --reviewed-plan`. |
| `chromo eval scaffold` | `evaluate.py` and `scaffold.py` | Ordered FASTA, assignment TSV, optional GFA, long-read PAF, and GAF | Prepare an editable review-event TSV for adjacent scaffold junctions, inferred gaps/overlaps, graph context, GAF path support, and read-bridge support. | Does not change FASTA. Accepted rows can later pin gap lengths through `chromo scaffold --reviewed-plan`. |
| `chromo eval gapfill` | `evaluate.py` and `gapfill.py` | Ordered FASTA, assignment TSV, GFA, optional GAF, Hi-C table, reference-placement PAF, and long-read PAF | Prepare an editable review-event TSV from the gapfill planner, including candidate path status, fillability, branch-support evidence, and read-bridge context. | Does not change FASTA. Accepted rows can later restrict `chromo gapfill --reviewed-plan --apply`. |
| `chromo fix` | `fix_contigs.py` | Assembly FASTA, coords or PAF, selected contigs or `--all`, optional GFA | Detect and split selected chimeric contigs into reference-labeled pieces using mode-specific breakpoint planning. | Writes a fixed FASTA and split report. Optionally writes only pieces and/or reverse-complements pieces. |
| `chromo cut` | `cut.py` | Assembly FASTA plus explicit cut positions | Apply reviewed 1-based cut positions and report emitted pieces. | Writes a cut FASTA. No cut is inferred from alignments. |
| `chromo manual` | `manual.py` | Reference FASTA, assembly FASTA, coords or PAF, optional GFA, long-read PAF, GAF, and review table | Build a self-contained review dashboard with alignment, contig, graph, initial piece information, and optional evidence panels. | Does not change FASTA. |
| `chromo manual fix/scaffold/gapfill` | `manual.py` and `review.py` | Same dashboard inputs plus task-specific `chromo eval` review table | Add a focused task-event queue and selected-event evidence panels for fix, scaffold, or gapfill review. | Does not change FASTA. Manual edits still require browser export plus `manual apply`; reviewed TSVs feed the matching executor. |
| `chromo manual apply` | `manual.py` | Assembly FASTA and exported recipe JSON | Apply a dashboard recipe containing slices, orientation choices, removals, and optional scaffold grouping. | Writes a manually edited FASTA and optional report. |
| `chromo scaffold` | `scaffold.py` | Ordered FASTA and `sort` assignment TSV, optional GFA, optional reviewed plan | Join ordered contigs per reference with inferred, fixed, or reviewed gaps; report reference overlaps and graph adjacency. | Writes scaffold FASTA. Overlap trimming occurs only under explicit overlap policies. |
| `chromo gapfill` | `gapfill.py` | Ordered FASTA, assignment TSV, GFA, optional GAF, Hi-C table, reference-placement PAF, reviewed plan | Plan graph-supported gap fills, annotate risk, resolve unique supported paths, and optionally apply fillable paths. | Writes gapfilled FASTA only with `--apply`; reviewed plans can further restrict which fillable gaps are applied. |
| `chromo plot` | `plot.py` | Reference FASTA, assembly FASTA, coords or PAF, optional assignment TSV, optional selected reference IDs, optional GFA overlay | Draw whole-genome, per-reference, and selected-reference dot plots from existing alignments, with optional projected query-axis graph boundaries. See the [dot-plot guide]({{ '/dot-plots/' | relative_url }}) for interpretation patterns. | Does not change FASTA. |
| `chromo graph-map` | `graph_map.py` and `graph.py` | Contig FASTA and unitig or contig GFA | Project GFA path/walk segment intervals into matching contig FASTA coordinates and write projection, summary, and warning reports. | Does not change FASTA. |

## Algorithm And Data Model Activation Map

This section makes explicit where each documented algorithm or data model is
used. "Activation" means the command path, mode, or parameter that calls the
logic; it does not imply that the evidence source is allowed to change
sequence by itself.

| Algorithm or data model | Activated by | Key parameter or mode controls | Main output fields or artifacts | Validation anchors |
| --- | --- | --- | --- | --- |
| Alignment normalization into `Segment` records | `sort`, `clean`, `fix`, `eval fix`, `manual`, and `plot` | `--coords` or `--paf`, `--min-segment-bp`, `--min-segment-idy`, `--min-mapq`, `--include-secondary-paf` | Match metrics, embedded alignment rows, dot-plot segments, fix blocks | `tests/test_reference_order.py`, `tests/test_fix_contigs.py`, `tests/test_plot.py`, `tests/test_manual.py` |
| Match aggregation and best-reference ranking | `sort`, `clean`, `manual`, and `plot --assignments` indirectly through sort tables | `--min-aligned-bp`, `--min-query-cov`, `--min-best-ref-share` | `match_report.tsv`, `contig_assignments.tsv`, manual query records | `tests/test_reference_order.py`, `tests/test_clean.py`, `tests/test_manual.py` |
| Duplicate-overlap filtering and terminal-extension rescue | `sort` and the sort/filter stage inside `clean` | `--no-overlap-filter`, `--overlap-mode`, `--min-novel-ref-bp`, `--min-novel-ref-frac`, `--novel-ref-criteria`, terminal-extension thresholds | Assignment statuses such as `duplicate_overlap`, `terminal_overlap`, and kept/rescued records | `tests/test_reference_order.py`, `tests/test_clean.py` |
| Split-candidate protection before fixing | `sort` and `clean` | `--split-candidate-min-aligned-bp`, `--split-candidate-min-query-frac`, `--split-candidate-max-best-share` | `kept_split_candidate` assignment status and split-candidate annotations | `tests/test_reference_order.py`, `tests/test_clean.py` |
| Fix breakpoint block planning | `fix`, `eval fix`, and the fix stage inside `clean` | `--mode sensitive/chromosome/comprehensive/conservative`, `--breakpoint-penalty-bp`, `--max-merge-gap`, piece support thresholds | `fix_report.tsv`, `fix_review.tsv`, cleaned FASTA pieces | `tests/test_fix_contigs.py`, `tests/test_eval.py`, `tests/test_clean.py` |
| Reviewed fix-row application | `fix --reviewed-plan` | Accepted `chromosort-review-event-v1` rows with task `fix` and action `split_piece` | Fixed FASTA and report rows generated from reviewed slices | `tests/test_eval.py`, `tests/test_fix_contigs.py` |
| Explicit cut planning | `cut` | `--cut`, `--cuts-file`, `--contig/--pos`, `--min-piece-bp` | Cut FASTA and cut report | `tests/test_cut.py` |
| GFA node and link evidence | `sort --gfa`, `fix --gfa`, `eval fix --gfa`, `manual --gfa`, `scaffold --gfa`, `eval scaffold --gfa`, `gapfill`, and `eval gapfill` | `--graph-guard`, `--graph-report`, `--graph-overlap-policy`, graph path depth controls | Graph assignment/report TSVs, graph panels, graph gap fields, fill-path nodes | `tests/test_graph.py`, `tests/test_reference_order.py`, `tests/test_fix_contigs.py`, `tests/test_scaffold.py`, `tests/test_gapfill.py`, `tests/test_manual.py` |
| Review-event TSV model | `eval fix/scaffold/gapfill`, `manual fix/scaffold/gapfill`, `fix --reviewed-plan`, `scaffold --reviewed-plan`, `gapfill --reviewed-plan` | `--review-table` in manual modes; `--reviewed-plan` in executors | `*.fix_review.tsv`, `*.scaffold_review.tsv`, `*.gapfill_review.tsv` using schema `chromosort-review-event-v1` | `tests/test_review.py`, `tests/test_eval.py`, executor tests |
| Long-read PAF evidence summaries | `eval fix --read-paf`, `eval scaffold --read-paf`, `eval gapfill --read-paf`, `manual --read-paf` panels | `--min-read-mapq`, `--read-window-bp`, `--read-terminal-window-bp`, `--read-min-anchor-bp` | `longread_*` review fields for breakpoint, bridge, orientation, read order, and median read-gap support | `tests/test_longreads.py`, `tests/test_eval.py`, `tests/test_manual.py` |
| GAF parsing and graph traversal summaries | `eval fix --gaf`, `eval scaffold --gaf`, `eval gapfill --gaf`, `manual --gaf` panels, and `gapfill --gaf` | `--min-gaf-mapq`, `--min-gaf-path-support`, candidate-path limits | `gaf_*` review fields, gapfill GAF support/status/read fields, selected-event GAF panel rows | `tests/test_gaf.py`, `tests/test_eval.py`, `tests/test_gapfill.py`, `tests/test_manual.py` |
| Scaffold gap and overlap construction | `scaffold` and `eval scaffold` | `--fixed-gap-bp`, `--overlap-policy`, `--trim-sequence-min-identity`, `--reviewed-plan` | Scaffold FASTA, `scaffold_gaps.tsv`, `scaffold_review.tsv` | `tests/test_scaffold.py`, `tests/test_eval.py` |
| Scaffold graph-junction reporting | `scaffold --gfa` and `eval scaffold --gfa` | `--graph-overlap-policy`, `--graph-max-path-edges`; `eval scaffold --gaf` adds path support context | `graph_gaps.tsv`, `graph_status`, `graph_path_nodes`, `gaf_*` review columns | `tests/test_scaffold.py`, `tests/test_eval.py` |
| Gapfill path enumeration and branch resolution | `gapfill` and `eval gapfill` | Required `--gfa`; `--gaf`, `--hic-pairs`, `--ref-paf`, `--max-path-edges`, `--max-candidate-paths`, support thresholds | `gapfill_plan.tsv`, `gapfill_review.tsv`, optional review HTML, optional gapfilled FASTA | `tests/test_gapfill.py`, `tests/test_eval.py` |
| Manual dashboard recipe model | `manual`, `manual fix/scaffold/gapfill`, and `manual apply` | Dashboard browser actions; `manual apply --scaffold/--no-scaffold`, `--gap-bp` | `chromosort-manual-v1` JSON recipe, manual FASTA, manual report | `tests/test_manual.py` |
| Dot-plot rendering | `plot` and the canvas view in `manual` | `--formats`, `--per-ref`, `--sel-ref`, `--assignments`, segment filters | SVG/PDF/PNG plots and manual canvas lines | `tests/test_plot.py`, `tests/test_manual.py` |

### Evidence Activation And Authority

| Evidence stream | Where it enters | Where it is report-only | Where it can affect sequence output |
| --- | --- | --- | --- |
| Whole-genome coords or PAF | `sort`, `clean`, `fix`, `eval fix`, `manual`, `plot` | Manual and plot views; eval tables before execution | `sort`, `clean`, and planner-driven `fix`; reviewed `fix` rows can replace planner inputs |
| GFA graph topology | `sort --gfa`, `fix --gfa`, `manual --gfa`, `scaffold --gfa`, `eval` modes, `gapfill` | Sort/fix/manual graph context; scaffold graph reports by default; eval advisory fields | `scaffold --graph-overlap-policy confirm` can confirm terminal overlap trimming; `gapfill --apply` can insert validated graph-path sequence |
| Long-read PAF to assembly | `eval fix/scaffold/gapfill --read-paf`, `manual --read-paf` | Review-event fields and dashboard panels | Does not directly edit sequence; users may accept or edit reviewed rows informed by the fields |
| Long-read GAF to graph | `eval fix/scaffold/gapfill --gaf`, `manual --gaf`, `gapfill --gaf` | Advisory in `eval fix`; support context in `eval scaffold` and manual panels | In `gapfill`, unique non-conflicting support can choose among candidate GFA paths before sequence validation |
| Hi-C-like graph contacts | `gapfill --hic-pairs` and `eval gapfill --hic-pairs` | Candidate path support and risk annotations | In `gapfill`, unique non-conflicting support can choose a candidate path before sequence validation |
| Reference-placement PAF for graph nodes | `gapfill --ref-paf` and `eval gapfill --ref-paf` | Candidate path support fields and review HTML | In `gapfill`, unique non-conflicting support can choose a candidate path before sequence validation |

## Formal Data Model

### Coordinate Conventions

ChromoSort accepts two alignment coordinate systems:

- MUMmer `show-coords` rows, which are read as 1-based closed endpoints.
- PAF rows, whose 0-based half-open endpoints are converted into equivalent
  1-based closed endpoints by adding one to the start coordinate and preserving
  the end coordinate.

Most interval arithmetic then converts endpoints to half-open intervals
$I = [a,b)$ with length

$$
|I| = b - a.
$$

For reverse-strand rows, interval endpoints are normalized with `min` and
`max` for arithmetic, while the original strand is retained as an orientation
field. This separation keeps overlap and coverage computations independent of
strand direction.

### Alignment Segment and Match Metrics

The shared `Segment` model records:

- reference name $r$ and query/contig name $q$,
- reference and query endpoints,
- aligned reference and query lengths,
- percent identity $p_s$,
- sequence lengths $L_r$ and $L_q$ when provided by the alignment format,
- orientation, MAPQ, and secondary-alignment status when available.

For a query $q$ and reference $r$, let $S_{q,r}$ be the set of passing
alignment segments after identity, MAPQ, secondary-PAF, and length filters. For
segment $s$, let $Q_s$ be its query interval, $R_s$ its reference interval, and
$l_s = |Q_s|$ its query-aligned length. ChromoSort aggregates those segments as:

$$
M^Q_{q,r} = \left|\bigcup_{s \in S_{q,r}} Q_s\right|
$$

$$
M^R_{q,r} = \left|\bigcup_{s \in S_{q,r}} R_s\right|
$$

$$
C^Q_{q,r} = \frac{M^Q_{q,r}}{L_q}, \qquad
C^R_{q,r} = \frac{M^R_{q,r}}{L_r}
$$

$$
\bar{p}_{q,r} =
\frac{\sum_{s \in S_{q,r}} p_s l_s}
     {\sum_{s \in S_{q,r}} l_s}.
$$

The resulting `MatchMetrics` object stores raw aligned bp, merged query bp,
merged reference bp, query coverage, reference coverage, weighted mean identity,
reference start/end, reference midpoint, dominant orientation, segment count,
and merged reference intervals.

The best-reference share used by `sort` is:

$$
B_q =
\frac{M^Q_{q,r^\*}}
     {\sum_{r'} M^Q_{q,r'}}
$$

where $r^\*$ is the best-ranked reference for query $q$.

### Graph Model

`graph.py` reads the subset of GFA needed by ChromoSort:

- `S` records become `GraphNode` objects with name, length, optional sequence,
  tags, and source line number.
- `L` records become `GraphEdge` objects with source node, target node,
  orientations, raw overlap CIGAR, parsed overlap bp when exact, tags, and line
  number.
- `P` and `W` records become `GraphPath` objects with ordered
  `GraphPathStep` entries, path/walk tags, and source line number.
- `AssemblyGraph` indexes outgoing and incoming edges by oriented key
  $(node, orientation)$ where $orientation \in \{+,-\}$.
- `GraphProjection` records store the conservative unitig-to-contig coordinate
  projection used by `chromo graph-map`, `chromo plot --gfa-overlay`, and
  projected graph context in `chromo eval fix`.

The parser is strict for missing linked segments and inconsistent `LN:i`
lengths when sequence is present. It only converts simple overlap CIGARs made
of `M`, `=`, and `X` operations into exact trim lengths. `*` or complex CIGARs
containing insertion, deletion, clipping, skipped, or padding operations return
unknown overlap (`None`) so downstream sequence-changing code cannot treat
them as validated trim lengths.

### Command-Level Records

The command modules layer additional records over these primitives:

- `Assignment` records store the selected reference, status, best and second
  match, overlap-filter annotations, split-candidate annotations, output name,
  and orientation outcome for `sort`.
- `QueryBlock`, `BlockGroup`, `SplitPiece`, and `ContigPlan` represent
  alignment blocks, smoothed breakpoint groups, emitted pieces, and planner
  status for `fix`.
- `ScaffoldMember`, `GapRecord`, `GraphGapRecord`, and `ScaffoldRecord`
  represent ordered contigs, inferred gaps or overlaps, graph evidence at
  adjacent contig junctions, and emitted scaffold records for `scaffold`.
- `FillPlan` records graph path status, support counts, risk flags, fill
  sequence, right-flank trim length, reviewed acceptance, and applied status
  for `gapfill`.
- `ReviewEvent` records use schema `chromosort-review-event-v1` and carry
  task type, candidate action, target, default acceptance, status, confidence,
  reason, notes, and task-specific evidence fields. They are written by
  `chromo eval fix/scaffold/gapfill`, embedded by `chromo manual
  fix/scaffold/gapfill`, and consumed by `fix`, `scaffold`, or `gapfill`
  through `--reviewed-plan`.
- Manual recipes use schema `chromosort-manual-v1` and store 1-based inclusive
  source slices, output names, strand choices, scaffold labels, and removal
  flags.

### Long-Read Evidence Records

`longreads.py` parses long-read-to-assembly PAF records into read-level
alignments keyed by read name and target contig. These records are not a
general-purpose aligner replacement. They are used only by `chromo eval` modes
to summarize breakpoint and contig-end support, and by `chromo manual` to label
the optional long-read PAF evidence panel when review-table fields are present.

`gaf.py` parses long-read-to-graph GAF rows into oriented node paths. The shared
summary functions are used in four places:

- `chromo eval fix --gaf` summarizes reads touching a candidate graph node.
- `chromo eval scaffold --gaf` compares support for selected and alternate GFA
  paths between adjacent contigs.
- `chromo eval gapfill --gaf` and `chromo gapfill --gaf` summarize candidate
  fill-path support from the same traversal logic.
- `chromo manual --gaf` displays `gaf_*` review-table fields for the selected
  event; it records the evidence input but does not independently recompute GAF
  support in the browser generator.

## Core Algorithms

### Alignment Normalization and Metric Construction

The parser layer streams coords or PAF rows, skips malformed rows, applies
filter thresholds, and yields normalized `Segment` objects. The aggregation
layer groups segments by `(query, reference)`, merges query and reference
intervals independently, and computes the formal metrics above. This design is
appropriate for fragmented whole-genome alignments because it avoids double
counting overlapping alignment rows while preserving raw aligned bp and segment
count for review.

### Reference Assignment and Duplicate Filtering

`chromo sort` evaluates every assembly FASTA record against its reference
matches:

1. For each query, matches are ranked by merged query bp, raw query bp, and
   mean identity.
2. The best match is rejected as `below_min_aligned_bp`,
   `below_min_query_cov`, or `ambiguous_ref_match` if it fails the configured
   support thresholds.
3. A strong multi-reference record can be protected as `kept_split_candidate`
   when at least two references meet the split-candidate bp and query-fraction
   thresholds and the best-reference share is not too high. This keeps likely
   chimeric contigs available for `chromo fix`.
4. A large alignment can be rescued as `kept_large_alignment` when it narrowly
   fails query coverage but satisfies the large-alignment thresholds.
5. Kept records are passed to the duplicate-overlap filter unless
   `--no-overlap-filter` is used.

For each reference, the duplicate-overlap pass ranks candidate assignments by
merged reference bp, merged query bp, query coverage, mean identity, and query
length. Higher-ranked contigs claim reference intervals first. Candidate
intervals are either:

- a broad first-to-last reference span (`--overlap-mode span`), or
- exact merged alignment intervals (`--overlap-mode alignment`).

Let $I$ be the candidate interval set and $C$ be the already covered interval
set on the same reference. Novel reference support is:

$$
N(I,C) = |I \setminus C|, \qquad
F(I,C) = \frac{N(I,C)}{|I|}.
$$

A lower-ranked contig is kept when the configured absolute and fractional novel
reference thresholds pass. If not, split-candidate protection and terminal
overlap rescue are considered before the contig is labeled
`terminal_overlap` or `duplicate_overlap`.

Terminal-overlap rescue is intentionally narrower than general overlap
acceptance. It applies only when the novel interval is a one-sided extension
of a single candidate interval:

$$
E = N(I,C), \qquad
E \ge E_{\min}, \qquad
\frac{E}{|I|} \ge f_{\min}.
$$

This distinguishes plausible dovetail extensions from contained duplicate
fragments.

### Breakpoint Planning for `chromo fix`

`chromo fix` starts from filtered alignment segments and converts them into
query-ordered `QueryBlock` records. Adjacent blocks with the same target
reference and orientation can be merged when they are separated by at most the
configured query gap. Candidate split signals are then planned differently by
mode:

| Mode | Candidate signal | Planner signature | Smoothing path | Consequence |
| --- | --- | --- | --- | --- |
| `chromosome` | Reference/chromosome changes only. | `reference` | Dynamic-programming breakpoint smoothing. | Same-reference orientation changes are absorbed into the dominant reference piece. |
| `conservative` | Reference/chromosome changes, plus same-reference orientation events only when the orientation-aware pre-plan meets the complex-inversion criteria. | Usually `reference`; switches to `(reference, orientation)` only for qualifying same-reference complex orientation events. | Dynamic-programming breakpoint smoothing. | Default breakpoint-averse behavior: strong inter-reference chimeras are eligible, simple same-reference inversions are usually left intact. |
| `comprehensive` | All reference/chromosome changes and all same-reference orientation changes. | `(reference, orientation)` | Dynamic-programming breakpoint smoothing. | Broader orientation-aware review mode. It can expose inversion candidates, but it is not a guaranteed superset of `conservative` because the target signature changes the smoothing model and piece-support checks. |
| `sensitive` | Every passing reference/orientation transition after adjacent same-target collapse. | `(reference, orientation)` | No breakpoint-penalty smoothing; only direct run collapse, minimum emitted-piece length, and the breakpoint guard apply. | Debugging or intentionally aggressive scans where every passing transition should be visible. |

The important implementation distinction is that `comprehensive` is not
implemented as "`conservative` plus extra calls." It runs the smoothed planner
with orientation included in every target signature. That can split the evidence
into different candidate pieces or reject a plan for weak-piece support even
when `conservative` would split the same contig using reference-only signatures.

The smoothed planner uses identity-weighted bp:

$$
w_i = l_i \frac{p_i}{100}
$$

for block $i$ with aligned length $l_i$ and identity $p_i$. For an interval of
blocks $G$, support is accumulated by signature $\sigma$ (reference alone or
reference plus orientation, depending on mode). The discordance cost of keeping
that interval as one piece is:

$$
D(G) = \sum_{i \in G} w_i - \max_{\sigma} \sum_{i \in G, \sigma_i=\sigma} w_i.
$$

Dynamic programming chooses a partition of the block list that minimizes:

$$
\operatorname{cost}(\mathcal{P}) =
\sum_{G \in \mathcal{P}} D(G) +
\lambda \left(|\mathcal{P}| - 1\right),
$$

where $\lambda$ is `--breakpoint-penalty-bp`. The planner reports a score
equivalent to the reduction in discordance after paying breakpoint penalties:

$$
\operatorname{score} =
D(\text{all blocks}) -
\sum_{G \in \mathcal{P}} D(G) -
\lambda \left(|\mathcal{P}| - 1\right).
$$

After partitioning, terminal groups that fail support thresholds can be merged
into neighbors, all remaining groups must pass minimum aligned-bp and
query-fraction thresholds, and a per-contig breakpoint guard can reject overly
fragmented plans. Accepted breakpoints are placed midway between adjacent
alignment blocks, then clamped to the contig sequence length. This produces
piece coordinates that are deterministic and easy to audit, while smoothing
small discordant blocks and local structural-variant-like noise by default.

### Explicit Cutting

`chromo cut` is deliberately not evidence-driven. It accepts repeatable
`CONTIG:POS[,POS...]` cut specifications, a cut file, or a single-contig
shorthand. Cut positions are 1-based and mean "cut after this base." The
planner rejects missing contigs, duplicate cut positions, terminal positions,
and any plan that would create pieces shorter than `--min-piece-bp`. It also
checks that emitted FASTA IDs will remain unique.

### Manual Review and Application

`chromo manual` computes the same best-match metrics used by `sort`, embeds
filtered alignment segments, optionally embeds graph node context, and writes a
browser-based dashboard. Initial pieces cover each source contig completely,
with the scaffold label initialized to the best reference or `unplaced`.
Task-specific `manual fix`, `manual scaffold`, and `manual gapfill` modes load
the shared review-event table with `--review-table`, open a focused event
queue, and render selected-event evidence panels for alignment, GFA,
long-read PAF, and GAF inputs when those data are available.

`chromo manual apply` is schema-driven. It validates the recipe schema, verifies
that requested source contigs and slices exist, reverse-complements pieces whose
strand is `-`, omits removed pieces, and optionally groups active pieces into
scaffolds joined by fixed `N` gaps. The recipe is therefore the sole authority
for manual edits.

### Eval Tables and Reviewed Execution

`chromo eval` is the table-first counterpart to task-specific manual review.
It does not introduce a separate decision algorithm; it reuses the algorithm
that the matching executor already exposes, then serializes candidate decisions
as `ReviewEvent` rows with enough evidence fields for spreadsheet or dashboard
review.

- `eval fix` activates the `fix` breakpoint planner and adds optional
  `graph_*`, `gaf_*`, and `longread_*` fields. Accepted rows are executable only
  through `fix --reviewed-plan`, which bypasses alignment re-planning and
  applies the reviewed source slices.
- `eval scaffold` activates scaffold gap/overlap calculation and optional graph
  junction evidence. Accepted rows are executable only through
  `scaffold --reviewed-plan`, where they override matching gap lengths but do
  not reorder, orient, or trim contigs by themselves.
- `eval gapfill` activates the gapfill planner and writes the same path,
  support, validation, and fillability concepts as a review-event table.
  Accepted rows are executable only through `gapfill --reviewed-plan --apply`,
  where the current plan must still match the reviewed path.

### Scaffold Construction

`chromo scaffold` groups kept assignment rows by assigned reference and uses
the ordered FASTA from the same `sort` run. For adjacent contigs $a$ and $b$ on
the same scaffold, the reference-inferred gap is:

$$
g(a,b) = start_b - end_a - 1.
$$

If $g \ge 0$, the scaffold inserts either $g$ Ns or a fixed configured gap. If
$g < 0$, the adjacent reference spans overlap by:

$$
o(a,b) = -g(a,b).
$$

Overlap handling is policy-controlled:

- `zero-gap` writes no Ns and keeps both sequences unchanged.
- `warn` does the same but emits warnings.
- `trim-reference` trims the right-contig prefix by the reference-inferred
  terminal overlap.
- `trim-sequence` trims only when the right-prefix/left-suffix sequence
  identity over the overlap meets the configured threshold.
- `graph-overlap-policy confirm` can trim a terminal overlap when an oriented
  direct GFA edge confirms the junction.

Nonterminal overlaps are reported but not trimmed by the terminal-only trimming
policies.

### Graph-Supported Gap Filling

`chromo gapfill` plans fills between adjacent sorted contigs using oriented GFA
paths. Path enumeration is breadth-first, orientation-aware, bounded by
`--max-path-edges`, bounded by `--max-candidate-paths`, and cycle-aware. Missing
graph nodes produce `missing_node`, absent paths produce `no_graph_path`, and
multiple candidate paths remain `ambiguous_paths` unless exactly one path is
selected by support evidence.

Support evidence is computed as follows:

- GAF support counts read paths that contain the candidate oriented path or its
  reverse-oriented counterpart as a contiguous subpath.
- Hi-C support sums contact counts across adjacent node pairs on the candidate
  path.
- Reference-placement PAF support sums the best identity-weighted query bp for
  intermediate graph nodes whose placements match the scaffold reference,
  orientation, and expected reference gap window.

If different evidence sources uniquely select different paths, ChromoSort
reports the conflict and refuses to fill automatically. When a path is selected,
sequence reconstruction still must pass validation:

- the GFA path must include sequence for both flank nodes and all intermediate
  nodes,
- the oriented flank node sequences must exactly match the adjacent FASTA
  records,
- every applied GFA link overlap must have a known simple bp length,
- overlap lengths must be valid relative to sequence lengths,
- the resulting inserted sequence must not exceed `--max-fill-bp`.

Only `fillable` plans can be applied. With a reviewed plan, the accepted row
must still exist in the current plan, the current path must match the reviewed
path exactly, and the current status must still be `fillable`.

### Plotting

`chromo plot` draws direct visualizations from existing alignment rows. It
constructs cumulative offsets for reference and query FASTA records:

$$
offset(x_i) = \sum_{j<i} L_j
$$

and maps normalized alignment endpoints into a whole-genome or per-reference
canvas. Forward and reverse alignments are rendered as line segments with
different colors. Optional assignment input reorders the query axis by kept
ChromoSort assignments without re-running placement. Optional selected-reference
input filters the reference axis, alignment rows, and `--per-ref` panel set
without changing the underlying FASTA or alignment file. The
[dot-plot guide]({{ '/dot-plots/' | relative_url }}) describes how these
segments map to common review patterns such as clean collinearity, reversed
contigs, internal inversions, chimeric candidates, duplicate signal, and
alignment gaps.

When `--gfa-overlay` is supplied, unitig-level graph intervals are first
projected through GFA `P` path or `W` walk records whose names match query FASTA
contigs. The overlay then draws query-axis graph boundaries in the same
coordinate system as the alignment rows; if projection is impossible,
ChromoSort reports that rather than comparing unitig and contig coordinates
directly.

## Evidence, Scoring, Ranking, or Decision Logic

ChromoSort does not currently implement a probabilistic confidence model.
Instead, it uses deterministic scores and thresholded evidence quantities that
are visible in reports.

| Decision point | Quantity | Interpretation |
| --- | --- | --- |
| Best reference in `sort` | `merged_query_bp`, `raw_query_bp`, `avg_identity` | Ranks each contig's candidate reference matches. |
| Ambiguous reference match | $B_q$ | Fraction of all merged query-aligned bp assigned to the best reference. |
| Duplicate-overlap filtering | $N(I,C)$ and $F(I,C)$ | Absolute and fractional novel reference contribution relative to stronger kept contigs. |
| Split-candidate protection | Number of supported references and best-reference share | Prevents strong multi-reference contigs from being discarded before `fix` review. |
| Breakpoint smoothing | $\operatorname{cost}(\mathcal{P})$ and planner score | Balances discordant alignment support against a fixed breakpoint penalty. |
| Long-read breakpoint support | Spanning, split, edge, and nearby read counts in `eval fix` | Reports support around candidate breakpoints without automatically changing the fix planner. |
| Scaffold overlap trimming | $g(a,b)$, $o(a,b)$, overlap class, optional sequence identity | Determines gap length and whether an explicitly allowed trim is safe to perform. |
| Long-read bridge support | Bridge count, orientation summary, read-order summary, and median read-space gap in `eval scaffold/gapfill` | Reports whether reads support adjacent contig ends in assembly-coordinate space. |
| GAF node or path support | Node traversal counts, selected-path support, best alternate support, and support status | Reports graph-read support in eval/manual; can resolve gapfill branches only under gapfill support thresholds. |
| Reviewed row acceptance | `accept` field in `chromosort-review-event-v1` rows | Allows human-reviewed rows to drive `fix`, `scaffold`, or `gapfill` after task-specific validation. |
| Gapfill branch selection | Unique GAF, Hi-C, or reference-placement support | Resolves graph branches only when one path clears threshold and is not contradicted. |
| Gapfill risk annotation | Branching, high-degree nodes, self-loops, unsequenced nodes, cycles, weak/tied support | Flags cases that need review even when a path is reported. |

The gapfill branch-complexity score is additive and intentionally simple:

$$
R =
\max(0, P-1) +
|H| + 2|S| + 2|U| + C + |X|,
$$

where $P$ is the number of candidate paths, $H$ is the set of high-degree nodes,
$S$ is the set of self-loop nodes, $U$ is the set of unsequenced nodes, $C$ is
the number of avoided cycles during path search, and $X$ is the set of extra
warning flags such as weak, tied, conflicting, or sequence-validation failures.
The score is a review aid, not a calibrated probability.

## Region, Window, Graph, or Table Construction Logic

### Interval Regions

All alignment-derived region logic is built from merged half-open intervals.
This applies to coverage, duplicate-overlap filtering, terminal-overlap
classification, and reference span reporting. A candidate's overlap with prior
claims is computed by two-pointer interval intersection after sorting; novel
intervals are computed by subtracting the covered set from the candidate set.

### Reference Gap Windows

Scaffold gaps and reference-placement gapfill support use adjacent assignment
coordinates. For same-reference gapfill candidates, the expected placement
window is bounded by the adjacent assignment edges:

$$
W(a,b) =
[\min(end_a, start_b), \max(end_a, start_b)].
$$

Intermediate graph-node placements from `--ref-paf` support a candidate path
only when their reference midpoint falls inside this window and the orientation
matches the candidate path.

### Graph Paths

Graph path construction uses oriented keys. Starting keys are determined from
the left assignment orientation when known, otherwise both orientations are
allowed. Target keys are built the same way from the right assignment. The BFS
state contains the current oriented node, the path of edges taken so far, and
the seen oriented nodes. Revisiting a seen oriented node is counted as an
avoided cycle and is not expanded.

### Output Tables

TSV reports are not secondary decorations; they are command contracts. For
example:

- `sort` assignment tables feed `scaffold`, `gapfill`, and ordered plots.
- `eval` review-event tables feed task-specific `manual` dashboards and the
  matching `fix`, `scaffold`, or `gapfill` reviewed execution path.
- `fix`, `cut`, and `manual apply` reports describe every emitted piece.
- `scaffold_gaps.tsv` records raw inferred gaps, final gap sizes, overlap
  class, policy, action, trim length, and optional sequence identity.
- `graph_gaps.tsv` records direct GFA adjacency, shortest reported path, and
  missing-node status for scaffold junctions.
- `gapfill_plan.tsv` records selected paths, alternative support, risk flags,
  fill status, fill length, right trim, acceptance, and applied status.

## Cross-source or Multi-stage Reconciliation

ChromoSort reconciles evidence sources through explicit names and stage-local
tables:

- Alignment evidence is keyed by reference and query FASTA names. PAF and
  coords inputs are interchangeable only after normalization into `Segment`
  records.
- Graph evidence resolves contig names conservatively. Scaffold and gapfill
  junctions try original contig names, ChromoSort `new_name` values, and FASTA
  record names when mapping records back to GFA segments.
- Review-event rows are task-scoped. `manual fix/scaffold/gapfill` rejects
  rows whose task does not match the selected mode, and executors consume only
  accepted rows for their own task.
- Manual-review recipes are validated against the original assembly FASTA; they
  do not inherit implicit state from the dashboard.
- Reviewed scaffold gaps are checked against the current scaffold, left contig,
  and right contig names before application.
- Reviewed gapfill plans are checked against the current graph path before
  application, preventing stale accepted paths from being silently applied.
- Sequence-editing stages can change contig IDs and coordinates. The intended
  workflow is to re-align edited FASTA records before using downstream
  alignment-dependent commands.

This reconciliation strategy favors reproducibility over convenience. It avoids
implicit cross-stage mutation and forces each downstream stage to consume the
current FASTA and the corresponding current reports.

## Implementation Architecture

The source tree is a small Python package under `src/chromosort/`.

| Module | Role |
| --- | --- |
| `cli.py` | Top-level `chromo` dispatcher and subcommand routing. |
| `reference_order.py` | FASTA length scanning, alignment parsing, interval utilities, match aggregation, assignment decisions, duplicate-overlap filtering, ordered FASTA writing, and sort reports. |
| `fix_contigs.py` | Chimeric-contig block collection, breakpoint planning, smoothing, split-piece writing, graph guard/report integration, and split reports. |
| `evaluate.py` | Table-only review-event generation for `eval fix`, `eval scaffold`, and `eval gapfill`, reusing executor planners and adding optional graph, long-read PAF, and GAF evidence fields. |
| `review.py` | Shared review-event schema, TSV column ordering, accept-value parsing, duplicate-id checks, and read/write helpers used by eval, manual, and reviewed executors. |
| `longreads.py` | Long-read-to-assembly PAF parsing and summaries for breakpoint support and contig-end bridge evidence. |
| `gaf.py` | Long-read-to-graph GAF parsing, oriented path matching, candidate traversal support summaries, and node-level advisory summaries. |
| `cut.py` | Explicit cut parsing, validation, FASTA splitting, and cut reports. |
| `manual.py` | Dashboard data construction, graph context embedding, task review-event queues, modular evidence panels, recipe validation, manual piece/scaffold application, and manual reports. |
| `plot.py` | Alignment filtering, axis construction, SVG/PDF/PNG rendering, and per-reference plots. |
| `scaffold.py` | Assignment-table reading, scaffold grouping, reviewed gap overrides, reference gap/overlap logic, graph junction reporting, and scaffold FASTA writing. |
| `gapfill.py` | GAF/Hi-C/reference-placement support integration, graph path enumeration, branch resolution, sequence reconstruction, review dashboard generation, reviewed-plan validation, and optional gapfilled FASTA writing. |
| `graph.py` | GFA `S`/`L` parser, oriented graph indexes, direct-edge queries, node evidence, and link evidence. |
| `paths.py` | Output directory creation and small shared filesystem helpers. |

Most dependencies are from the Python standard library. Pillow is used for PNG
plot output. Jekyll-related files support the documentation site and are
separate from runtime package behavior.

## Complexity and Scaling

The following bounds describe the implemented algorithms, not biological
limits. Let $S$ be the number of passing alignment segments, $V$ graph nodes,
$E$ graph edges, $Q$ query contigs, $B_q$ fix blocks for contig $q$, and $G$
adjacent scaffold gaps.

| Component | Time complexity | Memory behavior |
| --- | --- | --- |
| Alignment parsing | $O(S)$ | Stores grouped intervals and segment summaries needed by the command. |
| Match aggregation | $O(S + \sum_{q,r} k_{q,r}\log k_{q,r})$ | Interval lists per query-reference pair before merging. |
| Assignment ranking | $O(\sum_q R_q \log R_q)$ where $R_q$ is matched references per query | Per-query match lists. |
| Duplicate-overlap filtering | $O(\sum_r A_r \log A_r + I_r)$ plus interval subtraction/intersection work | Covered interval sets per reference; worst cases grow with fragmentation. |
| Fix block merging | $O(B_q \log B_q)$ per contig | Query blocks for selected contigs. |
| Fix smoothing DP | $O(B_q^2)$ per selected contig | DP table and interval-cost cache over block ranges. |
| Cut planning | $O(C \log C)$ per cut contig plus FASTA streaming | Stores cut positions and FASTA length records. |
| GFA parsing | $O(V+E)$ | Stores all parsed nodes, edges, and orientation indexes. |
| Review-event TSV read/write | $O(R)$ for $R$ rows | Stores reviewed rows while checking task and event-id uniqueness. |
| Long-read PAF support summaries | Proportional to read alignments for the queried contig or contig pair | Stores filtered read alignments grouped by read and target. |
| GAF traversal summaries | Proportional to candidate paths times filtered GAF read-path length | Stores filtered GAF records and per-path read-name summaries. |
| Scaffold construction | $O(Q+G)$ without graph path reporting | Stores ordered FASTA records and scaffold gap reports. |
| Scaffold graph path report | $O(G \cdot d^D)$ in the bounded worst case | Bounded by `--graph-max-path-edges`; reports one shortest path per gap. |
| Gapfill path enumeration | $O(G \cdot \min(P,d^D) \cdot D)$ | Bounded by `--max-candidate-paths` $P$ and `--max-path-edges` $D$. |
| Gapfill GAF support | Proportional to candidate paths times read-path scan length | Stores filtered GAF path records. |
| Plot rendering | $O(S)$ draw elements | PNG additionally requires $O(width \cdot height)$ image memory. |

The most performance-sensitive code paths are interval merging, the fix
planner's quadratic DP over alignment blocks, graph BFS at scaffold/gapfill
junctions, and large plot rendering. The exposed depth, candidate, segment,
and format filters exist to bound those costs for larger assemblies.

## Validation Strategy

The repository contains a pytest/unittest suite configured in `pyproject.toml`.
The tests are synthetic and deterministic; they are best interpreted as unit
and integration evidence for the implemented decision rules, not as biological
benchmarking.

Validation coverage currently includes:

- CLI dispatch and mutually exclusive argument checks.
- coords and PAF parity for sorting and fixing.
- best-reference assignment, duplicate-overlap rejection, split-candidate
  protection, terminal-overlap classification, terminal-extension rescue, and
  large-alignment rescue.
- graph assignment reports and graph guard warnings for sort/fix cases.
- breakpoint smoothing, sensitive/comprehensive/conservative fix modes,
  inversion handling, noisy small-event rejection, contig-subset parity, and
  breakpoint guardrails.
- review-event schema reading/writing, accepted/rejected value parsing,
  duplicate-id rejection, and task mismatch rejection.
- `eval fix`, `eval scaffold`, and `eval gapfill` table generation, including
  reviewed executor handoff cases.
- long-read PAF breakpoint and contig-bridge summaries.
- GAF parsing, MAPQ filtering, node summaries, selected-path support, alternate
  support, and mixed GFA/PAF/GAF scaffold evidence.
- explicit cut validation for duplicate, terminal, and file-based cut inputs.
- manual dashboard ordering, optional sequence embedding, graph context,
  task-event queues, modular evidence panels, recipe slicing,
  reverse-complementing, removal, and scaffold grouping.
- plot generation for PDF, SVG, PNG, per-reference plots, and assignment-based
  query ordering.
- scaffold inferred/fixed gaps, negative reference-gap overlap reporting,
  reference trimming, sequence-confirmed trimming, and graph-confirmed trimming.
- GFA parsing of orientations, simple versus complex overlap CIGARs, missing
  links, and length mismatches.
- gapfill unique path application, reviewed-plan acceptance/rejection, stale
  path rejection, ambiguous branch refusal, GAF/Hi-C/reference-placement branch
  resolution, conflict refusal, risk flags, and review HTML embedding.

The `tests/data/graph_gotchas/` fixture documents intentionally difficult graph
cases: an unambiguous path, an ambiguous branch, a cycle, orientation-specific
links, disconnected mapped nodes, unsequenced segments, and repeat/duplicate
warning contexts.

Runtime validation is also part of the architecture. Many user-facing errors
are raised as `ValueError` and converted to command exit status 2 with a clear
message. Examples include mutually exclusive alignment inputs, malformed GFA
records, missing assignment-table columns, stale reviewed gapfill paths, invalid
cut positions, missing source FASTA records, and unsupported graph safety flags
without `--gfa`.

The current GitHub Actions workflow deploys the documentation site from
`docs/`. At the time this document was written, test execution is configured
locally but is not represented as a repository CI workflow.

## Non-goals and Limitations

- ChromoSort does not run MUMmer, minimap2, GraphAligner, Hi-C processing, or
  assembly tools internally. It consumes their outputs.
- ChromoSort does not polish base calls, call variants, estimate haplotypes, or
  validate biological correctness beyond the supplied evidence.
- The ranking and risk scores are deterministic heuristics, not calibrated
  probabilities or statistical confidence intervals.
- Reference-guided decisions inherit reference bias. Contigs absent from, highly
  diverged from, or structurally incompatible with the reference can be
  under-supported or classified for review.
- Duplicate-overlap filtering is per reference and based on reference interval
  novelty. It is not a global assembly graph optimization.
- GFA support covers segment (`S`), link (`L`), path (`P`), and walk (`W`)
  records. Path/walk records are used for coordinate projection when their
  names match contig FASTA records; GFA files without those records can still
  provide topology but cannot project unitig intervals onto contig coordinates.
- Sequence-changing graph fills require segment sequences and simple exact link
  overlap lengths. Missing sequences, complex overlap CIGARs, flank mismatches,
  invalid overlaps, overly long fills, or ambiguous/conflicting branches prevent
  automatic filling.
- Scaffold overlap trimming is policy-limited and terminal-focused. Nonterminal
  overlaps are reported for review rather than broadly edited.
- After FASTA-changing steps, downstream alignment-dependent commands should be
  run with fresh alignments against the edited FASTA.
- The package metadata currently classifies the project as alpha, so reviewers
  should treat command behavior as implemented and tested, but still evolving.

## Reviewer's Reading Guide

For code review, start with the command boundary and then follow the evidence
models into the algorithms:

1. `src/chromosort/cli.py` for public command routing.
2. `src/chromosort/reference_order.py` for alignment normalization, interval
   math, match aggregation, sort assignment, and duplicate-overlap logic.
3. `src/chromosort/fix_contigs.py` for breakpoint planning, smoothing, mode
   behavior, and emitted split pieces.
4. `src/chromosort/evaluate.py` and `src/chromosort/review.py` for
   review-event table generation, schema constraints, and reviewed handoff
   fields.
5. `src/chromosort/graph.py`, `src/chromosort/gaf.py`, and
   `src/chromosort/longreads.py` for graph and read evidence constraints.
6. `src/chromosort/scaffold.py` and `src/chromosort/gapfill.py` for
   junction-level gap, overlap, graph path, and fill behavior.
7. `src/chromosort/manual.py`, `src/chromosort/cut.py`, and
   `src/chromosort/plot.py` for review, explicit edit, and visualization
   surfaces.
8. `tests/test_reference_order.py`, `tests/test_fix_contigs.py`,
   `tests/test_eval.py`, `tests/test_review.py`, `tests/test_longreads.py`,
   `tests/test_gaf.py`, `tests/test_scaffold.py`, `tests/test_gapfill.py`,
   `tests/test_graph.py`, `tests/test_manual.py`, `tests/test_cut.py`,
   `tests/test_plot.py`, and `tests/test_cli.py` for executable examples of
   the expected behavior.

For a methods review, the most important invariants to check are:

- interval normalization and merged-bp accounting,
- conservative handling of ambiguous or duplicate alignments,
- breakpoint smoothing and support thresholds,
- explicitness of every sequence-changing operation,
- review-event task scoping and accepted-row validation,
- graph path validation before sequence insertion,
- distinction between report-only evidence and evidence allowed to select or
  apply a sequence-changing path,
- stale-plan prevention for reviewed gap fills,
- report schemas used as contracts between stages.
