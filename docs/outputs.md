---
title: Output Files
description: Guide to ChromoSort FASTA, TSV, HTML, plot, and run-summary outputs.
---

# Output Files

Most ChromoSort commands write both sequence files and audit tables. FASTA
outputs are intended for downstream assembly workflows; TSV, HTML, plot, and
text reports are intended to make every keep, reject, split, cut, scaffold, or
graph-fill decision inspectable.

## AGP And Provenance Sidecars

Every ChromoSort command that writes a modified FASTA also writes a fresh AGP
2.1 map for that FASTA, a component-provenance TSV mirroring the AGP rows, and
a submission-oriented checklist. Downstream commands should compose a new AGP
for their own output rather than appending rows to an upstream AGP, because each
FASTA-changing step creates a new object coordinate system.

The AGP file is the structural map: output objects, component slices, inserted
N gaps, orientations, and part numbers. The component TSV and command-specific
reports carry the richer audit trail: whether a row was unchanged, algorithmic,
reviewed, manual, graph-filled, or a fallback gap. For submission workflows,
keep the final FASTA, final AGP, component TSV, command reports, and
`*.submission_checklist.tsv` together.

## Output Index

**Table 1. ChromoSort command output map.** Use this as a quick index to the
command-specific output sections below.

| Command | Primary outputs |
| --- | --- |
| [`chromo sort`](#chromo-sort-outputs) | `<prefix>.ordered.fa`, `<prefix>.ordered.agp`, `<prefix>.ordered_components.tsv`, assignment/match/chromosome reports, `<prefix>.submission_checklist.tsv`, optional `<prefix>.graph_assignments.tsv`, and `<prefix>.run_summary.txt`. |
| [`chromo clean`](#chromo-clean-outputs) | `<prefix>.clean.fa`, `<prefix>.clean.agp`, `<prefix>.clean_components.tsv`, initial-sort reports, `<prefix>.fix_targets.txt`, `<prefix>.fix_report.tsv`, `<prefix>.clean_contigs.tsv`, `<prefix>.clean_chromosome_summary.tsv`, optional discarded FASTA, `<prefix>.submission_checklist.tsv`, and `<prefix>.run_summary.txt`. |
| [`chromo eval`](#chromo-eval-outputs) | Task-specific review-event TSVs: `<prefix>.fix_review.tsv`, `<prefix>.scaffold_review.tsv`, or `<prefix>.gapfill_review.tsv`; `chromo eval all` writes all three plus `<prefix>.eval_all_outputs.tsv`. |
| [`chromo fix`](#chromo-fix-outputs) | Fixed FASTA at `--output-fasta`, split report at `--report`, default sidecars `<output-fasta>.agp`, `<output-fasta>.components.tsv`, `<output-fasta>.submission_checklist.tsv`, and optional graph report. |
| [`chromo cut`](#chromo-cut-outputs) | Cut FASTA at `--output-fasta`, cut-piece report at `--report`, and default sidecars `<output-fasta>.agp`, `<output-fasta>.components.tsv`, and `<output-fasta>.submission_checklist.tsv`. |
| [`chromo manual`](#chromo-manual-outputs) | Self-contained HTML dashboard, browser FASTA download, recipe JSON download, and reproducible `manual apply` FASTA/report outputs with AGP/component/checklist sidecars. |
| [`chromo gafprep`](#chromo-gafprep-outputs) | `<prefix>.targets.tsv`, `<prefix>.selected_reads.tsv`, `<prefix>.selected_read_review_links.tsv`, `<prefix>.selected.fastq.gz`, `<prefix>.graphaligner.gfa`, `<prefix>.graphaligner.sh`, sanitization reports, and `<prefix>.summary.tsv`. |
| [`chromo plot`](#chromo-plot-outputs) | Whole-genome and optional per-reference dot plots in PDF, SVG, or PNG, with interpretation examples in the [dot-plot guide]({{ '/dot-plots/' | relative_url }}). |
| [`chromo graph-map`](#chromo-graph-map-outputs) | `<prefix>.utg_to_ctg.tsv`, `<prefix>.path_summary.tsv`, and `<prefix>.warnings.tsv`. |
| [`chromo scaffold`](#chromo-scaffold-outputs) | `<prefix>.scaffold.fa`, `<prefix>.scaffold.agp`, `<prefix>.scaffold_components.tsv`, `<prefix>.scaffold_gaps.tsv`, optional `<prefix>.graph_gaps.tsv`, `<prefix>.scaffold_summary.tsv`, `<prefix>.submission_checklist.tsv`, and `<prefix>.run_summary.txt`. |
| [`chromo gapfill`](#chromo-gapfill-outputs) | `<prefix>.gapfill_plan.tsv`, `<prefix>.gapfilled.agp`, `<prefix>.gapfilled_components.tsv`, optional review HTML, optional `<prefix>.gapfilled.fa`, `<prefix>.submission_checklist.tsv`, and `<prefix>.run_summary.txt`. |

## `chromo sort` Outputs

Use `chromo sort` outputs to inspect contig assignment, duplicate-overlap
filtering, reference coverage, and the final reference-ordered FASTA.

| Output | Description |
| --- | --- |
| `<prefix>.ordered.fa` | Retained contigs ordered by reference chromosome and position. |
| `<prefix>.ordered.agp` | AGP 2.1 map for the ordered FASTA. Each row maps an emitted record to its original contig and orientation. |
| `<prefix>.ordered_components.tsv` | Human-readable provenance table mirroring the ordered AGP rows. |
| `<prefix>.contig_assignments.tsv` | One row per assembly contig with final status, assignment metrics, overlap classification, and split-candidate flags. |
| `<prefix>.contig_ref_matches.tsv` | One row per contig-reference match before final assignment. |
| `<prefix>.chromosome_summary.tsv` | One row per reference sequence with retained contig lists and covered reference bp. |
| `<prefix>.graph_assignments.tsv` | Optional report-only graph evidence for assignment and duplicate-overlap decisions when `--gfa` is provided. |
| `<prefix>.submission_checklist.tsv` | FASTA/AGP consistency checks and file manifest for downstream submission packaging. |
| `<prefix>.run_summary.txt` | Inputs, thresholds, output paths, status counts, and PAF diagnostics when `--paf` is used. |

**Table 2. Example `contig_assignments.tsv` rows.** Selected columns from the
synthetic fixture show common status classes.

| contig | status | kept | new_name | assigned_ref | query_cov | avg_identity | overlap_class |
| --- | --- | --- | --- | --- | ---: | ---: | --- |
| `contigA` | `kept` | `yes` | `chr1_contigA` | `chr1` | `1.0000` | `99.000` | `.` |
| `contigB` | `kept` | `yes` | `chr1_contigB` | `chr1` | `1.0000` | `98.500` | `.` |
| `contigDup` | `duplicate_overlap` | `no` | `.` | `chr1` | `1.0000` | `99.000` | `contained_overlap` |
| `contigNo` | `no_alignment` | `no` | `.` | `.` | `0.0000` | `0.000` | `.` |

**Table 3. Example `chromosome_summary.tsv` rows.** The summary table groups
retained contigs by reference.

| ref | ref_length | kept_contigs | covered_ref_bp | ref_cov | ordered_new_names |
| --- | ---: | ---: | ---: | ---: | --- |
| `chr1` | `120` | `2` | `120` | `1.0000` | `chr1_contigA,chr1_contigB` |
| `chr2` | `100` | `1` | `60` | `0.6000` | `chr2_contigC` |

**Listing 1. Example ordered FASTA headers.** Headers retain assignment and
orientation metadata for audit.

```text
>chr1_contigA original=contigA ref=chr1 ref_start=1 ref_end=80 orientation=+ reverse_complemented=no query_cov=1.0000 avg_identity=99.000
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>chr2_contigC original=contigC ref=chr2 ref_start=10 ref_end=69 orientation=- reverse_complemented=yes query_cov=1.0000 avg_identity=97.000
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
```

## `chromo clean` Outputs

Use `chromo clean` outputs to inspect the full conservative cleanup decision:
which raw contigs were discarded by sort logic, which retained contigs were
selected for fix planning, which split plans were accepted or rejected, and what
records were emitted into the cleaned FASTA.

| Output | Description |
| --- | --- |
| `<prefix>.clean.fa` | Cleaned FASTA with retained unsplit contigs and accepted split pieces, oriented if requested and ordered by reference placement. |
| `<prefix>.clean.agp` | AGP 2.1 map for the cleaned FASTA, mapping each emitted record to its raw source-contig slice and orientation. |
| `<prefix>.clean_components.tsv` | Human-readable provenance table mirroring the clean AGP rows, with clean status for each component. |
| `<prefix>.initial_sort.contig_assignments.tsv` | Initial raw-contig sort assignment report. |
| `<prefix>.initial_sort.contig_ref_matches.tsv` | Initial raw-contig per-reference match report. |
| `<prefix>.initial_sort.chromosome_summary.tsv` | Initial sort chromosome summary. |
| `<prefix>.fix_targets.txt` | Original raw contig IDs selected for fix planning. |
| `<prefix>.fix_report.tsv` | `chromo fix`-style report for selected retained contigs. |
| `<prefix>.clean_contigs.tsv` | Unified row-level audit table for discarded contigs, retained unsplit contigs, and split pieces. |
| `<prefix>.clean_chromosome_summary.tsv` | Final cleaned-record summary grouped by reference sequence. |
| `<prefix>.submission_checklist.tsv` | FASTA/AGP consistency checks and file manifest for downstream submission packaging. |
| `<prefix>.run_summary.txt` | Inputs, outputs, sort/fix/clean status counts, and a reminder to re-align the cleaned FASTA for validation. |
| `--discarded-fasta` path | Optional FASTA containing raw contigs discarded by sort filtering. |

**Example `clean_contigs.tsv` rows.** Selected columns show a discarded raw
contig, two emitted split pieces, and one retained unsplit contig.

| source_contig | clean_status | clean_name | sort_status | fix_selected | fix_status | dominant_ref |
| --- | --- | --- | --- | --- | --- | --- |
| `contig_01` | `discarded_no_alignment` | `.` | `no_alignment` | `no` | `.` | `.` |
| `contig_04` | `kept_split_piece` | `chrom02_contig_04_a` | `kept_split_candidate` | `yes` | `split` | `chrom02` |
| `contig_04` | `kept_split_piece` | `chrom07_contig_04_b` | `kept_split_candidate` | `yes` | `split` | `chrom07` |
| `contig_inv_mid` | `not_split_single_target` | `chrom06_contig_inv_mid` | `kept` | `yes` | `not_split_single_target` | `chrom06` |

## `chromo fix` Outputs

Use `chromo fix` outputs to inspect which requested or automatically detected
contigs were split, which were copied unchanged, and why each breakpoint plan
was accepted or rejected.

| Output | Description |
| --- | --- |
| `--output-fasta` | Full fixed assembly FASTA by default, with split pieces replacing fixed contigs. |
| `--report` | TSV report describing split pieces and unsplit requested contigs. |
| `--agp` / default `<output-fasta>.agp` | AGP 2.1 map for the fixed FASTA. Unchanged records are identity rows; split pieces map back to source-contig slices. |
| `--components` / default `<output-fasta>.components.tsv` | Human-readable provenance table mirroring the fixed AGP rows. |
| `--submission-checklist` / default `<output-fasta>.submission_checklist.tsv` | FASTA/AGP consistency checks and file manifest for downstream submission packaging. |
| `--graph-report` | Optional graph context TSV when `--gfa` is provided. Defaults to the `--report` path with a `.graph.tsv` suffix. |

**Table 4. Example split report rows.** Selected columns from a sensitive-mode
fixture show how one source contig can become multiple reference-labeled pieces.

| original_contig | status | new_contig | part_index | dominant_ref | slice_start | slice_end | piece_bp | orientation |
| --- | --- | --- | ---: | --- | ---: | ---: | ---: | --- |
| `contig_04` | `split` | `chrom02-contig_04-a` | `1` | `chrom02` | `1` | `20` | `20` | `+` |
| `contig_04` | `split` | `chrom07-contig_04-b` | `2` | `chrom07` | `21` | `40` | `20` | `+` |
| `contig_12` | `split` | `chrom04-contig_12-a` | `1` | `chrom04` | `1` | `5` | `5` | `+` |
| `contig_12` | `split` | `chrom05-contig_12-b` | `2` | `chrom05` | `6` | `35` | `30` | `+` |

**Listing 2. Example fixed FASTA records.** Split-piece headers carry the
original contig, slice interval, alignment interval, and orientation.

```text
>chrom02-contig_04-a original=contig_04 ref=chrom02 slice=1-20 alignment=1-20 orientation=+ reverse_complemented=no avg_identity=100.000
AAAAAAAAAAAAAAAAAAAA
>chrom07-contig_04-b original=contig_04 ref=chrom07 slice=21-40 alignment=21-40 orientation=+ reverse_complemented=no avg_identity=100.000
CCCCCCCCCCCCCCCCCCCC
```

## `chromo eval` Outputs

Use `chromo eval` outputs when a task should be reviewed in a spreadsheet or in
task-specific `chromo manual` dashboards before a sequence-changing command
runs.

| Output | Description |
| --- | --- |
| `<prefix>.fix_review.tsv` | Shared review-event table for candidate `split_piece` rows. Optional `graph_*`, `gaf_*`, and `longread_*` columns describe provenance; accepted rows can drive `chromo fix --reviewed-plan`. |
| `<prefix>.scaffold_review.tsv` | Shared review-event table for adjacent scaffold junctions. Accepted `scaffold_gap` rows can override matching gap lengths through `chromo scaffold --reviewed-plan`. |
| `<prefix>.gapfill_review.tsv` | Shared review-event table for graph fill candidates. Accepted `fill_path` rows can restrict `chromo gapfill --reviewed-plan --apply` after the current path is revalidated. |
| `<prefix>.eval_all_outputs.tsv` | Manifest written by `chromo eval all` listing the three review tables and their recommended `--eval-review-table` arguments for `chromo gafprep`. |

These tables use schema `chromosort-review-event-v1`. Rejected or deleted rows
do not apply sequence changes. Task-specific `manual` dashboards load the same
tables with `--review-table` and render their evidence fields as event panels.

## `chromo cut` Outputs

Use `chromo cut` outputs when exact reviewed breakpoints should be applied
without invoking the alignment-based `chromo fix` planner.

| Output | Description |
| --- | --- |
| `--output-fasta` | Full assembly FASTA with requested contigs replaced by cut pieces. Uncut contigs are copied unchanged. |
| `--report` | TSV report describing every emitted cut piece, including original contig, new contig, slice coordinates, piece length, and cut positions. |
| `--agp` / default `<output-fasta>.agp` | AGP 2.1 map for the cut FASTA. Cut pieces map back to original contig slices; uncut records are identity rows. |
| `--components` / default `<output-fasta>.components.tsv` | Human-readable provenance table mirroring the cut AGP rows. |
| `--submission-checklist` / default `<output-fasta>.submission_checklist.tsv` | FASTA/AGP consistency checks and file manifest for downstream submission packaging. |

**Table 5. Example cut report rows.** Cutting `contigA` after bases 20 and 50
emits three adjacent slices and records the complete cut set on each row.

| original_contig | new_contig | part_index | slice_start | slice_end | piece_bp | cut_after_positions |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| `contigA` | `contigA_cut001` | `1` | `1` | `20` | `20` | `20,50` |
| `contigA` | `contigA_cut002` | `2` | `21` | `50` | `30` | `20,50` |
| `contigA` | `contigA_cut003` | `3` | `51` | `80` | `30` | `20,50` |

**Listing 3. Example cut FASTA records.** Requested contigs are replaced by cut
pieces; unrequested contigs are copied through unchanged.

```text
>contigA_cut001
AAAAAAAAAAAAAAAAAAAA
>contigA_cut002
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>contigA_cut003
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>contigB
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
```

## `chromo manual` Outputs

Use `chromo manual` outputs when reviewed browser curation should be exported
and reproduced from the command line.

| Output | Description |
| --- | --- |
| `--output-html` | Self-contained manual-edit dashboard with embedded alignment/order metadata, optional graph context, task review queues, modular evidence panels, and optional embedded sequences. |
| Browser FASTA download | Edited FASTA exported from the dashboard. The source FASTA is not modified. |
| Browser recipe download | JSON recipe describing the current manual edits. |
| `chromo manual apply --output-fasta` | FASTA generated by applying a recipe from the command line. |
| `chromo manual apply --report` | Optional TSV report of emitted pieces and removed pieces. |
| `chromo manual apply --agp` / default `<output-fasta>.agp` | AGP 2.1 map for the manual-apply FASTA. In scaffold mode, active pieces become component rows and inserted N gaps become AGP gap rows. |
| `chromo manual apply --components` / default `<output-fasta>.components.tsv` | Human-readable provenance table mirroring the manual AGP rows. |
| `chromo manual apply --submission-checklist` / default `<output-fasta>.submission_checklist.tsv` | FASTA/AGP consistency checks and file manifest for downstream submission packaging. |

Manual scaffold-mode gaps default to AGP `contig/no/na`, which records an
unlinked N gap without claiming evidence. When a reviewed manual scaffold join
does have support, set `--agp-gap-type`, `--agp-linkage`, and
`--agp-linkage-evidence` explicitly during `chromo manual apply`.

<figure>
  <img src="https://rotheconrad.github.io/chromosort/assets/chromo_manual_graph_review.png" alt="Screenshot-style view of chromo manual showing contig graph filters, a selected-contig dot plot, and a graph-neighborhood panel." width="900">
  <figcaption><strong>Figure 1. Example <code>chromo manual</code> dashboard output.</strong> The self-contained HTML dashboard combines a contig list, dot-plot review panel, graph-neighborhood filters, and read-only GFA context for the selected contig.</figcaption>
</figure>

**Table 6. Example `chromo manual apply` report rows.** The command-line apply
step writes piece-level audit information for a browser-reviewed recipe.

| piece_id | source | output_name | scaffold | slice_start | slice_end | piece_bp | strand | removed | export_mode |
| --- | --- | --- | --- | ---: | ---: | ---: | --- | --- | --- |
| `contigA` | `contigA` | `contigA_manual001` | `chr1` | `1` | `80` | `80` | `+` | `no` | `pieces` |
| `contigB` | `contigB` | `contigB_manual002` | `chr1` | `1` | `40` | `40` | `+` | `no` | `pieces` |
| `contigNo` | `contigNo` | `contigNo_manual003` | `unplaced` | `1` | `30` | `30` | `+` | `no` | `pieces` |

**Listing 4. Example manual-apply FASTA records.** Browser edits become ordinary
FASTA records with the original source contig, slice interval, strand, and
scaffold label in the header.

```text
>contigA_manual001 original=contigA slice=1-80 strand=+ scaffold=chr1
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>contigB_manual002 original=contigB slice=1-40 strand=+ scaffold=chr1
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
```

## `chromo gafprep` Outputs

Use `chromo gafprep` outputs when you want targeted read-to-graph evidence
without aligning every read to the full graph. The command writes inputs for
GraphAligner and audit tables explaining why each read was selected.

| Output | Description |
| --- | --- |
| `<prefix>.targets.tsv` | One row per review-derived target interval, including source review table, review type, row ID, contig, coordinates, target type, reason, and source fields. |
| `<prefix>.selected_reads.tsv` | One row per deduplicated selected read with best read-to-assembly PAF alignment, selection reasons, linked target IDs, linked review row IDs, and target/background flags. |
| `<prefix>.selected_read_review_links.tsv` | One row per selected-read-to-review-row association with event type, target ID, distance, overlap, and selection reason. |
| `<prefix>.selected_read_ids.txt` | Selected read IDs, one per line. |
| `<prefix>.selected.fastq.gz` | FASTQ subset streamed from `.fastq`, `.fq`, `.fastq.gz`, or `.fq.gz` input. |
| `<prefix>.graphaligner.gfa` | GraphAligner-oriented GFA with sequences preserved, hifiasm `A` records removed, and pathological full-consuming links dropped. |
| `<prefix>.graphaligner.sh` | Executable shell script that runs GraphAligner on the selected FASTQ and sanitized GFA. |
| `<prefix>.summary.tsv` | PAF filter counts, selected/extracted read counts, output paths, and graph-sanitization warnings. |
| `<prefix>.dropped_gfa_links.tsv` | One row per dropped full-consuming `L` link with node names, orientations, lengths, CIGAR, overlap length, and reason. |
| `<prefix>.gfa_sanitize_summary.tsv` | GFA record counters and warnings, including whether sanitization touched a selected target contig. |

If summary files include `GAF evidence limited by graph sanitization near
target`, treat GAF results for that region as limited review evidence.

## `chromo plot` Outputs

Use `chromo plot` outputs to visually check MUMmer coords or minimap2 PAF
alignments without re-running the aligner. For examples of how to read clean,
reversed, chimeric, inverted, duplicated, gapped, and noisy patterns, see
[How to Interpret Dot Plots]({{ '/dot-plots/' | relative_url }}). When
`--sel-ref` is used, the whole-genome output is restricted to the selected
references and `--per-ref` writes only those selected reference panels.

| Output | Description |
| --- | --- |
| `<prefix>.pdf` | Whole-genome PDF dot plot by default. |
| `<prefix>.svg` | Whole-genome SVG dot plot when `--formats svg` is set. |
| `<prefix>.png` | Whole-genome PNG dot plot when `--formats png` is set. |
| `<prefix>.<ref>.<format>` | Per-reference plots when `--per-ref` is set; restricted to selected references when `--sel-ref` is also set. |
| `<prefix>.gfa_overlay.tsv` | Projected GFA interval report when `--gfa-overlay` is used and `--gfa-overlay-report` is not provided. |

**Table 7. Example plot file set.** A single run can write a whole-genome plot
and, with `--per-ref`, one plot per reference sequence in each requested format.

| Output path | Meaning |
| --- | --- |
| `plots/sample.pdf` | Whole-genome PDF dot plot. |
| `plots/sample.svg` | Whole-genome SVG dot plot from the same alignment rows. |
| `plots/sample.chr1.svg` | Per-reference SVG for `chr1`. |
| `plots/sample.chr2.svg` | Per-reference SVG for `chr2`. |

<figure>
  <img src="https://rotheconrad.github.io/chromosort/assets/chromo_plot_example.png" alt="ChromoSort whole-genome dot plot with forward and reverse alignments." width="900">
  <figcaption><strong>Figure 2. Example whole-genome dot-plot output.</strong> PNG output from the synthetic coords fixture showing all reference sequences on the x-axis and assembly contigs on the y-axis; forward alignments are blue and reverse alignments are red.</figcaption>
</figure>

<figure>
  <img src="https://rotheconrad.github.io/chromosort/assets/chromo_plot_example.chr1.png" alt="ChromoSort per-reference dot plot for chr1." width="900">
  <figcaption><strong>Figure 3. Example per-reference dot-plot output.</strong> Per-reference output for <code>chr1</code> generated with <code>--per-ref</code>, useful for inspecting one chromosome-level slice without re-running an aligner.</figcaption>
</figure>

## `chromo graph-map` Outputs

Use `chromo graph-map` outputs when a unitig-level GFA must be projected onto
contig FASTA coordinates before graph evidence can be compared to dot plots or
breakpoint candidates.

| Output | Description |
| --- | --- |
| `<prefix>.utg_to_ctg.tsv` | One row per projected path/walk step with contig coordinates, unitig name, orientation, segment length, overlap fields, source GFA, sequence-present flag, and reuse flag. |
| `<prefix>.path_summary.tsv` | One row per requested path/contig with projected bp, FASTA length, length-difference status, missing/zero-length step counts, and reused segment count. |
| `<prefix>.warnings.tsv` | Structured warnings for absent paths, missing segments, no-sequence/no-`LN:i` lengths, and path-vs-FASTA length mismatches. |

## `chromo scaffold` Outputs

Use `chromo scaffold` outputs after final sorting and review to inspect scaffold
FASTA construction, inferred gaps, overlap handling, and optional graph
junction evidence.

| Output | Description |
| --- | --- |
| `<prefix>.scaffold.fa` | One FASTA record per assigned reference sequence, with ordered contigs joined by Ns. |
| `<prefix>.scaffold.agp` | AGP 2.1 map of emitted scaffold objects, ordered contig components, and written N gaps. |
| `<prefix>.scaffold_components.tsv` | Human-readable provenance table mirroring the AGP rows, including component coordinates, gap lengths, source labels, and status fields. |
| `<prefix>.scaffold_gaps.tsv` | One row per inserted gap with flanking contigs, inferred gap, written gap, overlap bp/class/fractions, overlap policy/action, trimmed bp, and sequence-overlap identity when checked. |
| `<prefix>.graph_gaps.tsv` | Optional report-only GFA evidence for adjacent scaffold junctions when `--gfa` is provided, including direct links, short paths, orientations, overlap bp, and missing-node statuses. |
| `<prefix>.scaffold_summary.tsv` | One row per scaffold with contig count, scaffold length, sequence bp, gap bp, overlap totals, trimming totals, and ordered contig list. |
| `<prefix>.submission_checklist.tsv` | Spreadsheet-friendly submission checklist with output-file manifest rows, FASTA character/terminal-N checks, AGP continuity and FASTA/AGP consistency checks, gap/component counts, unresolved gap totals, and graph-fill span totals. |
| `<prefix>.run_summary.txt` | Inputs, gap model, output paths, and total scaffold counts. |

**Table 8. Example `scaffold_gaps.tsv` row.** The gap report records the flanking
contigs, inferred reference-space gap, FASTA gap actually written, overlap
classification, and overlap policy action.

| scaffold | left_contig | right_contig | raw_inferred_gap_bp | gap_bp | gap_mode | overlap_class | overlap_action |
| --- | --- | --- | ---: | ---: | --- | --- | --- |
| `chr1` | `chr1_contigA` | `chr1_contigB` | `5` | `5` | `inferred` | `no_overlap` | `none` |

**Table 9. Example `scaffold_summary.tsv` rows.** The summary table gives one
row per emitted scaffold record.

| scaffold | contigs | scaffold_bp | sequence_bp | gap_bp | gaps | ordered_contigs |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| `chr1` | `2` | `12` | `7` | `5` | `1` | `chr1_contigA,chr1_contigB` |
| `chr2` | `1` | `2` | `2` | `0` | `0` | `chr2_contigC` |

**Listing 5. Example scaffold FASTA output.** Scaffold headers summarize the
number of source contigs, sequence bases, gap bases, and gap mode.

```text
>chr1 contigs=2 sequence_bp=7 gap_bp=5 gap_mode=inferred
AAAANNNNNTTT
>chr2 contigs=1 sequence_bp=2 gap_bp=0 gap_mode=inferred
GG
```

## `chromo gapfill` Outputs

Use `chromo gapfill` outputs to review graph-supported candidate fills before
optional sequence-changing application.

| Output | Description |
| --- | --- |
| `<prefix>.gapfill_plan.tsv` | One row per adjacent sorted contig pair, or per adjacent AGP component pair in scaffold/AGP mode, with graph status, path nodes, GAF, Hi-C, reference-placement, optional long-read bridge support counts, optional external patch comparison fields, risk flags, branch-complexity score, high-degree/self-loop/unsequenced node lists, fill status, inserted bp, right-trim bp, fallback gap bp, editable `accept_fill`, and whether the fill was applied. With `--project-gfa-paths`, `left_graph_node` and `right_graph_node` are projected terminal unitigs; sequence-bearing projected rows can be `fillable`, while no-sequence rows remain `projected_path_planning_only`. |
| `<prefix>.gapfilled.agp` | AGP 2.1 map of ordered or AGP-derived components, fallback N gaps, and applied graph-fill component slices. In planning mode it describes the all-fallback map; in apply mode it describes `gapfilled.fa`. |
| `<prefix>.gapfilled_components.tsv` | Human-readable provenance table mirroring the AGP rows, including graph-fill segment coordinates, trimmed right-flank coordinates, gap lengths, source labels, and status fields. |
| `--review-html` path | Optional self-contained HTML table for reviewing gapfill-plan rows, comparing candidate paths, and exporting a reviewed-plan TSV. |
| `<prefix>.gapfilled.fa` | Optional FASTA written only with `--apply`, containing one record per assigned reference plus unassigned records. |
| `<prefix>.submission_checklist.tsv` | Spreadsheet-friendly submission checklist with output-file manifest rows, FASTA/AGP checks when a final FASTA exists, planning-mode warnings when it does not, unresolved gap totals, graph-filled span totals, non-ACGTN counts, and external-validation reminder. |
| `<prefix>.run_summary.txt` | Inputs, `graph_mode` (`direct_gfa_nodes` or `projected_gfa_paths`), parameters, output paths, and fill-status counts. |

Patch comparison columns are populated when `--patch-table` is provided:
`patch_candidate_count`, `patch_best_id`, `patch_best_source`, `patch_best_bp`,
`patch_graph_status`, `patch_best_notes`, and, with `--include-fill-sequences`,
`patch_best_sequence`. `patch_graph_status` values include
`exact_graph_match`, `same_length_graph_mismatch`, `graph_mismatch`, and
`patch_only_no_graph_sequence`.

**Table 10. Example `gapfill_plan.tsv` row.** Selected columns from a graph
fixture show a junction resolved by reference-placement PAF support. The default
`accept_fill=no` makes planning review explicit before strict reviewed
application.

| scaffold | left_contig | right_contig | graph_status | path_nodes | candidate_paths | ref_path_support | ref_best_alt_support | risk_flags | fill_status | fill_bp | right_trim_bp | accept_fill | applied |
| --- | --- | --- | --- | --- | ---: | ---: | ---: | --- | --- | ---: | ---: | --- | --- |
| `chr1` | `chr1_left` | `chr1_right` | `ref_paf_resolved_paths` | `left+,bridge_good+,right+` | `2` | `8` | `6` | `branching,high_degree` | `fillable` | `4` | `4` | `no` | `no` |

**Listing 6. Example applied gapfilled FASTA output.** With
`--apply --reviewed-plan`, accepted fillable paths insert graph sequence and
trim the right flank by the terminal GFA overlap; unresolved or unaccepted
junctions use fallback N gaps.

```text
>chr1 contigs=2 filled_gaps=1 fallback_gaps=0 fill_bp=4 fallback_gap_bp=0 trimmed_bp=4
AAAACCCCGGGGTTTT
```

## Audit Tables

The report tables are deliberately redundant. They keep original contig IDs,
new contig IDs, reference assignments, coordinates, status labels, and decision
metrics close to the sequence output they explain. Prefer using these tables
instead of parsing FASTA headers when downstream scripts need assignment,
split, cut, scaffold, or gap metadata.

## Sequence-Changing Outputs

Only a subset of commands change sequence:

- `chromo fix` replaces selected contigs with split pieces.
- `chromo cut` replaces selected contigs with pieces cut at exact reviewed positions.
- `chromo manual apply` reproduces a reviewed browser recipe.
- `chromo scaffold` joins sorted contigs with inferred or fixed N gaps and only trims overlaps when an explicit overlap policy asks it to.
- `chromo gapfill --apply --reviewed-plan` inserts graph sequence only for fillable and accepted paths; `--apply --apply-all-fillable` is the explicit exploratory all-fillable mode.

`chromo eval`, `chromo manual`, `chromo plot`, and report-only graph or
long-read evidence modes do not trim, polish, or rewrite sequence. `chromo sort`
does rewrite retained records into an ordered FASTA and can optionally orient
them with `--orient-to-reference`.
