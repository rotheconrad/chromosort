---
title: Current Status and Roadmap
description: Current ChromoSort command status and version history.
---

# Current Status and Roadmap

ChromoSort is under active development. Current version: `0.2.25`.

## What Works Now

- `chromo sort`
- `chromo clean`
- `chromo fix`
- `chromo cut`
- `chromo manual`
- `chromo plot`
- `chromo scaffold`
- `chromo gapfill`

For the plotting and manual review commands, the
[dot-plot guide]({{ '/dot-plots/' | relative_url }}) explains how to interpret
the visual patterns those tools show.

The repository includes synthetic tests and fixtures under `tests/`, including small FASTA, MUMmer coords, minimap2 PAF, GFA, GAF, and Hi-C-like files for graph-aware workflows.

## Production Upgrade Roadmap

The next planned production upgrade is captured in the
[Production Upgrade Roadmap]({{ '/roadmap/' | relative_url }}). The roadmap
describes planned `eval` table workflows and task-specific `manual` dashboards
for reviewed `fix`, `scaffold`, and `gapfill` decisions.

## Development Checks

```bash
pytest
python -m unittest discover -s tests -v
pixi run test
```

## Version History

| Version | Notes |
| --- | --- |
| Unreleased | No changes yet. |
| `0.2.25` | Synchronized package, citation, Pixi, conda recipe, README, and docs version metadata; added the production-upgrade roadmap for paired `eval` table workflows and task-specific `manual` dashboards feeding reviewed `fix`, `scaffold`, and `gapfill` execution paths. |
| `0.2.24` | Added `chromo clean`, a conservative cleanup command for mostly-correct assemblies that combines sort-style filtering with fix-style conservative splitting on retained raw contigs, then writes `<prefix>.clean.fa` plus initial-sort, fix, clean, and run-summary reports. Clarified README, command docs, and workflows around when FASTA-changing steps require fresh MUMmer or minimap2 alignments before downstream steps or final plots. |
| `0.2.23` | Renamed the graph gap-filling command from `chromo fill` to `chromo gapfill`, moved the package entry point to `chromosort.gapfill`, replaced the package script with `chromosort-gapfill`, and updated gapfill output names to `<prefix>.gapfill_plan.tsv` and `<prefix>.gapfilled.fa`. |
| `0.2.22` | Added Pixi installation support with `pixi.toml`, plus README figure assets and captions for `chromo manual` graph review and `chromo plot` whole-genome/per-reference examples. |
| `0.2.21` | Added graph-aware safety policies. `chromo sort` and `chromo fix` now have warning-only `--graph-guard` checks, while `chromo scaffold --graph-overlap-policy report|warn|confirm` keeps graph evidence report-only by default and only lets direct oriented GFA links confirm narrow terminal-overlap trimming when explicitly requested. |
| `0.2.20` | Added an end-to-end synthetic graph workflow to the README and shipped focused gapfill walkthrough inputs. The tutorial runs sort/manual/scaffold/gapfill with the graph-gotcha GFA, PAF, GAF, Hi-C-like contacts, review HTML, reviewed-plan TSV, and reviewed gapfill application. |
| `0.2.19` | Improved `chromo gapfill --review-html` candidate comparison. Review dashboards now embed per-candidate path rows with path nodes, support scores, validation status, fill length, trim length, risk flags, and optional fill sequence so reviewers can compare ambiguous branches directly before exporting a reviewed plan. |
| `0.2.18` | Added richer path-risk annotations to `chromo gapfill`. Gapfill plans and review HTML now report risk flags, branch-complexity score, high-degree graph nodes, self-loop nodes, unsequenced nodes, and cycle-guard counts so ambiguous or risky candidate paths are easier to triage. |
| `0.2.17` | Added reference-placement PAF evidence to `chromo gapfill`. The new `--ref-paf` path scorer reports selected and best-alternate reference support, can conservatively resolve ambiguous branches when one candidate has unique expected-gap placement support, and conflicts with GAF or Hi-C support leave the gap unresolved. |
| `0.2.16` | Expanded `chromo manual --gfa` review. Manual dashboards now include graph-neighborhood filtering, a selected-contig upstream/downstream neighbor panel, overlap/orientation details, and same-reference neighbor flags so branching graph context is easier to compare during manual curation. |
| `0.2.15` | Added `chromo manual --gfa` graph context. Manual dashboards now embed per-contig GFA node evidence, graph complexity labels, degree/neighbor counts, coverage tags such as `RC:i`, and oriented neighbor summaries so manual breakpoint and ordering review can consider local assembly-graph structure. |
| `0.2.14` | Added `chromo gapfill --review-html`, a self-contained HTML review table for gapfill plans. It embeds the same TSV columns, supports filtering and accepted-fill toggles, and exports a reviewed-plan TSV for `--reviewed-plan`; the TSV and HTML writers now share one row-generation path. |
| `0.2.13` | Added reviewed gapfill-plan application for `chromo gapfill`. Planning output now includes an editable `accept_fill` column, and `--reviewed-plan` makes `--apply` fill only accepted rows after rechecking the current scaffold, contig pair, path nodes, and fillability; rejected or unaccepted rows fall back to N gaps. |
| `0.2.12` | Added optional Hi-C pair support to `chromo gapfill`. Gapfill plans now report Hi-C path support and best alternate support, and otherwise ambiguous graph branches can be resolved when one candidate has unique summed contact support at or above `--min-hic-path-support`; conflicting GAF and Hi-C support leaves the junction unresolved. |
| `0.2.11` | Expanded the input-file documentation with a dedicated graph-input section describing where to find matching GFA files, which reference-to-assembly PAF files to keep for raw and fixed FASTAs, and how optional GAF read-to-graph alignments are used by `chromo gapfill`. |
| `0.2.10` | Added optional GAF read-path evidence to `chromo gapfill`. Gapfill plans now report GAF support counts, and otherwise ambiguous graph branches can be resolved when one candidate path has unique support after `--min-gaf-mapq` filtering and meets `--min-gaf-path-support`; weak, tied, or missing support still leaves the junction unresolved. |
| `0.2.9` | Added `chromo gapfill`, a conservative graph-gap planning and optional application command. It writes `<prefix>.gapfill_plan.tsv`, refuses ambiguous or unverifiable GFA paths, applies sequence only with `--apply`, trims the right flank by the final graph overlap when filling, and falls back to inferred or fixed N gaps for unresolved junctions. |
| `0.2.8` | Added report-only `--gfa` graph context to `chromo sort` and `chromo fix`. Sorting now writes `<prefix>.graph_assignments.tsv` with resolved graph nodes, node degree/self-loop evidence, and direct links to overlap-best contigs; fixing now writes a graph context table beside the split report so reviewed contigs can be checked against the assembly graph before gapfill workflows. |
| `0.2.7` | Added `chromo scaffold --gfa` report-only graph evidence. When a GFA is provided, scaffolding now writes `<prefix>.graph_gaps.tsv` with resolved graph nodes, orientation-aware direct links, link overlap bp, short explicit GFA paths up to `--graph-max-path-edges`, intermediate candidate nodes, and missing/no-path statuses without changing FASTA output. |
| `0.2.6` | Added the first graph-evidence foundation: a tested GFA parser for segment/link records, orientation-aware edge lookup helpers, overlap-CIGAR handling that preserves complex overlaps as non-trim lengths, and synthetic graph-gotcha fixtures with GFA, PAF, GAF, Hi-C-like, and expected-path files for future roadmap development. |
| `0.2.5` | Added `chromo manual`, a self-contained HTML dashboard for manual dot-plot review, contig removal/restoration, order changes, breakpoints, inversions, scaffold labeling/export, FASTA downloads, recipe JSON export, and reproducible `chromo manual apply` recipe execution. |
| `0.2.4` | Added `chromo cut` for exact reviewed breakpoint cuts, with repeatable `--cut CONTIG:POS[,POS...]`, single-contig `--contig/--pos`, batch `--cuts-file`, cut-piece FASTA output, and an audit TSV report. |
| `0.2.3` | Added explicit terminal-overlap classification/rescue in `chromo sort`, richer scaffold overlap reporting, and `chromo scaffold --overlap-policy` modes for warn-only, reference-coordinate trimming, and sequence-confirmed trimming. |
| `0.2.2` | Reworked `chromo fix` so `--contigs`/`--contigs-file` only select the inspection subset, `--all` scans every candidate contig, `--mode` controls planner behavior for both scopes, and breakpoint limits apply per contig. |
| `0.2.1` | Tightened `chromo sort` duplicate filtering for contaminated/alternate-fragment assemblies by using span-based overlap by default, requiring both novel coverage thresholds, rescuing very large near-threshold alignments, and letting split candidates protect their secondary reference spans. |
| `0.2.0` | Added minimap2 PAF input for `chromo sort` and `chromo fix`, plus `chromo plot` PDF/SVG/PNG dot plots for coords/PAF with optional assignment-report query ordering. |
| `0.1.2` | Raised the default auto-split query-span support threshold to 5% so small terminal off-target blocks are reported for review instead of being cut automatically. |
| `0.1.1` | Tightened `chromo fix` breakpoint placement by collapsing adjacent same-reference/orientation runs, added complex same-reference orientation detection, added a run-level auto breakpoint budget, protected strong multi-reference split candidates during `chromo sort`, and documented the fix-before-sort workflow for suspected misjoins. |
| `0.1.0` | Initial public package with `chromo sort`, `chromo fix`, `chromo scaffold`, duplicate-overlap filtering, user-nominated contig splitting, conservative auto smoothing, inferred/fixed-gap scaffolding, and synthetic tests. |
