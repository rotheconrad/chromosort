# ChromoSort

Reference-guided genome assembly utilities for sorting contigs, conservatively cleaning mostly-correct assemblies, splitting reviewed chimeric contigs, cutting exact breakpoints, manual dot-plot review, plotting alignments, scaffolding final ordered contigs, and applying reviewed graph-supported gap fills.

ChromoSort provides one command, `chromo`, with nine subcommands:

| Command | Purpose |
| --- | --- |
| `chromo sort` | Assign contigs to the best-supported reference sequence from MUMmer coords or minimap2 PAF, merge alignment evidence, filter contained or low-value duplicate overlaps, protect likely split candidates, and write a reference-ordered FASTA with TSV decision reports ([sort docs](https://rotheconrad.github.io/chromosort/commands/sort/)). |
| `chromo clean` | Apply sort-style filtering to raw contigs, conservatively fix retained contigs, orient/order the emitted records, and write a cleaned FASTA plus audit reports for mostly-correct assemblies ([clean docs](https://rotheconrad.github.io/chromosort/commands/clean/)). |
| `chromo eval` | Prepare editable TSV review tables for algorithm-assisted, human-reviewed decisions, including `fix`, `scaffold`, and `gapfill` tables for the matching `--reviewed-plan` execution paths ([eval docs](https://rotheconrad.github.io/chromosort/commands/eval/)). |
| `chromo fix` | Split chimeric or structurally inconsistent contigs into reference-labeled pieces by scanning query-ordered alignment blocks, smoothing ordinary gaps, selecting eligible reference/orientation transitions, and writing a fixed full-assembly FASTA plus an audit report ([fix docs](https://rotheconrad.github.io/chromosort/commands/fix/)). |
| `chromo cut` | Apply exact reviewed breakpoint edits when you already know the cut positions, replacing each requested contig with numbered pieces while copying uncut contigs unchanged and recording every emitted slice ([cut docs](https://rotheconrad.github.io/chromosort/commands/cut/)). |
| `chromo manual` | Build a self-contained browser dashboard for dot-plot curation, contig removal/restoration, reordering, inversion, breakpoint staging, optional GFA neighborhood review, FASTA export, and reproducible recipe application ([manual docs](https://rotheconrad.github.io/chromosort/commands/manual/), [dot-plot guide](https://rotheconrad.github.io/chromosort/dot-plots/)). |
| `chromo plot` | Draw whole-genome, per-reference, or selected-reference dot plots from existing MUMmer coords or minimap2 PAF alignments, optionally ordered by a `chromo sort` assignment report, without re-running an aligner ([plot docs](https://rotheconrad.github.io/chromosort/commands/plot/), [dot-plot guide](https://rotheconrad.github.io/chromosort/dot-plots/)). |
| `chromo scaffold` | Join the final sorted contigs into one scaffold FASTA record per assigned reference, infer or fix N-gap lengths, report overlaps and gap decisions, and optionally add report-only GFA junction evidence ([scaffold docs](https://rotheconrad.github.io/chromosort/commands/scaffold/)). |
| `chromo gapfill` | Plan graph-supported fills between adjacent sorted contigs using GFA paths plus optional GAF, Hi-C-like, or reference-placement PAF evidence, then apply only fillable and reviewed paths while unresolved junctions fall back to N gaps ([gapfill docs](https://rotheconrad.github.io/chromosort/commands/gapfill/)). |

## Documentation

Full documentation is available at <https://rotheconrad.github.io/chromosort/>.

New users should start with [Installation](https://rotheconrad.github.io/chromosort/installation/), then use [Input Files](https://rotheconrad.github.io/chromosort/input-files/) to prepare MUMmer, minimap2, GFA, GAF, or Hi-C-like evidence. The [Workflows](https://rotheconrad.github.io/chromosort/workflows/) page shows the recommended order for fixing, sorting, plotting, scaffolding, and graph-aware review. The [dot-plot guide](https://rotheconrad.github.io/chromosort/dot-plots/) is a mini tutorial for reading whole-genome and per-reference dot plots. The [Production Upgrade Roadmap](https://rotheconrad.github.io/chromosort/roadmap/) captures the planned `eval` and task-specific `manual` review modes. Command-specific pages are linked in the table above.

For interpreting results, see [Output Files](https://rotheconrad.github.io/chromosort/outputs/) and [Troubleshooting](https://rotheconrad.github.io/chromosort/troubleshooting/).

## Alignment Evidence Matches One FASTA

MUMmer coords and minimap2 PAF files describe one exact reference FASTA and one
exact assembly FASTA. If a ChromoSort step writes a changed FASTA by removing
records, splitting contigs, cutting contigs, reverse-complementing records,
renaming records, or scaffolding records, re-run MUMmer or minimap2 before using
that changed FASTA as the assembly input to another alignment-dependent command.

You can reuse the original coords or PAF to make decisions about the original
assembly. For example, run `chromo sort` on `raw.fa`, inspect
`split_candidate=yes` rows, then run `chromo fix` on that same `raw.fa` with the
same raw alignment file. You should not run `chromo fix` on
`sample.ordered.fa` from `chromo sort` with coords that were generated from
`raw.fa`.

`chromo plot --assignments` is also an important special case: it plots the
original alignment rows while ordering the query axis by a `chromo sort`
assignment report. This is useful for reviewing sort decisions without
re-aligning, but it is not a new alignment of the edited FASTA. To validate
`ordered.fa`, `fixed.fa`, or a manual-export FASTA, generate fresh coords or PAF
for that exact FASTA. For help reading the resulting visual patterns, use the
[dot-plot guide](https://rotheconrad.github.io/chromosort/dot-plots/).

## Quick Start

```bash
git clone https://github.com/rotheconrad/chromosort.git
cd chromosort

mamba env create -f environment.yml
mamba activate chromosort

chromo --help
chromo sort --help
chromo clean --help
```

Typical conservative cleanup workflow for a mostly-correct assembly:

```bash
chromo clean \
  --ref-fasta reference.fa \
  --assembly-fasta assembly.fa \
  --coords mummer/raw.coords \
  --output-prefix results/sample \
  --orient-to-reference \
  --discarded-fasta results/sample.discarded.fa

# Re-align results/sample.clean.fa before final validation plots.
```

Typical reviewed workflow, with re-alignment after FASTA edits:

```bash
# 1. Fix reviewed/suspect raw contigs.
chromo fix \
  --assembly-fasta assembly.fa \
  --coords mummer/raw.coords \
  --contigs suspect_contig_1 suspect_contig_2 \
  --output-fasta results/sample.fixed.fa \
  --report results/sample.fixed_contigs.tsv

# 2. Re-align results/sample.fixed.fa with MUMmer or minimap2.

# 3. Sort the fixed FASTA with the fixed-FASTA alignment.
chromo sort \
  --ref-fasta reference.fa \
  --assembly-fasta results/sample.fixed.fa \
  --coords mummer/fixed.coords \
  --output-prefix results/sample.fixed \
  --orient-to-reference

# 4. Plot from the same fixed-FASTA alignment for visual review.
chromo plot \
  --ref-fasta reference.fa \
  --assembly-fasta results/sample.fixed.fa \
  --coords mummer/fixed.coords \
  --assignments results/sample.fixed.contig_assignments.tsv \
  --output-prefix plots/sample.fixed \
  --per-ref

# Add --sel-ref Gm6 Gm12 Gm15 to redraw only selected references.
```

## Install With Pixi

```bash
git clone https://github.com/rotheconrad/chromosort.git
cd chromosort

pixi install
pixi run help
pixi run test
```

## Current Status

Current version: `0.2.25`. Operational commands are `sort`, `clean`, `eval`, `fix`, `cut`, `manual`, `plot`, `scaffold`, and `gapfill`. See [`docs/status.md`](docs/status.md) or [`CHANGELOG.md`](CHANGELOG.md) for version history. See [`docs/roadmap.md`](docs/roadmap.md) for the planned task-specific `manual` review-mode upgrade.

## Citation

If you use ChromoSort, cite this repository and cite MUMmer or minimap2 for the alignment files used by the workflow. See [`CITATION.cff`](CITATION.cff).

## Contact

Please use the GitHub issue tracker for bug reports, feature requests, and questions: <https://github.com/rotheconrad/chromosort/issues>.

## License

ChromoSort is released under the MIT License. See [`LICENSE`](LICENSE).

## Funding Support

This project is supported by the U.S. Department of Agriculture - Agricultural
Research Service (USDA-ARS) - Genomics and Bioinformatics Research Unit (GBRU)
through CRIS Project No. 6066-21310-006-000-D.

## Acknowledgements

ChromoSort can consume MUMmer and minimap2 whole-genome alignments. Thanks to
the genome assembly and comparative genomics communities whose workflows
motivated transparent reference-guided contig sorting, splitting, plotting, and
scaffolding tools.

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
