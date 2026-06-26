---
title: Educational Guide Roadmap
description: Roadmap for building ChromoSort concept guides modeled on the dot-plot tutorial.
---

# Educational Guide Roadmap

ChromoSort already has strong command references and architecture notes. The
[dot-plot guide]({{ '/dot-plots/' | relative_url }}) works because it teaches a
mental model before it asks readers to choose a command: axes and segments
first, pattern gallery next, then review workflow, common traps, and the next
ChromoSort action.

This roadmap turns that style into a guides section. Each guide should help a
reader decide what they are seeing, which evidence files are trustworthy for
that decision, and which ChromoSort command should or should not change
sequence.

## Available Guides

Start with the foundation, core decision, review-interface, evidence, and
worked-example guides when you are setting up or checking a workflow:

- [Alignment evidence and the exact FASTA rule]({{ '/guides/alignment-evidence/' | relative_url }})
- [Choosing PAF or MUMmer coords]({{ '/guides/paf-vs-coords/' | relative_url }})
- [FASTA and evidence name matching]({{ '/guides/name-matching/' | relative_url }})
- [Reading ChromoSort audit tables]({{ '/guides/audit-tables/' | relative_url }})
- [Sorting decisions and duplicate-overlap filtering]({{ '/guides/sorting-decisions/' | relative_url }})
- [Sort, clean, fix, cut, or manual?]({{ '/guides/choosing-commands/' | relative_url }})
- [Chimeric contig and breakpoint review]({{ '/guides/breakpoint-review/' | relative_url }})
- [Inversions and orientation changes]({{ '/guides/inversions-orientation/' | relative_url }})
- [Manual dashboard review]({{ '/guides/manual-dashboard-review/' | relative_url }})
- [Spreadsheet review tables]({{ '/guides/spreadsheet-review-tables/' | relative_url }})
- [Scaffolding, gaps, and overlaps]({{ '/guides/scaffolding-gaps-overlaps/' | relative_url }})
- [Assembly graph evidence]({{ '/guides/assembly-graph-evidence/' | relative_url }})
- [hifiasm unitig-to-contig projection]({{ '/guides/hifiasm-graph-projection/' | relative_url }})
- [Long-read PAF and GAF support]({{ '/guides/long-read-paf-gaf/' | relative_url }})
- [Graph-supported gap filling]({{ '/guides/graph-gap-filling/' | relative_url }})
- [Mostly-correct assembly cleanup]({{ '/guides/mostly-correct-cleanup/' | relative_url }})
- [Suspected chimeric contig review]({{ '/guides/chimeric-contig-walkthrough/' | relative_url }})
- [Graph-aware scaffold and gapfill review]({{ '/guides/graph-scaffold-gapfill-walkthrough/' | relative_url }})
- [Dataset handoff checklist]({{ '/guides/dataset-handoff-checklist/' | relative_url }})
- [How to interpret dot plots]({{ '/dot-plots/' | relative_url }})

## Dot-Plot Guide Review

The dot-plot guide should be the template for future guides because it:

- Defines the visual grammar before showing commands.
- Uses a pattern gallery for common biological and technical situations.
- Repeats the rule that plots describe one exact reference/assembly FASTA pair.
- Separates observation from interpretation, then interpretation from action.
- Links each pattern to the relevant ChromoSort command and report fields.
- Includes common traps so readers do not overcall repeats, inversions,
  reference differences, or stale alignments.

Future guides should use the same shape:

1. Core idea.
2. What ChromoSort reads or writes.
3. Pattern, status, or evidence gallery.
4. Practical review workflow.
5. Cheat sheet.
6. Common traps.
7. What to look at next in ChromoSort.

## Concept Inventory

The current docs, command pages, source modules, and tests point to these
guide-worthy concepts.

| Guide candidate | Main source material | Why it needs a guide |
| --- | --- | --- |
| Alignment evidence and the exact FASTA rule | README, input files, workflows, command docs | Users must know when coords/PAF evidence is still valid and when edited FASTA outputs require fresh alignments. |
| Choosing PAF or MUMmer coords | README, input files, `reference_order.py`, `fix_contigs.py`, plot docs | Users need a practical way to choose one primary alignment, tune minimap2 presets, and interpret expected PAF-vs-coords differences. |
| FASTA, PAF, GFA, GAF, and name matching | Input files, troubleshooting, `graph.py`, `gaf.py`, `graph_map.py` | Most failures come from IDs and coordinate systems drifting between files. |
| Reading ChromoSort audit tables | Outputs, sort/fix/scaffold/gapfill docs, `review.py` | TSV reports are central to trust, but users need help translating status labels into decisions. |
| Sorting decisions and duplicate-overlap filtering | Sort docs, `reference_order.py`, clean docs | `kept`, `duplicate_overlap`, `terminal_overlap`, and `kept_split_candidate` need a concept-level explanation separate from parameter reference. |
| Sort, clean, fix, cut, or manual? | Workflows, clean/fix/cut/manual docs | Users need a decision guide for choosing the least invasive command that fits the evidence. |
| Chimeric contigs and breakpoint review | Fix docs, eval docs, `fix_contigs.py`, tests/data/noisy_fix | The planner modes, smoothing, breakpoint penalties, and manual-review boundaries need an educational walkthrough. |
| Inversions and orientation changes | Dot plots, review playbook, fix docs | Users need to distinguish reverse-complemented contigs, true inversions, reference differences, and events that should not be split. |
| Manual dashboards and reviewed recipes | Manual docs, eval docs, `manual.py`, `review.py` | Browser review is powerful enough to deserve a user-centered guide from dashboard to reproducible recipe. |
| Spreadsheet review tables | Eval docs, reviewed-plan paths in fix/scaffold/gapfill | Users need a table-first guide for accept/reject fields, stale-row protection, and executor validation. |
| Scaffolding, gaps, and overlaps | Scaffold docs, outputs, `scaffold.py` | Negative inferred gaps, terminal overlaps, fixed gaps, and trim policies are common decision points. |
| Assembly graph evidence | Input files, graph-map docs, graph-aware command docs, `graph.py` | GFA segments, links, paths, walks, noseq graphs, and graph evidence policies need a gentle concept guide. |
| hifiasm unitig-to-contig projection | Graph-map docs, workflows, `graph_map.py`, plot overlay logic | Unitig GFA coordinates are not contig FASTA coordinates; this deserves a focused guide with diagrams. |
| Long-read PAF and GAF evidence | Eval/manual/gapfill docs, `longreads.py`, `gaf.py` | Readers need to understand bridge support, breakpoint support, traversal support, MAPQ filters, and report-only versus sequence-changing evidence. |
| Graph-supported gap filling | Gapfill docs, synthetic graph workflow, `gapfill.py` | Fillable, ambiguous, conflicting, unsequenced, and stale paths need a pattern gallery before users apply sequence. |
| End-to-end review handoffs | Review playbook, workflows | Users and assistants need a guide for manifests, evidence bundles, and reproducible handoffs between projects or chats. |

## Guide Build Roadmap

### Phase 1: Foundations

Build the guides that prevent the most expensive mistakes.

| Priority | Guide | Deliverable |
| --- | --- | --- |
| 1 | [Alignment evidence and the exact FASTA rule]({{ '/guides/alignment-evidence/' | relative_url }}) | Available. Explains raw, fixed, sorted, manual, scaffolded, and gapfilled FASTA stages; when old coords/PAF can be reused; and when re-alignment is mandatory. |
| 2 | [Choosing PAF or MUMmer coords]({{ '/guides/paf-vs-coords/' | relative_url }}) | Available. Compares minimap2 PAF and MUMmer coords, including `asm5`/`asm10`/`asm20`, `-c`, secondary rows, MAPQ, identity, and expected disagreement. |
| 3 | [FASTA and evidence name matching]({{ '/guides/name-matching/' | relative_url }}) | Available. Gives quick checks for FASTA IDs, PAF columns, coords query/reference names, GFA `S` names, GAF path nodes, and Hi-C node names. |
| 4 | [Reading ChromoSort audit tables]({{ '/guides/audit-tables/' | relative_url }}) | Available. Covers `contig_assignments.tsv`, fix reports, review-event TSVs, scaffold gap reports, graph reports, and gapfill plans. |

### Phase 2: Core Decision Guides

Build guides that help users pick the right ChromoSort action.

| Priority | Guide | Deliverable |
| --- | --- | --- |
| 5 | [Sorting and duplicate-overlap decisions]({{ '/guides/sorting-decisions/' | relative_url }}) | Available. Explains best-reference share, query coverage, novel reference span, terminal overlap rescue, split-candidate protection, and graph guardrails. |
| 6 | [Sort, clean, fix, cut, or manual?]({{ '/guides/choosing-commands/' | relative_url }}) | Available. Provides a decision tree for choosing conservative cleanup, automatic sorting, reviewed splitting, exact cuts, or browser curation. |
| 7 | [Chimeric contig and breakpoint review]({{ '/guides/breakpoint-review/' | relative_url }}) | Available. Covers `chromo fix`, planner modes, smoothing, breakpoint penalties, and selected-contig versus `--all` workflows. |
| 8 | [Inversions and orientation changes]({{ '/guides/inversions-orientation/' | relative_url }}) | Available. Distinguishes reverse orientation, internal inversion, complex same-reference orientation events, and biological/reference differences. |

### Phase 3: Review Interfaces

Build guides for human-in-the-loop workflows.

| Priority | Guide | Deliverable |
| --- | --- | --- |
| 9 | [Manual dashboard review]({{ '/guides/manual-dashboard-review/' | relative_url }}) | Available. Covers browsing dot plots, staging breakpoints, removing/restoring contigs, labeling scaffolds, exporting FASTA, and applying recipes. |
| 10 | [Spreadsheet review tables]({{ '/guides/spreadsheet-review-tables/' | relative_url }}) | Available. Covers `chromo eval fix/scaffold/gapfill`, accepted rows, stale-row validation, reviewed executor paths, and task-specific manual queues. |
| 11 | [Scaffolding, gaps, and overlaps]({{ '/guides/scaffolding-gaps-overlaps/' | relative_url }}) | Available. Covers inferred gaps, fixed gaps, negative gaps, terminal versus internal overlaps, sequence-confirmed trimming, and reviewed gap overrides. |

### Phase 4: Graph And Long-Read Evidence

Build guides for optional evidence streams that are powerful but easy to
overinterpret.

| Priority | Guide | Deliverable |
| --- | --- | --- |
| 12 | [Assembly graph evidence]({{ '/guides/assembly-graph-evidence/' | relative_url }}) | Available. Explains GFA `S`, `L`, `P`, and `W` records, noseq graphs, graph complexity labels, report-only evidence, and graph guardrails. |
| 13 | [hifiasm unitig-to-contig projection]({{ '/guides/hifiasm-graph-projection/' | relative_url }}) | Available. Explains why unitig coordinates need projection, how `chromo graph-map` and plot overlays use path/walk records, and how to review projection warnings. |
| 14 | [Long-read PAF and GAF support]({{ '/guides/long-read-paf-gaf/' | relative_url }}) | Available. Explains breakpoint support, contig-end bridge support, read traversal support, MAPQ thresholds, and advisory versus executor-changing evidence. |
| 15 | [Graph-supported gap filling]({{ '/guides/graph-gap-filling/' | relative_url }}) | Available. Explains path enumeration, candidate paths, GAF/Hi-C/reference-placement support, risk flags, flank validation, reviewed application, and unresolved fallbacks. |

### Phase 5: Worked Examples

Finish the section with example-driven walkthroughs that combine concepts.

| Priority | Guide | Deliverable |
| --- | --- | --- |
| 16 | [Mostly-correct assembly cleanup]({{ '/guides/mostly-correct-cleanup/' | relative_url }}) | Available. Guides `chromo clean` from raw FASTA to validation plots, including when to re-align. |
| 17 | [Suspected chimeric contig review]({{ '/guides/chimeric-contig-walkthrough/' | relative_url }}) | Available. Walks from dot-plot signal to `eval fix`, manual review, reviewed split, re-alignment, sort, and plot. |
| 18 | [Graph-aware scaffold and gapfill review]({{ '/guides/graph-scaffold-gapfill-walkthrough/' | relative_url }}) | Available. Uses the synthetic graph fixture to compare GFA, GAF, Hi-C, reference-placement PAF, and final reviewed application. |
| 19 | [Dataset handoff checklist]({{ '/guides/dataset-handoff-checklist/' | relative_url }}) | Available. Packages FASTA, coords/PAF, GFA/GAF, review tables, plots, run summaries, and notes for a new review session. |

## Definition Of Done

Each guide should be considered ready when it has:

- A short opening that names the decision the reader is trying to make.
- At least one concrete example table, plot, diagram, or command sequence.
- Links to the command reference pages that implement the action.
- Links to the source report fields or status labels used in the decision.
- A common-traps section.
- A final "what to look at next" section.
- Successful `bundle exec jekyll build --source docs --destination docs/_site --config docs/_config.yml`
  and `git diff --check` checks.

Guides that include diagrams should use the same rendering safeguards as the
architecture docs: set page front matter for any required rendering, keep labels
short, and verify the built site output rather than relying only on Markdown
preview.
