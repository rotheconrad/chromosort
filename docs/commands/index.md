---
title: Command Reference
description: Overview of ChromoSort subcommands.
---

# Command Reference

ChromoSort exposes one top-level command, `chromo`, with eleven subcommands. Each subcommand page includes the purpose, examples, outputs, key parameters, and reasoning behind the command.

| Command | Use it when you need to |
| --- | --- |
| [`chromo sort`](https://rotheconrad.github.io/chromosort/commands/sort/) | assign contigs to reference sequences, filter duplicate overlaps, and write an ordered FASTA. |
| [`chromo clean`](https://rotheconrad.github.io/chromosort/commands/clean/) | sort-filter raw contigs, conservatively fix retained contigs, and write a cleaned FASTA. |
| [`chromo eval`](https://rotheconrad.github.io/chromosort/commands/eval/) | prepare editable TSV review tables for command-line or spreadsheet-first `fix`, `scaffold`, and `gapfill` curation, including `eval all` for a three-table review bundle. |
| [`chromo fix`](https://rotheconrad.github.io/chromosort/commands/fix/) | split reviewed or automatically detected chimeric contigs into reference-labeled pieces. |
| [`chromo cut`](https://rotheconrad.github.io/chromosort/commands/cut/) | cut contigs at exact reviewed coordinates. |
| [`chromo manual`](https://rotheconrad.github.io/chromosort/commands/manual/) | generate a browser dashboard for manual dot-plot review, task-specific review-event queues, modular evidence panels, and reproducible recipe export. See the [dot-plot guide]({{ '/dot-plots/' | relative_url }}) if the visual patterns are unfamiliar. |
| [`chromo gafprep`](https://rotheconrad.github.io/chromosort/commands/gafprep/) | prepare targeted GraphAligner inputs from read-to-assembly PAF and ChromoSort review tables, commonly the three tables from `chromo eval all`. |
| [`chromo plot`](https://rotheconrad.github.io/chromosort/commands/plot/) | draw whole-genome, per-reference, or selected-reference dot plots from existing MUMmer coords or minimap2 PAF alignments. Use the [dot-plot guide]({{ '/dot-plots/' | relative_url }}) to interpret the patterns. |
| [`chromo graph-map`](https://rotheconrad.github.io/chromosort/commands/graph-map/) | project GFA unitig path/walk coordinates onto matching contig FASTA coordinates and write graph-map reports. |
| [`chromo scaffold`](https://rotheconrad.github.io/chromosort/commands/scaffold/) | join final sorted contigs into one scaffold per reference sequence, with optional reviewed gap overrides and report-only graph junction evidence. |
| [`chromo gapfill`](https://rotheconrad.github.io/chromosort/commands/gapfill/) | plan and optionally apply reviewed graph-supported fills between adjacent sorted contigs using guarded GFA paths and optional support evidence. |
