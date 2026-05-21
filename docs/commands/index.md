---
title: Command Reference
description: Overview of ChromoSort subcommands.
---

# Command Reference

ChromoSort exposes one top-level command, `chromo`, with seven subcommands. Each subcommand page includes the purpose, examples, outputs, key parameters, and reasoning behind the command.

| Command | Use it when you need to |
| --- | --- |
| [`chromo sort`](https://rotheconrad.github.io/chromosort/commands/sort/) | assign contigs to reference sequences, filter duplicate overlaps, and write an ordered FASTA. |
| [`chromo fix`](https://rotheconrad.github.io/chromosort/commands/fix/) | split reviewed or automatically detected chimeric contigs into reference-labeled pieces. |
| [`chromo cut`](https://rotheconrad.github.io/chromosort/commands/cut/) | cut contigs at exact reviewed coordinates. |
| [`chromo manual`](https://rotheconrad.github.io/chromosort/commands/manual/) | generate a browser dashboard for manual dot-plot review and reproducible recipe export. |
| [`chromo plot`](https://rotheconrad.github.io/chromosort/commands/plot/) | draw dot plots from existing MUMmer coords or minimap2 PAF alignments. |
| [`chromo scaffold`](https://rotheconrad.github.io/chromosort/commands/scaffold/) | join final sorted contigs into one scaffold per reference sequence. |
| [`chromo gapfill`](https://rotheconrad.github.io/chromosort/commands/gapfill/) | plan and optionally apply reviewed graph-supported fills between adjacent sorted contigs. |
