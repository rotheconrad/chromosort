---
title: Command Reference
description: Overview of ChromoSort subcommands.
---

# Command Reference

ChromoSort exposes one top-level command, `chromo`, with seven subcommands. Each subcommand page includes the purpose, examples, outputs, key parameters, and reasoning copied from the longform README.

| Command | Use it when you need to |
| --- | --- |
| [`chromo sort`]({{ '/commands/sort/' | relative_url }}) | assign contigs to reference sequences, filter duplicate overlaps, and write an ordered FASTA. |
| [`chromo fix`]({{ '/commands/fix/' | relative_url }}) | split reviewed or automatically detected chimeric contigs into reference-labeled pieces. |
| [`chromo cut`]({{ '/commands/cut/' | relative_url }}) | cut contigs at exact reviewed coordinates. |
| [`chromo manual`]({{ '/commands/manual/' | relative_url }}) | generate a browser dashboard for manual dot-plot review and reproducible recipe export. |
| [`chromo plot`]({{ '/commands/plot/' | relative_url }}) | draw dot plots from existing MUMmer coords or minimap2 PAF alignments. |
| [`chromo scaffold`]({{ '/commands/scaffold/' | relative_url }}) | join final sorted contigs into one scaffold per reference sequence. |
| [`chromo gapfill`]({{ '/commands/gapfill/' | relative_url }}) | plan and optionally apply reviewed graph-supported fills between adjacent sorted contigs. |
