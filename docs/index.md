---
title: Documentation
description: User guide for ChromoSort.
---

# ChromoSort Documentation

ChromoSort is a reference-guided assembly curation toolkit for sorting contigs, splitting reviewed chimeric contigs, making manual cuts, reviewing dot plots, scaffolding ordered contigs, applying reviewed graph-supported gap fills, and plotting existing whole-genome alignments.

This documentation is organized as a user guide rather than one very long README. The archived longform README remains available for historical reference.


## Start Here

- [Installation](https://rotheconrad.github.io/chromosort/installation/)
- [Input files](https://rotheconrad.github.io/chromosort/input-files/)
- [Workflows](https://rotheconrad.github.io/chromosort/workflows/)
- [Command reference](https://rotheconrad.github.io/chromosort/commands/)

## Commands

- [chromo sort](https://rotheconrad.github.io/chromosort/commands/sort/)
- [chromo fix](https://rotheconrad.github.io/chromosort/commands/fix/)
- [chromo cut](https://rotheconrad.github.io/chromosort/commands/cut/)
- [chromo manual](https://rotheconrad.github.io/chromosort/commands/manual/)
- [chromo plot](https://rotheconrad.github.io/chromosort/commands/plot/)
- [chromo scaffold](https://rotheconrad.github.io/chromosort/commands/scaffold/)
- [chromo gapfill](https://rotheconrad.github.io/chromosort/commands/gapfill/)

## Running Analyses

- [Output files](https://rotheconrad.github.io/chromosort/outputs/)
- [Architecture](https://rotheconrad.github.io/chromosort/architecture/)
- [Troubleshooting](https://rotheconrad.github.io/chromosort/troubleshooting/)

## Project Information

- [Current status and version history](https://rotheconrad.github.io/chromosort/status/)
- [Contributing](https://rotheconrad.github.io/chromosort/contributing/)
- [Archived longform README](https://rotheconrad.github.io/chromosort/archive/longform-readme/)

## Main Workflow

1. Generate MUMmer `show-coords` or minimap2 PAF alignments between the reference and assembly.
2. Use `chromo plot` and, when helpful, `chromo manual` to inspect suspicious contigs.
3. Use `chromo fix` or `chromo cut` for reviewed sequence edits.
4. Re-align the fixed assembly.
5. Use `chromo sort` to assign, filter, orient, and order contigs.
6. Use `chromo scaffold` to build one scaffold per reference sequence.
7. Use `chromo gapfill` only when reviewed GFA paths should replace N gaps.
