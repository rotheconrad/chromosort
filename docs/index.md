---
title: Documentation
description: User guide for ChromoSort.
---

# ChromoSort Documentation

ChromoSort is a reference-guided assembly curation toolkit for sorting contigs, splitting reviewed chimeric contigs, making manual cuts, reviewing dot plots, scaffolding ordered contigs, applying reviewed graph-supported gap fills, and plotting existing whole-genome alignments.

This documentation is organized as a user guide rather than one very long README. The archived longform README remains available for historical reference.


## Start Here

- [Installation]({{ '/installation/' | relative_url }})
- [Input files]({{ '/input-files/' | relative_url }})
- [Workflows]({{ '/workflows/' | relative_url }})
- [Command reference]({{ '/commands/' | relative_url }})

## Commands

- [chromo sort]({{ '/commands/sort/' | relative_url }})
- [chromo fix]({{ '/commands/fix/' | relative_url }})
- [chromo cut]({{ '/commands/cut/' | relative_url }})
- [chromo manual]({{ '/commands/manual/' | relative_url }})
- [chromo plot]({{ '/commands/plot/' | relative_url }})
- [chromo scaffold]({{ '/commands/scaffold/' | relative_url }})
- [chromo gapfill]({{ '/commands/gapfill/' | relative_url }})

## Running Analyses

- [Output files]({{ '/outputs/' | relative_url }})
- [Architecture]({{ '/architecture/' | relative_url }})
- [Troubleshooting]({{ '/troubleshooting/' | relative_url }})

## Project Information

- [Current status and version history]({{ '/status/' | relative_url }})
- [Contributing]({{ '/contributing/' | relative_url }})
- [Archived longform README]({{ '/archive/longform-readme/' | relative_url }})

## Main Workflow

1. Generate MUMmer `show-coords` or minimap2 PAF alignments between the reference and assembly.
2. Use `chromo plot` and, when helpful, `chromo manual` to inspect suspicious contigs.
3. Use `chromo fix` or `chromo cut` for reviewed sequence edits.
4. Re-align the fixed assembly.
5. Use `chromo sort` to assign, filter, orient, and order contigs.
6. Use `chromo scaffold` to build one scaffold per reference sequence.
7. Use `chromo gapfill` only when reviewed GFA paths should replace N gaps.
