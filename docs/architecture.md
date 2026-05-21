---
title: Architecture
description: How ChromoSort commands separate evidence, review, and sequence-changing steps.
---

# Architecture

ChromoSort is organized around one rule: sequence-changing steps should be explicit and auditable. Alignment evidence, graph evidence, dot plots, and manual review artifacts are kept as reports until a command is asked to write a new FASTA.

## Command Responsibilities

| Stage | Command | Responsibility |
| --- | --- | --- |
| Alignment review | `chromo plot`, `chromo manual` | Draw or review alignments without changing the source FASTA. |
| Placement and filtering | `chromo sort` | Assign contigs to reference sequences, classify overlaps, and write a reference-ordered FASTA. |
| Reviewed sequence edits | `chromo fix`, `chromo cut`, `chromo manual apply` | Split or cut contigs after explicit review or configured planning. |
| Scaffold construction | `chromo scaffold` | Join final ordered contigs with inferred or fixed gaps, with optional explicit overlap trimming. |
| Graph-supported gap fill | `chromo gapfill` | Plan graph paths first, then apply only fillable and accepted sequence paths. |

## Evidence Layers

ChromoSort accepts standard whole-genome alignment formats for placement and plotting:

- MUMmer `show-coords`
- minimap2 PAF

Graph-aware review and gap filling can add:

- GFA segment/link records
- GAF read-to-graph alignments
- Hi-C-like graph-node contact tables
- reference-placement PAF for graph intermediate nodes

Graph evidence is report-only in `sort`, `fix`, `manual`, and `scaffold` unless a command-specific explicit policy says otherwise. Graph sequence insertion is limited to `chromo gapfill --apply`.

## Data Boundaries

The command modules in `src/chromosort/` keep parsing, planning, and writing local to each workflow. Shared concepts such as FASTA reading, alignment segment normalization, graph parsing, and output reports are reused where they reduce duplication, but the command boundary remains visible to users through separate reports and subcommand documentation.
