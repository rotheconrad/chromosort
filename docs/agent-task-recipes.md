---
title: Agent Task Recipes
description: File, test, and documentation recipes for common ChromoSort coding-agent tasks.
---

# Agent Task Recipes

These recipes help coding agents make focused ChromoSort changes without
guessing which files, tests, or docs should move together. Use them with
`AGENTS.md` and the [Agent Readiness Roadmap]({{ '/agent-readiness-roadmap/' | relative_url }}).

## Shared First Steps

Before editing:

1. Read `AGENTS.md`.
2. Read the command page under `docs/commands/` that matches the change.
3. Search the relevant source module and tests with `rg`.
4. Use synthetic fixtures under `tests/data/` unless the task explicitly
   provides real data.

After editing, run the smallest focused check, then broaden to:

```bash
make agent-check
```

Pixi equivalent:

```bash
pixi run agent-check
```

## Add Or Change A CLI Flag

Inspect:

- `src/chromosort/cli.py` for top-level dispatch changes.
- The command module that owns the flag, such as `reference_order.py`,
  `fix_contigs.py`, `evaluate.py`, `scaffold.py`, or `gapfill.py`.
- `tests/test_cli.py` for help text and parser guardrails.

Update:

- command-specific tests for behavior;
- `docs/commands/<command>.md`;
- `README.md`, `docs/workflows.md`, or `docs/input-files.md` if the flag is
  part of a recommended workflow.

Verify:

```bash
python3 tests/test_cli.py
python3 -m unittest discover -s tests -v
```

## Add A TSV Output Column

Inspect:

- the writer function in the owning module;
- downstream readers in `review.py`, `evaluate.py`, `manual.py`, `scaffold.py`,
  or `gapfill.py`;
- output documentation in `docs/outputs.md` and command docs.

Update:

- fieldname lists and row dictionaries together;
- tests that assert headers or row values;
- docs that describe the table schema.

Guardrail:

For review-event TSVs, preserve schema compatibility. Add new columns as
explicit fields or discovered extras without breaking existing accepted rows.

## Change A FASTA-Writing Command

Inspect:

- FASTA-writing command module;
- `agp.py` and `submission.py`;
- tests for AGP, component-provenance, and submission-checklist sidecars.

Update:

- FASTA output tests;
- AGP/component/checklist tests;
- `docs/outputs.md`;
- command docs and workflow docs if command semantics change.

Guardrail:

Every FASTA-changing command creates a new coordinate system. Preserve the rule
that downstream alignment-dependent commands need fresh coords or PAF for that
exact FASTA.

## Change The Review-Event Schema

Inspect:

- `src/chromosort/review.py`;
- `src/chromosort/evaluate.py`;
- reviewed-plan application paths in `fix_contigs.py`, `scaffold.py`, and
  `gapfill.py`;
- manual task dashboards in `manual.py`.

Update:

- `tests/test_review.py`;
- `tests/test_eval.py`;
- matching executor tests;
- `docs/guides/spreadsheet-review-tables.md`;
- `docs/commands/eval.md` and command docs for reviewed execution.

Guardrail:

Accepted reviewed rows must be revalidated against the current FASTA,
assignment, graph, and command settings before sequence changes are applied.

## Change Plotting Or Manual Dashboard Behavior

Inspect:

- `plot.py` for static plots;
- `manual.py` for browser dashboards;
- `docs/dot-plots.md`;
- `docs/guides/manual-dashboard-review.md`.

Update:

- `tests/test_plot.py` or `tests/test_manual.py`;
- docs assets only when a visual example changes materially.

Verify:

```bash
python3 tests/test_plot.py
python3 tests/test_manual.py
```

## Add Graph, GAF, Or Long-Read Evidence Fields

Inspect:

- `graph.py`, `gaf.py`, and `longreads.py`;
- `evaluate.py`, `gapfill.py`, and `manual.py`;
- mixed-evidence fixtures under `tests/data/mixed_evidence/`.

Update:

- focused parser or summary tests;
- `tests/test_eval.py`, `tests/test_gapfill.py`, or `tests/test_manual.py`;
- `docs/input-files.md`, `docs/architecture.md`, and relevant guides.

Guardrail:

Keep evidence authority explicit. GFA, GAF, long-read PAF, Hi-C-like contacts,
reference-placement PAF, and patch tables are not interchangeable, and most are
review or branch-support evidence rather than direct sequence-edit authority.

## Docs-Only Change

Inspect:

- `docs/_config.yml` for navigation;
- `docs/index.md`;
- `README.md`, `docs/status.md`, and `CHANGELOG.md` when public behavior or
  roadmap status changes.

Verify:

```bash
git diff --check
```

If Ruby/Bundler is available:

```bash
make docs-check
```
