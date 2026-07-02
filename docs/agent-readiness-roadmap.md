---
title: Agent Readiness Roadmap
description: Roadmap for making ChromoSort easy and safe for coding agents and open-source automation harnesses.
---

# Agent Readiness Roadmap

This roadmap tracks the work needed to make ChromoSort straightforward for
coding agents, assistant chats, and open-source automation harnesses to inspect,
modify, test, and hand off safely. The target audience includes Codex, Claude,
Gemini, Jewel/Jules, Copilot, Aider, OpenHands, and any model-plus-harness
workflow that can read repository files and run commands.

The goal is not to let agents make biological calls automatically. The goal is
to give every agent the same repository map, install choices, check commands,
fixtures, and evidence boundaries so implementation work is reproducible.

## Design Stance

ChromoSort should support multiple install and execution styles:

- **Mamba/conda first-class support.** `environment.yml` remains the preferred
  scientific-stack setup path for users who want conda-forge and bioconda tools
  such as minimap2 and MUMmer in one environment.
- **Pixi first-class support.** `pixi.toml` remains a parallel lockable workflow
  for collaborators who prefer Pixi task execution and reproducible environments.
- **Plain Python fallback.** Lightweight code and test work should remain
  possible with `python -m pip install -e . pytest Pillow`, because generic
  harnesses often start from a Python image rather than a full conda stack.

Agent-facing commands should avoid requiring external aligners for routine test
work. The synthetic fixtures under `tests/data/` are the default local contract.

## Readiness Phases

| Phase | Status | Outcome |
| --- | --- | --- |
| 1. Save canonical roadmap and expose it in docs. | Done | This page records the plan and links from the public docs. |
| 2. Add a harness-neutral command surface. | Done | `Makefile` and Pixi tasks expose `test`, `smoke`, `diff-check`, and `agent-check`. |
| 3. Add agent adapter files. | Done | Common agent entry files point back to `AGENTS.md` as the canonical instruction source. |
| 4. Add code CI beyond documentation deployment. | Done | GitHub Actions runs unit tests, CLI smoke checks, and whitespace checks on code changes. |
| 5. Add contribution templates for agent-ready changes. | Done | Issue and PR templates ask for fixtures, checks, docs, and evidence-boundary notes. |
| 6. Add richer smoke or golden workflow fixtures. | Done | `tests/test_agent_smoke_workflow.py` exercises sort, plot, fix, eval-all, scaffold, and gapfill apply on synthetic fixtures. |
| 7. Add conservative linting. | Done | Ruff is configured with a narrow high-signal rule set and wired into Make, Pixi, and CI. |
| 8. Add optional packaging checks. | Done | `make package-check` and `pixi run package-check` build sdist/wheel without build isolation. |
| 9. Decide Pixi lockfile policy. | Done | Commit `pixi.lock` when generated or updated so Pixi users get a reproducible first-class environment. |

## Canonical Agent Contract

The canonical instruction file is `AGENTS.md`. Adapter files for individual
tools should stay short and should not fork the rules. If behavior changes,
update `AGENTS.md` first, then refresh adapters only when needed.

Agents should:

- read `AGENTS.md`, `README.md`, `docs/architecture.md`, and the command page
  matching the area being edited;
- work from source files, not generated documentation output;
- ignore `_site/`, `docs/_site/`, `.pytest_cache/`, `.venv/`, `vendor/`,
  `temp/`, `data/`, and `results/`;
- use synthetic fixtures under `tests/data/` for routine code checks;
- preserve the FASTA/alignment compatibility rule;
- add or update focused tests for behavior changes;
- update command docs, output docs, and changelog/status pointers when public
  behavior changes.

## Shared Checks

Routine local verification:

```bash
make agent-check
```

Equivalent pieces:

```bash
make diff-check
make lint
make smoke
make test
make package-check
```

Pixi equivalents:

```bash
pixi run agent-check
pixi run smoke
pixi run test
```

Plain Python equivalent:

```bash
git diff --check
ruff check .
PYTHONPATH=src python -m chromosort.cli --help
python -m unittest discover -s tests -v
python -m build --sdist --wheel --no-isolation
```

Docs-only changes should at least run `git diff --check`. Public documentation
page additions should also update `docs/_config.yml`, `docs/index.md`, and any
README/status/changelog pointers that should expose the new page.

## Agent Task Recipes To Add

Agent-facing recipes live in
[Agent Task Recipes]({{ '/agent-task-recipes/' | relative_url }}). They cover:

- adding or changing a CLI flag;
- adding a TSV output column;
- changing a FASTA-writing command;
- changing the review-event schema;
- changing plotting behavior;
- adding a graph, GAF, or long-read evidence field;
- making docs-only changes.

Each recipe should name the likely source files, docs, fixtures, and tests to
touch, plus the smallest useful verification command.

## Pixi Lockfile Policy

Commit `pixi.lock` whenever `pixi.toml` changes and Pixi can solve the
environment. The lockfile makes the Pixi path as concrete as the mamba/conda
path while still keeping `environment.yml` as a readable conda-compatible
environment specification.

If a contributor cannot run Pixi locally, they should update `pixi.toml`, run
the Make/Python checks they can run, and note that the lockfile needs refresh in
the PR checklist.

## Biological Decision Boundaries

Agents may change code, tests, docs, synthetic fixtures, CI, and packaging
metadata. They should not silently make dataset-specific biological decisions.
Any sequence-changing recommendation for real data should identify the exact
FASTA, alignment, graph/read evidence, review table, and command path involved.

When evidence is ambiguous, ChromoSort should keep the ambiguity visible through
reports, review tables, dashboards, or explicit user decisions rather than
guessing.
