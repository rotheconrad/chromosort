# Copilot Instructions

Use `AGENTS.md` as the canonical repository instructions. This file is only an
adapter for GitHub Copilot.

Important defaults:

- Work from source files, not generated `_site/` or `docs/_site/` output.
- Preserve the FASTA/alignment compatibility rule.
- Add or update focused tests for behavior changes.
- Update command docs when public CLI behavior or output schemas change.
- Prefer `make agent-check` for local verification.
- Pixi users can run `pixi run agent-check`.
