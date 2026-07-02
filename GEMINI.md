# Gemini Notes

Use `AGENTS.md` as the canonical repository instructions. This file is only an
adapter for Gemini-based tools.

Before editing, read:

- `AGENTS.md`
- `README.md`
- `docs/agent-readiness-roadmap.md`
- the command docs under `docs/commands/` that match the requested change

Prefer the shared checks:

```bash
make agent-check
```

For Pixi environments:

```bash
pixi run agent-check
```

Do not make dataset-specific biological decisions silently. Preserve ambiguity
in reports, review tables, docs, or explicit user-facing notes.
