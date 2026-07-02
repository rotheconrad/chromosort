# Jewel Notes

Use `AGENTS.md` as the canonical repository instructions. This file is only an
adapter for Jewel/Jules-style agent workflows.

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

Synthetic fixtures under `tests/data/` are the default test contract. Routine
agent checks should not require running MUMmer, minimap2, or GraphAligner.
