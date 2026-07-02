# Claude Notes

Use `AGENTS.md` as the canonical repository instructions. This file is only an
adapter for Claude-based tools.

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

Do not work from generated output directories such as `_site/`, `docs/_site/`,
`temp/`, `data/`, or `results/`.
