## Summary

-

## Change Type

- [ ] Code behavior
- [ ] CLI or output schema
- [ ] Documentation
- [ ] Tests or fixtures
- [ ] Packaging, CI, or agent-readiness infrastructure

## Checks Run

- [ ] `make agent-check`
- [ ] `pixi run agent-check`
- [ ] Focused test:
- [ ] Docs check:
- [ ] Not run:

## Docs And Fixtures

- [ ] Updated relevant command docs under `docs/commands/`
- [ ] Updated README, status, roadmap, or output docs if needed
- [ ] Added or updated synthetic fixtures under `tests/data/` if needed
- [ ] Refreshed `pixi.lock` if `pixi.toml` changed
- [ ] No generated `_site/`, `docs/_site/`, `temp/`, `data/`, or `results/` files committed

## Evidence Boundary

- [ ] This change preserves the FASTA/alignment compatibility rule
- [ ] This change does not make silent dataset-specific biological decisions
- [ ] Ambiguous biological evidence remains explicit in reports, review tables, dashboards, or docs
