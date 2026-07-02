---
title: Contributing
description: How to contribute issues, docs, tests, fixtures, and code to ChromoSort.
---

# Contributing

Contributions are welcome through GitHub Issues and pull requests. For the repository-level contributor guide, see [`CONTRIBUTING.md`](https://github.com/rotheconrad/chromosort/blob/main/CONTRIBUTING.md).

## Helpful Contributions

- Real-world regression fixtures for sorting, fixing, scaffolding, and graph-aware gap filling.
- Documentation improvements and worked examples for new users.
- Reports from assemblies where MUMmer and minimap2 produce meaningfully different placement evidence.
- Edge cases for GFA/GAF/Hi-C-like support and ambiguous graph paths.

## Development Checks

```bash
mamba env create -f environment.yml
mamba activate chromosort
make agent-check
```

or:

```bash
pixi install
pixi run agent-check
```

For a focused code check, `python -m unittest discover -s tests -v` and
`pytest` both exercise the synthetic fixture suite.
