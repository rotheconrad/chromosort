---
title: Installation
description: Install ChromoSort with mamba or Pixi and verify the command-line tools.
---

# Installation

## Installation With Mamba

Install Mambaforge, Miniforge, or another conda-compatible distribution with
`mamba`, then create the environment:

```bash
git clone https://github.com/rotheconrad/chromosort.git
cd chromosort

mamba env create -f environment.yml
mamba activate chromosort
```

The environment installs:

- Python
- minimap2
- MUMmer 4 (`nucmer`, `delta-filter`, `show-coords`)
- Pillow for optional PNG plot output
- pytest for the test suite
- ChromoSort in editable mode

Additional console entry points are retained for compatibility:

- `chromosort` is equivalent to `chromo sort`
- `chromosort-fix-contigs` is equivalent to `chromo fix`
- `chromosort-scaffold` is equivalent to `chromo scaffold`
- `chromosort-gapfill` is equivalent to `chromo gapfill`

New workflows should use `chromo sort`, `chromo fix`, `chromo cut`,
`chromo manual`, `chromo scaffold`, `chromo gapfill`, and `chromo plot`.

## Installation With Pixi

If you prefer [Pixi](https://pixi.sh), ChromoSort also ships a `pixi.toml`
environment file:

```bash
git clone https://github.com/rotheconrad/chromosort.git
cd chromosort

pixi install
pixi run help
pixi run test
```

The Pixi environment uses the same conda-forge/bioconda stack as
`environment.yml`, including Python, minimap2, MUMmer 4, Pillow, pytest, and an
editable install of the local ChromoSort package.

## Verify the Install

```bash
chromo --help
chromo sort --help
chromo fix --help
chromo cut --help
chromo manual --help
chromo plot --help
chromo scaffold --help
chromo gapfill --help
```

For development checks, run:

```bash
pytest
# or
pixi run test
```
