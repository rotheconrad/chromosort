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
- ChromoSort command-line tools

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
`environment.yml`, including Python, minimap2, MUMmer 4, Pillow, pytest, and
the ChromoSort command-line tools.

## Why Pillow Is Included

`chromo plot` writes PDF and SVG plots directly. PNG output is different because
it is a raster image format, so ChromoSort uses Pillow when `--formats png` is
requested. This keeps PNG export inside Python and avoids requiring separate
system renderers such as Cairo, ImageMagick, browser engines, or SVG-to-PNG
conversion tools. Pillow is not used by the sorting, fixing, scaffolding, or
gap-filling algorithms.

## Verify the Install

```bash
chromo --help
chromo sort --help
chromo clean --help
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
