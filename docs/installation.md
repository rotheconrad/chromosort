---
title: Installation
description: Install ChromoSort with mamba or Pixi and verify the command-line tools.
---

# Installation

## Install Format Policy

ChromoSort intentionally supports more than one install style. Use mamba or
conda when you want the familiar scientific Python plus bioconda stack. Use
Pixi when you prefer Pixi's project task runner and lockable environment
workflow. Use a plain Python editable install for lightweight code work in
generic automation harnesses.

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
pixi run agent-check
```

The Pixi environment uses the same conda-forge/bioconda stack as
`environment.yml`, including Python, minimap2, MUMmer 4, Pillow, pytest, and
the ChromoSort command-line tools.

When `pixi.toml` changes, the matching `pixi.lock` should be refreshed and
committed so Pixi users get a reproducible environment.

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
chromo eval all --help
chromo fix --help
chromo cut --help
chromo manual --help
chromo gafprep --help
chromo plot --help
chromo scaffold --help
chromo gapfill --help
```

For development checks, run:

```bash
make agent-check
# or, from a Pixi environment
pixi run agent-check
# or the underlying test suite
python -m unittest discover -s tests -v
# pytest is also fine when available
pytest
# legacy Pixi test task
pixi run test
```
