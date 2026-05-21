# Contributing to ChromoSort

Thank you for your interest in contributing. This document covers how to report issues, set up a development environment, run tests, update docs, and submit changes.

## Reporting Bugs and Requesting Features

Please use GitHub Issues: <https://github.com/rotheconrad/chromosort/issues>.

When reporting a bug, include:

- The exact command you ran, with private sample names redacted if needed.
- The full error message or traceback.
- Your operating system and Python version (`python --version`).
- The ChromoSort version (`python -c "import chromosort; print(chromosort.__version__)"`).
- A small FASTA, coords, PAF, GFA, GAF, or TSV fixture if possible.

When requesting a feature, describe:

- The assembly-curation problem you are trying to solve.
- Which inputs you have and which output you expect.
- Whether the request affects sorting, splitting, plotting, scaffolding, manual review, or graph gap filling.

## Development Setup

Recommended mamba setup:

```bash
git clone https://github.com/rotheconrad/chromosort.git
cd chromosort
mamba env create -f environment.yml
mamba activate chromosort
```

Pixi setup:

```bash
git clone https://github.com/rotheconrad/chromosort.git
cd chromosort
pixi install
pixi run help
```

Plain Python editable install for lightweight code work:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e . pytest Pillow
```

## Running Tests

```bash
pytest
python -m unittest discover -s tests -v
pixi run test
```

The tests use tiny synthetic fixtures under `tests/data` and do not require running MUMmer or minimap2. New behavior should include focused tests or fixtures when practical.

## Documentation Changes

The public documentation is built from `docs/` with Jekyll and GitHub Pages. The archived longform README is kept at `docs/archive/longform-readme.md` as a source reference, while command docs live under `docs/commands/`.

For docs-only changes, update the relevant Markdown page and check links to local files and assets. If command behavior changes, update both the command page and any affected output, workflow, troubleshooting, or status page.

## Pull Request Workflow

1. Fork the repository or create a branch from `main`.
2. Make a focused change.
3. Add or update tests for changed behavior.
4. Run the relevant tests before opening the pull request.
5. Write a clear PR description explaining what changed and why.

For large behavior changes, open an issue first so the design can be discussed before implementation.

## Code of Conduct

Be respectful and constructive. ChromoSort is affiliated with scientific software development; interactions should meet the standards expected in a professional research environment. See [`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md).
