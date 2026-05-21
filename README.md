# ChromoSort

Reference-guided genome assembly utilities for sorting contigs, splitting reviewed chimeric contigs, cutting exact breakpoints, manual dot-plot review, plotting alignments, scaffolding final ordered contigs, and applying reviewed graph-supported gap fills.

ChromoSort provides one command, `chromo`, with seven subcommands:

| Command | Purpose |
| --- | --- |
| `chromo sort` | Assign assembly contigs to reference sequences, filter low-value overlaps, and write a reference-ordered FASTA. |
| `chromo fix` | Split selected contigs, or all detected candidates, into reference-labeled pieces. |
| `chromo cut` | Cut contigs at exact reviewed positions. |
| `chromo manual` | Write a self-contained browser dashboard for dot-plot review and reproducible manual edits. |
| `chromo plot` | Draw PDF/SVG/PNG dot plots from existing MUMmer coords or minimap2 PAF alignments. |
| `chromo scaffold` | Join final sorted contigs into per-reference scaffold records with auditable gaps and optional overlap trimming. |
| `chromo gapfill` | Plan and optionally apply reviewed graph-supported fills between adjacent sorted contigs. |

## Documentation

The longform README has been split into a documentation site under [`docs/`](docs/index.md), using the same Jekyll/GitHub Pages mechanics as Panex Privus.

- Public docs: <https://rotheconrad.github.io/chromosort/>
- Local docs entry point: [`docs/index.md`](docs/index.md)
- Exact longform README snapshot: [`README.longform.md`](README.longform.md)
- Rendered longform archive: [`docs/archive/longform-readme.md`](docs/archive/longform-readme.md)
- Command reference: [`docs/commands/`](docs/commands/index.md)

## Quick Start

```bash
git clone https://github.com/rotheconrad/chromosort.git
cd chromosort

mamba env create -f environment.yml
mamba activate chromosort

chromo --help
chromo sort --help
```

Typical reviewed workflow:

```bash
# 1. Fix reviewed/suspect raw contigs.
chromo fix \
  --assembly-fasta assembly.fa \
  --coords mummer/raw.coords \
  --contigs suspect_contig_1 suspect_contig_2 \
  --output-fasta results/sample.fixed.fa \
  --report results/sample.fixed_contigs.tsv

# 2. Re-align the fixed FASTA, then sort.
chromo sort \
  --ref-fasta reference.fa \
  --assembly-fasta results/sample.fixed.fa \
  --coords mummer/fixed.coords \
  --output-prefix results/sample.fixed \
  --orient-to-reference

# 3. Plot from the same alignment file for visual review.
chromo plot \
  --ref-fasta reference.fa \
  --assembly-fasta results/sample.fixed.fa \
  --coords mummer/fixed.coords \
  --assignments results/sample.fixed.contig_assignments.tsv \
  --output-prefix plots/sample.fixed \
  --per-ref
```

## Install With Pixi

```bash
pixi install
pixi run help
pixi run test
```

## Project Files

- [`docs/installation.md`](docs/installation.md): mamba and Pixi setup.
- [`docs/input-files.md`](docs/input-files.md): MUMmer, minimap2, GFA, GAF, and Hi-C-like inputs.
- [`docs/workflows.md`](docs/workflows.md): quick start, overlap handling, and synthetic graph workflow.
- [`docs/outputs.md`](docs/outputs.md): FASTA, TSV, HTML, plot, and run-summary outputs.
- [`CONTRIBUTING.md`](CONTRIBUTING.md): issue reports, development setup, tests, and PR expectations.
- [`CITATION.cff`](CITATION.cff): citation metadata.

## Current Status

Current version: `0.2.23`. Operational commands are `sort`, `fix`, `cut`, `manual`, `plot`, `scaffold`, and `gapfill`. See [`docs/status.md`](docs/status.md) or [`CHANGELOG.md`](CHANGELOG.md) for version history.

## Citation

If you use ChromoSort, cite this repository and cite MUMmer or minimap2 for the alignment files used by the workflow. See [`CITATION.cff`](CITATION.cff).

## Contact

Please use the GitHub issue tracker for bug reports, feature requests, and questions: <https://github.com/rotheconrad/chromosort/issues>.

## License

ChromoSort is released under the MIT License. See [`LICENSE`](LICENSE).
