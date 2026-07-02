# Agent Notes For ChromoSort

These notes are for coding agents, assistant chats, and future project threads
working in this repository. They summarize the repo layout, safe verification
commands, and the assembly-curation judgment rules that are easy to lose when
moving between chats or datasets.

## Repository Orientation

- Public CLI entry point: `chromo`, dispatched from `src/chromosort/cli.py`.
- Core modules:
  - `reference_order.py`: `chromo sort` assignment, filtering, reports.
  - `fix_contigs.py`: `chromo fix` breakpoint planner and reviewed-plan apply.
  - `evaluate.py`: `chromo eval fix/scaffold/gapfill` review TSVs.
  - `plot.py`: `chromo plot` dot plots from existing coords or PAF.
  - `scaffold.py` and `gapfill.py`: ordered-contig scaffolding and graph fills.
  - `graph.py`, `gaf.py`, `longreads.py`: GFA, GAF, and read-PAF evidence.
- User docs live in `docs/`; command pages live in `docs/commands/`.
- Synthetic fixtures live in `tests/data/`. Tests do not require MUMmer or
  minimap2.
- Agent-readiness planning lives in `docs/agent-readiness-roadmap.md`.
- Change-type recipes live in `docs/agent-task-recipes.md`.
- Ignore generated or local-output directories unless a task explicitly targets
  them: `_site/`, `docs/_site/`, `.pytest_cache/`, `.venv/`, `vendor/`,
  `temp/`, `data/`, and `results/`.

## Install And Command Surface

Mamba/conda and Pixi are both first-class development paths. Use
`environment.yml` for the mamba/conda stack and `pixi.toml` for Pixi users.
Plain Python editable installs are also valid for lightweight code work:

```bash
python -m pip install -e . pytest Pillow
```

The shared harness-neutral checks are:

```bash
make agent-check
make smoke
make test
```

Pixi users can run the same contract through Pixi:

```bash
pixi run agent-check
pixi run smoke
pixi run test
```

When `pixi.toml` changes, refresh and commit `pixi.lock` if Pixi is available.

## Development Checks

Use the smallest relevant check while editing, then broaden when behavior
changes:

```bash
make agent-check
python3 -m unittest discover -s tests
python3 tests/test_reference_order.py
python3 tests/test_fix_contigs.py
python3 tests/test_eval.py
python3 tests/test_plot.py
python3 tests/test_longreads.py
git diff --check
```

If the environment has `pytest`, this is also valid:

```bash
pytest
```

For docs-only changes, at minimum run:

```bash
git diff --check
```

When adding a public docs page, update `docs/_config.yml`, `docs/index.md`,
and any README/status/changelog pointers that should expose the page.

## Evidence Discipline

Every MUMmer coords or minimap2 PAF file describes one exact reference FASTA
and one exact assembly FASTA. If a ChromoSort command writes a changed FASTA by
filtering, splitting, cutting, reverse-complementing, renaming, manually
exporting, scaffolding, or gap-filling records, re-align that changed FASTA
before using it as the assembly input to another alignment-dependent command.

`chromo plot --assignments` is a review convenience: it plots original
alignment rows and orders the query axis by a `chromo sort` assignment report.
It is not a fresh alignment of `ordered.fa`, `fixed.fa`, or a scaffold FASTA.

## Recommended PAF Recipe

For reference-vs-assembly PAF, keep `-c --secondary=no`. PAF is the recommended
primary alignment input for most new ChromoSort runs because it is fast and
supports MAPQ filtering:

```bash
minimap2 -x asm5 -c -t 16 --secondary=no reference.fa assembly.fa > sample.paf
```

Use `asm10` or `asm20` only when stricter presets miss expected chromosome-scale
alignments. With permissive presets, review dot plots and consider `--min-mapq`.
Avoid setting `--min-segment-idy` on PAF until the PAF identity distribution is
checked, because ChromoSort computes PAF identity from the match and block
length columns.

Recent `chromo sort --paf` run summaries include PAF diagnostics such as row
counts, `tp:A` classes, and `cg:Z` presence. If `cg:Z` is absent, check whether
minimap2 was run without `-c`.

## Review Decision Rules

Choose one primary whole-genome alignment source for ordinary ChromoSort runs:
MUMmer coords or minimap2 PAF. These are alternative representations of
reference-to-assembly alignment evidence, not independent biological
validation of a candidate event. Running both can be useful for benchmarking,
parser checks, or aligner-parameter tuning; in soybean fix testing, split
counts differed by about 5-10% and marginal split-contig sets by about 20-30%.
Treat those as alignment/output differences unless plots or tests indicate a
parser bug. Read-to-assembly PAF, GFA, GAF, and independent assemblies are
stronger support for real biological decisions.

- Multi-reference `kept_split_candidate` calls in the primary alignment are
  fix-review targets. Run targeted `chromo fix --contigs ...`, usually with
  `--mode conservative`, not `--all`, after plot and evidence review.
- PAF-only or coords-only candidates from an optional comparison are review
  targets, not automatic production fixes. Inspect per-reference plots and,
  when available, graph and long-read evidence with `chromo eval fix`.
- A whole contig aligning to one reference in the opposite orientation is
  usually an orientation issue. Use `chromo sort --orient-to-reference`; do not
  split it.
- A simple same-reference internal inversion is not automatically an assembly
  error. For pangenome graph inputs, preserve a real inversion as haplotype
  structure instead of reference-normalizing it. Use `chromo eval fix --mode
  comprehensive` plus read and graph evidence to decide whether it is real.
  Comprehensive mode is orientation-aware after smoothing and may differ from
  conservative mode; use it as review evidence, not as a guaranteed superset.
- `chromo fix` does not resolve cross-contig terminal overlaps. Use sort and
  scaffold overlap reports, then choose an explicit scaffold overlap policy
  only after review.

For pangenome graph preparation, especially reference-based minigraph-cactus
workflows, keep the distinction between biological structure and
reference-normalized pseudoassemblies clear. Whole-contig orientation to the
reference is often useful. Reorienting a validated internal inversion erases
the allele from the assembly path and should be reserved for intentional
reference-normalized products or confirmed assembly errors.

## Reusable Review Commands

Targeted conservative fix for reviewed multi-reference split candidates:

```bash
chromo fix \
  --assembly-fasta assembly.fa \
  --paf sample.ref_vs_asm.paf \
  --contigs contig_a contig_b \
  --mode conservative \
  --min-mapq 20 \
  --orient-to-reference \
  --output-fasta results/sample.fixed.fa \
  --report results/sample.fix_report.tsv
```

Same-reference inversion review without committing to a fix:

```bash
chromo eval fix \
  --assembly-fasta assembly.fa \
  --paf sample.ref_vs_asm.paf \
  --contigs contig_with_inversion \
  --mode comprehensive \
  --min-mapq 20 \
  --gfa assembly_graph.gfa \
  --gaf reads_to_graph.gaf \
  --read-paf reads_to_assembly.paf \
  --read-window-bp 10000 \
  --read-min-anchor-bp 2000 \
  --orient-to-reference \
  --output-prefix results/sample.inversion_eval
```

Long reads mapped back to the assembly for breakpoint and bridge evidence:

```bash
minimap2 -x map-hifi -c -t 16 --secondary=no assembly.fa reads.fastq.gz \
  > reads_to_assembly.paf
```

GAF read-to-graph evidence can come from GraphAligner or a compatible graph
mapper. GFA segment names should match assembly contig names for `eval fix`
node context, and graph-path names for scaffold/gapfill evidence.

## New Dataset Handoff Checklist

Before asking another agent or chat to continue an analysis, provide or create:

- Reference FASTA path and assembly FASTA path for each sample.
- The exact MUMmer and/or minimap2 commands used, including preset and filters.
- Output folders for the primary sort and plots, plus any optional comparison
  run.
- `*.run_summary.txt`, `*.contig_assignments.tsv`, and
  `*.contig_ref_matches.tsv` for every sample under review.
- Focused plot PDFs or PNGs for any contig-level decision.
- GFA, read-to-assembly PAF, and read-to-graph GAF when deciding whether a
  same-reference inversion or breakpoint is real.
- A short decision table: fix now, review before fix, leave native, or
  reference-normalized experimental output.

Prefer general statements over crop-specific assumptions. Reference names such
as `Gm18` are examples; the same logic applies to other species, assemblies,
and chromosome naming schemes.
