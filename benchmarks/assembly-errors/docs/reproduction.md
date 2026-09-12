# Reproduce the public benchmark

Commands below run from `benchmarks/assembly-errors/` in a repository clone.
No unpublished file or historical candidate archive is needed to verify the
fixtures, regenerate all 14 cases, summarize historical score records, or
evaluate an explicitly identified software installation.

## Dependencies and identities

| Work | Dependencies |
| --- | --- |
| Fixture/truth checks, retrieval from fixtures, generation, scoring, score summaries, tests | Python 3.9–3.12, standard library |
| Exact recorded mapping reproduction | minimap2 2.31-r1302 available on PATH |
| Optional BAM generation | samtools 1.23.1 in addition to minimap2 |
| ChromoSort runs | Identified installation; package-supported Python and Pillow dependencies |
| RagTag baseline runs | RagTag v2.1.0 at commit e993fb2787346345662c240967a3ec80bb2aac64, its Python dependencies and minimap2 |

Original generation used Python 3.9.6/zlib 1.2.12; rc1/rc2 used Python 3.12.14
and Pillow 12.3.0. Original RagTag ran with Python 3.13.12, intervaltree 3.1.0,
networkx 3.4.2, numpy 2.4.4 and pysam 0.23.3. These are observed environments,
not a solved cross-platform lock. Python/zlib differences can change compressed
FASTQ bytes; [expected identities](../sources/development_expected.json) compare
decoded FASTQ and exact FASTA/PAF/provenance. The minimap2 binary digest in
[software identities](../sources/software_identities.json) describes the original
platform binary, not a requirement that all operating systems use identical
executable bytes.

## Offline quick check

```sh
python3 scripts/verify.py --fixtures
python3 -m unittest discover -s tests -v
python3 scripts/summarize_history.py
python3 scripts/run_benchmark.py --tool unedited --output-root work/unedited
```

Expected: three fixture bundles, 1,853 PAF rows, exact sequence-to-truth
provenance, ten analytic tests and 294 validated historical development score
identities. The unedited command scores all three included examples without
installing ChromoSort or invoking any mapper.

For an isolated-copy check that excludes ignored generated data and caches:

```sh
python3 scripts/check_clean_checkout.py
python3 scripts/check_clean_checkout.py --with-mapping
```

The second command regenerates/remaps all 14 cases in a temporary copy using
only the selected public files and installed minimap2. It performs no network
request or moving-ChromoSort benchmark.

## Full development regeneration

```sh
python3 scripts/retrieve.py --output-root work
python3 scripts/generate.py --output-root work
python3 scripts/evidence.py work/development/D* --threads 2
python3 scripts/verify.py --regenerated work
```

The expected full result is 14 cases and 6,254 PAF rows. Inputs, target FASTA,
input-to-target maps, read origins and decoded FASTQ are checked against the
original development identities. Only the explanatory D001 construction prose
was corrected; see [the source erratum](sources.md). The public generator
rejects non-development case IDs. It makes no held-out dataset.

The three compact fixtures supply every canonical public source interval, so
the retrieval step is offline by default. To independently audit NCBI retrieval
instead, use a **fresh output root**, explicitly request downloads and supply
your contact email for E-utilities:

```sh
python3 scripts/retrieve.py --output-root work/ncbi-audit \
  --download --email YOUR_CONTACT_EMAIL
```

Five versioned accession URLs are fetched with rate limiting, normalized and
checked before their canonical files are installed. The script does not save
the email. Existing cache files are verified and reused; use a fresh root when
you want a network audit. Network availability is not needed for ordinary use.

Generate a subset with `--cases D003 D010`. Optional `scripts/evidence.py ...
--bam` adds BAM/index files and records commands using paths relative to each
generated case. The portable tracked examples retain PAF only; raw reads are
available after generation, not as dangling fixture paths.

## Evaluate a fixed ChromoSort installation

After choosing a release or clean source commit, install it in an isolated
environment. For a repository checkout, an example setup is:

```sh
python3.12 -m venv work/tool-env
work/tool-env/bin/python -m pip install ../..
```

Record the selected commit/archive SHA, inspect the checkout state, and replace
`COMMIT_OR_ARCHIVE_SHA256` below with that identity. The runner records the
installed version and a digest of the installed Python source files in addition
to the supplied identity. It preserves virtual-environment interpreter paths
and does not import the repository's moving source through PYTHONPATH.

```sh
python3 scripts/run_benchmark.py --tool chromosort \
  --python work/tool-env/bin/python --software-id COMMIT_OR_ARCHIVE_SHA256 \
  --output-root work/chromosort-fixtures
```

The default pairs `review` and `report` read-continuity policies in conservative,
comprehensive and sensitive modes across the three included fixtures. To use
all 14 regenerated cases:

```sh
python3 scripts/run_benchmark.py --tool chromosort \
  --python work/tool-env/bin/python --software-id COMMIT_OR_ARCHIVE_SHA256 \
  --case-root work/development --truth-root work/evaluator/development \
  --output-root work/chromosort-all14
```

The saved rc2 protocol selects only `reads_sim01`, MAPQ ≥20, 1 kb continuous
CIGAR anchors and two distinct spanning molecules. There is no read pooling,
BAM fallback or orientation normalization. `--modes` and `--policies` select
explicit subsets for a predeclared run; choosing settings after viewing scores
must be labeled development tuning, not independent validation.

For review-only inspection, run the installed command directly, for example:

```sh
work/tool-env/bin/chromo workflow scan \
  --manifest fixtures/development/D003/manifest.json --mode comprehensive \
  --read-evidence-id reads_sim01 --inspect-min-gap-bp 1000 \
  --output-dir work/D003-review
```

Pending decisions and inspection panels are not scored as corrected output.
Any later changed assembly must be freshly aligned before the next stage.

## RagTag comparison

Install the pinned source in a separate environment; do not copy its source or
environment into tracked benchmark content:

```sh
python3 -m venv work/ragtag-env
work/ragtag-env/bin/python -m pip install \
  'RagTag @ https://github.com/malonge/RagTag/archive/e993fb2787346345662c240967a3ec80bb2aac64.tar.gz'
```

RagTag's package declares intervaltree, numpy, pysam and networkx dependencies.
Use the recorded compatible versions when replaying the historical environment;
record any replacement environment as a new execution. Keep minimap2 on PATH.

```sh
python3 scripts/run_benchmark.py --tool ragtag \
  --ragtag-bin work/ragtag-env/bin \
  --software-id e993fb2787346345662c240967a3ec80bb2aac64 \
  --case-root work/development --truth-root work/evaluator/development \
  --output-root work/ragtag-all14
```

This runs five correction and two separate scaffold arms at two threads. For
a fast fixture-only run without raw FASTQ, select `--arms ragtag_default`.
Read-validation arms explicitly require regenerated FASTQ. Combined reference
FASTA is created in the result directory, preserving reference namespaces.

The runner calls official `ragtag_correct.py` / `ragtag_scaffold.py` directly
because the v2.1.0 dispatcher can suppress child failures. It does not patch
RagTag. Read [comparison fairness](results.md) before interpreting the arms.

## Output provenance and historical replay limits

Every new run writes software/input identity and parameter records, executes
without evaluator truth, commits output/run hashes, and only then opens the
independent scorer. Exact FASTA/AGP reconstruction is required. Logs and bulk
outputs remain in the chosen work root; original recorded scores are unchanged.

The repository contains the original 294 development score records and their
checksums, not all original bulk output FASTAs or software archives. rc1/rc2
archive SHA values identify historical snapshots, but no public download URL
is asserted for those archives. Exact historic tool replay requires the
matching snapshot if available; a newly installed repository commit produces
new results with its own identity. Data regeneration and historical summary
recalculation have no dependency on those snapshots or private artifacts.
