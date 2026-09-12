---
title: Publication release preparation
---

# Publication release preparation

The new interfaces are staged as **0.4.0rc2**. This candidate identifier keeps
development validation distinct from the later stable public release. No
registry upload, Git tag, DOI deposit or preprint submission is performed by
the preparation commands.

```bash
make agent-check DIST_DIR=build/package-check/dist
python3 scripts/prepare_release.py --dist-dir build/package-check/dist \
  --output-dir build/package-check/release
# Run in an environment with conda-build installed:
conda build build/package-check/release/conda-recipe
```

`build/package-check/dist/` contains this check's sdist and wheel. The explicit
`DIST_DIR` preserves any frozen archives already in `dist/`. The preparation script checks metadata
versions and records exact artifact SHA-256 values in
`build/package-check/release/release.json`. It generates a concrete conda recipe
against the local sdist, with a matching checksum. `conda-recipe/meta.yaml`
remains the historical 0.3.0 source recipe used as a template; its old source
checksum is never relabeled as the candidate's checksum.

Bioconda is the intended staged channel for this genomics CLI. The pure Python
wheel needs Pillow; BAM access additionally needs samtools, and the optional
alignment executor needs minimap2. The conda environment supplies these tools.
The Pixi lock covers the existing three platforms; routine CI keeps its
existing Ubuntu Python matrix and does not run the public biological benchmark
or upload new build artifacts.

Before stable publication, freeze an immutable software snapshot, public
simulation input versions, benchmark/comparator settings and claims. Validate
the public tutorial from fresh installs; archive the matching FASTA, AGP,
component ledgers, decision tables, figure settings, commands, software versions
and resource logs. Verify metadata/authorship and add the minted DOI to citation
metadata only after deposit. Use durable official release assets for the
approved version; do not publish working candidates as stable releases.

Read figure geometry and deterministic replay are implementation checks.
Comparative accuracy, review effort and biological error recovery require the
independent frozen benchmark. RagTag parameter and alignment-representation
differences must be reported, not treated as universal tool failures.

The frozen rc1 and rc2 source/wheel/conda artifacts retain their original
identities. A newer Git snapshot with the same candidate version is a separate
artifact: identify it by commit and archive checksum, and do not replace a
frozen candidate with a rebuild. Public-readiness changes validate read target
names and lengths against manifest assemblies and improve documentation and
package boundaries; they do not change valid-input correction policy. The rc2 refinement is described in [review policy]({{ '/review-refinement/' | relative_url }}); each candidate has its own immutable archive identity.
