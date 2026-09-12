# Conda Recipe Staging

This directory tracks a draft Bioconda recipe for ChromoSort so the package is
ready for future `conda install chromosort` distribution.

The recipe is not submitted to Bioconda yet. When ChromoSort is ready for
submission:

1. Cut and push the intended release tag.
2. Download the GitHub source archive for that tag and update `sha256`.
3. Copy `meta.yaml` into `bioconda-recipes/recipes/chromosort/meta.yaml`.
4. Run the Bioconda lint/build checks and open a Bioconda pull request.


## Publication candidate

For 0.4.0rc1 validation, run `python3 scripts/prepare_release.py` after building the sdist and wheel. It generates `build/release-candidate/conda-recipe/meta.yaml` with the exact local source archive and checksum. This directory's `meta.yaml` is the historical 0.3.0 source template, not a published candidate recipe.
