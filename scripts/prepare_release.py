#!/usr/bin/env python3
"""Prepare exact local package identities and a candidate conda recipe; never publish."""

import argparse
import hashlib
import json
import re
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--dist-dir", type=Path, help="Built sdist/wheel directory; defaults to ROOT/dist. Use a separate directory to preserve frozen candidates.")
    parser.add_argument("--output-dir", type=Path, default=Path("build/release-candidate"))
    args = parser.parse_args()
    root = args.root.resolve()
    version = re.search(r'^version = "([^"]+)"', (root / "pyproject.toml").read_text(), re.M).group(1)
    for path in [root / "CITATION.cff", root / "src/chromosort/__init__.py", root / "pixi.toml"]:
        if f'"{version}"' not in path.read_text():
            raise ValueError(f"Version mismatch: {path}")
    dist = args.dist_dir.resolve() if args.dist_dir else root / "dist"
    sdist = dist / f"chromosort-{version}.tar.gz"
    wheel = dist / f"chromosort-{version}-py3-none-any.whl"
    if not sdist.is_file() or not wheel.is_file():
        raise ValueError("Build the sdist and wheel before preparing the candidate")
    out = args.output_dir.resolve()
    recipe_dir = out / "conda-recipe"
    recipe_dir.mkdir(parents=True, exist_ok=True)
    recipe = (root / "conda-recipe/meta.yaml").read_text()
    recipe = re.sub(r'{% set version = "[^"]+" %}', '{% set version = "' + version + '" %}', recipe)
    recipe = re.sub(r'^  url: .*$', '  url: "' + sdist.as_uri() + '"', recipe, flags=re.M)
    recipe = re.sub(r'^  sha256: .*$', '  sha256: ' + digest(sdist), recipe, flags=re.M)
    recipe = recipe.replace('    - mummer4\n', '    - mummer4\n    - samtools\n')
    recipe = recipe.replace('    - chromo --help\n', '    - chromo --help\n    - chromo --version\n    - chromo reads --help\n    - chromo workflow --help\n')
    (recipe_dir / "meta.yaml").write_text(recipe)
    record = {
        "schema": "chromosort-release-preparation-v1", "version": version,
        "state": "local_candidate_not_published", "channel_target": "Bioconda (staged)",
        "artifacts": [{"path": str(p), "sha256": digest(p), "bytes": p.stat().st_size} for p in (sdist, wheel)],
        "conda_recipe": {"path": str(recipe_dir / "meta.yaml"), "sha256": digest(recipe_dir / "meta.yaml")},
        "publication_pending": ["stable version approval", "frozen benchmark and tutorial validation", "package registry upload",
                                "immutable public archive and DOI", "preprint submission"],
    }
    (out / "release.json").write_text(json.dumps(record, indent=2) + "\n")
    print(out / "release.json")


if __name__ == "__main__":
    main()
