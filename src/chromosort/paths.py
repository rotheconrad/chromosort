"""Filesystem helpers for ChromoSort outputs."""

from pathlib import Path


def ensure_parent_dir(path):
    """Create the parent directory for an output path when one was supplied."""
    if path is None:
        return None
    path = Path(path)
    parent = path.parent
    if parent and str(parent) != ".":
        parent.mkdir(parents=True, exist_ok=True)
    return path


def ensure_output_dirs(paths):
    """Create parent directories for a mapping or iterable of output paths."""
    values = paths.values() if hasattr(paths, "values") else paths
    for path in values:
        ensure_parent_dir(path)
