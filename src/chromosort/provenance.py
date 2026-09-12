"""Portable byte identities and audit records for analysis inputs and outputs."""

import hashlib
import json
from pathlib import Path

from . import __version__
from .paths import ensure_parent_dir


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def check_digest(path, expected, label="input"):
    observed = sha256_file(path)
    if expected and observed != expected:
        raise ValueError(f"Stale {label}: SHA-256 mismatch for {path}; re-create evidence for this FASTA stage.")
    return observed


def stable_id(kind, *identity):
    payload = json.dumps(identity, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return f"{kind}:{hashlib.sha256(payload.encode()).hexdigest()[:24]}"


def write_json(path, value):
    ensure_parent_dir(path)
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def file_record(path, relative_to=None):
    path = Path(path).resolve()
    name = str(path)
    if relative_to:
        import os
        name = os.path.relpath(path, Path(relative_to).resolve())
    return {"path": name, "sha256": sha256_file(path), "bytes": path.stat().st_size}


def audit_record(command, inputs, outputs, settings):
    return {
        "schema": "chromosort-audit-v1", "software": {"chromosort": __version__},
        "command": command, "inputs": [file_record(p) for p in inputs],
        "outputs": [file_record(p) for p in outputs], "settings": settings,
        "validation_state": "requires_fresh_alignment_after_sequence_change",
    }
