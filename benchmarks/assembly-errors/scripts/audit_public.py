#!/usr/bin/env python3
"""Audit the public selection and write a path/size/SHA-256 manifest."""
import argparse
import collections
import hashlib
import json
import re
from pathlib import Path
from urllib.parse import unquote, urlsplit
ROOT = Path(__file__).resolve().parents[1]
PRIVATE_MARKERS = ['/' + 'Users' + '/', 'One' + 'Drive-', 'file:' + '///', 'Library/' + 'CloudStorage']
TEXT = {'.py', '.md', '.json', '.jsonl', '.tsv', '.txt', '.yaml', '.yml', '.paf', '.fa'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    files, groups = {}, collections.defaultdict(lambda: dict(files=0, bytes=0))
    for path in sorted(ROOT.rglob('*')):
        rel = path.relative_to(ROOT)
        if rel.parts[0] == 'work' or '__pycache__' in rel.parts or path.suffix in {'.pyc', '.fai'} or '.ruff_cache' in rel.parts or path.name == '.DS_Store' or rel.as_posix() == 'PUBLIC_MANIFEST.json':
            continue
        if path.is_symlink():
            raise ValueError('Symlink in selection: ' + str(rel))
        if not path.is_file():
            continue
        if any(p in {'heldout', 'vendor', 'downloads', 'dist', '.venv', '.git', '.codex'} for p in rel.parts):
            raise ValueError('Excluded directory in selection: ' + str(rel))
        if re.search(r'(^|/)H[0-9]{3}(/|$)', str(rel)) or path.suffix in {'.bam', '.bai', '.enc', '.key', '.whl', '.gz', '.zip'}:
            raise ValueError('Excluded data/archive in selection: ' + str(rel))
        raw = path.read_bytes()
        if path.suffix in TEXT:
            text = raw.decode('utf-8')
            if any(marker in text for marker in PRIVATE_MARKERS):
                raise ValueError('Personal/local path in selection: ' + str(rel))
            if path.suffix == '.md':
                for target in re.findall(r'\]\(([^\s)]+)\)', text):
                    url = urlsplit(target.strip('<>'))
                    if url.scheme or not url.path:
                        continue
                    if not (path.parent / unquote(url.path)).exists():
                        raise ValueError('Broken local link: ' + str(rel) + ' -> ' + target)
        files[str(rel)] = dict(bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())
        groups[rel.parts[0]]['files'] += 1
        groups[rel.parts[0]]['bytes'] += len(raw)
    record = dict(schema='public-benchmark-content-manifest-v1', files=files, groups=dict(groups),
                  file_count=len(files), total_bytes=sum(r['bytes'] for r in files.values()),
                  excludes='Ignored work, Python caches and this self-referential manifest; no staging/commit/push is performed.')
    if args.write:
        (ROOT / 'PUBLIC_MANIFEST.json').write_text(json.dumps(record, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k: record[k] for k in ['file_count', 'total_bytes', 'groups']}, indent=2))


if __name__ == '__main__':
    main()
