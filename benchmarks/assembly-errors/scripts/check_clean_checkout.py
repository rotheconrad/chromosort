#!/usr/bin/env python3
"""Verify an isolated copy using only selected public content and explicit tools."""
import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--with-mapping', action='store_true')
    parser.add_argument('--out', type=Path, help='Save compact verification results')
    args = parser.parse_args()
    results = []
    with tempfile.TemporaryDirectory(prefix='chromosort-public-clean-') as temp:
        clean = Path(temp) / 'assembly-errors'
        for path in sorted(ROOT.rglob('*')):
            rel = path.relative_to(ROOT)
            if rel.parts[0] == 'work' or '__pycache__' in rel.parts or path.suffix in {'.pyc', '.fai'} or '.ruff_cache' in rel.parts or path.name == '.DS_Store':
                continue
            if path.is_symlink():
                raise ValueError('Unexpected symlink in benchmark selection')
            if path.is_file():
                target = clean / rel
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(path, target)
        commands = [['scripts/verify.py', '--fixtures'], ['-m', 'unittest', 'discover', '-s', 'tests', '-v'],
                    ['scripts/summarize_history.py'], ['scripts/retrieve.py', '--output-root', 'work'],
                    ['scripts/run_benchmark.py', '--tool', 'unedited', '--output-root', 'work/unedited']]
        if args.with_mapping:
            commands += [['scripts/generate.py', '--output-root', 'work'],
                         ['scripts/evidence.py', *['work/development/D%03d' % i for i in range(1, 15)], '--threads', '2'],
                         ['scripts/verify.py', '--regenerated', 'work']]
        for arguments in commands:
            process = subprocess.run([sys.executable, *arguments], cwd=clean, text=True,
                                     stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            results.append(dict(command=['python', *arguments], exit_code=process.returncode))
            print(process.stdout, end='', flush=True)
            if process.returncode:
                raise SystemExit(process.returncode)
    record = dict(schema='public-clean-copy-verification-v1', status='pass', python=sys.version.split()[0],
                  isolated_copy=True, network_used=False, with_mapping=args.with_mapping,
                  commands=results, dependency_note='Only Python plus minimap2 for full mapping; no sibling data/tool archive.')
    if args.with_mapping:
        record['minimap2'] = subprocess.check_output(['minimap2', '--version'], text=True).strip()
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(record, indent=2, sort_keys=True) + '\n')
    print('Isolated clean-copy verification passed.', flush=True)


if __name__ == '__main__':
    main()
