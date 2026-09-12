#!/usr/bin/env python3
"""Recompute compact development metrics from original, checksummed score records."""
import collections
import hashlib
import json
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from benchio import write_tsv


def main():
    groups = collections.defaultdict(collections.Counter)
    cases = collections.defaultdict(set)
    seen = set()
    for line in (ROOT / 'reports/historical_scores.jsonl').read_text().splitlines():
        record = json.loads(line)
        run, score = record['run'], record['score']
        key = (record['revision'], run['arm'], run['operation'])
        identity = key + (run['case_id'],)
        if identity in seen or not run['case_id'].startswith('D'):
            raise ValueError('Duplicate or non-development result')
        seen.add(identity)
        original = (json.dumps(score, indent=2, sort_keys=True) + '\n').encode()
        if hashlib.sha256(original).hexdigest() != record['original_score_json_sha256']:
            raise ValueError('Historical score identity changed')
        if run['output_sha256'] != score['output_sha256'] or run['status'] != 'complete':
            raise ValueError('Incomplete or mismatched original run')
        cases[key].add(run['case_id'])
        groups[key]['runs'] += 1
        groups[key]['retain_all_exact_runs'] += score['retain_all_exact']
        for name in ['new_breakpoints', 'expected_error_breakpoints', 'matched_error_breakpoints',
                     'false_breakpoints', 'input_unaccounted_bp', 'input_duplicated_bp',
                     'target_relation_discrepancy_count']:
            groups[key][name] += score[name]
    if len(seen) != 294 or any(len(ids) != 14 for ids in cases.values()):
        raise ValueError('Expected 294 development scores, 14 cases per arm')
    rows = [dict(revision=revision, arm=arm, operation=operation, **groups[(revision, arm, operation)])
            for revision, arm, operation in sorted(groups)]
    write_tsv(ROOT / 'reports/summary.tsv', rows, list(rows[0]))
    print('Verified 294 original score identities and summarized', len(rows), 'development arms')


if __name__ == '__main__':
    main()
