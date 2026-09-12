#!/usr/bin/env python3
"""Add full target-path completeness to the unchanged interval scorer."""
import argparse
import collections
import json
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'src'), str(ROOT / 'scripts')]
from benchio import read_fasta, read_tsv, save_json
from score import score


def complete_paths(target_lengths, output_lengths, projected):
    """Require one gap-free, full-length oriented output per target identity.

    Input is a projection ledger, not an alignment; identical repeat sequences
    cannot substitute for a different target copy. Records may be reordered and
    whole target records may be reverse-complemented.
    """
    matches = collections.defaultdict(list)
    details = []
    for oid, length in output_lengths.items():
        rows = sorted(projected.get(oid, []), key=lambda r: r['output_start'])
        valid = bool(rows)
        tid = rows[0]['target_id'] if rows else None
        strand = rows[0]['strand'] if rows else None
        cursor = 0 if strand == '+' else target_lengths.get(tid, 0)
        out_cursor = 0
        for row in rows:
            start = row['start'] if strand == '+' else row['end']
            end = row['end'] if strand == '+' else row['start']
            valid = valid and row['target_id'] == tid and row['strand'] == strand
            valid = valid and row['output_start'] == out_cursor and start == cursor
            valid = valid and row['output_end']-row['output_start'] == row['end']-row['start']
            cursor, out_cursor = end, row['output_end']
        valid = bool(valid and out_cursor == length and length == target_lengths.get(tid)
                     and cursor == (target_lengths[tid] if strand == '+' else 0))
        if valid:
            matches[tid].append(oid)
        details.append(dict(output_id=oid, complete_target_id=tid if valid else None,
                            complete_target_path=valid))
    exact = {tid: len(matches[tid]) == 1 for tid in target_lengths}
    return dict(schema='benchmark-complete-structure-v1', target_count=len(target_lengths),
        complete_target_count=sum(exact.values()), targets=exact, outputs=details,
        complete_structure=all(exact.values()) and len(output_lengths) == len(target_lengths)
            and all(row['complete_target_path'] for row in details),
        interpretation='Exact full target paths by provenance, including masked bases; fragmentation is incomplete. Whole-record reverse complement and output record order are free. This is stricter than cut localization.')


def score_structure(truth_dir, assembly_path, output_path, ledger_path):
    truth_dir = Path(truth_dir)
    # The frozen scorer validates exact FASTA and both ledgers before projection.
    base = score(truth_dir, assembly_path, output_path, ledger_path)
    truth = json.loads((truth_dir / 'truth.json').read_text())
    maps = collections.defaultdict(list)
    for row in read_tsv(truth_dir / truth['mapping']):
        maps[row['input_id']].append(row)
    projected = collections.defaultdict(list)
    for component in read_tsv(ledger_path):
        if component['input_id'] == '.':
            continue
        s, e = int(component['input_start']), int(component['input_end'])
        pos = int(component['output_start'])
        forward = component['strand'] == '+'
        for row in maps[component['input_id']]:
            ms, me = int(row['input_start']), int(row['input_end'])
            a, b = max(s, ms), min(e, me)
            if a >= b:
                continue
            ts, te = int(row['target_start']), int(row['target_end'])
            x, y = (ts+a-ms, ts+b-ms) if row['strand'] == '+' else (te-b+ms, te-a+ms)
            out = pos + (a-s if forward else e-b)
            projected[component['output_id']].append(dict(target_id=row['target_id'], start=x,
                end=y, strand='+' if component['strand'] == row['strand'] else '-',
                output_start=out, output_end=out+b-a))
    result = complete_paths({k: len(v) for k, v in read_fasta(truth_dir / 'target.fa').items()},
        {k: len(v) for k, v in read_fasta(output_path).items()}, projected)
    result.update(case_id=truth['case_id'], frozen_score=base)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--truth', type=Path, required=True)
    p.add_argument('--assembly', type=Path, required=True)
    p.add_argument('--output-fasta', type=Path, required=True)
    p.add_argument('--ledger', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    save_json(args.out, score_structure(args.truth, args.assembly, args.output_fasta, args.ledger))
