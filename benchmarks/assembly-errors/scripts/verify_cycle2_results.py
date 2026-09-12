#!/usr/bin/env python3
"""Recompute published metrics from portable input provenance ledgers, without tools."""
import argparse
import json
import sys
import tempfile
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'src'), str(ROOT / 'scripts')]
from benchio import digest, read_fasta, read_tsv, rc, save_json, write_fasta
from score_structure import score_structure


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--panel-root', type=Path, required=True)
    p.add_argument('--prior-root', type=Path, required=True)
    p.add_argument('--out', type=Path)
    args = p.parse_args()
    results = ROOT / 'cycle2/results'
    index = json.loads((results / 'artifact-index.json').read_text())
    for rel, record in index['files'].items():
        assert digest(results / rel) == record['sha256'], rel
    primary = json.loads((ROOT / 'cycle2/protocol.json').read_text())['case_ids']
    total = 0
    for line in (results / 'scores.jsonl').read_text().splitlines():
        record = json.loads(line)
        if record['status'] != 'complete':
            continue
        cid = record['case_id']
        panel = args.panel_root if cid in primary else args.prior_root
        case, truth = panel / 'development' / cid, panel / 'evaluator/development' / cid
        assembly_path = case / 'inputs/assembly.fa'
        assembly = read_fasta(assembly_path)
        ledger = results / 'artifacts' / record['role'] / cid / record['arm'] / 'components.tsv'
        outputs = {}
        for row in read_tsv(ledger):
            if row['input_id'] == '.':
                sequence = 'N' * (int(row['output_end']) - int(row['output_start']))
            else:
                sequence = assembly[row['input_id']][int(row['input_start']):int(row['input_end'])]
                if row['strand'] == '-':
                    sequence = rc(sequence)
            previous = outputs.get(row['output_id'], '')
            assert len(previous) == int(row['output_start'])
            outputs[row['output_id']] = previous + sequence
        with tempfile.TemporaryDirectory() as temporary:
            fasta = Path(temporary) / 'output.fa'
            write_fasta(fasta, outputs)
            actual = score_structure(truth, assembly_path, fasta, ledger)
        expected = record['structure_score']
        # Tool FASTA headers/wrapping are not recoverable from a base ledger.
        # Every sequence/provenance/structural metric and ledger hash must match.
        for item in [actual, expected]:
            item['frozen_score'].pop('output_sha256')
        assert actual == expected, (record['role'], cid, record['arm'])
        total += 1
    result = dict(schema='cycle2-portable-score-verification-v1', reconstructed_outputs=total,
        exact_artifacts=len(index['files']), metric_matches=total,
        note='Reconstruct sequence from public component ledgers; all metrics and ledger hashes match. Original FASTA byte hashes remain recorded separately; header and wrapping bytes are not reconstructed.')
    if args.out:
        save_json(args.out, result)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
