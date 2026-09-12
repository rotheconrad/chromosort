#!/usr/bin/env python3
"""Verify public fixture bytes, exact truth provenance and mapped PAF geometry."""
import argparse
import gzip
import hashlib
import json
import re
import sys
import tempfile
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'src'), str(ROOT / 'scripts')]
from benchio import digest, read_fasta, read_tsv, write_tsv
from adapt_agp import FIELDS
from score import score


def check_paf(path, queries, targets):
    count = 0
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        f = line.split('\t')
        assert len(f) >= 12, (path, 'short PAF')
        ql, qs, qe, tl, ts, te = map(int, [f[1], f[2], f[3], f[6], f[7], f[8]])
        assert queries[f[0]] == ql and targets[f[5]] == tl, (path, 'sequence lengths')
        assert 0 <= qs < qe <= ql and 0 <= ts < te <= tl and f[4] in ['+', '-']
        assert 0 <= int(f[9]) <= int(f[10]) and 0 <= int(f[11]) <= 255
        cigar = next((x[5:] for x in f[12:] if x.startswith('cg:Z:')), None)
        assert cigar and re.fullmatch(r'(?:[1-9][0-9]*[MIDNSHP=X])+', cigar), (path, 'CIGAR')
        ops = [(int(n), op) for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar)]
        assert sum(n for n, op in ops if op in 'MI=X') == qe - qs
        assert sum(n for n, op in ops if op in 'MDN=X') == te - ts
        count += 1
    return count


def read_lengths(case, manifest):
    result = {}
    for entry in manifest.get('readsets', []):
        if entry.get('path'):
            path = case / entry['path']
            assert digest(path) == entry['sha256']
            with gzip.open(path, 'rt') as handle:
                while True:
                    name = handle.readline().strip()
                    if not name:
                        break
                    sequence, plus, quality = [handle.readline().strip() for _ in range(3)]
                    assert name.startswith('@') and plus == '+' and len(sequence) == len(quality)
                    assert name[1:] not in result
                    result[name[1:]] = len(sequence)
        else:
            result.update({f'read{i:07d}': entry['read_length'] for i in range(1, entry['read_count'] + 1)})
    return result


def check_case(case, truth):
    manifest = json.loads((case / 'manifest.json').read_text())
    assert manifest['split'] == 'development' and case.name.startswith('D')
    assembly_path = case / manifest['assembly']['path']
    assert digest(assembly_path) == manifest['assembly']['sha256']
    assembly = read_fasta(assembly_path)
    al = {k: len(v) for k, v in assembly.items()}
    refs = {r['reference_id']: r for r in read_tsv(case / manifest['references'])}
    for ref in refs.values():
        assert digest(case / ref['fasta']) == ref['sha256']
    reads, count = read_lengths(case, manifest), 0
    for entry in manifest['evidence']:
        assert entry['assembly_sha256'] == manifest['assembly']['sha256']
        assert digest(case / entry['path']) == entry['sha256']
        if entry.get('index_path'):
            assert digest(case / entry['index_path']) == entry['index_sha256']
        if entry['kind'] == 'assembly_paf':
            ref = refs[entry['reference_id']]
            assert entry['reference_sha256'] == ref['sha256']
            lengths = {k: len(v) for k, v in read_fasta(case / ref['fasta']).items()}
            count += check_paf(case / entry['path'], al, lengths)
        elif entry['kind'] == 'read_paf':
            count += check_paf(case / entry['path'], reads, al)
    with tempfile.TemporaryDirectory() as temp:
        ledger = Path(temp) / 'identity.tsv'
        write_tsv(ledger, [dict(output_id=k, output_start=0, output_end=len(v), input_id=k,
                  input_start=0, input_end=len(v), strand='+') for k, v in assembly.items()], FIELDS)
        value = score(truth, assembly_path, assembly_path, ledger)
        assert value['retain_all_exact']
    return dict(case_id=case.name, input_bp=sum(al.values()), paf_rows=count, truth_provenance='exact')


def verify_fixtures():
    index = json.loads((ROOT / 'fixtures/index.json').read_text())
    for path, record in index['files'].items():
        assert digest(ROOT / path) == record['sha256'], path
        assert (ROOT / path).stat().st_size == record['bytes'], path
    return [check_case(ROOT / 'fixtures/development' / c['case_id'], ROOT / 'fixtures/truth' / c['case_id'])
            for c in index['cases']]


def verify_regenerated(root):
    expected = json.loads((ROOT / 'sources/development_expected.json').read_text())['cases']
    cases = sorted((root / 'development').glob('D*/manifest.json'))
    assert cases, 'No generated development cases'
    out = []
    for path in cases:
        cid = path.parent.name
        for rel, sha in expected[cid].items():
            if rel == 'decoded_read_fastq_sha256':
                data = gzip.decompress((path.parent / 'inputs/reads.fastq.gz').read_bytes())
                assert hashlib.sha256(data).hexdigest() == sha, (cid, 'decoded FASTQ')
            else:
                assert digest(root / rel) == sha, rel
        out.append(check_case(path.parent, root / 'evaluator/development' / cid))
    return out


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--fixtures', action='store_true')
    group.add_argument('--regenerated', type=Path)
    args = parser.parse_args()
    result = verify_fixtures() if args.fixtures else verify_regenerated(args.regenerated.resolve())
    print(json.dumps({'status': 'pass', 'cases': result}, indent=2))
