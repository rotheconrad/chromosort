#!/usr/bin/env python3
"""Verify cycle2 source/freeze identities and regenerated inputs without scoring."""
import argparse
import gzip
import hashlib
import json
import re
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'src'), str(ROOT / 'scripts')]
from benchio import digest, read_fasta, read_tsv, save_json
from verify import read_lengths
from run_benchmark import input_hashes


def check_paf(path, queries, targets, allow_missing_cigar=False):
    count = missing = 0
    for line in path.read_text().splitlines():
        fields = line.split('\t')
        q, qlen, qs, qe, strand, t, tlen, ts, te, matches, block, mapq = fields[:12]
        qlen, qs, qe, tlen, ts, te, matches, block, mapq = map(int, [qlen, qs, qe, tlen, ts, te, matches, block, mapq])
        assert queries[q] == qlen and targets[t] == tlen
        assert 0 <= qs < qe <= qlen and 0 <= ts < te <= tlen
        assert strand in '+-' and 0 <= matches <= block and 0 <= mapq <= 255
        cigar = next((f[5:] for f in fields[12:] if f.startswith('cg:Z:')), None)
        if cigar is None:
            assert allow_missing_cigar, (path, 'unexpected missing CIGAR')
            missing += 1
        else:
            assert re.fullmatch(r'(?:[1-9][0-9]*[MIDNSHP=X])+', cigar)
            ops = [(int(n), op) for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar)]
            assert sum(n for n, op in ops if op in 'MI=X') == qe-qs
            assert sum(n for n, op in ops if op in 'MDN=X') == te-ts
        count += 1
    return count, missing


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--generated-root', type=Path, required=True)
    p.add_argument('--out', type=Path)
    args = p.parse_args()
    cycle = ROOT / 'cycle2'
    for rel, sha in json.loads((cycle / 'prospective-freeze.json').read_text())['files'].items():
        assert digest(cycle / rel) == sha, rel
    protocol = json.loads((cycle / 'protocol.json').read_text())
    assert digest(cycle / protocol['scorer']['path']) == protocol['scorer']['sha256']
    for item in json.loads((cycle / 'sources/sequence_registry.json').read_text())['sequences']:
        assert digest(cycle / item['path']) == item['sha256']
    expected = json.loads((cycle / 'generated_expected.json').read_text())
    root = args.generated_root.resolve()
    for rel, item in expected['files'].items():
        assert digest(root / rel) == item['sha256'], rel
        assert (root / rel).stat().st_size == item['bytes'], rel
    cases, total_rows, total_missing = [], 0, 0
    for cid in protocol['case_ids']:
        case = root / 'development' / cid
        manifest = json.loads((case / 'manifest.json').read_text())
        input_hashes(case, manifest)
        assembly = {k: len(v) for k, v in read_fasta(case / manifest['assembly']['path']).items()}
        reads = read_lengths(case, manifest)
        if cid == 'D106':
            reads['read9000002'] = reads['read9000001']
        refs = {r['reference_id']: r for r in read_tsv(case / manifest['references'])}
        for entry in manifest['evidence']:
            is_read = entry['kind'] == 'read_paf'
            targets = assembly if is_read else {k: len(v) for k, v in read_fasta(case / refs[entry['reference_id']]['fasta']).items()}
            rows, missing = check_paf(case / entry['path'], reads if is_read else assembly, targets, cid == 'D104' and is_read)
            total_rows += rows
            total_missing += missing
        fastq = case / manifest['readsets'][0]['path']
        assert hashlib.sha256(gzip.decompress(fastq.read_bytes())).hexdigest() == expected['decoded_fastq_sha256'][cid]
        cases.append(dict(case_id=cid, input_bp=sum(assembly.values())))
    def read_rows(cid):
        return (root / f'development/{cid}/evidence/reads.assembly.paf').read_text().splitlines()
    stripped = ['\t'.join(f for f in row.split('\t') if not f.startswith('cg:Z:')) for row in read_rows('D101')]
    assert stripped == read_rows('D104')
    repeated = [r for r in read_rows('D105') if r.startswith('read9000001\t')]
    renamed = [r for r in read_rows('D106') if r.startswith(('read9000001\t', 'read9000002\t'))]
    assert len(repeated) == 20 and len(set(repeated)) == 1
    assert len(renamed) == 2 and renamed[0].split('\t')[1:] == renamed[1].split('\t')[1:] == repeated[0].split('\t')[1:]
    result = dict(schema='cycle2-input-verification-v1', cases=cases, exact_file_count=len(expected['files']),
        mapped_paf_rows=total_rows, deliberately_missing_cigar_rows=total_missing,
        derivation_checks='D104 exact stripped D101; D105 twenty identical rows; D106 same parent row under two names',
        evaluator_access='Byte hashes only; no output scoring')
    if args.out:
        save_json(args.out, result)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
