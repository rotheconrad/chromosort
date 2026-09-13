#!/usr/bin/env python3
"""Verify frozen cycle3 case inputs and native truth geometry without tool outputs."""
import argparse
import gzip
import hashlib
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'scripts')]
from benchio import digest,read_fasta,save_json
from verify import check_case


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--generated-root',type=Path,required=True)
    p.add_argument('--out',type=Path)
    args=p.parse_args();root=args.generated_root.resolve()
    freeze=json.loads((ROOT/'cycle3/prospective-freeze.json').read_text())
    for rel,sha in freeze['files'].items():assert digest(ROOT/'cycle3'/rel)==sha,rel
    for rel,sha in freeze['code'].items():assert digest(ROOT/rel)==sha,rel
    expected=json.loads((ROOT/'cycle3/generated_expected.json').read_text())
    for rel,item in expected['files'].items():
        assert digest(root/rel)==item['sha256'],rel
        assert (root/rel).stat().st_size==item['bytes']
    records=[]
    for cid,sha in expected['decoded_fastq_sha256'].items():
        assert hashlib.sha256(gzip.decompress((root/f'development/{cid}/inputs/reads.fastq.gz').read_bytes())).hexdigest()==sha
        # This verifies unedited input-to-target provenance, not candidate outputs.
        record=check_case(root/'development'/cid,root/'evaluator/development'/cid)
        truth=json.loads((root/'evaluator/development'/cid/'review_truth.json').read_text())
        assert truth['assembly_sha256']==digest(root/'development'/cid/'inputs/assembly.fa')
        record['review_errors']=len(truth['errors']);record['biological_controls']=len(truth['biological_controls'])
        records.append(record)
    result=dict(schema='cycle3-input-verification-v1',case_count=len(records),exact_files=len(expected['files']),cases=records,
        purpose='Input construction and truth provenance verification only; no candidate outputs opened')
    if args.out:save_json(args.out,result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
