#!/usr/bin/env python3
"""Reproduce public cycle3 proposal and sequence metrics without ChromoSort."""
import argparse
import json
import sys
import tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'scripts')]
from benchio import digest,read_fasta,read_tsv,rc,save_json,write_fasta
from score_review_proposals import score_review
from score_structure import score_structure


def reconstruct_and_score(assembly_path,truth,ledger):
    assembly=read_fasta(assembly_path);outputs={}
    for row in read_tsv(ledger):
        s,e=int(row['input_start']),int(row['input_end'])
        if row['input_id']=='.':sequence='N'*(int(row['output_end'])-int(row['output_start']))
        else:
            sequence=assembly[row['input_id']][s:e]
            if row['strand']=='-':sequence=rc(sequence)
        previous=outputs.get(row['output_id'],'')
        assert len(previous)==int(row['output_start'])
        outputs[row['output_id']]=previous+sequence
    with tempfile.TemporaryDirectory() as temporary:
        fasta=Path(temporary)/'output.fa';write_fasta(fasta,outputs)
        return score_structure(truth,assembly_path,fasta,ledger)


def assert_metrics(actual,expected):
    expected=json.loads(json.dumps(expected))
    actual['frozen_score'].pop('output_sha256');expected['frozen_score'].pop('output_sha256')
    assert actual==expected


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--panel-root',type=Path,required=True);p.add_argument('--known-root',type=Path,required=True)
    p.add_argument('--results-root',type=Path,default=ROOT/'cycle3/results');p.add_argument('--out',type=Path)
    args=p.parse_args();results=args.results_root
    index=json.loads((results/'artifact-index.json').read_text())
    for rel,item in index['files'].items():assert digest(results/rel)==item['sha256']
    reviews=autos=replays=unavailable_reviews=unavailable_autos=0
    for line in (results/'proposal-scores.jsonl').read_text().splitlines():
        record=json.loads(line)
        if record['status']!='complete':
            assert record['score'] is None;unavailable_reviews+=1;continue
        panel=args.panel_root if record['cohort']=='primary' else args.known_root
        truth=json.loads((panel/'evaluator/development'/record['case_id']/'review_truth.json').read_text())
        normalized=json.loads((results/record['normalized_events']).read_text())
        assert truth['assembly_sha256']==normalized['assembly_sha256']
        assert score_review(truth,normalized['events'])==record['score'];reviews+=1
    for line in (results/'automatic-scores.jsonl').read_text().splitlines():
        record=json.loads(line)
        if record['status']!='complete':
            assert record['score'] is None;unavailable_autos+=1;continue
        panel=args.panel_root if record['cohort']=='primary' else args.known_root;cid=record['case_id']
        ledger=results/'artifacts'/record['role']/cid/('chromosort_'+record['policy']+'_comprehensive')/'components.tsv'
        actual=reconstruct_and_score(panel/'development'/cid/'inputs/assembly.fa',panel/'evaluator/development'/cid,ledger)
        assert_metrics(actual,record['score']);autos+=1
    for record in json.loads((results/'scripted-selection.json').read_text()):
        if record['status']!='complete':continue
        ledger=results/'artifacts'/record['role']/'D201/scripted_selection/components.tsv'
        actual=reconstruct_and_score(args.panel_root/'development/D201/inputs/assembly.fa',args.panel_root/'evaluator/development/D201',ledger)
        assert_metrics(actual,record['score']);replays+=1
    value=dict(schema='cycle3-portable-metric-verification-v1',proposal_scores_reproduced=reviews,automatic_scores_reproduced=autos,
        scripted_scores_reproduced=replays,unavailable_reviews=unavailable_reviews,unavailable_automatic=unavailable_autos,
        exact_small_artifacts=len(index['files']),note='Public event geometries reproduce proposal metrics; input-provenance ledgers reproduce sequence metrics. Original FASTA header/wrapping byte hashes remain recorded separately.')
    if args.out:save_json(args.out,value)
    print(json.dumps(value,indent=2))


if __name__=='__main__':main()
