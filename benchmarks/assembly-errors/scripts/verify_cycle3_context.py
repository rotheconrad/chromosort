#!/usr/bin/env python3
"""Verify public rc4 input, mapping and coordinate-evidence identities independently."""
import argparse
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from benchio import digest,save_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--panel-root',type=Path,required=True);p.add_argument('--known-root',type=Path,required=True)
    p.add_argument('--results-root',type=Path,default=ROOT/'cycle3/results');p.add_argument('--out',type=Path)
    args=p.parse_args();cases=inputs=blocks=events=coordinates=0
    for line in (args.results_root/'proposal-scores.jsonl').read_text().splitlines():
        record=json.loads(line)
        if record['role']!='candidate' or record['status']!='complete':continue
        panel=args.panel_root if record['cohort']=='primary' else args.known_root
        case=panel/'development'/record['case_id'];manifest=json.loads((case/'manifest.json').read_text())
        normalized=json.loads((args.results_root/record['normalized_events']).read_text());context=normalized['producer_context']
        assembly_sha=digest(case/'inputs/assembly.fa')
        assert context['assembly_sha256']==manifest['assembly']['sha256']==assembly_sha
        assert context['manifest_sha256']==digest(case/'manifest.json')
        input_shas=set()
        for item in context['inputs']:
            assert item['path'].startswith('CASE/'),item['path']
            path=case/item['path'][5:]
            assert digest(path)==item['sha256'] and path.stat().st_size==item['bytes']
            input_shas.add(item['sha256']);inputs+=1
        for block in context['blocks']:
            assert block['reference']['sha256'] in input_shas
            for source in block['sources']:
                assert source['assembly_sha256']==assembly_sha
                assert source['sha256'] in input_shas and source['reference_sha256'] in input_shas
            blocks+=1
        for event in normalized['events']:
            detail=event['detail']
            if not detail:continue
            assert detail['assembly_sha256']==assembly_sha and detail['implicit_edits'] is False
            for coordinate in detail['coordinate_evidence']:
                reads=coordinate['reads'];graph=coordinate['graph']
                assert reads['position']==coordinate['position'] and reads['assembly_sha256']==graph['assembly_sha256']==assembly_sha
                assert reads['evidence_sha256'] in input_shas
                assert graph['status']=='unavailable_no_graph' and graph['evidence_id'] is None and graph['evidence_sha256'] is None
                coordinates+=1
            events+=1
        cases+=1
    value=dict(schema='cycle3-public-context-verification-v1',status='pass',cases=cases,input_records=inputs,
        mapping_blocks=blocks,structured_events=events,coordinate_evidence_records=coordinates,
        interpretation='Exact stage and evidence byte identities checked against regenerated public inputs; no biological correctness inferred from explanatory labels.')
    if args.out:save_json(args.out,value)
    print(json.dumps(value,indent=2))


if __name__=='__main__':main()
