#!/usr/bin/env python3
"""Normalize baseline scan rows and rc4 review details to independent native events."""
import argparse
import json
import re
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from benchio import digest,save_json


def portable(value,case,scan):
    if isinstance(value,dict):return {k:portable(v,case,scan) for k,v in value.items()}
    if isinstance(value,list):return [portable(v,case,scan) for v in value]
    if isinstance(value,str) and value.startswith('/'):
        path=Path(value)
        for name,root in [('CASE',case),('SCAN',scan)]:
            try:return name+'/'+str(path.relative_to(root))
            except ValueError:pass
        return 'EXTERNAL/'+path.name
    return value


def normalize(case,scan):
    case,scan=Path(case).resolve(),Path(scan).resolve()
    manifest=json.loads((case/'manifest.json').read_text())
    record=json.loads((scan/'scan.json').read_text())
    assert record['assembly']['sha256']==manifest['assembly']['sha256']
    assert record['manifest']['sha256']==digest(case/'manifest.json')
    details,context={},None
    path=scan/'review_proposals.json'
    if path.exists():
        context=json.loads(path.read_text())
        assert context['schema']=='chromosort-review-proposals-v1'
        assert context['assembly_sha256']==manifest['assembly']['sha256']
        assert context['manifest_sha256']==digest(case/'manifest.json')
        for field in ['inputs','parameters','blocks','contigs','events']:assert field in context
        block_ids={block['block_id'] for block in context['blocks']}
        assert len(block_ids)==len(context['blocks'])
        details={event['event_id']:event for event in context['events']}
        assert len(details)==len(context['events'])
        for detail in details.values():
            for field in ['supporting_block_ids','opposing_block_ids','overlapping_block_ids']:
                assert set(detail[field])<=block_ids
            assert [e['position'] for e in detail['coordinate_evidence']]==detail['positions']
    events=[]
    for row in record['proposals']:
        assert row['assembly_sha256']==manifest['assembly']['sha256']
        eid=row['event_id'];detail=details.get(eid)
        if detail:
            assert detail['target']==row['target'] and detail['decision']==row['decision']=='pending'
            assert detail['assembly_sha256']==manifest['assembly']['sha256'] and detail['action']==row['action']=='inspect'
            assert detail['implicit_edits'] is False
            kind=detail['kind'];positions=list(detail['positions'])
            if kind=='orientation_discrepancy':
                shape='interval';s,e=int(detail['start0']),int(detail['end0'])
                assert positions==[s,e]
            elif kind=='alternative_segmentation':
                assert positions==sorted(set(positions)) and positions
                shape='point' if len(positions)==1 else 'plan'
                s,e=positions[0],positions[-1]
            else:raise ValueError('Unknown review-proposal kind: '+kind)
            explanation='; '.join(map(str,detail.get('loss_stages',[])+detail.get('reasons',[])))
        elif row['action']=='cut':
            kind,shape,s,e,positions,explanation='automatic_cut_proposal','point',int(row['position']),int(row['position']),[int(row['position'])],row['reason']
        elif row['action']=='keep':
            kind,shape,s,e,positions,explanation='keep_acknowledgment','point',0,0,[],row['reason']
        else:
            match=re.search(r'internal interval \[(\d+), (\d+)\)',row['reason'])
            if row['action']!='inspect' or not match:raise ValueError('Unsupported baseline inspection geometry: '+eid)
            kind,shape='unaligned_gap','interval';s,e=map(int,match.groups());positions=[s,e];explanation=row['reason']
        event=dict(event_id=eid,contig=row['target'],kind=kind,shape=shape,start=s,end=e,positions=positions,
            queue=row['action']!='keep',decision=row['decision'],stage_explanation=explanation,
            assembly_sha256=manifest['assembly']['sha256'],detail=portable(detail,case,scan) if detail else None)
        events.append(event)
    assert set(details)<=set(e['event_id'] for e in events)
    return dict(schema='cycle3-normalized-review-events-v1',case_id=case.name,assembly_sha256=manifest['assembly']['sha256'],
        manifest_sha256=digest(case/'manifest.json'),scan_sha256=digest(scan/'scan.json'),
        structured_proposals_sha256=digest(path) if path.exists() else None,events=events,
        producer_context=portable({k:v for k,v in context.items() if k!='events'},case,scan) if context else None,
        normalization='One row per scan decision ID; structured details supply geometry, never duplicated as extra events. Whole plan counts once. Paths converted to portable labels; original artifact hashes retained.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',type=Path,required=True)
    parser.add_argument('--scan',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();save_json(args.out,normalize(args.case,args.scan))
