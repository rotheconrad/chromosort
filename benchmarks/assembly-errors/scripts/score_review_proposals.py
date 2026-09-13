#!/usr/bin/env python3
"""Independent native-geometry proposal coverage and review-burden scoring."""
import argparse
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from benchio import save_json


def close(a,b,tolerance):
    if a['contig']!=b['contig'] or a['shape']!=b['shape']:return False
    if a['shape']=='plan':
        return len(a['positions'])==len(b['positions']) and all(abs(x-y)<=tolerance for x,y in zip(a['positions'],b['positions']))
    return abs(a['start']-b['start'])<=tolerance and abs(a['end']-b['end'])<=tolerance


def points(event):
    return event['positions'] if event['shape']=='plan' else sorted({event['start'],event['end']})


def match_distance(event,error,tolerance):
    if event['contig']!=error['contig']:return None
    if error['shape']=='point' and event['shape'] in ['point','plan']:
        distance=min(abs(p-error['start']) for p in points(event))
        return distance if distance<=tolerance else None
    if close(event,error,tolerance):
        return abs(event['start']-error['start'])+abs(event['end']-error['end'])
    return None


def overlap(event,start,end):
    if event['shape'] in ['point','plan']:
        return any(start<=p<end for p in points(event))
    return max(start,event['start'])<min(end,event['end'])


def merge(intervals):
    result=[]
    for a,b in sorted(intervals):
        if result and a<=result[-1][1]:result[-1]=(result[-1][0],max(b,result[-1][1]))
        else:result.append((a,b))
    return result


def score_review(truth,events,tolerance=500):
    ids=[e['event_id'] for e in events]
    if len(ids)!=len(set(ids)):raise ValueError('Duplicate event IDs in normalized queue')
    queue=[e for e in events if e['queue']]
    for event in queue:
        length=truth['contig_lengths'][event['contig']]
        s,e=event['start'],event['end']
        if event['shape'] not in ['point','interval','plan'] or not 0<=s<=e<=length:
            raise ValueError('Invalid native event geometry')
        if (event['shape']=='point')!=(s==e):raise ValueError('Point/interval mismatch')
        if event['shape']=='plan' and (len(event['positions'])<2 or event['positions']!=sorted(set(event['positions'])) or points(event)[0]!=s or points(event)[-1]!=e):
            raise ValueError('Invalid ordered plan cut positions')
    errors=truth['errors']
    candidates=sorted((match_distance(e,t,tolerance),e['event_id'],t['error_id'],i,j)
        for i,e in enumerate(queue) for j,t in enumerate(errors) if match_distance(e,t,tolerance) is not None)
    matched_events,matched_errors,matches=set(),set(),[]
    for distance,eid,tid,i,j in candidates:
        if i not in matched_events and j not in matched_errors:
            matched_events.add(i);matched_errors.add(j)
            matches.append(dict(event_id=eid,error_id=tid,total_edge_distance_bp=distance))
    edge_total=edge_covered=0
    for error in errors:
        error_points=[error['start']] if error['shape']=='point' else [error['start'],error['end']]
        for point in error_points:
            edge_total+=1
            edge_covered+=any(e['contig']==error['contig'] and any(abs(x-point)<=tolerance for x in points(e)) for e in queue)
    # Greedy deterministic geometric deduplication, across kind and reference frame.
    representatives=[]
    for event in sorted(queue,key=lambda e:(e['contig'],e['shape'],e['start'],e['end'],e['event_id'])):
        if not any(close(event,r,tolerance) for r in representatives):representatives.append(event)
    bio=truth['biological_controls']
    bio_events=[e['event_id'] for e in queue if any(e['contig']==b['contig'] and overlap(e,b['start'],b['end']) for b in bio)]
    bio_hit=[b['control_id'] for b in bio if any(e['contig']==b['contig'] and overlap(e,b['start'],b['end']) for e in queue)]
    clean={}
    for contig,length in truth['contig_lengths'].items():
        excluded=[(max(0,e['start']-tolerance),min(length,e['end']+tolerance+1)) for e in errors if e['contig']==contig]
        excluded += [(b['start'],b['end']) for b in bio if b['contig']==contig]
        cursor=0;spans=[]
        for a,b in merge(excluded):
            if cursor<a:spans.append((cursor,a))
            cursor=max(cursor,b)
        if cursor<length:spans.append((cursor,length))
        clean[contig]=spans
    clean_events=[e['event_id'] for e in queue if any(overlap(e,a,b) for a,b in clean[e['contig']])]
    kinds={kind:sum(e['kind']==kind for e in queue) for kind in sorted({e['kind'] for e in queue})}
    return dict(schema='cycle3-review-proposal-score-v1',case_id=truth['case_id'],
        normalized_rows=len(events),queue_size=len(queue),acknowledgment_rows=len(events)-len(queue),
        pending_queue=sum(e['decision']=='pending' for e in queue),by_kind=kinds,
        expected_errors=len(errors),covered_errors=len(matched_errors),
        error_proposal_coverage=len(matched_errors)/len(errors) if errors else None,
        expected_error_edges=edge_total,covered_error_edges=edge_covered,matches=matches,
        unique_geometric_events=len(representatives),redundant_geometric_events=len(queue)-len(representatives),
        biological_control_count=len(bio),biological_controls_inspected=len(bio_hit),
        biological_burden_events=len(bio_events),biological_event_ids=bio_events,
        clean_burden_events=len(clean_events),clean_event_ids=clean_events,
        clean_control=truth.get('clean_control',False),whole_record_orientation_control=truth.get('whole_record_orientation_control',False),
        stage_explanations_available=sum(bool(e.get('stage_explanation')) for e in queue),
        interpretation='Proposal coverage and burden, not biological classification or false correction. Error matching requires one event with matching shape and both native edges; edge-only coverage is separate. Redundancy ignores event kind/reference frame.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--truth',type=Path,required=True)
    p.add_argument('--events',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    save_json(args.out,score_review(json.loads(args.truth.read_text()),json.loads(args.events.read_text())['events']))
