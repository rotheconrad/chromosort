#!/usr/bin/env python3
"""Summarize committed cycle3 scores and producer context without rerunning software."""
import argparse
import collections
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from benchio import save_json,write_tsv


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results-root',type=Path,default=ROOT/'cycle3/results')
    args=p.parse_args();root=args.results_root
    reviews=[json.loads(line) for line in (root/'proposal-scores.jsonl').read_text().splitlines()]
    autos=[json.loads(line) for line in (root/'automatic-scores.jsonl').read_text().splitlines()]
    keys=['expected_errors','covered_errors','expected_error_edges','covered_error_edges','queue_size','pending_queue',
        'biological_control_count','biological_burden_events','biological_controls_inspected','clean_burden_events',
        'redundant_geometric_events','unique_geometric_events','stage_explanations_available','acknowledgment_rows']
    summaries=[];automatic=[];events=[];bounds=[];coordinates=[]
    for role in ['baseline','candidate']:
        for cohort in ['primary','known']:
            selected=[r for r in reviews if r['role']==role and r['cohort']==cohort]
            complete=[r for r in selected if r['score'] is not None];kinds=collections.Counter()
            for r in complete:kinds.update(r['score']['by_kind'])
            summaries.append(dict(role=role,cohort=cohort,attempted_cases=len(selected),completed_cases=len(complete),
                unavailable_cases=[r['case_id'] for r in selected if r['score'] is None],by_kind=dict(kinds),
                **{k:sum(r['score'][k] for r in complete) for k in keys}))
            for policy in ['report','review']:
                selected=[r for r in autos if r['role']==role and r['cohort']==cohort and r['policy']==policy]
                complete=[r for r in selected if r['score'] is not None]
                errors=[r for r in complete if r['case_has_error']];controls=[r for r in complete if not r['case_has_error']]
                automatic.append(dict(role=role,cohort=cohort,policy=policy,attempted_cases=len(selected),completed_cases=len(complete),
                    error_cases=len(errors),complete_error_repairs=sum(r['score']['complete_structure'] for r in errors),
                    already_correct_controls=len(controls),complete_controls=sum(r['score']['complete_structure'] for r in controls),
                    exact_once_input_cases=sum(r['score']['frozen_score']['retain_all_exact'] for r in complete),
                    expected_error_edges=sum(r['score']['frozen_score']['expected_error_breakpoints'] for r in complete),
                    matched_cuts=sum(r['score']['frozen_score']['matched_error_breakpoints'] for r in complete),
                    false_cuts=sum(r['score']['frozen_score']['false_breakpoints'] for r in complete),
                    broken_controls=sum(sum(not c['continuous'] for c in r['score']['frozen_score']['preservation_controls']) for r in complete)))
    for record in reviews:
        if record['role']!='candidate' or record['score'] is None:continue
        normalized=json.loads((root/record['normalized_events']).read_text());context=normalized['producer_context']
        for contig in context['contigs']:
            orientation=contig['orientation'];alternative=contig['alternatives']
            bounds.append(dict(cohort=record['cohort'],case_id=record['case_id'],contig=contig['target'],
                canonical_blocks=contig['canonical_alignment_blocks'],merged_correction_blocks=contig['merged_correction_blocks'],
                orientation_status=orientation['status'],orientation_events=orientation.get('eligible_events'),
                orientation_omitted=orientation.get('omitted_by_event_cap'),alternative_status=alternative['status'],
                alternative_paths_examined=alternative.get('raw_paths_examined'),unique_alternatives=alternative.get('unique_alternatives'),
                emitted_alternatives=alternative.get('emitted'),
                frame_contexts=';'.join(frame['frame']+':'+frame['status'] for frame in orientation.get('frames',[]))))
        for event in normalized['events']:
            if not event['queue']:continue
            detail=event['detail'] or {};objective=detail.get('objective') or {}
            events.append(dict(cohort=record['cohort'],case_id=record['case_id'],event_id=event['event_id'],kind=event['kind'],
                shape=event['shape'],positions=','.join(map(str,event['positions'])),
                matched_error_ids=','.join(m['error_id'] for m in record['score']['matches'] if m['event_id']==event['event_id']),
                biological_burden=event['event_id'] in record['score']['biological_event_ids'],
                clean_burden=event['event_id'] in record['score']['clean_event_ids'],
                loss_stages=';'.join(detail.get('loss_stages',[])),objective_rank=objective.get('rank'),
                objective_total=objective.get('total'),objective_delta=objective.get('delta_to_best_partition'),
                structured_coordinate_count=len(detail.get('coordinate_evidence',[]))))
            for coordinate in detail.get('coordinate_evidence',[]):
                reads=coordinate['reads'];graph=coordinate['graph']
                coordinates.append(dict(cohort=record['cohort'],case_id=record['case_id'],event_id=event['event_id'],
                    position=coordinate['position'],assembly_sha256=detail['assembly_sha256'],
                    read_evidence_id=reads['evidence_id'],read_evidence_sha256=reads['evidence_sha256'],
                    read_dependency_group=reads['dependency_group'],spanning_molecules=reads['spanning_molecules'],
                    read_status=reads['continuity_status'],graph_status=graph['status'],graph_evidence_id=graph['evidence_id']))
    save_json(root/'cohort-summary.json',dict(schema='cycle3-cohort-summary-v1',
        interpretation='Sums of frozen per-case metrics over completed cases only. Primary and known cohorts remain separate. No software rerun or metric change.',
        proposal=summaries,automatic=automatic))
    for name,rows in [('candidate-events.tsv',events),('candidate-bounds.tsv',bounds),('candidate-coordinate-evidence.tsv',coordinates)]:
        write_tsv(root/name,rows,list(rows[0]))
    print('Summarized',len(reviews),'review statuses and',len(autos),'automatic statuses')


if __name__=='__main__':main()
