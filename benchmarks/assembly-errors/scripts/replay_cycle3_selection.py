#!/usr/bin/env python3
"""Supplementary D201 oracle event selection using emitted interval endpoints."""
import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'scripts')]
from benchio import digest,read_tsv,save_json,write_tsv
from adapt_agp import read_agp,FIELDS
from score_review_proposals import match_distance


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--python',type=Path,required=True)
    p.add_argument('--role',choices=['baseline','candidate'],required=True)
    p.add_argument('--software-root',type=Path,required=True)
    p.add_argument('--panel-root',type=Path,required=True)
    p.add_argument('--paired-commitment',type=Path,required=True)
    p.add_argument('--output-root',type=Path,required=True)
    args=p.parse_args();root=args.software_root.resolve();out=args.output_root.resolve()
    if out.exists():p.error('Choose a new replay directory')
    pair=json.loads(args.paired_commitment.read_text())
    for rel,sha in pair['software'][args.role]['files'].items():assert digest(root/rel)==sha
    case=args.panel_root.resolve()/'development/D201';scan=root/'review/D201/scan'
    normalized=json.loads((root/'review/D201/normalized-events.json').read_text())
    geometry=next(c for c in json.loads((ROOT/'cycle3/cases.json').read_text())['cases'] if c['case_id']=='D201')['error_interval']
    error=dict(contig='c001',shape='interval',start=geometry[0],end=geometry[1])
    candidates=sorted((match_distance(e,error,500),e['event_id'],e) for e in normalized['events']
        if e['queue'] and e['shape']=='interval' and match_distance(e,error,500) is not None)
    out.mkdir(parents=True)
    record=dict(schema='cycle3-scripted-selection-v1',role=args.role,case_id='D201',
        classification='Oracle event selection against frozen error interval, using emitted endpoints without refinement; not blinded judgment or automatic discovery',
        normalized_events_sha256=digest(root/'review/D201/normalized-events.json'),paired_commitment_sha256=digest(args.paired_commitment),
        status='unavailable_no_matching_interval')
    if not candidates:
        save_json(out/'run.json',record);save_json(out/'outputs.commitment.json',{'run.json':digest(out/'run.json')})
        print(args.role,record['status']);return
    event=candidates[0][2];record['selected_event']=event
    env=dict(os.environ);env.pop('PYTHONPATH',None);env['PYTHONDONTWRITEBYTECODE']='1'
    commands=[]
    def run(options,label):
        command=[str(args.python.absolute()),'-m','chromosort.cli',*map(str,options)]
        with (out/(label+'.log')).open('w') as log:
            result=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,env=env)
        commands.append(dict(command=['selected-python',*[s.replace(str(case),'CASE').replace(str(scan),'SCAN').replace(str(out),'OUTPUT') for s in command[1:]]],returncode=result.returncode))
        if result.returncode:raise RuntimeError(label+' failed')
    try:
        run(['workflow','edit','--scan-dir',scan,'--event-id',event['event_id'],'--action','reverse',
            '--start',event['start'],'--end',event['end'],'--notes','Oracle event selection; emitted interval endpoints unchanged; scripted feasibility diagnostic.',
            '--output-dir',out/'edit'],'edit')
        edits=read_tsv(out/'edit/reviewed_edits.tsv')
        for row in edits:row.update(decision='accept',reviewer='scripted-oracle-event-selection')
        write_tsv(out/'reviewed-edits.tsv',edits,list(edits[0]))
        decisions=read_tsv(scan/'decisions.tsv')
        for row in decisions:
            row.update(decision='accept' if row['event_id']==event['event_id'] else 'reject',reviewer='scripted-oracle-event-selection',
                notes='Acknowledge selected inspect item; explicit emitted interval edit only; no other edits.')
        write_tsv(out/'decisions.tsv',decisions,list(decisions[0]))
        common=['workflow','apply','--scan-dir',scan,'--decisions',out/'decisions.tsv','--reviewed-edits',out/'reviewed-edits.tsv']
        run([*common,'--plan-only','--output-dir',out/'preview'],'preview')
        assert not (out/'preview/reviewed.fa').exists()
        run([*common,'--output-dir',out/'applied'],'apply')
        run(['manual','apply','--assembly-fasta',case/'inputs/assembly.fa','--recipe',out/'applied/recipe.json',
            '--output-fasta',out/'replay.fa'],'replay')
        output=out/'applied/reviewed.fa'
        record['replay_identical']=output.read_bytes()==(out/'replay.fa').read_bytes()
        assert record['replay_identical']
        run(['workflow','align','--manifest',case/'manifest.json','--assembly-fasta',output,
            '--output-dir',out/'realigned','--stage','scripted-emitted-interval','--threads',2],'align')
        run(['workflow','validate','--manifest',out/'realigned/manifest.json','--apply-audit',out/'applied/apply.json',
            '--output',out/'validation.json'],'validate')
        record['validation_state']=json.loads((out/'validation.json').read_text())['state']
        shutil.copyfile(output,out/'output.fa');shutil.copyfile(out/'applied/reviewed.fa.agp',out/'output.agp')
        write_tsv(out/'components.tsv',read_agp(out/'output.agp'),FIELDS)
        record.update(status='complete',output_sha256=digest(out/'output.fa'),ledger_sha256=digest(out/'components.tsv'))
    except Exception as failure:record.update(status='failed',error=str(failure))
    record['commands']=commands;save_json(out/'run.json',record)
    save_json(out/'outputs.commitment.json',{str(p.relative_to(out)):digest(p) for p in sorted(out.rglob('*')) if p.is_file()})
    print(args.role,'D201 emitted-coordinate replay',record['status'])


if __name__=='__main__':main()
