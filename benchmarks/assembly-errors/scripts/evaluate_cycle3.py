#!/usr/bin/env python3
"""Commit the paired review outputs, run the one declared replay, then independently score."""
import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'scripts')]
from benchio import digest,save_json,write_tsv


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['baseline-root','candidate-root','baseline-python','candidate-python','panel-root','known-root','runtime-root','output-root']:
        p.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args();out=args.output_root.resolve();runtime=args.runtime_root.resolve()
    if out.exists() or runtime.exists():p.error('Choose fresh result and replay directories')
    roots={'baseline':args.baseline_root.resolve(),'candidate':args.candidate_root.resolve()}
    protocol=json.loads((ROOT/'cycle3/protocol.json').read_text())
    freeze=json.loads((ROOT/'cycle3/prospective-freeze.json').read_text())
    for rel,sha in freeze['files'].items():assert digest(ROOT/'cycle3'/rel)==sha
    for rel,sha in freeze['code'].items():assert digest(ROOT/rel)==sha
    for rel,sha in json.loads((ROOT/'cycle3/runner-freeze.json').read_text())['files'].items():assert digest(ROOT/rel)==sha
    commits={}
    for role,root in roots.items():
        value=json.loads((root/'outputs.commitment.json').read_text())
        for rel,sha in value.items():assert digest(root/rel)==sha,role+'/'+rel
        commits[role]=dict(commitment_sha256=digest(root/'outputs.commitment.json'),files=value)
    out.mkdir(parents=True);runtime.mkdir(parents=True)
    paired=out/'paired-output-commitment.json'
    save_json(paired,dict(schema='cycle3-paired-output-commitment-v1',phase='Written before selection/replay or evaluator truth access',software=commits,
        prospective_freeze_sha256=digest(ROOT/'cycle3/prospective-freeze.json'),evaluator_script_sha256=digest(__file__),
        replay_script_sha256=digest(ROOT/'scripts/replay_cycle3_selection.py')))
    for role,python in [('baseline',args.baseline_python),('candidate',args.candidate_python)]:
        subprocess.run([sys.executable,str(ROOT/'scripts/replay_cycle3_selection.py'),'--python',str(python.absolute()),
            '--role',role,'--software-root',str(roots[role]),'--panel-root',str(args.panel_root.resolve()),
            '--paired-commitment',str(paired),'--output-root',str(runtime/role)],check=True)
    replay_commits={}
    for role in roots:
        value=json.loads((runtime/role/'outputs.commitment.json').read_text())
        for rel,sha in value.items():assert digest(runtime/role/rel)==sha
        replay_commits[role]=value
    save_json(out/'replay-output-commitment.json',dict(schema='cycle3-replay-output-commitment-v1',phase='Before output structure scoring',software=replay_commits))
    from score_review_proposals import score_review
    from score_structure import score_structure
    reviews,autos,replays,review_table,auto_table,failures=[],[],[],[],[],[]
    identities={};artifacts={}
    def copy_file(source,target):
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target)
        artifacts[str(target.relative_to(out))]=dict(sha256=digest(target),original_sha256=digest(source),transformation='exact copy')
    def failure_record(role,cid,kind,directory,status,returncode=None):
        logs={}
        for name in ['stdout.log','stderr.log','scan.log']:
            path=directory/name
            if path.exists():
                content=path.read_text().replace(str(ROOT),'BENCHMARK')
                logs[name]=dict(original_sha256=digest(path),content=content,capture='combined stdout/stderr' if name=='scan.log' else name)
        failures.append(dict(role=role,case_id=cid,operation=kind,status=status,returncode=returncode,logs=logs,
            exit_status_note='Captured scan exit code' if returncode is not None else 'Legacy automatic runner records failed status and separate stdout/stderr, but not the numeric subcommand exit code'))
    for role,root in roots.items():
        identities[role]=json.loads((root/'software.json').read_text())
        for cid in protocol['primary_cases']+protocol['known_diagnostics']:
            cohort='primary' if cid in protocol['primary_cases'] else 'known'
            panel=args.panel_root.resolve() if cohort=='primary' else args.known_root.resolve()
            case=panel/'development'/cid;truth=panel/'evaluator/development'/cid
            directory=root/'review'/cid;run=json.loads((directory/'run.json').read_text())
            record=dict(role=role,case_id=cid,cohort=cohort,status=run['status'],returncode=run['returncode'],
                run_sha256=digest(directory/'run.json'),elapsed_seconds=run['elapsed_seconds'])
            if run['status']=='complete':
                normalized=json.loads((directory/'normalized-events.json').read_text())
                review_truth=json.loads((truth/'review_truth.json').read_text())
                record['score']=score_review(review_truth,normalized['events'])
                destination=out/'artifacts'/role/cid/'normalized-events.json';copy_file(directory/'normalized-events.json',destination)
                copy_file(directory/'scan/decisions.tsv',destination.parent/'pending-decisions.tsv')
                record['normalized_events']=str(destination.relative_to(out))
                score=record['score']
                review_table.append(dict(role=role,cohort=cohort,case_id=cid,status=record['status'],
                    expected_errors=score['expected_errors'],covered_errors=score['covered_errors'],expected_edges=score['expected_error_edges'],covered_edges=score['covered_error_edges'],
                    queue_size=score['queue_size'],biological_burden_events=score['biological_burden_events'],biological_controls_inspected=score['biological_controls_inspected'],
                    clean_burden_events=score['clean_burden_events'],redundant_events=score['redundant_geometric_events']))
            else:
                record['score']=None;record['interpretation']='Unavailable; excluded from numeric proposal/burden aggregates, not a zero queue or biological miss'
                failure_record(role,cid,'review_scan',directory,run['status'],run['returncode'])
            reviews.append(record)
            for policy in ['report','review']:
                arm='chromosort_'+policy+'_comprehensive';directory=root/'automatic'/cohort/cid/arm
                run=json.loads((directory/'run.json').read_text())
                record=dict(role=role,case_id=cid,cohort=cohort,policy=policy,status=run['status'],run_sha256=digest(directory/'run.json'),elapsed_seconds=run['elapsed_seconds'])
                if run['status']=='complete':
                    record['score']=score_structure(truth,case/'inputs/assembly.fa',directory/'output.fa',directory/'components.tsv')
                    record['case_has_error']=bool(json.loads((truth/'truth.json').read_text())['events'])
                    for filename in ['components.tsv','output.agp','fix.tsv']:
                        copy_file(directory/filename,out/'artifacts'/role/cid/arm/filename)
                    value=record['score'];base=value['frozen_score']
                    auto_table.append(dict(role=role,cohort=cohort,case_id=cid,policy=policy,status='complete',case_has_error=record['case_has_error'],
                        complete_structure=value['complete_structure'],matched_cuts=base['matched_error_breakpoints'],false_cuts=base['false_breakpoints'],
                        exact_once_input=base['retain_all_exact'],broken_controls=sum(not c['continuous'] for c in base['preservation_controls'])))
                else:
                    record['score']=None;failure_record(role,cid,'automatic_'+policy,directory,run['status'])
                autos.append(record)
        directory=runtime/role;run=json.loads((directory/'run.json').read_text())
        record=dict(role=role,case_id='D201',status=run['status'],classification=run['classification'],run_sha256=digest(directory/'run.json'),
            selected_event=run.get('selected_event'),replay_identical=run.get('replay_identical'),validation_state=run.get('validation_state'))
        if run['status']=='complete':
            record['score']=score_structure(args.panel_root/'evaluator/development/D201',args.panel_root/'development/D201/inputs/assembly.fa',directory/'output.fa',directory/'components.tsv')
            for name in ['components.tsv','output.agp','decisions.tsv','reviewed-edits.tsv']:
                copy_file(directory/name,out/'artifacts'/role/'D201'/'scripted_selection'/name)
            original=directory/'applied/recipe.json';recipe=json.loads(original.read_text());recipe.pop('sourceAssembly',None)
            path=out/'artifacts'/role/'D201'/'scripted_selection/recipe.json';save_json(path,recipe)
            artifacts[str(path.relative_to(out))]=dict(sha256=digest(path),original_sha256=digest(original),transformation='Remove advisory absolute sourceAssembly path; digest and emitted edit geometry unchanged')
        replays.append(record)
    for name,records in [('proposal-scores.jsonl',reviews),('automatic-scores.jsonl',autos)]:
        (out/name).write_text(''.join(json.dumps(r,sort_keys=True)+'\n' for r in records))
    write_tsv(out/'proposal-summary.tsv',review_table,list(review_table[0]));write_tsv(out/'automatic-summary.tsv',auto_table,list(auto_table[0]))
    save_json(out/'scripted-selection.json',replays);save_json(out/'failures.json',failures);save_json(out/'software-identities.json',identities)
    save_json(out/'artifact-index.json',dict(schema='cycle3-small-artifact-index-v1',files=artifacts))
    by_key={(r['role'],r['case_id'],r['policy']):r for r in autos};pairs=[]
    for cid in protocol['primary_cases']+protocol['known_diagnostics']:
        for policy in ['report','review']:
            a,b=[by_key[(role,cid,policy)] for role in roots];complete=a['status']==b['status']=='complete'
            pairs.append(dict(case_id=cid,cohort=a['cohort'],policy=policy,both_complete=complete,
                baseline_status=a['status'],candidate_status=b['status'],
                exact_fasta_ledger_equal=all(a['score']['frozen_score'][key]==b['score']['frozen_score'][key] for key in ['output_sha256','ledger_sha256']) if complete else None))
    save_json(out/'paired-automatic-equality.json',pairs)
    save_json(out/'resource-use.json',dict(review_elapsed_seconds={role:sum(r['elapsed_seconds'] for r in reviews if r['role']==role) for role in roots},
        automatic_elapsed_seconds={role:sum(r['elapsed_seconds'] for r in autos if r['role']==role) for role in roots},
        note='Local descriptive subprocess elapsed sums; candidate scans perform additional work, no equal-workload speed claim'))
    print(json.dumps(dict(review_runs=len(reviews),complete_reviews=sum(r['status']=='complete' for r in reviews),
        automatic_runs=len(autos),complete_automatic=sum(r['status']=='complete' for r in autos),
        paired_successes=sum(r['both_complete'] for r in pairs),identical_successful_pairs=sum(r['exact_fasta_ledger_equal'] is True for r in pairs)),indent=2))


if __name__=='__main__':main()
