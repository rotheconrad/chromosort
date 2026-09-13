#!/usr/bin/env python3
"""Run fixed cycle3 automatic arms and review scans, committing outputs before truth."""
import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'scripts')]
from benchio import digest,save_json
from run_benchmark import chromosort_identity,input_hashes
from normalize_cycle3_events import normalize


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--python',type=Path,required=True)
    p.add_argument('--software-id',required=True)
    p.add_argument('--role',choices=['baseline','candidate'],required=True)
    p.add_argument('--panel-root',type=Path,required=True)
    p.add_argument('--known-root',type=Path,required=True)
    p.add_argument('--output-root',type=Path,required=True)
    args=p.parse_args();out=args.output_root.resolve()
    if out.exists():p.error('Choose a fresh output directory')
    if out==ROOT or (ROOT in out.parents and out.relative_to(ROOT).parts[0]!='work'):p.error('Use work/ or an external output directory')
    protocol=json.loads((ROOT/'cycle3/protocol.json').read_text())
    env=dict(os.environ);env.pop('PYTHONPATH',None);env['PYTHONDONTWRITEBYTECODE']='1'
    python=args.python.absolute()
    identity=chromosort_identity(python,env);identity.update(role=args.role,declared_identity=args.software_id)
    out.mkdir(parents=True);records=[]
    save_json(out/'software.json',identity)
    for cohort,panel,ids in [('primary',args.panel_root.resolve(),protocol['primary_cases']),('known',args.known_root.resolve(),protocol['known_diagnostics'])]:
        automatic=[sys.executable,str(ROOT/'scripts/run_benchmark.py'),'--tool','chromosort','--python',str(python),
            '--software-id',args.software_id,'--case-root',str(panel/'development'),'--output-root',str(out/'automatic'/cohort),
            '--cases',*ids,'--modes','comprehensive','--policies','report','review','--defer-scoring']
        result=subprocess.run(automatic,env=env)
        if result.returncode:print('Automatic failures retained in',cohort,flush=True)
        for cid in ids:
            case=panel/'development'/cid;directory=out/'review'/cid;directory.mkdir(parents=True)
            manifest=json.loads((case/'manifest.json').read_text());hashes=input_hashes(case,manifest)
            scan=directory/'scan'
            replacement={'CASE/manifest.json':str(case/'manifest.json'),'OUTPUT/scan':str(scan)}
            command=[str(python),'-m','chromosort.cli',*[replacement.get(x,x) for x in protocol['commands']['scan_common']]]
            if args.role=='candidate':command+=protocol['commands']['candidate_scan_additions']
            start=time.monotonic()
            with (directory/'scan.log').open('w') as log:
                process=subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT)
            record=dict(case_id=cid,cohort=cohort,role=args.role,input_files_sha256=hashes,
                command=['selected-python',*[s.replace(str(case),'CASE').replace(str(directory),'OUTPUT') for s in command[1:]]],
                elapsed_seconds=time.monotonic()-start,returncode=process.returncode,status='complete' if process.returncode==0 else 'failed')
            if process.returncode==0:
                try:
                    value=normalize(case,scan)
                    save_json(directory/'normalized-events.json',value)
                    record['normalized_sha256']=digest(directory/'normalized-events.json')
                except Exception as error:
                    record.update(status='invalid_event_contract',error=str(error))
            save_json(directory/'run.json',record);records.append(record)
            print(args.role,cid,'scan',record['status'],flush=True)
    save_json(out/'execution.json',dict(schema='cycle3-output-execution-v1',software=identity,
        protocol_sha256=digest(ROOT/'cycle3/protocol.json'),runs=records,
        scorer_access='No truth opened or scorer imported; await both output commitments'))
    save_json(out/'outputs.commitment.json',{str(p.relative_to(out)):digest(p) for p in sorted(out.rglob('*')) if p.is_file()})


if __name__=='__main__':main()
