#!/usr/bin/env python3
"""Regenerate three known diagnostics and derive native review-region truth."""
import argparse
import json
import subprocess
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from benchio import digest,read_fasta,read_tsv,save_json


def derive(root):
    for cid in ['D101','D012','D108']:
        case=root/'development'/cid;directory=root/'evaluator/development'/cid
        truth=json.loads((directory/'truth.json').read_text())
        mappings=read_tsv(directory/truth['mapping']);controls=[]
        for control in truth['preservation_controls']:
            by_contig={}
            for row in mappings:
                if row['target_id']!=control['target_id']:continue
                s,e=max(control['start'],int(row['target_start'])),min(control['end'],int(row['target_end']))
                if s>=e:continue
                q=int(row['input_start']);t=int(row['target_start']);te=int(row['target_end'])
                a,b=(q+s-t,q+e-t) if row['strand']=='+' else (q+te-e,q+te-s)
                by_contig.setdefault(row['input_id'],[]).append((a,b))
            for contig,spans in by_contig.items():
                controls.append(dict(control_id=control['variant_id']+'_'+contig,kind='preserved_biology_or_empirical_copy',contig=contig,
                    start=min(s for s,e in spans),end=max(e for s,e in spans)))
        errors=[dict(error_id=e['event_id'],contig=e['input_id'],shape='interval' if e['kind']=='orientation' else 'point',
            start=e['breakpoints'][0],end=e['breakpoints'][-1],kind=e['kind']) for e in truth['events']]
        save_json(directory/'review_truth.json',dict(schema='cycle3-review-truth-v1',case_id=cid,
            assembly_sha256=truth['assembly_sha256'],original_truth_sha256=digest(directory/'truth.json'),
            contig_lengths={k:len(v) for k,v in read_fasta(case/'inputs/assembly.fa').items()},errors=errors,
            biological_controls=controls,scope='Known development diagnostic, excluded from prospective primary cohort'))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-root',type=Path,required=True)
    parser.add_argument('--derive-only',action='store_true')
    args=parser.parse_args();out=args.output_root.resolve()
    if out==ROOT or (ROOT in out.parents and out.relative_to(ROOT).parts[0]!='work'):
        parser.error('Use work/ or an external directory')
    if not args.derive_only:
        if out.exists():parser.error('Choose a new root or derive only from pre-existing known diagnostics')
        commands=[['retrieve.py','--output-root',out],['generate.py','--output-root',out,'--cases','D012'],
            ['evidence.py',out/'development/D012','--threads','2'],
            ['generate_cycle2.py','--output-root',out,'--cases','D101','D108']]
        for script,*options in commands:
            subprocess.run([sys.executable,str(ROOT/'scripts'/script),*map(str,options)],check=True)
    derive(out)
    print('Derived native review truth for three known diagnostics')
