#!/usr/bin/env python3
"""Regenerate and verify cycle3 using an isolated copy of public benchmark content."""
import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path)
    args=p.parse_args();records=[]
    with tempfile.TemporaryDirectory(prefix='chromosort-cycle3-public-') as temporary:
        root=Path(temporary)/'assembly-errors'
        def ignore(path,names):
            return [name for name in names if name in ['work','__pycache__','.ruff_cache','.DS_Store'] or name.endswith(('.pyc','.fai'))]
        shutil.copytree(ROOT,root,ignore=ignore)
        commands=[['-m','unittest','discover','-s','tests'],['scripts/generate_cycle3.py','--output-root','work/panel'],
            ['scripts/verify_cycle3_inputs.py','--generated-root','work/panel'],
            ['scripts/prepare_cycle3_known.py','--output-root','work/known'],
            ['scripts/verify_cycle3_results.py','--panel-root','work/panel','--known-root','work/known'],
            ['scripts/verify_cycle3_context.py','--panel-root','work/panel','--known-root','work/known'],
            ['scripts/summarize_cycle3.py'],
            ['scripts/diagnose_cycle3_sequence_equality.py','--panel-root','work/panel','--out','work/post-score-sequence-diagnostic.json'],
            ['scripts/summarize_history.py'],['scripts/verify.py','--fixtures']]
        summary_names=['cohort-summary.json','candidate-events.tsv','candidate-bounds.tsv','candidate-coordinate-evidence.tsv']
        expected_summaries={name:(root/'cycle3/results'/name).read_bytes() for name in summary_names}
        for command in commands:
            result=subprocess.run([sys.executable,*command],cwd=root,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            records.append(dict(command=['python',*command],returncode=result.returncode));print(result.stdout,end='',flush=True)
            if result.returncode:raise SystemExit(result.returncode)
        for name,raw in expected_summaries.items():assert (root/'cycle3/results'/name).read_bytes()==raw,name
        assert (root/'work/post-score-sequence-diagnostic.json').read_bytes()==(root/'cycle3/post-score-sequence-diagnostic.json').read_bytes()
    record=dict(schema='cycle3-clean-copy-verification-v1',status='pass',isolated_public_copy=True,network_used=False,
        python=sys.version.split()[0],minimap2=subprocess.check_output(['minimap2','--version'],text=True).strip(),commands=records,
        scope='Regeneration and portable metric verification; no installed ChromoSort, software-task outputs or private source files')
    if args.out:args.out.write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
    print('Cycle3 isolated public verification passed')


if __name__=='__main__':main()
