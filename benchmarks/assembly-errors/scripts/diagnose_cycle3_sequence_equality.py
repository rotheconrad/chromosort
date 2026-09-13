#!/usr/bin/env python3
"""Post-score D201 sequence comparison; never changes frozen provenance metrics."""
import argparse
import hashlib
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from benchio import digest,read_fasta,read_tsv,rc,save_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--panel-root',type=Path,required=True)
    p.add_argument('--results-root',type=Path,default=ROOT/'cycle3/results')
    p.add_argument('--actual-output',type=Path,help='Optional original committed replay FASTA for direct verification')
    p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    input_path=args.panel_root/'development/D201/inputs/assembly.fa'
    target_path=args.panel_root/'evaluator/development/D201/target.fa'
    ledger=args.results_root/'artifacts/candidate/D201/scripted_selection/components.tsv'
    selection=next(r for r in json.loads((args.results_root/'scripted-selection.json').read_text()) if r['role']=='candidate')
    assert selection['status']=='complete'
    assembly=read_fasta(input_path);target=read_fasta(target_path);outputs={}
    for row in read_tsv(ledger):
        assert row['input_id'] in assembly
        sequence=assembly[row['input_id']][int(row['input_start']):int(row['input_end'])]
        if row['strand']=='-':sequence=rc(sequence)
        previous=outputs.get(row['output_id'],'');assert len(previous)==int(row['output_start'])
        outputs[row['output_id']]=previous+sequence
    assert len(outputs)==len(target)==1
    observed=next(iter(outputs.values()));expected=next(iter(target.values()))
    if args.actual_output:
        assert digest(args.actual_output)==selection['score']['frozen_score']['output_sha256']
        assert read_fasta(args.actual_output)==outputs
        print('Original committed replay FASTA matches public-ledger reconstruction exactly')
    event=selection['selected_event'];s,e=70000,78000
    value=dict(schema='cycle3-post-score-sequence-diagnostic-v1',case_id='D201',role='candidate',
        classification='Requested after primary scores were known; supplementary interpretation only, not a revised success criterion or another candidate run',
        input_fasta_sha256=digest(input_path),target_fasta_sha256=digest(target_path),components_sha256=digest(ledger),
        recorded_original_output_fasta_sha256=selection['score']['frozen_score']['output_sha256'],
        output_sequence_sha256=hashlib.sha256(observed.encode('ascii')).hexdigest(),
        target_sequence_sha256=hashlib.sha256(expected.encode('ascii')).hexdigest(),
        output_bp=len(observed),target_bp=len(expected),exact_nucleotide_sequence_equal=observed==expected,
        differing_aligned_positions=sum(a!=b for a,b in zip(observed,expected)),length_difference_bp=len(observed)-len(expected),
        frozen_complete_structure=selection['score']['complete_structure'],
        frozen_target_relation_discrepancies=selection['score']['frozen_score']['target_relation_discrepancy_count'],
        declared_error_interval=[s,e],emitted_interval=[event['start'],event['end']],
        target_error_first_base=expected[s],target_error_last_base=expected[e-1],
        target_edge_bases_are_reverse_complements=expected[s]==rc(expected[e-1]),
        interpretation='Sequence equality and exact declared provenance are different outcomes. Equivalent edge bases leave exact target DNA after the emitted interior reversal while the frozen provenance path still fails. Oracle event selection and all original scores remain unchanged.')
    save_json(args.out,value);print(json.dumps(value,indent=2))


if __name__=='__main__':main()
