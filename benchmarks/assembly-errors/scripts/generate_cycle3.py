#!/usr/bin/env python3
"""Generate the fixed cycle3 visibility scenarios; no curation software imported."""
import argparse
import json
import random
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'src'), str(ROOT / 'scripts')]
from benchio import digest, read_fasta, rc, save_json, write_fasta, write_tsv
from generate import simulate_reads, MAP_FIELDS, REF_FIELDS, LABEL_FIELDS
from evidence import evidence


def generate(config, out):
    cid = config['case_id']
    case, truth_dir = out / 'development' / cid, out / 'evaluator/development' / cid
    if case.exists() or truth_dir.exists():
        raise ValueError('Refuse overwrite: ' + cid)
    case.mkdir(parents=True)
    truth_dir.mkdir(parents=True)
    def source(key):
        return next(iter(read_fasta(ROOT / 'cycle2/sources/sequences' / (key+'.fa')).values()))
    controls, events, native_controls = [], [], []
    if config['family'] == 'mixed_exchange_misjoin':
        a, d = [source(key) for key in config['source_ids']]
        s, e = config['exchange_interval']
        targets = {'t01': a[:s]+d[s:e]+a[e:], 't02': d[:s]+a[s:e]+d[e:]}
        references = [('CS_A', a, 'wheat_A', '1A'), ('CS_D', d, 'wheat_D', '1D')]
        parts = [('t01', 0, len(a), '+'), ('t02', 0, len(d), '+')]
        for offset, tid in [(0, 't01'), (len(a), 't02')]:
            for j, point in enumerate([s, e]):
                vid = f'{tid}_exchange_{j}'
                controls.append(dict(variant_id=vid, kind='simulated_biological_exchange', target_id=tid,
                    start=point-2000, end=point+2000, empirical_event=False))
                native_controls.append(dict(control_id=vid, kind='biological', contig='c001', start=offset+point-2000, end=offset+point+2000))
        events = [dict(event_id='E001', kind='interchromosomal_misjoin', input_id='c001',
            breakpoints=[len(a)], target_intervals=[['t01', len(a)-2000, len(a)], ['t02', 0, 2000]],
            historical_fidelity='documented_pattern_simulation')]
    else:
        s, e = config['source_interval']
        guide = source(config['source_id'])[s:e]
        assert len(guide) == 300000 and set(guide) <= set('ACGT')
        target = guide
        if config.get('biological_inversion'):
            s, e = config['biological_inversion']
            target = target[:s]+rc(target[s:e])+target[e:]
            controls = [dict(variant_id='V001', kind='simulated_biological_inversion', target_id='t01',
                start=s-2000, end=e+2000, empirical_event=False)]
            native_controls = [dict(control_id='V001', kind='biological', contig='c001', start=s-2000, end=e+2000)]
        targets = {'t01': target}
        references = [('CS_A', guide, 'wheat_A', '1A')]
        if config.get('second_reference'):
            references.append(('CS_A_reverse_frame', rc(guide), 'wheat_A', '1A'))
        parts = [('t01', 0, len(target), '-' if config['family']=='whole_record_reverse' else '+')]
        if config.get('error_interval'):
            s, e = config['error_interval']
            parts = [('t01', 0, s, '+'), ('t01', s, e, '-'), ('t01', e, len(target), '+')]
            events = [dict(event_id='E001', kind='orientation', input_id='c001', breakpoints=[s, e],
                target_intervals=[['t01', s, e]], historical_fidelity='documented_pattern_simulation')]
    chunks, mapping, position = [], [], 0
    for tid, s, e, strand in parts:
        seq = targets[tid][s:e]
        chunks.append(rc(seq) if strand == '-' else seq)
        mapping.append(dict(input_id='c001', input_start=position, input_end=position+e-s,
            target_id=tid, target_start=s, target_end=e, strand=strand))
        position += e-s
    assembly = {'c001': ''.join(chunks)}
    write_fasta(case/'inputs/assembly.fa', assembly)
    write_fasta(truth_dir/'target.fa', targets)
    write_tsv(truth_dir/'input_to_target.tsv', mapping, MAP_FIELDS)
    refrows, labelrows = [], []
    for rid, seq, group, chrom in references:
        name = f'inputs/ref_{rid}.fa'
        write_fasta(case/name, {'Chr1': seq})
        refrows.append(dict(reference_id=rid, fasta=name, sha256=digest(case/name),
            assembly_accession='GCA_018294505.1', output_group=group, role='backbone'))
        labelrows.append(dict(reference_id=rid, sequence_id='Chr1', chromosome=chrom,
            subgenome=chrom[-1], homoeolog_group='', copy_label='', output_group=group, backbone='yes'))
    write_tsv(case/'inputs/references.tsv', refrows, REF_FIELDS)
    write_tsv(case/'inputs/reference_sequences.tsv', labelrows, LABEL_FIELDS)
    fastq = case/'inputs/reads.fastq.gz'
    model = simulate_reads(targets, fastq, truth_dir, random.Random(config['seed']+1000003),
        config['coverage'], config['read_length'], config['read_error'])
    manifest = dict(schema='chromosort-input-v1', case_id=cid, split='development',
        assembly=dict(path='inputs/assembly.fa', sha256=digest(case/'inputs/assembly.fa'), sample_id=cid, stage='input'),
        references='inputs/references.tsv', reference_sequences='inputs/reference_sequences.tsv', evidence=[],
        readsets=[dict(readset_id='sim01', path='inputs/reads.fastq.gz', sha256=digest(fastq),
            evidence_class='mapped_simulated_reads', dependency_group='sim01')], coordinate_system='0-based-half-open')
    save_json(case/'manifest.json', manifest)
    labels = {tid:dict(output_group='wheat_A' if tid=='t01' else 'wheat_D', chromosome='1A' if tid=='t01' else '1D',subgenome='A' if tid=='t01' else 'D') for tid in targets}
    save_json(truth_dir/'truth.json', dict(schema='benchmark-truth-v1', case_id=cid, config=config,
        assembly_sha256=manifest['assembly']['sha256'], target_sha256=digest(truth_dir/'target.fa'),
        mapping='input_to_target.tsv', events=events, preservation_controls=controls, target_labels=labels,
        read_model=model, breakpoint_tolerance_bp=500, correct_endpoint_anchor_bp=100))
    save_json(truth_dir/'review_truth.json', dict(schema='cycle3-review-truth-v1',case_id=cid,
        assembly_sha256=manifest['assembly']['sha256'],contig_lengths={'c001':position},
        errors=[dict(error_id=e['event_id'],contig=e['input_id'],shape='interval' if e['kind']=='orientation' else 'point',
            start=e['breakpoints'][0],end=e['breakpoints'][-1],kind=e['kind']) for e in events],
        biological_controls=native_controls,whole_record_orientation_control=config['family']=='whole_record_reverse',
        clean_control=config['family']=='clean',scope='Known simulation geometry; inspection burden is not a false correction'))
    save_json(truth_dir/'config.json',config)
    evidence(case,threads=2)
    if config.get('assembly_paf_derivation'):
        manifest=json.loads((case/'manifest.json').read_text())
        entry=next(e for e in manifest['evidence'] if e['kind']=='assembly_paf')
        paf=case/entry['path'];original=paf.read_text().splitlines();raw_sha=digest(paf)
        matches=[r.split('\t') for r in original if r.split('\t')[4]=='-' and int(r.split('\t')[2])<=70500 and int(r.split('\t')[3])>=77500]
        assert len(matches)==1 and matches[0][2:4]==matches[0][7:9]
        assert 'cg:Z:'+str(int(matches[0][3])-int(matches[0][2]))+'M' in matches[0], 'Fixed nested derivation requires one exact covering reverse row'
        nested=matches[0][:12]
        nested[2:4]=['70500','77500'];nested[7:9]=['70500','77500'];nested[9:11]=['7000','7000']
        nested+=['tp:A:P','cg:Z:7000M']
        paf.write_text('\n'.join(original*3+['\t'.join(nested)])+'\n')
        entry.update(sha256=digest(paf),derivation=dict(mapped_paf_sha256=raw_sha,
            description=config['assembly_paf_derivation'],dependency='same original mapping; no independent source'))
        # Remove convenience union rather than leave it inconsistent with the selected derived PAF.
        for key in ['combined_reference_path','combined_paf_path']:
            (case/manifest['convenience'][key]).unlink()
        manifest.pop('convenience')
        save_json(case/'manifest.json',manifest)
    print(cid,position,'input bp',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-root',type=Path,required=True)
    args=parser.parse_args();out=args.output_root.resolve()
    if out==ROOT or (ROOT in out.parents and out.relative_to(ROOT).parts[0]!='work'):
        parser.error('Use work/ or an external directory')
    config=json.loads((ROOT/'cycle3/cases.json').read_text())
    assert digest(ROOT/'cycle3'/config['source_registry'])==config['source_registry_sha256']
    for source in json.loads((ROOT/'cycle3'/config['source_registry']).read_text())['sequences']:
        assert digest(ROOT/'cycle2'/source['path'])==source['sha256']
    for case in config['cases']:generate(case,out)
