#!/usr/bin/env python3
"""Regenerate the frozen second development panel from bundled public sequences."""
import argparse
import gzip
import json
import random
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'src'), str(ROOT / 'scripts')]
from benchio import digest, read_fasta, write_fasta, rc, save_json, write_tsv
from generate import simulate_reads, MAP_FIELDS, REF_FIELDS, LABEL_FIELDS
from evidence import evidence


def generate_case(config, output_root, map_reads=True):
    cid = config['case_id']
    case = output_root / 'development' / cid
    truth_dir = output_root / 'evaluator/development' / cid
    if case.exists() or truth_dir.exists():
        raise ValueError('Refuse overwrite: ' + cid)
    case.mkdir(parents=True)
    truth_dir.mkdir(parents=True)
    def source(key):
        return next(iter(read_fasta(ROOT / 'cycle2/sources/sequences' / (key + '.fa')).values()))
    controls, events = [], []
    if config['family'] == 'mixed_exchange_misjoin':
        a, d = source('wheat_A_new'), source('wheat_D_new')
        start, end = config['exchange_interval']
        targets = {'t01': a[:start] + d[start:end] + a[end:],
                   't02': d[:start] + a[start:end] + d[end:]}
        references = [('CS_A', a, 'wheat_A', '1A'), ('CS_D', d, 'wheat_D', '1D')]
        parts = [[('t01', 0, len(a), '+'), ('t02', 0, len(d), '+')]]
        labels = {'t01': {'output_group': 'wheat_A', 'chromosome': '1A', 'subgenome': 'A'},
                  't02': {'output_group': 'wheat_D', 'chromosome': '1D', 'subgenome': 'D'}}
        for tid in targets:
            for j, point in enumerate([start, end]):
                controls.append(dict(variant_id=f'{tid}_exchange_edge{j + 1}',
                    kind='simulated_biological_reciprocal_exchange', target_id=tid,
                    start=point - 5000, end=point + 5000, empirical_event=False))
        events = [dict(event_id='E001', kind='interchromosomal_misjoin', input_id='c001',
            breakpoints=[len(a)], target_intervals=[['t01', len(a)-5000, len(a)], ['t02', 0, 5000]],
            historical_fidelity='documented_pattern_simulation_not_historical_edit')]
    else:
        if config['family'] == 'documented_scf139_orientation':
            native = source('wheat_v1_scf139')
            a, b = config['error_interval']
            target = native[:a] + rc(native[a:b]) + native[b:]
            guide = source('wheat_v2_scf139')
            fidelity = 'documented_orientation_reconstruction_with_gap_boundary_uncertainty'
        else:
            s, e = config['source_interval']
            guide = source(config['source'])[s:e]
            assert len(guide) == 300000 and set(guide) <= set('ACGT')
            s, e = config['biological_inversion']
            target = guide[:s] + rc(guide[s:e]) + guide[e:]
            controls.append(dict(variant_id='V001', kind='simulated_biological_inversion',
                target_id='t01', start=s - 2000, end=e + 2000, empirical_event=False))
            fidelity = 'documented_pattern_simulation_not_historical_edit'
        a, b = config['error_interval']
        targets = {'t01': target}
        references = [('CS_A', guide, 'wheat_A', '1A')]
        labels = {'t01': {'output_group': 'wheat_A', 'chromosome': '1A', 'subgenome': 'A'}}
        parts = [[('t01', 0, a, '+'), ('t01', a, b, '-'), ('t01', b, len(target), '+')]]
        events = [dict(event_id='E001', kind='orientation', input_id='c001', breakpoints=[a, b],
                      target_intervals=[['t01', a, b]], historical_fidelity=fidelity)]
    assembly, mapping = {}, []
    for i, group in enumerate(parts, 1):
        name, position, chunks = f'c{i:03d}', 0, []
        for tid, start, end, strand in group:
            sequence = targets[tid][start:end]
            chunks.append(rc(sequence) if strand == '-' else sequence)
            mapping.append(dict(input_id=name, input_start=position, input_end=position+len(sequence),
                target_id=tid, target_start=start, target_end=end, strand=strand))
            position += len(sequence)
        assembly[name] = ''.join(chunks)
    write_fasta(case / 'inputs/assembly.fa', assembly)
    write_fasta(truth_dir / 'target.fa', targets)
    write_tsv(truth_dir / 'input_to_target.tsv', mapping, MAP_FIELDS)
    refrows, labelrows = [], []
    for rid, sequence, group, chrom in references:
        filename = f'inputs/ref_{rid}.fa'
        write_fasta(case / filename, {'Chr1': sequence})
        refrows.append(dict(reference_id=rid, fasta=filename, sha256=digest(case / filename),
            assembly_accession='GCA_018294505.1', output_group=group, role='backbone'))
        labelrows.append(dict(reference_id=rid, sequence_id='Chr1', chromosome=chrom,
            subgenome=chrom[-1], homoeolog_group='', copy_label='', output_group=group, backbone='yes'))
    write_tsv(case / 'inputs/references.tsv', refrows, REF_FIELDS)
    write_tsv(case / 'inputs/reference_sequences.tsv', labelrows, LABEL_FIELDS)
    fastq = case / 'inputs/reads.fastq.gz'
    read_model = simulate_reads(targets, fastq, truth_dir, random.Random(config['seed']+1000003),
        config['coverage'], config['read_length'], config['read_error'])
    condition = config['evidence_condition']
    if condition.startswith('one_molecule'):
        point = config['misjoin_position']
        chimera = assembly['c001'][point - 5000:point + 5000]
        assert len(chimera) == 10000 and set(chimera) <= set('ACGT')
        with open(fastq, 'ab') as raw, gzip.GzipFile(fileobj=raw, mode='wb', filename='', mtime=0) as stream:
            stream.write(('@read9000001\n'+chimera+'\n+\n'+'?'*10000+'\n').encode())
        save_json(truth_dir / 'adversarial_read_origin.json', dict(parent_molecule='read9000001',
            source='corrupted_input_not_biological_target', input_id='c001', start=point-5000,
            end=point+5000, condition=condition, molecular_origin_count=1,
            exposed_names=['read9000001'] if condition.endswith('20_rows') else ['read9000001','read9000002']))
        read_model['added_artificial_chimeric_molecules'] = 1
    manifest = dict(schema='chromosort-input-v1', case_id=cid, split='development',
        assembly=dict(path='inputs/assembly.fa', sha256=digest(case / 'inputs/assembly.fa'), sample_id=cid, stage='input'),
        references='inputs/references.tsv', reference_sequences='inputs/reference_sequences.tsv', evidence=[],
        readsets=[dict(readset_id='sim01', path='inputs/reads.fastq.gz', sha256=digest(fastq),
            evidence_class='mapped_simulated_reads', dependency_group='sim01')], coordinate_system='0-based-half-open')
    save_json(case / 'manifest.json', manifest)
    save_json(truth_dir / 'truth.json', dict(schema='benchmark-truth-v1', case_id=cid, config=config,
        assembly_sha256=manifest['assembly']['sha256'], target_sha256=digest(truth_dir / 'target.fa'),
        mapping='input_to_target.tsv', events=events, preservation_controls=controls, target_labels=labels,
        read_model=read_model, breakpoint_tolerance_bp=500, correct_endpoint_anchor_bp=100,
        source_constructions=[dict(construction='See frozen cycle2/cases.json and source admission')],
        baseline_mask='All non-ACGT target positions', limitations=['Development mechanism panel; not independent genomes or empirical raw reads']))
    save_json(truth_dir / 'config.json', config)
    if map_reads:
        evidence(case, threads=2)
        manifest = json.loads((case / 'manifest.json').read_text())
        entry = next(e for e in manifest['evidence'] if e['evidence_id'] == 'reads_sim01')
        paf = case / entry['path']
        mapped_sha = digest(paf)
        rows = [line.split('\t') for line in paf.read_text().splitlines()]
        transform = 'identity'
        if condition == 'missing_cigar':
            rows = [[f for f in row if not f.startswith('cg:Z:')] for row in rows]
            transform = 'remove all cg:Z CIGAR tags'
        elif condition.startswith('one_molecule'):
            chimera_rows = [r for r in rows if r[0] == 'read9000001']
            assert len(chimera_rows) == 1, 'Prespecified single chimera did not map as one row; report construction failure'
            rows = [r for r in rows if r[0] != 'read9000001']
            if condition.endswith('20_rows'):
                rows += [list(chimera_rows[0]) for _ in range(20)]
                transform = 'repeat one chimeric-molecule PAF row 20 times, same read ID'
            else:
                rows += [[name] + chimera_rows[0][1:] for name in ['read9000001', 'read9000002']]
                transform = 'two copies of one chimeric-molecule PAF row, renamed IDs; raw FASTQ contains its one parent'
        paf.write_text(''.join('\t'.join(r)+'\n' for r in rows))
        entry.update(sha256=digest(paf), derivation=dict(mapped_paf_sha256=mapped_sha,
            transformation=transform, condition=condition, script='scripts/generate_cycle2.py'))
        save_json(case / 'manifest.json', manifest)
    print(cid, sum(map(len, assembly.values())), 'input bp')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-root', type=Path, required=True)
    p.add_argument('--cases', nargs='+')
    p.add_argument('--without-mapping', action='store_true')
    args = p.parse_args()
    out = args.output_root.resolve()
    if out == ROOT or (ROOT in out.parents and out.relative_to(ROOT).parts[0] != 'work'):
        p.error('Use work/ or a directory outside the benchmark')
    config = json.loads((ROOT / 'cycle2/cases.json').read_text())
    selected = args.cases or [c['case_id'] for c in config['cases']]
    if set(selected) - {c['case_id'] for c in config['cases']}:
        p.error('Unknown case')
    for item in json.loads((ROOT / 'cycle2/sources/sequence_registry.json').read_text())['sequences']:
        assert digest(ROOT / 'cycle2' / item['path']) == item['sha256']
    for case in config['cases']:
        if case['case_id'] in selected:
            generate_case(case, out, not args.without_mapping)


if __name__ == '__main__':
    main()
