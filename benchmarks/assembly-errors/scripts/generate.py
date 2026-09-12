#!/usr/bin/env python3
"""Create public-sequence cases and separate evaluator truth. No tool logic."""
import argparse
import gzip
import json
import math
import random
import re
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from benchio import digest, read_fasta, write_fasta, rc, save_json, write_tsv

ROOT = Path(__file__).resolve().parents[1]
MAP_FIELDS = ['input_id', 'input_start', 'input_end', 'target_id', 'target_start', 'target_end', 'strand']
REF_FIELDS = ['reference_id', 'fasta', 'sha256', 'assembly_accession', 'output_group', 'role']
LABEL_FIELDS = ['reference_id', 'sequence_id', 'chromosome', 'subgenome', 'homoeolog_group', 'copy_label', 'output_group', 'backbone']


def read_source(key, source_dir):
    return next(iter(read_fasta(Path(source_dir) / f'{key}.fa').values()))


def mutate(seq, rate, rng):
    seq = list(seq)
    for i in rng.sample(range(len(seq)), round(len(seq) * rate)):
        if seq[i] in 'ACGT':
            seq[i] = rng.choice([b for b in 'ACGT' if b != seq[i]])
    return ''.join(seq)


def simulate_reads(targets, out, truth_dir, rng, coverage, read_length, error_rate):
    """Uniform starts within eligible ACGT runs; no reads across baseline Ns.

    This is a transparent i.i.d. substitution mechanism model, not an empirical
    PacBio/ONT error profile. Provenance stays evaluator-side; opaque read IDs.
    """
    origins, rows = [], []
    runs = []
    for tid, seq in targets.items():
        for m in re.finditer('[ACGT]+', seq):
            if m.end() - m.start() >= read_length:
                runs.append((tid, m.start(), m.end()))
    eligible = sum(e - s for _, s, e in runs)
    nreads = math.ceil(coverage * eligible / read_length)
    weights = [e - s - read_length + 1 for _, s, e in runs]
    if not runs:
        raise ValueError('No interval can generate the requested read length')
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, 'wb') as raw, gzip.GzipFile(fileobj=raw, mode='wb', filename='', mtime=0) as f:
        for i in range(nreads):
            tid, s, e = rng.choices(runs, weights=weights, k=1)[0]
            start = rng.randint(s, e - read_length)
            seq = targets[tid][start:start + read_length]
            strand = rng.choice('+-')
            if strand == '-':
                seq = rc(seq)
            seq = mutate(seq, error_rate, rng)
            name = f'read{i + 1:07d}'
            f.write(f'@{name}\n{seq}\n+\n'.encode() + b'?' * len(seq) + b'\n')
            origins.append(dict(read_id=name, target_id=tid, start=start, end=start + read_length, strand=strand))
    write_tsv(truth_dir / 'read_origins.tsv', origins, ['read_id', 'target_id', 'start', 'end', 'strand'])
    return dict(model='iid_substitution_v1', coverage_requested=coverage,
                coverage_denominator='ACGT runs at least one read long',
                read_count=nreads, read_length=read_length, eligible_bp=eligible,
                total_target_bp=sum(map(len, targets.values())), substitution_rate=error_rate,
                indel_rate=0, quality_char='?', quality_semantics='constant Q30 model value, not measured',
                source='biological_target_before_corruption', strand_probability=0.5)


def generate(config, output_root, source_dir=None):
    output_root = Path(output_root).resolve()
    if output_root == ROOT or ROOT / 'fixtures' in [output_root, *output_root.parents]:
        raise ValueError('Use a work directory outside tracked fixtures')
    source_dir = Path(source_dir) if source_dir else output_root / 'sources/sequences'
    def source(key):
        return read_source(key, source_dir)
    cid, family, seed = config['case_id'], config['family'], config['seed']
    split = config.get('split', 'development')
    if split != 'development' or not re.fullmatch(r'D[0-9]{3}', cid):
        raise ValueError('This public generator accepts development cases only')
    case = output_root / split / cid
    truth_dir = output_root / 'evaluator' / split / cid
    if (case / 'manifest.json').exists():
        raise ValueError('Refusing to overwrite an existing case: ' + cid)
    rng = random.Random(seed)
    targets, refs, labels, variants, events, parts, source_maps = {}, [], {}, [], [], [], []
    def add_target(tid, seq, key, group, chrom, sub='', note='unmodified public sequence'):
        targets[tid] = seq
        labels[tid] = dict(output_group=group, chromosome=chrom, subgenome=sub)
        source_maps.append(dict(target_id=tid, source_id=key, construction=note))
    def add_ref(rid, seqs, group, chrom, sub='', accession='', backbone=True):
        refs.append(dict(reference_id=rid, sequences=seqs, output_group=group,
                         chromosome=chrom, subgenome=sub, accession=accession, backbone=backbone))
    def piece(tid, start, end, strand='+'):
        return (tid, start, end, strand)
    def whole(tid):
        return [piece(tid, 0, len(targets[tid]))]

    if family == 'wheat_documented_orientation':
        native = source('wheat_v1_1b')
        a, b = 250000, 393811
        target = native[:a] + rc(native[a:b]) + native[b:]
        add_target('t01', target, 'wheat_v1_1b', 'wheat_B', '1B', 'B',
                   'Reverse-complement the Table S3 scaffold198 reported interval; retain v1 sequence and gap Ns. Reported inclusive endpoints define operational truth despite the 1 bp table-length discrepancy; only non-ACGT bases are masked, not the adjacent called base.')
        add_ref('CS_v2', {'Chr1B': source('wheat_v2_1b')}, 'wheat_B', '1B', 'B', 'GCA_018294505.1')
        parts = [[piece('t01', 0, a), piece('t01', a, b, '-'), piece('t01', b, len(target))]]
        events = [dict(event_id='E001', kind='orientation', source_id='wheat_scf198',
                       historical_fidelity='documented_orientation_reconstruction_with_gap_boundary_uncertainty',
                       input_id='c001', breakpoints=[a, b], target_intervals=[['t01', a, b]])]
    elif family.startswith('multi_'):
        public = [('A', 'wheat_v2_1a', '1A'), ('B', 'wheat_v2_1b', '1B'), ('D', 'wheat_v2_1d', '1D')]
        for sub, key, chrom in public:
            tid = 't0' + str(len(targets) + 1)
            seq = source(key)[:250000]
            add_target(tid, seq, key, 'wheat_' + sub, chrom, sub,
                       'First 250000 bp of registered interval; loci are chromosome-labeled, not asserted syntenic homoeologous intervals')
            add_ref('CS_' + sub, {'Chr1': seq}, 'wheat_' + sub, chrom, sub, 'GCA_018294505.1')
        if family == 'multi_exchange':
            # Simulated reciprocal exchange is target biology, not an assembly error.
            a, b = targets['t01'], targets['t02']
            targets['t01'], targets['t02'] = a[:125000] + b[125000:], b[:125000] + a[125000:]
            for tid in ['t01', 't02']:
                variants.append(dict(variant_id='V' + tid, kind='simulated_reciprocal_exchange', target_id=tid,
                                     start=120000, end=130000, source_id='design_preservation_control',
                                     empirical_event=False))
                source_maps.append(dict(target_id=tid, construction='Reciprocal A/B target-tail exchange at 125000; chromosome membership remains supplied; ancestry changes.'))
        parts = [whole(tid) for tid in targets]
        if family == 'multi_ambiguous':
            # An alternative label for exactly the same sequence is an ambiguity
            # fixture. It is not evidence that B and D are biologically identical.
            add_ref('label_conflict', {'Chr1': targets['t02']}, 'candidate_other', '1B', '',
                    'GCA_018294505.1:label_conflict_fixture', True)
            labels['t02']['allowed_groups'] = ['wheat_B', 'candidate_other']
            labels['t02']['require_ambiguous'] = True
        if family == 'multi_misjoin':
            parts = [whole('t01') + whole('t03'), whole('t02')]
            events = [dict(event_id='E001', kind='interchromosomal_misjoin', source_id='wheat_chimeras',
                           historical_fidelity='documented_pattern_simulation', input_id='c001',
                           breakpoints=[250000], target_intervals=[['t01', 245000, 250000], ['t03', 0, 5000]])]
    else:
        native = source('barley_mla')
        target = mutate(native, config.get('target_divergence', 0.0), rng)
        # A target allele independent of the corrupted copy. The actual tandem
        # duplication in the public BAC remains an empirical preservation control.
        if config.get('biological_inversion', True):
            start, end = 190000, 205000
            target = target[:start] + rc(target[start:end]) + target[end:]
            variants.append(dict(variant_id='V001', kind='simulated_biological_inversion', target_id='t01',
                                 start=start - 2000, end=end + 2000, empirical_event=False,
                                 source_id='design_preservation_control'))
        add_target('t01', target, 'barley_mla', 'barley', '1H', note='BAC plus declared target-allele simulation, before assembly corruption')
        variants.extend([dict(variant_id='V002', kind='empirical_tandem_copy_1', target_id='t01',
                              start=59567, end=99282, empirical_event=True, source_id='barley_mla'),
                         dict(variant_id='V003', kind='empirical_tandem_copy_2', target_id='t01',
                              start=99282, end=138992, empirical_event=True, source_id='barley_mla')])
        guide = native
        if family == 'reference_gap':
            guide = native[:150000] + native[155000:]
            variants.append(dict(variant_id='V004', kind='simulated_target_specific_sequence', target_id='t01',
                                 start=150000, end=155000, empirical_event=False, source_id='design_preservation_control'))
        if family == 'reference_overlap_only':
            guide = native[:59567] + native[99282:]
        if family == 'novel_retention':
            guide = native[:170000] + native[185000:]
            variants.append(dict(variant_id='V004', kind='simulated_target_specific_sequence', target_id='t01',
                                 start=170000, end=185000, empirical_event=False, source_id='design_preservation_control'))
        add_ref('Mla_BAC', {'Chr1': guide}, 'barley', '1H', accession='AF427791.1')
        parts = [whole('t01')]
        a = config.get('error_start', 151000)
        b = a + config.get('error_size', 8000)
        if family in ['barley_orientation', 'barley_two_errors']:
            assert 139000 <= a < b < 189000
            parts = [[piece('t01', 0, a), piece('t01', a, b, '-'), piece('t01', b, len(target))]]
            events = [dict(event_id='E001', kind='orientation', source_id='barley_mla_error',
                           historical_fidelity='documented_pattern_simulation_not_exact_historical_edit',
                           input_id='c001', breakpoints=[a, b], target_intervals=[['t01', a, b]])]
        if family in ['barley_collapse', 'barley_two_errors']:
            c, d = 99282, 138992
            if family == 'barley_two_errors':
                parts[0] = [piece('t01', 0, c), piece('t01', d, a)] + parts[0][1:]
                events[0]['breakpoints'] = [a - (d - c), b - (d - c)]
            else:
                parts = [[piece('t01', 0, c), piece('t01', d, len(target))]]
            events.append(dict(event_id='E002', kind='collapse_missing_sequence', source_id='barley_mla_error',
                               historical_fidelity='documented_pattern_simulation_on_empirical_repeat',
                               input_id='c001', breakpoints=[c], target_intervals=[['t01', c, d]],
                               missing_bp=d - c, expected_capability='unresolved_without_sequence_source'))
        if family == 'true_overlap':
            parts = [[piece('t01', 0, 150000)], [piece('t01', 140000, len(target))]]
        if family in ['repeat_overlap', 'reference_overlap_only']:
            parts = [[piece('t01', 0, 99282)], [piece('t01', 99282, len(target))]]
        if family == 'reference_gap':
            parts = [[piece('t01', 0, 150000)], [piece('t01', 155000, len(target))]]
        if family == 'oversegmentation':
            # Clean input, repeated alignment fragments produced by real repeats.
            parts = [whole('t01')]
        if family == 'novel_retention':
            parts = [[piece('t01',0,170000)], [piece('t01',185000,len(target))], [piece('t01',170000,185000)]]
    case.mkdir(parents=True, exist_ok=True)
    truth_dir.mkdir(parents=True, exist_ok=True)
    assembly, mapping = {}, []
    for i, group in enumerate(parts, 1):
        name, chunks, pos = f'c{i:03d}', [], 0
        for tid, start, end, strand in group:
            seq = targets[tid][start:end]
            if strand == '-':
                seq = rc(seq)
            chunks.append(seq)
            mapping.append(dict(input_id=name, input_start=pos, input_end=pos + len(seq),
                                target_id=tid, target_start=start, target_end=end, strand=strand))
            pos += len(seq)
        assembly[name] = ''.join(chunks)
    write_fasta(case / 'inputs/assembly.fa', assembly)
    write_fasta(truth_dir / 'target.fa', targets)
    write_tsv(truth_dir / 'input_to_target.tsv', mapping, MAP_FIELDS)
    refrows, labelrows = [], []
    for r in refs:
        filename = 'inputs/ref_' + r['reference_id'] + '.fa'
        write_fasta(case / filename, r['sequences'])
        refrows.append(dict(reference_id=r['reference_id'], fasta=filename, sha256=digest(case / filename),
                            assembly_accession=r['accession'], output_group=r['output_group'], role='backbone'))
        for name in r['sequences']:
            labelrows.append(dict(reference_id=r['reference_id'], sequence_id=name, chromosome=r['chromosome'],
                                  subgenome=r['subgenome'], homoeolog_group='', copy_label='',
                                  output_group=r['output_group'], backbone='yes' if r['backbone'] else 'no'))
    write_tsv(case / 'inputs/references.tsv', refrows, REF_FIELDS)
    write_tsv(case / 'inputs/reference_sequences.tsv', labelrows, LABEL_FIELDS)
    read_model = simulate_reads(targets, case / 'inputs/reads.fastq.gz', truth_dir,
                                random.Random(seed + 1000003), config.get('coverage', 12),
                                config.get('read_length', 10000), config.get('read_error', 0.001))
    manifest = dict(schema='chromosort-input-v1', case_id=cid,
                    assembly=dict(path='inputs/assembly.fa', sha256=digest(case / 'inputs/assembly.fa'),
                                  sample_id=cid, stage='input'),
                    references='inputs/references.tsv', reference_sequences='inputs/reference_sequences.tsv',
                    evidence=[], readsets=[dict(readset_id='sim01', path='inputs/reads.fastq.gz',
                                               sha256=digest(case / 'inputs/reads.fastq.gz'),
                                               evidence_class='mapped_simulated_reads', dependency_group='sim01')],
                    coordinate_system='0-based-half-open', split=split)
    save_json(case / 'manifest.json', manifest)
    save_json(truth_dir / 'truth.json', dict(schema='benchmark-truth-v1', case_id=cid, config=config,
              assembly_sha256=digest(case / 'inputs/assembly.fa'), target_sha256=digest(truth_dir / 'target.fa'),
              mapping='input_to_target.tsv', events=events, preservation_controls=variants,
              target_labels=labels, source_constructions=source_maps, read_model=read_model,
              baseline_mask='All non-ACGT target positions; reported separately from known sequence',
              breakpoint_tolerance_bp=500, correct_endpoint_anchor_bp=100,
              limitations=['No platform-calibrated read model', 'No de novo phasing or inferred subgenome labels',
                           'Locus cases are not whole-genome accuracy estimates']))
    save_json(truth_dir / 'config.json', config)
    print(cid, sum(map(len, assembly.values())), 'input bp;', read_model['read_count'], 'reads')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--config', type=Path, default=ROOT / 'configs/development.json')
    p.add_argument('--output-root', type=Path, default=ROOT / 'work')
    p.add_argument('--source-dir', type=Path)
    p.add_argument('--cases', nargs='+', help='Subset of configured development IDs; default all 14')
    a = p.parse_args()
    configs = json.loads(a.config.read_text())
    selected = set(a.cases or [c['case_id'] for c in configs['cases']])
    if selected - {c['case_id'] for c in configs['cases']}:
        p.error('Unknown development case ID')
    for config in configs['cases']:
        if config['case_id'] in selected:
            generate(config, a.output_root, a.source_dir)
