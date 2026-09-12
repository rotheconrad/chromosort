#!/usr/bin/env python3
"""Independent interval/provenance scorer, with exact FASTA validation.

The scorer imports only generic file helpers, never generator or curation code.
Truth identifies *input* bases, so repeated copies are not reassigned by alignment.
"""
import argparse
import collections
import json
import re
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from benchio import digest, read_fasta, read_tsv, rc, save_json

FIELDS = ['output_id', 'output_start', 'output_end', 'input_id', 'input_start', 'input_end', 'strand']


def union_intervals(ranges):
    out = []
    for a, b in sorted(ranges):
        if a == b:
            continue
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def length(ranges):
    return sum(b - a for a, b in ranges)


def intersect_length(a, b):
    a, b = union_intervals(a), union_intervals(b)
    i = j = total = 0
    while i < len(a) and j < len(b):
        total += max(0, min(a[i][1], b[j][1]) - max(a[i][0], b[j][0]))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return total


def numeric(rows, names):
    for r in rows:
        for k in names:
            r[k] = int(r[k])
    return rows


def score(truth_dir, assembly_path, output_path, ledger_path, assignments_path=None):
    truth_dir = Path(truth_dir)
    truth = json.loads((truth_dir / 'truth.json').read_text())
    if digest(assembly_path) != truth['assembly_sha256']:
        raise ValueError('Input FASTA does not match frozen truth')
    if digest(truth_dir / 'target.fa') != truth['target_sha256']:
        raise ValueError('Target FASTA does not match frozen truth')
    inputs, outputs, targets = read_fasta(assembly_path), read_fasta(output_path), read_fasta(truth_dir / 'target.fa')
    maps = numeric(read_tsv(truth_dir / truth['mapping']), ['input_start', 'input_end', 'target_start', 'target_end'])
    by_input = collections.defaultdict(list)
    for r in maps:
        by_input[r['input_id']].append(r)
    for iid, seq in inputs.items():
        rows = sorted(by_input[iid], key=lambda r:r['input_start'])
        pos = 0
        for r in rows:
            if r['input_start'] != pos or r['input_end'] - pos != r['target_end'] - r['target_start']:
                raise ValueError('Invalid frozen input map')
            t = targets[r['target_id']][r['target_start']:r['target_end']]
            if r['strand'] == '-':
                t = rc(t)
            if seq[r['input_start']:r['input_end']] != t:
                raise ValueError('Frozen input map is not sequence-exact')
            pos = r['input_end']
        if pos != len(seq):
            raise ValueError('Frozen map misses input bases')
    ledger = numeric(read_tsv(ledger_path), ['output_start','output_end','input_start','input_end'])
    by_out = collections.defaultdict(list)
    for r in ledger:
        by_out[r['output_id']].append(r)
    if set(by_out) != set(outputs):
        raise ValueError('Ledger/output FASTA record sets differ')
    projected = collections.defaultdict(list)
    retained_inputs = collections.defaultdict(list)
    retained_targets = collections.defaultdict(list)
    output_gaps = 0
    for oid, seq in outputs.items():
        rows = sorted(by_out[oid], key=lambda r:r['output_start'])
        pos = 0
        for r in rows:
            if r['output_start'] != pos or r['output_end'] <= pos or r['output_end'] > len(seq):
                raise ValueError('Ledger has overlapping, missing, or out-of-bounds output bases')
            n = r['output_end'] - pos
            if r['input_id'] == '.':
                if seq[pos:pos + n] != 'N' * n:
                    raise ValueError('Gap ledger includes sequence other than N')
                output_gaps += n
                pos += n
                continue
            iid, s, e, strand = r['input_id'], r['input_start'], r['input_end'], r['strand']
            if iid not in inputs or not 0 <= s < e <= len(inputs[iid]) or e - s != n or strand not in ['+', '-']:
                raise ValueError('Invalid input ledger interval')
            expected = inputs[iid][s:e]
            if strand == '-':
                expected = rc(expected)
            if seq[pos:pos + n] != expected:
                raise ValueError('Output FASTA disagrees with the supplied provenance ledger')
            retained_inputs[iid].append((s, e))
            pieces = []
            for m in by_input[iid]:
                a, b = max(s, m['input_start']), min(e, m['input_end'])
                if a >= b:
                    continue
                if m['strand'] == '+':
                    ts, te = m['target_start'] + a - m['input_start'], m['target_start'] + b - m['input_start']
                else:
                    ts, te = m['target_end'] - b + m['input_start'], m['target_end'] - a + m['input_start']
                outstart = pos + (a - s if strand == '+' else e - b)
                p = dict(target_id=m['target_id'], start=ts, end=te,
                         strand='+' if strand == m['strand'] else '-',
                         output_start=outstart, output_end=outstart + b - a)
                pieces.append(p)
                retained_targets[m['target_id']].append((ts, te))
            projected[oid].extend(sorted(pieces, key=lambda p:p['output_start']))
            pos += n
        if pos != len(seq):
            raise ValueError('Ledger misses trailing output sequence')
    # New cuts: interior input boundaries retaining sequence on both sides.
    # Terminal trimming is recorded as loss, not a correction breakpoint.
    cuts = []
    for iid, ranges in retained_inputs.items():
        points = sorted({x for r in ranges for x in r if 0 < x < len(inputs[iid])})
        for point in points:
            if any(a < point <= b for a,b in ranges) and any(a <= point < b for a,b in ranges):
                cuts.append((iid, point))
    expected_cuts = [(e['input_id'], p) for e in truth['events'] if e['kind'] != 'collapse_missing_sequence' for p in e['breakpoints']]
    # One-to-one matching prevents many nearby cuts from receiving one credit.
    candidates = sorted((abs(p-q), i, j) for i,(iid,p) in enumerate(cuts)
                        for j,(eid,q) in enumerate(expected_cuts)
                        if iid == eid and abs(p-q) <= truth['breakpoint_tolerance_bp'])
    matched_i, matched_j, localization = set(), set(), []
    for distance, i, j in candidates:
        if i not in matched_i and j not in matched_j:
            matched_i.add(i); matched_j.add(j); localization.append(distance)
    relation_errors, paths = [], []
    for oid, pieces in projected.items():
        current = None
        for p in pieces:
            if current is None:
                current = dict(p)
                continue
            left = current['end'] if current['strand'] == '+' else current['start']
            right = p['start'] if p['strand'] == '+' else p['end']
            outgap = p['output_start'] - current['output_end']
            same = current['target_id'] == p['target_id'] and current['strand'] == p['strand']
            tg = (right - left) * (1 if p['strand'] == '+' else -1) if same else None
            if same and tg == 0 and outgap == 0:
                current['start'], current['end'] = min(current['start'],p['start']), max(current['end'],p['end'])
                current['output_end'] = p['output_end']
            else:
                if not same or tg != outgap:
                    relation_errors.append(dict(output_id=oid, output_position=current['output_end'],
                                                kind='orientation_or_interchromosomal' if not same else ('overlap' if tg < 0 else 'gap_length'),
                                                target_gap=tg, output_gap=outgap,
                                                gap_error_bp=None if tg is None else outgap-tg))
                paths.append(dict(current, output_id=oid))
                current = dict(p)
        if current:
            paths.append(dict(current, output_id=oid))
    input_total = sum(map(len,inputs.values()))
    input_union = sum(length(union_intervals(r)) for r in retained_inputs.values())
    input_mult = sum(length(r) for r in retained_inputs.values())
    callable_runs = {tid:[(m.start(),m.end()) for m in re.finditer('[ACGT]+',seq)] for tid,seq in targets.items()}
    callable_total = sum(length(r) for r in callable_runs.values())
    callable_retained = sum(intersect_length(retained_targets[tid],runs) for tid,runs in callable_runs.items())
    callable_mult = sum(sum(intersect_length([r],callable_runs[tid]) for r in rs) for tid,rs in retained_targets.items())
    controls = []
    for c in truth['preservation_controls']:
        interval = [(c['start'], c['end'])]
        denom = intersect_length(interval, callable_runs[c['target_id']])
        keep = intersect_length(union_intervals(retained_targets[c['target_id']]), interval)
        continuous = any(p['target_id'] == c['target_id'] and p['start'] <= c['start'] and p['end'] >= c['end'] for p in paths)
        controls.append(dict(variant_id=c['variant_id'],kind=c['kind'],retained_bp=keep,
                             interval_bp=c['end']-c['start'],callable_bp=denom,continuous=continuous))
    event_scores = []
    for e in truth['events']:
        expected = [(e['input_id'], p) for p in e['breakpoints']]
        nmatch = sum(1 for j in matched_j if expected_cuts[j] in expected)
        if e['kind'] == 'collapse_missing_sequence':
            missing_recovered = sum(intersect_length(retained_targets[tid],[(s,t)]) for tid,s,t in e['target_intervals'])
            event_scores.append(dict(event_id=e['event_id'],kind=e['kind'],missing_bp=e['missing_bp'],
                                     missing_bp_recovered=missing_recovered,break_recall=None,
                                     interpretation='Placement cannot recover unavailable sequence; report unresolved'))
        else:
            event_scores.append(dict(event_id=e['event_id'],kind=e['kind'],expected_breakpoints=len(expected),
                                     matched_breakpoints=nmatch,break_recall=nmatch/len(expected),
                                     all_expected_breakpoints_localized=nmatch==len(expected)))
    assignment = dict(status='not_supplied', assigned_bp=None, correct_assigned_bp=None, ambiguous_required_bp=None)
    if assignments_path:
        assignments = {r['output_id']:r for r in read_tsv(assignments_path)}
        if set(assignments) != set(outputs):
            raise ValueError('Assignments must explicitly account for every output record')
        assigned = correct = required = wrongly_forced = 0
        for oid, pieces in projected.items():
            ar = assignments[oid]
            for piece in pieces:
                n = piece['end']-piece['start']; label = truth['target_labels'][piece['target_id']]
                ambiguous = ar['status'] in ['ambiguous','unplaced']
                if label.get('require_ambiguous'):
                    required += n
                    if not ambiguous: wrongly_forced += n
                if not ambiguous:
                    assigned += n
                    if ar['output_group'] == label['output_group']: correct += n
        assignment = dict(status='scored', assigned_bp=assigned, correct_assigned_bp=correct,
                          accuracy_among_assigned=correct/assigned if assigned else None,
                          ambiguous_required_bp=required, inappropriate_forced_assignment_bp=wrongly_forced)
    return dict(schema='benchmark-score-v1', case_id=truth['case_id'],
                output_sha256=digest(output_path), ledger_sha256=digest(ledger_path),
                truth_sha256=digest(truth_dir / 'truth.json'), sequence_provenance_valid=True,
                input_bp=input_total, input_retained_bp=input_union, input_unaccounted_bp=input_total-input_union,
                input_duplicated_bp=input_mult-input_union, retain_all_exact=input_union==input_total and input_mult==input_total,
                target_callable_bp=callable_total, target_retained_callable_bp=callable_retained,
                target_retention=callable_retained/callable_total if callable_total else None,
                target_missing_callable_bp=callable_total-callable_retained,
                target_extra_copy_callable_bp=callable_mult-callable_retained,
                target_masked_bp=sum(map(len,targets.values()))-callable_total,
                new_breakpoints=len(cuts), expected_error_breakpoints=len(expected_cuts),
                matched_error_breakpoints=len(matched_j), false_breakpoints=len(cuts)-len(matched_i),
                false_breakpoints_per_callable_mb=(len(cuts)-len(matched_i))*1e6/callable_total,
                breakpoint_localization_error_bp=localization, cut_precision=len(matched_i)/len(cuts) if cuts else None,
                cut_recall=len(matched_j)/len(expected_cuts) if expected_cuts else None,
                output_gap_bp=output_gaps, target_relation_discrepancies=relation_errors,
                target_relation_discrepancy_count=len(relation_errors), events=event_scores,
                preservation_controls=controls, assignments=assignment,
                interpretation='Cut localization is not full repair. Read retention, copy counts, structural discrepancies and control continuity jointly. Missing/discarded input needs an external policy ledger before claiming accounted retain-all.')


if __name__ == '__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--truth',type=Path,required=True)
    p.add_argument('--assembly',type=Path,required=True)
    p.add_argument('--output-fasta',type=Path,required=True)
    p.add_argument('--ledger',type=Path,required=True)
    p.add_argument('--assignments',type=Path)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    save_json(a.out,score(a.truth,a.assembly,a.output_fasta,a.ledger,a.assignments))
