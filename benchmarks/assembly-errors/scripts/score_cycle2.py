#!/usr/bin/env python3
"""Commit both completed software output trees, then score the fixed cycle2 pass."""
import argparse
import collections
import json
import re
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'src'), str(ROOT / 'scripts')]
from benchio import digest, read_fasta, save_json, write_tsv


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--baseline-root', type=Path, required=True)
    p.add_argument('--candidate-root', type=Path, required=True)
    p.add_argument('--panel-root', type=Path, required=True)
    p.add_argument('--prior-root', type=Path, required=True)
    p.add_argument('--output-root', type=Path, required=True)
    args = p.parse_args()
    out = args.output_root.resolve()
    if out.exists():
        p.error('Choose a fresh score directory')
    roots = {'baseline': args.baseline_root.resolve(), 'candidate': args.candidate_root.resolve()}
    freeze = json.loads((ROOT / 'cycle2/execution-freeze.json').read_text())
    for rel, sha in freeze['files'].items():
        assert digest(ROOT / rel) == sha, 'Execution-frozen file changed: '+rel
    commitments = {}
    for role, root in roots.items():
        commitment = json.loads((root / 'outputs.commitment.json').read_text())
        for rel, sha in commitment.items():
            assert digest(root / rel) == sha, 'Changed committed output: '+role+'/'+rel
        commitments[role] = dict(commitment_sha256=digest(root / 'outputs.commitment.json'), files=commitment)
    out.mkdir(parents=True)
    # This public ledger is written before importing the evaluator or opening truth.
    save_json(out / 'paired-output-commitment.json', dict(schema='cycle2-paired-output-commitment-v1',
        phase='Both completed software trees verified and recorded before evaluator import',
        execution_freeze_sha256=digest(ROOT / 'cycle2/execution-freeze.json'), software=commitments))
    from score_structure import score_structure
    records, inspection, summary = [], [], []
    primary = json.loads((ROOT / 'cycle2/protocol.json').read_text())['case_ids']
    identities = {}
    for role, root in roots.items():
        identities[role] = json.loads((root / 'execution.json').read_text())['software']
        for run_path in sorted(root.rglob('run.json')):
            run = json.loads(run_path.read_text())
            cid = run['case_id']
            panel = args.panel_root.resolve() if cid in primary else args.prior_root.resolve()
            case = panel / 'development' / cid
            truth = panel / 'evaluator/development' / cid
            manifest = json.loads((case / 'manifest.json').read_text())
            directory = run_path.parent
            automatic = 'automatic' in run_path.relative_to(root).parts
            unedited = 'unedited' in run_path.relative_to(root).parts
            arm = run['arm'] if automatic or unedited else 'scripted_reviewed_repair'
            record = dict(role=role, case_id=cid, cohort='cycle2_primary' if cid in primary else 'prior_diagnostic',
                arm=arm, status=run['status'], workflow_class=run['workflow_class'],
                run_sha256=digest(run_path), declared_software_identity=identities[role]['declared_identity'])
            if run['status'] == 'complete':
                assert digest(directory / 'output.fa') == run['output_sha256']
                assert digest(directory / 'components.tsv') == run['ledger_sha256']
                value = score_structure(truth, case / manifest['assembly']['path'], directory / 'output.fa', directory / 'components.tsv')
                record['structure_score'] = value
                base = value['frozen_score']
                summary.append(dict(role=role, cohort=record['cohort'], case_id=cid, arm=arm, status='complete',
                    matched_cuts=base['matched_error_breakpoints'], expected_cuts=base['expected_error_breakpoints'],
                    false_cuts=base['false_breakpoints'], complete_structure=value['complete_structure'],
                    exact_once_input=base['retain_all_exact'], target_missing_bp=base['target_missing_callable_bp'],
                    target_extra_bp=base['target_extra_copy_callable_bp'], broken_controls=sum(not c['continuous'] for c in base['preservation_controls']),
                    remaining_adjacency_errors=base['target_relation_discrepancy_count']))
            else:
                record['error'] = run.get('error', run['status'])
            records.append(record)
            if not automatic and not unedited:
                proposals = run.get('inspection', {}).get('proposals', [])
                truth_value = json.loads((truth / 'truth.json').read_text())
                bounds = [e['breakpoints'] for e in truth_value['events'] if e['kind'] == 'orientation']
                intervals = []
                for row in proposals:
                    match = re.search(r'internal interval \[(\d+), (\d+)\)', row['reason'])
                    if row['action'] == 'inspect' and match:
                        intervals.append([int(match[1]), int(match[2])])
                panels = [{k: panel[k] for k in ['breakpoint', 'spanning_reads', 'split_reads', 'clipped_near_boundary',
                    'molecules_in_window', 'alignments_in_window', 'evidence_role', 'coverage_detail'] if k in panel}
                    for panel in run.get('supplied_boundary_panels', [])]
                inspection.append(dict(role=role, case_id=cid, cohort=record['cohort'], status=run['status'],
                    scan_state=run.get('inspection', {}).get('state'), proposals=proposals,
                    pending_proposals=sum(r['decision'] == 'pending' for r in proposals),
                    inspect_intervals=intervals, orientation_both_edges_within_500bp=[any(abs(a-x)<=500 and abs(b-y)<=500 for x,y in intervals) for a,b in bounds],
                    supplied_boundary_panels=panels, replay_identical=run.get('replay_identical'),
                    fresh_validation_state=run.get('fresh_stage_validation', {}).get('state'),
                    biological_validation=run.get('fresh_stage_validation', {}).get('biological_validation')))
    with (out / 'scores.jsonl').open('w') as handle:
        for record in records:
            handle.write(json.dumps(record, sort_keys=True)+'\n')
    write_tsv(out / 'summary.tsv', summary, list(summary[0]))
    save_json(out / 'inspection-and-replay.json', inspection)
    save_json(out / 'software-identities.json', identities)
    pair_index = {(r['role'], r['case_id'], r['arm']): r for r in records}
    paired = []
    for cid in primary:
        for policy in ['report', 'review']:
            arm = f'chromosort_{policy}_comprehensive'
            a, b = [pair_index[(role, cid, arm)] for role in ['baseline', 'candidate']]
            valid = a['status'] == b['status'] == 'complete'
            paired.append(dict(case_id=cid, policy=policy, both_complete=valid,
                exact_fasta_and_ledger_equal=valid and all(a['structure_score']['frozen_score'][k] == b['structure_score']['frozen_score'][k] for k in ['output_sha256', 'ledger_sha256']),
                structural_score_equal=valid and a['structure_score'] == b['structure_score']))
    totals = collections.Counter((r['role'], r['cohort'], r['arm'], r['status']) for r in records)
    save_json(out / 'evaluation-summary.json', dict(schema='cycle2-evaluation-summary-v1', scored_outputs=len(summary),
        recorded_runs=len(records), arms=[dict(role=k[0], cohort=k[1], arm=k[2], status=k[3], count=v) for k,v in totals.items()],
        automatic_pairs=paired, source_task_constraint='No candidate tuning after this independent pass',
        interpretation='Scripted repairs use supplied coordinates. Automatic correction, pending inspection and complete target structure are separate outcomes.'))
    print(json.dumps(dict(scored_outputs=len(summary), recorded_runs=len(records),
        automatic_pairs=len(paired), identical_automatic_pairs=sum(r['exact_fasta_and_ledger_equal'] for r in paired)), indent=2))


if __name__ == '__main__':
    main()
