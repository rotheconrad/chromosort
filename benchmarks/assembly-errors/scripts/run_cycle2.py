#!/usr/bin/env python3
"""Run one immutable cycle2 software arm without opening evaluator truth."""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'src'), str(ROOT / 'scripts')]
from benchio import digest, read_fasta, read_tsv, save_json, write_tsv
from adapt_agp import read_agp, FIELDS
from run_benchmark import chromosort_identity, input_hashes


def manual_recipe(assembly, edit):
    pieces = []
    for name, sequence in read_fasta(assembly).items():
        if name != edit['contig']:
            spans = [(0, len(sequence), '+', name)]
        elif edit['action'] == 'cut':
            cut = edit['position']
            spans = [(0, cut, '+', name+'_1'), (cut, len(sequence), '+', name+'_2')]
        else:
            a, b = edit['start'], edit['end']
            spans = [(0, a, '+', name), (a, b, '-', name), (b, len(sequence), '+', name)]
        for s, e, strand, group in spans:
            if s < e:
                pieces.append(dict(id=f'p{len(pieces)+1}', source=name, name=group, start=s+1,
                    end=e, strand=strand, scaffold=group, removed=False))
    return dict(schema='chromosort-manual-v1', sourceAssemblySha256=digest(assembly),
        requireCompleteCoverage=True, scaffoldEnabled=True, gapBp=0, pieces=pieces,
        benchmark_class='supplied-oracle-coordinate scripted development replay')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--python', type=Path, required=True)
    p.add_argument('--software-id', required=True)
    p.add_argument('--role', choices=['baseline', 'candidate'], required=True)
    p.add_argument('--case-root', type=Path, required=True)
    p.add_argument('--prior-case-root', type=Path, required=True)
    p.add_argument('--output-root', type=Path, required=True)
    args = p.parse_args()
    out, case_root, prior = args.output_root.resolve(), args.case_root.resolve(), args.prior_case_root.resolve()
    if out.exists():
        p.error('Choose a fresh output root')
    if out == ROOT or (ROOT in out.parents and out.relative_to(ROOT).parts[0] != 'work'):
        p.error('Use work/ or an external output root')
    out.mkdir(parents=True)
    env = dict(os.environ)
    env.pop('PYTHONPATH', None)
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    python = str(args.python.absolute())
    identity = chromosort_identity(Path(python), env)
    identity.update(declared_identity=args.software_id, role=args.role)
    protocol = json.loads((ROOT / 'cycle2/protocol.json').read_text())
    edits = {k: v[0] for k, v in protocol['scripted_review_arm']['edits'].items()}
    for cid, record in json.loads((ROOT / 'cycle2/prior-case-diagnostics.json').read_text())['cases'].items():
        edits[cid] = dict(action=record['action'], contig=record['input_id'], start=record['start'], end=record['end'])
    commands, runs = [], []
    def run(arguments, label):
        command = [python, '-m', 'chromosort.cli', *map(str, arguments)]
        start = time.monotonic()
        log = out / (label + '.log')
        log.parent.mkdir(parents=True, exist_ok=True)
        with log.open('w') as handle:
            result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, env=env)
        commands.append(dict(label=label, command=['selected-python', *[s.replace(str(out), '{output_root}').replace(str(case_root), '{case_root}').replace(str(prior), '{prior_case_root}') for s in command[1:]]],
            returncode=result.returncode, elapsed_seconds=round(time.monotonic()-start, 6)))
        if result.returncode:
            raise RuntimeError(f'{label} failed; see its log')
    automatic = [sys.executable, str(ROOT / 'scripts/run_benchmark.py'), '--tool', 'chromosort',
        '--python', python, '--software-id', args.software_id, '--case-root', str(case_root),
        '--output-root', str(out / 'automatic'), '--modes', 'comprehensive', '--policies', 'report', 'review', '--defer-scoring']
    auto = subprocess.run(automatic, env=env)
    if auto.returncode:
        print('Automatic failures retained in automatic run records', flush=True)
    if args.role == 'baseline':
        subprocess.run([sys.executable, str(ROOT / 'scripts/run_benchmark.py'), '--tool', 'unedited',
            '--case-root', str(case_root), '--output-root', str(out / 'unedited'), '--defer-scoring'], check=True, env=env)
    for cid in sorted(edits):
        case = (case_root if cid in protocol['case_ids'] else prior) / cid
        manifest = json.loads((case / 'manifest.json').read_text())
        hashes = input_hashes(case, manifest)
        assembly = case / manifest['assembly']['path']
        directory = out / cid
        directory.mkdir()
        scan, edit = directory / 'scan', edits[cid]
        record = dict(case_id=cid, software=identity, input_files_sha256=hashes,
            workflow_class='scripted_oracle_coordinate_development_repair', status='pending')
        try:
            run(['workflow', 'scan', '--manifest', case / 'manifest.json', '--mode', 'comprehensive',
                '--read-evidence-id', 'reads_sim01', '--inspect-min-gap-bp', '1000', '--output-dir', scan], cid+'/scan')
            scan_record = json.loads((scan / 'scan.json').read_text())
            record['inspection'] = dict(state=scan_record['state'], proposals=scan_record['proposals'], figures=scan_record['figures'])
            # Diagnostic panels use the frozen supplied positions, never inferred truth.
            points = [edit['position']] if edit['action'] == 'cut' else [edit['start'], edit['end']]
            record['supplied_boundary_panels'] = []
            for index, point in enumerate(points):
                prefix = directory / ('supplied_edge'+str(index))
                run(['reads', '--manifest', case / 'manifest.json', '--evidence-id', 'reads_sim01',
                    '--contig', edit['contig'], '--start', max(0, point-10000), '--end', min(len(read_fasta(assembly)[edit['contig']]), point+10000),
                    '--breakpoint', point, '--min-mapq', 20, '--min-anchor-bp', 1000, '--formats', 'svg',
                    '--action', 'supplied-coordinate development diagnostic', '-o', prefix], cid+'/panel'+str(index))
                data = json.loads(Path(str(prefix)+'.figure.json').read_text())
                record['supplied_boundary_panels'].append({k: v for k, v in data.items() if k not in ['outputs', 'evidence', 'depth_bins', 'input_files', 'parameters']})
            output, agp = directory / 'output.fa', directory / 'output.agp'
            if args.role == 'baseline':
                recipe = directory / 'recipe.json'
                save_json(recipe, manual_recipe(assembly, edit))
                run(['manual', 'apply', '--assembly-fasta', assembly, '--recipe', recipe,
                    '--output-fasta', output, '--agp', agp], cid+'/manual')
            else:
                decisions = read_tsv(scan / 'decisions.tsv')
                links = [r for r in decisions if r['target'] == edit['contig']]
                if not links:
                    record['status'] = 'unavailable_no_scan_event'
                    save_json(directory / 'run.json', record)
                    runs.append(record)
                    continue
                # Resolve an exact proposed cut in its original decision row. Other
                # requested cuts/reversals get an explicit linked edit record.
                existing = next((r for r in links if edit['action'] == 'cut' and r['action'] == 'cut'
                    and int(r['position']) == edit['position']), None)
                for row in decisions:
                    row.update(decision='accept' if row is existing else 'reject', reviewer='scripted-development-oracle',
                        notes='Supplied prespecified edit only; preserve other native structure. Not a blinded review.')
                decision_path = directory / 'scripted-decisions.tsv'
                write_tsv(decision_path, decisions, list(decisions[0]))
                extra = []
                if existing is None:
                    point = edit.get('position')
                    run(['workflow', 'edit', '--scan-dir', scan, '--event-id', links[0]['event_id'],
                        '--action', edit['action'], '--start', point if point is not None else edit['start'],
                        '--end', point if point is not None else edit['end'], '--notes',
                        'Supplied prespecified oracle-coordinate development edit; not independently inferred.',
                        '--output-dir', directory / 'explicit'], cid+'/edit')
                    edit_path = directory / 'explicit/reviewed_edits.tsv'
                    rows = read_tsv(edit_path)
                    for row in rows:
                        row.update(decision='accept', reviewer='scripted-development-oracle')
                    reviewed = directory / 'scripted-edits.tsv'
                    write_tsv(reviewed, rows, list(rows[0]))
                    extra = ['--reviewed-edits', reviewed]
                common = ['workflow', 'apply', '--scan-dir', scan, '--decisions', decision_path, *extra]
                run([*common, '--plan-only', '--output-dir', directory / 'preview'], cid+'/preview')
                assert not (directory / 'preview/reviewed.fa').exists()
                run([*common, '--output-dir', directory / 'applied'], cid+'/apply')
                shutil.copyfile(directory / 'applied/reviewed.fa', output)
                shutil.copyfile(directory / 'applied/reviewed.fa.agp', agp)
                recipe = directory / 'applied/recipe.json'
                run(['workflow', 'align', '--manifest', case / 'manifest.json', '--assembly-fasta', directory / 'applied/reviewed.fa',
                    '--output-dir', directory / 'realigned', '--stage', 'reviewed', '--threads', 2], cid+'/realign')
                run(['workflow', 'validate', '--manifest', directory / 'realigned/manifest.json',
                    '--apply-audit', directory / 'applied/apply.json', '--output', directory / 'validation.json'], cid+'/validate')
                record['fresh_stage_validation'] = json.loads((directory / 'validation.json').read_text())
            replay = directory / 'replay.fa'
            run(['manual', 'apply', '--assembly-fasta', assembly, '--recipe', recipe, '--output-fasta', replay], cid+'/replay')
            record['replay_identical'] = output.read_bytes() == replay.read_bytes()
            if not record['replay_identical']:
                raise ValueError('Manual recipe replay bytes differ')
            write_tsv(directory / 'components.tsv', read_agp(agp), FIELDS)
            record.update(status='complete', output_sha256=digest(output), ledger_sha256=digest(directory / 'components.tsv'))
        except Exception as error:
            record.update(status='failed', error=str(error))
        save_json(directory / 'run.json', record)
        runs.append(record)
        print(cid, args.role, record['status'], flush=True)
    save_json(out / 'execution.json', dict(software=identity, commands=commands, cases=[dict(case_id=r['case_id'], status=r['status']) for r in runs],
        scoring='No evaluator truth opened; await both software output commitments'))
    files = {str(p.relative_to(out)): digest(p) for p in sorted(out.rglob('*')) if p.is_file()}
    save_json(out / 'outputs.commitment.json', files)


if __name__ == '__main__':
    main()
