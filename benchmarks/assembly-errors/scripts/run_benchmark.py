#!/usr/bin/env python3
"""Run explicit development arms, commit output hashes, then independently score.

Install the chosen software separately and supply its commit/archive identity.
The benchmark does not download, patch, tune, or publish third-party software.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from urllib.parse import quote
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'src'), str(ROOT / 'scripts')]
from benchio import digest, read_fasta, read_tsv, save_json, write_fasta, write_tsv
from adapt_agp import read_agp, FIELDS


def input_hashes(case, manifest):
    paths = ['manifest.json', manifest['assembly']['path'], manifest['references'], manifest['reference_sequences']]
    assert digest(case / manifest['assembly']['path']) == manifest['assembly']['sha256']
    for ref in read_tsv(case / manifest['references']):
        assert digest(case / ref['fasta']) == ref['sha256']
        paths.append(ref['fasta'])
    for evidence in manifest['evidence']:
        assert evidence['assembly_sha256'] == manifest['assembly']['sha256']
        assert digest(case / evidence['path']) == evidence['sha256']
        paths.append(evidence['path'])
        if evidence.get('index_path'):
            assert digest(case / evidence['index_path']) == evidence['index_sha256']
            paths.append(evidence['index_path'])
    for reads in manifest.get('readsets', []):
        if reads.get('path'):
            assert digest(case / reads['path']) == reads['sha256']
            paths.append(reads['path'])
    return {p: digest(case / p) for p in sorted(set(paths))}


def chromosort_identity(python, env):
    code = '''import hashlib,json,sys,chromosort
from pathlib import Path
root=Path(chromosort.__file__).parent
files={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(root.rglob('*.py'))}
print(json.dumps({'version':chromosort.__version__,'python':sys.version.split()[0],
 'source_tree_sha256':hashlib.sha256(json.dumps(files,sort_keys=True).encode()).hexdigest(),'source_files':files}))
'''
    return json.loads(subprocess.check_output([str(python), '-c', code], env=env, text=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tool', choices=['unedited', 'chromosort', 'ragtag'], required=True)
    parser.add_argument('--case-root', type=Path, default=ROOT / 'fixtures/development')
    parser.add_argument('--truth-root', type=Path, default=ROOT / 'fixtures/truth')
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--cases', nargs='+')
    parser.add_argument('--software-id', help='Exact software git commit or source archive SHA-256')
    parser.add_argument('--python', type=Path, default=Path(sys.executable), help='Python with the selected ChromoSort installed')
    parser.add_argument('--modes', nargs='+', choices=['conservative', 'comprehensive', 'sensitive'],
                        default=['conservative', 'comprehensive', 'sensitive'])
    parser.add_argument('--policies', nargs='+', choices=['review', 'report'], default=['review', 'report'])
    parser.add_argument('--ragtag-bin', type=Path, help='Directory containing the installed RagTag subcommand scripts')
    parser.add_argument('--arms', nargs='+', help='RagTag protocol arms; default all seven')
    parser.add_argument('--defer-scoring', action='store_true',
                        help='Commit outputs without opening truth, for paired multi-version evaluation')
    args = parser.parse_args()
    case_root, truth_root, output_root = [p.resolve() for p in (args.case_root, args.truth_root, args.output_root)]
    if output_root == ROOT or (ROOT in output_root.parents and output_root.relative_to(ROOT).parts[0] != 'work'):
        parser.error('Write generated outputs to work/ or outside this benchmark directory')
    if output_root.exists():
        parser.error('Output root already exists; use a fresh directory')
    if args.tool != 'unedited' and not args.software_id:
        parser.error('--software-id must identify the selected immutable software')
    selected = args.cases or [p.parent.name for p in sorted(case_root.glob('D*/manifest.json'))]
    if not selected or len(selected) != len(set(selected)):
        parser.error('Select one or more unique development cases')
    env = dict(os.environ)
    env.pop('PYTHONPATH', None)
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    software = {'declared_identity': args.software_id or 'identity-baseline', 'tool': args.tool}
    protocol = json.loads((ROOT / 'configs/candidate-0.4.0rc2.json').read_text())
    if args.tool == 'chromosort':
        software.update(chromosort_identity(args.python.absolute(), env))
        arms = [(f'chromosort_{policy}_{mode}', 'correct', mode, policy)
                for policy in args.policies for mode in args.modes]
    elif args.tool == 'ragtag':
        if args.ragtag_bin:
            tool_dir = args.ragtag_bin.resolve()
        else:
            located = shutil.which('ragtag.py')
            if not located:
                parser.error('Install pinned RagTag or provide --ragtag-bin')
            tool_dir = Path(located).resolve().parent
        env['PATH'] = str(tool_dir) + os.pathsep + env.get('PATH', '')
        version = subprocess.check_output([str(tool_dir / 'ragtag.py'), '--version'], env=env, text=True).strip()
        if version != 'v2.1.0':
            parser.error('The saved comparator protocol requires RagTag v2.1.0')
        software.update(version=version, entrypoint_sha256={p.name: digest(p) for p in
                        [tool_dir / 'ragtag_correct.py', tool_dir / 'ragtag_scaffold.py']})
        v1 = json.loads((ROOT / 'configs/comparators.json').read_text())
        v2 = json.loads((ROOT / 'configs/comparators.v2.json').read_text())
        corrections = dict(v1['correction_arms'], **v2['correction_arms'])
        options = {k: ('correct', v) for k, v in corrections.items()}
        options.update({k: ('scaffold', v) for k, v in v1['scaffold_arms'].items()})
        if set(args.arms or options) - set(options):
            parser.error('Unknown RagTag arm; see configs/comparators*.json')
        arms = [(name, *options[name], None) for name in (args.arms or options)]
    else:
        arms = [('unedited', 'identity', None, None)]
    inputs = {}
    for cid in selected:
        case = case_root / cid
        manifest = json.loads((case / 'manifest.json').read_text())
        if manifest['split'] != 'development' or manifest['case_id'] != cid or not cid.startswith('D'):
            parser.error('Only development manifests are accepted')
        inputs[cid] = (manifest, input_hashes(case, manifest))
        if args.tool == 'ragtag' and any('{reads}' in arm[2] for arm in arms):
            if not any(r.get('path') for r in manifest.get('readsets', [])):
                parser.error('Read-validation arms require regenerated FASTQ; tracked fixtures contain PAF only')
    output_root.mkdir(parents=True)
    save_json(output_root / 'protocol.json', dict(tool=software, arms=arms, input_files={k: v[1] for k, v in inputs.items()},
        evidence='Exactly reads_sim01 for ChromoSort; RagTag uses its separate declared protocol.',
        saved_protocol_sha256={p.name: digest(p) for p in (ROOT / 'configs').glob('*.json')},
        scorer_sha256=digest(ROOT / 'scripts/score.py')))
    pending, committed, failures = [], {}, []
    for cid in selected:
        case = case_root / cid
        manifest, hashes = inputs[cid]
        assembly = case / manifest['assembly']['path']
        for name, operation, mode_or_params, policy in arms:
            out = output_root / cid / name
            out.mkdir(parents=True)
            output, ledger = out / 'output.fa', out / 'components.tsv'
            command = []
            start = time.monotonic()
            status = 'complete'
            if args.tool == 'unedited':
                shutil.copyfile(assembly, output)
                write_tsv(ledger, [dict(output_id=k, output_start=0, output_end=len(v), input_id=k,
                    input_start=0, input_end=len(v), strand='+') for k, v in read_fasta(assembly).items()], FIELDS)
            elif args.tool == 'chromosort':
                replacements = {'CASE/manifest.json': case / 'manifest.json', 'MODE': mode_or_params,
                    'POLICY': policy, 'OUTPUT/output.fa': output, 'OUTPUT/fix.tsv': out / 'fix.tsv',
                    'OUTPUT/output.agp': out / 'output.agp'}
                command = [str(args.python.absolute()), '-m', 'chromosort.cli'] + [str(replacements.get(x, x))
                           for x in protocol['fix_arguments']]
            else:
                combined = {}
                for ref in read_tsv(case / manifest['references']):
                    for seqid, seq in read_fasta(case / ref['fasta']).items():
                        combined[quote(ref['reference_id'], safe='') + '::' + quote(seqid, safe='')] = seq
                reference = out / 'references.combined.fa'
                write_fasta(reference, combined)
                reads = next((case / r['path'] for r in manifest.get('readsets', []) if r.get('path')), None)
                command = [str(tool_dir / f'ragtag_{operation}.py'), str(reference), str(assembly),
                           '-t', '2', '-o', str(out / 'tool')]
                command += [str(reads) if x == '{reads}' else x for x in mode_or_params]
            if command:
                with (out / 'stdout.log').open('w') as stdout, (out / 'stderr.log').open('w') as stderr:
                    process = subprocess.run(command, env=env, stdout=stdout, stderr=stderr)
                if process.returncode:
                    status = 'failed'
                else:
                    agp = out / 'output.agp'
                    if args.tool == 'ragtag':
                        shutil.copyfile(out / 'tool' / f'ragtag.{operation}.fasta', output)
                        shutil.copyfile(out / 'tool' / f'ragtag.{operation}.agp', agp)
                    write_tsv(ledger, read_agp(agp, reverse=args.tool == 'ragtag' and operation == 'correct'), FIELDS)
            record = dict(schema='public-benchmark-run-v1', case_id=cid, arm=name, operation=operation,
                workflow_class='unedited_baseline' if args.tool == 'unedited' else 'automatic', status=status,
                software=software, input_files_sha256=hashes, elapsed_seconds=time.monotonic() - start,
                command=[x.replace(str(case_root), '{case_root}').replace(str(output_root), '{output_root}') for x in command],
                output_sha256=digest(output) if output.exists() else None,
                ledger_sha256=digest(ledger) if ledger.exists() else None)
            if command:
                record['command'][0] = 'selected-python' if args.tool == 'chromosort' else Path(command[0]).name
            save_json(out / 'run.json', record)
            committed[str((out / 'run.json').relative_to(output_root))] = digest(out / 'run.json')
            if status == 'complete':
                pending.append((cid, assembly, out))
            else:
                failures.append(cid + '/' + name)
            print(cid, name, status, flush=True)
    save_json(output_root / 'outputs.commitment.json', committed)
    if args.defer_scoring:
        print('Committed outputs without opening evaluator truth.', flush=True)
        if failures:
            raise RuntimeError('Failed arms: ' + ', '.join(failures))
        return
    # Outputs are immutable before evaluator truth is opened.
    from score import score
    for cid, assembly, out in pending:
        save_json(out / 'score.json', score(truth_root / cid, assembly, out / 'output.fa', out / 'components.tsv'))
    if failures:
        raise RuntimeError('Failed arms: ' + ', '.join(failures))
    print('Independently scored', len(pending), 'development outputs.', flush=True)


if __name__ == '__main__':
    main()
