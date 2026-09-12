#!/usr/bin/env python3
"""Export small exact ledgers and portable recipes from committed cycle2 outputs."""
import argparse
import json
import shutil
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from benchio import digest, read_tsv, save_json, write_tsv


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--baseline-root', type=Path, required=True)
    p.add_argument('--candidate-root', type=Path, required=True)
    p.add_argument('--results-root', type=Path, required=True)
    args = p.parse_args()
    results = args.results_root.resolve()
    destination = results / 'artifacts'
    if destination.exists():
        p.error('Artifacts already exist')
    commitment = json.loads((results / 'paired-output-commitment.json').read_text())
    artifacts, continuity, panels, resource = {}, [], [], []
    for role, root in [('baseline', args.baseline_root.resolve()), ('candidate', args.candidate_root.resolve())]:
        for rel, sha in commitment['software'][role]['files'].items():
            assert digest(root / rel) == sha
        execution = json.loads((root / 'execution.json').read_text())
        resource.append(dict(role=role, workflow_commands=len(execution['commands']),
            workflow_elapsed_seconds=sum(r['elapsed_seconds'] for r in execution['commands']),
            interpretation='Local descriptive subprocess elapsed sum; candidate includes more workflow stages, no speed benchmark'))
        for path in sorted(root.rglob('run.json')):
            record = json.loads(path.read_text())
            cid, directory = record['case_id'], path.parent
            arm = record.get('arm', 'scripted_reviewed_repair')
            output = destination / role / cid / arm
            output.mkdir(parents=True)
            for filename in ['components.tsv', 'output.agp', 'fix.tsv', 'fix.tsv.breakpoint_evidence.tsv']:
                if (directory / filename).exists():
                    shutil.copyfile(directory / filename, output / filename)
            recipe = directory / ('recipe.json' if role == 'baseline' else 'applied/recipe.json')
            if arm == 'scripted_reviewed_repair' and recipe.exists():
                data = json.loads(recipe.read_text())
                data.pop('sourceAssembly', None)
                save_json(output / 'recipe.json', data)
                artifacts[str((output / 'recipe.json').relative_to(results))] = dict(sha256=digest(output / 'recipe.json'),
                    original_sha256=digest(recipe), transformation='Remove advisory absolute sourceAssembly path if present; assembly SHA and all edit geometry preserved')
                if role == 'candidate':
                    for filename in ['decisions.tsv', 'decision_outcomes.tsv', 'reviewed_edits.tsv']:
                        source = directory / 'applied' / filename
                        if source.exists():
                            shutil.copyfile(source, output / filename)
                for panel in record.get('supplied_boundary_panels', []):
                    panels.append(dict(role=role, case_id=cid, supplied_position=panel['settings']['breakpoint'],
                        spanning_read_names=panel['spanning_reads'], split_read_names=panel['split_reads'],
                        molecules_in_window=panel['molecules_in_window'], alignments_in_window=panel['alignments_in_window'],
                        cigar_unavailable_rows=panel['filters']['cigar_unavailable'],
                        coverage_detail=panel['coverage_detail'], inference='Supplied-coordinate diagnostic; common molecule origin is not inferred from read IDs'))
            if arm.startswith('chromosort_'):
                for row in read_tsv(directory / 'fix.tsv.breakpoint_evidence.tsv'):
                    continuity.append(dict(role=role, case_id=cid, arm=arm, **row))
        auto_times = [json.loads(p.read_text())['elapsed_seconds'] for p in (root/'automatic').rglob('run.json')]
        resource[-1].update(automatic_commands=len(auto_times), automatic_elapsed_seconds=sum(auto_times))
    for path in sorted(destination.rglob('*')):
        if path.is_file():
            artifacts.setdefault(str(path.relative_to(results)), dict(sha256=digest(path), transformation='exact copy of committed tool output'))
    save_json(results / 'artifact-index.json', dict(schema='cycle2-small-artifacts-v1', files=artifacts,
        exclusions='Full FASTA, reads, figures, software archives and runtime logs remain regenerable local outputs'))
    write_tsv(results / 'continuity.tsv', continuity, list(continuity[0]))
    write_tsv(results / 'supplied-boundary-support.tsv', panels, list(panels[0]))
    save_json(results / 'resource-use.json', resource)
    print('Exported', len(artifacts), 'small artifacts')


if __name__ == '__main__':
    main()
