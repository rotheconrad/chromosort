#!/usr/bin/env python3
"""Map references and biological-target reads to each exact input assembly."""
import argparse
import json
import shlex
import subprocess
import sys
import time
from pathlib import Path
from urllib.parse import quote
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from benchio import digest, read_tsv, read_fasta, write_fasta, save_json

ROOT = Path(__file__).resolve().parents[1]


def run(cmd, cwd, output, log):
    start = time.monotonic()
    with open(output, 'wb') as out, open(log, 'wb') as err:
        p = subprocess.run(cmd, cwd=cwd, stdout=out, stderr=err)
    if p.returncode:
        raise RuntimeError(f'Command failed ({p.returncode}): {shlex.join(cmd)}; {log}')
    return dict(command=cmd, elapsed_seconds=round(time.monotonic() - start, 6),
                log=str(Path(log).relative_to(cwd)), output_sha256=digest(output),
                cost_note='minimap2 stderr includes process CPU and peak RSS; no GitHub Actions used')


def evidence(case, threads=2, bam=False):
    case = case.resolve()
    if ROOT / 'fixtures' in [case, *case.parents]:
        raise ValueError('Tracked fixtures are immutable; regenerate into a work directory')
    if threads < 1:
        raise ValueError('threads must be positive')
    mpath = case / 'manifest.json'
    manifest = json.loads(mpath.read_text())
    if manifest['evidence']:
        raise ValueError('Evidence already exists; refuse stale overwrite: ' + str(case))
    asm = case / manifest['assembly']['path']
    if digest(asm) != manifest['assembly']['sha256']:
        raise ValueError('Assembly checksum mismatch')
    (case / 'evidence').mkdir(exist_ok=True)
    (case / 'logs').mkdir(exist_ok=True)
    refs = read_tsv(case / manifest['references'])
    version = subprocess.check_output(['minimap2', '--version'], text=True).strip()
    history, entries, combined_seqs, combined_paf = [], [], {}, []
    for row in refs:
        rid = row['reference_id']
        ref = case / row['fasta']
        if digest(ref) != row['sha256']:
            raise ValueError('Reference checksum mismatch: ' + rid)
        paf = case / f'evidence/{rid}.assembly.paf'
        cmd = ['minimap2', '-x', 'asm5', '-c', '--secondary=no', '-t', str(threads), row['fasta'], manifest['assembly']['path']]
        rec = run(cmd, case, paf, case / f'logs/{rid}.assembly.log')
        history.append(rec)
        entries.append(dict(evidence_id='asm_' + rid, kind='assembly_paf', path=str(paf.relative_to(case)),
                            sha256=digest(paf), assembly_sha256=digest(asm), reference_id=rid,
                            reference_sha256=digest(ref), software={'minimap2': version}, command=cmd,
                            evidence_class='mapped_assembly', dependency_group='reference_' + rid,
                            coordinate_system='PAF-0-based-half-open'))
        for seqid, seq in read_fasta(ref).items():
            combined_seqs[quote(rid, safe='') + '::' + quote(seqid, safe='')] = seq
        for line in paf.read_text().splitlines():
            f = line.split('\t')
            f[5] = quote(rid, safe='') + '::' + quote(f[5], safe='')
            combined_paf.append('\t'.join(f))
    # Convenience inputs for single-reference comparators. PAF is the union of
    # independent mappings, not an independent evidence source. RagTag remaps
    # the same concatenated FASTA itself; separate-ref MAPQ is not comparable.
    write_fasta(case / 'inputs/references.combined.fa', combined_seqs)
    (case / 'evidence/references.combined.paf').write_text('\n'.join(combined_paf) + '\n')
    reads = manifest['readsets'][0]
    if digest(case / reads['path']) != reads['sha256']:
        raise ValueError('Read checksum mismatch')
    paf = case / 'evidence/reads.assembly.paf'
    cmd = ['minimap2', '-x', 'map-hifi', '-c', '--secondary=no', '-t', str(threads), manifest['assembly']['path'], reads['path']]
    history.append(run(cmd, case, paf, case / 'logs/reads.assembly.log'))
    entries.append(dict(evidence_id='reads_sim01', kind='read_paf', path=str(paf.relative_to(case)),
                        sha256=digest(paf), assembly_sha256=digest(asm), readset_id='sim01',
                        software={'minimap2': version}, command=cmd, coordinate_system='PAF-0-based-half-open',
                        evidence_class='mapped_simulated_reads', dependency_group='sim01',
                        limitation='iid substitution reads; not platform-calibrated; supplementary mappings preserved'))
    if bam:
        sam = case / 'evidence/reads.assembly.sam'
        cmd = ['minimap2', '-ax', 'map-hifi', '--secondary=no', '-t', str(threads), manifest['assembly']['path'], reads['path']]
        history.append(run(cmd, case, sam, case / 'logs/reads.sam.log'))
        bpath = case / 'evidence/reads.assembly.bam'
        subprocess.run(['samtools', 'sort', '-@', str(threads), '-o', str(bpath.relative_to(case)), str(sam.relative_to(case))], check=True, cwd=case)
        subprocess.run(['samtools', 'index', str(bpath.relative_to(case))], check=True, cwd=case)
        sam.unlink()
        entries.append(dict(evidence_id='reads_bam_sim01', kind='bam', path=str(bpath.relative_to(case)),
                            sha256=digest(bpath), assembly_sha256=digest(asm), readset_id='sim01',
                            index_path=str(bpath.relative_to(case)) + '.bai', index_sha256=digest(str(bpath) + '.bai'),
                            evidence_class='mapped_simulated_reads', dependency_group='sim01',
                            command=cmd, software={'minimap2':version,'samtools':subprocess.check_output(['samtools','--version'],text=True).splitlines()[0]}))
    manifest['evidence'] = entries
    manifest['convenience'] = dict(combined_reference_path='inputs/references.combined.fa',
                                  combined_reference_sha256=digest(case / 'inputs/references.combined.fa'),
                                  combined_paf_path='evidence/references.combined.paf',
                                  combined_paf_sha256=digest(case / 'evidence/references.combined.paf'),
                                  target_namespace='percent-encoded reference_id::sequence_id',
                                  note='Combined PAF is per-reference union, not joint mapping. Comparator remaps combined FASTA.')
    save_json(case / 'logs/evidence_commands.json', history)
    save_json(mpath, manifest)
    print(case.name, len(entries), 'mapped evidence entries')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('case', nargs='+', type=Path)
    p.add_argument('--threads', type=int, default=2)
    p.add_argument('--bam', action='store_true')
    a = p.parse_args()
    for case in a.case:
        evidence(case, a.threads, a.bam)
