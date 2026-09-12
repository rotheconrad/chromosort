#!/usr/bin/env python3
"""Recover canonical public sources from fixtures, or audit by NCBI download."""
import argparse
import json
import sys
import tempfile
import time
import urllib.parse
import urllib.request
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from benchio import digest, read_fasta, write_fasta


def retrieve(output_root, download=False, email=None):
    if download and not email:
        raise ValueError('Provide --email for the NCBI E-utilities contact parameter')
    registry = json.loads((ROOT / 'sources/sequence_registry.json').read_text())
    embedded = json.loads((ROOT / 'sources/embedded_sources.json').read_text())
    destination = Path(output_root).resolve() / 'sources/sequences'
    if ROOT / 'fixtures' in [destination, *destination.parents] or Path(output_root).resolve() == ROOT:
        raise ValueError('Use work/ or an external output root')
    destination.mkdir(parents=True, exist_ok=True)
    for row in registry['sequences']:
        key = row['source_id']
        target = destination / (key + '.fa')
        if target.exists():
            if digest(target) != row['sha256']:
                raise ValueError('Existing source checksum mismatch: ' + key)
            continue
        with tempfile.TemporaryDirectory(dir=destination) as temp:
            staged = Path(temp) / 'canonical.fa'
            if download:
                url = row['retrieval_url'] + '&' + urllib.parse.urlencode(
                    {'tool': 'chromosort-assembly-errors-benchmark', 'email': email})
                raw = Path(temp) / 'source.fa'
                for attempt in range(3):
                    try:
                        with urllib.request.urlopen(url, timeout=45) as response:
                            raw.write_bytes(response.read())
                        records = read_fasta(raw)
                        if len(records) != 1:
                            raise ValueError('Expected exactly one NCBI FASTA record')
                        sequence = next(iter(records.values()))
                        break
                    except Exception:
                        if attempt == 2:
                            raise
                        time.sleep(attempt + 1)
                time.sleep(0.4)
            else:
                source = embedded[key]
                sequence = read_fasta(ROOT / source['path'])[source['sequence_id']]
            if len(sequence) != row['length']:
                raise ValueError('Source length mismatch: ' + key)
            write_fasta(staged, {key: sequence})
            if digest(staged) != row['sha256']:
                raise ValueError('Canonical source checksum mismatch: ' + key)
            staged.replace(target)
    print('Verified', len(registry['sequences']), 'canonical public sources')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-root', type=Path, default=ROOT / 'work')
    parser.add_argument('--download', action='store_true', help='Fetch from NCBI instead of embedded public intervals')
    parser.add_argument('--email', help='NCBI contact parameter, used only for download requests and never saved')
    args = parser.parse_args()
    retrieve(args.output_root, args.download, args.email)
