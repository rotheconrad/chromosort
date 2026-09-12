"""Small file helpers; no ChromoSort or RagTag code is imported."""
import csv
import gzip
import hashlib
import json
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def read_fasta(path):
    records = {}
    name = None
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                name = line[1:].split()[0]
                if name in records:
                    raise ValueError('Duplicate FASTA ID: ' + name)
                records[name] = []
            elif line.strip():
                if name is None:
                    raise ValueError('Not FASTA: ' + str(path))
                records[name].append(line.strip().upper())
    if not records:
        raise ValueError('Empty FASTA: ' + str(path))
    return {k: ''.join(v) for k, v in records.items()}


def write_fasta(path, records):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w') as f:
        for name, seq in records.items():
            if not seq or set(seq) - set('ACGTRYSWKMBDHVN'):
                raise ValueError('Empty/invalid DNA: ' + name)
            f.write('>' + name + '\n')
            for start in range(0, len(seq), 80):
                f.write(seq[start:start + 80] + '\n')


def rc(seq):
    return seq.translate(str.maketrans('ACGTRYSWKMBDHVN', 'TGCAYRSWMKVHDBN'))[::-1]


def save_json(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def write_tsv(path, rows, fields):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))
