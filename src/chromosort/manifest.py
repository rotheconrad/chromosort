"""User-labeled references with separate coordinate frames and evidence identities.

All internal alignment coordinates retain the existing Segment convention
(1-based closed). The manifest's annotation intervals are 0-based half-open.
No reference multiplicity is interpreted as independent biological support.
"""

import csv
import json
from collections import defaultdict
from dataclasses import dataclass, replace
from pathlib import Path
from urllib.parse import quote
import re

from .provenance import check_digest, sha256_file, file_record


SCHEMA = "chromosort-input-v1"
ALIGNMENT_KINDS = {"assembly_paf": "paf", "assembly_coords": "coords"}


def checked_input(path, digest, label):
    if not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest):
        raise ValueError(f"A lowercase SHA-256 digest is required for {label}: {path}")
    return check_digest(path, digest, label)


def namespaced(*labels):
    return "::".join(quote(str(label), safe="").replace(".", "%2E") or "." for label in labels)


def table(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@dataclass(frozen=True)
class ReferenceSequence:
    reference_id: str
    sequence_id: str
    chromosome: str
    output_group: str
    subgenome: str = ""
    homoeolog_group: str = ""
    copy_label: str = ""
    backbone: bool = False

    @property
    def key(self):
        return namespaced(self.reference_id, self.sequence_id)

    @property
    def group(self):
        return (self.output_group, self.chromosome, self.subgenome, self.copy_label)

    @property
    def output_id(self):
        return namespaced(*self.group)


class InputBundle:
    def __init__(self, path, assembly_fasta=None):
        from .reference_order import FastaIndexRecord, scan_fasta_lengths

        self.path = Path(path).resolve()
        self.read_cache = {}
        self._input_records = None
        self.root = self.path.parent
        self.data = json.loads(self.path.read_text())
        if self.data.get("schema") != SCHEMA:
            raise ValueError(f"Expected manifest schema {SCHEMA}")
        assembly = self.data["assembly"]
        self.assembly = self.resolve(assembly["path"])
        if assembly_fasta and Path(assembly_fasta).resolve() != self.assembly:
            if sha256_file(assembly_fasta) != sha256_file(self.assembly):
                raise ValueError("--assembly-fasta differs from the manifest assembly")
            self.assembly = Path(assembly_fasta).resolve()
        self.assembly_digest = checked_input(self.assembly, assembly["sha256"], "assembly")
        self.query_records, self.query_by_name = scan_fasta_lengths(self.assembly)
        if not self.query_records or len(self.query_records) != len(self.query_by_name) or any(r.length < 1 for r in self.query_records):
            raise ValueError("Assembly FASTA must have unique nonempty records")
        self.references = self.rows("references")
        self.sequences = {}
        self.reference_files = {}
        self.reference_digests = {}
        self.native_records = {}
        for row in self.references:
            ref_id = row["reference_id"]
            if not ref_id or ref_id in self.reference_files:
                raise ValueError(f"Missing or duplicate reference_id: {ref_id!r}")
            if row.get("role", "supporting") not in {"backbone", "supporting", "alternative"}:
                raise ValueError(f"Unknown reference role for {ref_id}")
            fasta = self.resolve(row["fasta"])
            self.reference_files[ref_id] = fasta
            self.reference_digests[ref_id] = checked_input(fasta, row["sha256"], "reference")
            records, _ = scan_fasta_lengths(fasta)
            if not records or len({r.name for r in records}) != len(records) or any(r.length < 1 for r in records):
                raise ValueError(f"Reference FASTA must have unique nonempty records: {ref_id}")
            for rec in records:
                self.native_records[(ref_id, rec.name)] = rec
        for row in self.rows("reference_sequences"):
            fields = {name: row.get(name, "") or "" for name in ReferenceSequence.__dataclass_fields__}
            if not all(fields[k] for k in ("reference_id", "sequence_id", "chromosome", "output_group")):
                raise ValueError("Reference sequence requires reference_id, sequence_id, chromosome and output_group")
            raw = str(fields["backbone"]).lower()
            if raw not in {"true", "false", "yes", "no", "1", "0", ""}:
                raise ValueError(f"Invalid backbone flag: {raw}")
            fields["backbone"] = raw in {"true", "yes", "1"}
            seq = ReferenceSequence(**fields)
            if seq.key in self.sequences:
                raise ValueError(f"Duplicate reference sequence label: {seq.key}")
            if (seq.reference_id, seq.sequence_id) not in self.native_records:
                raise ValueError(f"Labeled sequence absent from reference FASTA: {seq.key}")
            self.sequences[seq.key] = seq
        if len(self.sequences) != len(self.native_records):
            raise ValueError("Every reference FASTA record must have exactly one reference_sequences row")
        self.backbones = {}
        for seq in self.sequences.values():
            if seq.backbone:
                if seq.group in self.backbones:
                    raise ValueError(f"Multiple coordinate backbones for {seq.output_id}")
                self.backbones[seq.group] = seq.key
        self.ref_records = [
            FastaIndexRecord(seq.key, self.native_records[(seq.reference_id, seq.sequence_id)].length, -1, 0, 0)
            for seq in self.sequences.values()
        ]
        self.ref_by_name = {rec.name: rec for rec in self.ref_records}
        self.evidence = self.data.get("evidence", [])
        seen = set()
        aligned_references = set()
        for entry in self.evidence:
            eid = entry["evidence_id"]
            if not eid or eid in seen:
                raise ValueError(f"Duplicate or empty evidence_id: {eid!r}")
            seen.add(eid)
            checked_input(self.resolve(entry["path"]), entry["sha256"], "evidence")
            if entry.get("index_path"):
                checked_input(self.resolve(entry["index_path"]), entry["index_sha256"], "evidence index")
            if entry["assembly_sha256"] != self.assembly_digest:
                raise ValueError(f"Evidence {eid} was generated against a different assembly stage")
            if entry["kind"] in ALIGNMENT_KINDS:
                ref_id = entry["reference_id"]
                if ref_id not in self.reference_files:
                    raise ValueError(f"Unknown reference_id in evidence: {ref_id}")
                if entry["reference_sha256"] != self.reference_digests[ref_id]:
                    raise ValueError(f"Evidence {eid} uses a different reference FASTA")
                if ref_id in aligned_references:
                    raise ValueError(f"Choose one primary alignment source for reference {ref_id}")
                aligned_references.add(ref_id)

    def resolve(self, value):
        path = Path(value)
        return (self.root / path).resolve() if not path.is_absolute() else path.resolve()

    def rows(self, field):
        value = self.data[field]
        return table(self.resolve(value)) if isinstance(value, str) else value

    def input_paths(self):
        paths = [self.path, self.assembly, *self.reference_files.values()]
        paths.extend(self.resolve(self.data[k]) for k in ("references", "reference_sequences") if isinstance(self.data[k], str))
        for entry in self.evidence:
            paths.append(self.resolve(entry["path"]))
            if entry.get("index_path"):
                paths.append(self.resolve(entry["index_path"]))
        return list(dict.fromkeys(paths))

    def input_records(self):
        if self._input_records is None:
            self._input_records = [file_record(path) for path in self.input_paths()]
        return self._input_records

    def iter_segments(self, min_identity=0, min_mapq=0, include_secondary_paf=False):
        from .reference_order import iter_alignments, inspect_paf_metadata
        for entry in self.evidence:
            fmt = ALIGNMENT_KINDS.get(entry["kind"])
            if not fmt:
                continue
            if fmt == "paf" and inspect_paf_metadata(self.resolve(entry["path"])).malformed_rows:
                raise ValueError(f"Malformed primary PAF evidence: {entry['evidence_id']}")
            for seg in iter_alignments(self.resolve(entry["path"]), fmt, min_identity, min_mapq, include_secondary_paf):
                ref = namespaced(entry["reference_id"], seg.ref)
                if ref not in self.ref_by_name or seg.query not in self.query_by_name:
                    raise ValueError(f"Alignment name mismatch in {entry['evidence_id']}: {seg.ref}/{seg.query}")
                qlen = self.query_by_name[seg.query].length
                rlen = self.ref_by_name[ref].length
                if seg.query_length != qlen or seg.ref_length != rlen:
                    raise ValueError(f"Alignment length mismatch for {seg.query}/{ref}")
                if not 0 <= seg.identity <= 100 or (seg.mapq is not None and not 0 <= seg.mapq <= 255):
                    raise ValueError(f"Invalid alignment identity or MAPQ for {seg.query}/{ref}")
                if not (1 <= min(seg.ref_start, seg.ref_end) <= max(seg.ref_start, seg.ref_end) <= rlen
                        and 1 <= min(seg.query_start, seg.query_end) <= max(seg.query_start, seg.query_end) <= qlen):
                    raise ValueError(f"Alignment interval outside FASTA for {seg.query}/{ref}")
                yield replace(seg, ref=ref)

    def correction_segments(self, **filters):
        """Use only explicit backbones; abstain at overlapping homoeolog signals.

        Supporting/alternative references remain in reports but never manufacture
        transitions by interleaving coordinate frames. An overlapping alignment
        to a different output group is withheld if identity differs by <= 2 points.
        Whole segments are withheld rather than inventing base-level projection.
        """
        grouped = defaultdict(list)
        for seg in self.iter_segments(**filters):
            if self.sequences[seg.ref].backbone:
                grouped[seg.query].append(seg)
        for segments in grouped.values():
            for seg in segments:
                lo, hi = sorted((seg.query_start, seg.query_end))
                competitors = [other for other in segments if other.ref != seg.ref
                               and max(0, min(hi, max(other.query_start, other.query_end))
                                       - max(lo, min(other.query_start, other.query_end)) + 1) >= (hi - lo + 1) * 0.5]
                if any(other.identity >= seg.identity - 2.0 for other in competitors):
                    continue
                yield seg


def add_manifest_argument(group):
    group.add_argument("--manifest", help="Checked chromosort-input-v1 JSON bundle with labeled references and evidence.")


def resolve_manifest_args(args, require_reference=False):
    if getattr(args, "manifest", None):
        bundle = InputBundle(args.manifest, getattr(args, "assembly_fasta", None))
        args.assembly_fasta = str(bundle.assembly)
        if getattr(args, "ref_fasta", None):
            raise ValueError("--ref-fasta cannot be combined with --manifest")
        return bundle
    if not getattr(args, "assembly_fasta", None):
        raise ValueError("--assembly-fasta is required without --manifest")
    if require_reference and not getattr(args, "ref_fasta", None):
        raise ValueError("--ref-fasta is required without --manifest")
    return None
