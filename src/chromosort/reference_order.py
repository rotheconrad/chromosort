#!/usr/bin/env python3
"""
Order assembly contigs by their best MUMmer match to a reference genome.

This script is intended for coords files produced from a filtered nucmer delta,
for example:

    nucmer -p sample ref.fa assembly.fa
    delta-filter -i 95 -l 10000 -1 sample.delta > sample.filter
    show-coords -r -c -l sample.filter > sample.coords

For each assembly contig, the script:
  1. merges overlapping alignment intervals from show-coords,
  2. chooses the best reference chromosome by merged query coverage,
  3. keeps confident contig-to-chromosome assignments,
  4. removes contigs that duplicate already-kept reference intervals,
  5. writes an ordered FASTA named as REF_CONTIG,
  6. writes TSV reports explaining the choices.

No show-tiling output is required.
"""

import argparse
import gzip
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class FastaIndexRecord:
    name: str
    length: int
    offset: int
    line_bases: int
    line_width: int


@dataclass
class Segment:
    ref: str
    query: str
    ref_start: int
    ref_end: int
    query_start: int
    query_end: int
    len_ref: int
    len_query: int
    identity: float
    ref_length: int
    query_length: int
    orientation: str


@dataclass
class MatchMetrics:
    query: str
    ref: str
    query_length: int
    ref_length: int
    raw_ref_bp: int
    raw_query_bp: int
    merged_query_bp: int
    merged_ref_bp: int
    query_cov: float
    ref_cov: float
    avg_identity: float
    ref_start: int
    ref_end: int
    ref_midpoint: float
    orientation: str
    segment_count: int
    merged_ref_intervals: list


@dataclass
class Assignment:
    query: str
    query_length: int
    status: str
    kept: bool
    best: Optional[MatchMetrics]
    second: Optional[MatchMetrics]
    best_ref_share: float
    total_refs_matched: int
    new_name: str = "."
    order_in_ref: Optional[int] = None
    reverse_complemented: bool = False
    novel_ref_bp: Optional[int] = None
    novel_ref_frac: Optional[float] = None
    overlap_best_contig: str = "."
    overlap_best_bp: Optional[int] = None


RC_TABLE = str.maketrans(
    "ACGTRYKMSWBDHVNacgtrykmswbdhvn",
    "TGCAYRMKSWVHDBNtgcayrmkswvhdbn",
)


def parse_args(argv=None, prog=None):
    ap = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Use MUMmer show-coords output to assign assembly contigs to their "
            "best reference chromosome and write a reference-ordered FASTA."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--ref-fasta", required=True, help="Reference FASTA.")
    ap.add_argument(
        "--ref-fai",
        default=None,
        help="Reference FASTA index. Defaults to <ref-fasta>.fai when present.",
    )
    ap.add_argument(
        "-f",
        "--assembly-fasta",
        required=True,
        help="Assembly FASTA whose contigs should be ordered.",
    )
    ap.add_argument(
        "--assembly-fai",
        default=None,
        help="Assembly FASTA index. Defaults to <assembly-fasta>.fai when present.",
    )
    ap.add_argument(
        "-c",
        "--coords",
        required=True,
        help="MUMmer show-coords file, preferably produced from delta-filter output.",
    )
    ap.add_argument(
        "-o",
        "--output-prefix",
        required=True,
        help=(
            "Output prefix. Writes <prefix>.ordered.fa, <prefix>.contig_assignments.tsv, "
            "<prefix>.contig_ref_matches.tsv, <prefix>.chromosome_summary.tsv, "
            "and <prefix>.run_summary.txt."
        ),
    )
    ap.add_argument(
        "--min-aligned-bp",
        type=int,
        default=100_000,
        help="Minimum merged query-aligned bp required to keep a contig.",
    )
    ap.add_argument(
        "--min-query-cov",
        type=float,
        default=0.50,
        help="Minimum merged query coverage required to keep a contig.",
    )
    ap.add_argument(
        "--min-best-ref-share",
        type=float,
        default=0.50,
        help=(
            "Minimum fraction of a contig's total matched bp that must belong "
            "to the best reference chromosome."
        ),
    )
    ap.add_argument(
        "--min-segment-idy",
        type=float,
        default=0.0,
        help="Ignore individual show-coords rows below this percent identity.",
    )
    ap.add_argument(
        "--min-novel-ref-bp",
        type=int,
        default=50_000,
        help=(
            "During overlap filtering, keep a contig if it adds at least this "
            "many previously-uncovered reference bp."
        ),
    )
    ap.add_argument(
        "--min-novel-ref-frac",
        type=float,
        default=0.20,
        help=(
            "During overlap filtering, keep a contig if this fraction of its "
            "merged reference match is previously uncovered."
        ),
    )
    ap.add_argument(
        "--no-overlap-filter",
        action="store_true",
        help=(
            "Disable the duplicate-overlap pass. With this set, all contigs "
            "passing the basic match thresholds are written."
        ),
    )
    ap.add_argument(
        "--name-separator",
        default="_",
        help="Separator used in new FASTA IDs: REF<sep>CONTIG.",
    )
    ap.add_argument(
        "--orient-to-reference",
        action="store_true",
        help="Reverse-complement contigs whose dominant alignment is reverse-strand.",
    )
    ap.add_argument(
        "--simple-headers",
        action="store_true",
        help="Write FASTA headers containing only the new sequence ID.",
    )
    ap.add_argument(
        "--discarded-fasta",
        default=None,
        help="Optional FASTA path for contigs that were not kept.",
    )
    ap.add_argument(
        "--reports-only",
        action="store_true",
        help="Write reports but skip ordered/discarded FASTA output.",
    )
    return ap.parse_args(argv)


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def default_fai_path(fasta_path, explicit_path=None):
    if explicit_path:
        return Path(explicit_path)
    candidate = Path(str(fasta_path) + ".fai")
    return candidate if candidate.exists() else None


def read_fai(path):
    records = []
    by_name = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                raise ValueError(f"Malformed FAI line in {path}: {line.rstrip()}")
            rec = FastaIndexRecord(
                name=cols[0],
                length=int(cols[1]),
                offset=int(cols[2]),
                line_bases=int(cols[3]),
                line_width=int(cols[4]),
            )
            records.append(rec)
            by_name[rec.name] = rec
    return records, by_name


def scan_fasta_lengths(path):
    records = []
    lengths = {}
    with open_text(path) as fh:
        current = None
        length = 0
        for line in fh:
            if line.startswith(">"):
                if current is not None:
                    records.append(FastaIndexRecord(current, length, -1, 0, 0))
                    lengths[current] = records[-1]
                current = line[1:].strip().split()[0]
                length = 0
            else:
                length += len(line.strip())
        if current is not None:
            records.append(FastaIndexRecord(current, length, -1, 0, 0))
            lengths[current] = records[-1]
    return records, lengths


def read_fasta_lengths(fasta_path, fai_path=None):
    fai = default_fai_path(fasta_path, fai_path)
    if fai is not None:
        return read_fai(fai)
    return scan_fasta_lengths(fasta_path)


def parse_coords_line(line):
    cols = line.split()
    if not cols:
        return None
    try:
        if "|" in cols:
            s1 = int(cols[0])
            e1 = int(cols[1])
            s2 = int(cols[3])
            e2 = int(cols[4])
            len1 = int(cols[6])
            len2 = int(cols[7])
            identity = float(cols[9])
            lenr = int(cols[11])
            lenq = int(cols[12])
            rname = cols[-2]
            qname = cols[-1]
        else:
            s1 = int(cols[0])
            e1 = int(cols[1])
            s2 = int(cols[2])
            e2 = int(cols[3])
            len1 = int(cols[4])
            len2 = int(cols[5])
            identity = float(cols[6])
            lenr = int(cols[7])
            lenq = int(cols[8])
            rname = cols[-2]
            qname = cols[-1]
    except (ValueError, IndexError):
        return None

    orientation = "+" if (e1 - s1) * (e2 - s2) >= 0 else "-"
    return Segment(
        ref=rname,
        query=qname,
        ref_start=s1,
        ref_end=e1,
        query_start=s2,
        query_end=e2,
        len_ref=len1,
        len_query=len2,
        identity=identity,
        ref_length=lenr,
        query_length=lenq,
        orientation=orientation,
    )


def iter_coords(path, min_identity=0.0):
    in_table = False
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.strip().startswith("="):
                in_table = True
                continue

            segment = parse_coords_line(line)
            if segment is None:
                if not in_table:
                    continue
                continue

            in_table = True
            if segment.identity < min_identity:
                continue
            yield segment


def closed_to_half_open(a, b):
    lo = min(a, b)
    hi = max(a, b)
    return lo, hi + 1


def merge_intervals(intervals):
    if not intervals:
        return []
    merged = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def interval_bp(intervals):
    return sum(end - start for start, end in intervals)


def intersect_bp(left_intervals, right_intervals):
    total = 0
    i = 0
    j = 0
    left = sorted(left_intervals)
    right = sorted(right_intervals)
    while i < len(left) and j < len(right):
        start = max(left[i][0], right[j][0])
        end = min(left[i][1], right[j][1])
        if end > start:
            total += end - start
        if left[i][1] < right[j][1]:
            i += 1
        else:
            j += 1
    return total


def subtract_intervals(intervals, covered_intervals):
    if not intervals:
        return []
    if not covered_intervals:
        return list(intervals)

    covered = merge_intervals(covered_intervals)
    novel = []
    for start, end in sorted(intervals):
        cursor = start
        for cov_start, cov_end in covered:
            if cov_end <= cursor:
                continue
            if cov_start >= end:
                break
            if cov_start > cursor:
                novel.append((cursor, min(cov_start, end)))
            cursor = max(cursor, cov_end)
            if cursor >= end:
                break
        if cursor < end:
            novel.append((cursor, end))
    return novel


def build_match_metrics(coords_path, ref_lengths, query_lengths, min_identity):
    groups = defaultdict(
        lambda: {
            "raw_ref_bp": 0,
            "raw_query_bp": 0,
            "identity_bp": 0.0,
            "plus_bp": 0,
            "minus_bp": 0,
            "ref_intervals": [],
            "query_intervals": [],
            "midpoint_bp": 0.0,
            "segment_count": 0,
            "ref_length": 0,
            "query_length": 0,
        }
    )
    skipped_unknown_query = 0

    for seg in iter_coords(coords_path, min_identity):
        if seg.query not in query_lengths:
            skipped_unknown_query += 1
            continue

        key = (seg.query, seg.ref)
        data = groups[key]
        ref_interval = closed_to_half_open(seg.ref_start, seg.ref_end)
        query_interval = closed_to_half_open(seg.query_start, seg.query_end)
        ref_midpoint = (ref_interval[0] + ref_interval[1] - 1) / 2

        data["raw_ref_bp"] += seg.len_ref
        data["raw_query_bp"] += seg.len_query
        data["identity_bp"] += seg.identity * seg.len_query
        data["ref_intervals"].append(ref_interval)
        data["query_intervals"].append(query_interval)
        data["midpoint_bp"] += ref_midpoint * seg.len_query
        data["segment_count"] += 1
        data["ref_length"] = ref_lengths.get(seg.ref, FastaIndexRecord(seg.ref, seg.ref_length, -1, 0, 0)).length
        data["query_length"] = query_lengths[seg.query].length
        if seg.orientation == "+":
            data["plus_bp"] += seg.len_query
        else:
            data["minus_bp"] += seg.len_query

    metrics = []
    by_query = defaultdict(list)
    for (query, ref), data in groups.items():
        merged_ref = merge_intervals(data["ref_intervals"])
        merged_query = merge_intervals(data["query_intervals"])
        merged_ref_bp = interval_bp(merged_ref)
        merged_query_bp = interval_bp(merged_query)
        query_length = data["query_length"]
        ref_length = data["ref_length"]
        raw_query_bp = data["raw_query_bp"]
        avg_identity = data["identity_bp"] / raw_query_bp if raw_query_bp else 0.0
        ref_start = min(start for start, _ in merged_ref) if merged_ref else 0
        ref_end = max(end for _, end in merged_ref) - 1 if merged_ref else 0
        ref_midpoint = data["midpoint_bp"] / raw_query_bp if raw_query_bp else 0.0
        orientation = "-" if data["minus_bp"] > data["plus_bp"] else "+"

        match = MatchMetrics(
            query=query,
            ref=ref,
            query_length=query_length,
            ref_length=ref_length,
            raw_ref_bp=data["raw_ref_bp"],
            raw_query_bp=raw_query_bp,
            merged_query_bp=merged_query_bp,
            merged_ref_bp=merged_ref_bp,
            query_cov=merged_query_bp / query_length if query_length else 0.0,
            ref_cov=merged_ref_bp / ref_length if ref_length else 0.0,
            avg_identity=avg_identity,
            ref_start=ref_start,
            ref_end=ref_end,
            ref_midpoint=ref_midpoint,
            orientation=orientation,
            segment_count=data["segment_count"],
            merged_ref_intervals=merged_ref,
        )
        metrics.append(match)
        by_query[query].append(match)

    return metrics, by_query, skipped_unknown_query


def choose_assignments(query_records, by_query, args):
    assignments = {}
    for rec in query_records:
        matches = sorted(
            by_query.get(rec.name, []),
            key=lambda m: (m.merged_query_bp, m.raw_query_bp, m.avg_identity),
            reverse=True,
        )
        if not matches:
            assignments[rec.name] = Assignment(
                query=rec.name,
                query_length=rec.length,
                status="no_alignment",
                kept=False,
                best=None,
                second=None,
                best_ref_share=0.0,
                total_refs_matched=0,
            )
            continue

        best = matches[0]
        second = matches[1] if len(matches) > 1 else None
        total_matched_bp = sum(m.merged_query_bp for m in matches)
        best_ref_share = best.merged_query_bp / total_matched_bp if total_matched_bp else 0.0

        status = "kept"
        if best.merged_query_bp < args.min_aligned_bp:
            status = "below_min_aligned_bp"
        elif best.query_cov < args.min_query_cov:
            status = "below_min_query_cov"
        elif best_ref_share < args.min_best_ref_share:
            status = "ambiguous_ref_match"

        assignments[rec.name] = Assignment(
            query=rec.name,
            query_length=rec.length,
            status=status,
            kept=status == "kept",
            best=best,
            second=second,
            best_ref_share=best_ref_share,
            total_refs_matched=len(matches),
        )
    return assignments


def best_overlap_contig(candidate, kept_assignments):
    best_contig = "."
    best_bp = 0
    candidate_intervals = candidate.best.merged_ref_intervals
    for kept in kept_assignments:
        overlap_bp = intersect_bp(candidate_intervals, kept.best.merged_ref_intervals)
        if overlap_bp > best_bp:
            best_contig = kept.query
            best_bp = overlap_bp
    return best_contig, best_bp


def resolve_duplicate_overlaps(assignments, args):
    """
    Let the strongest contigs on each reference claim intervals first.

    A lower-ranked contig is marked duplicate_overlap when it contributes
    neither enough absolute novel reference bp nor enough fractional novel
    reference coverage compared with already-kept contigs on the same
    chromosome.
    """
    if args.no_overlap_filter:
        for assignment in assignments.values():
            if assignment.kept and assignment.best is not None:
                assignment.novel_ref_bp = assignment.best.merged_ref_bp
                assignment.novel_ref_frac = 1.0
        return

    by_ref = defaultdict(list)
    for assignment in assignments.values():
        if assignment.kept and assignment.best is not None:
            by_ref[assignment.best.ref].append(assignment)

    for ref_assignments in by_ref.values():
        ranked = sorted(
            ref_assignments,
            key=lambda a: (
                a.best.merged_ref_bp,
                a.best.merged_query_bp,
                a.best.query_cov,
                a.best.avg_identity,
                a.query_length,
            ),
            reverse=True,
        )
        covered = []
        kept_for_ref = []
        for assignment in ranked:
            ref_bp = assignment.best.merged_ref_bp
            novel_intervals = subtract_intervals(assignment.best.merged_ref_intervals, covered)
            novel_bp = interval_bp(novel_intervals)
            novel_frac = novel_bp / ref_bp if ref_bp else 0.0

            assignment.novel_ref_bp = novel_bp
            assignment.novel_ref_frac = novel_frac
            if kept_for_ref:
                overlap_contig, overlap_bp = best_overlap_contig(assignment, kept_for_ref)
                assignment.overlap_best_contig = overlap_contig
                assignment.overlap_best_bp = overlap_bp
            else:
                assignment.overlap_best_bp = 0

            adds_enough_ref_bp = novel_bp >= args.min_novel_ref_bp
            adds_enough_ref_frac = novel_frac >= args.min_novel_ref_frac
            if novel_bp == 0 or (not adds_enough_ref_bp and not adds_enough_ref_frac):
                assignment.status = "duplicate_overlap"
                assignment.kept = False
                continue

            kept_for_ref.append(assignment)
            covered = merge_intervals(covered + assignment.best.merged_ref_intervals)


def order_assignments(assignments, ref_order, separator, orient_to_reference):
    ref_rank = {name: idx for idx, name in enumerate(ref_order)}
    kept = [a for a in assignments.values() if a.kept and a.best is not None]
    kept.sort(
        key=lambda a: (
            ref_rank.get(a.best.ref, float("inf")),
            a.best.ref_start,
            a.best.ref_end,
            a.query,
        )
    )

    per_ref = Counter()
    used_names = Counter()
    for assignment in kept:
        ref = assignment.best.ref
        per_ref[ref] += 1
        assignment.order_in_ref = per_ref[ref]
        base_name = f"{ref}{separator}{assignment.query}"
        used_names[base_name] += 1
        assignment.new_name = (
            base_name
            if used_names[base_name] == 1
            else f"{base_name}{separator}{used_names[base_name]}"
        )
        assignment.reverse_complemented = (
            orient_to_reference and assignment.best.orientation == "-"
        )
    return kept


def reverse_complement(seq):
    return seq.translate(RC_TABLE)[::-1]


def iter_fasta_records(path):
    name = None
    header = None
    parts = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, header, "".join(parts)
                header = line.rstrip("\n")
                name = line[1:].strip().split()[0]
                parts = []
            else:
                parts.append(line.strip())
        if name is not None:
            yield name, header, "".join(parts)


class FastaReader:
    def __init__(self, fasta_path, fai_path=None):
        self.fasta_path = Path(fasta_path)
        self.fai_path = default_fai_path(fasta_path, fai_path)
        self.index = None
        self.seqs = None
        self.handle = None
        if self.fai_path is not None and not str(self.fasta_path).endswith(".gz"):
            _, self.index = read_fai(self.fai_path)
            self.handle = open(self.fasta_path, "rb")
        else:
            self.seqs = {name: seq for name, _, seq in iter_fasta_records(self.fasta_path)}

    def close(self):
        if self.handle is not None:
            self.handle.close()

    def fetch(self, name):
        if self.index is None:
            try:
                return self.seqs[name]
            except KeyError as exc:
                raise KeyError(f"Sequence {name!r} not found in {self.fasta_path}") from exc

        rec = self.index[name]
        self.handle.seek(rec.offset)
        remaining = rec.length
        pieces = []
        while remaining > 0:
            line = self.handle.readline()
            if not line:
                break
            seq_line = line.rstrip(b"\r\n")
            if not seq_line:
                continue
            pieces.append(seq_line.decode("ascii"))
            remaining -= len(seq_line)
        seq = "".join(pieces)
        if len(seq) != rec.length:
            raise ValueError(
                f"Expected {rec.length} bp for {name}, read {len(seq)} bp from {self.fasta_path}"
            )
        return seq


def write_wrapped(out, seq, width=80):
    for i in range(0, len(seq), width):
        out.write(seq[i : i + width] + "\n")


def fasta_header(assignment, simple_headers):
    if simple_headers:
        return assignment.new_name
    best = assignment.best
    fields = [
        assignment.new_name,
        f"original={assignment.query}",
        f"ref={best.ref}",
        f"ref_start={best.ref_start}",
        f"ref_end={best.ref_end}",
        f"orientation={best.orientation}",
        f"reverse_complemented={'yes' if assignment.reverse_complemented else 'no'}",
        f"query_cov={best.query_cov:.4f}",
        f"avg_identity={best.avg_identity:.3f}",
    ]
    return " ".join(fields)


def write_ordered_fasta(path, fasta_path, kept_assignments, assembly_fai, simple_headers):
    reader = FastaReader(fasta_path, assembly_fai)
    try:
        with open(path, "w") as out:
            for assignment in kept_assignments:
                seq = reader.fetch(assignment.query)
                if assignment.reverse_complemented:
                    seq = reverse_complement(seq)
                out.write(f">{fasta_header(assignment, simple_headers)}\n")
                write_wrapped(out, seq)
    finally:
        reader.close()


def write_discarded_fasta(path, fasta_path, assignments, assembly_fai):
    discard_ids = [name for name, assignment in assignments.items() if not assignment.kept]
    reader = FastaReader(fasta_path, assembly_fai)
    try:
        with open(path, "w") as out:
            for name in discard_ids:
                seq = reader.fetch(name)
                out.write(f">{name} status={assignments[name].status}\n")
                write_wrapped(out, seq)
    finally:
        reader.close()


def fmt(value, digits=4):
    if value is None:
        return "."
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def write_assignment_report(path, query_records, assignments):
    header = [
        "contig",
        "length",
        "status",
        "kept",
        "new_name",
        "assigned_ref",
        "order_in_ref",
        "ref_start",
        "ref_end",
        "ref_midpoint",
        "orientation",
        "reverse_complemented",
        "raw_aligned_bp",
        "merged_query_bp",
        "query_cov",
        "merged_ref_bp",
        "ref_cov",
        "avg_identity",
        "best_ref_share",
        "novel_ref_bp",
        "novel_ref_frac",
        "overlap_best_contig",
        "overlap_best_new_name",
        "overlap_best_bp",
        "second_ref",
        "second_merged_query_bp",
        "second_query_cov",
        "total_refs_matched",
    ]
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for rec in query_records:
            assignment = assignments[rec.name]
            best = assignment.best
            second = assignment.second
            overlap_assignment = assignments.get(assignment.overlap_best_contig)
            overlap_new_name = (
                overlap_assignment.new_name
                if overlap_assignment is not None and overlap_assignment.new_name != "."
                else "."
            )
            row = [
                rec.name,
                rec.length,
                assignment.status,
                "yes" if assignment.kept else "no",
                assignment.new_name,
                best.ref if best else ".",
                assignment.order_in_ref if assignment.order_in_ref is not None else ".",
                best.ref_start if best else ".",
                best.ref_end if best else ".",
                fmt(best.ref_midpoint, 1) if best else ".",
                best.orientation if best else ".",
                "yes" if assignment.reverse_complemented else "no",
                best.raw_query_bp if best else 0,
                best.merged_query_bp if best else 0,
                fmt(best.query_cov) if best else "0.0000",
                best.merged_ref_bp if best else 0,
                fmt(best.ref_cov) if best else "0.0000",
                fmt(best.avg_identity, 3) if best else "0.000",
                fmt(assignment.best_ref_share),
                assignment.novel_ref_bp if assignment.novel_ref_bp is not None else ".",
                fmt(assignment.novel_ref_frac) if assignment.novel_ref_frac is not None else ".",
                assignment.overlap_best_contig,
                overlap_new_name,
                assignment.overlap_best_bp if assignment.overlap_best_bp is not None else ".",
                second.ref if second else ".",
                second.merged_query_bp if second else ".",
                fmt(second.query_cov) if second else ".",
                assignment.total_refs_matched,
            ]
            out.write("\t".join(str(x) for x in row) + "\n")


def write_match_report(path, matches):
    header = [
        "contig",
        "ref",
        "query_length",
        "ref_length",
        "segment_count",
        "raw_query_aligned_bp",
        "raw_ref_aligned_bp",
        "merged_query_bp",
        "query_cov",
        "merged_ref_bp",
        "ref_cov",
        "ref_start",
        "ref_end",
        "ref_midpoint",
        "orientation",
        "avg_identity",
    ]
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for match in sorted(
            matches,
            key=lambda m: (m.query, -m.merged_query_bp, m.ref),
        ):
            row = [
                match.query,
                match.ref,
                match.query_length,
                match.ref_length,
                match.segment_count,
                match.raw_query_bp,
                match.raw_ref_bp,
                match.merged_query_bp,
                fmt(match.query_cov),
                match.merged_ref_bp,
                fmt(match.ref_cov),
                match.ref_start,
                match.ref_end,
                fmt(match.ref_midpoint, 1),
                match.orientation,
                fmt(match.avg_identity, 3),
            ]
            out.write("\t".join(str(x) for x in row) + "\n")


def write_chromosome_summary(path, ref_records, kept_assignments):
    by_ref = defaultdict(list)
    for assignment in kept_assignments:
        by_ref[assignment.best.ref].append(assignment)

    header = [
        "ref",
        "ref_length",
        "kept_contigs",
        "covered_ref_bp",
        "ref_cov",
        "first_ref_start",
        "last_ref_end",
        "ordered_new_names",
        "ordered_original_contigs",
    ]
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for rec in ref_records:
            assignments = by_ref.get(rec.name, [])
            intervals = []
            for assignment in assignments:
                intervals.extend(assignment.best.merged_ref_intervals)
            merged = merge_intervals(intervals)
            covered = interval_bp(merged)
            first = min((a.best.ref_start for a in assignments), default=".")
            last = max((a.best.ref_end for a in assignments), default=".")
            row = [
                rec.name,
                rec.length,
                len(assignments),
                covered,
                fmt(covered / rec.length if rec.length else 0.0),
                first,
                last,
                ",".join(a.new_name for a in assignments) if assignments else ".",
                ",".join(a.query for a in assignments) if assignments else ".",
            ]
            out.write("\t".join(str(x) for x in row) + "\n")


def write_run_summary(path, args, output_paths, ref_records, query_records, assignments, skipped_unknown_query):
    status_counts = Counter(a.status for a in assignments.values())
    kept = sum(1 for a in assignments.values() if a.kept)
    with open(path, "w") as out:
        out.write("fasta_mummer_reference_order.py\n")
        out.write("\nInputs\n")
        out.write(f"ref_fasta\t{args.ref_fasta}\n")
        out.write(f"assembly_fasta\t{args.assembly_fasta}\n")
        out.write(f"coords\t{args.coords}\n")
        out.write("\nThresholds\n")
        out.write(f"min_aligned_bp\t{args.min_aligned_bp}\n")
        out.write(f"min_query_cov\t{args.min_query_cov}\n")
        out.write(f"min_best_ref_share\t{args.min_best_ref_share}\n")
        out.write(f"min_segment_idy\t{args.min_segment_idy}\n")
        out.write(f"overlap_filter_enabled\t{not args.no_overlap_filter}\n")
        out.write(f"min_novel_ref_bp\t{args.min_novel_ref_bp}\n")
        out.write(f"min_novel_ref_frac\t{args.min_novel_ref_frac}\n")
        out.write(f"orient_to_reference\t{args.orient_to_reference}\n")
        out.write("\nOutputs\n")
        for label, path_value in output_paths.items():
            out.write(f"{label}\t{path_value}\n")
        out.write("\nSummary\n")
        out.write(f"reference_sequences\t{len(ref_records)}\n")
        out.write(f"assembly_contigs\t{len(query_records)}\n")
        out.write(f"kept_contigs\t{kept}\n")
        out.write(f"discarded_contigs\t{len(query_records) - kept}\n")
        out.write(f"coords_rows_skipped_unknown_query\t{skipped_unknown_query}\n")
        for status, count in sorted(status_counts.items()):
            out.write(f"status_{status}\t{count}\n")
        out.write(
            "\nNote\tCoverage and ordering are based on merged intervals from show-coords. "
            "This avoids counting overlapping MUMmer rows more than once. "
            "The duplicate-overlap filter rejects otherwise-good contigs that add "
            "little or no new reference coverage after better contigs claim intervals.\n"
        )


def main(argv=None, prog=None):
    args = parse_args(argv, prog=prog)
    prefix = Path(args.output_prefix)
    if prefix.parent and str(prefix.parent) != ".":
        prefix.parent.mkdir(parents=True, exist_ok=True)

    output_paths = {
        "ordered_fasta": Path(str(prefix) + ".ordered.fa"),
        "contig_assignments": Path(str(prefix) + ".contig_assignments.tsv"),
        "contig_ref_matches": Path(str(prefix) + ".contig_ref_matches.tsv"),
        "chromosome_summary": Path(str(prefix) + ".chromosome_summary.tsv"),
        "run_summary": Path(str(prefix) + ".run_summary.txt"),
    }
    if args.discarded_fasta:
        output_paths["discarded_fasta"] = Path(args.discarded_fasta)
    for output_path in output_paths.values():
        if output_path.parent and str(output_path.parent) != ".":
            output_path.parent.mkdir(parents=True, exist_ok=True)

    ref_records, ref_by_name = read_fasta_lengths(args.ref_fasta, args.ref_fai)
    query_records, query_by_name = read_fasta_lengths(args.assembly_fasta, args.assembly_fai)
    ref_lengths = {name: rec for name, rec in ref_by_name.items()}
    query_lengths = {name: rec for name, rec in query_by_name.items()}

    matches, by_query, skipped_unknown_query = build_match_metrics(
        args.coords,
        ref_lengths,
        query_lengths,
        args.min_segment_idy,
    )
    assignments = choose_assignments(query_records, by_query, args)
    resolve_duplicate_overlaps(assignments, args)
    kept_assignments = order_assignments(
        assignments,
        [rec.name for rec in ref_records],
        args.name_separator,
        args.orient_to_reference,
    )

    write_assignment_report(output_paths["contig_assignments"], query_records, assignments)
    write_match_report(output_paths["contig_ref_matches"], matches)
    write_chromosome_summary(output_paths["chromosome_summary"], ref_records, kept_assignments)
    write_run_summary(
        output_paths["run_summary"],
        args,
        output_paths,
        ref_records,
        query_records,
        assignments,
        skipped_unknown_query,
    )

    if not args.reports_only:
        write_ordered_fasta(
            output_paths["ordered_fasta"],
            args.assembly_fasta,
            kept_assignments,
            args.assembly_fai,
            args.simple_headers,
        )
        if args.discarded_fasta:
            write_discarded_fasta(
                output_paths["discarded_fasta"],
                args.assembly_fasta,
                assignments,
                args.assembly_fai,
            )

    status_counts = Counter(a.status for a in assignments.values())
    sys.stderr.write(
        f"Kept {len(kept_assignments)}/{len(query_records)} contigs across "
        f"{len(ref_records)} reference sequences.\n"
    )
    for status, count in sorted(status_counts.items()):
        sys.stderr.write(f"  {status}: {count}\n")
    sys.stderr.write(f"Wrote reports with prefix: {prefix}\n")
    if args.reports_only:
        sys.stderr.write("Skipped FASTA output because --reports-only was set.\n")
    else:
        sys.stderr.write(f"Wrote ordered FASTA: {output_paths['ordered_fasta']}\n")


if __name__ == "__main__":
    main()
