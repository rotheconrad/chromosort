"""Format-aware read adapters; SAM detail is never invented for PAF input."""

import re
import math
import subprocess
from collections import Counter
from dataclasses import dataclass

from .longreads import LongReadAlignment, LongReadEvidence
from .reference_order import open_text, parse_paf_line, merge_intervals


CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


@dataclass
class ReadInput:
    evidence: LongReadEvidence
    counts: dict
    format: str
    detail: str


def cigar_parts(cigar):
    parts = [(int(n), op) for n, op in CIGAR.findall(cigar)]
    if not parts or "".join(f"{n}{op}" for n, op in parts) != cigar or any(n < 1 for n, _ in parts):
        raise ValueError(f"Invalid CIGAR: {cigar!r}")
    return parts


def cigar_geometry(cigar, start):
    parts = cigar_parts(cigar)
    position = start
    blocks = []
    for length, op in parts:
        if op in "M=X":
            blocks.append((position - 1, position + length - 1))
        if op in "MDN=X":
            position += length
    left = sum(n for n, op in parts[:next((i for i, (_, op) in enumerate(parts) if op not in "SH"), len(parts))])
    right = 0
    for n, op in reversed(parts):
        if op not in "SH":
            break
        right += n
    return position - 1, tuple((lo + 1, hi) for lo, hi in merge_intervals(blocks)), left, right


def paf_alignment(line):
    seg = parse_paf_line(line)
    if seg is None:
        raise ValueError("Malformed read PAF row")
    cols = line.rstrip().split("\t")
    tags = {c.split(":", 2)[0]: c.split(":", 2)[2] for c in cols[12:] if c.count(":") >= 2}
    cigar = tags.get("cg")
    blocks = None
    if cigar:
        end, blocks, _, _ = cigar_geometry(cigar, seg.ref_start)
        if end != seg.ref_end:
            raise ValueError(f"PAF CIGAR reference span disagrees for read {seg.query}")
        if sum(n for n, op in cigar_parts(cigar) if op in "MI=X") != seg.query_end - seg.query_start + 1:
            raise ValueError(f"PAF CIGAR query span disagrees for read {seg.query}")
    return LongReadAlignment(seg.query, seg.query_length, seg.query_start, seg.query_end,
                             seg.ref, seg.ref_length, seg.ref_start, seg.ref_end,
                             seg.orientation, seg.identity, seg.mapq, seg.is_secondary,
                             cigar=cigar, aligned_blocks=blocks)


def sam_alignment(line, lengths):
    cols = line.rstrip("\n").split("\t")
    if len(cols) < 11:
        raise ValueError("Malformed SAM alignment")
    flag = int(cols[1])
    if flag & 4:
        return None
    if cols[2] not in lengths:
        raise ValueError(f"SAM contig lacks @SQ length: {cols[2]}")
    start = int(cols[3])
    end, blocks, left_clip, right_clip = cigar_geometry(cols[5], start)
    parts = cigar_parts(cols[5])
    aligned_query = sum(n for n, op in parts if op in "MI=X")
    read_len = aligned_query + left_clip + right_clip
    qstart, qend = left_clip + 1, left_clip + aligned_query
    if flag & 16:
        qstart, qend = read_len - qend + 1, read_len - qstart + 1
    tags = {c.split(":", 2)[0]: c.split(":", 2)[2] for c in cols[11:] if c.count(":") >= 2}
    denominator = sum(n for n, op in parts if op in "MID=X")
    identity = max(0, 100 * (1 - int(tags["NM"]) / denominator)) if denominator and "NM" in tags else float("nan")
    return LongReadAlignment(cols[0], read_len, qstart, qend, cols[2], lengths[cols[2]],
                             start, end, "-" if flag & 16 else "+", identity, int(cols[4]),
                             bool(flag & 256), bool(flag & 2048), cols[5], blocks,
                             left_clip, right_clip, "SA" in tags)


def validate_read_target(contig, length, targets):
    """Check the declared assembly before filtering away incompatible evidence."""
    if targets is None:
        return
    if contig not in targets:
        raise ValueError(f"Read evidence target absent from manifest assembly: {contig}")
    if length != targets[contig].length:
        raise ValueError(f"Read evidence target length differs from manifest assembly: {contig}")


def load_reads(path, fmt="paf", min_mapq=0, min_identity=0, include_secondary=False,
               samtools="samtools", region=None, targets=None):
    if min_mapq < 0 or min_mapq > 254 or not 0 <= min_identity <= 100:
        raise ValueError("Read filters require MAPQ 0–254 and identity 0–100")
    if fmt not in {"paf", "sam", "bam"}:
        raise ValueError("Supported read formats: paf, sam, bam; CRAM requires a separate verified reference adapter")
    counts = Counter()
    lengths, alignments = {}, []
    process = None
    if fmt == "bam":
        command = [samtools, "view", "-h", str(path)]
        if region:
            command.append(region)
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        handle = process.stdout
    else:
        handle = open_text(path)
    try:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            if fmt != "paf" and line.startswith("@"):
                if line.startswith("@SQ\t"):
                    fields = dict(c.split(":", 1) for c in line.rstrip().split("\t")[1:])
                    lengths[fields["SN"]] = int(fields["LN"])
                    validate_read_target(fields["SN"], lengths[fields["SN"]], targets)
                continue
            counts["rows_seen"] += 1
            aln = paf_alignment(line) if fmt == "paf" else sam_alignment(line, lengths)
            if aln is None:
                counts["unmapped_excluded"] += 1
                continue
            validate_read_target(aln.contig, aln.contig_length, targets)
            if not (1 <= aln.contig_start <= aln.contig_end <= aln.contig_length and
                    1 <= aln.read_start <= aln.read_end <= aln.read_length):
                raise ValueError(f"Read alignment outside sequence bounds: {aln.read}")
            if aln.is_secondary and not include_secondary:
                counts["secondary_excluded"] += 1
                continue
            if aln.mapq is None or aln.mapq == 255:
                counts["mapq_unavailable"] += 1
                if min_mapq:
                    counts["mapq_excluded"] += 1
                    continue
            elif aln.mapq < min_mapq:
                counts["mapq_excluded"] += 1
                continue
            if math.isnan(aln.identity):
                counts["identity_unavailable"] += 1
            if aln.identity < min_identity or (min_identity and math.isnan(aln.identity)):
                counts["identity_excluded"] += 1
                continue
            alignments.append(aln)
            counts["alignments_retained"] += 1
            counts["supplementary_retained"] += bool(aln.is_supplementary)
            counts["sa_tag_retained"] += aln.has_sa
            counts["cigar_unavailable"] += aln.cigar is None
        if process:
            error = process.stderr.read()
            if process.wait() != 0:
                raise ValueError(f"samtools failed: {error.strip()}")
    finally:
        handle.close()
        if process:
            process.stderr.close()
            if process.poll() is None:
                process.terminate()
                process.wait()
    counts["molecules_retained"] = len({a.read for a in alignments})
    return ReadInput(LongReadEvidence.from_alignments(alignments), dict(counts), fmt,
                     "alignment-span" if any(a.cigar is None for a in alignments) else "aligned-block")
