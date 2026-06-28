"""AGP and component-provenance helpers for FASTA-changing outputs."""

from collections import OrderedDict
from dataclasses import dataclass
from typing import Optional

from .paths import ensure_parent_dir


@dataclass(frozen=True)
class AgpPart:
    object_name: str
    object_start: int
    object_end: int
    part_number: int
    part_type: str
    component_type: str
    component_id: str = "."
    component_start: Optional[int] = None
    component_end: Optional[int] = None
    orientation: str = "+"
    gap_length: Optional[int] = None
    gap_type: str = "."
    linkage: str = "."
    linkage_evidence: str = "."
    source: str = "."
    status: str = "."
    notes: str = "."


@dataclass(frozen=True)
class AgpFastaRecord:
    name: str
    seq: str


GAP_COMPONENT_TYPES = {"N", "U"}


def default_sidecar_path(path, suffix):
    return f"{path}{suffix}"


def parse_int(value, field, line_number):
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Malformed AGP line {line_number}: expected integer {field}, found {value!r}"
        ) from exc


def read_agp(path):
    parts = []
    with open(path) as fh:
        for line_number, line in enumerate(fh, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                raise ValueError(
                    f"Malformed AGP line {line_number}: expected 9 columns, found {len(cols)}"
                )
            object_name = cols[0]
            object_start = parse_int(cols[1], "object_start", line_number)
            object_end = parse_int(cols[2], "object_end", line_number)
            part_number = parse_int(cols[3], "part_number", line_number)
            component_type = cols[4]
            if object_start < 1 or object_end < object_start:
                raise ValueError(
                    f"Malformed AGP line {line_number}: invalid object coordinates "
                    f"{object_start}-{object_end}"
                )
            if component_type in GAP_COMPONENT_TYPES:
                gap_length = parse_int(cols[5], "gap_length", line_number)
                if gap_length < 1:
                    raise ValueError(
                        f"Malformed AGP line {line_number}: gap_length must be positive"
                    )
                parts.append(
                    AgpPart(
                        object_name=object_name,
                        object_start=object_start,
                        object_end=object_end,
                        part_number=part_number,
                        part_type="gap",
                        component_type=component_type,
                        gap_length=gap_length,
                        gap_type=cols[6],
                        linkage=cols[7],
                        linkage_evidence=cols[8],
                    )
                )
            else:
                component_start = parse_int(cols[6], "component_start", line_number)
                component_end = parse_int(cols[7], "component_end", line_number)
                if component_start < 1 or component_end < component_start:
                    raise ValueError(
                        f"Malformed AGP line {line_number}: invalid component coordinates "
                        f"{component_start}-{component_end}"
                    )
                parts.append(
                    AgpPart(
                        object_name=object_name,
                        object_start=object_start,
                        object_end=object_end,
                        part_number=part_number,
                        part_type="component",
                        component_type=component_type,
                        component_id=cols[5],
                        component_start=component_start,
                        component_end=component_end,
                        orientation=cols[8],
                    )
                )
    return parts


def group_parts_by_object(parts):
    grouped = OrderedDict()
    for part in parts:
        grouped.setdefault(part.object_name, []).append(part)
    for object_name in grouped:
        grouped[object_name].sort(key=lambda item: (item.object_start, item.part_number))
    return grouped


def component_part(
    object_name,
    object_start,
    object_end,
    part_number,
    component_id,
    component_start,
    component_end,
    orientation="+",
    component_type="W",
    source="contig",
    status=".",
    notes=".",
):
    return AgpPart(
        object_name=object_name,
        object_start=object_start,
        object_end=object_end,
        part_number=part_number,
        part_type="component",
        component_type=component_type,
        component_id=component_id,
        component_start=component_start,
        component_end=component_end,
        orientation=orientation,
        source=source,
        status=status,
        notes=notes,
    )


def gap_part(
    object_name,
    object_start,
    object_end,
    part_number,
    gap_length,
    gap_type="scaffold",
    linkage="yes",
    linkage_evidence="align_genus",
    source="gap",
    status=".",
    notes=".",
):
    return AgpPart(
        object_name=object_name,
        object_start=object_start,
        object_end=object_end,
        part_number=part_number,
        part_type="gap",
        component_type="N",
        gap_length=gap_length,
        gap_type=gap_type,
        linkage=linkage,
        linkage_evidence=linkage_evidence,
        source=source,
        status=status,
        notes=notes,
    )


def write_agp(path, parts):
    ensure_parent_dir(path)
    with open(path, "w") as out:
        out.write("##agp-version\t2.1\n")
        for part in parts:
            if part.part_type == "gap":
                row = [
                    part.object_name,
                    part.object_start,
                    part.object_end,
                    part.part_number,
                    part.component_type,
                    part.gap_length,
                    part.gap_type,
                    part.linkage,
                    part.linkage_evidence,
                ]
            else:
                row = [
                    part.object_name,
                    part.object_start,
                    part.object_end,
                    part.part_number,
                    part.component_type,
                    part.component_id,
                    part.component_start,
                    part.component_end,
                    part.orientation,
                ]
            out.write("\t".join(str(item) for item in row) + "\n")


def write_component_tsv(path, parts):
    header = [
        "object",
        "object_start",
        "object_end",
        "part_number",
        "part_type",
        "source",
        "status",
        "component_type",
        "component_id",
        "component_start",
        "component_end",
        "orientation",
        "gap_length",
        "gap_type",
        "linkage",
        "linkage_evidence",
        "notes",
    ]
    ensure_parent_dir(path)
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for part in parts:
            row = [
                part.object_name,
                part.object_start,
                part.object_end,
                part.part_number,
                part.part_type,
                part.source,
                part.status,
                part.component_type,
                part.component_id,
                "." if part.component_start is None else part.component_start,
                "." if part.component_end is None else part.component_end,
                part.orientation,
                "." if part.gap_length is None else part.gap_length,
                part.gap_type,
                part.linkage,
                part.linkage_evidence,
                part.notes,
            ]
            out.write("\t".join(str(item) for item in row) + "\n")
