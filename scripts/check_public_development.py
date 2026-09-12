#!/usr/bin/env python3
"""Check frozen public development inputs without consulting benchmark truth.

Run with the installed candidate Python, not a mutable source PYTHONPATH.
"""

import argparse
import json
import time
from pathlib import Path

import chromosort
from chromosort.manifest import InputBundle
from chromosort.provenance import file_record, write_json
from chromosort.reference_order import main as sort_main
from chromosort.workflow import verify_outputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--development-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--software-archive", type=Path, required=True)
    args = parser.parse_args()
    if "heldout" in str(args.development_root).lower():
        raise ValueError("This diagnostic runner accepts development inputs only")
    result = {"schema": "chromosort-public-development-check-v1", "version": chromosort.__version__,
              "import_path": str(Path(chromosort.__file__).resolve()),
              "software_archive": file_record(args.software_archive), "cases": []}
    for manifest in sorted(args.development_root.glob("D*/manifest.json")):
        started = time.perf_counter()
        record = {"case": manifest.parent.name, "manifest": file_record(manifest)}
        try:
            bundle = InputBundle(manifest)
            prefix = args.output_dir / manifest.parent.name / "sort"
            command = ["--manifest", str(manifest), "--retain-all", "--orient-to-reference",
                       "--min-aligned-bp", "1000", "-o", str(prefix)]
            sort_main(command)
            fasta, agp = Path(str(prefix) + ".ordered.fa"), Path(str(prefix) + ".ordered.agp")
            record.update(command=["chromo", "sort", *command], accounting=verify_outputs(bundle.assembly, fasta, agp, True),
                          fasta=file_record(fasta), agp=file_record(agp), outcome="pass")
        except Exception as exc:
            record.update(outcome="failure", error=f"{type(exc).__name__}: {exc}")
        record["seconds"] = time.perf_counter() - started
        result["cases"].append(record)
    if not result["cases"]:
        raise ValueError("No development manifests found")
    write_json(args.output_dir / "validation.json", result)
    print(json.dumps({"cases": len(result["cases"]), "passed": sum(r["outcome"] == "pass" for r in result["cases"])}))
    if any(r["outcome"] != "pass" for r in result["cases"]):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
