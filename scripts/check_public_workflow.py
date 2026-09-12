#!/usr/bin/env python3
"""Scripted replay diagnostic; decisions are mechanism checks, not biological review."""

import argparse
import json
from pathlib import Path

import chromosort
from chromosort import manual, workflow
from chromosort.manifest import InputBundle
from chromosort.multireference import write_tsv
from chromosort.provenance import file_record, write_json
from chromosort.reference_order import iter_fasta_records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--read-evidence-id", required=True)
    parser.add_argument("--mode", default="sensitive")
    parser.add_argument("--realign", action="store_true")
    args = parser.parse_args()
    if "heldout" in str(args.manifest).lower():
        raise ValueError("Use development inputs only for this diagnostic")
    bundle = InputBundle(args.manifest)
    out = args.output_dir.resolve()
    commands = []

    def run(module, command):
        commands.append(["chromo", module, *command])
        return {"workflow": workflow, "manual": manual}[module].main(command)

    scan = out / "review"
    run("workflow", ["scan", "--manifest", str(bundle.path), "--output-dir", str(scan),
                     "--mode", args.mode, "--read-evidence-id", args.read_evidence_id])
    original = {name: seq for name, _, seq in iter_fasta_records(bundle.assembly)}
    results = {}
    for decision in ("accept", "reject"):
        rows = workflow.read_decisions(scan / "decisions.tsv")
        for row in rows:
            row.update(decision=decision, reviewer="scripted replay diagnostic",
                       notes="Mechanism check only; no biological adjudication or benchmark truth used.")
        table = out / (decision + ".decisions.tsv")
        write_tsv(table, rows, workflow.DECISION_COLUMNS)
        applied = out / decision
        fasta = run("workflow", ["apply", "--scan-dir", str(scan), "--decisions", str(table),
                                 "--output-dir", str(applied)])
        replay = applied / "replayed.fa"
        run("manual", ["apply", "--assembly-fasta", str(bundle.assembly), "--recipe", str(applied / "recipe.json"),
                       "-o", str(replay)])
        if fasta.read_bytes() != replay.read_bytes():
            raise ValueError("Recipe replay differs from applied FASTA")
        records = {name: seq for name, _, seq in iter_fasta_records(fasta)}
        if decision == "reject" and records != original:
            raise ValueError("Rejecting every cut must preserve the original contigs")
        results[decision] = {"fasta": file_record(fasta), "replay": file_record(replay),
                             "output_records": len(records), "output_bp": sum(map(len, records.values())),
                             "byte_identical_replay": True,
                             "accounting": json.loads((applied / "apply.json").read_text())["accounting"]}
    if args.realign:
        aligned = out / "realigned"
        run("workflow", ["align", "--manifest", str(bundle.path), "--assembly-fasta",
                         str(out / "accept/reviewed.fa"), "--output-dir", str(aligned),
                         "--stage", "scripted-replay-diagnostic", "--threads", "1"])
        run("workflow", ["validate", "--manifest", str(aligned / "manifest.json"),
                         "--apply-audit", str(out / "accept/apply.json"), "--output", str(out / "validation.json")])
    write_json(out / "diagnostic.json", {"schema": "chromosort-workflow-diagnostic-v1", "chromosort_version": chromosort.__version__,
               "manifest": file_record(bundle.path), "commands": commands, "results": results,
               "interpretation": "Scripted replay/identity checks only; biological correctness and review effort not evaluated."})


if __name__ == "__main__":
    main()
