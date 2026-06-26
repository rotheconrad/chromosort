---
title: FASTA And Evidence Name Matching
description: Practical checks for keeping FASTA, coords, PAF, GFA, GAF, and review-table names aligned in ChromoSort workflows.
---

# FASTA And Evidence Name Matching

Use this guide when a ChromoSort command reports missing alignments, missing
graph nodes, empty plots, stale reviewed rows, or unexpected unplaced contigs.

Most confusing input problems come from names drifting apart. ChromoSort is
conservative about this: it will not guess that two similar-looking identifiers
mean the same record if the file formats do not actually agree.

## The Core Idea

Names are part of the evidence. A PAF row whose query is `contigA` describes
`contigA` in the assembly FASTA that was aligned. A GFA link between `utg42+`
and `utg43-` describes graph segments named `utg42` and `utg43`. A reviewed
gapfill row applies only when its current scaffold, flanking contigs, and path
nodes still match.

If names drift, coordinates drift with them.

## What Must Match

| File or field | Must match | Why |
| --- | --- | --- |
| Reference FASTA IDs | Coords reference names or PAF target names | Reference axes and reference-space intervals depend on these IDs. |
| Assembly FASTA IDs | Coords query names or PAF query names | Query coordinates and contig records depend on these IDs. |
| FASTA `.fai` IDs and lengths | The exact FASTA file being used | Stale indexes can preserve old lengths after edits. |
| GFA `S` segment names | Assembly contig IDs or assignment-report names, depending on command | Graph node lookup uses these identifiers. |
| GFA `P` path or `W` walk names | Contig FASTA names for unitig-to-contig projection | Unitig coordinates can be projected only through matching paths or walks. |
| GAF path nodes | GFA `S` segment names | GAF support is counted on graph-node traversals. |
| Hi-C pair nodes | GFA `S` segment names | Contact support is summed between graph nodes. |
| `--ref-paf` query names for gapfill | GFA intermediate node names | Reference-placement support is scored for graph nodes. |
| Reviewed scaffold rows | Current scaffold, left contig, and right contig | Stale reviewed gap overrides are rejected. |
| Reviewed gapfill rows | Current scaffold, flanks, and `path_nodes` | Stale graph fills are rejected before sequence application. |

FASTA record IDs are read as the first whitespace-delimited token after `>`.
Avoid spaces in identifiers unless you are certain every upstream and downstream
tool treats them the same way.

## Preflight Checks

Run these before a long workflow:

```bash
# FASTA IDs.
grep '^>' reference.fa | head
grep '^>' assembly.fa | head

# MUMmer coords should end each data row with reference and query names.
head -40 mummer/sample.coords

# PAF query names are column 1; reference/target names are column 6.
awk 'BEGIN{FS="\t"} !/^#/ && NF>=12 {print $1, $6; exit}' paf/sample.paf

# GFA graph nodes are S records.
awk 'BEGIN{FS="\t"} $1=="S" {print $2; count++} count==10 {exit}' assembly_graph.gfa

# GAF path strings are column 6; MAPQ is column 12.
awk 'BEGIN{FS="\t"} !/^#/ && NF>=12 {print $1, $6, $12; exit}' reads_to_graph.gaf

# Hi-C pair files should have node names and integer counts.
head graph_contacts.tsv
```

For ChromoSort outputs, also inspect the assignment report:

```bash
head results/sample.contig_assignments.tsv
```

The `contig` column records original assembly IDs. The `new_name` column records
the names written to `<prefix>.ordered.fa`.

## Pattern Gallery

### Many `no_alignment` Rows

Likely cause: the alignment query names do not match the assembly FASTA passed
to ChromoSort, or the wrong assembly FASTA was aligned.

Check PAF column 1 or coords query names against:

```bash
grep '^>' assembly.fa | head
```

### Reference Names Are Missing Or Swapped

Likely cause: the reference and query were swapped during alignment, or the
reference FASTA has different IDs from the alignment target/reference names.

For normal reference ordering with minimap2, use:

```bash
minimap2 -x asm5 -c --secondary=no reference.fa assembly.fa > sample.paf
```

The assembly should be the PAF query, and the reference should be the PAF
target.

### Graph Reports Show Missing Nodes

Likely cause: GFA `S` segment names do not match the assembly names visible to
the command. This often happens after polishing, renaming, splitting, or
scaffolding the FASTA after graph export.

For scaffold and gapfill, ChromoSort tries both original contig IDs and
ChromoSort `new_name` values from the assignment report. If neither matches,
use graph evidence from the same assembly stage or keep a name map.

### Unitig GFA Overlay Is Empty

Likely cause: the dot plot is in contig coordinates, but the GFA is in
unitig-local coordinates and has no matching `P` path or `W` walk records.

Run:

```bash
chromo graph-map \
  --ctg-fasta assembly.p_ctg.fa \
  --utg-gfa assembly.p_utg.noseq.gfa \
  --output-prefix review/sample.graphmap
```

If the warning table says paths are missing, the graph can still be useful for
topology, but it cannot define query-axis unitig intervals for a contig dot
plot.

### Reviewed Plan Is Rejected As Stale

Likely cause: the FASTA, assignment report, GFA, flanking contigs, or graph path
changed after the review table was exported.

Regenerate the eval table or gapfill plan from the current inputs. Stale-row
rejection is a safety feature, not a formatting nuisance.

## Practical Review Workflow

1. Confirm the reference FASTA IDs.
2. Confirm the assembly FASTA IDs.
3. Confirm the alignment file points to those exact IDs.
4. Confirm the command is reading the FASTA stage that produced the evidence.
5. Confirm graph node names match the assembly or assignment-report names used
   by the graph-aware command.
6. For unitig evidence on contig plots, confirm GFA `P` or `W` records support
   projection.
7. For reviewed plans, regenerate the table after changing FASTA, assignments,
   graph, or path-search settings.

## Cheat Sheet

| Symptom | First check |
| --- | --- |
| Empty or sparse plot | FASTA IDs versus coords/PAF names. |
| Many `no_alignment` rows | PAF query names or coords query names versus assembly FASTA IDs. |
| Everything maps to unexpected references | Reference/query order in the aligner command. |
| Ordered FASTA missing kept assignments | `ordered.fa` and `contig_assignments.tsv` came from different sort runs. |
| Graph context missing | GFA `S` names versus original contig IDs and `new_name` values. |
| GFA overlay empty | Unitig-vs-contig coordinate mismatch or missing GFA paths/walks. |
| `--ref-paf` support is zero | PAF query names do not match GFA node names. |
| Reviewed gapfill row rejected | Current `path_nodes` no longer matches the reviewed row. |

## Common Traps

Do not assume visually similar names are equivalent. `chr01`, `Chr01`, and
`Gm01` are distinct identifiers unless you explicitly made a mapping.

Do not use unitig coordinates as contig coordinates. A unitig graph feature
needs projection before it can be compared with a contig FASTA dot plot.

Do not keep a stale `.fai` beside a rewritten FASTA. Regenerate indexes after
sequence edits.

Do not use a graph exported before renaming or polishing unless names still
match the current FASTA or assignment report.

Do not edit a reviewed TSV and then apply it after changing the graph or sorted
FASTA. Regenerate review rows for the current inputs.

## What To Look At Next In ChromoSort

- Use [Input Files]({{ '/input-files/' | relative_url }}#name-matching-rules)
  for formal name-matching rules.
- Use [Alignment Evidence And The Exact FASTA Rule]({{ '/guides/alignment-evidence/' | relative_url }})
  when the names match but the FASTA stage may not.
- Use [chromo graph-map]({{ '/commands/graph-map/' | relative_url }}) when
  hifiasm unitig evidence must be projected onto contig coordinates.
- Use [Troubleshooting]({{ '/troubleshooting/' | relative_url }}) for symptom
  checks.
