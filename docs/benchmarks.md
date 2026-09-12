---
title: Public benchmark inputs and validation
---

# Public benchmark inputs and validation

The [assembly-error benchmark](https://github.com/rotheconrad/chromosort/tree/main/benchmarks/assembly-errors)
contains compact D001, D003 and D010 development fixtures, all 14 development
case definitions, retrieval/generation scripts and independent scoring methods.
The fixtures use public reference sequence and controlled modifications; they
contain no unpublished project assemblies. Read the benchmark's source and
reproduction documentation before interpreting the examples.

Use `fixtures/development/D003/manifest.json` from that directory with
`chromo workflow scan`; select `--read-evidence-id reads_sim01`. These manifests
resolve relative paths from their own directory. The adjacent truth records
are for independent scoring and must be withheld from curation and assisted
review. Larger generated inputs and results belong in a separate working
directory. Benchmark data is available from the repository, separately from
the installation packages.

## What the development checks establish

The frozen rc1 implementation passed 151 synthetic unit tests. The frozen rc2
implementation passed 160, including continuity deferral and two-edge gap
inspection. The current public-readiness snapshot passes 162 tests, adding
invalid read-target contract cases. Source/wheel consistency, fresh installation, noarch conda builds,
CLI entry points, static figures and deterministic recipe replay were checked.
These establish implementation and replay behavior, not biological accuracy.
Interactive dashboard operation still requires a separate browser review;
static exports and JavaScript syntax checks do not establish that behavior.

Independent development evaluation preserved every input base exactly once in
all 56 rc1 outputs and all 84 paired rc2 outputs. Under rc2's default review
policy, conservative/comprehensive/sensitive modes localized 1/3/3 of nine
scored boundaries, with no unmatched new cuts in this development set. Most
errors therefore remained unresolved. Report-policy FASTAs were byte-identical
to rc1 in all 42 paired correction runs. These development observations do not
estimate held-out accuracy or guarantee safe correction.

The [full results and limitations](https://github.com/rotheconrad/chromosort/blob/main/benchmarks/assembly-errors/docs/results.md)
separate sequence accounting, boundary localization, biological preservation
and unresolved orientation events. The [review policy]({{ '/review-refinement/' | relative_url }})
explains why a continuity conflict defers a whole contig plan and why absent
spanning reads are neutral. All refinement choices used development evidence;
held-out inputs and truth were excluded from software refinement.

Frozen candidate archives retain their original checksums. The current Git
snapshot includes later input-contract validation and portability changes;
identify a run by commit and artifact checksum as well as version. See
[release preparation]({{ '/release-preparation/' | relative_url }}) for builds
that preserve the frozen archives, and [workflow diagnostics]({{ '/commands/workflow/' | relative_url }})
for input-only accounting and replay checks.
