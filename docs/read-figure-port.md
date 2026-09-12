---
title: Read figure design and provenance
---

# Read figure design and provenance

The generic read renderer adapts the user-authorized, user-authored improved
read-snapshot design from `round07_ssr_read_snapshot_matrix.py`, specifically
the result-bundle copy with SHA-256
`0bcad6c14f13c856c1cedc21140926ccf73ab601c6875029f5c1079b07bed898`.
The earlier copy was not used as the layout reference. The inspected scripts
had no separate license header; the port follows the repository's MIT license
and the author's explicit instruction to adapt this work.

Retained design elements: region/feature/read models, native and reversed axis
transforms, sorting segments in display coordinates before packing, feature,
coverage, orientation and optional variant lanes, and static export. The port
uses ChromoSort drawing primitives for consistent SVG/PDF/PNG output instead
of requiring an external SVG converter.

The inspected structural-metrics companion motivated supplementary, SA and
terminal-clip summaries. The generic implementation changes those semantics:
supplementary segments are retained; unique molecules are distinct from
alignments; full endpoints determine clipping; CIGAR D/N gaps do not count as
aligned coverage; viewport cropping does not create a breakpoint. It shares
the breakpoint summary model with `longreads.py`.

No unpublished source sequences, figures, genotype lists, marker constants,
native intervals, HPC paths, or all-sample mapping jobs are distributed in this
port. Controlled fixtures exercise geometry and repeatability. Public example
figures must be regenerated from the separately versioned publication benchmark.
See [read figure semantics]({{ '/commands/reads/' | relative_url }}).
