# Graph Gotchas Synthetic Fixture

This fixture is intentionally tiny and hand-authored so ChromoSort graph-aware
features can be developed without external aligners.

Files:

- `ref.fa`: toy two-chromosome reference.
- `assembly.fa`: toy contig/unitig sequences.
- `unitigs.gfa`: assembly graph with sequences, orientations, branches, a cycle,
  one segment with only `LN:i`, and link overlaps.
- `unitig_to_ref.paf`: minimap2-like unitig-to-reference placements.
- `reads_to_graph.gaf`: GraphAligner-like long-read paths through the graph.
- `hic_pairs.tsv`: simple node-pair contact counts.
- `fill_ordered.fa`: tiny two-flank ordered FASTA for fill walkthroughs.
- `fill_assignments.tsv`: matching two-row assignment table for
  `fill_ordered.fa`.
- `expected_gap_paths.tsv`: expected high-level path labels for future candidate
  ranking tests.

Scenarios covered:

- `confident_gap_path`: `left+ -> bridge_good+ -> right+` has clean GFA
  adjacency plus stronger PAF, GAF, and Hi-C support.
- `ambiguous_branch`: `bridge_alt` also connects left to right but has weaker
  placement and support.
- `cycle_guard`: `bridge_alt+ -> bridge_alt-` creates a cycle for path-search
  depth-limit tests.
- `orientation_specific`: `right- -> reverse_only-` should only match in that
  orientation.
- `disconnected_mapped_node`: `isolated` has a PAF placement but no graph path
  from the chr1 gap.
- `repeat_or_duplicate_warning`: `repeat_shared` and `duplicate_fragment` create
  branches that future ranking/reporting code should mark cautiously.
