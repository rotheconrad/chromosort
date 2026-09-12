# Controlled labeled-reference and read geometry

This fixture is synthetic mechanism-test data, not a realistic genome simulation
or biological performance benchmark. It has two user-labeled subgenomes, an
alternative reference for one group, repeated chromosome IDs, two preserved
copies, an ambiguous contig, an unplaced contig, and a two-chromosome chimera.
The short SAM fixture contains a spanning read, distinct split/supplementary
segments, a CIGAR deletion, clipping, secondary/low-MAPQ exclusions and one tiny
alignment that must not count as two split segments.

All paths in `manifest.json` are relative to this folder. Regenerate geometry
with the `make_bundle` / `add_read_fixture` functions in the corresponding test
modules. No aligner, proprietary sample data or held-out benchmark truth was
used to construct this fixture.

```bash
chromo sort --manifest tests/data/labeled_references/manifest.json \
  --min-aligned-bp 1 --retain-all -o results/labeled/sort
chromo reads --manifest tests/data/labeled_references/manifest.json \
  --evidence-id reads --contig copyA --start 20 --end 80 --breakpoint 50 \
  --min-anchor-bp 10 --read-window-bp 10 --bin-bp 5 --reverse \
  -o results/labeled/reads
```
