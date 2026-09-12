# Development case construction

Definitions and random seeds are in [the configuration](../configs/development.json).
They are public development parameters. All cases are built by the same
[generator](../scripts/generate.py), which creates target sequence/read origins
separately from the corrupted assembly and its input-to-target provenance map.

| Case | Construction | What correct handling must preserve or acknowledge |
| --- | --- | --- |
| D001 | Wheat scaffold198 orientation reconstructed from the reported table | Original v1 sequence/N gaps and optical-map boundary uncertainty; not complete v2 reconstruction |
| D002 | Clean variation-bearing barley target | Simulated biological inversion and empirical tandem copies; fragmented alignments need not be errors |
| D003 | Injected 8 kb orientation error | Separate simulated target inversion must remain intact |
| D004 | Deleted native-copy interval plus orientation error | Two different defects; absent sequence cannot be recovered by placement |
| D005 | Deleted native-copy interval only | Report missing copy content without inventing recovered bases |
| D006 | Two input fragments overlap by 10 kb | Input retain-all and correct biological de-duplication are different metrics |
| D007 | Native repeat copies in separate fragments | Similar terminal sequence does not justify deleting a distinct copy |
| D008 | 5 kb interval absent from guide and input fragments | Biological target contains sequence unavailable to placement-only output |
| D009 | Identical candidate reference with incompatible group label | Explicit ambiguity and retained copies; not a biological equivalence claim |
| D010 | Two reciprocal A/B target-tail exchanges | Simulated target chromosome membership differs from local reference ancestry |
| D011 | Concatenated A/D input plus separate B | One interchromosomal misjoin, all original copies accounted for |
| D012 | 25 kb orientation error, 0.2% target divergence, 6× requested reads | Orientation-only injected error on a native-repeat background, not an injected collapse |
| D013 | Guiding reference lacks one labeled repeat copy | Reference overlap cannot erase a real target copy |
| D014 | 15 kb target-specific fragment absent from guide | Retain unplaced sequence and biological continuity |

The barley orientation/deletion coordinates are documented-pattern simulations,
not recovered historical TRITEX edits. Target inversions, exchanges, divergence
and reference deletions are declared simulations. The BAC tandem-repeat
background is empirical, with conservative labeled copy intervals rather than
proven molecular unit endpoints. Wheat A/B/D labels are supplied chromosome
identities; the selected windows are not asserted syntenic homoeologous loci.

Reads have fixed length, iid substitutions, equal strand probability and a
constant Q30 model character. Starts are sampled uniformly within eligible
ACGT runs; no reads cross baseline Ns. Requested depth is measured over eligible
runs and is not uniform observed base coverage. Counts/settings are stored in
each truth record. Gzip output bytes can vary with zlib; decoded FASTQ identity
is checked separately.

Source loci are shared across scenarios. Cases, boundaries and reads are not
independent sampled genomes. These definitions permit reproducible mechanism
tests; they do not support population-wide confidence intervals.
