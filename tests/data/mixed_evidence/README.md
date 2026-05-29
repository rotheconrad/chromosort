Mixed GFA/PAF/GAF review fixture.

The scaffold inputs in `tests/data/scaffold` place `contigA` before `contigB`.
This fixture gives their graph junction two possible intermediate nodes:
`bridge_good` and `bridge_alt`. The long-read PAF supports the contig-to-contig
junction generally, while the GAF alignments support the alternate graph branch
more strongly than the first enumerated graph path.
