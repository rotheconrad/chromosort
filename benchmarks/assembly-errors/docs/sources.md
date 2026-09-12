# Public source audit

Verified 12 September 2026. The [sequence registry](../sources/sequence_registry.json)
contains exact accession versions, inclusive intervals, organism/cultivar,
assembly and chromosome identities, retrieval URLs and canonical FASTA hashes.
The [error registry](../sources/source_error_registry.json) records admission
status, independent source evidence and uncertainty. No unpublished sequence
or unpublished experimental read data are used.

## Wheat scaffold-orientation reconstruction

[Zhu et al. 2021](https://doi.org/10.1111/tpj.15289),
[full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC8360199/), report DLS and
independent NLRS optical-map reconciliation of scaffold placement/orientation
and chimeras. Table S3, sheet `v2.1 vs v1.0`, row 219 assigns scaffold198 /
chr1B_scf_198 to 1B in reverse orientation; flanking scaffolds197/199 remain
forward. Its v1 interval is 90,394,774–90,538,584. The stated 143,810 bp length
differs by one from inclusive endpoint arithmetic.

Original LS992081.1 (GCA_900519105.1), positions 90,144,774–90,788,584,
provides 250 kb flanks around the interval. Corrected CM031179.1
(GCA_018294505.1), positions 96,012,942–96,994,776, supplies the guide. D001
reverse-complements the reported central interval to form its target, then
restores the original orientation for its erroneous assembly. Original bases
and gap Ns remain, so complete corrected-v2 sequence reconstruction is not
claimed. The optical maps do not establish base-precise breakpoint truth.

The compact D001 fixture retains its original frozen truth JSON. Its construction
note incorrectly says the one-base discrepancy is "masked at the gap boundary":
only non-ACGT positions are masked, not the adjacent called base. This public
generator corrects that prose; generated sequence, mapping and operational
boundaries remain identical. Use the reported inclusive endpoints and the
documented 500 bp scoring tolerance, not a one-base accuracy claim.

The [supplement download](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8360199/supplementaryFiles)
contains `TPJ-107-303-s002.xlsx`; its member SHA-256 is in the error registry.
The workbook and paper figures are not redistributed. Mapping the registered
original and guiding FASTAs with `minimap2 -x asm5 -c --secondary=no` regenerates
the supporting sequence geometry. It is not independent optical-map evidence.

## Barley Mla patterns and actual tandem copies

[Mascher et al. 2021](https://doi.org/10.1093/plcell/koab077),
[full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC8290290/), Fig. 1 describe
TRITEX collapse and a small erroneous inversion at Mla versus CCS_Canu.
The figure names scaffold_10x_111 and tig00014565; supplementary tables identify
GCA_904848325.1 / GCA_904848295.1 and independent Sanger BAC AF427791.1.
The BAC was published by [Wei et al. 2002](https://doi.org/10.1105/tpc.002238).

AF427791.1 provides 261,265 bp. Its self-alignment supports an approximately
39.7 kb tandem repeat. We use nonoverlapping labels [59,567, 99,282) and
[99,282, 138,992). These are conservative labels derived from alignments,
not experimentally localized molecular tandem-unit boundaries.

Exact historical TRITEX edits were not recovered from public contig identifiers.
The benchmark therefore injects declared orientation/deletion patterns on the
BAC background. It does not call those sequences the published erroneous
assembly. The actual repeat copies are distinct from simulated target biology.
Missing sequence remains unresolved without an available sequence source.

## Labeled wheat backgrounds and compact source reconstruction

CM031178.1/CM031179.1/CM031180.1 supply Chinese Spring 1A/1B/1D intervals.
These are actual chromosome labels, not asserted syntenic homoeologous windows.
D009's duplicate label and D010's reciprocal exchanges are controlled simulations.
They do not recover a published natural exchange or infer phasing.

The three committed fixtures contain all five canonical public intervals.
[embedded_sources.json](../sources/embedded_sources.json) maps each accession
source to its included FASTA record; the retriever restores the canonical source
header/line wrapping and verifies the original checksum. D001's input assembly
is the original v1 source sequence, so it safely supplies that canonical source.
No damaged or simulated sequence is substituted for an unmodified source.

Data use follows the [NCBI molecular-data policy](https://www.ncbi.nlm.nih.gov/home/about/policies/)
and the [distribution notice](../NOTICE.md). Additional exact historical errors
or empirical biological variants require their own accession/interval recovery
and independent evidence before admission.
