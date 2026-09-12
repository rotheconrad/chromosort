# Paired development evaluation: published baseline and 0.4.0rc3

The candidate completed explicit, supplied-coordinate repairs through the
supervised workflow. **Automatic correction was unchanged:** all sixteen paired
FASTA and component-ledger outputs were byte-identical to public baseline
`dd66e3d987e775d99e612f4cdc19adec039b43e6`. Neither version completely repaired an
eight-case primary-panel assembly under either automatic policy. The baseline's
existing manual recipe API already completed the same supplied-coordinate
repairs; this result supports an improved explicit workflow, not an automatic
accuracy advantage.

The [design](README.md), [prospective protocol](protocol.json), and
[execution freeze](execution-freeze.json) precede the evaluated outputs. The
candidate source SHA-256 is
`6e090a5174ff32b301666aba5f16c0c6d3789ffbe6e65984052c7d914c96f173`;
the independently installed wheel SHA-256 is
`ff5770c61aa4b5eab4aa45cbe6028cf824674488fdb4663db6d42036f0f0bdf2`.
[Installed file identities](results/software-identities.json) distinguish the
two immutable installations. Both completed output trees were recorded in the
[paired hash ledger](results/paired-output-commitment.json) before evaluator
import. No candidate changes followed scoring.

## Primary eight-case results

Each arm processes 3,646,070 input bases, ten expected error boundaries and 25
biological preservation intervals. The following totals are the same for both
software versions. Scripted repairs are a separate supplied-coordinate arm.

| Outcome | Unedited | Automatic report | Automatic review | Scripted repair |
|---|---:|---:|---:|---:|
| Expected boundaries localized | 0/10 | 1/10 | 1/10 | 10/10 |
| Extra cut boundaries | 0 | 7 | 4 | 0 |
| Complete target structures | 0/8 | 0/8 | 0/8 | 8/8 |
| Broken preservation intervals | 0/25 | 6/25 | 3/25 | 0/25 |
| Remaining target-path adjacency discrepancies | 10 | 10 | 10 | 0 |
| Every input base retained exactly once | 8/8 | 8/8 | 8/8 | 8/8 |

All sixty recorded outputs completed: 32 automatic, eight unedited, and twenty
scripted repairs (eight primary plus two prior diagnostics per version). Every
output retained all input bases exactly once, with no missing callable target
sequence or extra target copies. Retention alone did not establish structural
correctness. See [per-case metrics](results/summary.tsv),
[full scores](results/scores.jsonl), and [paired equality checks](results/evaluation-summary.json).

## Mixed contigs and evidence limits

For all D101–D106 inputs, smoothing proposed a cut at **320010**, near a simulated
biological exchange edge, and did not propose the true misjoin at **200000**.
The report policy applies that extra cut. Review defers the plan when sufficient
read continuity crosses the proposed biological-edge cut.

| Case | Evidence condition | Read names spanning proposed cut 320010 | Review outcome | Read names spanning supplied true join 200000 |
|---|---|---:|---|---:|
| D101 | 12×, 10 kb reads | 11 | Defers plan; preserves the interval | 0 |
| D102 | 1×, 10 kb reads | 1 | Applies extra cut | 0 |
| D103 | 12×, 600 bp reads | 0 | Applies extra cut | 0 |
| D104 | CIGAR removed | Unavailable (480 rows) | Applies extra cut | Unavailable |
| D105 | One chimera, 20 repeated rows | 11 | Defers plan | 1 |
| D106 | One chimera, two renamed rows | 11 | Defers plan | 2 |

The 1× condition has one spanning name, below the two-name threshold. The
600 bp reads cannot provide two 1 kb anchors. Missing CIGAR is explicitly
`unavailable_cigar` in the continuity audit, rather than evidence of an error;
the generic supplied-position figure counter still reports zero and must be
interpreted with its CIGAR-availability field. Read absence remains neutral,
while the reference-based proposal can still produce an incorrect cut.

The added D105/D106 chimera crosses the **true join**, not the automatically
proposed cut. Twenty identical rows count as one read name; renaming a shared
parent gives two names. This exposes the identifier-based independence limit,
but it does **not** change the automatic plan in this panel. Likewise, these
cases cannot demonstrate suppression of an otherwise useful mixed-contig cut:
the true cut was never proposed. Independent software geometry tests address
selective-plan accounting; they do not replace this observed proposal failure.

Exact [continuity audits](results/continuity.tsv) and
[supplied-position counters](results/supplied-boundary-support.tsv) retain the
positions, filters, evidence hashes and availability fields.

## Orientation and inspection

D107's reported source interval has operational boundaries 150000 and 796070.
Both automatic versions cut at 153089 and 796109. The latter is within the
frozen 500 bp tolerance; the former is 3089 bp away and receives no localization
credit. The original orientation remains unrepaired. The original sequence has
21,460 masked bases; no simulated reads appear at the first supplied boundary.
That absence cannot validate the historical orientation decision.

D108 receives no automatic cut and no matching inspection interval around its
injected 8 kb error. Supplying the declared reversal explicitly repairs its
sequence while preserving the separate simulated biological inversion. Linking
that edit to an existing `keep` event does not mean the scan discovered the
error. The mapped 8 kb reverse block is present in its raw PAF; it falls
below the 10 kb correction-segment floor. D012's mapped 25 kb reverse block
passes that floor but is smoothed out by the fixed breakpoint penalty. Gap
inspection unions alignment coverage regardless of strand, so a mapped reverse
block is not itself an uncovered interval. These are post-pass stage
interpretations; no parameters or candidate outputs were changed.

The prior D003 diagnostic has an inspection interval `[151002,158998)` covering
both error edges within tolerance. D012 has no matching error inspection
interval; its visible inspection interval concerns the simulated biological
inversion instead. These misses are retained in
[the inspection record](results/inspection-and-replay.json). All scan proposals
start pending, and inspection alone changes no sequence.

## Explicit repairs and portability

Both versions restore the complete target structure in **8/8 primary cases** and
**2/2 prior diagnostics** using the supplied edit coordinates. The candidate
successfully previews ten plans without sequence output, applies the reviewed
edits, reproduces each output byte-for-byte through manual recipe replay, and
validates exact input accounting against fresh reference alignments in all ten.
That validation confirms identity and reconstruction; biological validation
remains unresolved by those checks alone.

The baseline uses its existing manual API; its arm does not perform the
candidate's plan-preview and apply-audit workflow. Therefore these runs do not
measure equal-workload speed, reviewer effort, or human error rates. Local
workflow subprocess elapsed sums were about 4.5 seconds for the baseline and
12.4 seconds for the candidate, which executes more stages. Mapping uses two
threads, no network during regeneration, no HPC or shared CI. Exact descriptive
records are in [resource use](results/resource-use.json).

The [small artifact index](results/artifact-index.json) contains exact output
AGP/component ledgers, automatic reports and continuity tables, plus portable
manual recipes. Exported candidate recipes omit only the advisory absolute
source path; source SHA-256 and edit geometry remain unchanged. Original and
exported recipe hashes are recorded. Full sequence outputs, figures, reads,
software archives and logs are regenerable local artifacts and are not bundled.

[Portable verification](results/portable-score-verification.json) reconstructed
all sixty output sequences from these public provenance ledgers and reproduced
every metric and ledger hash. It does not recreate tool-specific FASTA header or
wrapping bytes; original byte hashes remain separately recorded. The first
fourteen case definitions, source identities, scorer, fixture bytes and
historical score records remain unchanged.

This bounded development panel adds new loci from one already represented wheat
genome. It is not independent species-level validation or an empirical read
study. A future development cycle would need prospectively specified proposal
coverage and independent evidence tests before claiming broader correction or
inspection performance.
