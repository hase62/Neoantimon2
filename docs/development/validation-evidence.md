# Versioned validation evidence

[README](../../README.md) · [Reproduction](../user/reproducibility.md) ·
[Scientific limitations](../user/methods-and-limitations.md)

This page distinguishes historical scientific/software evidence from fresh
integration testing. The historical inventory table is transcribed from machine-readable primary
receipts whose exact byte hashes and revision anchors are in the
[compact inventory](validation-evidence.json). They are not new runs performed
during README/repository integration. A package test assertion is not an
independent molecular case, patient or experimentally tested neoantigen.

## Evidence inventory

| Evidence and revision | Measured scope/result | Interpretation and limitations |
| --- | --- | --- |
| Independent synthetic molecular oracles, `35eeed3` | 54 cases attempted, 51 literal expected-outcome matches, three qualified cases, zero execution errors | Includes legitimate zero-output/rejection cases and low-level molecular calls; not 54 patient samples or 54 public end-to-end analyses |
| HCC1395 saved scientific graph, `11c6816` | 41,072 source variants/Events, 244,191 Origins, 110,030 pMHC candidates and 949,015 candidate–Origin links | Genuine SEQC2 HCC1395 cancer-cell-line data; every source key accounted for; generated candidates are not all evaluated binders |
| Full saved-result recovery | All Origins, candidates, links, VAF/gene and Grade/rank projection checks passed | Original driver SIGTERM -15, controller 143 remains historical fact. Scientific checkpoint completion and later successful read-only recovery are separate. No second full independent biological rerun is claimed |
| Real MHCflurry 2.0.6 BA pilot | 75 genuine HCC1395 peptide–HLA queries; 150 imported measurements (affinity and percentile) | Genuine inference and registered Converter/import verification; public numeric projection has a documented 15-significant-digit boundary; original score strings retained |
| Real BigMHC EL/IM, source `199e890fa69a63646730bcecea0a520e08ca4f02` | 75 evaluated/imported rows per channel, five of five exact HCC1395 alleles, no fuzzy substitution, no failed/unsupported/missing rows | Academic-use pretrained inference; separate EL and IM tasks; demonstrates execution/interoperability, not biological predictive accuracy |
| Extended genuine pVACseq comparison, pVACtools 7.1.5 / VEP 116.3 | 24 variants selected (10 missense SNVs, six in-frame deletions, eight frameshifts); 164 exact shared transcript contexts, equal mutant proteins in all 164; four mutant contexts unavailable | No eligible authentic in-frame insertions in the matched GT/input intersection. Not complete mutation-class coverage. Eligibility/window and WT-mapping policies remain separate from molecular disagreement |
| Peptide-output comparison in that cohort | 18,001 Neo context-peptide observations, 17,826 pVACseq observations, 17,826 shared; 135 self-identical-output-policy and 40 eligibility differences; 2,363 WT-pair policy differences | Counts are context-scoped, not independent neoantigens. A larger emitted set is not better accuracy; no unexplained mutant-protein discrepancy was observed in the exact shared scope |
| Experimental Yield/reporting saved-output demonstration, `11c6816` plus explicitly hashed uncommitted features | 15 unique peptides, 75 pMHCs, 36 unique Origins, 180 candidate–Origin links, 300 preserved measurements; expected positives 0.5346 and independence probability 0.4208215262598308 | Reuses prior real predictions. IMPROVE-derived scenario rates on BigMHC-IM ordering are cross-model assumptions, not calibration or immunogenicity validation. All five alternative HLAs retained |
| Unfiltered report of the same small five-variant project | Five Events, 373 sequence units, 3,696 Origins, 1,450 pMHCs, 14,400 links | Small-project reporting completeness, including Origins without candidates. Not a rerun of full HCC1395 or a new prediction experiment |

The three qualified oracle cases remain explicit: an invalid cross-intron genomic
adapter (`indel_exon_boundary`), unknown actual germline background
(`germline_somatic_unknown`), and unsupported stop-loss extension boundary
(`stoploss_extension_snv`). They are not silently promoted to unconditional PASS.
The underlying historical inventory contains 49 primary, two secondary and three
fusion cases, plus eight auxiliary comparison rows; these denominators are distinct.

## Current molecular-repair integration

[PR #9](https://github.com/hase62/Neoantimon2/pull/9) integrated the bounded repairs
at merge `0fa9af969e3a0e812c2b4ff5efbbc9c1822d858a`, with tested candidate
`2aa579df884b3f5a5141cef4a562823748614bda`. Its public PR records 2,992 related
source assertions, 12,491 installed assertions with one expected repository-only
example skip, and normal `R CMD check --no-manual` with zero errors, warnings and
notes. The unchanged repository-only example separately passed 12 assertions.
The recorded R environment was 4.5.1; a testthat build-under-R-4.5.2 startup warning
was retained separately from assertion/check warning counts. These are historical
integration receipts, not newly executed tests during a documentation update,
CI success or independent experimental validation.

The same integration recorded 2,004 fixed-value checks in an 18-condition
comparison, six additional conditions / 946 checks and R09 initial/fresh query
checks of 45 / 49. The current-source side executed 18 conditions (17 completed
projects and one intentional invalid-Sec input rejection); the earlier baseline
side reused its separately recorded executions. These are condition/check units,
not independent patients. The first integration harness rejected the correctly
installed 22-API source because it required 20 exports; a source-specific API
preflight adapter resolved that technical mismatch without changing biological
inputs, control values or expected sequences. Initial failures were retained.

### Retention and later fixed-source replay

The earlier 18-condition / 2,004-check figures above are historical receipt
claims. Their complete original replay inputs, independent expectations and driver
have not been recovered; this is not a presently executable public replay scope.
The separately restored six-condition candidate (946 recorded checks) and R09
candidate retain their fixed 20-API source identities and prior execution records.

MR11 v002 is a separate study. Its preserved candidate and pre-run freeze record
nine sequence expectations, one unresolved-contract condition and one planned
input rejection (ten exit-zero cases, one expected exit-one case), 2,500 checks
and four separate R readers. Initial harness failures remain in the historical
record. The [fixed paper guide](https://github.com/hase62/Neoantimon2/blob/4234f653e9572ffbe5db6237f30fd5a89fc03e4f/docs/reproducibility/paper/README.md)
provides the public MR11, six-condition and R09 replay procedures. The
development-version documentation update did not rerun those studies.

R09 used controlled downstream input, native pVACtools binding/aggregate outputs
and external queries. At that tested boundary, the available pVACtools full table
and detailed JSON supported equivalent information queries. This does not show
that another tool loses provenance, nor a full pVACseq reconstruction, predictor
inference or clinical advantage. See [repair behavior](../developer/bounded-repair.md)
and [public replay scope](../user/reproducibility.md#current-main-and-additional-replay-materials).
The historical evidence inventory above has not been reassigned to this revision.

## Source and experimental boundaries

The HCC1395 callset is SEQC2 v1.2.1 (39,447 SNV plus 1,625 indel source records),
GRCh38.p14 with GENCODE v50. Source URLs:
[SEQC2 release](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest/)
and [HCC1395 identity/HLA](https://www.cellosaurus.org/CVCL_1249).
Archival filenames, actual input/reference hashes and full recipes remain in the
original retained scientific evidence. A mutable latest URL alone is not an
archival identifier. Truth-VCF pooled INFO/TVAF is a source-reported estimate,
not AD/BAM counts or fabricated GT.

The full run had explicit Class I 8–12mers and five verified Class I alleles.
Class II 15-mer Origins were retained, but no Class II HLA candidates or genuine
Class II predictor evaluation were claimed for that sample. All generated full-run
candidates were unevaluated in that workflow; model execution was a separate
bounded pilot. Later BA/EL/IM successes supersede earlier environment/license
blockers; those old failure logs remain historical evidence.

The original Neoantimon is a historical software baseline, not an independent
biological oracle. See [Hasegawa et al. (2020)](https://doi.org/10.1093/bioinformatics/btaa616).
No reconstruction match, predictor execution, converter pass or software regression
establishes sensitivity, specificity, clinical utility or experimentally confirmed
T-cell recognition. Native A–E Grades, A–D/B1/B2 evidence and numerical assumptions
have separate meanings.

## Reproducing and extending the evidence

The public source includes automated synthetic fixtures and installed examples.
Use the [package validation commands](../user/reproducibility.md), not a new full
HCC1395 run, to check an integration. Historical test counts are versioned receipts; see the
[recorded repair integration](#current-molecular-repair-integration) for its
source/installed/check scope, not a new run during documentation changes.

Large references, licensed models, original scientific projects and experimental
reference databases are intentionally outside Git. The local audit locators in
the inventory identify owner-held receipts, not publicly downloadable archives.
The compact inventory is not a substitute for those originals. The [public fixed snapshot](../user/reproducibility.md#snapshot-saved-result-reproduction)
provides the scoped eight-table reaggregation. The [fixed paper replay guide](https://github.com/hase62/Neoantimon2/blob/4234f653e9572ffbe5db6237f30fd5a89fc03e4f/docs/reproducibility/paper/README.md) and the [fixed candidate workflow guide](https://github.com/hase62/Neoantimon2/blob/e88934784ae8247a3ac3be2ca3e43f9ad9d4a969/docs/user/candidate-review-and-reference-percentiles.md#reproduce-from-a-fixed-checkout)
describe separate execution routes. The original full scientific projects and
private audit receipts remain distinct from those replays.
