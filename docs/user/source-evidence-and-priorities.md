# Source-file evidence and numerical priorities (unreleased)

These additions exist on the local development branch, **not immutable v0.3.0**.
The 16 released function signatures and somatic Grade/rank algorithms are retained;
`neo_priorities()` is one additional read-only accessor. Research use only.

## Supply original files once

Declare samples, patient HLA, VCFs and references with the existing
`neo_config()` → `neo_project()` → `p <- neo_run(p)` workflow.
`neo_run()` returns the initialized/completed **project**, not a results table.
`neo_results(p)` retrieves the existing native results.

Single-sample somatic VCF FORMAT evidence is acquired automatically. No BAM,
`dna_reads` table, preliminary run or Event ID is needed. Optional evidence uses
the existing `somatic_evidence$sources` record array in `neo_config()`:

```r
# Insert into your existing neo_config(...). Replace paths and sample/build
# declarations with verified values; this is configuration, not a full example.
somatic_evidence = list(sources = list(
  list(kind = "rna_bam", sample_id = "S1", path = "S1.rna.bam",
    reference_build = "GRCh38.p14"),
  list(kind = "expression", format = "salmon", sample_id = "S1",
    path = "salmon/S1/quant.sf"),
  list(kind = "ascat", format = "ascat_rds", sample_id = "S1",
    path = "S1.ascat-output.rds", reference_build = "GRCh38.p14",
    options = list(coordinate_system = "1-based-inclusive", normal_cn = 2))
))
```

HLA remains keyed by **patient**, inputs/evidence by **sample**. Each VCF has one
sample column matching its declaration. Multiple separately declared files are
supported; multi-sample VCF interpretation has not changed. Existing generic
DNA/RNA/expression/purity/FACETS/CCF/LOH tables remain supported.

`normal_cn = 2` is an explicit example assumption, **not an inferred measurement**;
use only for scientifically justified loci. Leave it absent when unknown. ASCAT
alone does not establish normal chromosome copy number or mutation multiplicity.

## What DNA VAF means

| Input | Result |
| --- | --- |
| Tumor `FORMAT/AD=80,20` | 20/(80+20) = 0.20 |
| Multiallelic `AD=50,20,30` | ALT1=0.20, ALT2=0.30; denominator includes all alleles |
| `FORMAT/AO=20; RO=80` (FORMAT colon fields) | ALT/(RO+sum(AO)); never ALT/DP alone |
| Declared `FORMAT/AF`, Number=A, describing allele fraction/frequency | Caller-estimated fraction, retained separately from count fraction |
| DP alone, INFO/population AF, ambiguous/undeclared sample AF | No invented VAF |
| AD=0,0 | Observed zero usable depth; VAF unavailable |
| AD=20,0 | Evaluated VAF=0 |

Negative, nonfinite, impossible depths/fractions fail. AD has precedence over
AO/RO. Reported AF disagreements remain in the VCF provenance alongside raw
FORMAT, ALT index and count-derived value. No averaging or silent replacement.
AD totals can exceed DP when their caller filters differ; this is retained as a
depth discrepancy, not used to change the denominator. See [GATK allele depths](https://gatk.broadinstitute.org/hc/en-us/articles/360037421851-DepthPerAlleleBySample).
Normal/germline sample values never become tumor VAF. When normalized Events
merge source records with conflicting count/fraction assertions, the evidence
is marked ambiguous with unavailable VAF, while supported reconstruction continues.
An explicit evidence table can resolve precedence without pooling callers.

Allele-bundle precedence is **complete explicit table/canonical allele evidence →
tumor VCF → DNA BAM fallback**. A depth-only row never hides VCF AD or VAF.
Depth is an independent measurement: an existing VCF/BAM depth stays with its
count bundle; otherwise a same-unit explicit depth may be attached with its own
provenance, without using it to recompute VAF. All explicit observations remain
in `alternatives`, and depth/count or source discrepancies are flagged.
Incomplete explicit allele assertions that contradict a complete source make
assay REF/ALT counts and VAF unavailable, with the original bundle retained for
audit; they are not averaged or made into negative evidence. Conflicting explicit canonical/table values retain their existing
error behavior. Alternatives and disagreements between selected/VCF/BAM
measurements are retained in explanation provenance; they are not pooled.
`neo_explain(p, candidate_id = id)$candidates[[1]]$input_evidence` exposes the
source-file evidence without a second analysis. `neo_priorities(p)` includes
DNA/RNA VAF and biological evidence even when prediction status is `NOT_RUN`.
It also exposes native Grade/rank, REF/ALT support, purity/ploidy and local
copy number directly; users do not need to join these tables by internal IDs.

## Optional BAM allele evidence

`kind="dna_bam"`, `"normal_dna_bam"` or `"rna_bam"` requires existing readable BAM
and exactly one unambiguous BAI/CSI index. Neoantimon2 does not sort, rewrite or
index the user's BAM. Optional Bioconductor packages Rsamtools, GenomicRanges
and IRanges are required only for BAM processing; no predictor is needed.

Defaults: mapping/base quality ≥20, exclude duplicates, secondary/supplementary,
unmapped and QC-failed alignments. `options=list(min_mapq=20, min_baseq=20,
include_duplicates=FALSE, mode="fallback")` is explicit and reproducible.
DNA fallback does not query loci with sufficient VCF/explicit VAF. Request
`mode="independent"` to retain an independent BAM measurement without overriding
higher-priority evidence. RNA BAM counts are obtained when supplied.

Read-group SM must match sample identity. An independently verified upstream
name can be declared with `options=list(sample_name="upstream_tumor_name")`;
this is an explicit identity mapping, never automatic tumor/normal inference.
Missing/unknown read groups, incompatible reference assembly, queried contig
names/lengths, REF mismatches and unreadable indexes fail with a corrective message.
Extra non-target ALT, decoy or unplaced contigs in a BAM header are permitted;
every actually queried contig must exist in both BAM and reference with exactly
the same name and length. This is locus compatibility, not genome-wide BAM identity.
The assembly declaration and reference checks do not prove genome-wide sequence
identity for BAMs lacking trustworthy reference metadata. Previously completed
runs retain their original serialized evidence and method provenance; reading
them does not retroactively repair v1 BAM counts. To apply these corrected
counting rules, run the same declared inputs into a new output project while
preserving the old run. The new source-evidence checkpoint version prevents
v1 results from being reused when this stage is explicitly executed.

Queries use indexed random access. SNVs compare the observed base. Anchored
insertions/deletions assign the **target CIGAR edit**, normalized by the same
reference-based algorithm as VCF Events. This recognizes repeat-equivalent
placements and ignores unrelated aligned substitutions, including nearby SNVs;
it does not require the entire read window to equal a pristine REF/ALT sequence.
Ten-base flanks and the repeat context are still required for coverage/quality.
A target edit with a conflicting anchor or multiple interacting local edits is
ambiguous. A single clearly different indel is `other`, not the target ALT.
`N` is a splice gap, never a deletion allele or REF base. Introns outside the
variant window are allowed. Reads missing necessary context/quality are
unusable. Repeat contexts exceeding 2,000 bases fail explicitly; complex/structural
events have unsupported status and no fabricated counts. This is not validation
of arbitrary complex haplotypes, long-read alleles or structural variants.

Fragments are keyed by read group and query name: agreeing overlapping mates
count once; discordant mates are ambiguous. At an indel locus with unresolved
fragments, primary assay REF/ALT counts and VAF are unavailable, so frozen Grade
rules cannot mistake partial counts for negative evidence. Definitive partial
counts remain in `resolved_ref_fragments`/`resolved_alt_fragments` for QC.
REF/ALT counts, other/ambiguous counts,
quality policy and unit `fragment` are retained. VAF uses REF+selected-ALT usable
fragments; other alleles are excluded from this BAM denominator, unlike VCF AD's
all-allele denominator. Preserve that method distinction when comparing sources.
Zero coverage has zero counts and unavailable VAF; absent BAM has no measurement.
REF means absence of the target indel edit under the supplied alignment; it does
not require absence of every nearby substitution. Mutation-specific read support
does not establish full peptide expression or multi-event haplotype phasing.

Strong file provenance is unchanged: every declared BAM is fully SHA-256 hashed
before checkpoint reuse, even if DNA fallback later skips all allele queries.
For a WGS BAM this is O(file bytes) I/O and can dominate startup; indexed counting
is not the same as avoiding that hash. Bound descriptors are reused by the next
evidence stage, not rehashed there. No stat-only cache or partial hash is passed
off as full integrity. If VCF evidence suffices and no independent measurement is
needed, omit the optional DNA BAM declaration to avoid the full read. No WGS
hash-throughput benchmark was performed for this correction.

## Expression from existing quantification

Choose an explicit `format`: `salmon`, `rsem_genes`, `rsem_transcripts`,
`featurecounts` or the existing `generic` table. Headers are validated against
upstream schemas. Salmon quant.sf and RSEM preserve TPM and estimated counts;
RSEM FPKM is retained separately. featureCounts preserves counts and requires
an exact sample count column (`options=list(column="S1.bam")` if necessary).
No counts/FPKM are silently converted to TPM. Sample assignment is explicit in
the source record because some native single-sample files contain no sample ID.

Transcript versions and gene IDs are retained. Gene-level data joins to the
selected reference-derived transcript models, rather than fuzzy gene symbols.
No BAM expression quantification, transcript assembly or expression inference is
performed. In the additional numerical view, all nonmissing gene/transcript TPM
assertions for an Origin must agree exactly. Different values yield unavailable
TPM with `gene_tpm_status`/`transcript_tpm_status = "conflict"`, never the maximum,
mean or first source. Original values and provenance remain in `expression_evidence`.
Paper scores requiring conflicting gene TPM are unavailable with an explicit
status; prediction-only and VAF-adjusted ordering still work. Existing native
Grade/rank source-selection rules are not replaced by this view. Unmatched identifiers remain unassigned; do not silently rename them.

## ASCAT and conditional CCF

Run ASCAT externally. Supported native schemas are:

- `format="ascat_segments"`: native segments TSV columns
  `sample, chr, startpos, endpos, nMajor, nMinor`.
- `format="ascat_rds"`: `saveRDS(ascat.output, ...)`, containing the native
  `segments` data frame and named `purity` (or older `aberrantcellfraction`)
  and `ploidy` vectors. ASCAT need not be installed to read this data.

The schema is verified against ASCAT's output-writing source at
[61ddf3b / ascat.runAscat.R](https://github.com/VanLoo-lab/ascat/blob/61ddf3b24453eea91134798cc41d4e828a82fa90/ASCAT/R/ascat.runAscat.R).
An explicitly verified upstream identity can be declared as
`options=list(coordinate_system="1-based-inclusive", sample_name="tumor_A")`
for Neoantimon2 sample `S1`. In a multi-sample ASCAT output, only that exact
upstream sample's segments and named purity/ploidy are selected; both identities
and the declared mapping are preserved. Without an explicit mapping, segment
sample identity must match exactly; nothing is guessed from filenames. A
multi-sample object with unnamed purity/ploidy is rejected as ambiguous.

Declare the assembly and `coordinate_system="1-based-inclusive"`; these are
SNP-position segments, not arbitrary BED intervals. Sample mismatch, overlapping
segments, noninteger rounded copy numbers, minor>major and wrong assemblies fail.
Total CN is nMajor+nMinor; raw LogR/BAF are never treated as absolute CN. Assay
origin (array/WGS/etc.) does not alter the import schema. A TSV without purity
metadata does not invent purity/ploidy; use a saved native output when available.

Existing conditional CCF uses
`VAF * (purity * tumor_CN + (1-purity) * normal_CN) / (purity * multiplicity)`.
It preserves supplied CCF, enumerates unresolved multiplicities, and reports
ambiguity/out-of-range instead of choosing a clone. ASCAT does not establish
multiplicity; the existing enumeration is not an allele-specific clonal model
or the original CCFP posterior. Missing normal CN, purity or VAF prevents inference.
For existing CCF/read tables, `contig,pos,ref,alt` plus declared reference_build
can now resolve internal Event IDs automatically; existing event_id tables work.

## Additional numerical priorities

```r
# Completed p; native predictions if available:
v <- neo_priorities(p)
# Completed external analysis a, without touching native results:
v <- neo_priorities(p, analysis = a)
v[, c("candidate_id", "origin_id", "category", "raw_score",
       "prediction_only", "prediction_rank", "evidence_adjusted",
       "dna_vaf", "vaf_missing_substituted")]
```

Rows retain **each candidate–Origin–model channel**, with all Event relationships,
sample, complete HLA, provider/model/run, task, metric, unit, direction and raw
prediction/source provenance. Native Grade/rank remain candidate-level decisions;
`native_priority_origin_id` identifies the original native selection while each
row retains its own Origin-specific numerical inputs. BA, EL and IM remain separate. Duplicate Origins
never inflate the prediction rank; ranks are computed on unique candidates per
sample/class/model/version/run/task/metric/unit/direction. No cross-model consensus.

`model-v1` prediction-only priority applies only direction: score for higher-is-
better, negative score for lower-is-better. Units remain original; this is not a
probability or cross-model scale. A tie-preserving within-channel order quantile
is `(N-rank)/(N-1)` (N=1 →1), based only on the submitted candidate set. It is
**exploratory, candidate-set-dependent ordering**, not model calibration.
The adjusted number is **not comparable across cohorts/candidate sets**, and is
not a biological probability or a BA+EL+IM composite. It is explicitly marked
`adjusted_priority_status="exploratory_candidate_set_dependent"`.
Evidence-adjusted priority is
this quantile multiplied by chosen DNA/RNA VAF. Missing VAF uses zero **only in
this calculation**, with an indicator. Biological missingness stays NA and is
never passed as zero to CCF or Grade eligibility. Compound Origins lacking one
well-defined VAF remain explicitly unavailable. Prediction-only ordering still
works without expression, VAF, ASCAT or CCF.

## Published 2020 equations and discrepancies

For `method="neoantimon-paper-2020-v1"`, the implemented numerical specification is
[Hasegawa et al., Section 2.4](https://academic.oup.com/bioinformatics/article/36/18/4813/5906504):

```
P_I = L_I(I_m) f(A,E) [1 - 2^(-M) L_I(I_w)] C(x)
P_R = L_R(R_m) f(A,E) [1 - 2^(-M) L_R(R_w)] C(x)
f(A,E) = tanh(g(A) E);  g(A) = tanh(5 A)
L_I(z) = 1 / (1 + exp(0.015 (z - 500)))    # IC50, nM
L_R(z) = 1 / (1 + exp(5 (z - 2)))         # binding percentile, percent
C(x) = 1 / (1 - exp(-30 (x - 0.8)))
```

I_m/I_w and R_m/R_w are mutant/patient-normal BA measurements, M is their peptide
mismatch count, E is corresponding gene expression, A selects DNA/RNA VAF, and x
is the supplied/conditional CCF summary used in this numerical evaluation.
Original CCFP is a posterior summary; a conditional point CCF is not that posterior.
The paper does not fix a universally interchangeable expression unit. This
implementation requires explicitly acknowledged gene TPM; counts/FPKM/transcript
TPM are not substituted for gene E. Constants are exactly the published ones,
not recalibrated cutoffs or modern EL/IM calibration.

Paper IC50 inputs require a recognized metric (`BA_nM`, `affinity_nM`, `IC50`,
`ic50`, `ic50_nM`), nM units and lower-is-better direction. Other nM measurements
such as Kd are not silently interpreted as IC50. Binding ranks require actual
binding-percentile semantics, not a raw percentage score.

Supply `wild_type_predictions` as a standard TSV; joins use exact patient-normal
peptide/HLA/provider/model/version/metric/unit/direction. These predictions must
be legitimately obtained separately; the accessor never invents missing WT
binding or executes a predictor. Missing prerequisites have diagnostic status.
`paper_equation` identifies P_I or P_R for each compatible channel.
`paper_without_ccf` explicitly omits C; `paper_with_ccf` needs CCF. At x=0.8 the
score is unavailable with singularity status. Finite negative or near-singular
values (absolute denominator <1e-6 is the diagnostic instability threshold) are
retained for mathematical auditing, labelled unsuitable for automatic
biological prioritization, and never used to replace native ranks.

Pinned original source [2363e88](https://github.com/hase62/Neoantimon/blob/2363e88a22abf031cb3de7386d377f087cf72ff7/R/FunctionsCalculatePriorityScores.R)
uses **A*tanh(5A)** and substitutes **1** for missing VAF. This implementation
uses the paper's tanh(5A), and the owner's **Neoantimon2 zero-substitution
convention**, which is not unchanged historical missing-input behavior.
Original `Mut_EL`/`Wt_EL` are affinity-like inputs to an IC50 equation, not modern
EL scores. The original comparison P uses binding-rank L_R, A*tanh(E), WT/mismatch
penalty and no CCF factor; it is the source's MuPeXI comparison, not P_I/P_R.

## Generate paired mutant/WT BA queries explicitly

`neo_priorities()` does not automatically execute a predictor while inspecting a
completed run. For an already installed/licensed **configured native predictor**,
the shipped version-pinned example generates both query sets, runs the same
provider/model in a separate sidecar, and writes both standard TSVs. No manual
WT table construction is needed:

```r
library(Neoantimon2)
p <- neo_project(path = "/path/to/completed-project")
source(system.file("examples", "paper-binding.R", package = "Neoantimon2"))
# request selects an existing project$config$predictors record explicitly.
plan <- run_paper_binding_queries(p, "/new/separate/paired-BA", request = 1L)
print(plan) # Version/capability inspection; no prediction or output write.
paired <- run_paper_binding_queries(p, "/new/separate/paired-BA",
  request = 1L, metric = "BA_nM", execute = TRUE) # Explicit execution permission.
a <- neo_analyze_external(p, paired$mutant, policy = list(
  mode = "external_rank", provider = paired$plan$provider,
  run_id = paired$run_id, task = paired$plan$task,
  metric = paired$plan$metric, score_direction = "lower_is_better"))
pi <- neo_priorities(p, analysis = a, method = "neoantimon-paper-2020-v1",
  wild_type_predictions = paired$wild_type, expression_unit = "TPM")
```

Use `metric="BA_rank"` and a different output directory for P_R. Both BA inputs
are genuinely obtained from the same explicitly requested provider execution;
no raw EL/IM score is converted into BA. The workflow uses the package's internal
native adapters and is a tested example, **not an additional public API**. It
supports configured NetMHCpan, NetMHCIIpan or MHCflurry only to the extent their
existing adapter/capability supports the requested class, sequence and metric.
The original Native predictions/Grades/ranks are untouched; the imported sidecar
is an independent External analysis. The native adapter deduplicates physical
peptide/HLA queries while per-sample, per-Origin query records and the saved
mapping retain every relationship. Raw sidecar paths resolve within its own
`intermediate/prediction_runs/`; the completion manifest hashes both TSVs, the
query-Origin mapping and the plan. Missing WT/equal-length contexts are listed
as excluded Origins; unavailable or ambiguous provider outputs fail and retain
raw sidecar evidence rather than fabricate scores. WT sequences come from the
Origin's normal context, which may be reference-derived; no matched-normal assay
is implied. Gene TPM and supported CCF are still separate prerequisites for the
complete paper score. No packages, predictors or model weights are installed.

The commands above were exercised with a clearly labelled synthetic NetMHCpan
stub, including both IC50 and binding-rank routes, strict import, explanation and
native/input hash preservation. Real licensed execution and this example's real
MHCflurry/NetMHCIIpan execution remain untested. Existing provider adapter tests
remain separate evidence; synthetic scores are not predictor validation.

## Actually executable offline demonstrations

With this development package installed and optional BAM dependencies already
available, run the shipped script into a new directory:

```r
library(Neoantimon2)
source(system.file("examples", "source-evidence.R", package = "Neoantimon2"))
demo <- run_source_evidence_demo(tempfile("neo-source-demo-"))
```

It creates labelled synthetic inputs and runs A (VCF/HLA), B (indexed RNA BAM),
C (expression/ASCAT), D (independent BA/EL/IM and paper calculations), and E
(missing VAF). It records native equality/reopening and structured provenance.
B is explicitly skipped if Rsamtools is unavailable. Flat CSVs exclude list
columns; the RDS retains full provenance. No licensed prediction engine, model
weights, real patient data, large download or expression quantifier is used.

No synthetic test establishes biological accuracy or clinical efficacy. See
[inputs](inputs.md) for existing technical contracts and [results](getting-started.md)
for the released workflow. This guide describes local unreleased additions.
