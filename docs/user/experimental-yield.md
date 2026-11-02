# Experimental Yield Estimation v0.1

Experimental Yield Estimation provides an assumption-based estimate of the
number of peptide candidates required to reach a specified probability of
detecting at least one positive result. Evidence categories A–D summarize prior
experimental information and **do not alter the probability calculation**.

This is an optional, read-only development feature, separate from the immutable
v0.3.0 release. It never runs predictors, downloads a database, changes Native
A–E Grades/ranking, or modifies BA/EL/IM observations or Event/Origin links.

## Fast executable example

Run the bundled, entirely **synthetic** offline example:

```r
library(Neoantimon2)
source(system.file("examples", "experimental-yield.R", package="Neoantimon2"))
plan <- experimental_yield_example()
plan$summary
plan$targets
plan$curve
plot(plan)
```

The example uses explicitly assumed rates and a fixed hypothetical score scale.
It creates no project, executes no predictor and makes no biological claim.
Targets beyond the available pool return `minimum_N = NA`, with a reason. Never
create extra candidates simply to make a target reachable.

## Existing native or imported external output

```r
# p is an existing completed neo_project, reopened with neo_project(path=...).
# These accessors read existing results; do not call neo_run().
x <- neo_priorities(p)                       # native score/provenance view
# Or: x <- neo_priorities(p, analysis=existing_external_analysis_id)

# Select exactly ONE sample and one predictor/model/run/task/metric.
# Replace these literals with identities actually present in x.
x <- x[x$sample_id == "YOUR_SAMPLE_ID", , drop=FALSE]
plan <- neo_experimental_yield(
  x, profile="improve-exploratory-v0.1", reference=NULL,
  channel=list(provider="YOUR_PROVIDER", metric="YOUR_METRIC"),
  unit="peptide"
)
```

The second block is an API template requiring your actual saved project and
identities, not a self-contained runnable dataset. A completed project can also
be supplied directly as the first argument, with `analysis=` when using an
existing external analysis. `channel=` accepts exact identity-field filters;
ambiguity between runs/models/metrics fails. BA, EL and IM are not pooled.
Unpredicted candidate rows with no predictor identity are retained as unavailable,
not silently dropped. Filter a **deliberately selected scored pilot subset** first
if that is the intended experimental pool; record that selection explicitly.

All original row fields remain in `input_provenance` and the per-pair
`source_rows` list column, including original scores, units, direction, candidate
IDs, sample, all Origins, Events, model/run provenance and existing Grade/rank.
Conflicting duplicate observations fail instead of being averaged or maximized.
Use bounded small saved subsets; a full-project priority view can be expensive.

## Independent A–D evidence lookup, including both B subtypes

`neo_yield_evidence()` works without a score or probability profile:

```r
pairs <- data.frame(peptide_sequence="ACDEFGHIK", mhc_key="HLA-A*02:01")
neo_yield_evidence(pairs) # explicit REFERENCE_UNAVAILABLE; category NA
# With a local audited reference:
# ev <- neo_yield_evidence(pairs, "/absolute/path/reference.json",
#   endpoint="T_cell_recognition", assay="DNA_barcoded_pMHC_multimer")
```

Priority is **A → B → C → D**:

| Category | Meaning |
| --- | --- |
| A | Exact normalized peptide and exact full HLA have an experimental result. |
| B1 | Exact peptide, different exact HLA key. |
| B2 | Exact HLA, equal-length peptide differing at exactly one residue (Hamming1). |
| C | No A/B; sufficient relevant tested units and independently identified patients for that HLA. |
| D | No A/B and insufficient relevant reference data; not “low immunogenicity.” |

The new label is stored in `evidence_category`; the existing BA/EL/IM `category`
column is preserved.

B1/B2 have **no validated ordering**. Both are retained when applicable, even
when A takes classification priority. Each A/B1/B2 group independently reports
`positive`, `negative`, `conflicting` or `none`; neither B subtype requires a
positive result. Every matched record retains its actual allele, sequence,
substitution (for B2), patient/study, assay and citation. Cross-patient conflicting
outcomes are not a claim that one patient had contradictory measurements.
Previously positive does not imply positive in a new patient; B2 does not
establish equivalent TCR recognition.

Select one explicit endpoint and assay. Incompatible assays are not silently
combined. A missing database or absent requested assay/endpoint scope produces
**NA category and an unavailable status**, not D or an experimentally negative
label. An available relevant database with little evidence for an HLA yields D.

Default C thresholds are **100 tested units and3 independent patients**. These
are configurable provisional engineering choices, not validated cutoffs. Units
are deduplicated by study/patient/peptide/exact HLA so technical replicates and
conflicting repeats do not inflate counts. Unknown patients do not count. For
study-local patient IDs, the maximum distinct-patient count within a single
study is used; cross-study identity independence is not guessed. A reference
may explicitly declare `global_verified` identities instead.

Peptides are uppercased and trimmed; I/L are not equated. HLA prefix and star
punctuation are normalized (`HLA-A02:01` → `HLA-A*02:01`), but allele fields and
expression suffixes are never truncated or fuzzily substituted. Class II pairs
require explicit compatible alpha/beta chains. No supertypes are inferred.

## Replaceable local reference schema and audited importer

A JSON/list reference has:

- `schema_version: "experimental-reference-0.1.0"`, `id`, `version`, `provenance`;
- `patient_identity: "study_local"` or `"global_verified"`;
- `records`, each with `record_id`, `peptide_sequence`, `mhc_key`, `endpoint`,
  `assay`, `result` (`positive` or `negative`), `patient_id` (can be missing),
  `study_id`, `citation`. Additional source fields are retained.

Never turn untested or uncertain records into negatives. Duplicate record IDs
fail. Files are read explicitly, hashed with SHA-256 and not modified.

An independently authored importer is bundled at
`inst/examples/import-improve-reference.R`. It accepts only the audited original
`01_Validated_neoepitopes.txt` member of the official IMPROVE paper archive at
commit `c67094276b52358b5202ad9e5597acd87888ea1d`, SHA-256
`b27854c98a34b6bb83af248979cc3ea121534fa2e8d9e97415654dc0c0a95c95`.

```r
source(system.file("examples", "import-improve-reference.R", package="Neoantimon2"))
# Supply the genuine local source; output must not already exist.
# import_improve_reference("/local/01_Validated_neoepitopes.txt",
#                         "/local/new-reference.json")
```

The source contains17,520 explicitly sufficiently tested records,467 positive
and17,053 negative, across70 source patients and36 HLA keys. The importer
preserves `yes/no/1/0`, original patient identity, source HLA/record and validation
status. Assay/endpoint annotations come from the publication's T-cell screening
methods: DNA-barcoded pMHC multimers, with patient-derived samples. It is **not**
HCC1395 experimental truth and does not calibrate BigMHC. No source records are
bundled: upstream repository redistribution terms were not established. Users
must obtain data under its source terms. A future IEDB/CEDAR adapter must retain
record-level endpoint/assay and citations; do not treat a simplified benchmark
`response` table without patient/assay metadata as a fully characterized reference.

## Four independent score bins and two probability modes

Bins are ordered by **prediction quality**, respecting score direction:
1 bottom30%,2 30–60%,3 60–90%,4 top10%.

### Scenario planning

An explicit user profile supplies four probabilities. Minimal schema:

```r
profile <- list(
  schema_version="experimental-yield-profile-0.1.0",
  id="my-assumptions", version="1", mode="scenario",
  endpoint="example-T-cell-endpoint", assay="example-assay",
  population="My exploratory experiment design",
  limitations="Assumed rates, not observed calibration",
  provenance=list(source="Investigator-selected scenario"),
  binning=list(basis="fixed_thresholds", thresholds=c(.3,.6,.9),
               score_direction="higher_is_better"),
  rates=c(.01,.02,.03,.06)
)
```

These illustrative cutpoints require a compatible score scale; they are not
universal. Ordinarily use frozen reference30th/60th/90th raw-score quantiles.
Increasing quality means increasing raw cutpoints for higher-is-better metrics,
and decreasing raw cutpoints for lower-is-better metrics. Equality stays in the
lower bin; ties can create empty bins.

For exploratory cohort-relative bins set `binning=list(basis="cohort_relative")`.
Mid-rank percentiles `(rank - 0.5)/n` are computed over **distinct scored
sample/peptide/exact-HLA pairs**, never repeated Origin rows. These bins are
candidate-set dependent and cannot be described as fixed empirical calibration.

The built-in `improve-exploratory-v0.1` rates are0.0139,0.0230,0.0306,0.0639.
They were supplied by the owner from a reanalysis of the published IMPROVE
cross-validation score distribution; exact CV-bin raw counts/thresholds were
not independently reconstructed here. They are shipped only as an **assumed
reference scenario**. Applying them to BigMHC, BA or EL ordering explicitly
reports `EXPLORATORY_CROSS_MODEL_ASSUMPTION_NOT_CALIBRATED`. Another model's raw
score thresholds cannot be transplanted; cross-model scenarios require the
exploratory cohort-relative approach.

### Empirical profile

A versioned empirical JSON/list uses the same common metadata and adds:

- `mode="empirical"`, fixed thresholds;
- `model` containing exact `provider`, `predictor_version`, `model_version`,
  `task`, `metric`, `metric_kind`, `metric_semantics`, `unit`, `score_direction`;
- `bins`: four ordered records with `bin`, integer `tested`, `positive`,
  `negative`, and `rate` equal to `positive/tested`;
- explicit dataset provenance, endpoint/assay, calibration population and limitations.

Missing model versions or incompatible providers/metrics/units/directions fail.
There must be tested observations in each bin; positives+negatives must equal
tested. This is a strict profile reader, not a training framework. A matching
user-supplied table is labelled `empirical_model_endpoint_matched`, with an
explicit transportability warning. The software cannot independently establish
that a user-supplied dataset is unbiased, held-out or clinically validated.
**No independently validated predictor-specific empirical profile is shipped.**
Frozen reference scores, preprocessing settings and population should be
archived alongside any future empirical profile, not only the rate table.

## Experimental units, ordering and probability calculation

Default `unit="peptide"` counts unique synthesized sequences. A deterministic
primary HLA is chosen by best direction-aware score (ties: lexical HLA then
candidate ID), or explicitly by `primary_hla="lexicographic"`. This is a labelled
single-HLA **proxy**, not a joint probability across a person's HLA genotype.
The best-score rule can be optimistic. All alternative associations, candidate
IDs and complete original Origin/source rows remain in `alternative_pairs`.

`unit="peptide_hla"` counts exact assays and separately displays distinct
synthesized peptides. Repeated Origins are not extra assays or syntheses. When
a peptide appears with multiple HLAs, the expectation counts positive assays;
joint success probability/targets become unavailable after that peptide repeats,
rather than silently treating its HLA outcomes as independent. No arbitrary
HLA probability sum/max is used as a peptide success probability.

The additional view orders by predictor quality, or `ordering="input"`; it
never overwrites existing Native rank. For probabilities `p_i`:

- expected positive units: `sum(p_i)`;
- probability of at least one: `-expm1(sum(log1p(-p_i)))`, assuming independence;
- minimum N for50/80/90/95% is found only within the actual available ordered pool.

The curve includes N=0 and cumulative distinct synthesis counts. Probabilities
0 and1 are valid. Missing predictions stay NA; biological VAF remains untouched.
Full-pool expected/probability values are NA when unavailable, and separately
named `*_scored_subset` summaries describe only known rates. `reachable=FALSE`
means the known pool cannot reach the target; `reachable=NA` means missing/joint
information prevents deciding. The reason is returned for every target.

Different candidates can be correlated through patient/tumor/T-cell biology;
independence is a planning assumption, not a finding. pMHC-multimer rates are
not automatically ELISpot, stimulation, tumor-recognition or clinical response
rates. Neither A/B/C/D, B1/B2, Native Grade, VAF nor expression multiplies the
four-bin probability.

## Evidence used for v0.1 verification

Synthetic fixtures independently exercise exact and conflicting A evidence,
both B subtypes including negatives, C/D thresholds, absent references, bins,
model mismatch, deduplication, missingness, stable analytical probabilities,
unreachable targets and read-only Native/External access. A separate real-output
demo reads the previously saved75 HCC1395 BigMHC-IM peptide/HLA predictions and
uses15 synthesized-peptide proxies. This demonstrates software semantics and
metric/provenance retention; it is **not immunogenicity validation**. The original
whole-sample HCC1395 analysis, predictor integrations and pVACseq benchmarks are
not rerun.

Sources: [Borch et al., IMPROVE,2024](https://doi.org/10.3389/fimmu.2024.1360281),
[official fixed data repository](https://github.com/SRHgroup/IMPROVE_paper/tree/c67094276b52358b5202ad9e5597acd87888ea1d),
[Wells et al., TESLA,2020](https://doi.org/10.1016/j.cell.2020.09.015),
[IEDB database export](https://www.iedb.org/database_export_v3.php).
TESLA/IEDB endpoints and source sampling differ; no automatic pooled calibration
or general-purpose patient-independent probability claim is made.
