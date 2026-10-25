# Neoantimon2 analysis playbook

[Shared rules](AI_INSTRUCTIONS.md) · [Manifest](INPUT_MANIFEST.example.yaml) ·
[Requests](EXAMPLE_REQUESTS.md) · [Inputs](../user/inputs.md) · [Setup](../user/setup.md)

This is the authoritative workflow for both Codex and Claude Code. Use R 4.1+
and the already installed v0.3.0 package/dependencies. The instructions work with
an installed package; repository utilities are explicitly labeled. No specific
home directory, development library or hidden environment variable is required.

## 1. Inspect and choose a mode

Locate the manifest, configuration, reference provenance, VCFs and any completed
project. Read the original declarations, not normalized internal project records.
Check `packageVersion("Neoantimon2")` and `getNamespaceExports("Neoantimon2")`.
Record the source commit if available, R/platform and `sessionInfo()`.
The 16 exports are listed in [API reference](../user/api-reference.md).

The YAML example is an **agent declaration**. The agent translates it into the
real API: samples data frame, patient-keyed HLA list, unnamed list of input
records and an explicit `neo_config()`. `neo_project()` does not accept this
manifest file. Do not invent `neo_run(samples=..., hla=..., inputs=...)`.

| Mode | Procedure |
| --- | --- |
| A: offline demonstration | Run the installed synthetic Quick Start below; no human reference download or predictor |
| B: one patient / one sample | One sample row, one patient-keyed HLA genotype, one input record |
| C: one patient / multiple samples | Repeat patient ID across samples; declare HLA once; separate single-sample VCF records |
| D: multiple patients | Distinct patient IDs and independent HLA list entries; preserve every sample occurrence |
| E: external prediction | Complete a project, export, separately run an authorized predictor or use a labeled fixture, convert and import as described below |
| F: explanation only | Open the existing completed project and inspect saved results; do not rerun without reason |

Longitudinal, vaccine-set and restricted indirect-HLA/allo workflows require
explicit settings and their own prerequisites in [workflows](../user/workflows.md).
Multiple samples, filenames and informal timepoint labels do not activate them.
Do not infer chronology or turn absent evidence into evaluated absence.

## 2. Offline demonstration (Mode A)

```r
library(Neoantimon2)
stopifnot(as.character(packageVersion("Neoantimon2")) == "0.3.0")
e <- new.env()
source(system.file("examples", "quickstart.R", package = "Neoantimon2"), local = e)
p <- e$p
stopifnot(nrow(neo_events(p)) > 0L, nrow(neo_origins(p)) > 0L,
  nrow(neo_results(p)) > 0L)
r <- neo_results(p)
x <- neo_explain(p, candidate_id = r$candidate_id[[1L]])
reopened <- neo_project(path = p$path)
stopifnot(identical(neo_results(reopened), r))
```

Tiny references and peptide lengths are technical fixtures, not a production
profile or biological validation. The example prints a temporary run path;
copy required outputs to a new authorized durable location before the R session
ends. The example must report NOT_RUN and NOT_GRADED when no predictions/evidence
were requested. Completion does not imply predictor success or immunogenicity.

## 3. Preflight real inputs (Modes B–D)

Verify each sample has exactly one patient and a unique sample ID. HLA must be
keyed by patient and explicitly supplied; ask about absent/ambiguous declarations.
The package can retain events/origins without HLA and yield no candidates; this
agent preflight intentionally stops rather than interpreting that as a negative.
Require explicit supported Class II alpha~beta pairs. Do not type or guess HLA.

Inspect VCF header sample columns and exact sample bindings; require one sample
column per supported VCF. Reject a multi-sample file or request an explicitly
provided, reproducibly extracted single-sample file. Do not rewrite the header
or select a column silently. Check declared assembly, contig names and REF against
the selected reference; parser errors must not be repaired by guessing an alias,
assembly lift-over or alternative allele. See [input contract](../user/inputs.md)
and [four executable cohort examples](../user/multi-sample-examples.md).

Use the existing reference configuration; assembly/release labels must match
actual bytes and indexes. A convenient initial check is:

```r
cfg <- neo_config(file = "/path/to/reviewed-config.yaml", env = character())
resource_report <- neo_check_resources(config = cfg, checksum = TRUE,
  validate_reference = FALSE)
stopifnot(resource_report$ok)
predictor_report <- neo_check_predictors(config = cfg, smoke_test = FALSE)
```

Replace placeholder paths first. Full reference validation/production runs can be
expensive; agree on cost before `validate_reference=TRUE` or execution. A check is
not proof that every reference is assembly-compatible or every predictor allele
works. Missing optional predictors are not a failed biological prediction.
When no predictors are configured, the package may inspect installed provider
availability; it does not enable them automatically.

The released setup signatures are:

```r
# Signatures, not commands to download or install:
# neo_download_resources(profile=NULL, resources=NULL, dest, overwrite=FALSE,
#   python=NULL, timeout=600)
# neo_check_resources(profile=NULL, resources=NULL, dest=NULL, config=NULL,
#   checksum=TRUE, validate_reference=FALSE, python=NULL)
# neo_check_predictors(config=NULL, providers=NULL, python=NULL, smoke_test=FALSE,
#   run_root=tempdir(), timeout=30)
# neo_setup_predictor(provider, executable=NULL, model_path=NULL, task=NULL,
#   options=list(), python=NULL, environment=NULL, install=FALSE, version=NULL,
#   dependencies=character(), smoke_test=FALSE)
```

Inspect existing resources first. Downloads need explicit destination/selection
and authorization for size; use `overwrite=FALSE`. Configure legitimately obtained
licensed tools, without redistribution. MHCflurry installation needs explicit
consent, a new environment and pinned versions; never mutate the active/global
Python environment silently. A capability/version probe is not prediction
execution. No setup API is required to continue using an existing valid config.

## 4. Declare, plan and execute

The complete executable declarations for all four cohorts are maintained in
[the cohort guide](../user/multi-sample-examples.md) and installed
`system.file("examples", "multi-sample.R", package="Neoantimon2")`.
Their API pattern is:

```r
samples <- data.frame(sample_id = c("S1", "S2"), patient_id = c("P1", "P1"))
hla <- list(P1 = c("HLA-A*02:01", "HLA-B*07:02"))
inputs <- list(
  list(input_id = "S1-somatic", input_role = "somatic_vcf", sample_id = "S1", path = "/path/to/S1.vcf"),
  list(input_id = "S2-somatic", input_role = "somatic_vcf", sample_id = "S2", path = "/path/to/S2.vcf"))
# Use validated cfg from preflight. Replace paths before execution.
p <- neo_project(samples = samples, hla = hla, inputs = inputs, config = cfg)
# Show the plan and obtain any missing expensive-run authorization before:
p <- neo_run(p)
```

Use a fresh run root outside original inputs/reference directories. Keep an
existing completed run immutable; reopen it instead of overwriting it.
Summarize patient/sample counts, resource identities, intended providers,
estimated cost, output destination and known prerequisites in the plan.

The optional repository-only helper [run-analysis.R](run-analysis.R) translates
this manifest for the conservative somatic single-sample-VCF route. It accepts
canonical two-field HLA/pairs, checks bindings/resources, rejects nonempty
optional evidence/external declarations for explicit review, and does not probe
or install predictors. It is not a general replacement for the full package
input contract. It performs a read-only preflight by default:

```sh
Rscript --vanilla docs/ai/run-analysis.R /absolute/path/manifest.yaml
# Only after the user has authorized the configured execution/cost:
Rscript --vanilla docs/ai/run-analysis.R /absolute/path/manifest.yaml --execute
```

The session directory must be new. The helper preserves input hashes, saves
configuration/results/explanation locally, verifies reopening, and writes a
counts/status/provenance report. Provider configuration alone does not establish
availability: inspect trusted tools with the setup API before a real prediction
run. It does not skip a missing requested provider or silently invent scores.
Optional RNA/CCF/MS inputs require their actual [evidence schemas](../user/inputs.md),
including source kind/scope/configuration. Do not append arbitrary TSV paths as
somatic_vcf records or claim they were analyzed when they were only declared.

## 5. Inspect, explain and report (also Mode F)

```r
# For explanation-only work, open the completed run first:
p <- neo_project(path = "/path/to/completed/run")
events <- neo_events(p, sample_id = "S1")
origins <- neo_origins(p, sample_id = "S1")
results <- neo_results(p, sample_id = "S1")
if (nrow(results)) explanation <- neo_explain(p, candidate_id = results$candidate_id[[1L]])
```

Inspect technical/evidence statuses, not just output counts. An empty candidate
set can reflect technical exclusions or missing HLA, not a negative tumor.
NOT_RUN is not a numerical zero/negative; NOT_GRADED is not grade E. Report why
those states occurred and which dependencies or observations were absent.
Trace selected candidates to all origins/events and preserve original IDs. Select
patients through their sample IDs; `neo_results()` has no patient_id argument.

Retain the exact R script/manifest/config, `sessionInfo()`, package/source and
reference/model versions when known, input/resource checksums, logs and local
outputs. Distinguish true observed negatives from missing evidence. Redact
identifiable details in any shared report; do not copy patient data into Git.
Use a short report: execution facts; sample-specific candidate counts; evidence
coverage; traced examples; limitations; reproducible local paths. Do not present
candidate counts as clinical efficacy or model accuracy.

## 6. External prediction (Mode E)

Retain the completed original project. Native scores/grades are not required:

```r
bundle <- neo_export_bundle(p)
# Separately execute an authorized/installed predictor, or use clearly labeled
# synthetic fixtures. Review the pinned converter and its license/checksum.
# "standard.tsv" below is an actual converted standard-format prediction file.
a <- neo_analyze_external(p, "/path/to/standard.tsv",
  policy = list(mode = "evidence_only"))
er <- neo_external_results(a)
if (nrow(er)) ex <- neo_external_explain(a, er$candidate_id[[1L]])
```

Actual converter selection adds `converter_id`, `registry` and explicit first-use
`trust=TRUE` after code/checksum review and authorization. Scripts run with user
permissions, not in a sandbox. Do not execute unreviewed code or upload queries.
Registry scripts adapt files; they do not execute the model or provide weights.
Use the [offline BA/EL/IM walkthrough](https://github.com/hase62/Neoantimon2-Converters/blob/main/EXTERNAL-QUICKSTART.md)
and [core external example](../../scripts/development/external_prediction_quickstart.R).
Strict peptide or peptide+presenter matching preserves every original Origin.

Select an explicit provider/run/task/metric/direction for external ranking.
BA nM, binding percentile, EL score/rank and IM score/rank remain different
channels. IM raw scores/percentiles are not EL percentiles and do not qualify
for somatic A–E grades. The released eligible-percentile/core-prerequisite route
is conditional; never request grading merely because a number exists. External
results/explanations are separate derived analyses; native evidence, grades and
ranks remain unchanged. No automatic consensus across incompatible scales.
