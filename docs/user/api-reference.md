# Public API reference

[Documentation index](../README.md) · [Quickstart](getting-started.md)

Current **main (0.3.0.9000) exports 22 APIs**, counted from NAMESPACE `export()`
entries; S3 registrations and imports are not additional public functions.
The immutable published **v0.3.0 has 16 APIs**. The frozen **reproducibility-v1 has
20 APIs**, adding `neo_priorities`, `neo_experimental_yield`, `neo_yield_evidence`
and `neo_report_bundle`. Main additionally includes
[`neo_indirect_autorun()`](indirect-autorun.md) and
[`neo_indirect_class2()`](indirect-class2.md). These do not change either tag.
See [which version to install](../../README.md#choose-a-version); record its
commit as well as version. Existing baseline signatures and scientific contracts
remain intact.

| Function | Main arguments and purpose | Return / side effects |
| --- | --- | --- |
| `neo_config(file = NULL, ..., env = Sys.getenv())` | Load JSON/YAML; overlay named settings and optional environment | `neo_config` S3 list; validates configuration. No analysis or tool execution |
| `neo_project(samples = NULL, hla = NULL, reference = NULL, config = NULL, path = NULL, mode = "somatic", project_id = "project", analysis_context = NULL, inputs = list(), initialize = FALSE, fixture_kind = "real")` | Bind metadata/resources/inputs, or reopen saved run with `path` | `neo_project` S3 list; observes file integrity and derives identities. New run writes only when initialized; allo opt-ins validate bound manifests and derive their separate evidence objects |
| `neo_run(project)` | Execute configured native stages and optional providers/overlays, or verify/reuse checkpoints | Initialized project; writes manifests, checkpoints, generated artifacts and logs under run root |
| `neo_results(project, grade = NULL, prediction_status = NULL, grading_status = NULL, reason = NULL, sample_id = NULL, candidate_id = NULL)` | Candidate data-frame view with exact filters | Fresh data frame; optional overlay columns/attributes; reads verify persisted integrity |
| `neo_events(project, sample_id = NULL, event_id = NULL)` | Inspect canonical input-event rows | Data frame, including empty schema-shaped tables when execution produced no events |
| `neo_origins(project, sample_id = NULL, origin_id = NULL, peptide = NULL)` | Inspect materialized sequence/provenance paths | Data frame; `peptide` selects literal peptide sequence |
| `neo_multisample(project, patient_id = NULL, peptide = NULL)` | Inspect configured subject-scoped multi-sample summaries | Data frame with typed/list cells when longitudinal integration applies; not a patient ranking |
| `neo_explain(project, candidate_id = NULL, reason = NULL, peptide = NULL, sample_id = NULL, format = c("structured", "text"))` | Trace selected candidates to events, origins and available evidence | `neo_explanation` list or character text; requires a selector except query-free configured allo explanations |

| Additional function | Purpose |
| --- | --- |
| `neo_export_bundle(project, dir = NULL)` | Verify/export a completed run's literal peptide/HLA queries and all origins/annotations |
| `neo_analyze_external(project, predictions, policy, ...)` | Strict standardized import, optional local/consented registry conversion, explicit external evidence/rank/eligible somatic grading |
| `neo_external_results(analysis, project = NULL)` | Read the separately checksummed derived result without changing native results |
| `neo_external_explain(analysis, candidate_id, project = NULL)` | Explain external evidence and every preserved candidate origin |

See [external prediction interoperability](external-predictions.md) and installed
help `?neo_export_bundle` / `?neo_analyze_external`. The canonical optional registry
is https://www.hase62.jp/Neoantimon2-Converters/; the default native workflow
requires neither its network availability nor consent to execute a converter.

Selectors accept nonempty character vectors, use exact matches, and combine
different fields as intersections. `reason` matches any retained warning,
downgrade, exclusion, not-evaluated or grade reason. Filtering never changes the
stored candidate universe. Unexecuted/missing output raises an error rather than
returning a successful empty result.

```r
neo_results(p, sample_id = "sample-1", grade = "NOT_GRADED")
neo_events(p, event_id = events$event_id[[1L]])
neo_origins(p, peptide = origins$peptide_sequence[[1L]])
neo_multisample(p)
neo_explain(p, candidate_id = results$candidate_id[[1L]], format = "text")
```

The example's single sample has no configured multi-tumor overlay, so its
multi-sample summary is empty. This does not mean a measured absence.


## R setup APIs

These functions are included in v0.3.0 and are absent from the immutable
v0.2.0 release. Existing
public API signatures and analysis contracts are unchanged. See [R setup](setup.md)
and installed help `?neo_download_resources`.

| Function | Purpose |
| --- | --- |
| `neo_download_resources(profile, resources, dest, ...)` | Explicit public selections; native R minimum-profile download/preparation; advanced selections use the installed existing helper |
| `neo_check_resources(..., config = NULL)` | Report missing/invalid files, integrity states, required tools and optionally full reference compatibility |
| `neo_check_predictors(config = NULL, ...)` | Inspect native tools and runtimes; optional explicit synthetic execution check |
| `neo_setup_predictor(provider, ...)` | Return ordinary predictor configuration; optionally prepare a new pinned MHCflurry environment with explicit consent |

## Development evidence, planning and reporting APIs

| Function | Actual signature and behavior |
| --- | --- |
| `neo_priorities(project, analysis = NULL, method = c("model-v1", "neoantimon-paper-2020-v1"), vaf = c("dna", "rna"), wild_type_predictions = NULL, expression_unit = NULL)` | Read existing Native/External numeric priorities with model and Origin provenance; paper methods have separate prerequisites |
| `neo_yield_evidence(x, reference = NULL, endpoint = NULL, assay = NULL, min_records = 100L, min_patients = 3L)` | Descriptive A–D/B1/B2 evidence; no probability adjustment; explicit unavailable status |
| `neo_experimental_yield(x, profile, reference = NULL, analysis = NULL, channel = NULL, unit = c("peptide", "peptide_hla"), ordering = c("score", "input"), primary_hla = c("best_score", "lexicographic"), min_records = 100L, min_patients = 3L)` | Existing predictor evidence or compatible score table; independent scenario/empirical profile; units, curve, targets and provenance |
| `neo_report_bundle(project, dir = NULL, analyses = list(), candidate_id = NULL, yield_plan = NULL, external_channels = NULL, reference = NULL, endpoint = NULL, assay = NULL)` | Read-only four-level report; separate channel representatives and all links; optional unchanged planning tables |

See installed help for defaults and validation, [source evidence/priorities](source-evidence-and-priorities.md),
[Experimental Yield](experimental-yield.md) and [complete reporting](result-reporting.md).
No function executes a predictor merely to read a saved result. The reporting
export has its own contract; `neo_export_bundle()` remains unchanged.

## Optional standalone indirect-HLA autorun

`neo_indirect_autorun(allo, presenters, reference, executable = NULL, model = NULL, output_dir, dry_run = FALSE, resume = FALSE, authorization = NULL, timeout = 600L, batch_size = 256L, predictor_version = c("4.3k", "3.2", "auto"), inversion = c("vendor_default", "all"))`

[Setup and exact examples](indirect-autorun.md). Standard 4.3k BA uses its own orientation-aware assumption profile; explicit 3.2 preserves the original strict baseline; completed-run resume is hash checked, and mock runs are explicitly synthetic. This API is not in immutable `v0.3.0` or frozen `reproducibility-v1`.

## Expanded indirect class II workflow

`neo_indirect_class2(donor, recipient, profile, installation, reference, output_dir, presenters=NULL, source_loci, self_loci, thresholds, inversion, dry_run=FALSE, resume=FALSE, batch_size=256L, timeout=600L, min_free_gib=25)` is the optional high-level research wrapper. [Full signatures, reference preparation, uncertainty and examples](indirect-class2.md). It does not change the legacy autorun contract.
