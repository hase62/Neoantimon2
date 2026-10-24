# Public API reference

[Documentation index](../README.md) · [Quickstart](getting-started.md)

The eight native functions below remain exported. Four additional external
prediction APIs are also exported and listed separately below; optional
workflows use their supported configuration and evidence inputs.
The v0.3.0 candidate exports all sixteen interfaces below. Its four
[R setup APIs](#r-setup-apis) also exist in the published development source,
but are absent from the immutable v0.2.0 release. The candidate is not yet a
published release; existing signatures and scientific contracts remain intact.

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

These functions are included in the v0.3.0 candidate and the earlier published
development source; they are absent from the immutable v0.2.0 release. Existing
public API signatures and analysis contracts are unchanged. See [R setup](setup.md)
and installed help `?neo_download_resources`.

| Function | Purpose |
| --- | --- |
| `neo_download_resources(profile, resources, dest, ...)` | Explicit public selections; native R minimum-profile download/preparation; advanced selections use the installed existing helper |
| `neo_check_resources(..., config = NULL)` | Report missing/invalid files, integrity states, required tools and optionally full reference compatibility |
| `neo_check_predictors(config = NULL, ...)` | Inspect native tools and runtimes; optional explicit synthetic execution check |
| `neo_setup_predictor(provider, ...)` | Return ordinary predictor configuration; optionally prepare a new pinned MHCflurry environment with explicit consent |
