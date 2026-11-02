# Complete result reporting

`neo_report_bundle()` is an optional, read-only projection of saved Event, peptide, pMHC and Origin records. It does not run reconstruction, predictors or imports, and does not change Native Grade, rank, External evidence or Experimental Yield probabilities. The established `neo_export_bundle()` query/interoperability contract remains unchanged.

## Saved project and explicitly selected analyses

```r
library(Neoantimon2)
p <- neo_project(path = saved_project_path)
r <- neo_report_bundle(p, analyses = list(existing_BA_analysis,
                                         existing_EL_analysis,
                                         existing_IM_analysis))
r$tables$event_summary
r$tables$event_best_by_channel
r$tables$pmhc_predictions
# A fresh directory; never inside the project, inputs, reference or repository:
r <- neo_report_bundle(p, dir = new_report_directory,
                       analyses = list(existing_IM_analysis))
```

The variables above refer to existing paths/analysis IDs, not a new convenience API. The compact executable synthetic example is `inst/examples/result-reporting.R`. The review package also contains an executed HCC1395 saved-output demonstration and its receipts. No full-sample analysis is needed for reporting.

For a bounded report, pass `candidate_id = exact_candidate_ids`. Every original Event node stays visible, but association counts and Origins/links are explicitly scoped to these candidates. Without that filter, all Origins remain, including those having no HLA candidate. Such Origins are not evaluated pMHCs. Ordinary multisample projects are supported with separate sample scopes. For longitudinal aggregate handles, report each underlying sample project individually; aggregate reporting is not supported in v0.1.

## Event representatives and complete alternatives

`event_summary.tsv` contains counts, all candidate/Origin references and status for every Event. `event_best_by_channel.tsv` provides separate representatives:

- Native uses the existing persisted rank only when its primary prediction is evaluated and rank status is RANKED, separately by sample, MHC class, policy and analysis context. A display rank on NOT_GRADED/unpredicted results is not evidence of a best evaluated neoantigen.
- External uses raw score direction within each exact analysis/provider/version/model/run/task/metric/kind/semantics/unit/prediction-scope/sample/class. Scores from BA, EL and IM are never combined. Identical scores break ties by lexical candidate ID.

`external_channels = list(list(provider = "bigmhc-im"))` restricts representatives, not the long prediction table. Every observation in every explicitly selected analysis remains. A different run or metric gets a different representative. Peptide-only predictions retain their model-input scope; HLA fan-out is not represented as HLA-conditioned inference.

If no evaluable candidate exists, the Event and scope remain with `UNAVAILABLE` and a reason. Representatives do not replace alternatives. The chosen Origin must link to the Event and candidate: use the Native priority Origin when associated, otherwise the lexical associated Origin. The reported Native Grade remains the global candidate Grade, not a new Event-specific Origin Grade.

`origin_event_links.tsv` distinguishes joint generation source-set members from germline/non-somatic contextual variants. A compound Origin is not evidence that each Event independently generates that peptide. Source extension and long-read IDs have explicit namespaces; exact annotation links are used, never gene-name or coordinate guesses. Reporting-only `event_ref` keys do not replace original Event IDs. Missing source metadata remains explicitly unavailable.

## Unique peptides and Experimental Yield

```r
# Same established optional planning API and unchanged score/HLA proxy:
plan <- neo_experimental_yield(saved_scores, "improve-exploratory-v0.1",
                              reference = explicit_reference_path)
r <- neo_report_bundle(p, analyses = list(existing_IM_analysis),
                       candidate_id = unique(saved_scores$candidate_id),
                       yield_plan = plan)
```

`peptide_summary.tsv` has one sequence per sample and includes all alternative HLA, candidate and Origin associations. Its planning fields copy the supplied unique-peptide plan's primary-HLA proxy, score bin and probability. Alternative HLAs do not become additional synthesized peptides or independent extra probabilities. An assay-unit plan remains separately identified; it does not fabricate a single-HLA peptide proxy.

All original planning units,  pair evidence, curve, targets and summary remain in `experimental_yield_*.tsv`; the curve/order and target mathematics are not recomputed by reporting. Prediction links identify which selected persisted channel supplied a planning score, or explicitly flag an unbound channel. A supplied reference must agree with the plan's identity, content and assay/endpoint. Existing configurable evidence thresholds are reused.

A--D and B1/B2 evidence is descriptive. It never changes probability. The IMPROVE-derived scenario applied to BigMHC-IM is an exploratory cross-model assumption, not BigMHC calibration or experimental validation. Shared peptide biology can violate the independence assumption.

## Complete long predictions and molecular provenance

- `events.tsv`: full original Events, original identity/location and source records.
- `pmhc_candidates.tsv`: every selected pMHC and its active public Native Grade/rank plus original saved candidate record.
- `pmhc_predictions.tsv`: every Native persisted prediction/status and selected External evidence row. Keep provider/model/run/analysis/task/metric, unit, direction, status, original score text and public numeric value. Explicit `coverage_status` rows mean no matching observation; they are not prediction measurements or experimental negatives.
- `origins.tsv`: every selected full molecular Origin, transcript/gene/protein, mutation context, WT comparison and source/provenance fields.
- `candidate_origin_links.tsv`, `origin_event_links.tsv`, `candidate_peptide_links.tsv`: explicit many-to-many joins. No display-string parsing is necessary.
- `molecular_event_evidence.tsv`, `expression_evidence.tsv`: available original evidence. Join molecular measurements by declared Event ID; a compound Origin does not collapse several VAFs into one invented scalar.
- `source_observations.tsv`, `long_read_observations.tsv`, `native_run_statuses.tsv`: additional saved source/run records where available.

Source-derived VAF retains its actual measurement type and provenance. Missing evidence stays missing, observed zero stays zero. Gene symbols are annotations, not unique identities. Complete original records and active Native projections are both retained because an original graph snapshot and a subsequently materialized policy may legitimately differ.

## Serialization and reproducibility

The independent reporting contract is `result-report-0.1.0`; it has its own manifest, identities and SHA-256 table hashes. It does not change export-bundle hashes. Output is staged and renamed to a new directory; existing directories cannot be overwritten. Saved scientific files are read only.

TSVs are quoted and list-valued cells are lossless JSON; `source_record_json` preserves the complete original record, including explicit null and empty arrays. The manifest describes column encodings. New numeric serialization preserves available R doubles with 17 significant digits. Original predictor decimal strings are retained separately because the existing public projection can already have a documented 15-significant-digit boundary; reporting cannot recover precision absent from that projection. Missing scalar cells are empty; use original JSON to distinguish original null/empty string where relevant. These exports do not retain arbitrary R attributes or serialize executable R objects.

For example:

```r
m <- jsonlite::read_json(file.path(new_report_directory, "manifest.json"))
candidates <- read.delim(file.path(new_report_directory, "pmhc_candidates.tsv"),
                        check.names = FALSE, colClasses = "character")
origin_ids <- jsonlite::fromJSON(candidates$origin_ids[[1]])
```

Only explicitly requested analyses/references are included. An absent reference gives evidence-unavailable rather than a fabricated D or negative result. Reporting does not download reference databases. Keep confidential or restricted experimental records out of shared bundles unless redistribution is permitted.
