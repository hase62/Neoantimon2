# Results, evaluability and explanations

[Getting started](getting-started.md) · [Workflows](workflows.md) ·
[Troubleshooting](troubleshooting.md)

## Trace a result through its evidence

An **event** is a canonical input observation. An **origin** is a particular
transcript/sequence/phase path that produces a peptide. A **candidate** combines
a ligand with a presenter within a sample. Several origins can support one
candidate, so their row counts need not agree. Provider observations, processing,
MS and receptor evidence have their own scopes; a provider row is not a new
biological origin.

After the [installed quick start](getting-started.md):

```r
events <- neo_events(p)
origins <- neo_origins(p)
results <- neo_results(p)
results[, c("sample_id", "ligand_sequence", "mhc_key",
  "prediction_status", "grading_status", "final_grade")]
```

There is one event, seven origins and four candidates. All four have
`prediction_status = "NOT_RUN"`, `grading_status = "NOT_GRADED"` and
`final_grade = "NOT_GRADED"`. No predictor was requested and no quantitative
somatic evidence was supplied. These rows are reconstructed possibilities,
not measurements of recognition or negative binding.

| Field family | Representative fields / purpose |
| --- | --- |
| Identity | `candidate_id`, `candidate_key`, `sample_id`, `ligand_sequence`, `mhc_key`, `hla_class`; exact scoped result and biological ligand/presenter keys |
| Sequence provenance | Origin/event IDs, transcript IDs, coordinates, phase and reconstruction links; inspect `neo_origins()` and explanations |
| Status/evaluability | `prediction_status`, `grading_status`, reason/list fields, self status; optional `technical_eligibility`, `primary_prediction_status`, `rna_status` |
| Grade/rank | `final_grade`; when somatic evidence is configured: `pmhc_band`, `primary_rank`, `rank_within_sample_class`, `rank_status`, selected priority origin |
| Provider provenance | `primary_status_records` and explanation observations, provider/version/run/model/resource identities; values belong to their original metric and scope |

Optional evidence can add columns and attributes. Inspect `names(results)` and
`attributes(results)` before assuming a field is always present. A table cell may
be scalar `NA`; nested array/object fields are list columns and may contain
`NULL`. Preserve these distinctions during export rather than flattening unknown
values to zero or FALSE.

## Select views without changing the analysis

```r
neo_results(p, sample_id = "sample-1", prediction_status = "NOT_RUN")
neo_results(p, grade = c("A", "B"))
neo_results(p, candidate_id = results$candidate_id[[1L]])
neo_events(p, event_id = events$event_id[[1L]])
neo_origins(p, peptide = origins$peptide_sequence[[1L]])
```

Filters are exact nonempty character vectors: alternatives within one filter,
intersection across different filters. An empty filtered view does not mean an
empty original analysis. `reason` searches warning, downgrade, exclusion,
not-evaluated and grade reasons; supply an actual retained code. There is no
`peptide` argument on `neo_results()`; select peptide with origins or explain.
Accessors read and verify persisted artifacts and return fresh views. Missing
unexecuted output raises an error, distinct from a successful zero-row view.

## Missing is not zero or a biological negative

| State | Practical interpretation |
| --- | --- |
| `NOT_RUN` | Candidate projection has no successful evaluated primary prediction. In the demo none was requested. In somatic priority output it can also reflect unavailable, unsupported or failed primary evidence: inspect `primary_prediction_status`, `primary_status_records` and explain |
| `SUCCESS` | A primary prediction was evaluated; it does not establish T-cell recognition or clinical benefit |
| `FAILED` | A retained provider/source result or execution failed. Run-manifest failure and candidate/provider failure are different levels; a completed run may retain failed provider evidence |
| `PARTIAL` | Typed prediction state for incomplete evidence. Read retained observations and reasons for scope; do not infer that every somatic priority projection uses this label |
| `NOT_GRADED` | Required grading evidence/eligibility is unavailable, unsupported or incomplete, or this source lacks a grading policy; it is not grade E |
| `evaluated_negative` | An explicit scoped assay/evaluability assertion with adequate opportunity and provenance; absence of a file does not create it |
| `not_evaluated`, `unavailable`, `unsupported`, `ambiguous` | Distinct evidence-specific limits, not interchangeable negatives; retain the accompanying reason and scope |

Measured zero expression/read support differs from an unmeasured value. A missing
MS ligand does not imply absence unless the assay/scope supports an explicit
negative. Unknown receptor pairing, missing RNA and unresolved self-search scope
remain unknown. A run's completion describes execution, not complete biological
evaluability across every optional evidence family.

## Somatic grades and ranking

Quantitative somatic evidence enables the existing technical eligibility and
priority policy. Defaults require tumor alternative reads >= 5 and normal
alternative reads <= 1 for technical eligibility; unavailable required counts
remain not evaluated. Default RNA support requires at least one alternative read.
Alternative origins are evaluated separately; the chosen priority origin must
support its required events. The policy does not pool incompatible paths to
manufacture eligibility or RNA support.

For the default primary EL percentile-rank metrics:

| Class | Primary provider/metric | Strong band | Binder band |
| --- | --- | --- | --- |
| I | NetMHCpan / `EL_rank` | <= 0.5 | <= 2 |
| II | NetMHCIIpan / `EL_rank` | <= 1 | <= 5 |

These are percentile ranks, not nM affinity or probability. The active policy and
`pmhc_comparisons` retain actual observed values, thresholds and comparison
operators. Explicit supported configuration can select the primary/thresholds;
read the active policy rather than interpreting every provider with these defaults.
The indirect-HLA BA threshold is a separate comparison.

| Grade | Current meaning |
| --- | --- |
| `INVALID` | Technical eligibility is invalid; excluded from rank |
| `NOT_GRADED` | Technical eligibility is not fully evaluable, or the primary prediction/band is not evaluated |
| A | Eligible, strong band, RNA support |
| B | Eligible, strong band, without supporting RNA evidence |
| C | Eligible, binder band, RNA support |
| D | Eligible, binder band, without supporting RNA evidence |
| E | Eligible, evaluated primary prediction outside the binder band |

“Without supporting RNA evidence” includes unavailable support; it is not a
claim of RNA absence. Fusion/splice and long-read candidates retain their
source-specific `NOT_GRADED` policy reasons rather than borrowing somatic grades.

Default ranking is within each sample and HLA class, lexicographically by grade,
primary rank, RNA state, transcript TPM, gene TPM, CCF, tumor alternative reads,
exact-self status, then candidate ID for deterministic ties. Favorable numeric
values sort in their defined direction; missing numeric values sort last.
Supported `ranking_fields` configuration can reorder the defined field set.
This is a scoped priority order, not a universal immunogenicity probability.
Presentation, recognition, HLA LOH, longitudinal, vaccine-set and native allo
overlays do not rewrite the frozen somatic grade/priority/rank.

## Explain a reported candidate

```r
explanation <- neo_explain(p, candidate_id = results$candidate_id[[1L]])
names(explanation)
explanation$candidates[[1L]]$result
explanation$candidates[[1L]]$origins
neo_explain(p, candidate_id = results$candidate_id[[1L]], format = "text")
```

The structured explanation traces the selected candidate to retained origins and
events, reasons, reference provenance, and available provider/optional evidence.
Use it to investigate why a candidate was ungraded or why different origins
support it. Text is a concise view; structured fields carry fuller provenance.
Selectors also accept literal `peptide`, `sample_id`, or retained `reason` codes.
Ordinary explanations require a selector. A configured allo project additionally
allows `neo_explain(project)` to inspect its separate allo graph without a
somatic candidate selector.

Source-extension fusion/splice graphs are available in explanations. Do not assume
the ordinary `neo_events()`/`neo_origins()` views expose every extension record:
the origin accessor adds source/long-read records when long-read inputs are bound.
Aggregate longitudinal projects refer to separate completed sample graphs;
inspect the original sample projects for their individual event/origin views.

Longitudinal summaries use `neo_multisample()`; vaccine sets are in the
`vaccine_set` result attribute and explanation; indirect results are separate
`indirect_hla_result`/`indirect_hla_baseline` attributes. Native allo evidence is
in `allo_native_evidence`. See [workflow examples](workflows.md) for their distinct
evaluability, counting and selection rules.
