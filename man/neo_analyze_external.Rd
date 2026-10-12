\name{neo_analyze_external}
\alias{neo_analyze_external}
\alias{neo_external_results}
\alias{neo_external_explain}
\title{Analyze explicitly selected external prediction evidence}
\usage{
neo_analyze_external(project, predictions, policy, converter = NULL,
  converter_id = NULL, registry = NULL, trust = FALSE, timeout = 60)
neo_external_results(analysis, project = NULL)
neo_external_explain(analysis, candidate_id, project = NULL)
}
\arguments{
\item{project}{Completed source project. Required with a character analysis ID.}
\item{predictions}{Standard external-prediction-1.0.0 TSV, or raw input when converting.}
\item{policy}{Explicit named list. Mode is evidence_only, external_rank or somatic_grade.
Ranking and grading require provider, run_id, task, metric and score_direction.
Optional predictor_version/model_version disambiguate channels. Grading also requires hla_class.}
\item{converter}{Local R or Python script path.}
\item{converter_id}{Immutable community converter ID, mutually exclusive with converter.}
\item{registry}{Explicit local or HTTPS registry index. The optional community
index is \url{https://www.hase62.jp/Neoantimon2-Converters/index.json}.
Use this canonical HTTPS URL; the github.io alias can redirect through HTTP.}
\item{trust}{Explicit consent to execute a first-use community converter as arbitrary user code.}
\item{timeout}{Finite converter timeout in seconds.}
\item{analysis}{Returned analysis handle or its fingerprint with project.}
\item{candidate_id}{Canonical candidate IDs to explain.}}
\description{Validates literal sequences and full canonical presenters, fans evidence
out to all matching origins and creates a separately checksummed derived result.
Native results, manifests, IDs and predictor sidecars are unchanged.}
\value{Analyze returns an explicit analysis handle. Results returns a data frame
including native annotations and external evidence, scores and ranks. Explain
returns structured full origins, source graph, evidence and reproducibility metadata.}
\details{Missing predictions remain NOT_PREDICTED and unranked. External ranks
order evaluated scores directionally within sample and HLA class with lexical
candidate ID ties. Percentile grading requires declared MHC percent rank semantics
and original somatic evidence, under external-percentile-grade-v1. It reuses the
existing eligibility, origin, RNA, threshold and rank rules. Extensions without
somatic prerequisites are NOT_APPLICABLE. Derived grades are explicit in
external_final_grade and external_grade_status; external_priority traces the
selected source. Native grade and primary predictor fields are retained.
No arbitrary raw score becomes a grade.
Vaccine, longitudinal, receptor, presentation and allo overlays continue to
consume native results. A completed local project is required; an export bundle
alone is not an executable project. Community feedback is a manual snapshot,
not scientific validation or execution consent. Local converters are arbitrary
user code, not a security sandbox. Native analysis and direct standardized
imports never require the registry. Synthetic interoperability verification does
not certify real licensed predictor accuracy, clinical use or other platforms.}
