\name{neo_results}
\alias{neo_results}
\title{Inspect candidate results and optional evidence overlays}
\usage{
neo_results(project, grade = NULL, prediction_status = NULL,
  grading_status = NULL, reason = NULL, sample_id = NULL,
  candidate_id = NULL)
}
\arguments{
\item{project}{Completed \code{neo_project} with materialized output.}
\item{grade}{Exact final-grade values: INVALID, NOT_GRADED, A, B, C, D, E.}
\item{prediction_status}{Exact prediction status, such as NOT_RUN, SUCCESS,
FAILED, PARTIAL; inspect primary/provider status for precise meaning.}
\item{grading_status}{Optional exact grading status.}
\item{reason}{Retained warning/downgrade/exclusion/not-evaluated/grade code.}
\item{sample_id, candidate_id}{Optional exact sample/candidate IDs.}
}
\description{Returns fresh candidate views with identity, evidence, status and
available grade/rank fields. Filters do not modify persisted results.}
\details{
Filters are nonempty character vectors: OR within a field, AND across fields.
Missing output raises an error; successful empty/filtered views have zero rows.
Nested object/array fields are list columns preserving NULL; scalar unknowns use
typed NA. Optional policies add columns/attributes, including vaccine_set,
indirect_hla_result and allo_native_evidence.

Somatic priority requires quantitative evidence and evaluated primary percentile
rank. Default Class I strong/binder thresholds are <=0.5/2, Class II <=1/5.
Eligible strong candidates with RNA support are A, without support B; binder
candidates C/D; evaluated outside-band E. Invalid technical eligibility is INVALID;
incomplete prerequisites are NOT_GRADED. Unavailable RNA does not establish absence.
Scoped lexicographic rank is not a probability or treatment recommendation.
Optional overlays preserve grade/rank.

Somatic priority projects unevaluated primary evidence to NOT_RUN even if the
underlying primary failed or was unsupported; inspect primary status records
and explanations. NOT_GRADED is not E; missing evidence is not negative.
}
\value{A data frame with core schema-shaped and optional columns/attributes.}
\examples{
source(system.file("examples", "quickstart.R", package = "Neoantimon2"), local = TRUE)
neo_results(p, prediction_status = "NOT_RUN", sample_id = "sample-1")
}
\seealso{\code{\link{neo_config}}, \code{\link{neo_project}},
\code{\link{neo_run}}, \code{\link{neo_results}}, \code{\link{neo_events}},
\code{\link{neo_origins}}, \code{\link{neo_multisample}}, \code{\link{neo_explain}}}
