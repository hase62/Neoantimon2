\name{neo_report_bundle}
\alias{neo_report_bundle}
\title{Complete Read-Only Event, Peptide, pMHC and Origin Reporting}
\usage{
neo_report_bundle(project, dir = NULL, analyses = list(),
  candidate_id = NULL, yield_plan = NULL, external_channels = NULL,
  reference = NULL, endpoint = NULL, assay = NULL)
}
\arguments{
\item{project}{A completed local Neoantimon2 project. Ordinary multisample projects and core, source and long-read Event namespaces are supported. Report underlying sample projects individually for longitudinal aggregate handles.}
\item{dir}{A new output directory outside the source project, repository, inputs and references. NULL returns tables in memory. Existing directories are never overwritten.}
\item{analyses}{Explicit list of existing external analysis handles or IDs belonging to this project. No predictor, import, download or reconstruction is executed.}
\item{candidate_id}{Optional exact candidate IDs bounding pMHC, Origin and linkage output. All original Event nodes remain visible with counts explicitly scoped to the selection. NULL retains every candidate and Origin, including Origins without a candidate.}
\item{yield_plan}{Optional existing \code{neo_experimental_yield()} result. Its ordered units, alternative HLA pairs, curve, targets and probability mathematics remain unchanged.}
\item{external_channels}{Optional list of named exact selectors, such as \code{list(list(provider="bigmhc-im"))}, restricting external representatives only. All prediction observations from selected analyses remain in long output. Each analysis, model, version, run, task, metric, scope, unit, direction, sample and MHC class remains a separate ranking scope. NULL reports each distinct channel separately.}
\item{reference, endpoint, assay}{Optional versioned experimental reference and exact endpoint/assay for descriptive A--D/B1/B2 evidence. The explicit reference in yield_plan is reused by default, with integrity checks. Missing reference evidence remains unavailable. Evidence never modifies probabilities.}
}
\value{A \code{neo_report_bundle} list containing \code{tables}, a versioned \code{manifest}, and optional \code{path}. Written bundles contain quoted TSV tables and a manifest with file SHA-256, row counts and column encodings. List cells use JSON. Original source records are retained in \code{source_record_json}.}
\description{Adds a complete reporting projection without changing saved scientific results or the existing \code{neo_export_bundle()} interoperability contract. Event representatives are convenience summaries; every alternative pMHC and Origin remains exported.}
\details{
Native representatives consume the existing evaluated rank within sample, class, policy and context. An unevaluated candidate with a numeric display rank is not a best evaluated candidate. External representatives use raw score direction within an exact compatible channel and analysis. Ties use lexical candidate ID. No BA/EL/IM composite or cross-class ranking is calculated.

Origin Event associations are authoritative declared generation source sets, not assertions that each member independently caused the peptide. Germline context is separately labeled. Structural and long-read namespaces remain distinct even when IDs have identical spelling. The Event-associated Origin is the existing priority Origin if associated, otherwise the lexical associated Origin. Native Grade remains the persisted global candidate Grade, not a newly inferred Event-specific Origin Grade.

The reporting contract is \code{result-report-0.1.0}. Reporting IDs are separate from biological identities. Missing observations remain missing; coverage-status records are explicitly distinguished from measurements. Original upstream decimal score strings are retained separately from the existing public numeric projection. New TSV numeric serialization uses 17 significant digits for round-trip preservation of the available R double. Experimental Yield remains a scenario or calibration according to its supplied profile; A--D categories and B subtypes do not alter it.
}
\examples{
# p is a completed saved project:
# p <- neo_project(path = "saved-project")
# r <- neo_report_bundle(p)  # read-only, no output directory
# r$tables$event_best_by_channel
# r$tables$pmhc_predictions
# r <- neo_report_bundle(p, "new-report", analyses = list(existing_analysis_id))
}
\seealso{\code{\link{neo_events}}, \code{\link{neo_results}}, \code{\link{neo_origins}}, \code{\link{neo_export_bundle}}, \code{\link{neo_experimental_yield}}}
