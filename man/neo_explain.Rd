\name{neo_explain}
\alias{neo_explain}
\title{Trace candidate results to evidence and provenance}
\usage{
neo_explain(project, candidate_id = NULL, reason = NULL, peptide = NULL,
  sample_id = NULL, format = c("structured", "text"))
}
\arguments{
\item{project}{A \code{neo_project}; candidate explanations need materialized output.
Opt-in allo graphs can be inspected after validated construction.}
\item{candidate_id}{Optional exact candidate ID(s).}
\item{reason}{Retained warning/downgrade/exclusion/not-evaluated/grade code(s).}
\item{peptide}{Optional literal peptide sequence(s).}
\item{sample_id}{Optional exact sample ID(s).}
\item{format}{Structured list (default) or concise character text.}
}
\description{Traces selected results to retained origins, source events, reasons,
reference/provider provenance and available optional evidence.}
\details{Ordinary explanations need a selector. Nonempty character vectors
select alternatives within a field and intersections across fields. Opt-in allo
also allows query-free inspection of its separate foundation/indirect/native
graph. Structured output retains more detail than text. Optional longitudinal,
vaccine-set, presentation and receptor audits appear where available.
Missing evidence is not fabricated. Reads verify integrity without regrading,
reranking or rewriting scientific output.}
\value{A \code{neo_explanation} list with query and candidate evidence records and
configured overlays, or character text for \code{format = "text"}.}
\examples{
source(system.file("examples", "quickstart.R", package = "Neoantimon2"), local = TRUE)
x <- neo_explain(p, candidate_id = results$candidate_id[[1L]])
x$candidates[[1L]]$result
neo_explain(p, candidate_id = results$candidate_id[[1L]], format = "text")
}
\seealso{\code{\link{neo_config}}, \code{\link{neo_project}},
\code{\link{neo_run}}, \code{\link{neo_results}}, \code{\link{neo_events}},
\code{\link{neo_origins}}, \code{\link{neo_multisample}}, \code{\link{neo_explain}}}
