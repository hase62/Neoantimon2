\name{neo_multisample}
\alias{neo_multisample}
\title{Inspect subject-scoped multi-sample evidence}
\usage{neo_multisample(project, patient_id = NULL, peptide = NULL)}
\arguments{
\item{project}{Completed \code{neo_project}.}
\item{patient_id}{Optional exact subject identity/identities.}
\item{peptide}{Optional exact literal peptide sequence(s).}
}
\description{Returns multi-sample summaries or configured longitudinal trajectories
with scoped evidence and evaluability.}
\details{Versioned settings declare subjects, contexts, relations, comparable
time and sample graphs from the current project or explicitly bound completed
sample runs. Informal timepoint labels alone do not activate this workflow.
Missing observations are not negative assays.
Gain/loss requires chronology and supported absence. Imported provider prevalence
is separate from event CCF; no de novo clustering, universal trajectory score or
patient ranking is computed. Somatic priority/grade/rank are preserved.
Single-sample/matched-normal-only projects do not create a multi-tumor overlay.
Exact-vector filters never modify saved results.}
\value{A data frame; multi-tumor integration uses typed/list cells for scoped
trajectories. Otherwise the core shape is retained, possibly with zero rows.}
\examples{
source(system.file("examples", "quickstart.R", package = "Neoantimon2"), local = TRUE)
neo_multisample(p, patient_id = "synthetic-subject")
}
\seealso{\code{\link{neo_config}}, \code{\link{neo_project}},
\code{\link{neo_run}}, \code{\link{neo_results}}, \code{\link{neo_events}},
\code{\link{neo_origins}}, \code{\link{neo_multisample}}, \code{\link{neo_explain}}}
