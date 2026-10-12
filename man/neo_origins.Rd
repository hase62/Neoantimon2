\name{neo_origins}
\alias{neo_origins}
\title{Inspect peptide sequence origins}
\usage{neo_origins(project, sample_id = NULL, origin_id = NULL, peptide = NULL)}
\arguments{
\item{project}{Completed \code{neo_project}.}
\item{sample_id, origin_id}{Optional nonempty character vectors of exact IDs.}
\item{peptide}{Optional exact literal peptide sequence(s).}
}
\description{Returns materialized peptide origin/provenance rows. Multiple
transcript/sequence/phase paths can support one peptide-presenter candidate.}
\details{Filters use alternatives within a field and intersections across fields.
Unexecuted output raises an error; successful empty views retain schema.
Uncertainty stays explicit. Source/long-read origins are added when long-read
inputs are bound; extension graphs are also available through explanations.
Aggregate longitudinal projects retain individual sample graphs; inspect those
original projects for their origin tables.}
\value{A data frame with origin identity, sequence/transcript/phase and evidence
fields, including list columns. No stored graph is mutated.}
\examples{
source(system.file("examples", "quickstart.R", package = "Neoantimon2"), local = TRUE)
neo_origins(p, peptide = origins$peptide_sequence[[1L]])
}
\seealso{\code{\link{neo_config}}, \code{\link{neo_project}},
\code{\link{neo_run}}, \code{\link{neo_results}}, \code{\link{neo_events}},
\code{\link{neo_origins}}, \code{\link{neo_multisample}}, \code{\link{neo_explain}}}
