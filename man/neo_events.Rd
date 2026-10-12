\name{neo_events}
\alias{neo_events}
\title{Inspect canonical input events}
\usage{neo_events(project, sample_id = NULL, event_id = NULL)}
\arguments{
\item{project}{Completed \code{neo_project}.}
\item{sample_id, event_id}{Optional nonempty character vectors of exact IDs.}
}
\description{Reads verified materialized event rows and returns an exact-filtered
view. An event is an input observation, not an origin or candidate.}
\details{Filters use alternatives within a field and intersections across fields.
Unexecuted output raises an error; successful empty views retain schema.
Noncoding/unresolved/SV events need not generate candidates. Source-extension
graphs and aggregate sample references can require explanations/original sample
projects rather than the core event view.}
\value{A data frame of event fields; nested objects/arrays remain list columns.}
\examples{
source(system.file("examples", "quickstart.R", package = "Neoantimon2"), local = TRUE)
neo_events(p, event_id = events$event_id[[1L]])
}
\seealso{\code{\link{neo_config}}, \code{\link{neo_project}},
\code{\link{neo_run}}, \code{\link{neo_results}}, \code{\link{neo_events}},
\code{\link{neo_origins}}, \code{\link{neo_multisample}}, \code{\link{neo_explain}}}
