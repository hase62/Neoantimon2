\name{neo_run}
\alias{neo_run}
\title{Execute or resume a configured analysis}
\usage{neo_run(project)}
\arguments{\item{project}{A \code{neo_project} object.}}
\description{Initializes an execution directory and runs configured reconstruction,
peptide/candidate generation, supported predictors and optional evidence workflows.
Verified checkpoints can be reused on resume.}
\details{
Ordinary somatic analysis needs supported event input and selected references;
HLA resolution is needed for peptide-presenter candidates. Optional quantitative,
source, presentation, recognition, longitudinal, vaccine-set and allo inputs retain
their own evaluability and provider restrictions. Not every imported provider
has an executable adapter.

Writes manifests, checkpoints, intermediates, caches, results and logs under run_root.
Required input/reference changes and corrupted artifacts are not silently accepted.
Run directories may be retained and reopened; regeneration uses original inputs
and configuration rather than accidental historical development runs.

No configured predictor yields NOT_RUN; absent grading evidence yields NOT_GRADED.
Provider failure and run execution failure are separate levels. Missing evidence
is not zero or a biological negative. Overlays do not rewrite somatic grade/rank.
Indirect-HLA is literature-derived, not official PIRCHE.
}
\value{The initialized \code{neo_project} with saved path and run identity.}
\examples{
# SYNTHETIC / TEST-SCALE / NON-CLINICAL; no network or licensed tools.
source(system.file("examples", "quickstart.R", package = "Neoantimon2"), local = TRUE)
reopened <- neo_project(path = p$path)
stopifnot(identical(neo_results(neo_run(reopened)), results))
}
\seealso{\code{\link{neo_config}}, \code{\link{neo_project}},
\code{\link{neo_run}}, \code{\link{neo_results}}, \code{\link{neo_events}},
\code{\link{neo_origins}}, \code{\link{neo_multisample}}, \code{\link{neo_explain}}}
