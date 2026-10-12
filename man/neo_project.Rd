\name{neo_project}
\alias{neo_project}
\title{Create an analysis project or reopen a saved run}
\usage{
neo_project(samples = NULL, hla = NULL, reference = NULL, config = NULL,
  path = NULL, mode = "somatic", project_id = "project",
  analysis_context = NULL, inputs = list(), initialize = FALSE,
  fixture_kind = "real")
}
\arguments{
\item{samples}{Data frame or array of records with unique sample_id, explicit
patient_id subject identity and optional timepoint/tissue/role strings.}
\item{hla}{Named list mapping subjects to declared allele vectors; declare a
subject once. Class II requires explicit alpha~beta pairs.}
\item{reference}{Named resource records overriding config references; role names
receive the reference. prefix.}
\item{config}{Validated \code{neo_config}; defaults are used if omitted.}
\item{path}{Existing run directory to verify/reopen, or a new path inside run_root
to initialize. Relative paths require explicit config.}
\item{mode}{Currently only somatic; allo workflows use configuration opt-ins.}
\item{project_id}{Descriptive label, excluded from scientific identity.}
\item{analysis_context}{Optional context records; defaults bind each subject's
samples to one somatic context.}
\item{inputs}{Declarations with input_id, input_role, sample_id, path, optional
version and checksum. File checksum declares algorithm sha256, scope file, and
a 64-character lowercase hexadecimal value. Roles choose supported parsers.}
\item{initialize}{Create an execution directory immediately; default FALSE.}
\item{fixture_kind}{Provenance category real/synthetic/mock. Default real enables
ordinary analysis; synthetic/mock runtime fixtures differ from running the
ordinary pipeline on explicitly synthetic input bytes.}
}
\description{Binds metadata, inputs/references and configuration to an analysis
identity, or reopens a verified run. No ordinary variant reconstruction or
external predictor execution occurs at construction.}
\details{
File integrity observations and opt-in allo manifest validation occur during
construction. Initialize/path creation writes a manifest, logs, checkpoints,
intermediates, cache and results directories. Original input/resource roots must
not overlap generated run trees.

Scientific identity includes metadata/context, input/resource identity and
scientific settings; runtime options, descriptive input labels and output
locations are excluded. Without content checksums, paths/stat metadata contribute
conservatively to input identity. Each execution has a distinct run ID. Reopening
uses saved config and verifies identity, safe paths and integrity. Changed
scientific inputs/settings require a new project. Locks are not automatically
stolen after interruption.
}
\value{A \code{neo_project} S3 list with resolved config, inputs/resources,
metadata, analysis identity and optional run identity/path. Manifest is authoritative.}
\examples{
p <- neo_project(samples = data.frame(sample_id = "s1", patient_id = "subject-1"),
  hla = list("subject-1" = "HLA-A*02:01"), config = neo_config(env = character()))
class(p)
}
\seealso{\code{\link{neo_config}}, \code{\link{neo_project}},
\code{\link{neo_run}}, \code{\link{neo_results}}, \code{\link{neo_events}},
\code{\link{neo_origins}}, \code{\link{neo_multisample}}, \code{\link{neo_explain}}}
