\name{neo_config}
\alias{neo_config}
\title{Load and validate analysis configuration}
\usage{neo_config(file = NULL, ..., env = Sys.getenv())}
\arguments{
\item{file}{Optional JSON/YAML configuration filename, absolute or home-relative.}
\item{...}{Named overrides: \code{paths}, \code{scientific}, \code{runtime},
\code{checksum}, \code{resources}, \code{peptides}, \code{predictors},
\code{somatic_evidence}, \code{somatic_priority}, \code{hla_loh},
\code{vaccine_set}, \code{allo}, \code{indirect_hla}, \code{allo_native_evidence}.
Unknown keys and invalid types are rejected.}
\item{env}{Named character environment vector. Use \code{character()} to disable
environment overrides in a saved recipe.}
}
\description{Combines defaults, a JSON/YAML file, environment settings and explicit
arguments, in that order. Validates configuration without running an analysis
or installing providers.}
\details{
Paths include data_root, input_root, tools_root, run_root and optional repository.
File-relative roots use the config directory; argument/environment roots must
resolve absolutely. Defaults use R's per-user data directory. Relative references
use data_root, inputs use input_root, logical tools/models use tools_root.
No developer directory layout is required.

Resources are named records with path, required, type, resource_id, version,
checksum and source. Types are file/directory/executable. Checksum modes are
none, metadata (default), sha256. Directory SHA256 hashing is rejected.
A checksum declaration and an enabled verification mode are distinct.

Scientific settings include species taxon, genome assembly and optional policy ID.
Runtime settings include positive integral threads and verbosity quiet/info/debug.
Environment prefix is NEOANTIMON2_ with suffix DATA_ROOT, RUN_ROOT, INPUT_ROOT,
TOOLS_ROOT, REPOSITORY, THREADS, VERBOSITY or CHECKSUM.

Optional policies enable workflows only when explicitly bound. Vaccine-set policy
needs version vaccine-set-policy-1.0.0 and explicit max_peptides; it may be a list
or JSON/YAML filename. Allo and indirect-HLA use separate versioned manifests.
Missing evidence is not a negative observation.
}
\value{A validated \code{neo_config} S3 list. Printing omits paths and values.}
\examples{
cfg <- neo_config(env = character(), runtime = list(threads = 2))
class(cfg)
}
\seealso{\code{\link{neo_config}}, \code{\link{neo_project}},
\code{\link{neo_run}}, \code{\link{neo_results}}, \code{\link{neo_events}},
\code{\link{neo_origins}}, \code{\link{neo_multisample}}, \code{\link{neo_explain}}}
