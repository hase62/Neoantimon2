\name{neo_download_resources}
\alias{neo_download_resources}
\alias{neo_check_resources}
\alias{neo_check_predictors}
\alias{neo_setup_predictor}
\title{Prepare resources and configure native predictors from R}
\usage{
neo_download_resources(profile = NULL, resources = NULL, dest,
  overwrite = FALSE, python = NULL, timeout = 600)

neo_check_resources(profile = NULL, resources = NULL, dest = NULL,
  config = NULL, checksum = TRUE, validate_reference = FALSE, python = NULL)

neo_check_predictors(config = NULL, providers = NULL, python = NULL,
  smoke_test = FALSE, run_root = tempdir(), timeout = 30)

neo_setup_predictor(provider, executable = NULL, model_path = NULL,
  task = NULL, options = list(), python = NULL, environment = NULL,
  install = FALSE, version = NULL, dependencies = character(), smoke_test = FALSE)
}
\arguments{
\item{profile}{An explicit \code{"minimum"} or \code{"full"} public-resource profile.}
\item{resources}{Character vector of individual resource IDs. May be combined
with a profile; duplicates are processed once. Not a configuration resource list.}
\item{dest}{Chosen data directory. No repository checkout or fixed layout is required.}
\item{overwrite}{Explicit permission to replace invalid existing artifacts after
new content has passed validation. The default retains existing content.}
\item{python}{Optional Python executable. Minimum-profile downloads use R and
do not require Python. Advanced resources use the installed Python helper.
For MHCflurry, use the interpreter belonging to its executable/environment.}
\item{timeout}{Positive seconds for a native download or predictor probe.
Advanced acquisition can take longer; its existing helper manages the operation.}
\item{config}{Existing \code{neo_config} object to inspect; it is not modified.}
\item{checksum}{Compute observed SHA256 values. Declared and receipt checksums
are verified even when this is false. An observed hash is not publisher authentication.}
\item{validate_reference}{Explicitly run the existing full cross-resource
reference validator. This can consume substantial time/memory for full catalogs.}
\item{providers}{Native predictor IDs to inspect. Default: configured providers,
or all three native external providers when no requests are configured.}
\item{smoke_test}{Explicitly execute a small synthetic native-adapter query.
Default checks versions/configuration only. No scientific accuracy claim is made.}
\item{run_root}{Location for owned synthetic predictor-check artifacts.}
\item{provider}{One of \code{"netmhcpan"}, \code{"netmhciipan"}, \code{"mhcflurry"}.}
\item{executable}{Existing legitimately obtained predictor executable, or NULL
to use the provider's command on PATH.}
\item{model_path}{Explicit model directory, required for MHCflurry readiness.}
\item{task}{Existing provider task; default is the provider's binding task.}
\item{options}{Existing native predictor options, validated by the existing policy.}
\item{environment}{New isolated Python environment directory, only with install=TRUE.}
\item{install}{Explicit consent to create a new MHCflurry environment and install
pinned packages. Never modifies an existing or active environment. Licensed DTU
software is not downloaded or installed.}
\item{version}{Explicit MHCflurry 2.x package release for a requested installation.}
\item{dependencies}{Additional explicit \code{package==version} pins for that new
environment, for example an appropriate TensorFlow distribution for older releases.}
}
\description{Opt-in setup convenience. Existing configurations and analyses do
not require these functions. These interfaces are development additions and are
not present in the immutable v0.2.0 release.}
\details{
Minimum selects genome, GENCODE and MANE. R downloads the selected public
resources, streams decompression and uses samtools when a genome FAI must be
prepared. Other existing full/individual selections reuse the installed Python
helper and their external validators. R and Python use one packaged catalog.
No licensed predictor or human reference is bundled.

Downloads use owned staging and validate before promotion. Invalid existing
files are retained unless overwrite=TRUE. Do not concurrently acquire into one
destination through different tools. A stale setup lock requires inspection;
no unknown staging or user data is automatically removed.

Resource reports distinguish missing files, missing dependencies, invalid
content, declared/receipt SHA256 verification and observed-only hashes.
Configured reference roles follow the existing R contract. Full catalog
relationships remain NOT_RUN unless validate_reference=TRUE. File checks alone
do not establish biological compatibility or source authenticity.

Predictor checks preserve existing version/task/model contracts. A detected
version produces CONFIGURED_UNTESTED, not a prediction-validation claim.
Optional smoke checks use the existing native adapters. Indirect-HLA imports
and other import-only providers are separate. MHCflurry models must be supplied;
they are not automatically fetched. Failed new environments remain marked
incomplete. Python package/backend availability does not prove model compatibility.
}
\value{Download/check functions return lists containing file/provider diagnostic
data frames and detailed verification states. Resource reports also contain
\code{resources} and \code{scientific} fragments accepted by \code{neo_config}.
Predictor setup returns an existing-contract \code{predictors} request list,
checks and installation metadata. Missing predictors return diagnostics rather
than failing unrelated operations. Malformed arguments and failed acquisitions
raise informative errors.}
\examples{
neo_check_resources(profile = "minimum",
  dest = file.path(tempdir(), "absent-neo-resources"))$files[, c("resource_id", "status")]
neo_check_predictors(providers = "netmhcpan")$predictors
x <- neo_setup_predictor("netmhcpan",
  executable = file.path(tempdir(), "absent-netmhcpan"))
x$checks$predictors
}
\seealso{\code{\link{neo_config}}, \code{\link{neo_project}}, \code{\link{neo_run}}}
