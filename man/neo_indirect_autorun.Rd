\name{neo_indirect_autorun}
\alias{neo_indirect_autorun}
\title{Run the literature-derived indirect-HLA Class II workflow}
\usage{
neo_indirect_autorun(allo, presenters, reference, executable = NULL,
  model = NULL, output_dir, dry_run = FALSE, resume = FALSE,
  authorization = NULL, timeout = 600L, batch_size = 256L)
}
\arguments{
\item{allo}{Existing allo-foundation configuration, or a list of donor and recipient typing lists. The compact form requires two explicit allele copies at A, B, C, DRB1 and DQB1, plus recipient DRA.}
\item{presenters}{Explicit recipient DRA*01:01~DRB1 pairs, with two-field DRB1 names, or existing presenter declaration records. No inferred alpha chains or fuzzy substitution.}
\item{reference}{List with path and SHA-256 of an explicit versioned literal precursor hla-reference-1.0.0 manifest.}
\item{executable}{Legitimately installed local NetMHCIIpan 3.2 command. Never downloaded automatically.}
\item{model}{List with path (installed -rdir root containing data/), identity, and optional expected tree SHA-256. Internal symlinks are rejected.}
\item{output_dir}{New output directory; existing files are preserved.}
\item{dry_run}{Generate query and provenance plans without executing any provider command.}
\item{resume}{Reopen an identical completed run or dry plan with strict request/tool/artifact checks. Interrupted or failed runs cannot resume in place.}
\item{authorization}{Explicit list with status and terms. Status is authorized_local_installation or mock_only. This records a user declaration; it does not obtain a license.}
\item{timeout}{Positive per-command timeout, in seconds.}
\item{batch_size}{Positive maximum number of unique peptides per presenter invocation.}
}
\description{Optional local orchestration reusing the existing allo foundation, 15mer generation, self catalog, imported-binding contract and strict literature-derived counting. It is not official PIRCHE reproduction or a clinical risk model.}
\details{First scope is human SOT donor-to-recipient response, source/self loci A/B/C/DRB1/DQB1 and explicit recipient DRA~DRB1. NetMHCIIpan 3.2 BA is required. Other versions, DQ/DP presentation and alternate source directions are not silently substituted. Primary nine-residue cores and zero-based offsets must be present and consistent in provider output. Binding uses the existing strictly-less-than-1000 nM threshold and same-presenter primary-core self policy. Missing/incomplete evidence is never converted into a definitive zero.

The model tree and executable are hashed before and after prediction. Help, supported allele list, exact quoted argument vectors, stdout/stderr, query links, normalized evidence, public result/explanation and SHA-256 lists are saved. The same output directory is never overwritten. A completed run can be reopened only with resume=TRUE, matching inputs/tool identity and unchanged artifacts. No weights, binaries or reference resources are acquired by the API.

The bundled mock is a synthetic software fixture. Its provenance explicitly states that no real predictor was executed, while recording the actual mock process. See the installed example and repository user guide for setup and real-provider limitations.}
\value{A list with status, output_dir, manifest, plan and (after completion) existing public result and explanation projections.}
\examples{
ref <- system.file("extdata", "indirect-autorun", "synthetic-reference.json",
                   package = "Neoantimon2")
a <- jsonlite::read_json(system.file("extdata", "indirect-autorun",
  "synthetic-foundation.json", package = "Neoantimon2"))
p <- "HLA-DRA*01:01~HLA-DRB1*01:01"
plan <- neo_indirect_autorun(a, p,
  list(path = ref, sha256 = digest::digest(file = ref, algo = "sha256",
    serialize = FALSE)), output_dir = tempfile("indirect-plan-"),
  dry_run = TRUE)
plan$status
}
