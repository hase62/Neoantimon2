\name{neo_experimental_yield}
\alias{neo_experimental_yield}
\alias{plot.neo_experimental_yield}
\title{Optional Experimental Yield Estimation v0.1}
\description{Assumption-based planning of the number of candidate experimental
units needed to detect at least one positive result. Independent A--D evidence
labels, including B1 and B2, never modify assigned probabilities. This optional
read-only view does not modify native Grades, ranking, scores or Origin links.}
\usage{
neo_experimental_yield(x, profile, reference = NULL, analysis = NULL,
  channel = NULL, unit = c("peptide", "peptide_hla"),
  ordering = c("score", "input"),
  primary_hla = c("best_score", "lexicographic"),
  min_records = 100L, min_patients = 3L)

\method{plot}{neo_experimental_yield}(x, y = NULL, ...)
}
\arguments{
\item{x}{Completed \code{neo_project}, or a data frame returned by
\code{neo_priorities}. For plotting, a returned yield object. A supplied table
must retain sample/candidate/peptide/HLA identity, raw score, status, source,
provider, predictor/model versions, run, task, metric, metric kind/semantics,
unit and direction. Exactly one sample and one selected channel are required.}
\item{profile}{Explicit versioned \code{experimental-yield-profile-0.1.0} JSON
path or list. The identifier \code{"improve-exploratory-v0.1"} loads an
illustrative scenario with assumed rates 0.0139, 0.0230, 0.0306, 0.0639.
Those IMPROVE-derived owner-supplied reanalysis rates are not a BigMHC calibration.}
\item{reference}{Optional versioned \code{experimental-reference-0.1.0} JSON
path or list. No implicit network download. Missing references produce an explicit
unavailable status, not a fabricated D category.}
\item{analysis}{Existing external analysis handle or ID, for project input only.}
\item{channel}{Named list of exact identity-field values selecting one
predictor/model/run/task/metric channel. BA, EL and IM are never combined.
Unpredicted candidates lacking channel identity remain unavailable pool members.}
\item{unit}{\code{peptide}: each distinct synthesized sequence uses one
explicit primary-HLA probability proxy. \code{peptide_hla}: exact assays;
repeated Origins are deduplicated, and synthesized sequences are separately counted.}
\item{ordering}{Direction-aware score order with deterministic sequence/HLA ties,
or the first input occurrence. Existing native ranking is never changed.}
\item{primary_hla}{Highest direction-aware score, then lexical full HLA and
candidate ID; or lexical HLA. All alternative associations and original records
are retained. Best-score selection is an optimistic, unvalidated peptide proxy.}
\item{min_records, min_patients}{Provisional engineering thresholds for C;
these do not alter any probability. See \code{neo_yield_evidence}.}
\item{y}{Unused, for the plot generic.}
\item{...}{Additional arguments passed to the base plot method.}
}
\details{
There are four bins ordered by prediction quality: bottom 30\%, 30--60\%,
60--90\%, top 10\%. Fixed profiles declare three raw-score thresholds in quality
order; equality stays in the lower bin. Lower-is-better metrics use decreasing
raw cutpoints. Thresholds ordinarily come from a frozen reference distribution.
Scenario cohort-relative mid-rank bins are explicitly exploratory and candidate-set
dependent; ties can leave bins empty. Origin duplicates do not affect bins.

Scenario rates are assumptions. Cross-model use is explicitly labelled and only
permitted using exploratory relative ranks, not another model's raw cutpoints.
Empirical profiles require exact provider, predictor/model versions, task, metric,
kind, semantics, units and direction, frozen cutpoints, endpoint/assay, population,
limitations, provenance and four tested/positive/negative counts with matching rates.
A matching user-supplied profile does not establish independent validation or
population transportability. No empirical profile is shipped as validated.

For the selected ordered units, expected positives equal \code{sum(p)}.
Success probability equals \code{-expm1(sum(log1p(-p)))} under independence.
No probability is changed by evidence class, VAF, expression or native Grade.
Missing predictions are not measured zero or negative; full-pool quantities become
unavailable and scored-subset quantities are separately named. No extra candidates
are created for unreachable targets. A missing probability in an input-ordered
prefix makes that prefix unavailable.

Multiple HLA assays sharing the same peptide are not silently treated as independent
successes. Expected positive assay counts remain additive; joint success and target
values become unavailable after the first repeated peptide. Use the single-HLA
peptide proxy for a simple synthesis plan. Even different peptides can have
correlated outcomes; the independence model is an explicit simplifying assumption.
pMHC-multimer response rates do not establish ELISpot or tumor recognition rates.
}
\value{A \code{neo_experimental_yield} list: \code{summary}, \code{units},
\code{pair_evidence}, \code{curve}, \code{targets} for 50/80/90/95\%,
\code{composition}, full profile/reference, original input provenance, assumptions
and warnings. List columns retain all original candidate/Origin/source records and
alternative HLA pairs. \code{reachable} is FALSE when known rates cannot reach a
target, or NA when missing/joint evidence prevents deciding; \code{minimum_N} is
NA in either unavailable case, with an explicit reason. The plot method returns
the curve invisibly.}
\seealso{\code{\link{neo_yield_evidence}}, \code{\link{neo_priorities}},
\code{\link{neo_external_results}}}
\examples{
source(system.file("examples", "experimental-yield.R", package = "Neoantimon2"))
plan <- experimental_yield_example() # synthetic offline example
plan$targets
}
