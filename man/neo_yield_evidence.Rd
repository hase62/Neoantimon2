\name{neo_yield_evidence}
\alias{neo_yield_evidence}
\title{Independent Experimental Evidence Categories A--D, Including B1/B2}
\description{Descriptive versioned reference lookup, independent of score bins
and probabilities. Retains all matching source records, patients, studies,
assay types, HLA identities, citations and positive/negative/conflicting results.}
\usage{
neo_yield_evidence(x, reference = NULL, endpoint = NULL, assay = NULL,
  min_records = 100L, min_patients = 3L)
}
\arguments{
\item{x}{Data frame with literal \code{peptide_sequence} and \code{mhc_key}.
Extra fields and original strings are retained. One output row per input row.}
\item{reference}{Explicit \code{experimental-reference-0.1.0} JSON path or
list with id, version, provenance, patient identity policy and experimental records.
NULL or an absent file is explicitly unavailable. Malformed files fail.}
\item{endpoint, assay}{One exact endpoint and assay. Required with an available
reference. Other assays/endpoints are excluded from classification and remain in
reference provenance; they are never silently pooled.}
\item{min_records}{Minimum independent tested peptide/patient/HLA units for C;
default100. Technical duplicates and conflicting replicates do not inflate this count.}
\item{min_patients}{Minimum independently identified patients for C; default3.
Unknown patients do not count. Study-local identities use the largest number of
distinct patients within any single study, avoiding cross-study identity assumptions.
A reference explicitly declaring verified global identity can count globally.}
}
\details{
Priority is A, then B, then C, then D. A is exact normalized peptide and full exact
HLA with a tested positive or negative result. B1 is exact peptide with a different
exact HLA key. B2 is exact HLA and exactly one amino-acid substitution at the same
length (Hamming distance1). Positive, negative and conflicting observations are
retained for both subtypes. B1 and B2 may both apply and have no validated biological
ordering. Even if A takes priority, related B matches remain available.

C requires the declared provisional tested-unit and patient thresholds for the
same endpoint/assay and HLA. D indicates limited reference evidence, not low
immunogenicity. A missing reference or missing endpoint/assay scope produces NA
category plus an unavailable status. An untested peptide is never experimentally
negative. Peptides are uppercased and trimmed without I/L substitution or wildcard
matching. HLA prefix/star punctuation is normalized conservatively, but fields and
expression suffixes are never truncated, inferred or fuzzily substituted. Class II
pairs require explicit compatible alpha/beta chains.

Every record requires record id, peptide, HLA, endpoint, assay, positive/negative
result, study id and citation. Patient id can be missing but cannot contribute to C.
Reference ids and versions are explicit. All matched original measurements are retained.
Assay-specific conflicts across people are not proof of a contradiction within a
person. Prior positivity is not a guarantee in another patient, and a one-residue
neighbor does not establish shared TCR recognition. These categories never assign
or multiply probability.
}
\value{Data frame with evidence_category, evidence status, both B subtypes, separate
A/B1/B2 result status, HLA tested/patient counts and full matched records in list
columns. Reference artifact hashes and the policy are retained as attributes.}
\seealso{\code{\link{neo_experimental_yield}}}
\examples{
neo_yield_evidence(data.frame(peptide_sequence = "ACDEFGHIK",
                            mhc_key = "HLA-A*02:01")) # reference unavailable
}
