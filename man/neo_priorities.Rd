\name{neo_priorities}
\alias{neo_priorities}
\title{Additional Model-Specific Numerical Priorities}
\usage{
neo_priorities(project, analysis = NULL,
  method = c("model-v1", "neoantimon-paper-2020-v1"),
  vaf = c("dna", "rna"), wild_type_predictions = NULL,
  expression_unit = NULL)
}
\arguments{
\item{project}{A completed Neoantimon2 project.}
\item{analysis}{Optional completed external analysis handle or ID.}
\item{method}{Versioned numerical scoring method.}
\item{vaf}{DNA or RNA measurement for adjusted priorities.}
\item{wild_type_predictions}{Optional standard TSV for exact patient-normal peptide/HLA predictions with matching nonempty run identity.}
\item{expression_unit}{Explicit \code{"TPM"} acknowledgement for paper scoring.}
}
\value{Data frame with one row per candidate, origin and prediction channel,
including raw scores, provenance list columns, prediction-only ordering,
and separately labelled evidence-adjusted values.}
\description{Provides additional read-only numerical views. It never changes
native Grades, native ranks or external analyses. See the source evidence user guide.}
\details{Missing VAF remains an unavailable biological observation and contributes
zero only to the adjusted calculation. Raw model scores are only comparable within
the same provider, model, run and metric. The adjusted empirical-order-times-VAF
value is exploratory and candidate-set dependent, not calibrated or comparable
across cohorts. Disagreeing expression assertions are retained with unavailable
TPM/conflict status, never resolved by selecting the largest value. Paper CCF singularities return unavailable
scores; finite negative/unstable values are retained with diagnostic status.
The paper method needs matching wild-type BA predictions, including the same
nonempty run ID and model/metric semantics, and expression. Cross-run pairing
is unavailable; imported IDs remain user provenance assertions, not proof of execution;
without CCF, its separately labelled no-CCF value is not a complete paper score.}
\references{Hasegawa et al. (2020), Bioinformatics 36:4813--4816,
\doi{10.1093/bioinformatics/btaa616}.}
