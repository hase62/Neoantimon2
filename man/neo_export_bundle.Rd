\name{neo_export_bundle}
\alias{neo_export_bundle}
\title{Export a verified literal peptide and annotation bundle}
\usage{neo_export_bundle(project, dir = NULL)}
\arguments{
\item{project}{A completed \code{neo_project}.}
\item{dir}{New export directory; defaults to results/export_bundle in the run.}}
\description{Exports distinct literal AA20+U peptide and canonical peptide-HLA queries,
every candidate, all origins and links, events and complete provenance JSON.
Native runs invoke this operation after completion without needing a predictor.}
\value{A bundle handle with path, manifest and annotation graph. Repeating the
operation verifies the existing immutable bundle and its source run checksums.}
\details{Annotation TSV array/object cells are JSON. Source and long-read
observations, events and peptide provenance are retained in annotations.json.
Export refuses existing incomplete targets and overlap with original inputs.
No biological ID or native result is changed.}
