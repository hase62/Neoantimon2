\name{neo_indirect_class2}
\alias{neo_indirect_class2}
\title{Run Expanded Research-Only Indirect HLA Class II Analysis}
\usage{
neo_indirect_class2(donor, recipient,
  profile = c("dr-expanded-v1", "dqdp-exploratory-v1"),
  installation = NULL, reference, output_dir, presenters = NULL,
  source_loci = ie_loci(), self_loci = source_loci,
  thresholds = list(DR = 1000, DQ = 1000, DP = 1000),
  inversion = c("vendor_default", "all"), dry_run = FALSE, resume = FALSE,
  batch_size = 256L, timeout = 600L, min_free_gib = 25)
}
\arguments{
\item{donor,recipient}{Named locus lists. Character values explicitly declare typed allele copies; records can state status typed/absent/untyped/uncertain, alleles, optional copy_number and evidence. No DRB gene or copy number is inferred.}
\item{profile}{Separate DR-expanded or exploratory DQ/DP research profile.}
\item{installation}{List with path to a legitimately installed 4.3k root and authorization=list(status,terms). Status is authorized_local_installation or mock_only. No model acquisition.}
\item{reference}{List(release,cache, optional download/commit) for official IPD, or pinned local list(path,sha256) manifest. References are locally generated and must not be redistributed contrary to upstream terms.}
\item{output_dir}{Fresh output directory. Existing results are never overwritten.}
\item{presenters}{List of pair, pairing and expression_confirmed records, optionally evidence. DQ/DP requires explicitly established pairs. DR auto-pairing is restricted to one typed DRA*01:01 and typed beta alleles.}
\item{source_loci,self_loci}{Selected donor and recipient precursor loci. Default is all 13 source loci; the historical five can be explicitly selected. Self scope is these selected proteins, not the human proteome.}
\item{thresholds}{Separate DR/DQ/DP BA nM assumptions; each defaults to strictly less than 1000 nM. No calibration or clinical interpretation.}
\item{inversion}{Vendor default, or explicit all-locus inversion.}
\item{dry_run}{Prepare reference, typing and windows without invoking a predictor. Offline operation requires download=FALSE and an existing cache.}
\item{resume}{Reopen completed/incomplete/dry results after exact request, code, tool, reference and artifact hash verification; no inference or window generation.}
\item{batch_size}{Maximum 1024 peptides per bounded local CPU invocation.}
\item{timeout}{Positive per-command timeout in seconds.}
\item{min_free_gib}{Free-space reserve, at least 25 GiB.}
}
\description{Optional research orchestration reusing existing 15mer enumeration, 4.3k BA parser and execution primitives. Legacy 3.2/DRB1, somatic reconstruction, Native Grade/rank, Origin identities and Yield remain unchanged.}
\details{Supports configured A/B/C, DRB1/3/4/5, DQA1/DQB1, DPA1/DPB1, MICA/MICB source precursors. Exact typed DRA--DRB and explicitly established DQ/DP pairs are distinct profiles/schemas. Missing, partial, untyped and absent observations stay distinct. Unknown group counts are NULL; conditional counts are provisional. Two-field protein-name reference mapping may use identical fully sequenced synonymous records of that exact protein group; no higher-resolution genotype or partial sequence is completed. All origins and query/candidate links remain after sequence deduplication. Affinity and BA rank remain separate from EL metrics. Inverted core offsets index the reversed ligand; self matches require the same core, orientation and presenter. No official PIRCHE/Frost equivalence or clinical performance is established. No DQ/DP pooled clinical risk is returned. IPD reference and vendor assets are prepared locally under their terms, not redistributed.}
\value{A saved run list with status, manifest, plan, result, explanation and output_dir. TSV list cells are lossless JSON; result/explanation JSON arrays are authoritative for empty tables.}
\seealso{\code{\link{neo_indirect_autorun}}}
