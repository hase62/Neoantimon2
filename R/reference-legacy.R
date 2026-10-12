reference_legacy <- function(paths) {
  x <- reference_table(paths[["refFlat"]], 11L, "refFlat")
  names(x) <- c("gene_name", "transcript_id", "seqname", "strand", "txStart", "txEnd",
    "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds")
  if (any(!nzchar(x$gene_name)) || any(!nzchar(x$transcript_id)) || any(!nzchar(x$seqname)) ||
      anyDuplicated(x$transcript_id)) stop("Duplicate/empty refFlat identity")
  if (any(!x$strand %in% c("+", "-"))) stop("Invalid refFlat strand")
  for (n in c("txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount"))
    x[[n]] <- reference_integer(x[[n]], paste("refFlat", n))
  if (any(x$txStart >= x$txEnd | x$cdsStart > x$cdsEnd | x$cdsStart < x$txStart |
    x$cdsEnd > x$txEnd | x$exonCount < 1)) stop("Invalid refFlat interval")
  for (n in c("exonStarts", "exonEnds")) x[[n]] <- lapply(x[[n]], function(s) {
    if (!grepl("^[0-9]+(,[0-9]+)*,?$", s)) stop("Malformed refFlat exon list")
    reference_integer(strsplit(sub(",$", "", s), ",", fixed = TRUE)[[1]], "refFlat exon")
  })
  for (i in seq_len(nrow(x))) {
    s <- x$exonStarts[[i]]; e <- x$exonEnds[[i]]
    if (length(s) != x$exonCount[i] || length(e) != x$exonCount[i]) stop("refFlat exonCount mismatch")
    if (any(s >= e | s < x$txStart[i] | e > x$txEnd[i]) ||
        (length(s) > 1L && any(s[-1] < utils::head(e, -1)))) stop("Invalid refFlat exon interval/order")
    if (s[1] != x$txStart[i] || utils::tail(e, 1) != x$txEnd[i]) stop("refFlat transcript/exon bounds mismatch")
  }
  tx <- reference_fasta(paths[["refMrna"]], "legacy")
  if (any(!x$transcript_id %in% tx$id)) stop("Missing matching refMrna transcript")
  expected <- vapply(seq_len(nrow(x)), function(i) sum(x$exonEnds[[i]] - x$exonStarts[[i]]), 0)
  mismatches <- x$transcript_id[expected != tx$length[match(x$transcript_id, tx$id)]]
  # RefSeq/genomic differences are evidence, never silently repaired.
  list(annotation = reference_sort(x), transcripts = tx,
    coordinate_convention = "0-based half-open", observed = list(assembly = NULL),
    mismatches = list(transcript_exonic_length = mismatches,
      fasta_without_annotation = setdiff(tx$id, x$transcript_id)),
    capabilities = list(genome_sequence_available = FALSE, transcript_sequence_available = TRUE,
      protein_reference_available = FALSE, mane_available = FALSE, sec_metadata_available = FALSE,
      assembly_verified = NA), mitochondrial_contigs = character())
}
