long_read_coordinates <- function(s) as.integer(unlist(lapply(s$exons, function(e)
  if (s$strand == "+") seq.int(e$start, e$end) else seq.int(e$end, e$start)),
  use.names = FALSE))
long_read_reconstruct <- function(s, reference) {
  pieces <- vapply(s$exons, function(e) {
    seq <- event_reference_fetch(reference, s$contig, e$start, e$end)
    if (s$strand == "-") reconstruction_reverse_complement(seq) else seq
  }, "")
  seq <- paste0(pieces, collapse = "")
  if (!grepl("^[ACGT]+$", seq)) stop("Long-read reference sequence contains unsupported bases")
  ends <- cumsum(nchar(pieces)); starts <- c(1L, utils::head(ends, -1L) + 1L)
  list(sequence = seq, coordinates = long_read_coordinates(s),
    contigs = rep(s$contig, nchar(seq)), strand = s$strand,
    blocks = lapply(seq_along(pieces), function(i) c(s$exons[[i]], list(
      start_nt = as.integer(starts[[i]]), end_nt = as.integer(ends[[i]]),
      contig = s$contig, strand = s$strand))),
    junction_nt = as.integer(utils::head(ends, -1L)),
    sequence_authority = "REFERENCE_GENOME", reconstruction_policy = "long-read-reconstruction-v1")
}

long_read_reference_models <- function(bundle, structures, reference) {
  if (!length(structures)) return(list())
  a <- bundle$annotation
  tx <- a[a$feature == "transcript" & !is.na(a$transcript_id), , drop = FALSE]
  selected <- unique(unlist(lapply(structures, function(s) tx$transcript_id[
    tx$seqname == s$contig & tx$strand == s$strand & tx$start <=
      max(unlist(lapply(s$exons, `[[`, "end"))) & tx$end >=
      min(unlist(lapply(s$exons, `[[`, "start")))])))
  f <- reference_transcript_features(bundle, sort(selected))
  lapply(split(f, f$transcript_id), function(rows) {
    exon <- rows[rows$feature == "exon", , drop = FALSE]
    if (!nrow(exon)) return(NULL)
    s <- long_read_structure(bundle$profile$assembly_release, exon$seqname[[1L]],
      exon$strand[[1L]], exon[, c("start", "end"), drop = FALSE])
    rec <- long_read_reconstruct(s, reference)
    start <- rows[rows$feature == "start_codon", c("start", "end"), drop = FALSE]
    start_coords <- if (nrow(start)) translation_coordinates(start, s$strand) else integer()
    p <- if (length(start_coords) == 3L) match(start_coords[[1L]], rec$coordinates) else NA_integer_
    cds <- rows[rows$feature == "CDS", , drop = FALSE]
    incomplete <- any(grepl("cds_start_NF", rows$raw_attributes, fixed = TRUE))
    valid <- !is.na(p) && p + 2L <= length(rec$coordinates) &&
      identical(rec$coordinates[p + 0:2], as.integer(start_coords)) &&
      substring(rec$sequence, p, p + 2L) == "ATG" && nrow(cds) > 0L && !incomplete
    list(transcript_id = rows$transcript_id[[1L]], gene_id = exon$gene_id[[1L]],
      structure = s, reconstructed = rec,
      coding_start = if (valid) as.integer(p) else NULL,
      coding_start_coordinates = if (valid) as.integer(start_coords) else integer(),
      sec_aa = as.integer(bundle$sec$position[bundle$sec$transcript_id == rows$transcript_id[[1L]]]))
  })
}

long_read_intervals <- function(positions, coordinates, s) {
  if (!length(positions)) return(list())
  groups <- split(positions, cumsum(c(TRUE, diff(positions) != 1L)))
  lapply(groups, function(p) list(start_nt = as.integer(min(p)), end_nt = as.integer(max(p)),
    reference_coordinates = as.list(coordinates[p]), contig = s$contig, strand = s$strand))
}
long_read_junction_key <- function(j) paste(j$donor, j$acceptor, sep = ":")

long_read_compare <- function(s, rec, references) {
  refs <- Filter(function(r) !is.null(r) && r$structure$contig == s$contig &&
    r$structure$strand == s$strand &&
    min(r$reconstructed$coordinates) <= max(rec$coordinates) &&
    max(r$reconstructed$coordinates) >= min(rec$coordinates), references)
  exact <- Filter(function(r) identical(s, r$structure), refs)
  # Exact structural match to any reference overrides relative differences from
  # other isoforms; alternative AUGs on an FSM cannot become novel sources.
  catalog <- unique(unlist(lapply(refs, function(r) vapply(r$structure$junctions,
    long_read_junction_key, ""))))
  donors <- unique(unlist(lapply(refs, function(r) vapply(r$structure$junctions, `[[`, 0L, "donor"))))
  acceptors <- unique(unlist(lapply(refs, function(r) vapply(r$structure$junctions, `[[`, 0L, "acceptor"))))
  novel_junction <- which(!vapply(s$junctions, long_read_junction_key, "") %in% catalog)
  evidence <- lapply(seq_along(s$junctions), function(i) c(s$junctions[[i]], list(
    junction_nt = rec$junction_nt[[i]], status = if (i %in% novel_junction)
      "novel_junction" else "known_junction",
    donor_status = if (s$junctions[[i]]$donor %in% donors) "known_donor" else "novel_donor",
    acceptor_status = if (s$junctions[[i]]$acceptor %in% acceptors) "known_acceptor" else "novel_acceptor")))
  reference_exonic_coordinates <- unique(unlist(lapply(refs, function(r) r$reconstructed$coordinates)))
  inserted_positions <- if (length(exact) || !length(refs)) integer() else
    which(!rec$coordinates %in% reference_exonic_coordinates)
  novel_intervals <- long_read_intervals(inserted_positions, rec$coordinates, s)
  comparisons <- lapply(refs, function(r) {
    mapped <- match(rec$coordinates, r$reconstructed$coordinates)
    inserted <- inserted_positions
    missing <- setdiff(r$reconstructed$coordinates, rec$coordinates)
    cls <- if (identical(s, r$structure)) "reference-compatible_structure" else character()
    if (!length(exact)) {
      if (length(novel_junction)) cls <- c(cls, "novel_junction")
      if (any(vapply(evidence, function(j) j$donor_status == "novel_donor", TRUE))) cls <- c(cls, "novel_donor")
      if (any(vapply(evidence, function(j) j$acceptor_status == "novel_acceptor", TRUE))) cls <- c(cls, "novel_acceptor")
      if (length(inserted)) {
        for (e in s$exons) {
          p <- which(rec$coordinates >= e$start & rec$coordinates <= e$end)
          overlap <- any(!is.na(mapped[p]))
          internal <- min(rec$coordinates[p]) > min(r$reconstructed$coordinates) &&
            max(rec$coordinates[p]) < max(r$reconstructed$coordinates)
          if (all(is.na(mapped[p])) && internal) cls <- c(cls, "exon_insertion")
          if (overlap && any(is.na(mapped[p])) && internal) cls <- c(cls, "retained_intronic_block")
        }
        # A merged exon spanning a full annotated intron is retained intron.
        for (j in r$structure$junctions) if (any(vapply(s$exons, function(e)
          e$start <= min(j$donor, j$acceptor) && e$end >= max(j$donor, j$acceptor), TRUE)))
          cls <- c(cls, "retained_intronic_block")
      }
      if (length(missing) && any(vapply(r$structure$exons, function(e)
          all(!seq.int(e$start, e$end) %in% rec$coordinates), TRUE))) cls <- c(cls, "exon_skipping")
      if (s$transcript_start != r$structure$transcript_start ||
          s$transcript_end != r$structure$transcript_end) cls <- c(cls, "alternative_terminal_structure")
      if (!length(cls)) cls <- "known_junction_combination"
    } else if (!length(cls)) cls <- "alternative_reference_comparison"
    list(transcript_id = r$transcript_id, gene_id = r$gene_id,
      exact_structure = identical(s, r$structure), classes = as.list(sort(unique(cls))),
      novel_nt_intervals = long_read_intervals(inserted, rec$coordinates, s),
      reference_coordinate_map = mapped, reference = r)
  })
  # Multi-gene models whose blocks have disjoint reference gene assignments are
  # deferred, never translated as an ordinary isoform or a second fusion model.
  block_genes <- lapply(s$exons, function(e) unique(vapply(Filter(function(r)
    any(r$reconstructed$coordinates >= e$start & r$reconstructed$coordinates <= e$end), refs),
    `[[`, "", "gene_id")))
  known <- Filter(function(g) length(g) > 0L, block_genes)
  fusion_like <- length(known) > 1L && !length(Reduce(intersect, known))
  list(status = if (fusion_like) "FUSION_MAPPING_DEFERRED" else if (length(exact))
    "REFERENCE_EQUIVALENT" else if (length(refs)) "COMPARED" else "UNRESOLVED_COMPARISON",
    comparisons = comparisons, junctions = evidence, novel_nt_intervals = novel_intervals,
    novel_junction_nt = if (length(exact) || !length(refs)) integer() else
      rec$junction_nt[novel_junction], reference_equivalent = length(exact) > 0L,
    comparison_policy = "long-read-reference-comparison-v1")
}

long_read_fasta_comparison <- function(sequence, provider_sequence = NULL, supplied = FALSE) {
  status <- if (!supplied) "not_supplied" else if (is.null(provider_sequence)) "unresolved" else
    if (nchar(sequence) != nchar(provider_sequence)) "length_difference" else
      if (identical(sequence, provider_sequence)) "exact_match" else "sequence_difference"
  list(status = status, canonical_sequence_authority = "REFERENCE_GENOME",
    provider_sequence_checksum = if (is.null(provider_sequence)) NULL else peptide_digest(provider_sequence),
    canonical_sequence_checksum = peptide_digest(sequence),
    mismatch_nt = if (status == "sequence_difference") as.list(which(
      strsplit(sequence, "", fixed = TRUE)[[1L]] != strsplit(provider_sequence, "", fixed = TRUE)[[1L]])) else list(),
    mutation_inference = FALSE)
}
