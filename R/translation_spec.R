# T009 consumes T008 nucleotide contexts. All positions below are one-based in
# the spliced transcript, with separate original genomic coordinate maps.
translation_policy <- function() "t009-translation-spec-v1"

translation_coordinates <- function(segments, strand) {
  if (!nrow(segments)) return(integer())
  order <- if (identical(strand, "-")) rev(seq_len(nrow(segments))) else seq_len(nrow(segments))
  unlist(lapply(order, function(i) if (identical(strand, "-"))
    seq.int(segments$end[[i]], segments$start[[i]]) else
      seq.int(segments$start[[i]], segments$end[[i]])), use.names = FALSE)
}

translation_specification <- function(model, context, bundle, protein_catalog,
                                      transcript_catalog, features) {
  basic <- list(status = "UNRESOLVED", reason = NULL,
    transcript_id = context$transcript_id, protein_id = context$protein_id,
    gene_id = context$gene_id, reference_source = "GENCODE_PRIMARY_GTF_ALL_FASTA",
    reference_identity = bundle$resource_identity, strand = context$strand,
    contig = context$contig, reconstruction_policy = context$reconstruction_policy,
    translation_policy = translation_policy(),
    translation_system = if (identical(context$contig, "chrM"))
      "vertebrate_mitochondrial" else "nuclear_standard",
    genetic_code_id = if (identical(context$contig, "chrM")) 2L else 1L)
  fail <- function(reason) { basic$reason <- reason; basic }
  if (is.null(model) || !is.null(model$model_reason) ||
      is.null(context$reference_transcript) || is.null(context$reference_cds))
    return(fail(if (is.null(model)) "REFERENCE_MODEL_MISSING" else
      if (!is.null(model$model_reason)) model$model_reason else "REFERENCE_SEQUENCE_MISSING"))
  tx <- context$reference_transcript
  if (is.null(transcript_catalog[[context$transcript_id]]) ||
      !identical(tx, transcript_catalog[[context$transcript_id]]))
    return(fail("REFERENCE_TRANSCRIPT_MISMATCH"))
  ex <- translation_coordinates(model$exons, model$strand)
  cds <- translation_coordinates(model$segments, model$strand)
  if (!length(ex) || !length(cds) || nchar(tx) != length(ex))
    return(fail("TRANSCRIPT_CDS_MAP_UNRESOLVED"))
  cds_tx <- match(cds, ex)
  if (anyNA(cds_tx) || !identical(as.integer(cds_tx),
    seq.int(cds_tx[[1L]], length.out = length(cds))) ||
      !identical(substring(tx, cds_tx[[1L]], utils::tail(cds_tx, 1L)),
                 context$reference_cds)) return(fail("TRANSCRIPT_CDS_MAP_UNRESOLVED"))
  stop_genome <- translation_coordinates(model$stop_segments, model$strand)
  stop_tx <- match(stop_genome, ex)
  if (length(stop_tx) && (anyNA(stop_tx) ||
      !identical(as.integer(stop_tx), seq.int(stop_tx[[1L]], length.out = length(stop_tx)))))
    return(fail("ANNOTATED_STOP_MAP_UNRESOLVED"))
  stop_relation <- if (!length(stop_tx)) "ABSENT" else
    if (all(stop_tx %in% cds_tx)) "INSIDE_CDS" else
      if (all(stop_tx > utils::tail(cds_tx, 1L))) "AFTER_CDS" else "AMBIGUOUS"
  if (identical(stop_relation, "AMBIGUOUS")) return(fail("ANNOTATED_STOP_CDS_RELATION_UNRESOLVED"))
  if (length(stop_tx) && length(stop_tx) != 3L)
    return(fail("ANNOTATED_STOP_LENGTH_UNRESOLVED"))
  coding_end <- if (identical(stop_relation, "INSIDE_CDS")) stop_tx[[1L]] - 1L else
    utils::tail(cds_tx, 1L)
  if (coding_end < cds_tx[[1L]]) return(fail("CODING_EXTENT_UNRESOLVED"))
  start_features <- features[features$feature == "start_codon", , drop = FALSE]
  start_genome <- translation_coordinates(start_features[, c("start", "end"), drop = FALSE],
                                          model$strand)
  annotated_start <- length(start_genome) == 3L &&
    identical(as.integer(start_genome), as.integer(utils::head(cds, 3L)))
  raw <- features$raw_attributes
  incomplete_5 <- any(grepl("cds_start_NF", raw, fixed = TRUE)) ||
    as.integer(model$segments$phase[[if (model$strand == "-") nrow(model$segments) else 1L]]) != 0L ||
    (!annotated_start && !identical(context$contig, "chrM"))
  incomplete_3 <- any(grepl("cds_end_NF", raw, fixed = TRUE))
  if (length(stop_tx) && stop_tx[[1L]] != coding_end + 1L)
    return(fail("ANNOTATED_STOP_NOT_ADJACENT_TO_CDS"))
  sec_aa <- as.integer(bundle$sec$position[bundle$sec$transcript_id == context$transcript_id])
  sec_tx <- as.integer(cds_tx[[1L]] + (sec_aa - 1L) * 3L)
  if (length(sec_tx) && any(sec_tx + 2L > coding_end))
    return(fail("SEC_POSITION_OUTSIDE_CODING_EXTENT"))
  catalog <- protein_catalog[[context$protein_id]]
  if (is.null(catalog)) return(fail("REFERENCE_PROTEIN_MISSING"))
  catalog <- sub("\\*$", "", catalog)
  partial_terminal <- if (coding_end < nchar(tx)) substring(tx, coding_end + 1L) else ""
  # An inferred mitochondrial completion is a boundary rule, never added DNA.
  mito_partial <- identical(context$contig, "chrM") && !length(stop_tx) &&
    coding_end == nchar(tx) && nchar(context$reference_cds) %% 3L %in% 1:2 &&
    (endsWith(context$reference_cds, "T") || endsWith(context$reference_cds, "TA"))
  terminal_completion <- if (mito_partial) "CATALOG_VALIDATED_PARTIAL_TAA_PENDING" else "NONE"
  basic$status <- "READY"
  basic$reason <- NULL
  basic$reference_transcript <- tx
  basic$reference_cds <- context$reference_cds
  basic$reference_protein_catalog <- catalog
  basic$exon_segments <- model$exons
  basic$cds_segments <- model$segments
  basic$cds_phase <- model$segments$phase
  basic$start_feature <- start_features[, c("start", "end"), drop = FALSE]
  basic$stop_feature <- model$stop_segments
  basic$stop_relation <- stop_relation
  basic$annotated_start <- annotated_start
  basic$start_source <- if (annotated_start) "GTF_START_CODON" else
    if (identical(context$contig, "chrM")) "CDS_BOUNDARY_CATALOG_PENDING" else "UNANNOTATED"
  basic$start_codon <- substring(tx, cds_tx[[1L]], cds_tx[[1L]] + 2L)
  basic$coding_start_tx <- as.integer(cds_tx[[1L]])
  basic$coding_end_tx <- as.integer(coding_end)
  basic$annotated_stop_tx <- if (length(stop_tx)) as.integer(stop_tx[[1L]]) else NULL
  basic$annotated_stop_codon <- if (length(stop_tx))
    substring(tx, stop_tx[[1L]], stop_tx[[1L]] + 2L) else NULL
  basic$transcript_length <- as.integer(nchar(tx))
  basic$transcript_end_tx <- as.integer(nchar(tx))
  basic$exon_coordinates <- ex
  basic$cds_coordinates <- cds
  basic$partial_5 <- incomplete_5
  basic$partial_3 <- incomplete_3 || (!length(stop_tx) && !mito_partial)
  basic$partial_reason <- c(if (incomplete_5) "FIVE_PRIME_INCOMPLETE",
    if (incomplete_3) "THREE_PRIME_INCOMPLETE",
    if (!length(stop_tx) && !mito_partial) "STOP_FEATURE_ABSENT")
  basic$terminal_completion <- terminal_completion
  basic$sec_aa_positions <- sec_aa
  basic$sec_tx_positions <- sec_tx
  basic$translation_exception_ids <- paste0("GENCODE_SEC:", context$transcript_id,
                                            ":", sec_aa)
  basic$spec_id <- paste0("translation_spec:", digest::digest(charToRaw(canonical_bytes(list(
    policy = translation_policy(), transcript_id = context$transcript_id,
    protein_id = context$protein_id, reference_identity = bundle$resource_identity,
    code = basic$genetic_code_id, sec = as.list(sec_aa),
    terminal_completion = terminal_completion))), algo = "sha256", serialize = FALSE))
  basic
}
