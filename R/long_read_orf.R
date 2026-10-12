long_read_orf_policy <- function() "noncanonical-orf-v1"
long_read_start_context <- function(sequence, start) list(
  context = substring(sequence, max(1L, start - 6L), min(nchar(sequence), start + 5L)),
  minus_3 = if (start > 3L) substring(sequence, start - 3L, start - 3L) else NULL,
  plus_4 = if (start + 3L <= nchar(sequence)) substring(sequence, start + 3L, start + 3L) else NULL,
  hard_filter = FALSE)

long_read_provider_orfs <- function(observation, rec) {
  out <- list()
  for (meta in observation$metadata) {
    r <- meta$raw_fields
    if (is.null(r$CDS_start) || r$CDS_start %in% c("", "NA", ".")) next
    out[[length(out) + 1L]] <- list(start = long_read_integer_evidence(r$CDS_start),
      end = if (is.null(r$CDS_end) || r$CDS_end %in% c("", "NA", ".")) NULL else
        long_read_integer_evidence(r$CDS_end),
      frame = if (is.null(r$frame)) NULL else long_read_integer_evidence(r$frame),
      strand = r$strand, source = list(metadata_input = meta$input_id, raw_fields = r))
  }
  cds <- observation$model$cds_features
  if (length(cds)) {
    cds <- cds[order(vapply(cds, `[[`, 0L, "start"))]
    if (rec$strand == "-") cds <- rev(cds)
    coords <- as.integer(unlist(lapply(cds, function(e) if (rec$strand == "+")
      seq.int(e$start, e$end) else seq.int(e$end, e$start))))
    mapped <- match(coords, rec$coordinates)
    contiguous <- !anyNA(mapped) && identical(mapped, seq.int(min(mapped), max(mapped)))
    sizes <- vapply(cds, function(e) e$end - e$start + 1L, 0L)
    phases <- suppressWarnings(as.integer(vapply(cds, `[[`, "", "phase")))
    expected <- (3L - c(0L, utils::head(cumsum(sizes), -1L)) %% 3L) %% 3L
    out[[length(out) + 1L]] <- list(start = if (contiguous) min(mapped) else NA_integer_,
      end = if (contiguous) max(mapped) else NA_integer_,
      invalid_mapping = !contiguous || anyNA(phases) || any(phases != expected),
      strand = rec$strand, source = list(gtf_cds = cds))
  }
  out
}
long_read_provider_orf_validate <- function(p, rec) {
  start <- p$start; end <- p$end
  reason <- if (isTRUE(p$invalid_mapping)) "CDS_MAPPING_OR_PHASE_INCONSISTENT" else
    if (length(start) != 1L || is.na(start) || start < 1L || start + 2L > nchar(rec$sequence))
      "START_OUT_OF_BOUNDS" else if (!is.null(p$strand) && p$strand != rec$strand)
      "STRAND_INCONSISTENT" else if (substring(rec$sequence, start, start + 2L) != "ATG")
      "START_CODON_INCONSISTENT" else if (!is.null(p$frame) &&
        (length(p$frame) != 1L || is.na(p$frame) || p$frame != (start - 1L) %% 3L))
      "FRAME_INCONSISTENT" else NULL
  if (is.null(reason) && !is.null(end)) {
    if (length(end) != 1L || is.na(end) || end < start + 2L || end > nchar(rec$sequence) ||
        (end - start + 1L) %% 3L != 0L) reason <- "CDS_END_INCONSISTENT" else {
      protein <- reconstruction_translate(substring(rec$sequence, start, end))
      stops <- which(strsplit(protein, "", fixed = TRUE)[[1L]] == "*")
      if (length(stops) && !identical(stops, as.integer(nchar(protein))))
        reason <- "PREMATURE_STOP_IN_CDS"
      if (!length(stops) && end + 3L <= nchar(rec$sequence) &&
          !substring(rec$sequence, end + 1L, end + 3L) %in% c("TAA", "TAG", "TGA"))
        reason <- "SUPPLIED_END_NOT_TERMINAL"
    }
  }
  list(status = if (is.null(reason)) "valid" else "inconsistent", reason = reason)
}

long_read_translate_orf <- function(structure_id, rec, comparison, start, type,
                                    policy, reference_ids = character(), evidence = list()) {
  spec <- list(sec_tx_positions = integer(), translation_exception_ids = character(),
    translation_system = "nuclear_standard", start_source = type, terminal_completion = "NONE")
  # Recoding is retained only for an exact authoritative transcript and its
  # annotated start. Novel context never acquires Sec just by genomic overlap.
  if (type == "reference_anchored_orf" && comparison$reference_equivalent) {
    for (cmp in comparison$comparisons) {
      r <- cmp$reference
      if (!cmp$exact_structure || !r$transcript_id %in% reference_ids ||
          is.null(r$coding_start) || start != r$coding_start) next
      for (aa in r$sec_aa) {
        pos <- as.integer(start + (aa - 1L) * 3L)
        spec$sec_tx_positions <- c(spec$sec_tx_positions, pos)
        spec$translation_exception_ids <- c(spec$translation_exception_ids,
          paste0("GENCODE_SEC:", r$transcript_id, ":", aa))
      }
    }
  }
  read <- translation_read(rec$sequence, seq_len(nchar(rec$sequence)), start, spec)
  id <- long_read_id("orf_hypothesis_id", list(isoform_structure_id = structure_id,
    start = start, frame = as.integer((start - 1L) %% 3L), policy = policy,
    evidence_type = type, reference_ids = as.list(sort(reference_ids))))
  if (is.null(read$protein) || !nzchar(read$protein)) return(list(
    orf_hypothesis_id = id, evidence_type = type, status = "UNAVAILABLE", reason = read$boundary))
  starts <- as.integer(start + (seq_len(nchar(read$protein)) - 1L) * 3L)
  inserted <- integer(); shifted <- integer(); novelty_paths <- list()
  if (!comparison$reference_equivalent && comparison$status == "COMPARED") {
    for (cmp in comparison$comparisons) {
      r <- cmp$reference
      # A reference-anchored origin uses only its own authoritative reference.
      if (type == "reference_anchored_orf" && !r$transcript_id %in% reference_ids) next
      aa_inserted <- integer()
      for (nt in cmp$novel_nt_intervals) aa_inserted <- union(aa_inserted,
        which(starts <= nt$end_nt & starts + 2L >= nt$start_nt))
      reference_start <- match(rec$coordinates[[start]], r$reconstructed$coordinates)
      same_start <- !is.na(reference_start) && reference_start + 2L <=
        length(r$reconstructed$coordinates) &&
        identical(rec$coordinates[start + 0:2], r$reconstructed$coordinates[reference_start + 0:2])
      aa_shifted <- integer()
      if (same_start) {
        # Per-codon phase comparison permits restoration after another junction;
        # unchanged in-frame downstream codons never become globally novel.
        for (i in seq_along(starts)) {
          mapped <- cmp$reference_coordinate_map[starts[[i]] + 0:2]
          if (!anyNA(mapped) && all(diff(mapped) == 1L) &&
              (mapped[[1L]] - reference_start) %% 3L != 0L) aa_shifted <- c(aa_shifted, i)
        }
      }
      inserted <- union(inserted, aa_inserted); shifted <- union(shifted, aa_shifted)
      novelty_paths[[length(novelty_paths) + 1L]] <- list(transcript_id = r$transcript_id,
        novel_nt_intervals = cmp$novel_nt_intervals, inserted_sequence_aa = as.list(aa_inserted),
        novel_frame_aa = as.list(aa_shifted))
    }
  }
  inserted_boundaries <- unlist(lapply(novelty_paths, function(path) unlist(lapply(
    path$novel_nt_intervals, function(nt) c(nt$start_nt - 1L, nt$end_nt)))))
  cuts <- sort(unique(c(comparison$novel_junction_nt, inserted_boundaries)))
  split <- integer(); junction_aa <- integer()
  for (j in cuts) {
    if (j < start || j >= start + nchar(read$protein) * 3L) next
    cut <- as.integer((j - start + 1L) %/% 3L)
    if (cut >= 1L && cut < nchar(read$protein)) junction_aa <- c(junction_aa, cut)
    if ((j - start + 1L) %% 3L != 0L) split <- c(split, cut + 1L)
  }
  novel <- sort(unique(as.integer(c(inserted, shifted, split))))
  list(status = "TRANSLATED", orf_hypothesis_id = id,
    protein_origin_id = long_read_id("protein_origin_id", list(orf = id, protein = read$protein)),
    evidence_type = type, evidence = evidence, reference_transcript_ids = as.list(sort(reference_ids)),
    policy = policy, start_nt = as.integer(start), start_codon = "ATG", start_type = "AUG",
    frame = as.integer((start - 1L) %% 3L), protein = read$protein, orf_length_aa = nchar(read$protein),
    termination = read$boundary, stop_nt = read$stop_nt, residual_nt = read$residual_nt,
    stop_status = if (is.null(read$stop_nt)) "no_stop" else "first_in_frame_stop",
    boundary_status = if (is.null(read$stop_nt)) "transcript_boundary" else "terminated",
    reference_orf_relationships = lapply(comparison$comparisons, function(cmp) {
      r <- cmp$reference
      anchor <- if (is.null(r$coding_start)) NA_integer_ else
        match(r$coding_start_coordinates[[1L]], rec$coordinates)
      list(transcript_id = r$transcript_id, relationship = if (is.na(anchor)) "unresolved" else
        if (start < anchor) "uORF" else if (start == anchor) "annotated_start" else
          if ((start - anchor) %% 3L == 0L) "downstream_ORF" else "overlapping_alternative_frame_ORF")
    }),
    start_context = long_read_start_context(rec$sequence, start),
    translation_exception_ids = as.list(read$sec_sites_used),
    frame_state = if (length(shifted)) "out_of_frame" else "in_frame",
    novel_aa = as.list(novel), inserted_sequence_aa = as.list(sort(as.integer(inserted))),
    novel_frame_aa = as.list(sort(as.integer(shifted))), junction_aa = as.list(sort(unique(junction_aa))),
    novelty_paths = novelty_paths, peptide_eligible = length(novel) > 0L || length(junction_aa) > 0L,
    reference_comparison = if (type == "de_novo_aug_orf") "noncanonical_orf_no_single_reference" else
      if (comparison$reference_equivalent) "reference_transcript_counterpart" else "source_transcript_mapping")
}

long_read_orfs <- function(structure_id, rec, comparison, observations) {
  out <- list(); unavailable <- list()
  if (comparison$status %in% c("FUSION_MAPPING_DEFERRED") ||
      observations[[1L]]$model$structure$contig %in% c("chrM", "MT"))
    return(list(hypotheses = list(), unavailable = list(list(status = "UNSUPPORTED_DEFERRED"))))
  for (cmp in comparison$comparisons) {
    r <- cmp$reference; start <- r$coding_start_coordinates
    p <- if (length(start)) match(start[[1L]], rec$coordinates) else NA_integer_
    retained <- length(start) == 3L && !is.na(p) && p + 2L <= length(rec$coordinates) &&
      identical(start, rec$coordinates[p + 0:2])
    if (!retained) unavailable[[length(unavailable) + 1L]] <- list(evidence_type = "reference_anchored_orf",
      transcript_id = r$transcript_id, status = "UNAVAILABLE", reason = "ANNOTATED_START_NOT_RETAINED") else
      out[[length(out) + 1L]] <- long_read_translate_orf(structure_id, rec, comparison,
        p, "reference_anchored_orf", "long-read-reference-orf-v1", r$transcript_id)
  }
  for (observation in observations) for (p in long_read_provider_orfs(observation, rec)) {
    validated <- long_read_provider_orf_validate(p, rec)
    if (validated$status != "valid") unavailable[[length(unavailable) + 1L]] <- list(
      evidence_type = "provider_supplied_orf", provider_orf_status = validated$status,
      reason = validated$reason, observation_id = observation$long_read_observation_id, evidence = p) else
      out[[length(out) + 1L]] <- long_read_translate_orf(structure_id, rec, comparison,
        p$start, "provider_supplied_orf", "long-read-provider-orf-v1",
        evidence = list(observation_id = observation$long_read_observation_id, coordinates = p))
  }
  if (any(vapply(observations, `[[`, TRUE, "de_novo_aug"))) {
    starts <- gregexpr("ATG", rec$sequence, fixed = TRUE)[[1L]]
    for (p in starts[starts > 0L]) out[[length(out) + 1L]] <- long_read_translate_orf(
      structure_id, rec, comparison, p, "de_novo_aug_orf", long_read_orf_policy())
  }
  if (length(out)) {
    groups <- split(out, vapply(out, `[[`, "", "orf_hypothesis_id"))
    out <- lapply(groups, function(g) {
      x <- g[[1L]]; x$evidence_records <- lapply(g, `[[`, "evidence"); x
    })
  }
  list(hypotheses = unname(out), unavailable = unavailable)
}
