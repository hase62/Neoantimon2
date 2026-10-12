translation_side <- function(spec, sequence, origin, coding_edits) {
  mapped <- match(spec$coding_start_tx + 0:2, origin)
  intact_map <- !anyNA(mapped) && identical(as.integer(mapped),
                                            as.integer(mapped[[1L]] + 0:2))
  start <- if (intact_map) as.integer(mapped[[1L]]) else NA_integer_
  codon <- if (intact_map && start + 2L <= nchar(sequence))
    substring(sequence, start, start + 2L) else NULL
  frames <- translation_frame_state(coding_edits)
  if (!intact_map || !identical(codon, spec$start_codon))
    return(list(status = "START_LOST", start_codon = codon,
      start_effect = "lost", protein = NULL, stop_aa = NULL, stop_nt = NULL,
      boundary = "START_UNRESOLVED", residual_nt = NULL,
      sec_sites_used = character(), codon_origins = list(), frames = frames,
      start_tx = start))
  end <- if (isTRUE(spec$partial_3)) {
    mapped_end <- match(spec$coding_end_tx, origin)
    if (is.na(mapped_end)) NA_integer_ else as.integer(mapped_end)
  } else nchar(sequence)
  read <- translation_read(sequence, origin, start, spec, end)
  c(read, list(start_codon = codon, start_effect = "intact",
               frames = frames, start_tx = start))
}

translation_coding_edits <- function(rows, spec) {
  if (!length(rows)) return(list())
  all <- list()
  for (row in rows) {
    projected <- translation_project_rows(list(row), spec$cds_coordinates, spec$strand)
    if (is.null(projected)) next
    all <- c(all, projected)
  }
  if (!length(all)) return(all)
  all[order(vapply(all, `[[`, 0L, "ref_start"),
            vapply(all, `[[`, "", "event_id"), method = "radix")]
}

translation_frame_reach <- function(frame, read) {
  if (!length(frame$shifted_intervals)) return(frame)
  frame$shifted_intervals <- lapply(frame$shifted_intervals, function(x) {
    restore <- x$restoration_aa
    reached <- if (is.null(restore)) FALSE else
      !is.null(read$protein) && (is.null(read$stop_aa) || read$stop_aa >= restore)
    x$restoration_reached_by_translation <- reached
    x$altered_aa_start <- x$start_aa
    x$altered_aa_end <- if (is.null(read$protein)) NULL else
      if (reached) min(nchar(read$protein), restore - 1L) else nchar(read$protein)
    x$altered_aa_segment <- if (!is.null(x$altered_aa_end) &&
      x$altered_aa_end >= x$altered_aa_start)
      substring(read$protein, x$altered_aa_start, x$altered_aa_end) else ""
    x
  })
  frame
}

translation_frame_effect <- function(frame) {
  if (!length(frame$changes)) return("unchanged")
  restored <- length(frame$shifted_intervals) > 0L &&
    all(vapply(frame$shifted_intervals, function(x)
      isTRUE(x$restoration_reached_by_translation), TRUE))
  if (restored) "shifted_restored" else "shifted_unrestored"
}

# A codon is an anchor only when the same three original transcript bases
# encode the same amino acid on both haplotypes. This yields separate blocks
# for separated changes and lets a restored frame rejoin downstream codons.
translation_difference_intervals <- function(normal, mutant) {
  if (is.null(normal$protein) || is.null(mutant$protein)) return(NULL)
  a <- strsplit(normal$protein, "", fixed = TRUE)[[1L]]
  b <- strsplit(mutant$protein, "", fixed = TRUE)[[1L]]
  signature <- function(x) vapply(x$codon_origins, function(y)
    if (length(y) != 3L || anyNA(y)) "" else paste(y, collapse = ":"), "")
  sa <- signature(normal); sb <- signature(mutant)
  lookup <- stats::setNames(seq_along(sb)[nzchar(sb)], sb[nzchar(sb)])
  anchors <- list(); last_b <- 0L
  for (i in seq_along(sa)) {
    if (!nzchar(sa[[i]]) || !sa[[i]] %in% names(lookup)) next
    j <- as.integer(lookup[[sa[[i]]]])
    if (j > last_b && identical(a[[i]], b[[j]])) {
      anchors[[length(anchors) + 1L]] <- c(i, j); last_b <- j
    }
  }
  anchors <- c(list(c(0L, 0L)), anchors, list(c(length(a) + 1L, length(b) + 1L)))
  out <- list()
  for (k in seq_len(length(anchors) - 1L)) {
    p <- anchors[[k]]; q <- anchors[[k + 1L]]
    a_start <- p[[1L]] + 1L; a_end <- q[[1L]] - 1L
    b_start <- p[[2L]] + 1L; b_end <- q[[2L]] - 1L
    if (a_start > a_end && b_start > b_end) next
    old <- if (a_start <= a_end) paste0(a[seq.int(a_start, a_end)], collapse = "") else ""
    new <- if (b_start <= b_end) paste0(b[seq.int(b_start, b_end)], collapse = "") else ""
    # A changed codon can still encode the same amino acid. Trim equal ends
    # within each mapping-bounded block, without joining separate blocks.
    while (nzchar(old) && nzchar(new) && substring(old, 1L, 1L) ==
           substring(new, 1L, 1L)) {
      old <- substring(old, 2L); new <- substring(new, 2L)
      a_start <- a_start + 1L; b_start <- b_start + 1L
    }
    while (nzchar(old) && nzchar(new) && substring(old, nchar(old)) ==
           substring(new, nchar(new))) {
      old <- substring(old, 1L, nchar(old) - 1L)
      new <- substring(new, 1L, nchar(new) - 1L)
      a_end <- a_end - 1L; b_end <- b_end - 1L
    }
    if (identical(old, new)) next
    out[[length(out) + 1L]] <- list(normal_aa_start = as.integer(a_start),
      normal_aa_end = as.integer(a_end), mutant_aa_start = as.integer(b_start),
      mutant_aa_end = as.integer(b_end), normal_sequence = old,
      mutant_sequence = new, left_anchor = as.integer(p), right_anchor = as.integer(q))
  }
  out
}

translation_termination_effect <- function(normal, tumor) {
  if (is.null(tumor$protein) || is.null(normal$protein)) return("unresolved")
  n <- normal$stop_aa; t <- tumor$stop_aa
  if (!is.null(normal$stop_origin) && !is.null(tumor$stop_origin) &&
      identical(normal$stop_origin, tumor$stop_origin)) return("unchanged")
  if (!is.null(n) && !is.null(t) && t < n) return("premature")
  if (!is.null(n) && is.null(t)) return("lost_open_end")
  if (!is.null(n) && !is.null(t)) return("lost_with_recovered_stop")
  if (is.null(t) && tumor$boundary == "INCOMPLETE_TERMINAL_CODON")
    return("incomplete_terminal_codon")
  if (is.null(t)) return("transcript_boundary")
  "unchanged"
}

translation_protein_effect <- function(normal, tumor, frame_effect,
                                       termination_effect, intervals) {
  if (is.null(normal$protein) || is.null(tumor$protein)) return("unavailable")
  if (identical(normal$protein, tumor$protein)) return("identical")
  if (frame_effect == "shifted_unrestored") return("altered_tail")
  if (frame_effect == "shifted_restored") return("altered_segment")
  if (termination_effect == "premature") return("truncated")
  if (termination_effect %in% c("lost_with_recovered_stop", "lost_open_end"))
    return("extended")
  if (nchar(normal$protein) != nchar(tumor$protein)) return("inframe_indel")
  "substitution"
}

translation_side_from_context <- function(spec, context, rows, expected, role) {
  edits <- translation_project_rows(rows, spec$exon_coordinates, spec$strand)
  if (is.null(edits)) return(list(status = "UNRESOLVED", reason = "EXON_EDIT_PROJECTION_UNRESOLVED"))
  applied <- translation_apply_map(spec$reference_transcript, edits)
  if (applied$status != "COMPLETE")
    return(list(status = "UNRESOLVED", reason = applied$status))
  if (is.null(expected) || !identical(applied$sequence, expected))
    return(list(status = "UNRESOLVED", reason = paste0("T008_", role,
      "_TRANSCRIPT_DISAGREEMENT")))
  coding <- translation_coding_edits(rows, spec)
  side <- translation_side(spec, applied$sequence, applied$origin, coding)
  side$sequence <- applied$sequence
  side$origin <- applied$origin
  side$edits <- applied$edits
  side$frames <- translation_frame_reach(side$frames, side)
  side
}

translation_termination_geometry <- function(spec, tumor, termination_effect) {
  if (!identical(termination_effect, "premature"))
    return(list(status = "NOT_APPLICABLE", ptc_tx = NULL,
      final_junction_tx = NULL, distance_to_final_junction_nt = NULL,
      ptc_in_last_exon = NULL, intronless = NULL, nmd_prediction = NULL))
  exons <- spec$exon_segments
  if (is.null(exons) || !nrow(exons) || is.null(tumor$stop_nt))
    return(list(status = "JUNCTION_INFORMATION_INSUFFICIENT",
      ptc_tx = tumor$stop_nt, final_junction_tx = NULL,
      distance_to_final_junction_nt = NULL, ptc_in_last_exon = NULL,
      intronless = NULL, nmd_prediction = NULL))
  order <- if (spec$strand == "-") rev(seq_len(nrow(exons))) else seq_len(nrow(exons))
  lengths <- as.integer(exons$end[order] - exons$start[order] + 1L)
  starts <- cumsum(c(1L, utils::head(lengths, -1L)))
  ends <- cumsum(lengths)
  for (e in tumor$edits) {
    hit <- which(starts <= e$ref_start & ends >= max(e$ref_start, e$ref_end))
    if (length(hit) != 1L)
      return(list(status = "JUNCTION_INFORMATION_INSUFFICIENT",
        ptc_tx = tumor$stop_nt, final_junction_tx = NULL,
        distance_to_final_junction_nt = NULL, ptc_in_last_exon = NULL,
        intronless = NULL, nmd_prediction = NULL))
    lengths[[hit]] <- lengths[[hit]] + e$delta
  }
  intronless <- length(lengths) == 1L
  junction <- if (intronless) NULL else as.integer(sum(utils::head(lengths, -1L)))
  list(status = "STRUCTURAL_ONLY_NO_BINARY_POLICY", ptc_tx = tumor$stop_nt,
    final_junction_tx = junction,
    distance_to_final_junction_nt = if (is.null(junction)) NULL else
      as.integer(junction - tumor$stop_nt),
    ptc_in_last_exon = if (intronless) TRUE else tumor$stop_nt > junction,
    intronless = intronless, nmd_prediction = NULL)
}

translation_extension_events <- function(context, rows, spec, tumor, termination_effect) {
  if (!termination_effect %in% c("lost_with_recovered_stop", "lost_open_end") ||
      is.null(spec$annotated_stop_tx))
    return(list(status = "NONE", event_ids = character()))
  included <- c(context$germline_event_ids, context$somatic_event_ids)
  first <- spec$annotated_stop_tx + 3L
  limit <- if (is.null(tumor$stop_nt)) spec$transcript_length else
    min(spec$transcript_length, tumor$stop_nt + 2L)
  if (first > limit) return(list(status = "NONE", event_ids = character()))
  touched <- spec$exon_coordinates[seq.int(first, limit)]
  found <- list()
  for (row in rows) {
    if (row$event_id %in% included || !identical(row$sample_id, context$sample_id) ||
        !identical(row$contig, context$contig) ||
        !row$event_role %in% c("somatic", "germline")) next
    core <- reconstruction_edit_core(row)
    if (is.null(core)) next
    e <- translation_project_rows(list(row), spec$exon_coordinates, spec$strand)
    affects <- if (!is.null(e) && length(e))
      e[[1L]]$ref_start <= limit && e[[1L]]$ref_end >= first - 1L else
        any(if (core$insertion) touched %in% c(core$start - 1L, core$start) else
          touched >= core$start & touched <= core$end)
    if (!affects) next
    relations <- vapply(context$somatic_event_ids, function(id)
      phase_relative_events(rows, id, row$event_id)$relation, "")
    if (length(relations) && all(relations == "OPPOSITE_HAPLOTYPE")) next
    found[[row$event_id]] <- if (is.null(e)) "PROJECTION_UNRESOLVED" else
      if (reconstruction_homozygous(row) ||
          (length(relations) && all(relations == "SAME_HAPLOTYPE")))
        "KNOWN_SAME_HAPLOTYPE" else "PHASE_UNRESOLVED"
  }
  if (!length(found)) return(list(status = "NONE", event_ids = character()))
  ids <- sort(names(found), method = "radix")
  kinds <- unlist(found[ids], use.names = FALSE)
  status <- if ("PROJECTION_UNRESOLVED" %in% kinds)
    "EXTENSION_PROJECTION_UNRESOLVED" else if ("PHASE_UNRESOLVED" %in% kinds)
      "EXTENSION_PHASE_UNRESOLVED" else "EXTENSION_EVENT_OUTSIDE_T008_HAPLOTYPE"
  list(status = status, event_ids = ids, event_statuses = stats::setNames(kinds, ids))
}

translation_context <- function(context, model, rows, bundle, proteins,
                                transcripts, features) {
  result <- context
  fallback <- list(status = "UNRESOLVED", reason = context$reason,
    reconstruction_id = context$reconstruction_id, translation_policy = translation_policy(),
    start_effect = "unresolved", frame_effect = "unresolved",
    termination_effect = "unresolved", protein_effect = "unresolved",
    translation_system = if (context$contig == "chrM")
      "vertebrate_mitochondrial" else "nuclear_standard")
  result$translation <- fallback
  spec <- translation_specification(model, context, bundle, proteins, transcripts, features)
  result$translation$spec <- spec
  if (spec$status != "READY") {
    result$translation$reason <- spec$reason
    return(result)
  }
  oracle <- translation_reference_oracle(spec)
  result$translation$reference_oracle <- oracle$status
  if (oracle$status == "UNRESOLVED") {
    result$translation$reason <- oracle$reason
    return(result)
  }
  spec <- oracle$spec
  result$translation$spec <- spec
  ids <- unique(c(context$germline_event_ids, context$somatic_event_ids))
  by_id <- stats::setNames(rows, vapply(rows, `[[`, "", "event_id"))
  if (any(!ids %in% names(by_id))) {
    result$translation$reason <- "CONTRIBUTING_EVENT_MISSING"
    return(result)
  }
  normal_rows <- unname(by_id[context$germline_event_ids])
  somatic_rows <- unname(by_id[context$somatic_event_ids])
  combined <- reconstruction_tumor_rows(normal_rows, somatic_rows)
  if (combined$status != "COMPLETE") {
    result$translation$reason <- combined$status
    return(result)
  }
  normal <- translation_side_from_context(spec, context, normal_rows,
    context$patient_normal_transcript, "NORMAL")
  tumor <- translation_side_from_context(spec, context, combined$rows,
    context$tumor_transcript, "TUMOR")
  if (identical(normal$status, "UNRESOLVED") || identical(tumor$status, "UNRESOLVED")) {
    result$translation$reason <- if (identical(normal$status, "UNRESOLVED"))
      normal$reason else tumor$reason
    return(result)
  }
  if (identical(context$status, "COMPLETE") &&
      (!identical(normal$protein, context$patient_normal_protein) ||
       !identical(tumor$protein, context$tumor_protein))) {
    result$translation$reason <- "T008_COMPLETE_PROTEIN_DISAGREEMENT"
    return(result)
  }
  frame_effect <- translation_frame_effect(tumor$frames)
  termination_effect <- translation_termination_effect(normal, tumor)
  outside <- translation_extension_events(context, rows, spec, tumor, termination_effect)
  if (outside$status != "NONE") {
    result$translation$status <- "UNRESOLVED"
    result$translation$reason <- outside$status
    result$translation$extension_event_ids <- outside$event_ids
    result$translation$extension_event_statuses <- outside$event_statuses
    result$translation$spec <- spec
    return(result)
  }
  intervals <- translation_difference_intervals(normal, tumor)
  reference_side <- if (oracle$status == "EXACT") oracle$read else NULL
  reference_intervals <- if (!is.null(reference_side))
    translation_difference_intervals(reference_side, tumor) else NULL
  protein_effect <- translation_protein_effect(normal, tumor, frame_effect,
    termination_effect, intervals)
  translation_status <- if (isTRUE(spec$partial_5)) "PARTIAL" else
    if (is.null(tumor$protein)) "UNRESOLVED" else
      if (isTRUE(spec$partial_3)) "PARTIAL" else
        if (is.null(tumor$stop_aa)) "BOUNDED" else "EVALUABLE"
  if (isTRUE(spec$partial_5)) {
    normal$protein <- NULL; tumor$protein <- NULL
    protein_effect <- "unavailable"; intervals <- NULL; reference_intervals <- NULL
  }
  first <- if (length(intervals)) intervals[[1L]]$mutant_aa_start else NULL
  last <- if (length(intervals)) utils::tail(intervals, 1L)[[1L]]$mutant_aa_end else NULL
  frameshift_start <- if (length(tumor$frames$shifted_intervals))
    tumor$frames$shifted_intervals[[1L]]$start_aa else NULL
  restored <- Filter(function(x) isTRUE(x$restoration_reached_by_translation),
                     tumor$frames$shifted_intervals)
  extension <- if (termination_effect %in% c("lost_with_recovered_stop", "lost_open_end") &&
    !is.null(normal$stop_aa) && !is.null(tumor$protein) &&
    nchar(tumor$protein) >= normal$stop_aa)
    substring(tumor$protein, normal$stop_aa) else NULL
  result$translation <- list(status = translation_status,
    reason = if (translation_status == "PARTIAL") paste(spec$partial_reason, collapse = ";") else
      if (translation_status == "UNRESOLVED") "START_BOUNDARY_UNRESOLVED" else NULL,
    reconstruction_id = context$reconstruction_id, translation_policy = translation_policy(),
    translation_id = paste0("translation:", digest::digest(charToRaw(canonical_bytes(list(
      reconstruction_id = context$reconstruction_id, policy = translation_policy(),
      spec_id = spec$spec_id))), algo = "sha256", serialize = FALSE)),
    spec = spec, reference_oracle = oracle$status,
    analysis_id = context$analysis_id, sample_id = context$sample_id,
    transcript_id = context$transcript_id, protein_id = context$protein_id,
    germline_event_ids = context$germline_event_ids,
    somatic_event_ids = context$somatic_event_ids,
    phase_state = context$phase_state, phase_context = context$phase_context,
    reference_identity = context$reference_resources,
    reconstruction_policy = context$reconstruction_policy,
    translation_system = spec$translation_system, genetic_code_id = spec$genetic_code_id,
    start_effect = if (isTRUE(spec$partial_5)) "unresolved" else
      if (identical(tumor$start_effect, "lost")) "lost" else "intact",
    normal_start_effect = if (isTRUE(spec$partial_5)) "unresolved" else normal$start_effect,
    frame_effect = frame_effect, normal_frame_effect = translation_frame_effect(normal$frames),
    termination_effect = termination_effect, protein_effect = protein_effect,
    reference_protein = if (oracle$status == "EXACT") oracle$read$protein else NULL,
    patient_normal_protein = normal$protein, mutant_protein = tumor$protein,
    patient_normal_boundary = normal$boundary, mutant_boundary = tumor$boundary,
    translation_boundary_reached = tumor$boundary,
    reference_stop_aa = if (oracle$status == "EXACT") oracle$read$stop_aa else NULL,
    normal_stop_aa = normal$stop_aa, mutant_stop_aa = tumor$stop_aa,
    reference_stop_origin = if (oracle$status == "EXACT")
      oracle$read$stop_origin else NULL,
    normal_stop_origin = normal$stop_origin,
    mutant_stop_origin = tumor$stop_origin,
    first_changed_aa = first, last_changed_aa = last,
    frameshift_start_aa = frameshift_start,
    frame_restore_aa = if (length(restored)) restored[[1L]]$restoration_aa else NULL,
    frame_restore_nt = if (length(restored)) restored[[1L]]$haplotype_restoration_nt else NULL,
    normal_tumor_difference_intervals = intervals,
    reference_tumor_difference_intervals = reference_intervals,
    novel_tail_sequence = if (frame_effect == "shifted_unrestored" &&
      !is.null(tumor$protein) && !is.null(first)) substring(tumor$protein, first) else NULL,
    extension_sequence = extension,
    sec_sites_used = tumor$sec_sites_used,
    translation_exception_ids = tumor$sec_sites_used,
    normal_frame = normal$frames, mutant_frame = tumor$frames,
    normal_stop_nt = normal$stop_nt, mutant_stop_nt = tumor$stop_nt,
    residual_terminal_nt = tumor$residual_nt,
    termination_geometry = translation_termination_geometry(spec, tumor,
      termination_effect))
  result
}

translation_contexts <- function(project, bundle, rows, contexts) {
  if (!length(contexts)) return(list())
  ids <- sort(unique(vapply(contexts, `[[`, "", "transcript_id")), method = "radix")
  models <- reconstruction_models(bundle, ids)
  features <- reference_transcript_features(bundle, ids)
  protein_ids <- unique(stats::na.omit(vapply(models, `[[`, "", "protein_id")))
  proteins <- reconstruction_fasta_select(bundle$resources$protein_fasta$resolved_path,
                                          protein_ids)
  transcripts <- reconstruction_fasta_select(bundle$resources$transcript_fasta$resolved_path,
                                             ids)
  lapply(contexts, function(x) translation_context(x, models[[x$transcript_id]], rows,
    bundle, proteins, transcripts,
    features[features$transcript_id == x$transcript_id, , drop = FALSE]))
}
