translation_codons <- function(system) {
  code <- reconstruction_codons
  if (identical(system, "vertebrate_mitochondrial"))
    code[c("TGA", "ATA", "AGA", "AGG")] <- c("W", "M", "*", "*")
  code
}

# Project original-reference assertions in transcript orientation. This is a
# coordinate audit of T008 output; it is not a second haplotype caller.
translation_project_rows <- function(rows, coordinates, strand) {
  out <- list()
  for (row in rows) {
    e <- reconstruction_edit_core(row)
    if (is.null(e)) return(NULL)
    if (e$insertion) {
      first <- match(if (strand == "-") e$start - 1L else e$start, coordinates)
      if (is.na(first)) return(NULL)
      last <- first - 1L
    } else {
      first <- match(if (strand == "-") e$end else e$start, coordinates)
      last <- match(if (strand == "-") e$start else e$end, coordinates)
      if (is.na(first) || is.na(last) || last - first + 1L != nchar(e$old)) return(NULL)
      covered <- coordinates[seq.int(first, last)]
      expected <- if (strand == "-") seq.int(e$end, e$start) else
        seq.int(e$start, e$end)
      if (!identical(as.integer(covered), as.integer(expected))) return(NULL)
    }
    old <- if (strand == "-") reconstruction_reverse_complement(e$old) else e$old
    new <- if (strand == "-") reconstruction_reverse_complement(e$new) else e$new
    out[[length(out) + 1L]] <- list(event_id = e$event_id,
      ref_start = as.integer(first), ref_end = as.integer(last), old = old,
      new = new, delta = as.integer(e$delta), new_length = nchar(new))
  }
  if (!length(out)) return(out)
  out[order(vapply(out, `[[`, 0L, "ref_start"),
            vapply(out, `[[`, "", "event_id"), method = "radix")]
}

translation_atoms <- function(edits) {
  out <- list()
  for (e in edits) {
    if (nchar(e$old) == nchar(e$new) && nchar(e$old) > 1L) {
      old <- strsplit(e$old, "", fixed = TRUE)[[1L]]
      new <- strsplit(e$new, "", fixed = TRUE)[[1L]]
      for (j in which(old != new)) {
        a <- e; a$ref_start <- e$ref_start + j - 1L
        a$ref_end <- a$ref_start; a$old <- old[[j]]; a$new <- new[[j]]
        a$new_length <- 1L; out[[length(out) + 1L]] <- a
      }
    } else out[[length(out) + 1L]] <- e
  }
  if (!length(out)) return(out)
  out[order(vapply(out, `[[`, 0L, "ref_start"),
    vapply(out, function(e) as.integer(nchar(e$old) > 0L), 0L),
    vapply(out, `[[`, "", "event_id"), method = "radix")]
}

translation_apply_map <- function(reference, edits) {
  atoms <- translation_atoms(edits)
  keep <- list()
  for (e in atoms) {
    observed <- if (nchar(e$old)) substring(reference, e$ref_start, e$ref_end) else ""
    if (!identical(observed, e$old)) return(list(status = "REFERENCE_ALLELE_MISMATCH"))
    if (length(keep)) {
      p <- keep[[length(keep)]]
      if (identical(c(e$ref_start, e$ref_end, e$old, e$new),
                    c(p$ref_start, p$ref_end, p$old, p$new))) next
      if (e$ref_start <= p$ref_end || (e$ref_start == p$ref_start &&
          !nchar(e$old) && !nchar(p$old)))
        return(list(status = "OVERLAPPING_EDIT_CONFLICT"))
    }
    keep[[length(keep) + 1L]] <- e
  }
  pieces <- character(); origin <- integer(); cursor <- 1L
  for (e in keep) {
    if (e$ref_start > cursor) {
      pieces <- c(pieces, substring(reference, cursor, e$ref_start - 1L))
      origin <- c(origin, seq.int(cursor, e$ref_start - 1L))
    }
    pieces <- c(pieces, e$new)
    if (nchar(e$new)) origin <- c(origin, if (nchar(e$old) == nchar(e$new))
      seq.int(e$ref_start, length.out = nchar(e$new)) else
        rep(NA_integer_, nchar(e$new)))
    cursor <- e$ref_end + 1L
  }
  if (cursor <= nchar(reference)) {
    pieces <- c(pieces, substring(reference, cursor))
    origin <- c(origin, seq.int(cursor, nchar(reference)))
  }
  list(status = "COMPLETE", sequence = paste0(pieces, collapse = ""),
       origin = as.integer(origin), edits = keep)
}

translation_frame_state <- function(edits) {
  frame <- 0L; offset <- 0L; changes <- list(); intervals <- list(); open <- NULL
  for (e in edits) {
    before <- frame
    frame <- as.integer((frame + e$delta) %% 3L)
    if (e$delta %% 3L) {
      hap_start <- as.integer(e$ref_start + offset)
      change <- list(event_id = e$event_id, reference_nt = e$ref_start,
        haplotype_nt = hap_start, delta = e$delta,
        frame_before = before, frame_after = frame)
      changes[[length(changes) + 1L]] <- change
      if (before == 0L && frame != 0L) open <- list(
        start_event_id = e$event_id, reference_start_nt = e$ref_start,
        haplotype_start_nt = hap_start,
        start_aa = as.integer((hap_start - 1L) %/% 3L + 1L))
      if (before != 0L && frame == 0L && !is.null(open)) {
        open$restoration_event_id <- e$event_id
        open$reference_restoration_nt <- as.integer(e$ref_end + 1L)
        open$haplotype_restoration_nt <- as.integer(e$ref_end + 1L + offset + e$delta)
        open$first_original_frame_codon_nt <- as.integer(1L + 3L *
          ceiling((open$haplotype_restoration_nt - 1L) / 3L))
        open$restoration_aa <- as.integer((open$first_original_frame_codon_nt - 1L) %/% 3L + 1L)
        intervals[[length(intervals) + 1L]] <- open; open <- NULL
      }
    }
    offset <- offset + e$delta
  }
  if (!is.null(open)) intervals[[length(intervals) + 1L]] <- open
  list(changes = changes, shifted_intervals = intervals,
       final_frame = frame, net_length_delta = as.integer(offset))
}

translation_sec_used <- function(spec, sequence, origin, start) {
  used <- integer(); ids <- character()
  for (i in seq_along(spec$sec_tx_positions)) {
    original <- spec$sec_tx_positions[[i]] + 0:2
    mapped <- match(original, origin)
    if (anyNA(mapped) || !identical(as.integer(mapped), as.integer(mapped[[1L]] + 0:2)) ||
        (mapped[[1L]] - start) %% 3L != 0L ||
        !identical(substring(sequence, mapped[[1L]], mapped[[1L]] + 2L), "TGA")) next
    used <- c(used, mapped[[1L]])
    ids <- c(ids, spec$translation_exception_ids[[i]])
  }
  list(positions = as.integer(used), ids = ids)
}

translation_read <- function(sequence, origin, start, spec, end = nchar(sequence)) {
  if (is.null(sequence) || is.na(start) || start < 1L || start > nchar(sequence))
    return(list(status = "UNRESOLVED", boundary = "START_BOUNDARY_UNAVAILABLE",
      protein = NULL, stop_aa = NULL, stop_nt = NULL, residual_nt = NULL,
      sec_sites_used = character(), codon_origins = list()))
  if (end < start || end > nchar(sequence))
    return(list(status = "UNRESOLVED", boundary = "CODING_BOUNDARY_UNAVAILABLE",
      protein = NULL, stop_aa = NULL, stop_nt = NULL, residual_nt = NULL,
      sec_sites_used = character(), codon_origins = list()))
  count <- as.integer((end - start + 1L) %/% 3L)
  residual <- as.integer((end - start + 1L) %% 3L)
  if (!count) return(list(status = "BOUNDED", boundary = "INCOMPLETE_TERMINAL_CODON",
    protein = "", stop_aa = NULL, stop_nt = NULL, residual_nt = residual,
    sec_sites_used = character(), codon_origins = list(), stop_origin = NULL))
  begins <- start + (seq_len(count) - 1L) * 3L
  codons <- substring(sequence, begins, begins + 2L)
  if (any(!grepl("^[ACGT]{3}$", codons)))
    return(list(status = "UNRESOLVED", boundary = "AMBIGUOUS_CODON",
      protein = NULL, stop_aa = NULL, stop_nt = NULL, residual_nt = residual,
      sec_sites_used = character(), codon_origins = list()))
  sec <- translation_sec_used(spec, sequence, origin, start)
  aa <- unname(translation_codons(spec$translation_system)[codons])
  aa[begins %in% sec$positions & codons == "TGA"] <- "U"
  if (length(aa) && !identical(spec$start_source, "UNANNOTATED")) aa[[1L]] <- "M"
  stop <- which(aa == "*")
  first <- if (length(stop)) stop[[1L]] else NULL
  last <- if (is.null(first)) length(aa) else first - 1L
  protein <- if (last) paste0(aa[seq_len(last)], collapse = "") else ""
  completed <- FALSE
  if (is.null(first) && identical(spec$terminal_completion,
      "CATALOG_VALIDATED_PARTIAL_TAA") && residual %in% 1:2 && end == nchar(sequence)) {
    tail <- substring(sequence, end - residual + 1L)
    original_tail <- seq.int(spec$transcript_length - residual + 1L,
                             spec$transcript_length)
    completed <- identical(tail, if (residual == 1L) "T" else "TA") &&
      identical(as.integer(utils::tail(origin, residual)), as.integer(original_tail))
    if (completed) first <- as.integer(length(aa) + 1L)
  }
  n_codons <- if (is.null(first)) length(aa) else first - 1L
  list(status = if (!is.null(first)) "TERMINATED" else "BOUNDED",
    boundary = if (completed) "CATALOG_COMPLETED_PARTIAL_STOP" else
      if (!is.null(first)) "STOP_CODON" else
        if (residual) "INCOMPLETE_TERMINAL_CODON" else "TRANSCRIPT_END_NO_STOP",
    protein = protein, stop_aa = first,
    stop_nt = if (completed) as.integer(end - residual + 1L) else
      if (is.null(first)) NULL else as.integer(begins[[first]]),
    stop_origin = if (completed) as.integer(utils::tail(origin, residual)) else
      if (is.null(first)) NULL else
        as.integer(origin[seq.int(begins[[first]], begins[[first]] + 2L)]),
    residual_nt = residual, sec_sites_used = sec$ids,
    codon_origins = lapply(seq_len(n_codons), function(i)
      as.integer(origin[seq.int(begins[[i]], begins[[i]] + 2L)])))
}

translation_reference_oracle <- function(spec, strict = FALSE) {
  failure <- function(reason) {
    if (strict) stop(reason, call. = FALSE)
    list(status = "UNRESOLVED", reason = reason, spec = spec, read = NULL)
  }
  if (!identical(spec$status, "READY")) return(failure(spec$reason))
  if (isTRUE(spec$partial_5) || isTRUE(spec$partial_3))
    return(list(status = "PARTIAL", reason = paste(spec$partial_reason, collapse = ";"),
      spec = spec, read = NULL))
  starts <- if (identical(spec$translation_system, "vertebrate_mitochondrial"))
    c("ATG", "ATA", "ATT", "ATC", "GTG") else "ATG"
  if (!spec$annotated_start && !identical(spec$translation_system, "vertebrate_mitochondrial"))
    return(failure("REFERENCE_START_ANNOTATION_MISSING"))
  if (!spec$annotated_start && !spec$start_codon %in% starts)
    return(failure("MITOCHONDRIAL_START_BOUNDARY_UNRESOLVED"))
  if (identical(spec$terminal_completion, "CATALOG_VALIDATED_PARTIAL_TAA_PENDING"))
    spec$terminal_completion <- "CATALOG_VALIDATED_PARTIAL_TAA"
  if (identical(spec$start_source, "CDS_BOUNDARY_CATALOG_PENDING"))
    spec$start_source <- "CDS_BOUNDARY_CATALOG_VALIDATED"
  read <- translation_read(spec$reference_transcript,
    seq_len(spec$transcript_length), spec$coding_start_tx, spec)
  if (!identical(read$protein, spec$reference_protein_catalog))
    return(failure("REFERENCE_PROTEIN_MISMATCH"))
  if (!is.null(spec$annotated_stop_tx) &&
      !identical(read$stop_nt, spec$annotated_stop_tx))
    return(failure("REFERENCE_TERMINATION_MISMATCH"))
  if (is.null(spec$annotated_stop_tx) &&
      !identical(read$boundary, "CATALOG_COMPLETED_PARTIAL_STOP"))
    return(failure("REFERENCE_TERMINAL_BOUNDARY_UNRESOLVED"))
  list(status = "EXACT", reason = NULL, spec = spec, read = read)
}
