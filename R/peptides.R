peptide_policy <- function(x = NULL) {
  defaults <- list(class_i_lengths = as.list(8:11), class_ii_lengths = as.list(13:25))
  if (!is.null(x)) {
    strict_keys(x, names(defaults), "peptides")
    for (name in names(x)) defaults[name] <- x[name]
  }
  for (name in names(defaults)) {
    z <- unlist(defaults[[name]], use.names = FALSE)
    if (!is.numeric(z) || !length(z) || anyNA(z) || any(!is.finite(z)) ||
        any(z < 1 | z > .Machine$integer.max | z != floor(z))) stop("Invalid peptide lengths")
    defaults[[name]] <- as.list(as.integer(sort(unique(z), method = "radix")))
  }
  defaults
}

peptide_digest <- function(sequence) digest::digest(charToRaw(sequence),
  algo = "sha256", serialize = FALSE)

peptide_map <- function(normal, mutant, intervals) {
  map <- rep(NA_integer_, nchar(mutant)); a <- 1L; b <- 1L
  text <- function(s, first, last) if (first > last) "" else substring(s, first, last)
  for (x in intervals) {
    ns <- x$normal_aa_start; ne <- x$normal_aa_end
    ms <- x$mutant_aa_start; me <- x$mutant_aa_end
    if (ns < a || ms < b || ne < ns - 1L || me < ms - 1L ||
        ne > nchar(normal) || me > nchar(mutant) ||
        !identical(text(normal, ns, ne), x$normal_sequence) ||
        !identical(text(mutant, ms, me), x$mutant_sequence) ||
        !identical(text(normal, a, ns - 1L), text(mutant, b, ms - 1L)))
      stop("Invalid T009 protein difference mapping")
    if (ms > b) map[seq.int(b, ms - 1L)] <- seq.int(a, ns - 1L)
    if (ne - ns == me - ms && ms <= me) map[seq.int(ms, me)] <- seq.int(ns, ne)
    a <- ne + 1L; b <- me + 1L
  }
  if (!identical(text(normal, a, nchar(normal)), text(mutant, b, nchar(mutant))))
    stop("Invalid T009 unchanged protein suffix")
  if (b <= nchar(mutant)) map[seq.int(b, nchar(mutant))] <- seq.int(a, nchar(normal))
  map
}

peptide_hits <- function(intervals, start, end) which(vapply(intervals, function(x) {
  if (x$mutant_aa_start <= x$mutant_aa_end)
    start <= x$mutant_aa_end && end >= x$mutant_aa_start
  else start < x$mutant_aa_start && end >= x$mutant_aa_start
}, TRUE))

peptide_correspondence <- function(normal, map, intervals, start, end,
                                   frame_effect, mutant_sequence, frameshift_start = NULL) {
  missing <- list(status = "NO_DIRECT_COUNTERPART", sequence = NULL,
    reference_start_aa = NULL, reference_end_aa = NULL, reference_span = NULL,
    sequence_comparison = "UNRESOLVED")
  if (is.null(normal) || is.null(map)) { missing$status <- "UNAVAILABLE"; return(missing) }
  if (identical(frame_effect, "shifted_unrestored") &&
      (is.null(frameshift_start) || end >= frameshift_start)) return(missing)
  left <- map[[start]]; right <- map[[end]]
  if (is.na(left) || is.na(right) || left > right) return(missing)
  span <- substring(normal, left, right)
  direct <- right - left == end - start
  # Equal-length endpoints are insufficient when multiple indels cancel.
  positions <- map[seq.int(start, end)]
  continuous <- direct && !anyNA(positions) &&
    identical(as.integer(positions), seq.int(left, right))
  list(status = if (!continuous) "JUNCTION_INTERVAL" else
      if (left == start) "COORDINATE_ALIGNED" else "INDEL_SHIFTED",
    sequence = if (continuous) span else NULL,
    reference_start_aa = left, reference_end_aa = right, reference_span = span,
    sequence_comparison = if (!continuous) "NO_SAME_LENGTH_COUNTERPART" else
      if (identical(span, mutant_sequence)) "IDENTICAL" else "DIFFERENT")
}

peptide_enumerate <- function(context, policy = peptide_policy()) {
  t <- context$translation
  provenance <- peptide_translation_provenance(context)
  answer <- function(status, peptides = list(), reason = NULL) list(
    status = status, reason = reason, provenance = provenance, peptides = peptides)
  if (is.null(t$mutant_protein) || is.null(t$patient_normal_protein) ||
      !t$status %in% c("EVALUABLE", "BOUNDED", "PARTIAL"))
    return(answer("UNAVAILABLE", reason = if (is.null(t$reason))
      "MUTANT_OR_NORMAL_PROTEIN_UNAVAILABLE" else t$reason))
  normal <- t$patient_normal_protein; mutant <- t$mutant_protein
  intervals <- t$normal_tumor_difference_intervals
  if (is.null(intervals)) return(answer("UNAVAILABLE", reason = "DIFFERENCE_MAPPING_UNAVAILABLE"))
  if (identical(normal, mutant)) return(answer("PROTEIN_IDENTICAL"))
  map <- peptide_map(normal, mutant, intervals)
  reference_map <- if (!is.null(t$reference_protein) &&
    !is.null(t$reference_tumor_difference_intervals))
    peptide_map(t$reference_protein, mutant, t$reference_tumor_difference_intervals) else NULL
  lengths <- sort(unique(unlist(policy)), method = "radix")
  peptides <- list(); n <- nchar(mutant)
  for (len in lengths) {
    if (len > n) next
    starts <- sort(unique(unlist(lapply(intervals, function(x) {
      lo <- max(1L, x$mutant_aa_start - len + 1L)
      hi <- min(n - len + 1L, if (x$mutant_aa_end < x$mutant_aa_start)
        x$mutant_aa_start - 1L else x$mutant_aa_end)
      if (lo > hi) integer() else seq.int(lo, hi)
    }))), method = "radix")
    for (start in starts) {
      end <- start + len - 1L
      hits <- peptide_hits(intervals, start, end)
      if (!length(hits)) next
      sequence <- substring(mutant, start, end)
      peptide_sequence_validate(sequence)
      if (grepl("U", sequence, fixed = TRUE) && !length(t$translation_exception_ids))
        stop("T009 Sec translation provenance missing")
      changed <- unique(unlist(lapply(intervals[hits], function(x) {
        a <- max(start, x$mutant_aa_start); b <- min(end, x$mutant_aa_end)
        if (a > b) integer() else seq.int(a, b)
      })))
      junctions <- vapply(Filter(function(x) x$mutant_aa_start > x$mutant_aa_end,
        intervals[hits]), `[[`, 0L, "mutant_aa_start")
      peptide <- list(peptide_sequence = sequence, peptide_start_aa = as.integer(start),
        peptide_end_aa = as.integer(end), peptide_length = as.integer(len),
        classes = as.list(c(if (len %in% unlist(policy$class_i_lengths)) "I",
          if (len %in% unlist(policy$class_ii_lengths)) "II")),
        mutant_reference_status = "MUTANT", protein_checksum = peptide_digest(mutant),
        difference_interval_indices = as.list(as.integer(hits)),
        altered_intervals = intervals[hits], altered_positions = as.list(as.integer(changed)),
        junction_positions = as.list(as.integer(junctions)),
        overlaps_altered_sequence = length(changed) > 0L,
        overlaps_altered_junction = length(junctions) > 0L,
        altered_at_n_terminus = start %in% changed,
        altered_at_c_terminus = end %in% changed,
        altered_internal = any(changed > start & changed < end),
        normal_correspondence = peptide_correspondence(normal, map, intervals,
          start, end, t$frame_effect, sequence, t$frameshift_start_aa),
        reference_correspondence = peptide_correspondence(t$reference_protein,
          reference_map, t$reference_tumor_difference_intervals, start, end,
          t$frame_effect, sequence, t$frameshift_start_aa), provenance = provenance)
      peptides[[length(peptides) + 1L]] <- peptide
    }
  }
  answer(if (length(peptides)) "GENERATED" else "NO_ALTERED_LIGAND", peptides)
}
