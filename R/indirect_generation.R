indirect_sequence_windows <- function(record) {
  # Never concatenate intervals or infer residues. Known partial windows are
  # retained as insufficient candidates; they are not eligible provider queries.
  seq <- record$protein_sequence
  if (is.null(seq)) return(list(windows = list(), blocked = list(), status = "unavailable",
    reason_codes = list("source_sequence_unavailable")))
  n <- nchar(seq)
  if (n < 15L) return(list(windows = list(), blocked = list(), status = "not_applicable",
    reason_codes = list("not_applicable_length")))
  coverage <- unique(unlist(lapply(record$covered_regions, function(d) seq.int(d$start, d$end))))
  windows <- list(); blocked <- list()
  for (start in seq_len(n - 14L)) {
    end <- start + 14L; peptide <- substring(seq, start, end)
    domains <- allo_records(Filter(function(d) d$start <= end && d$end >= start, record$domains))
    valid <- all(seq.int(start, end) %in% coverage) && !grepl("X", peptide)
    if (!valid) {
      blocked[[length(blocked) + 1L]] <- list(start = start, end = end,
        status = "insufficient", reason_codes = list("window_coverage_unavailable"))
    } else windows[[length(windows) + 1L]] <- list(start = start, end = end,
      peptide_sequence = peptide, domains = domains,
      status = if (record$sequence_status == "complete") "evaluated" else "insufficient")
  }
  list(windows = windows, blocked = blocked,
    status = if (record$sequence_status == "complete" && !length(blocked)) "evaluated" else "insufficient",
    reason_codes = if (record$sequence_status == "complete" && !length(blocked)) list() else list("source_sequence_partial"))
}
indirect_generate <- function(foundation, context) {
  policy <- indirect_generation_policy()
  policy_id <- allo_id("policy", policy)
  subject <- context$semantics$source_subject_id
  typing <- Filter(function(t) t$subject_id == subject && t$locus %in% unlist(indirect_loci()), foundation$typing)
  missing <- setdiff(unlist(indirect_loci()), vapply(typing, `[[`, "", "locus"))
  occurrences <- list(); ledger <- list()
  for (locus in missing) ledger[[length(ledger) + 1L]] <- list(locus = locus,
    status = "insufficient", reason_codes = list("hla_resolution_insufficient"), blocked_windows = list())
  for (t in typing) if (t$evaluability$status != "evaluated") ledger[[length(ledger) + 1L]] <- list(
    locus = t$locus, status = t$evaluability$status, reason_codes = t$evaluability$reason_codes, blocked_windows = list())
  source <- Filter(function(o) o$semantics$source_subject_id == subject &&
    any(vapply(typing, function(t) o$semantics$allele_hypothesis %in% unlist(t$compatible_alleles), TRUE)), foundation$source_occurrences)
  for (o in source) {
    rr <- Filter(function(r) identical(r$hla_source_sequence_id, o$semantics$sequence_id) &&
      r$allele == o$semantics$allele_hypothesis, foundation$reference$records)
    r <- rr[[1L]]
    # Ambiguous and unsupported hypotheses are not chosen. Partial known segments
    # are inspectable, but remain insufficient and cannot enter the counted set.
    if (!o$evaluability$status %in% c("evaluated", "insufficient") ||
      (o$evaluability$status == "insufficient" && !"source_sequence_partial" %in% unlist(o$evaluability$reason_codes))) next
    g <- indirect_sequence_windows(r)
    if (o$evaluability$status != "evaluated") {
      g$status <- "insufficient"
      g$reason_codes <- allo_set(c(g$reason_codes, o$evaluability$reason_codes))
      g$windows <- lapply(g$windows, function(w) { w$status <- "insufficient"; w })
    }
    ledger[[length(ledger) + 1L]] <- list(source_occurrence_id = o$hla_source_occurrence_id,
      status = g$status, reason_codes = g$reason_codes, blocked_windows = g$blocked)
    for (w in g$windows) {
      key <- biological_id("peptide_key", list(ligand_sequence = w$peptide_sequence), peptide_contract())
      payload <- indirect_occurrence_payload(o$hla_source_occurrence_id, w$start, w$end, key, policy_id)
      occurrences[[length(occurrences) + 1L]] <- list(schema_version = "indirect-peptide-occurrence-1.0.0",
        indirect_peptide_occurrence_id = indirect_id("indirect-peptide-occurrence-id-v1", payload),
        semantics = payload, hla_source_sequence_id = r$hla_source_sequence_id,
        source_subject_id = subject, source_allele_hypothesis = o$semantics$allele_hypothesis,
        sequence_coordinates = "precursor", peptide_sequence = w$peptide_sequence,
        domains = w$domains, status = w$status,
        mismatch_annotation = list(status = "not_evaluated", reason_codes = list("alignment_not_performed")),
        provenance = list(reference_policy_id = foundation$reference$reference_policy_id,
          source_occurrence_id = o$hla_source_occurrence_id))
    }
  }
  complete <- !length(missing) && length(typing) > 0L &&
    all(vapply(typing, function(t) t$evaluability$status == "evaluated", TRUE)) &&
    all(vapply(ledger, function(l) l$status %in% c("evaluated", "not_applicable"), TRUE))
  list(policy = policy, generation_policy_id = policy_id, occurrences = allo_records(occurrences),
    ledger = allo_records(ledger), complete = complete)
}
