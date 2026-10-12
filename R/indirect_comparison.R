indirect_self_catalog <- function(foundation, context, provider) {
  found <- Filter(function(u) u$semantics$scenario == context$allo_context_id, foundation$self_universes)
  if (!length(found)) return(list(universe = NULL, windows = list(), complete = FALSE,
    self_subject_id = context$semantics$self_subject_id, reason_codes = list("self_search_not_evaluated")))
  u <- found[[1L]]
  contract <- u$semantics$scope == "all_self_hla" &&
    identical(allo_set(u$semantics$required_loci), indirect_loci()) &&
    u$semantics$sequence_policy == "precursor"
  if (!contract) stop("Strict profile requires all-self-HLA A/B/C/DRB1/DQB1 precursor search scope")
  windows <- list(); complete <- u$coverage_complete_for_requested_scope
  reasons <- u$evaluability$reason_codes
  for (t in u$typing_evidence) {
    if (t$evaluability$status != "evaluated") { complete <- FALSE; next }
    for (r in Filter(function(r) !is.null(r$hla_source_sequence_id) &&
      r$hla_source_sequence_id %in% unlist(t$sequence_ids) &&
      r$allele %in% unlist(t$compatible_alleles), foundation$reference$records)) {
      g <- indirect_sequence_windows(r)
      if (!g$status %in% c("evaluated", "not_applicable")) complete <- FALSE
      reasons <- c(reasons, g$reason_codes)
      for (w in g$windows) if (w$status == "evaluated") windows[[length(windows) + 1L]] <- list(
        peptide_sequence = w$peptide_sequence, sequence_id = r$hla_source_sequence_id,
        allele = r$allele, start = w$start, end = w$end)
    }
  }
  list(universe = u, self_subject_id = u$semantics$self_subject_id,
    windows = allo_records(windows), complete = complete, reason_codes = allo_set(reasons))
}
indirect_self_prediction <- function(catalog, presenter, provider) {
  hits <- list(); ledger <- list(); complete <- catalog$complete
  for (peptide in allo_set(lapply(catalog$windows, `[[`, "peptide_sequence"))) {
    row <- indirect_query(provider, peptide, presenter)
    state <- row$status
    primary <- NULL
    if (state == "evaluated" && row$affinity_nM < 1000) {
      primary <- indirect_primary(row)
      state <- primary$status
      if (state == "evaluated") hits[[length(hits) + 1L]] <- list(
        core_sequence = primary$core$core_sequence, presenter = presenter,
        peptide_sequence = peptide, evidence_id = row$evidence_id,
        matches = allo_records(Filter(function(w) w$peptide_sequence == peptide, catalog$windows)))
    }
    if (state != "evaluated") complete <- FALSE
    ledger[[length(ledger) + 1L]] <- list(peptide_sequence = peptide, presenter = presenter,
      status = state, reason_codes = allo_set(c(row$reason_codes, if (!is.null(primary)) primary$reason_codes else list())),
      evidence_id = row$evidence_id)
  }
  list(hits = allo_records(hits), ledger = allo_records(ledger), complete = complete)
}
indirect_compare <- function(query_identity, type, sequence, presenter, catalog, search,
                             provider, occurrence_ids = list()) {
  algorithm <- list(version = "exact-literal-same-presenter-primary-core-1.0.0",
    matching = "literal_exact", wildcard = FALSE, il_equivalence = FALSE,
    comparison_level = type, presenter = presenter)
  self_policy <- if (type == "core_pmhc") list(provider_configuration_id = provider$provider_configuration_id,
    threshold = indirect_threshold(), registers = "designated_primary_forward_only") else NULL
  payload <- indirect_comparison_payload(query_identity, type,
    if (is.null(catalog$universe)) NULL else catalog$universe$self_universe_id, algorithm, self_policy)
  matches <- if (type == "core_pmhc") Filter(function(h) h$core_sequence == sequence && h$presenter == presenter, search$hits) else
    Filter(function(w) grepl(sequence, w$peptide_sequence, fixed = TRUE), catalog$windows)
  complete <- if (type == "core_pmhc") search$complete else catalog$complete
  ambiguous <- !is.null(catalog$universe) && any(vapply(catalog$universe$typing_evidence,
    function(t) t$evaluability$status == "ambiguous", TRUE))
  status <- if (length(matches)) "self_match" else if (complete) "supported_nonself_within_scope" else
    if (ambiguous) "ambiguous" else if (is.null(catalog$universe)) "not_evaluated" else "insufficient_search"
  reasons <- allo_set(c(catalog$reason_codes,
    if (!complete) list("self_search_incomplete") else list(),
    if (type == "core_pmhc") unlist(lapply(search$ledger, `[[`, "reason_codes"), recursive = FALSE) else list()))
  list(schema_version = "self-comparison-1.0.0",
    self_comparison_id = if (is.null(catalog$universe)) NULL else indirect_id("self-comparison-id-v1", payload), semantics = payload,
    self_subject_id = catalog$self_subject_id,
    scope = "all_self_hla", status = status, search_complete = complete,
    biological_self_completeness_established = FALSE, reason_codes = reasons,
    matches = allo_records(matches), occurrence_ids = allo_set(occurrence_ids),
    search_ledger = if (type == "core_pmhc") search$ledger else list(),
    provenance = list(reference_policy_id = if (is.null(catalog$universe)) NULL else catalog$universe$reference_policy_id))
}
