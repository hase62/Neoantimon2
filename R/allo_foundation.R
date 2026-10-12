allo_foundation <- function(config) {
  contexts <- allo_contexts(config)
  reference <- allo_reference(config$reference)
  typing <- allo_resolve_typing(config, reference)
  universes <- allo_self_universes(config, contexts, reference, typing)
  sources <- unique(vapply(contexts, function(c) c$semantics$source_subject_id, ""))
  for (subject in sources) {
    scenarios <- unique(vapply(Filter(function(t) t$subject_id == subject, typing),
      `[[`, "", "typing_scenario_id"))
    if (length(scenarios) > 1L)
      stop("Multiple source typing scenarios unsupported in Phase A: ", subject)
  }
  occurrences <- allo_records(unlist(lapply(Filter(function(t) t$subject_id %in% sources, typing), function(t) {
    # All compatible mappings are retained; an ambiguous occurrence is a
    # hypothesis, not a choice of genotype or a license to generate peptides.
    lapply(Filter(function(r) !is.null(r$hla_source_sequence_id) && r$allele %in% unlist(t$compatible_alleles), reference$records), function(r) {
      payload <- list(source_subject_id = t$subject_id, allele_hypothesis = r$allele,
        sequence_id = r$hla_source_sequence_id, domains = r$domains,
        coordinate_map = r$coordinate_map, typing_scenario_id = t$typing_scenario_id)
      list(hla_source_occurrence_id = allo_id("occurrence", payload), semantics = payload,
        evaluability = t$evaluability, typing_evidence = list(t), resolution_provenance = list(r$resolution_provenance))
    })
  }), recursive = FALSE))
  occurrences <- lapply(allo_set(lapply(occurrences, `[[`, "hla_source_occurrence_id")), function(id) {
    matches <- Filter(function(o) o$hla_source_occurrence_id == id, occurrences)
    states <- vapply(matches, function(o) canonical_bytes(o$evaluability), "")
    if (length(unique(states)) != 1L) stop("Conflicting source occurrence evaluability in typing scenario")
    out <- matches[[1L]]
    out$typing_evidence <- allo_records(unlist(lapply(matches, `[[`, "typing_evidence"), recursive = FALSE))
    out$resolution_provenance <- allo_set(unlist(lapply(matches, `[[`, "resolution_provenance"), recursive = FALSE))
    out
  })
  missing_source <- setdiff(sources, vapply(typing, `[[`, "", "subject_id"))
  missing_self <- setdiff(vapply(contexts, `[[`, "", "allo_context_id"), vapply(universes, function(u) u$semantics$scenario, ""))
  states <- c(reference$status, vapply(typing, function(t) t$evaluability$status, ""),
    vapply(universes, function(u) u$evaluability$status, ""))
  reasons <- allo_set(c(reference$reason_codes, unlist(lapply(typing, function(t) t$evaluability$reason_codes), recursive = FALSE),
    unlist(lapply(universes, function(u) u$evaluability$reason_codes), recursive = FALSE),
    if (length(missing_source)) list("source_typing_not_evaluated") else list(),
    if (length(missing_self)) list("self_search_not_evaluated") else list()))
  readiness <- if (length(missing_source) || length(missing_self)) "not_evaluated" else if ("failed" %in% states) "failed" else
    if ("ambiguous" %in% states) "ambiguous" else if ("unavailable" %in% states) "unavailable" else
      if ("unsupported" %in% states) "unsupported" else if (any(states != "evaluated")) "insufficient" else "evaluated"
  payload <- list(phase = "A", context_ids = allo_set(lapply(contexts, `[[`, "allo_context_id")),
    sequence_ids = allo_set(lapply(reference$records, `[[`, "hla_source_sequence_id")),
    universe_ids = allo_set(lapply(universes, `[[`, "self_universe_id")),
    reference_policy_id = reference$reference_policy_id,
    self_search_evidence = allo_records(lapply(universes, function(u) list(
      self_universe_id = u$self_universe_id, typing_scenario_id = u$typing_scenario_id,
      catalog_complete = u$catalog_complete, evaluability = u$evaluability))),
    typing = allo_records(lapply(typing, function(t) t[setdiff(names(t), "raw_alleles")])),
    occurrence_ids = allo_set(lapply(occurrences, `[[`, "hla_source_occurrence_id")))
  list(schema_version = allo_version(), specification_id = "neo-allo-indirect-spec-0.1", phase = "A",
    allo_analysis_id = allo_id("analysis", payload), contexts = contexts, reference = reference,
    typing = typing, source_occurrences = occurrences, self_universes = universes,
    readiness = allo_state(readiness, reasons, as.list(c(sources, vapply(contexts, `[[`, "", "allo_context_id"))),
      as.list(c(vapply(typing, `[[`, "", "subject_id"), vapply(universes, function(u) u$semantics$scenario, ""))), "phase_a_foundation"),
    analysis_status = "not_evaluated",
    limitations = list("Foundation readiness is not indirect-HLA analysis or clinical risk", "No peptide generation, prediction, self query search or count"))
}
allo_read <- function(project) {
  if (!allo_enabled(project)) return(NULL)
  bound <- allo_config(project$config$allo)
  expected <- allo_foundation(bound)
  if (is.null(project$allo_foundation) || !identical(canonical_bytes(project$allo_foundation), canonical_bytes(expected)))
    stop("Allo foundation derived contract mismatch")
  project$allo_foundation
}
allo_summary <- function(x) list(schema_version = x$schema_version, phase = "A",
  allo_analysis_id = x$allo_analysis_id, context_ids = lapply(x$contexts, `[[`, "allo_context_id"),
  reference_policy_id = x$reference$reference_policy_id, readiness = x$readiness,
  analysis_status = "not_evaluated", limitations = x$limitations)
allo_results_append <- function(out, project) {
  if (allo_enabled(project)) attr(out, "allo_foundation") <- allo_summary(allo_read(project))
  if (indirect_enabled(project)) {
    b <- indirect_read(project)
    attr(out, "indirect_hla_baseline") <- indirect_summary(b)
    attr(out, "indirect_hla_result") <- indirect_public_result(b, project$allo_foundation)
  }
  if (native_allo_enabled(project)) attr(out, "allo_native_evidence") <- native_allo_summary(native_allo_read(project))
  out
}
# Query-less explanations are allowed only for opt-in Phase A metadata. The
# historical no-allo explanation function and its errors stay byte-equivalent.
neo_explain <- function(project, candidate_id = NULL, reason = NULL,
                        peptide = NULL, sample_id = NULL,
                        format = c("structured", "text")) {
  if (!allo_enabled(project)) return(neo_explain_core(project, candidate_id, reason, peptide, sample_id, format))
  format <- match.arg(format)
  if (allo_enabled(project) && all(vapply(list(candidate_id, reason, peptide, sample_id), is.null, TRUE))) {
    out <- if (format == "text") "" else structure(list(query = list(candidate_id = NULL,
      reason = NULL, peptide = NULL, sample_id = NULL), candidates = list()), class = "neo_explanation")
  } else out <- neo_explain_core(project, candidate_id, reason, peptide, sample_id, format)
  if (!allo_enabled(project)) return(out)
  x <- allo_read(project)
  if (format == "structured") {
    out$allo_foundation <- x
    if (indirect_enabled(project)) {
      b <- indirect_read(project)
      out$indirect_hla_baseline <- b
      out$indirect_hla_explanation <- indirect_public_explanation(b, x)
    }
    if (native_allo_enabled(project)) out$allo_native_evidence <- native_allo_read(project)
    return(out)
  }
  if (indirect_enabled(project)) {
    explanation <- indirect_public_explanation(indirect_read(project), x)
    text <- paste0(out, if (nzchar(out)) "\n" else "", indirect_public_text(explanation))
    if (native_allo_enabled(project)) text <- paste0(text, "\n", native_allo_text(native_allo_read(project)))
    return(text)
  }
  paste0(out, if (nzchar(out)) "\n" else "", "ALLO PHASE A: ", x$allo_analysis_id,
    "; readiness=", x$readiness$status, "; analysis=not_evaluated; no indirect peptide count or clinical interpretation")
}
