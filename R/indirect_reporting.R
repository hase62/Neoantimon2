# T049-T051: read-only projections of the frozen Phase B computation.
# No namespace, prediction, self search, or scientific selection is defined here.
indirect_report_order <- function(rows, field) {
  if (!length(rows)) return(list())
  rows[order(vapply(rows, function(r) if (length(field) == 1L) r[[field]] else canonical_bytes(r[field]), ""), method = "radix")]
}
indirect_report_one <- function(rows, field, id) {
  found <- Filter(function(r) identical(r[[field]], id), rows)
  if (length(found) != 1L) stop("Indirect reporting reference must resolve exactly once: ", field, " ", id)
  found[[1L]]
}
indirect_result_context <- function(x, foundation, context) {
  a <- indirect_report_one(foundation$contexts, "allo_context_id", context$allo_context_id)
  source <- Filter(function(t) t$subject_id == a$semantics$source_subject_id &&
    t$locus %in% unlist(x$policy$generation$loci), foundation$typing)
  source <- lapply(source, function(t) t[c("subject_id", "typing_scenario_id", "locus",
    "resolved_allele", "compatible_alleles", "sequence_ids", "evaluability")])
  source <- allo_records(source)
  s <- context$summary
  evaluated <- allo_set(lapply(Filter(function(p) p$status == "evaluated", context$presenters), `[[`, "presenter"))
  list(schema_version = "indirect-hla-result-1.0.0", allo_context_id = context$allo_context_id,
    profile = x$profile, profile_version = x$profile_version,
    claim_boundary = x$claim, official_reproduction = FALSE,
    evidence_kind = x$provider$semantics$evidence_kind,
    analysis_status = s$status, completeness_status = if (s$analysis_complete) "complete" else "incomplete",
    complete_count = s$complete_count, conditional_count = s$conditional_included_count,
    conditional_semantics = "inspectable_conditional_result_only",
    included_unit_count = length(s$member_ids), self_excluded_count = s$excluded_self_count,
    binding_policy_excluded_count = s$excluded_nonbinder_count, unevaluable_count = s$unevaluable_count,
    evaluated_presenter_count = length(evaluated), evaluated_presenters = evaluated,
    declared_presenters = indirect_report_order(context$presenters, "presenter"),
    presenter_scope = context$presenter_scope, source_hla_summary = source,
    source_subject_id = a$semantics$source_subject_id, self_subject_id = a$semantics$self_subject_id,
    self_universe_id = context$self_universe_id, allo_policy_id = x$allo_policy_id,
    allo_analysis_id = x$allo_analysis_id, reason_codes = allo_set(c(s$reason_codes,
      foundation$reference$reason_codes,
      unlist(lapply(context$generation$ledger, `[[`, "reason_codes"), recursive = FALSE))))
}
indirect_result_collection <- function(contexts) {
  contexts <- indirect_report_order(contexts, "allo_context_id")
  if (anyDuplicated(vapply(contexts, `[[`, "", "allo_context_id"))) stop("Duplicate public allo context")
  list(schema_version = "indirect-hla-result-1.0.0", specification_id = "neo-allo-indirect-spec-0.1",
    official_reproduction = FALSE, aggregation = "separate_contexts_no_patient_aggregate", contexts = contexts)
}
indirect_public_result <- function(x, foundation) indirect_result_collection(lapply(x$contexts,
  function(context) indirect_result_context(x, foundation, context)))

indirect_reconciliation <- function(context) {
  units <- context$baseline_units
  selected <- function(state) Filter(function(u) u$state == state, units)
  included <- selected("included"); self <- selected("excluded_self")
  included_keys <- allo_set(lapply(included, `[[`, "epitope_key"))
  self_keys <- allo_set(lapply(self, `[[`, "epitope_key"))
  s <- context$summary
  valid <- identical(included_keys, allo_set(s$member_ids)) &&
    identical(length(included_keys), s$conditional_included_count) &&
    identical(length(self_keys), s$excluded_self_count) &&
    identical(length(selected("excluded_nonbinder")), as.integer(s$excluded_nonbinder_count)) &&
    identical(length(selected("unevaluable")), as.integer(s$unevaluable_count)) &&
    if (s$analysis_complete) identical(s$complete_count, length(included_keys)) else is.null(s$complete_count)
  if (!valid) stop("Frozen Phase B count reconciliation failed; scientific contract must be reviewed")
  list(status = "reconciled", included_unique_core_presenter_units = length(included_keys),
    included_pmhc_records = length(included), duplicate_included_union_members = length(included) - length(included_keys),
    self_excluded_unique_core_presenter_units = length(self_keys), self_excluded_pmhc_records = length(self),
    binding_policy_excluded_pmhc_records = length(selected("excluded_nonbinder")),
    unevaluable_pmhc_records = length(selected("unevaluable")),
    complete_count = s$complete_count, conditional_count = s$conditional_included_count,
    semantics = "included and self counts are separate unique epitope-key sets; nonbinder and unevaluable counts are pMHC records; sets may overlap, no conservation equation")
}

indirect_explanation_context <- function(x, foundation, context) {
  # Graph node keys are existing semantic IDs, qualified by type. Decisions use
  # their pMHC key and union nodes use their epitope key: no new biological ID.
  types <- c("context", "source_sequence", "source_occurrence", "peptide_occurrence", "peptide",
    "candidate", "pmhc", "core_assignment", "epitope", "self_comparison", "evidence", "self_universe",
    "policy", "analysis", "decision", "aggregate_unit")
  nodes <- stats::setNames(lapply(types, function(t) list()), types)
  edges <- list()
  add <- function(type, id, data) {
    if (is.null(id)) return(invisible(NULL))
    previous <- nodes[[type]][[id]]
    record <- list(id = id, data = data)
    if (!is.null(previous) && !identical(canonical_bytes(previous), canonical_bytes(record)))
      stop("Conflicting reporting node: ", type, " ", id)
    nodes[[type]][[id]] <<- record
    invisible(NULL)
  }
  link <- function(from_type, from_id, relation, to_type, to_id) {
    if (is.null(from_id) || is.null(to_id)) return(invisible(NULL))
    edges[[length(edges) + 1L]] <<- list(from_type = from_type, from_id = from_id,
      relation = relation, to_type = to_type, to_id = to_id)
    invisible(NULL)
  }
  cid <- context$allo_context_id
  a <- indirect_report_one(foundation$contexts, "allo_context_id", cid)
  add("context", cid, a[c("schema_version", "semantics", "evaluability")])
  add("analysis", x$allo_analysis_id, list(phase = "B", profile = x$profile, profile_version = x$profile_version))
  add("analysis", foundation$allo_analysis_id, list(phase = "A", readiness = foundation$readiness))
  link("context", cid, "reported_by", "analysis", x$allo_analysis_id)
  link("analysis", x$allo_analysis_id, "foundation", "analysis", foundation$allo_analysis_id)
  add("policy", x$allo_policy_id, list(kind = "baseline", semantics = x$policy))
  add("policy", context$generation$generation_policy_id, list(kind = "generation", semantics = context$generation$policy))
  add("policy", foundation$reference$reference_policy_id, list(kind = "reference", semantics = foundation$reference$semantics))
  add("policy", x$provider$provider_configuration_id, list(kind = "provider_configuration", semantics = x$provider$semantics))
  for (pid in c(x$allo_policy_id, foundation$reference$reference_policy_id, x$provider$provider_configuration_id))
    link("analysis", x$allo_analysis_id, "uses_policy", "policy", pid)
  universe <- Filter(function(u) identical(u$self_universe_id, context$self_universe_id), foundation$self_universes)
  for (u in universe) {
    add("self_universe", u$self_universe_id, u)
    link("context", cid, "self_universe", "self_universe", u$self_universe_id)
    link("self_universe", u$self_universe_id, "uses_reference", "policy", u$reference_policy_id)
  }
  # One molecular sequence node may have several allele/reference mappings.
  seqids <- allo_set(lapply(foundation$reference$records, `[[`, "hla_source_sequence_id"))
  for (id in seqids) {
    rr <- Filter(function(r) identical(r$hla_source_sequence_id, id), foundation$reference$records)
    add("source_sequence", id, list(semantics = allo_sequence_identity_payload(foundation$reference$semantics$species, rr[[1L]]),
      mappings = allo_records(lapply(rr, function(r) r[c("allele", "sequence_identifier", "sequence_status", "supported")]))))
    link("source_sequence", id, "reference_policy", "policy", foundation$reference$reference_policy_id)
  }
  for (o in Filter(function(o) o$semantics$source_subject_id == a$semantics$source_subject_id, foundation$source_occurrences)) {
    add("source_occurrence", o$hla_source_occurrence_id, o[c("semantics", "evaluability")])
    link("context", cid, "source_origin", "source_occurrence", o$hla_source_occurrence_id)
    link("source_occurrence", o$hla_source_occurrence_id, "sequence", "source_sequence", o$semantics$sequence_id)
  }
  for (o in context$generation$occurrences) {
    id <- o$indirect_peptide_occurrence_id
    add("peptide_occurrence", id, o[setdiff(names(o), "provenance")])
    add("peptide", o$semantics$peptide_key, list(ligand_sequence = o$peptide_sequence))
    link("source_occurrence", o$semantics$source_occurrence_id, "generates", "peptide_occurrence", id)
    link("peptide_occurrence", id, "sequence", "source_sequence", o$hla_source_sequence_id)
    link("peptide_occurrence", id, "molecular_ligand", "peptide", o$semantics$peptide_key)
    link("peptide_occurrence", id, "generation_policy", "policy", o$semantics$generation_policy_id)
  }
  for (r in x$provider$rows) {
    add("evidence", r$evidence_id, r[setdiff(names(r), c("raw_output", "raw_output_alternatives", "provenance"))])
    link("evidence", r$evidence_id, "provider_configuration", "policy", x$provider$provider_configuration_id)
  }
  for (q in context$comparisons) {
    add("self_comparison", q$self_comparison_id, q[setdiff(names(q), "provenance")])
    link("self_comparison", q$self_comparison_id, "searches", "self_universe", q$semantics$self_universe_id)
    if (q$semantics$query_type == "full_peptide")
      link("self_comparison", q$self_comparison_id, "query", "peptide", q$semantics$query_identity)
    for (oid in q$occurrence_ids) link("self_comparison", q$self_comparison_id, "query_origin", "peptide_occurrence", oid)
    for (h in q$matches) {
      if (!is.null(h$evidence_id)) link("self_comparison", q$self_comparison_id, "supported_by", "evidence", h$evidence_id)
      matches <- if (is.null(h$matches)) list(h) else h$matches
      for (w in matches) link("self_comparison", q$self_comparison_id, "self_sequence", "source_sequence", w$sequence_id)
    }
    for (l in q$search_ledger) link("self_comparison", q$self_comparison_id, "search_evidence", "evidence", l$evidence_id)
  }
  for (p in context$pmhcs) {
    id <- p$indirect_pmhc_id
    add("pmhc", id, p[setdiff(names(p), c("binding_evidence", "core_assignments"))])
    add("candidate", p$candidate_key, list(peptide_key = p$peptide_key, presenter = p$presenter, hla_class = p$hla_class))
    link("context", cid, "presents", "pmhc", id)
    link("pmhc", id, "molecular_candidate", "candidate", p$candidate_key)
    link("candidate", p$candidate_key, "molecular_ligand", "peptide", p$peptide_key)
    link("pmhc", id, "molecular_ligand", "peptide", p$peptide_key)
    link("peptide", p$peptide_key, "presented_as", "pmhc", id)
    link("pmhc", id, "binding_evidence", "evidence", p$binding_evidence$evidence_id)
    link("pmhc", id, "provider_configuration", "policy", p$provider_configuration_id)
    for (oid in p$occurrence_ids) link("pmhc", id, "source_support", "peptide_occurrence", oid)
    for (ca in p$core_assignments) {
      add("core_assignment", ca$indirect_core_assignment_id, ca)
      link("pmhc", id, "core_assignment", "core_assignment", ca$indirect_core_assignment_id)
      link("core_assignment", ca$indirect_core_assignment_id, "provider_configuration", "policy", ca$semantics$provider_configuration_id)
      if (!is.null(ca$epitope_key)) {
        add("epitope", ca$epitope_key, list(core_sequence = ca$semantics$core_sequence, presenter = p$presenter))
        link("core_assignment", ca$indirect_core_assignment_id, "molecular_core_presenter", "epitope", ca$epitope_key)
      }
    }
  }
  groups <- lapply(context$summary$member_ids, function(key) {
    uu <- Filter(function(u) u$state == "included" && identical(u$epitope_key, key), context$baseline_units)
    ids <- allo_set(lapply(uu, `[[`, "indirect_pmhc_id"))
    pp <- lapply(ids, function(id) indirect_report_one(context$pmhcs, "indirect_pmhc_id", id))
    occurrence_ids <- allo_set(unlist(lapply(uu, `[[`, "occurrence_ids"), recursive = FALSE))
    oo <- lapply(occurrence_ids, function(id) indirect_report_one(context$generation$occurrences, "indirect_peptide_occurrence_id", id))
    ep <- nodes$epitope[[key]]$data
    record <- list(counted_unit_key = key, presenter = ep$presenter, core_sequence = ep$core_sequence,
      supporting_pmhc_ids = ids, peptide_keys = allo_set(lapply(pp, `[[`, "peptide_key")),
      full_ligands = allo_set(lapply(pp, `[[`, "peptide_sequence")), occurrence_ids = occurrence_ids,
      source_occurrence_ids = allo_set(lapply(oo, function(o) o$semantics$source_occurrence_id)),
      supporting_evidence_ids = allo_set(lapply(pp, function(p) p$binding_evidence$evidence_id)),
      threshold_policy = x$policy$threshold, inclusion_reason = "supported_nonself_within_scope_under_frozen_policy")
    add("aggregate_unit", key, record)
    link("aggregate_unit", key, "counted_by", "analysis", x$allo_analysis_id)
    link("aggregate_unit", key, "within_context", "context", cid)
    link("aggregate_unit", key, "counted_core_presenter", "epitope", key)
    for (id in ids) link("aggregate_unit", key, "union_member", "pmhc", id)
    record
  })
  groups <- indirect_report_order(groups, "counted_unit_key")
  for (u in context$baseline_units) {
    p <- indirect_report_one(context$pmhcs, "indirect_pmhc_id", u$indirect_pmhc_id)
    r <- p$binding_evidence
    aggregate <- if (u$state == "included") indirect_report_one(groups, "counted_unit_key", u$epitope_key) else NULL
    duplicate <- !is.null(aggregate) && u$indirect_pmhc_id != aggregate$supporting_pmhc_ids[[1L]]
    comparison <- if (is.null(u$self_comparison_id)) NULL else
      indirect_report_one(context$comparisons, "self_comparison_id", u$self_comparison_id)
    outcome <- if (u$state == "excluded_self") "supported_self_match" else
      if (u$state == "excluded_nonbinder") "provider_not_binding_under_policy" else
        if (u$state == "included") "supported_nonself_within_scope" else
          if (foundation$reference$status == "failed") "reference_failure" else
          if ("presenter_out_of_profile_scope" %in% unlist(u$reason_codes)) "presenter_out_of_profile" else
            if (r$status == "failed") "provider_failure" else
              if (r$status != "evaluated") "provider_unevaluable" else
                if (any(c("primary_core_ambiguous", "primary_core_unresolved") %in% unlist(u$reason_codes))) "core_assignment_unevaluable" else
                  if (foundation$reference$status == "failed") "reference_failure" else
                    if (is.null(context$self_universe_id)) "missing_self_universe" else
                      if (!is.null(comparison) && comparison$status == "ambiguous") "ambiguous_self_typing" else
                        if (!is.null(comparison) && comparison$status == "insufficient_search") "insufficient_self_search" else "presenter_resolution_insufficient"
    explanation_codes <- allo_set(c(u$reason_codes, list(outcome),
      if (duplicate) list("duplicate_core_presenter_union_member") else list()))
    decision <- c(u, list(outcome_kind = outcome, explanation_codes = explanation_codes,
      affinity_nM = r$affinity_nM, provider_status = r$status, provider_evidence_id = r$evidence_id,
      self_comparison_status = if (is.null(comparison)) NULL else comparison$status,
      threshold_policy = p$threshold_policy,
      affinity_under_frozen_policy = if (r$status == "evaluated") r$affinity_nM < p$threshold_policy$value else NULL,
      aggregation_role = if (is.null(aggregate)) "uncounted" else if (duplicate) "duplicate_union_member" else "union_representative",
      counted_unit_key = if (u$state == "included") u$epitope_key else NULL))
    add("decision", u$indirect_pmhc_id, decision)
    link("pmhc", u$indirect_pmhc_id, "decision", "decision", u$indirect_pmhc_id)
    link("evidence", r$evidence_id, "informs_decision", "decision", u$indirect_pmhc_id)
    link("decision", u$indirect_pmhc_id, "self_comparison", "self_comparison", u$self_comparison_id)
    link("decision", u$indirect_pmhc_id, "assigned_core_presenter", "epitope", u$epitope_key)
    if (u$state == "included") link("decision", u$indirect_pmhc_id, "aggregation_membership", "aggregate_unit", u$epitope_key)
  }
  for (q in context$comparisons) if (q$semantics$query_type == "core_pmhc") {
    link("self_comparison", q$self_comparison_id, "query", "epitope", q$semantics$query_identity)
    link("epitope", q$semantics$query_identity, "compared_by", "self_comparison", q$self_comparison_id)
  }
  # Convert keyed maps to arrays in semantic-ID order, independent of contents.
  nodes <- lapply(nodes, function(table) {
    if (!length(table)) return(list())
    unname(table[sort(names(table), method = "radix")])
  })
  edges <- allo_records(edges)
  graph <- list(nodes = nodes, edges = edges)
  indirect_graph_integrity(graph)
  list(allo_context_id = cid, summary = indirect_result_context(x, foundation, context),
    graph = graph, included_units = groups,
    decisions = lapply(nodes$decision, `[[`, "data"),
    generation_coverage = allo_records(context$generation$ledger), self_searches = indirect_report_order(context$self_searches, "presenter"),
    unidentified_self_comparisons = allo_records(Filter(function(q) is.null(q$self_comparison_id), context$comparisons)),
    reference_status = foundation$reference$status, reference_reason_codes = foundation$reference$reason_codes,
    reconciliation = indirect_reconciliation(context))
}
indirect_graph_integrity <- function(graph) {
  keys <- lapply(graph$nodes, function(table) vapply(table, `[[`, "", "id"))
  if (any(vapply(keys, anyDuplicated, 0L) > 0L)) stop("Duplicate evidence graph node")
  for (edge in graph$edges) {
    if (!edge$from_id %in% keys[[edge$from_type]] || !edge$to_id %in% keys[[edge$to_type]])
      stop("Orphan indirect evidence graph reference: ", edge$relation)
  }
  # Check embedded semantic references as well as materialized graph edges.
  namespaces <- c("allo-context-id-v1", "hla-source-sequence-id-v1", "hla-source-occurrence-id-v1",
    "indirect-peptide-occurrence-id-v1", "peptide_key", "candidate_key", "indirect-pmhc-id-v1",
    "indirect-core-assignment-id-v1", "epitope_key", "self-comparison-id-v1", "indirect-evidence-id-v1",
    "self-universe-id-v1", "allo-policy-id-v1", "allo-analysis-id-v1")
  values <- unlist(graph$nodes, recursive = TRUE, use.names = FALSE)
  references <- unique(values[grepl(paste0("^(", paste(namespaces, collapse = "|"), "):[0-9a-f]{64}$"), values)])
  if (length(setdiff(references, unlist(keys, use.names = FALSE)))) stop("Orphan embedded indirect semantic reference")
  invisible(TRUE)
}
indirect_public_explanation <- function(x, foundation) {
  contexts <- indirect_explanation_contexts(lapply(x$contexts, function(c) indirect_explanation_context(x, foundation, c)))
  rows <- indirect_report_order(x$provider$rows, "evidence_id")
  list(schema_version = "indirect-hla-explanation-1.0.0", specification_id = x$specification_id,
    profile = x$profile, profile_version = x$profile_version, claim_boundary = x$claim,
    official_reproduction = FALSE, contexts = contexts, compatibility = indirect_compatibility(),
    provenance = list(reference = foundation$reference$provenance, provider = x$provider$provenance,
      provider_artifact = x$provenance,
      contexts = lapply(indirect_report_order(foundation$contexts, "allo_context_id"),
        function(c) list(allo_context_id = c$allo_context_id, provenance = c$provenance)),
      sequence_mappings = lapply(indirect_report_order(foundation$reference$records, c("allele", "sequence_identifier")),
        function(r) r[c("allele", "sequence_identifier", "hla_source_sequence_id", "resolution_provenance")]),
      reference_resolution = lapply(indirect_report_order(foundation$source_occurrences, "hla_source_occurrence_id"),
        function(o) list(hla_source_occurrence_id = o$hla_source_occurrence_id, resolution_provenance = o$resolution_provenance)),
      evidence = lapply(rows, function(r) r[c("evidence_id", "raw_output", "raw_output_alternatives", "provenance")])),
    limitations = x$limitations)
}
indirect_explanation_contexts <- function(contexts) {
  contexts <- indirect_report_order(contexts, "allo_context_id")
  if (anyDuplicated(vapply(contexts, `[[`, "", "allo_context_id"))) stop("Duplicate explanation allo context")
  contexts
}
indirect_public_text <- function(explanation) {
  lines <- vapply(explanation$contexts, function(c) paste0(c$allo_context_id,
    "; completeness=", c$summary$completeness_status,
    "; complete_count=", if (is.null(c$summary$complete_count)) "unavailable (NULL)" else c$summary$complete_count,
    "; conditional_count=", c$summary$conditional_count, " (inspectable conditional result only)",
    "; self_excluded=", c$summary$self_excluded_count, "; unevaluable=", c$summary$unevaluable_count,
    "; reasons=", paste(unlist(c$summary$reason_codes), collapse = ",")), "")
  paste(c(paste0("INDIRECT HLA: ", explanation$profile, " version ", explanation$profile_version,
    "; ", explanation$claim_boundary), lines), collapse = "\n")
}
