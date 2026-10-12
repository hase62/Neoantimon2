# This assesses the declared search space. No query search is performed in A.
allo_self_universes <- function(config, contexts, reference, typing) {
  out <- lapply(config$self_search, function(x) {
    allo_keys(x, c("scenario_id", "self_subject_id", "typing_scenario_id", "scope",
      "required_loci", "required_domains", "catalog_complete", "isoform_policy", "sequence_policy"), label = "self search envelope")
    for (k in c("scenario_id", "self_subject_id", "typing_scenario_id", "scope", "isoform_policy", "sequence_policy"))
      x[[k]] <- allo_text(x[[k]], k)
    context <- Filter(function(c) c$semantics$scenario_id == x$scenario_id, contexts)
    if (length(context) != 1L || x$self_subject_id != context[[1L]]$semantics$self_subject_id) stop("Self-search subject/context binding conflict")
    if (!x$scope %in% c("corresponding_hla", "all_self_hla") ||
      x$isoform_policy != "explicit_mapped_isoforms" || x$sequence_policy != reference$semantics$sequence_policy)
      stop("Unsupported Phase A self-search policy; proteome/pMHC searches deferred")
    x$required_loci <- allo_strings(x$required_loci, "required loci", TRUE)
    x$required_domains <- allo_strings(x$required_domains, "required domains", TRUE)
    if (!all(unlist(x$required_loci) %in% c("A", "B", "C", "E", "F", "G", "DRA", "DRB1", "DRB3", "DRB4", "DRB5", "DQA1", "DQB1", "DPA1", "DPB1"))) stop("Unsupported required self locus")
    if (x$scope == "corresponding_hla" && length(x$required_loci) != 1L) stop("Corresponding-HLA scope requires exactly one declared locus")
    x$catalog_complete <- allo_boolean(x$catalog_complete, "catalog complete")
    own <- Filter(function(t) t$subject_id == x$self_subject_id && t$typing_scenario_id == x$typing_scenario_id && t$locus %in% unlist(x$required_loci), typing)
    examined <- allo_set(lapply(own, `[[`, "locus"))
    missing <- setdiff(unlist(x$required_loci), unlist(examined))
    seqids <- allo_set(unlist(lapply(own, `[[`, "sequence_ids"), recursive = FALSE))
    records <- Filter(function(r) !is.null(r$hla_source_sequence_id) && r$hla_source_sequence_id %in% unlist(seqids), reference$records)
    coverage <- allo_records(lapply(records, function(r) list(sequence_id = r$hla_source_sequence_id,
      domains = r$domains, covered_regions = r$covered_regions)))
    domains_missing <- !length(records) || any(vapply(records, function(r)
      !all(unlist(x$required_domains) %in% vapply(r$domains, `[[`, "", "name")), TRUE))
    reasons <- allo_set(c(if (length(missing)) list("hla_resolution_insufficient") else list(),
      if (!x$catalog_complete) list("self_catalog_incomplete") else list(),
      if (domains_missing) list("self_domain_coverage_incomplete") else list(),
      unlist(lapply(own, function(t) t$evaluability$reason_codes), recursive = FALSE),
      reference$reason_codes))
    ready <- x$catalog_complete && !length(missing) && !domains_missing && reference$status == "evaluated" &&
      length(own) > 0L && all(vapply(own, function(t) t$evaluability$status == "evaluated", TRUE))
    if (!ready) reasons <- allo_set(c(reasons, list("self_search_incomplete")))
    payload <- list(self_subject_id = x$self_subject_id, scenario = context[[1L]]$allo_context_id,
      scope = x$scope,
      required_loci = x$required_loci, required_domains = x$required_domains,
      isoform_policy = x$isoform_policy,
      sequence_policy = x$sequence_policy, sequences = coverage)
    list(schema_version = "allo-self-universe-1.0.0", self_universe_id = allo_id("universe", payload),
      semantics = payload, typing_scenario_id = x$typing_scenario_id,
      catalog_complete = x$catalog_complete, typing_evidence = allo_records(own),
      reference_policy_id = reference$reference_policy_id,
      coverage_complete_for_requested_scope = ready, biological_self_completeness_established = FALSE,
      comparison_status = "not_evaluated", search_performed = FALSE,
      evaluability = allo_state(if (ready) "evaluated" else "insufficient_search", reasons,
        x$required_loci, examined, x$scope),
      limitations = list("No query search performed; no supported absence or immune tolerance claim"))
  })
  if (anyDuplicated(vapply(out, function(u) u$semantics$scenario, ""))) stop("Duplicate self-search envelope for context")
  allo_records(out)
}
