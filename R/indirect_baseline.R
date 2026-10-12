indirect_context <- function(foundation, config, context, provider) {
  c <- context$semantics
  if (c$transplant_type != "SOT" || c$scenario_id != "sot_recipient_response" || c$direction != "donor_to_recipient")
    stop("Strict literature-derived profile supports SOT recipient response only")
  if (foundation$reference$semantics$sequence_policy != "precursor") stop("Strict profile requires precursor source/self proteins")
  generation <- indirect_generate(foundation, context)
  catalog <- indirect_self_catalog(foundation, context, provider)
  presenters <- indirect_presenters(config, foundation, context)
  presenter_scope <- indirect_presenter_scope(foundation, context, presenters)
  pmhcs <- list(); comparisons <- list(); units <- list(); searches <- list()
  for (p in presenters) {
    search <- indirect_self_prediction(catalog, p$presenter, provider)
    if (p$status != "evaluated") search$complete <- FALSE
    searches[[length(searches) + 1L]] <- list(presenter = p$presenter, complete = search$complete,
      hits = search$hits, ledger = search$ledger)
    ligands <- allo_set(lapply(Filter(function(o) o$status == "evaluated", generation$occurrences), `[[`, "peptide_sequence"))
    for (ligand in ligands) {
      origins <- Filter(function(o) o$peptide_sequence == ligand && o$status == "evaluated", generation$occurrences)
      peptide_key <- origins[[1L]]$semantics$peptide_key
      presenter <- peptide_presenters(list(p$presenter))[[1L]]
      # Out-of-scope/invalid molecules remain an explicit failed stratum, never a
      # fabricated Class II identity or a zero-count presenter.
      if (presenter$status != "RESOLVED") next
      key <- biological_id("candidate_key", list(species = foundation$reference$semantics$species,
        hla_class = presenter$hla_class, ligand_sequence = ligand, mhc_key = presenter$mhc_key), peptide_contract())
      payload <- indirect_pmhc_payload(context$allo_context_id, key, c$apc_compartment, c$scenario_id)
      id <- indirect_id("indirect-pmhc-id-v1", payload)
      row <- indirect_query(provider, ligand, p$presenter)
      evidence <- row
      if (p$status != "evaluated") row$status <- p$status
      assignments <- lapply(evidence$cores, function(core) {
        cp <- list(indirect_pmhc_id = id, core_sequence = core$core_sequence, offset = core$offset,
          register = core$register, orientation = core$orientation,
          provider_configuration_id = provider$provider_configuration_id,
          assignment_semantics = list(primary = core$primary, status = core$status))
        list(indirect_core_assignment_id = indirect_id("indirect-core-assignment-id-v1", cp),
          mismatch_annotation = list(status = "not_evaluated", reason_codes = list("alignment_not_performed")),
          semantics = cp, epitope_key = if (core$orientation == "forward" && core$status == "evaluated")
            biological_id("epitope_key", list(species = foundation$reference$semantics$species,
              hla_class = "II", binding_core = core$core_sequence, mhc_key = p$presenter), peptide_contract()) else NULL)
      })
      links <- allo_set(lapply(origins, `[[`, "indirect_peptide_occurrence_id"))
      pmhcs[[length(pmhcs) + 1L]] <- list(schema_version = "indirect-pmhc-evidence-1.0.0",
        indirect_pmhc_id = id, semantics = payload, candidate_key = key, peptide_key = peptide_key,
        peptide_sequence = ligand, presenter = p$presenter, hla_class = presenter$hla_class,
        occurrence_ids = links, provider_configuration_id = provider$provider_configuration_id,
        binding_evidence = evidence, threshold_policy = indirect_threshold(),
        core_assignments = allo_records(assignments), status = row$status)
      full <- indirect_compare(peptide_key, "full_peptide", ligand, NULL, catalog, search, provider, links)
      comparisons[[length(comparisons) + 1L]] <- full
      state <- "unevaluable"; reason <- row$reason_codes; core_key <- NULL; comparison_id <- NULL
      if (row$status == "evaluated") {
        if (row$affinity_nM >= 1000) { state <- "excluded_nonbinder"; reason <- list("strict_affinity_not_met") }
        else {
          primary <- indirect_primary(row)
          reason <- primary$reason_codes
          if (primary$status == "evaluated") {
            core_key <- biological_id("epitope_key", list(species = foundation$reference$semantics$species,
              hla_class = "II", binding_core = primary$core$core_sequence, mhc_key = p$presenter), peptide_contract())
            comparison <- indirect_compare(core_key, "core_pmhc", primary$core$core_sequence, p$presenter,
              catalog, search, provider, links)
            comparison_id <- comparison$self_comparison_id
            comparisons[[length(comparisons) + 1L]] <- comparison
            state <- if (comparison$status == "self_match") "excluded_self" else
              if (comparison$status == "supported_nonself_within_scope") "included" else "unevaluable"
            reason <- comparison$reason_codes
          }
        }
      } else reason <- allo_set(c(reason, p$reason_codes, list("binding_unevaluable")))
      units[[length(units) + 1L]] <- list(indirect_pmhc_id = id, epitope_key = core_key,
        self_comparison_id = comparison_id, occurrence_ids = links, presenter = p$presenter,
        state = state, reason_codes = allo_set(reason))
    }
  }
  # A comparison belongs to a molecular query/universe, not an occurrence. Union
  # relational links after normalizing equivalent queries from overlapping windows.
  comparison_keys <- vapply(comparisons, function(r) canonical_bytes(list(
    self_subject_id = r$self_subject_id, semantics = r$semantics)), "")
  comparisons <- lapply(sort(unique(comparison_keys), method = "radix"), function(key) {
    rr <- comparisons[comparison_keys == key]
    out <- rr[[1L]]; out$occurrence_ids <- allo_set(unlist(lapply(rr, `[[`, "occurrence_ids"), recursive = FALSE)); out
  })
  included <- allo_set(lapply(Filter(function(u) u$state == "included", units), `[[`, "epitope_key"))
  self <- allo_set(lapply(Filter(function(u) u$state == "excluded_self", units), `[[`, "epitope_key"))
  failures <- Filter(function(u) u$state == "unevaluable", units)
  complete <- generation$complete && catalog$complete && presenter_scope$complete && length(presenters) > 0L &&
    all(vapply(presenters, function(p) p$status == "evaluated", TRUE)) &&
    all(vapply(searches, `[[`, TRUE, "complete")) && !length(failures)
  summary <- function(selected, is_complete) {
    members <- allo_set(lapply(Filter(function(u) u$state == "included", selected), `[[`, "epitope_key"))
    list(complete_count = if (is_complete) length(members) else NULL, conditional_included_count = length(members),
      member_ids = members, excluded_self_count = length(allo_set(lapply(Filter(function(u) u$state == "excluded_self", selected), `[[`, "epitope_key"))),
      excluded_nonbinder_count = sum(vapply(selected, function(u) u$state == "excluded_nonbinder", TRUE)),
      unevaluable_count = sum(vapply(selected, function(u) u$state == "unevaluable", TRUE)),
      ligand_diagnostic_count = length(allo_set(lapply(Filter(function(u) u$state == "included", selected), `[[`, "indirect_pmhc_id"))),
      analysis_complete = is_complete, status = if (is_complete) "evaluated" else "insufficient",
      conditional_count_is_not_a_lower_bound = TRUE)
  }
  attributed <- list()
  for (allele in allo_set(lapply(generation$occurrences, `[[`, "source_allele_hypothesis"))) for (p in presenters) {
    ids <- unlist(lapply(Filter(function(o) o$source_allele_hypothesis == allele, generation$occurrences), `[[`, "indirect_peptide_occurrence_id"))
    selected <- Filter(function(u) u$presenter == p$presenter && any(unlist(u$occurrence_ids) %in% ids), units)
    attributed[[length(attributed) + 1L]] <- c(list(source_allele = allele, presenter = p$presenter), summary(selected, complete))
  }
  list(allo_context_id = context$allo_context_id, generation = generation, presenters = presenters,
    presenter_scope = presenter_scope,
    self_universe_id = if (is.null(catalog$universe)) NULL else catalog$universe$self_universe_id,
    self_searches = allo_records(searches), pmhcs = allo_records(pmhcs), comparisons = allo_records(comparisons),
    baseline_units = allo_records(units), source_presenter_summaries = allo_records(attributed),
    summary = c(summary(units, complete), list(reason_codes = allo_set(c(
      if (!generation$complete) list("source_generation_incomplete") else list(),
      if (!catalog$complete) list("self_search_incomplete") else list(),
      if (!length(presenters)) list("presenter_not_evaluated") else list(),
      presenter_scope$reason_codes,
      unlist(lapply(presenters, `[[`, "reason_codes"), recursive = FALSE),
      if (any(!vapply(searches, `[[`, TRUE, "complete"))) list("self_prediction_incomplete") else list(),
      unlist(lapply(failures, `[[`, "reason_codes"), recursive = FALSE))),
      coverage = list(source_generation = generation$ledger, self_search = lapply(searches, `[[`, "ledger")),
      source_attributed_counts_may_overlap = TRUE)))
}
indirect_baseline <- function(foundation, config) {
  # Preserve raw declarations in bound configuration; identity uses the frozen
  # canonical alpha/beta molecule rather than input spelling or chain order.
  config$presenters <- allo_records(lapply(config$presenters, indirect_presenter_semantics))
  if (any(!vapply(config$presenters, function(p) p$scenario_id %in%
    vapply(foundation$contexts, function(c) c$semantics$scenario_id, ""), TRUE)))
    stop("Unknown Phase B presenter scenario")
  provider <- indirect_provider(config$provider$manifest)
  policy <- list(profile = indirect_profile(), version = "1.0.0", generation = indirect_generation_policy(),
    reference_policy_id = foundation$reference$reference_policy_id, presenters = allo_records(config$presenters),
    provider_configuration_id = provider$provider_configuration_id, threshold = indirect_threshold(),
    self = "primary_forward_core_same_presenter_from_qualifying_self_15mers",
    counting = "unique_core_presenter_union", uncertainty = "fail_closed_no_imputation")
  policy_id <- allo_id("policy", policy)
  contexts <- allo_records(lapply(foundation$contexts, function(c) indirect_context(foundation, config, c, provider)))
  payload <- list(context_ids = allo_set(lapply(contexts, `[[`, "allo_context_id")),
    sequence_ids = allo_set(lapply(foundation$reference$records, `[[`, "hla_source_sequence_id")),
    universe_ids = allo_set(lapply(foundation$self_universes, `[[`, "self_universe_id")),
    policy_id = policy_id, provider_configuration_ids = list(provider$provider_configuration_id),
    input_evidence_ids = allo_set(c(list(foundation$allo_analysis_id), lapply(provider$rows, `[[`, "evidence_id"))))
  list(schema_version = "indirect-hla-baseline-1.0.0", phase = "B", profile = indirect_profile(),
    profile_version = "1.0.0", specification_id = "neo-allo-indirect-spec-0.1",
    claim = "Neoantimon2 literature-derived molecular baseline; not official PIRCHE reproduction or a clinical risk score",
    allo_analysis_id = allo_id("analysis", payload), allo_policy_id = policy_id, policy = policy,
    provider = provider, contexts = contexts,
    summaries = lapply(contexts, function(c) c(list(allo_context_id = c$allo_context_id), c$summary)),
    provenance = list(artifact = config$provider[c("path", "sha256")]),
    limitations = list("Predicted binding is not observed presentation or T-cell recognition",
      "Scoped self comparison does not establish immune tolerance", "No official PIRCHE or clinical reproduction"))
}
