longitudinal_cell <- function(settings, sample, domain, key, observations = list()) {
  a <- Filter(function(a) a$sample_id == sample && a$domain == domain && a$feature_key == key, settings$evaluability)
  assay <- if (length(a)) a[[1L]] else NULL
  state <- if (is.null(assay)) "not_evaluated" else switch(assay$state,
    evaluated_positive = "observed", evaluated_negative = "not_observed_evaluable",
    insufficient_evidence = "not_observed_insufficient", not_evaluated = "not_evaluated",
    not_applicable = "not_applicable", ambiguous = "ambiguous")
  if (length(observations)) state <- if (state == "not_observed_evaluable" || state == "not_applicable") "ambiguous" else "observed"
  list(sample_id = sample, state = state, evaluability = assay, observations = observations)
}
longitudinal_transition <- function(a, b, contexts, relations) {
  ca <- Filter(function(s) s$sample_id == a$sample_id, contexts)[[1L]]
  cb <- Filter(function(s) s$sample_id == b$sample_id, contexts)[[1L]]
  ta <- longitudinal_time(ca); tb <- longitudinal_time(cb)
  replicate_group <- a$sample_id
  repeat {
    expanded <- unique(c(replicate_group, unlist(lapply(Filter(function(r)
      r$type == "technical_replicate" && any(c(r$from, r$to) %in% replicate_group), relations), function(r) c(r$from, r$to)))))
    if (setequal(expanded, replicate_group)) break
    replicate_group <- expanded
  }
  replicate <- b$sample_id %in% replicate_group ||
    (!is.null(ca$biological_context_id) && identical(ca$biological_context_id, cb$biological_context_id))
  if (replicate || is.null(ta$value) || is.null(tb$value) || !identical(ta$axis, tb$axis) || ta$value >= tb$value) return(NULL)
  state <- if (a$state == "observed" && b$state == "not_observed_evaluable") "loss_with_supported_absence" else
    if (a$state == "not_observed_evaluable" && b$state == "observed") "gain_with_prior_supported_absence" else
    if (a$state == "observed" && b$state == "observed") "observed_at_both_timepoints" else "unresolved"
  list(from = a$sample_id, to = b$sample_id, state = state)
}
longitudinal_summary <- function(cells, contexts, relations) {
  states <- vapply(cells, `[[`, "", "state")
  observed <- sum(states == "observed"); evaluated <- sum(states %in% c("observed", "not_observed_evaluable"))
  transitions <- list()
  if (length(cells) > 1L) for (i in seq_len(length(cells) - 1L)) for (j in seq.int(i + 1L, length(cells))) {
    t <- longitudinal_transition(cells[[i]], cells[[j]], contexts, relations)
    if (!is.null(t)) transitions[[length(transitions) + 1L]] <- t
  }
  decisive <- Filter(function(t) t$state %in% c("loss_with_supported_absence", "gain_with_prior_supported_absence"), transitions)
  complete <- all(states %in% c("observed", "not_observed_evaluable", "not_applicable"))
  state <- if (!complete || !observed) "insufficient_to_classify" else
    if (length(decisive)) if (length(unique(vapply(decisive, `[[`, "", "state"))) == 1L) decisive[[1L]]$state else "mixed_observation" else
    if (observed == evaluated) "observed_all_evaluable_samples" else if (observed > 1L) "observed_multiple_samples" else "observed_single_sample"
  regions <- Filter(function(s) !is.null(s$region_id), contexts)
  spatial <- if (length(unique(vapply(regions, `[[`, "", "region_id"))) < 2L) "not_applicable" else if (!complete) "insufficient_evaluability" else
    if (observed == evaluated) "observed_all_regions" else if (observed == 1L) "region_specific_observation" else "observed_subset_regions"
  list(sample_count_observed = observed, evaluable_sample_count = evaluated,
    trajectory_state = state, distribution = if (observed == 1L) "observed_in_one_sample_only" else
      if (observed == length(cells)) "observed_in_all_samples" else "observed_in_subset",
    spatial_state = spatial, transitions = transitions)
}
longitudinal_hla_consistency <- function(contexts, subjects) {
  lapply(subjects, function(s) {
    declarations <- unlist(lapply(Filter(function(c) c$subject_id == s$subject_id, contexts), function(c)
      Filter(Negate(is.null), list(c$hla_genotype, c$observed_hla_genotype))), recursive = FALSE)
    genotypes <- lapply(declarations, function(g) sort(vapply(peptide_presenters(g), `[[`, "", "mhc_key")))
    list(subject_id = s$subject_id, status = if (length(genotypes) < 2L || isTRUE(s$transplant_context)) "not_evaluated" else
      if (length(unique(vapply(genotypes, function(g) infra_hash(as.list(g)), ""))) == 1L) "hla_genotype_consistent" else "hla_genotype_conflict",
      observed_genotypes = genotypes)
  })
}
longitudinal_escape <- function(cells, contexts, relations) {
  measurements <- lapply(cells, function(c) list(sample_id = c$sample_id,
    evidence = lapply(c$observations, function(o) list(candidate_id = o$candidate_id,
      hla_loh = o$evidence$hla_loh, presentation = o$evidence$presentation_evidence))))
  loh <- lapply(cells, function(c) {
    states <- unique(unlist(lapply(c$observations, function(o) {
      l <- o$evidence$hla_loh$candidate_link
      if (identical(l$presenter_hla_loh_status, "evaluated")) l$presenter_hla_loh_call else "not_evaluated"
    })))
    list(sample_id = c$sample_id, state = if (length(states) == 1L) states else "not_evaluated")
  })
  transitions <- list()
  apm_transitions <- list()
  if (length(loh) > 1L) for (i in seq_len(length(loh) - 1L)) for (j in seq.int(i + 1L, length(loh))) {
    a <- loh[[i]]; b <- loh[[j]]
    if (all(c(a$state, b$state) %in% c("loss_supported", "loss_not_supported")) &&
        !is.null(longitudinal_transition(a, b, contexts, relations)))
      transitions[[length(transitions) + 1L]] <- list(from = a$sample_id, to = b$sample_id,
        earlier_state = a$state, later_state = b$state,
        description = if (a$state == "loss_not_supported" && b$state == "loss_supported") "HLA loss newly supported in later sample" else "sample-specific HLA LOH states")
    if (!is.null(longitudinal_transition(a, b, contexts, relations))) {
      components <- function(cell) unlist(lapply(cell$observations, function(o)
        Filter(function(l) l$applicability == "applicable", o$evidence$presentation_evidence$antigen_presentation_machinery$links)), recursive = FALSE)
      before <- components(cells[[i]]); after <- components(cells[[j]])
      for (gene in intersect(vapply(before, `[[`, "", "gene"), vapply(after, `[[`, "", "gene"))) {
        earlier <- Filter(function(l) l$gene == gene, before); later <- Filter(function(l) l$gene == gene, after)
        apm_transitions[[length(apm_transitions) + 1L]] <- list(from = a$sample_id, to = b$sample_id, gene = gene,
          earlier_evidence = earlier, later_evidence = later,
          description = "Ordered component observations; expression or alteration alone does not establish acquired escape")
      }
    }
  }
  list(status = if (length(transitions)) "ordered_evidence_available" else "sample_specific_evidence_only",
    hla_loh_states = loh, hla_loh_transitions = transitions, apm_transitions = apm_transitions, sample_measurements = measurements,
    interpretation = "HLA LOH, gene/allele expression and component APM measurements remain separate; no universal escape score.")
}
