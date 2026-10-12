# T041: ordered objectives, never a weighted or summed biological score.
vaccine_covered <- function(selected) {
  atoms <- vaccine_atoms_empty()
  for (f in names(atoms)) atoms[[f]] <- vaccine_strings(unlist(lapply(selected, function(u) u$atoms[[f]])), f)
  atoms
}
vaccine_unmet <- function(selected, policy) {
  covered <- vaccine_covered(selected); out <- list()
  for (p in setdiff(unlist(policy$required_presenters), unlist(covered$presenters)))
    out <- c(out, list(list(constraint = "required_presenter", required = p, observed = FALSE)))
  for (cl in c("I", "II")) {
    n <- sum(vapply(selected, function(u) u$hla_class == cl, TRUE)); min <- policy[[if (cl == "I") "min_class_i" else "min_class_ii"]]
    if (n < min) out <- c(out, list(list(constraint = paste0("min_class_", tolower(cl)), required = min, observed = n)))
  }
  out
}
vaccine_admission_checks <- function(unit, selected, policy) {
  checks <- list()
  add <- function(k, count, limit, scope = NULL) checks[[length(checks) + 1L]] <<- list(
    constraint = k, status = if (count <= limit) "satisfied" else "failed",
    detail = list(count_after_admission = count, limit = limit, scope = scope))
  add("max_peptides", length(selected) + 1L, policy$max_peptides)
  max <- policy[[if (unit$hla_class == "I") "max_class_i" else "max_class_ii"]]
  if (!is.null(max)) add("class_maximum", sum(vapply(selected, function(u) u$hla_class == unit$hla_class, TRUE)) + 1L, max, unit$hla_class)
  if (!is.null(policy$max_per_presenter)) add("max_per_presenter", sum(vapply(selected, function(u) u$presenter == unit$presenter, TRUE)) + 1L,
    policy$max_per_presenter, unit$presenter)
  if (!is.null(policy$max_per_redundancy_group)) for (g in unit$strong_groups) {
    n <- sum(vapply(selected, function(u) g %in% unlist(u$strong_groups), TRUE))
    add("max_per_redundancy_group", n + 1L, policy$max_per_redundancy_group, g)
  }
  checks
}
vaccine_admission <- function(unit, selected, policy) {
  lapply(Filter(function(c) c$status == "failed", vaccine_admission_checks(unit, selected, policy)),
    function(c) vaccine_reason("set_constraint_failure", c$constraint, c$detail))
}
# Exact feasibility witness for declarative set constraints, not an objective optimizer.
# Branch on an unmet requirement; positive coverage fillers cannot repair deficits.
vaccine_feasibility <- function(units, selected, policy) {
  memo <- new.env(parent = emptyenv())
  search <- function(chosen) {
    unmet <- vaccine_unmet(chosen, policy)
    if (!length(unmet)) return(list())
    if (length(chosen) >= policy$max_peptides) return(NULL)
    key <- paste0("state:", paste(sort(vapply(chosen, `[[`, "", "candidate_key"), method = "radix"), collapse = ";"))
    if (exists(key, envir = memo, inherits = FALSE)) return(NULL)
    own <- vapply(chosen, `[[`, "", "candidate_key")
    remaining <- Filter(function(u) !u$candidate_key %in% own && !length(vaccine_admission(u, chosen, policy)), units)
    deficit <- unmet[[1L]]
    candidates <- Filter(function(u) if (deficit$constraint == "required_presenter") u$presenter == deficit$required else
      u$hla_class == if (deficit$constraint == "min_class_i") "I" else "II", remaining)
    for (u in candidates) {
      witness <- search(c(chosen, list(u)))
      if (!is.null(witness)) return(c(list(u$candidate_key), witness))
    }
    assign(key, TRUE, envir = memo)
    NULL
  }
  search(selected)
}
vaccine_marginal <- function(unit, selected, policy, relations) {
  covered <- vaccine_covered(selected)
  atoms <- lapply(names(covered), function(f) as.list(setdiff(unlist(unit$atoms[[f]]), unlist(covered[[f]]))))
  names(atoms) <- names(covered)
  required <- intersect(unlist(atoms$presenters), unlist(policy$required_presenters))
  min <- policy[[if (unit$hla_class == "I") "min_class_i" else "min_class_ii"]]
  class_gain <- as.integer(sum(vapply(selected, function(u) u$hla_class == unit$hla_class, TRUE)) < min)
  selected_keys <- vapply(selected, `[[`, "", "candidate_key")
  related <- if (!length(selected_keys) || !length(relations)) list() else Filter(function(r) unit$candidate_key %in% unlist(r$candidate_keys) && any(selected_keys %in% unlist(r$candidate_keys)), relations)
  strong <- length(unique(vapply(Filter(function(r) r$strong && r$type == "same_epitope_key", related), `[[`, "", "group_id")))
  overlap <- if (policy$reduce_overlap) length(unique(vapply(Filter(function(r) r$type == "source_region_overlap", related),
    function(r) paste(unlist(r$candidate_keys), collapse = "~"), ""))) else 0L
  # Every dimension remains separately named in both objective vector and audit.
  vector <- c(required_presenter_gain = length(required), class_quota_gain = class_gain,
    new_presenter_gain = length(atoms$presenters), new_context_gain = length(atoms$contexts),
    new_region_gain = length(atoms$regions), shared_gain = length(atoms$shared), persistent_gain = length(atoms$persistent),
    clonal_gain = length(atoms$clonal), supported_hla_loss = -length(unit$escape$hla_loh),
    supported_apm_loss = -length(unit$escape$apm_loss), observed_apm_abnormality = -length(unit$escape$apm_abnormality),
    strong_redundancy = -strong, policy_overlap = -overlap)
  list(vector = as.list(vector), newly_contributed = atoms, required_presenters = as.list(required),
    class_contribution = unit$hla_class, escape_evidence = unit$escape, redundancy_relations = related,
    has_incremental_value = any(vector[seq_len(8L)] > 0))
}
vaccine_choose <- function(considered, units) {
  alive <- seq_along(considered); decisive <- "stable_biological_key"
  for (f in names(considered[[1L]]$marginal$vector)) {
    values <- vapply(considered[alive], function(c) c$marginal$vector[[f]], 0)
    keep <- alive[values == max(values)]
    if (length(keep) < length(alive)) decisive <- f
    alive <- keep
    if (length(alive) == 1L) break
  }
  source <- if (length(alive) > 1L) "stable_biological_key" else "not_needed"; tied <- alive
  if (length(alive) > 1L) {
    obs <- lapply(units[alive], function(u) if (length(u$observations) == 1L) u$observations[[1L]] else NULL)
    comparable <- all(vapply(obs, function(o) !is.null(o) && !is.null(o$frozen_result$policy_version) &&
      !identical(o$frozen_result$final_grade, "NOT_GRADED") && is.numeric(o$frozen_result$rank_within_sample_class) &&
      length(o$frozen_result$rank_within_sample_class) == 1L && is.finite(o$frozen_result$rank_within_sample_class), TRUE))
    if (comparable) comparable <- length(unique(vapply(obs, function(o) paste(o$sample_id, o$hla_class, o$frozen_result$policy_version), ""))) == 1L
    if (comparable) {
      rank <- vapply(obs, function(o) o$frozen_result$rank_within_sample_class, 0)
      alive <- alive[rank == min(rank)]; source <- "frozen_rank_same_sample_class_policy"
    }
  }
  if (length(alive) > 1L && source != "stable_biological_key") source <- paste0(source, "+stable_biological_key")
  keys <- vapply(units[alive], `[[`, "", "candidate_key")
  chosen <- alive[order(keys, method = "radix")][[1L]]
  list(index = chosen, decisive_dimension = decisive, tie_break_source = source,
    tied_candidate_keys = as.list(sort(vapply(units[tied], `[[`, "", "candidate_key"), method = "radix")))
}
vaccine_select_subject <- function(projection, subject, binding) {
  policy <- binding$normalized; units <- vaccine_units(projection, subject)
  relations <- vaccine_redundancy(projection, subject)
  index_relations <- function(field) {
    keys <- lapply(relations, `[[`, field)
    indices <- split(rep(seq_along(relations), lengths(keys)), unlist(keys, use.names = FALSE))
    lapply(indices, function(i) relations[unique(i)])
  }
  by_key <- index_relations("candidate_keys"); by_id <- index_relations("candidate_ids")
  for (i in seq_along(units)) units[[i]]$relations <- by_key[[units[[i]]$candidate_key]]
  witness <- vaccine_feasibility(units, list(), policy)
  feasible <- !is.null(witness); selected <- list(); trace <- list(); stop_reason <- "no_incremental_coverage"
  repeat {
    selected_keys <- vapply(selected, `[[`, "", "candidate_key")
    remaining <- Filter(function(u) !u$candidate_key %in% selected_keys, units)
    if (!length(remaining)) { stop_reason <- "pool_exhausted"; break }
    considered <- lapply(remaining, function(u) {
      reasons <- vaccine_admission(u, selected, policy)
      if (!length(reasons) && feasible && length(vaccine_unmet(selected, policy)) &&
          is.null(vaccine_feasibility(units, c(selected, list(u)), policy)))
        reasons <- list(vaccine_reason("prevents_required_completion", "required_coverage"))
      list(candidate_key = u$candidate_key, marginal = vaccine_marginal(u, selected, policy, u$relations),
        set_constraint_checks = vaccine_admission_checks(u, selected, policy), constraint_failures = reasons, decision = "considered")
    })
    available <- which(vapply(considered, function(c) !length(c$constraint_failures) &&
      (c$marginal$has_incremental_value || policy$fill_to_capacity), TRUE))
    if (!length(available)) {
      trace <- c(trace, list(list(iteration = length(selected) + 1L, selected_candidate_key = NULL,
        considered = considered, tie_break = NULL, fill_to_capacity = policy$fill_to_capacity)))
      stop_reason <- if (length(selected) >= policy$max_peptides) "capacity_reached" else "no_incremental_coverage"
      break
    }
    choice <- vaccine_choose(considered[available], remaining[available]); idx <- available[[choice$index]]
    for (i in seq_along(considered)) considered[[i]]$decision <- if (i == idx) "selected" else if (length(considered[[i]]$constraint_failures))
      "constraint_blocked" else if (considered[[i]]$candidate_key %in% unlist(choice$tied_candidate_keys)) "lost_deterministic_tie_break" else "lower_lexicographic_contribution"
    trace <- c(trace, list(list(iteration = length(selected) + 1L, selected_candidate_key = remaining[[idx]]$candidate_key,
      considered = considered, tie_break = choice[setdiff(names(choice), "index")],
      fill_to_capacity = policy$fill_to_capacity && !considered[[idx]]$marginal$has_incremental_value)))
    selected <- c(selected, remaining[idx])
  }
  unmet <- vaccine_unmet(selected, policy)
  if (feasible && length(unmet)) stop("Selector failed a proven feasible policy")
  own <- Filter(function(p) p$subject_id == subject, projection$candidates)
  selected_keys <- vapply(selected, `[[`, "", "candidate_key")
  audit <- lapply(own, function(p) {
    in_scope <- any(vapply(projection$contexts, function(c) c$subject_id == p$subject_id && c$sample_id == p$sample_id, TRUE))
    supported <- vaccine_observation_supported(p)
    picked <- in_scope && supported && p$candidate_key %in% selected_keys && p$eligibility$status == "eligible"
    history <- unlist(lapply(trace, function(step) lapply(Filter(function(c) c$candidate_key == p$candidate_key, step$considered),
      function(c) c(list(iteration = step$iteration, tie_break = step$tie_break), c))), recursive = FALSE)
    if (!length(history)) history <- list()
    reasons <- p$eligibility$reasons
    if (!in_scope) reasons <- c(list(vaccine_reason("outside_tumor_scope", "subject_context_scope")), reasons)
    else if (!supported) reasons <- c(list(vaccine_reason("authoritative_observation_not_supported", "t039_observation_support",
      list(cells = Filter(function(cell) identical(cell$sample_id, p$sample_id), p$trajectory$cells),
        interpretation = "Unsupported selection evidence; not negative biology"))), reasons)
    else if (picked) reasons <- list(vaccine_reason("selected", "lexicographic_marginal_contribution",
      Filter(function(h) h$decision == "selected", history)))
    else if (!length(reasons)) {
      u <- Filter(function(u) u$candidate_key == p$candidate_key, units)
      if (!length(u)) reasons <- list(vaccine_reason("outside_tumor_scope", "subject_context_scope")) else {
        failures <- vaccine_admission(u[[1L]], selected, policy)
        m <- vaccine_marginal(u[[1L]], selected, policy, u[[1L]]$relations)
        redundant <- any(vapply(m$redundancy_relations, function(r) r$strong ||
          (policy$reduce_overlap && r$type == "source_region_overlap"), TRUE))
        reasons <- if (length(failures)) failures else list(vaccine_reason(
          if (any(vapply(history, function(h) h$decision == "lost_deterministic_tie_break", TRUE))) "lost_deterministic_tie_break" else
          if (redundant && !m$has_incremental_value) "redundant_with_selected_coverage" else "no_incremental_coverage", stop_reason))
      }
    }
    list(candidate_id = p$candidate_id, candidate_key = p$candidate_key, sample_id = p$sample_id,
      status = if (picked) "selected" else "unselected", eligibility = p$eligibility, reasons = reasons,
      projection = p, selection_history = history,
      redundancy_relations = if (is.null(by_id[[p$candidate_id]])) list() else by_id[[p$candidate_id]])
  })
  scope <- list(subject_id = subject, sample_ids = as.list(sort(vapply(Filter(function(c) c$subject_id == subject, projection$contexts), `[[`, "", "sample_id"), method = "radix")))
  semantic <- list(scope = scope, policy_id = binding$policy_id, algorithm = vaccine_algorithm(), selected_candidate_keys = as.list(sort(selected_keys, method = "radix")))
  list(vaccine_set_id = vaccine_id("set", semantic), scope = scope,
    status = if (length(unmet)) "infeasible" else "feasible", selected_count = length(selected),
    max_peptides = policy$max_peptides, class_i_count = sum(vapply(selected, function(u) u$hla_class == "I", TRUE)),
    class_ii_count = sum(vapply(selected, function(u) u$hla_class == "II", TRUE)),
    selected_candidates = lapply(selected, function(u) u[c("candidate_key", "candidate_ids", "ligand_sequence", "presenter", "hla_class")]),
    candidate_audit = audit, coverage = vaccine_covered(selected),
    evaluability = lapply(own, function(p) list(candidate_id = p$candidate_id,
      trajectory = p$trajectory, single_sample_observation = if (is.null(p$trajectory)) p$sample_id else NULL)),
    unmet_constraints = unmet, redundancy_relations = relations, selection_trace = trace,
    feasibility = list(method = "exact-declarative-constraint-search-v1", initial_feasible = feasible,
      witness_candidate_keys = witness), stop_reason = stop_reason)
}
vaccine_construct <- function(projection, binding) {
  subjects <- sort(unique(c(unlist(projection$subject_ids), vapply(projection$contexts, `[[`, "", "subject_id"),
    vapply(projection$candidates, `[[`, "", "subject_id"))), method = "radix")
  list(schema_version = vaccine_selection_version(), policy_version = vaccine_policy_version(),
    policy_id = binding$policy_id, normalized_policy = binding$normalized,
    algorithm = list(version = vaccine_algorithm(), kind = "deterministic_lexicographic_heuristic", global_optimality_claim = FALSE),
    identity_namespaces = vaccine_namespaces(), selections = lapply(subjects, function(s) vaccine_select_subject(projection, s, binding)),
    provenance = list(semantic = list(policy_id = binding$policy_id, normalized_policy = binding$normalized,
      algorithm = vaccine_algorithm()), policy_input = binding$input, physical = binding$physical, raw_policy_text = binding$raw_text),
    interpretation = "Selection overlay only. Missing is not negative. No response probability, causal claim or universal immunogenicity score.")
}
