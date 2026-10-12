test_that("policy requires opt-in explicit capacity and rejects hidden constraints", {
  expect_null(neo_config(env = character())$vaccine_set)
  expect_error(Neoantimon2:::vaccine_policy_bind(list(version = "vaccine-set-policy-1.0.0")), "max_peptides")
  for (v in list(-1, 1.5, NA_real_, Inf, "4", TRUE)) expect_error(vaccine_binding(max_class_i = v), "count")
  expect_error(vaccine_binding(immunogenicity_score = .4), "Unknown")
  expect_error(vaccine_binding(fill_to_capacity = 1), "boolean")
  expect_error(vaccine_binding(forbidden_motifs = list("A.*")), "literal")
  expect_error(vaccine_binding(required_presenters = list("A0201")), "canonical")
  expect_identical(vaccine_binding()$normalized$max_peptides, 4L)
})

test_that("semantic policy identity ignores order paths bytes formats and defaults", {
  x <- vaccine_fixture("policy.json"); a <- Neoantimon2:::vaccine_policy_bind(x)
  y <- x[rev(names(x))]; y$required_presenters <- rev(y$required_presenters)
  expect_identical(a$policy_id, Neoantimon2:::vaccine_policy_bind(y)$policy_id)
  f <- tempfile(fileext = ".json"); g <- tempfile(fileext = ".yaml")
  jsonlite::write_json(x, f, pretty = TRUE, auto_unbox = TRUE); yaml::write_yaml(y, g)
  b <- Neoantimon2:::vaccine_policy_bind(f); c <- Neoantimon2:::vaccine_policy_bind(g)
  expect_identical(a$policy_id, b$policy_id); expect_identical(b$policy_id, c$policy_id)
  expect_false(identical(b$physical$sha256, c$physical$sha256))
  expect_false(identical(b$physical$path, c$physical$path))
  expect_identical(b$raw_text, paste0(paste(readLines(f), collapse = "\n"), "\n"))
  y$fill_to_capacity <- FALSE; y$min_class_i <- 2.0
  expect_identical(a$policy_id, Neoantimon2:::vaccine_policy_bind(y)$policy_id)
  y$max_peptides <- 5L
  expect_false(identical(a$policy_id, Neoantimon2:::vaccine_policy_bind(y)$policy_id))
})

test_that("capacity is an upper bound and fill requires explicit policy", {
  f <- vaccine_fixture(); f$candidates <- f$candidates[c(2L, 10L, 11L)]
  b <- vaccine_binding(); s <- vaccine_fixture_select(b, f)$selections[[1L]]
  expect_identical(s$selected_count, 1L); expect_lt(s$selected_count, length(f$candidates))
  expect_identical(s$stop_reason, "no_incremental_coverage")
  b <- vaccine_binding(fill_to_capacity = TRUE, max_class_i = 2L)
  s <- vaccine_fixture_select(b, f)$selections[[1L]]
  expect_identical(s$selected_count, 2L); expect_true(any(vapply(s$selection_trace, `[[`, TRUE, "fill_to_capacity")))
  for (limit in 0:5) {
    b <- Neoantimon2:::vaccine_policy_bind(list(version = "vaccine-set-policy-1.0.0", max_peptides = limit, fill_to_capacity = TRUE))
    s <- vaccine_fixture_select(b)$selections[[1L]]
    expect_lte(s$selected_count, limit)
    expect_true(all(vapply(s$candidate_audit, function(a) length(a$reasons) > 0L, TRUE)))
  }
})

test_that("required presenters and heterodimer coverage precede optional coverage", {
  b <- Neoantimon2:::vaccine_policy_bind(vaccine_fixture("policy.json")); s <- vaccine_fixture_select(b)$selections[[1L]]
  expect_identical(s$status, "feasible"); expect_length(s$unmet_constraints, 0L)
  expect_true(all(unlist(b$normalized$required_presenters) %in% unlist(s$coverage$presenters)))
  expect_gte(s$class_i_count, 2); expect_gte(s$class_ii_count, 1)
  expect_lte(s$selected_count, 4)
  expect_true(any(vapply(s$selection_trace, function(step) any(vapply(step$considered, function(c) c$marginal$vector$required_presenter_gain > 0L, TRUE)), TRUE)))
})

test_that("impossible presenter and class constraints produce structured infeasibility", {
  policies <- list(vaccine_binding(required_presenters = list("HLA-A*01:01")),
    vaccine_binding(min_class_ii = 2L), vaccine_binding(min_class_i = 4L, min_class_ii = 1L),
    vaccine_binding(min_class_i = 2L, max_class_i = 1L),
    vaccine_binding(required_presenters = list("HLA-A*02:01"), max_per_presenter = 0L))
  for (b in policies) {
    s <- vaccine_fixture_select(b)$selections[[1L]]
    expect_identical(s$status, "infeasible"); expect_gt(length(s$unmet_constraints), 0L)
    expect_false(s$feasibility$initial_feasible); expect_lte(s$selected_count, 4L)
    expect_true(all(vapply(s$unmet_constraints, function(u) all(c("constraint", "required", "observed") %in% names(u)), TRUE)))
  }
})

test_that("class minimum maximum and presenter maximum constraints remain separate", {
  for (b in list(vaccine_binding(min_class_i = 2L), vaccine_binding(min_class_ii = 1L),
      vaccine_binding(min_class_i = 2L, min_class_ii = 1L, max_class_i = 2L, max_class_ii = 1L, max_per_presenter = 1L))) {
    s <- vaccine_fixture_select(b)$selections[[1L]]
    expect_identical(s$status, "feasible")
    expect_gte(s$class_i_count, b$normalized$min_class_i); expect_gte(s$class_ii_count, b$normalized$min_class_ii)
    if (!is.null(b$normalized$max_class_i)) expect_lte(s$class_i_count, b$normalized$max_class_i)
    if (!is.null(b$normalized$max_per_presenter)) {
      expect_true(all(table(unlist(s$coverage$presenters)) <= 1L))
      expect_true(all(table(vapply(s$selected_candidates, `[[`, "", "presenter")) <= b$normalized$max_per_presenter))
    }
  }
})

test_that("constraint feasibility guards a greedy choice that consumes a required slot", {
  f <- vaccine_fixture(); f$candidates <- f$candidates[c(1L, 3L, 4L)]
  b <- Neoantimon2:::vaccine_policy_bind(list(version = "vaccine-set-policy-1.0.0", max_peptides = 2L,
    required_presenters = list("HLA-A*02:01"), min_class_ii = 1L))
  s <- vaccine_fixture_select(b, f)$selections[[1L]]
  expect_identical(s$status, "feasible"); expect_identical(s$selected_count, 2L)
  expect_identical(s$class_ii_count, 1L)
  expect_true(any(vapply(s$selection_trace, function(step) any(vapply(step$considered, function(c)
    any(vapply(c$constraint_failures, function(r) r$code == "prevents_required_completion", TRUE)), TRUE)), TRUE)))
})

test_that("candidate and context input ordering cannot alter set identity or trace", {
  b <- Neoantimon2:::vaccine_policy_bind(vaccine_fixture("policy.json")); f <- vaccine_fixture()
  records <- vaccine_records(f); a <- Neoantimon2:::vaccine_construct(Neoantimon2:::vaccine_project_records(records, f$contexts, list(), b$normalized), b)
  for (i in 1:4) {
    order <- if (i == 1L) rev(seq_along(records)) else c(seq.int(i, length(records)), seq_len(i - 1L))
    r <- lapply(records[order], function(r) { r$trajectory$cells <- rev(r$trajectory$cells); r })
    z <- Neoantimon2:::vaccine_construct(Neoantimon2:::vaccine_project_records(r, rev(f$contexts), list(), b$normalized), b)
    expect_identical(vaccine_keys(a), vaccine_keys(z))
    expect_identical(a$selections[[1L]]$vaccine_set_id, z$selections[[1L]]$vaccine_set_id)
    expect_identical(a$selections[[1L]]$selection_trace, z$selections[[1L]]$selection_trace)
  }
  expect_identical(a, Neoantimon2:::vaccine_construct(Neoantimon2:::vaccine_project_records(records, f$contexts, list(), b$normalized), b))
})

test_that("set identity responds only to semantic policy scope and chosen frozen keys", {
  a <- vaccine_fixture_select(); b <- vaccine_fixture_select(vaccine_binding(blocklist = list(vaccine_keys(a)[[1L]])))
  expect_false(identical(a$selections[[1L]]$vaccine_set_id, b$selections[[1L]]$vaccine_set_id))
  f <- vaccine_fixture(); f$contexts[[1L]]$subject_id <- "other-subject"; f$subject_id <- "other-subject"
  for (i in seq_along(f$contexts)) f$contexts[[i]]$subject_id <- "other-subject"
  expect_false(identical(a$selections[[1L]]$vaccine_set_id, vaccine_fixture_select(fixture = f)$selections[[1L]]$vaccine_set_id))
  binding <- vaccine_binding(); binding$physical <- list(path = "/different/file", sha256 = "different-sha")
  b <- Neoantimon2:::vaccine_construct(vaccine_fixture_projection(binding), binding)
  expect_identical(a$selections[[1L]]$vaccine_set_id, b$selections[[1L]]$vaccine_set_id)
  expect_identical(a$selections[[1L]]$selection_trace, b$selections[[1L]]$selection_trace)
})

test_that("frozen INVALID is excluded while NOT_GRADED is independently eligible", {
  p <- vaccine_fixture_projection()
  invalid <- Filter(function(p) p$frozen_result$final_grade == "INVALID", p$candidates)
  ungraded <- Filter(function(p) p$frozen_result$final_grade == "NOT_GRADED", p$candidates)
  expect_gt(length(invalid), 0); expect_gt(length(ungraded), 0)
  for (p in invalid) { expect_identical(p$eligibility$status, "ineligible"); expect_identical(p$frozen_result$final_grade, "INVALID") }
  for (p in ungraded) { expect_identical(p$eligibility$status, "eligible"); expect_identical(p$frozen_result$final_grade, "NOT_GRADED") }
  s <- vaccine_fixture_select(vaccine_binding(min_class_ii = 1L))$selections[[1L]]
  expect_true(any(vapply(s$candidate_audit, function(a) a$status == "selected" && a$projection$frozen_result$final_grade == "NOT_GRADED", TRUE)))
  expect_false(any(vapply(s$candidate_audit, function(a) a$status == "selected" && a$projection$frozen_result$final_grade == "INVALID", TRUE)))
})

test_that("NOT_GRADED fusion splice and long-read sources never acquire a synthetic grade", {
  s <- vaccine_fixture_select(vaccine_binding(fill_to_capacity = TRUE))$selections[[1L]]
  for (source in c("fusion", "splice", "long_read")) {
    audits <- Filter(function(a) any(vapply(a$projection$origins, function(o) identical(o$source_type, source), TRUE)), s$candidate_audit)
    expect_gt(length(audits), 0L)
    for (a in audits) {
      expect_identical(a$projection$frozen_result$final_grade, "NOT_GRADED")
      expect_identical(a$projection$frozen_result$policy_version, "source-not-graded-v1")
      expect_null(a$projection$frozen_result$rank_within_sample_class)
    }
  }
})

test_that("missing insufficient and evaluable nonobservation remain distinct", {
  p <- vaccine_fixture_projection()
  candidate <- Filter(function(p) p$presenter == "HLA-C*07:02", p$candidates)[[1L]]
  expect_identical(vapply(candidate$trajectory$cells, `[[`, "", "state"), c("observed", "not_observed_insufficient", "not_evaluated"))
  expect_identical(candidate$trajectory$summary$trajectory_state, "insufficient_to_classify")
  expect_null(candidate$frozen_result$self_exact_match); expect_null(candidate$evidence$hla_loh)
  expect_length(candidate$escape$hla_loh, 0L); expect_length(candidate$escape$apm_loss, 0L)
  u <- Filter(function(u) u$presenter == "HLA-C*07:02", Neoantimon2:::vaccine_units(p, "patient-1"))[[1L]]
  expect_identical(u$atoms$contexts, list("S1")); expect_length(u$atoms$persistent, 0L); expect_length(u$atoms$shared, 0L)
  lost <- Filter(function(p) p$presenter == "HLA-A*24:02", p$candidates)[[1L]]
  expect_identical(lost$trajectory$cells[[3L]]$state, "not_observed_evaluable")
  expect_identical(lost$trajectory$summary$trajectory_state, "loss_with_supported_absence")
  gained <- Filter(function(p) p$presenter == "HLA-B*08:01", p$candidates)[[1L]]
  expect_identical(gained$trajectory$summary$trajectory_state, "gain_with_prior_supported_absence")
})

test_that("unordered shared regions have spatial coverage without ordered persistence", {
  f <- vaccine_fixture(); for (i in seq_along(f$contexts)) f$contexts[[i]]$ordinal_time <- NULL
  f$candidates <- f$candidates[1L]
  p <- vaccine_fixture_projection(fixture = f); u <- Neoantimon2:::vaccine_units(p, "patient-1")[[1L]]
  expect_identical(p$candidates[[1L]]$trajectory$ordering, "unordered")
  expect_length(p$candidates[[1L]]$trajectory$summary$transitions, 0L)
  expect_length(u$atoms$shared, 3L); expect_length(u$atoms$persistent, 0L)
  expect_false(identical(Neoantimon2:::vaccine_pair_atom(c("a~b", "c")), Neoantimon2:::vaccine_pair_atom(c("a", "b~c"))))
  expect_identical(jsonlite::fromJSON(Neoantimon2:::vaccine_pair_atom(c("a~b", "c"))), c("a~b", "c"))
})

test_that("technical replicate observations cannot inflate longitudinal coverage", {
  f <- vaccine_fixture(); f$candidates <- f$candidates[1L]
  f$relations <- list(list(from = "S1", to = "S2", type = "technical_replicate"), list(from = "S2", to = "S3", type = "technical_replicate"))
  p <- vaccine_fixture_projection(fixture = f); u <- Neoantimon2:::vaccine_units(p, "patient-1")[[1L]]
  expect_length(u$atoms$contexts, 1L); expect_length(u$atoms$regions, 1L)
  expect_length(u$atoms$shared, 0L); expect_length(u$atoms$persistent, 0L)
  expect_length(p$candidates, 3L)
  bridge <- list(list(from = "S1", to = "normal-bridge", type = "technical_replicate"),
    list(from = "normal-bridge", to = "S2", type = "technical_replicate"), list(from = "S2", to = "S3", type = "technical_replicate"))
  expect_length(unique(Neoantimon2:::vaccine_context_units(f$contexts, bridge)), 1L)
  f$relations <- list(); for (i in seq_along(f$contexts)) f$contexts[[i]]$biological_context_id <- "biopsy:1"
  u <- Neoantimon2:::vaccine_units(vaccine_fixture_projection(fixture = f), "patient-1")[[1L]]
  expect_length(u$atoms$contexts, 1L); expect_length(u$atoms$persistent, 0L)
})

test_that("linked clonal quantities are retained separately and multiple origins survive", {
  f <- vaccine_fixture(); records <- vaccine_records(f)
  records[[1L]]$evidence$origins <- c(records[[1L]]$evidence$origins, list(list(origin_id = "origin:extra", source_event_ids = list("E9"), origin_validity = "VALID")))
  p <- Neoantimon2:::vaccine_project_records(records, f$contexts, list(), vaccine_binding()$normalized)
  a <- Filter(function(p) length(p$origins) == 2L, p$candidates)[[1L]]
  expect_length(a$origins, 2L); expect_identical(a$frozen_result$ccf_value, .31)
  expect_identical(a$clonal_evidence[[1L]]$cellular_prevalence, .72)
  expect_identical(a$clonal_evidence[[1L]]$cellular_prevalence_std, .08)
  expect_false(identical(a$frozen_result$ccf_value, a$clonal_evidence[[1L]]$cellular_prevalence))
  b <- Filter(function(p) p$presenter == "HLA-C*07:02", p$candidates)[[1L]]
  expect_length(b$clonal_evidence, 0L)
  units <- Neoantimon2:::vaccine_units(p, "patient-1")
  before <- lapply(units, function(u) u$atoms$clonal)
  for (i in seq_along(records)) for (j in seq_along(records[[i]]$trajectory$cells))
    for (k in seq_along(records[[i]]$trajectory$cells[[j]]$clonal_evidence)) records[[i]]$trajectory$cells[[j]]$clonal_evidence[[k]]$cellular_prevalence <- .01
  after <- Neoantimon2:::vaccine_units(Neoantimon2:::vaccine_project_records(records, f$contexts, list(), vaccine_binding()$normalized), "patient-1")
  expect_identical(before, lapply(after, function(u) u$atoms$clonal))
})

test_that("allele loss remains allele-specific and B2M remains Class I-specific", {
  f <- vaccine_fixture(); f$candidates[[4L]]$apm_samples <- list("S3")
  p <- vaccine_fixture_projection(fixture = f)
  a <- Filter(function(p) p$presenter == "HLA-A*02:01" && p$sample_id == "S3" && p$ligand_sequence == "GILGFVFTL", p$candidates)[[1L]]
  b <- Filter(function(p) p$presenter == "HLA-B*07:02" && p$sample_id == "S3", p$candidates)[[1L]]
  ii <- Filter(function(p) p$hla_class == "II" && p$sample_id == "S3", p$candidates)[[1L]]
  expect_length(a$escape$hla_loh, 1L); expect_length(b$escape$hla_loh, 0L)
  expect_length(a$escape$apm_loss, 1L); expect_length(b$escape$apm_loss, 1L)
  expect_length(ii$escape$hla_loh, 0L); expect_length(ii$escape$apm_loss, 0L)
  expect_length(ii$evidence$presentation_evidence$antigen_presentation_machinery$observations, 1L)
  expect_identical(a$evidence$hla_loh$candidate_link$presenter_hla_loh_call, "loss_supported")
  expect_identical(a$escape$hla_loh[[1L]]$presenter, "HLA-A*02:01")
})

test_that("escape changes a tied marginal decision without any intact claim", {
  f <- vaccine_fixture(); f$candidates <- f$candidates[1:2]
  b <- vaccine_binding(); s <- vaccine_fixture_select(b, f)$selections[[1L]]
  expect_identical(s$selected_candidates[[1L]]$ligand_sequence, "ILGFVFTLA")
  selected <- Filter(function(c) c$decision == "selected", s$selection_trace[[1L]]$considered)[[1L]]
  other <- Filter(function(c) c$decision != "selected", s$selection_trace[[1L]]$considered)[[1L]]
  expect_identical(selected$marginal$vector$supported_hla_loss, 0L)
  expect_identical(other$marginal$vector$supported_hla_loss, -1L)
  expect_length(other$marginal$escape_evidence$hla_loh, 1L)
  expect_false(grepl("intact|no_escape", jsonlite::toJSON(Neoantimon2:::longitudinal_plain(selected))))
})

test_that("redundancy preserves duplicate identities presenters epitopes overlap and events", {
  f <- vaccine_fixture(); records <- vaccine_records(f); b <- vaccine_binding()
  a <- Neoantimon2:::vaccine_project_records(records, f$contexts, list(), b$normalized)
  d <- Neoantimon2:::vaccine_project_records(c(records, records[1L]), f$contexts, list(), b$normalized)
  expect_identical(a, d)
  changed <- records[[1L]]; changed$result$final_grade <- "A"
  expect_error(Neoantimon2:::vaccine_project_records(c(records, list(changed)), f$contexts, list(), b$normalized), "Conflicting duplicate")
  relations <- Neoantimon2:::vaccine_redundancy(a, "patient-1")
  for (type in c("same_frozen_candidate_key", "same_epitope_key", "source_region_overlap", "same_ligand_different_presenter", "same_source_event"))
    expect_true(any(vapply(relations, function(r) r$type == type, TRUE)))
  expect_false(any(vapply(Filter(function(r) r$type == "same_ligand_different_presenter", relations), `[[`, TRUE, "strong")))
  expect_false(any(vapply(Filter(function(r) r$type == "source_region_overlap", relations), `[[`, TRUE, "strong")))
  expect_identical(length(a$candidates), length(records))
  # Coordinates from distinct protein/ORF axes do not establish overlap.
  for (i in seq_along(records)) records[[i]]$evidence$origins[[1L]]$sequence_context_id <- records[[i]]$candidate_key
  axes <- Neoantimon2:::vaccine_redundancy(Neoantimon2:::vaccine_project_records(records, f$contexts, list(), b$normalized), "patient-1")
  expect_false(any(vapply(axes, function(r) r$type == "source_region_overlap" && length(r$candidate_keys) > 1L, TRUE)))
})

test_that("redundancy limits and optional overlap affect only explicit rules", {
  f <- vaccine_fixture(); f$candidates <- f$candidates[1:2]
  b <- vaccine_binding(fill_to_capacity = TRUE, max_per_redundancy_group = 1L)
  s <- vaccine_fixture_select(b, f)$selections[[1L]]
  expect_identical(s$selected_count, 1L)
  expect_true(any(vapply(s$candidate_audit, function(a) any(vapply(a$reasons, function(r) r$constraint == "max_per_redundancy_group", TRUE)), TRUE)))
  p <- vaccine_fixture_projection(fixture = f); units <- Neoantimon2:::vaccine_units(p, "patient-1")
  rel <- Neoantimon2:::vaccine_redundancy(p, "patient-1")
  m <- Neoantimon2:::vaccine_marginal(units[[2L]], units[1L], vaccine_binding()$normalized, rel)
  expect_identical(m$vector$policy_overlap, 0L)
  m <- Neoantimon2:::vaccine_marginal(units[[2L]], units[1L], vaccine_binding(reduce_overlap = TRUE)$normalized, rel)
  expect_lt(m$vector$policy_overlap, 0L)
  expect_gt(length(m$redundancy_relations), 0L)
})

test_that("frozen rank is used only for a comparable final tie", {
  b <- vaccine_binding(); p <- vaccine_fixture_projection(b)
  units <- Neoantimon2:::vaccine_units(p, "patient-1")[1:2]
  for (i in seq_along(units)) {
    units[[i]]$observations <- units[[i]]$observations[1L]
    units[[i]]$observations[[1L]]$frozen_result$rank_within_sample_class <- as.integer(3L - i)
    units[[i]]$observations[[1L]]$frozen_result$final_grade <- "B"
    units[[i]]$observations[[1L]]$frozen_result$policy_version <- "somatic-priority-v1"
    units[[i]]$observations[[1L]]$hla_class <- "I"
    units[[i]]$observations[[1L]]$sample_id <- "S1"
  }
  considered <- lapply(units, function(u) list(marginal = list(vector = list(new_context_gain = 1L))))
  chosen <- Neoantimon2:::vaccine_choose(considered, units)
  expect_identical(chosen$index, 2L); expect_identical(chosen$tie_break_source, "frozen_rank_same_sample_class_policy")
  units[[2L]]$observations[[1L]]$frozen_result$final_grade <- "NOT_GRADED"
  chosen <- Neoantimon2:::vaccine_choose(considered, units)
  expect_identical(chosen$tie_break_source, "stable_biological_key")
  expect_identical(units[[chosen$index]]$candidate_key, sort(vapply(units, `[[`, "", "candidate_key"))[[1L]])
})

test_that("declarative length residues motifs allowlists and blocklists are auditable", {
  b <- vaccine_binding(min_length = 9L, max_length = 10L, forbidden_residues = list("P"), forbidden_motifs = list("RRR"), forbid_confirmed_exact_self = TRUE)
  p <- vaccine_fixture_projection(b)
  for (candidate in p$candidates) {
    if (nchar(candidate$ligand_sequence) > 10L || grepl("P|RRR", candidate$ligand_sequence) || isTRUE(candidate$frozen_result$self_exact_match))
      expect_identical(candidate$eligibility$status, "ineligible")
    expect_gt(length(candidate$eligibility$checks), 0L)
  }
  key <- p$candidates[[1L]]$candidate_key
  a <- vaccine_fixture_projection(vaccine_binding(allowlist = list(key)))
  expect_true(all(vapply(Filter(function(p) p$eligibility$status == "eligible", a$candidates), function(p) p$candidate_key == key, TRUE)))
  a <- vaccine_fixture_projection(vaccine_binding(blocklist = list(key)))
  expect_true(all(vapply(Filter(function(p) p$candidate_key == key, a$candidates), function(p) p$eligibility$status == "ineligible", TRUE)))
  no <- vaccine_fixture_projection(); expect_true(all(vapply(no$candidates, function(p) identical(vapply(p$eligibility$checks, `[[`, "", "constraint"), "frozen_technical_validity"), TRUE)))
})

test_that("explicit evidence requirements report insufficient evidence independently", {
  for (req in c("exact_self_evaluated", "longitudinal_evaluable", "linked_clonal", "ms_observed")) {
    p <- vaccine_fixture_projection(vaccine_binding(require_evidence = list(req)))
    not_met <- Filter(function(p) p$eligibility$status == "ineligible" && p$frozen_result$final_grade != "INVALID", p$candidates)
    expect_gt(length(not_met), 0L)
    expect_true(all(vapply(not_met, function(p) any(vapply(p$eligibility$reasons, function(r) r$code == "insufficient_evidence_for_requirement", TRUE)), TRUE)))
    expect_true(all(vapply(not_met, function(p) p$frozen_result$final_grade %in% c("B", "NOT_GRADED"), TRUE)))
  }
})

test_that("MS and receptor observations stay sample-specific without negative conversion", {
  p <- vaccine_fixture_projection(); observed <- p$candidates[[1L]]
  for (a in p$candidates) {
    ms <- a$evidence$presentation_evidence$immunopeptidomics$observations
    expect_true(all(vapply(ms, function(m) m$context$sample_id == a$sample_id, TRUE)))
    expect_identical(a$evidence$receptor_evidence$status, "no_observation")
    if (a$sample_id != "S1") expect_length(ms, 0L)
  }
})

test_that("every selected and unselected observation has independent structured rationale", {
  s <- vaccine_fixture_select(Neoantimon2:::vaccine_policy_bind(vaccine_fixture("policy.json")))$selections[[1L]]
  expect_true(any(vapply(s$candidate_audit, function(a) a$status == "selected", TRUE)))
  expect_true(any(vapply(s$candidate_audit, function(a) a$status == "unselected", TRUE)))
  for (a in s$candidate_audit) {
    expect_gt(length(a$reasons), 0L); expect_true(all(c("frozen_result", "evidence", "eligibility", "trajectory", "origins") %in% names(a$projection)))
    if (a$status == "selected") expect_true(any(vapply(a$selection_history, function(h) h$decision == "selected", TRUE)))
  }
  for (step in s$selection_trace) for (c in step$considered) {
    expect_true(any(vapply(c$set_constraint_checks, function(check) check$constraint == "max_peptides", TRUE)))
    if (c$decision == "selected") expect_true(all(vapply(c$set_constraint_checks, function(check) check$status == "satisfied", TRUE)))
    expect_true(all(c("required_presenter_gain", "class_quota_gain", "new_context_gain", "shared_gain", "persistent_gain", "clonal_gain", "supported_hla_loss", "supported_apm_loss", "strong_redundancy") %in% names(c$marginal$vector)))
    expect_false(any(grepl("score|probability|fitness", names(c$marginal$vector))))
  }
})

test_that("no universal score is present in runtime schema or results", {
  x <- vaccine_fixture_select()
  names_recursive <- function(x) if (is.list(x)) c(names(x), unlist(lapply(x, names_recursive))) else character()
  expect_false(any(grepl("immunogenicity_score|response_probability|presentation_fitness|combined_score|weighted_score", names_recursive(x))))
  for (f in c("vaccine-set-policy-1.0.0.json", "vaccine-set-selection-1.0.0.json")) {
    schema <- jsonlite::read_json(system.file("schema", f, package = "Neoantimon2"))
    expect_false(any(grepl("score|probability|fitness", names_recursive(schema))))
  }
  expect_identical(x$algorithm$kind, "deterministic_lexicographic_heuristic")
  expect_false(x$algorithm$global_optimality_claim)
})

test_that("empty pools zero capacity normal contexts and subject isolation are explicit", {
  f <- vaccine_fixture(); f$candidates <- list()
  s <- vaccine_fixture_select(fixture = f)$selections[[1L]]
  expect_identical(s$selected_count, 0L); expect_identical(s$status, "feasible")
  s <- vaccine_fixture_select(vaccine_binding(min_class_i = 1L), f)$selections[[1L]]
  expect_identical(s$status, "infeasible")
  f <- vaccine_fixture(); f$contexts[[1L]]$role <- "normal"
  p <- vaccine_fixture_projection(fixture = f)
  expect_false("S1" %in% unlist(Neoantimon2:::vaccine_units(p, "patient-1")[[1L]]$atoms$contexts))
  for (i in seq_along(f$contexts)) f$contexts[[i]]$role <- "normal"
  normal <- vaccine_fixture_select(vaccine_binding(min_class_i = 1L), f)$selections[[1L]]
  expect_identical(normal$status, "infeasible"); expect_identical(normal$selected_count, 0L)
  expect_true(all(vapply(normal$candidate_audit, function(a) a$status == "unselected" && length(a$reasons) > 0L, TRUE)))
  f <- vaccine_fixture(); records <- vaccine_records(f)
  other <- records; for (i in seq_along(other)) { other[[i]]$subject_id <- "patient-2"; other[[i]]$result$candidate_id <- paste0(other[[i]]$result$candidate_id, ":P2") }
  ctx <- lapply(f$contexts, function(c) { c$subject_id <- "patient-2"; c })
  x <- Neoantimon2:::vaccine_construct(Neoantimon2:::vaccine_project_records(c(records, other), c(f$contexts, ctx), list(), vaccine_binding()$normalized), vaccine_binding())
  expect_length(x$selections, 2L); expect_false(identical(x$selections[[1L]]$vaccine_set_id, x$selections[[2L]]$vaccine_set_id))
  expect_true(all(vapply(x$selections, function(s) all(vapply(s$candidate_audit, function(a) a$projection$subject_id == s$scope$subject_id, TRUE)), TRUE)))
})

test_that("real multi-context vaccine E2E preserves all frozen results and trajectories", {
  e <- vaccine_real_e2e(); a <- neo_results(e$before); b <- neo_results(e$after)
  summary <- attr(b, "vaccine_set"); attr(b, "vaccine_set") <- NULL
  expect_identical(a, b); expect_identical(summary$algorithm$version, "lexicographic-greedy-v1")
  expect_identical(Neoantimon2:::longitudinal_read(e$before), Neoantimon2:::longitudinal_read(e$after))
  expect_identical(neo_events(e$before), neo_events(e$after)); expect_identical(neo_origins(e$before), neo_origins(e$after))
  expect_identical(neo_multisample(e$before), neo_multisample(e$after))
  s <- e$selection$selections[[1L]]; expect_identical(s$status, "feasible")
  expect_lte(s$selected_count, 4L); expect_gte(s$class_i_count, 2L); expect_gte(s$class_ii_count, 1L)
  expect_length(s$coverage$presenters, 3L); expect_length(s$coverage$contexts, 3L)
  expect_gt(length(s$coverage$persistent), 0L); expect_gt(length(s$coverage$shared), 0L)
  expect_true(any(vapply(s$candidate_audit, function(a) length(a$projection$escape$hla_loh) > 0L, TRUE)))
  expect_true(any(vapply(s$candidate_audit, function(a) length(a$projection$escape$apm_loss) > 0L, TRUE)))
  expect_true(any(vapply(s$redundancy_relations, function(r) r$type == "source_region_overlap", TRUE)))
  expect_true(any(b$final_grade == "NOT_GRADED"))
  expect_true(all(vapply(s$candidate_audit, function(a) length(a$reasons) > 0L, TRUE)))
  z <- neo_explain(e$after, candidate_id = b$candidate_id)
  expect_length(z$candidates, nrow(b)); expect_identical(z$vaccine_set, e$selection)
  expect_true(all(vapply(z$candidates, function(i) !is.null(i$vaccine_selection), TRUE)))
  expect_match(neo_explain(e$after, candidate_id = b$candidate_id[[1L]], format = "text"), "deterministic lexicographic heuristic", fixed = TRUE)
  expect_identical(e$selection, Neoantimon2:::vaccine_read(neo_run(e$after)))
})

test_that("single sample no-policy output and enabled candidate invariance are exact", {
  before <- neo_run(recognition_gold_project()); r <- neo_results(before); events <- neo_events(before)
  expect_null(attr(r, "vaccine_set")); expect_identical(neo_results(neo_run(before)), r)
  after <- vaccine_annotate(before, vaccine_binding(min_class_ii = 1L)); x <- neo_results(after)
  expect_true(!is.null(attr(x, "vaccine_set"))); attr(x, "vaccine_set") <- NULL
  expect_identical(x, r); expect_identical(neo_events(after), events)
  expect_identical(neo_multisample(before), neo_multisample(after))
  expect_identical(neo_results(before), r)
  expect_true(all(vapply(neo_explain(after, candidate_id = r$candidate_id)$candidates, function(i) !is.null(i$vaccine_selection), TRUE)))
})

test_that("genuine fusion splice long-read candidates stay NOT_GRADED when selected", {
  before <- presentation_source_pair()$after; rows <- neo_results(before)
  after <- vaccine_annotate(before, vaccine_binding(fill_to_capacity = TRUE)); a <- neo_results(after)
  attr(a, "vaccine_set") <- NULL; expect_identical(a, rows)
  expect_true(all(a$final_grade == "NOT_GRADED"))
  s <- Neoantimon2:::vaccine_read(after)$selections[[1L]]
  expect_gt(s$selected_count, 0L)
  expect_true(all(vapply(s$candidate_audit, function(a) a$projection$frozen_result$final_grade == "NOT_GRADED", TRUE)))
})

test_that("sidecar corruption and policy file changes fail closed", {
  p <- neo_run(recognition_gold_project()); path <- tempfile(fileext = ".json")
  jsonlite::write_json(vaccine_binding()$input, path, auto_unbox = TRUE)
  p <- vaccine_annotate(p, Neoantimon2:::vaccine_policy_bind(path))
  root <- Neoantimon2:::vaccine_root(p, Neoantimon2:::vaccine_inputs(p))
  writeLines("corrupt", file.path(root, "selection.rds"))
  expect_error(neo_results(p), "integrity")
  writeLines("{}", path); expect_error(neo_run(p), "changed since binding")
  b <- vaccine_binding(); b$policy_id <- "tampered"
  expect_error(Neoantimon2:::vaccine_policy_bind(b), "mismatch")
})

test_that("derived validator rejects invented coverage grades and selection identity", {
  b <- vaccine_binding(); p <- vaccine_fixture_projection(b); x <- Neoantimon2:::vaccine_construct(p, b)
  expect_silent(Neoantimon2:::vaccine_validate(x, list(binding = b, projection = p)))
  for (field in c("coverage", "grade", "set_id")) {
    bad <- x
    if (field == "coverage") bad$selections[[1L]]$coverage$persistent <- list("invented")
    if (field == "grade") bad$selections[[1L]]$candidate_audit[[1L]]$projection$frozen_result$final_grade <- "A"
    if (field == "set_id") bad$selections[[1L]]$vaccine_set_id <- "invented"
    expect_error(Neoantimon2:::vaccine_validate(bad, list(binding = b, projection = p)), "derived contract")
  }
})

test_that("public exports match the exact approved milestone API", {
  expect_length(getNamespaceExports("Neoantimon2"), length(current_public_api))
  expect_setequal(getNamespaceExports("Neoantimon2"), current_public_api)
})

test_that("authoritative ambiguous cells cannot become supported coverage", {
  f <- vaccine_fixture(); f$candidates <- f$candidates[1L]; records <- vaccine_records(f)
  for (i in seq_along(records)) {
    records[[i]]$trajectory$cells[[1L]]$state <- "ambiguous"
    records[[i]]$trajectory$summary$trajectory_state <- "insufficient_to_classify"
  }
  p <- Neoantimon2:::vaccine_project_records(records, f$contexts, list(), vaccine_binding()$normalized)
  u <- Neoantimon2:::vaccine_units(p, "patient-1")[[1L]]
  expect_false("S1" %in% unlist(u$atoms$contexts))
  expect_length(u$atoms$shared, 0L); expect_length(u$atoms$persistent, 0L)
  expect_identical(p$candidates[[1L]]$trajectory$cells[[1L]]$state, "ambiguous")
})

test_that("provider-derived epitope keys are consumed without rederiving cores", {
  f <- vaccine_fixture(); f$candidates <- f$candidates[c(1L, 2L)]; records <- vaccine_records(f)
  for (i in seq_along(records)) {
    records[[i]]$result$epitope_key <- NULL
    records[[i]]$evidence$predictor_core_evidence <- list(list(candidate_id = records[[i]]$result$candidate_id,
      presenter = records[[i]]$result$mhc_key, status = "evaluated", epitope_key_status = "derived", epitope_key = "frozen:core"))
  }
  p <- Neoantimon2:::vaccine_project_records(records, f$contexts, list(), vaccine_binding()$normalized)
  expect_true(all(vapply(p$candidates, function(p) identical(p$epitope_keys, list("frozen:core")), TRUE)))
  rel <- Neoantimon2:::vaccine_redundancy(p, "patient-1")
  strong <- Filter(function(r) r$type == "same_epitope_key" && r$strong, rel)
  expect_gt(length(strong), 0L)
  groups <- Neoantimon2:::vaccine_units(p, "patient-1")[[1L]]$strong_groups
  expect_true(all(vapply(strong, function(r) r$group_id %in% unlist(groups), TRUE)))
  for (i in seq_along(records)) records[[i]]$evidence$predictor_core_evidence[[1L]]$epitope_key_status <- "deferred_orientation"
  p <- Neoantimon2:::vaccine_project_records(records, f$contexts, list(), vaccine_binding()$normalized)
  expect_true(all(vapply(p$candidates, function(p) length(p$epitope_keys) == 0L, TRUE)))
})

test_that("treatment metadata never becomes causal or changes set selection", {
  f <- vaccine_fixture(); a <- vaccine_fixture_select(fixture = f)
  f$treatments <- list(list(treatment_id = "therapy-1", before_samples = list("S1"), after_samples = list("S3")))
  b <- vaccine_fixture_select(fixture = f)
  expect_identical(a$selections[[1L]]$selection_trace, b$selections[[1L]]$selection_trace)
  expect_identical(a$selections[[1L]]$vaccine_set_id, b$selections[[1L]]$vaccine_set_id)
})

test_that("policy-relative files exact input and config roundtrips are retained", {
  root <- tempfile("vaccine-config-"); dir.create(root)
  jsonlite::write_json(vaccine_binding()$input, file.path(root, "policy.json"), auto_unbox = TRUE)
  yaml::write_yaml(list(vaccine_set = "policy.json"), file.path(root, "config.yaml"))
  cfg <- neo_config(file.path(root, "config.yaml"), env = character())
  expect_identical(cfg$vaccine_set$physical$path, normalizePath(file.path(root, "policy.json")))
  expect_identical(cfg$vaccine_set$input, jsonlite::read_json(file.path(root, "policy.json")))
  roundtrip <- do.call(neo_config, c(unclass(cfg), list(env = character())))
  expect_identical(roundtrip, cfg)
})

test_that("all unsupported authoritative states block every positive atom quota and fill", {
  f <- vaccine_fixture(); f$candidates <- f$candidates[c(1L, 4L)]
  original <- vaccine_records(f)
  for (state in vaccine_fixture("review-support.json")$unsupported_states) {
    records <- lapply(original, function(r) {
      for (i in seq_along(r$trajectory$cells)) r$trajectory$cells[[i]]$state <- state
      r
    })
    for (fill in c(FALSE, TRUE)) {
      b <- vaccine_binding(required_presenters = list("HLA-A*02:01", "HLA-DRA*01:01~HLA-DRB1*04:01"),
        min_class_i = 1L, min_class_ii = 1L, fill_to_capacity = fill)
      p <- Neoantimon2:::vaccine_project_records(records, f$contexts, list(), b$normalized)
      expect_length(Neoantimon2:::vaccine_units(p, "patient-1"), 0L)
      atoms <- Neoantimon2:::vaccine_atoms(p$candidates, Neoantimon2:::vaccine_context_units(f$contexts, list()), f$contexts)
      expect_identical(atoms, Neoantimon2:::vaccine_atoms_empty())
      s <- Neoantimon2:::vaccine_construct(p, b)$selections[[1L]]
      expect_identical(s$selected_count, 0L); expect_identical(s$class_i_count, 0L); expect_identical(s$class_ii_count, 0L)
      expect_identical(s$coverage, Neoantimon2:::vaccine_atoms_empty())
      expect_identical(s$status, "infeasible"); expect_false(s$feasibility$initial_feasible)
      expect_setequal(vapply(s$unmet_constraints, `[[`, "", "constraint"), c("required_presenter", "min_class_i", "min_class_ii"))
      expect_length(s$unmet_constraints, 4L)
      for (a in s$candidate_audit) {
        expect_identical(a$status, "unselected")
        expect_identical(a$reasons[[1L]]$code, "authoritative_observation_not_supported")
        expect_identical(a$reasons[[1L]]$detail$cells[[1L]]$state, state)
        r <- Filter(function(r) r$result$candidate_id == a$candidate_id, records)[[1L]]
        expect_identical(a$projection$frozen_result, r$result)
        expect_identical(a$projection$trajectory, r$trajectory)
        expect_identical(a$projection$evidence, r$evidence)
        expect_identical(a$eligibility$status, "eligible")
      }
    }
  }
})

test_that("one observed sample supports a key while ambiguous samples retain unselected evidence", {
  f <- vaccine_fixture(); f$candidates <- f$candidates[1L]
  records <- vaccine_records(f)
  for (i in seq_along(records)) {
    records[[i]]$trajectory$cells[[1L]]$state <- "ambiguous"
    records[[i]]$trajectory$cells[[3L]]$state <- "not_applicable"
  }
  b <- vaccine_binding(required_presenters = list("HLA-A*02:01"), min_class_i = 1L, fill_to_capacity = TRUE)
  p <- Neoantimon2:::vaccine_project_records(records, f$contexts, list(), b$normalized)
  u <- Neoantimon2:::vaccine_units(p, "patient-1")[[1L]]
  expect_identical(vapply(u$observations, `[[`, "", "sample_id"), "S2")
  expect_length(u$atoms$presenters, 1L); expect_identical(u$atoms$contexts, list("S2"))
  expect_identical(u$atoms$regions, list("primary-right")); expect_length(u$atoms$clonal, 1L)
  expect_length(u$atoms$shared, 0L); expect_length(u$atoms$persistent, 0L)
  s <- Neoantimon2:::vaccine_construct(p, b)$selections[[1L]]
  expect_identical(s$status, "feasible"); expect_identical(s$selected_count, 1L); expect_identical(s$class_i_count, 1L)
  expect_identical(s$selected_candidates[[1L]]$candidate_ids, list(u$observations[[1L]]$candidate_id))
  for (sample in c("S1", "S2", "S3")) {
    a <- Filter(function(a) a$sample_id == sample, s$candidate_audit)[[1L]]
    expect_identical(a$status, if (sample == "S2") "selected" else "unselected")
    expect_identical(a$reasons[[1L]]$code, if (sample == "S2") "selected" else "authoritative_observation_not_supported")
    r <- Filter(function(r) r$result$sample_id == sample, records)[[1L]]
    expect_identical(a$projection$trajectory, r$trajectory)
    expect_identical(a$projection$frozen_result, r$result)
    expect_identical(a$projection$evidence, r$evidence)
  }
})

test_that("mixed normal tumor observations sharing a key have sample local scope reasons", {
  f <- vaccine_fixture(); f$candidates <- f$candidates[1L]; f$contexts[[1L]]$role <- "normal"
  b <- vaccine_binding(fill_to_capacity = TRUE, min_class_i = 1L)
  s <- vaccine_fixture_select(b, f)$selections[[1L]]
  expect_identical(s$status, "feasible"); expect_identical(s$selected_count, 1L)
  expect_setequal(unlist(s$coverage$contexts), c("S2", "S3"))
  expect_false("primary-left" %in% unlist(s$coverage$regions)); expect_length(s$coverage$clonal, 2L)
  normal <- Filter(function(a) a$sample_id == "S1", s$candidate_audit)[[1L]]
  expect_identical(normal$status, "unselected"); expect_identical(normal$reasons[[1L]]$code, "outside_tumor_scope")
  expect_identical(normal$reasons[[1L]]$constraint, "subject_context_scope")
  expect_true(all(vapply(Filter(function(a) a$sample_id != "S1", s$candidate_audit), function(a) a$status == "selected", TRUE)))
  expect_false(normal$candidate_id %in% unlist(s$selected_candidates[[1L]]$candidate_ids))
  expect_true(normal$candidate_key %in% vaccine_keys(list(selections = list(s))))
})

test_that("presenter admission limits count selected candidates rather than unique coverage", {
  f <- vaccine_fixture(); f$candidates <- f$candidates[1:2]
  unlimited <- vaccine_fixture_select(vaccine_binding(fill_to_capacity = TRUE), f)$selections[[1L]]
  expect_identical(unlimited$selected_count, 2L)
  expect_identical(vapply(unlimited$selected_candidates, `[[`, "", "presenter"), rep("HLA-A*02:01", 2L))
  b <- vaccine_binding(fill_to_capacity = TRUE, max_per_presenter = 1L)
  s <- vaccine_fixture_select(b, f)$selections[[1L]]
  counts <- table(vapply(s$selected_candidates, `[[`, "", "presenter"))
  expect_identical(s$selected_count, 1L); expect_true(all(counts <= b$normalized$max_per_presenter))
  expect_true(any(vapply(s$candidate_audit, function(a) any(vapply(a$reasons, function(r) r$constraint == "max_per_presenter", TRUE)), TRUE)))
})

test_that("raw duplicate sets remain exact while normalized policy and semantic identity deduplicate", {
  raw <- vaccine_review_duplicate_policy()
  b <- Neoantimon2:::vaccine_policy_bind(raw)
  unique <- raw; for (k in names(raw)[vapply(raw, is.list, TRUE)]) unique[[k]] <- as.list(sort(unique(unlist(raw[[k]])), method = "radix"))
  a <- Neoantimon2:::vaccine_policy_bind(unique)
  expect_identical(b$input, raw); expect_identical(a$normalized, b$normalized); expect_identical(a$policy_id, b$policy_id)
  x <- vaccine_fixture_select(b); y <- vaccine_fixture_select(a)
  expect_gt(x$selections[[1L]]$selected_count, 0L)
  expect_identical(x$provenance$policy_input, raw); expect_identical(x$normalized_policy, b$normalized)
  expect_identical(vaccine_keys(x), vaccine_keys(y))
  expect_identical(x$selections[[1L]]$vaccine_set_id, y$selections[[1L]]$vaccine_set_id)
  expect_identical(x$selections[[1L]]$selection_trace, y$selections[[1L]]$selection_trace)
  for (k in names(raw)[vapply(raw, is.list, TRUE)]) {
    expect_gt(length(raw[[k]]), length(b$normalized[[k]]))
    expect_identical(b$normalized[[k]], as.list(sort(unique(unlist(raw[[k]])), method = "radix")))
  }
})

test_that("explicit policy nulls fail for every field and omitted optional fields default", {
  raw_schema <- jsonlite::read_json(system.file("schema", "vaccine-set-policy-1.0.0.json", package = "Neoantimon2"))
  base <- vaccine_binding()$input
  for (k in names(raw_schema$properties)) {
    raw <- base; raw[k] <- list(NULL)
    expect_error(Neoantimon2:::vaccine_policy_bind(raw), paste0("Explicit null.*", k))
  }
  for (k in c("max_class_i", "fill_to_capacity", "allowlist")) {
    raw <- base; raw[k] <- list(NULL); path <- tempfile(fileext = ".json")
    jsonlite::write_json(raw, path, auto_unbox = TRUE, null = "null")
    expect_error(Neoantimon2:::vaccine_policy_bind(path), paste0("Explicit null.*", k))
  }
  n <- vaccine_binding()$normalized
  expect_identical(n$min_class_i, 0L); expect_identical(n$min_class_ii, 0L)
  for (k in c("fill_to_capacity", "reduce_overlap", "forbid_confirmed_exact_self")) expect_identical(n[[k]], FALSE)
  for (k in c("required_presenters", "allowlist", "blocklist", "forbidden_residues", "forbidden_motifs", "require_evidence")) expect_identical(n[[k]], list())
  for (k in c("max_class_i", "max_class_ii", "max_per_presenter", "max_per_redundancy_group", "min_length", "max_length")) expect_false(k %in% names(n))
  expect_false(any(vapply(n, is.null, TRUE)))
})
