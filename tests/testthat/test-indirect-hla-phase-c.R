test_that('complete public neo_results adds a compact profile-aware summary', {
  p <- phase_c_project(); r <- neo_results(p); s <- attr(r, 'indirect_hla_result')
  expect_identical(s, phase_c_case()$public_result)
  c <- s$contexts[[1L]]
  expect_identical(c$complete_count, 5L); expect_identical(c$conditional_count, 5L)
  expect_identical(c$completeness_status, 'complete'); expect_identical(c$analysis_status, 'evaluated')
  expect_identical(c$included_unit_count, 5L); expect_identical(c$self_excluded_count, 1L)
  expect_identical(c$unevaluable_count, 0L); expect_identical(c$evaluated_presenter_count, 1L)
  expect_identical(c$profile_version, '1.0.0'); expect_false(c$official_reproduction)
  expect_match(c$claim_boundary, 'not official PIRCHE reproduction')
  expect_false(any(c('graph', 'member_ids', 'provider', 'coverage') %in% names(c)))
})
test_that('incomplete self summary preserves null and conditional semantics', {
  s <- attr(neo_results(phase_c_project('incomplete_self')), 'indirect_hla_result')$contexts[[1L]]
  expect_null(s$complete_count); expect_identical(s$completeness_status, 'incomplete')
  expect_identical(s$conditional_count, 0L); expect_identical(s$self_excluded_count, 1L)
  expect_true('self_search_incomplete' %in% unlist(s$reason_codes))
  expect_identical(s$conditional_semantics, 'inspectable_conditional_result_only')
})
test_that('incomplete presenter scope exposes omitted DRB1 without a false zero', {
  s <- attr(neo_results(phase_c_project('incomplete_scope')), 'indirect_hla_result')$contexts[[1L]]
  expect_null(s$complete_count); expect_identical(s$conditional_count, 5L)
  expect_false(s$presenter_scope$complete)
  expect_identical(s$presenter_scope$undeclared_beta_alleles, list('HLA-DRB1*03:01'))
  expect_true('presenter_scope_incomplete' %in% unlist(s$reason_codes))
})
test_that('unsupported accessory evidence remains visible and contributes no count', {
  s <- attr(neo_results(phase_c_project('accessory')), 'indirect_hla_result')$contexts[[1L]]
  e <- phase_c_detail('accessory')
  expect_null(s$complete_count); expect_identical(s$conditional_count, 5L)
  unsupported <- Filter(function(d) grepl('DRB3', d$presenter), e$decisions)
  expect_length(unsupported, 6L)
  expect_true(all(vapply(unsupported, function(d) d$outcome_kind == 'presenter_out_of_profile' && is.null(d$counted_unit_key), TRUE)))
  expect_false(any(vapply(e$included_units, function(u) grepl('DRB3', u$presenter), TRUE)))
})
test_that('provider failure summary is unavailable evidence rather than nonbinding', {
  s <- attr(neo_results(phase_c_project('provider_failure')), 'indirect_hla_result')$contexts[[1L]]
  expect_null(s$complete_count); expect_identical(s$unevaluable_count, 1L)
  expect_identical(s$binding_policy_excluded_count, 0L)
  expect_true('provider_execution_failed' %in% unlist(s$reason_codes))
  d <- Filter(function(d) d$outcome_kind == 'provider_failure', phase_c_detail('provider_failure')$decisions)
  expect_length(d, 1L); expect_null(d[[1L]]$affinity_nM); expect_null(d[[1L]]$affinity_under_frozen_policy)
})
test_that('provider integrity and pinning failures retain frozen fail-closed errors', {
  x <- indirect_input(); x$b$provider$sha256 <- strrep('0', 64)
  expect_error(indirect_test(x), 'integrity mismatch')
  m <- indirect_precomputed_contract(); m$provenance$model_artifact <- NULL
  expect_error(indirect_test(indirect_with_binding(m)), 'missing required field')
  p <- phase_c_project(); p$config$indirect_hla$provider$sha256 <- strrep('0',64)
  expect_error(neo_results(p), 'integrity mismatch'); expect_error(neo_explain(p), 'integrity mismatch')
})
test_that('included-unit chain retains affinity threshold provenance and source links', {
  e <- neo_explain(phase_c_project())$indirect_hla_explanation
  expect_identical(e, phase_c_case()$public_explanation)
  c <- e$contexts[[1L]]
  expect_length(c$included_units, 5L)
  for (u in c$included_units) {
    expect_length(u$source_occurrence_ids, 5L); expect_length(u$occurrence_ids, 5L)
    expect_identical(u$threshold_policy$comparator, '<'); expect_identical(u$threshold_policy$value, 1000)
    expect_identical(u$threshold_policy$unit, 'nM'); expect_length(u$supporting_evidence_ids, 1L)
    d <- Filter(function(d) d$epitope_key == u$counted_unit_key, c$decisions)
    expect_true(d[[1L]]$affinity_under_frozen_policy); expect_identical(d[[1L]]$affinity_nM, 500)
    expect_identical(d[[1L]]$self_comparison_status, 'supported_nonself_within_scope')
    expect_match(u$inclusion_reason, 'supported_nonself')
  }
  expect_identical(e$provenance$provider$kind, 'synthetic_fixture')
  expect_identical(e$provenance$provider$execution, 'not_executed')
})
test_that('self-excluded chain preserves positive self evidence and universe', {
  e <- phase_c_detail(); d <- Filter(function(d) d$state == 'excluded_self', e$decisions)
  expect_length(d, 1L); expect_identical(d[[1L]]$outcome_kind, 'supported_self_match')
  expect_identical(d[[1L]]$self_comparison_status, 'self_match'); expect_null(d[[1L]]$counted_unit_key)
  q <- Neoantimon2:::indirect_report_one(e$graph$nodes$self_comparison, 'id', d[[1L]]$self_comparison_id)$data
  expect_gt(length(q$matches), 0L); expect_identical(q$semantics$self_universe_id, e$summary$self_universe_id)
  expect_false(q$biological_self_completeness_established)
})
test_that('insufficient self search exposes a distinct unevaluable chain', {
  e <- phase_c_detail('incomplete_self'); d <- Filter(function(d) d$state == 'unevaluable', e$decisions)
  expect_length(d, 5L)
  expect_true(all(vapply(d, function(d) d$outcome_kind == 'insufficient_self_search' && d$self_comparison_status == 'insufficient_search', TRUE)))
  expect_true(all(vapply(d, function(d) 'self_search_incomplete' %in% unlist(d$reason_codes), TRUE)))
})
test_that('missing provider, reference failure, core ambiguity and absent universe are distinct', {
  expect_true(all(vapply(phase_c_detail('missing_provider')$decisions, function(d) d$outcome_kind == 'insufficient_self_search', TRUE)))
  failed <- phase_c_detail('reference_failure')
  expect_identical(failed$reference_status, 'failed')
  expect_true('reference_provider_failure' %in% unlist(failed$summary$reason_codes))
  expect_true(all(vapply(failed$generation_coverage, function(l) l$status == 'failed', TRUE)))
  expect_null(failed$summary$complete_count); expect_length(failed$decisions, 0L)
  expect_true(any(vapply(phase_c_detail('core_ambiguous')$decisions, function(d) d$outcome_kind == 'core_assignment_unevaluable', TRUE)))
  e <- phase_c_detail('missing_universe')
  expect_null(e$summary$self_universe_id); expect_gt(length(e$unidentified_self_comparisons), 0L)
  expect_true(any(vapply(e$decisions, function(d) d$outcome_kind == 'missing_self_universe', TRUE)))
})
test_that('multiple source origins survive one counted molecular unit', {
  e <- phase_c_detail()
  expect_length(e$graph$nodes$source_occurrence, 5L); expect_length(e$graph$nodes$peptide_occurrence, 30L)
  expect_length(e$graph$nodes$peptide, 6L); expect_length(e$graph$nodes$pmhc, 6L)
  expect_length(e$graph$nodes$aggregate_unit, 5L)
  for (u in e$included_units) expect_length(u$source_occurrence_ids, 5L)
})
test_that('same core presenter union explains duplicate full ligands without counting twice', {
  e <- phase_c_detail('union'); expect_identical(e$summary$complete_count, 5L)
  expect_identical(e$reconciliation$included_pmhc_records, 6L)
  expect_identical(e$reconciliation$duplicate_included_union_members, 1L)
  u <- Filter(function(u) length(u$supporting_pmhc_ids) == 2L, e$included_units)
  expect_length(u, 1L); expect_length(u[[1L]]$full_ligands, 2L); expect_length(u[[1L]]$occurrence_ids, 10L)
  d <- Filter(function(d) d$aggregation_role == 'duplicate_union_member', e$decisions)
  expect_length(d, 1L); expect_true('duplicate_core_presenter_union_member' %in% unlist(d[[1L]]$explanation_codes))
  expect_identical(d[[1L]]$state, 'included')
})
test_that('different presenters retain separate core presenter union members', {
  e <- phase_c_detail('two_presenters')
  expect_identical(e$summary$complete_count, 10L); expect_length(e$included_units, 10L)
  expect_length(unique(vapply(e$included_units, `[[`, '', 'presenter')), 2L)
  expect_length(unique(vapply(e$included_units, `[[`, '', 'core_sequence')), 5L)
})
test_that('independently computed contexts stay separate in schema-valid collection order', {
  a <- phase_c_case(); input <- indirect_input(); input$a$contexts[[1L]]$occasion_id <- 'occasion:independent-2'
  subjects <- c('donor-1' = 'donor-2', 'patient-1' = 'patient-2')
  for (i in seq_along(input$a$subjects)) input$a$subjects[[i]]$subject_id <- unname(subjects[input$a$subjects[[i]]$subject_id])
  for (i in seq_along(input$a$typing)) input$a$typing[[i]]$subject_id <- unname(subjects[input$a$typing[[i]]$subject_id])
  for (field in names(input$a$contexts[[1L]])) if (endsWith(field, 'subject_id')) input$a$contexts[[1L]][[field]] <- unname(subjects[input$a$contexts[[1L]][[field]]])
  input$a$self_search[[1L]]$self_subject_id <- 'patient-2'; input$b$presenters[[1L]]$subject_id <- 'patient-2'
  b <- indirect_test(input)
  first <- a$public_result$contexts[[1L]]
  second <- Neoantimon2:::indirect_result_context(b$result, b$foundation, b$result$contexts[[1L]])
  x <- Neoantimon2:::indirect_result_collection(list(first, second))
  y <- Neoantimon2:::indirect_result_collection(list(second, first))
  expect_identical(x, y); expect_length(x$contexts, 2L)
  expect_identical(vapply(x$contexts, `[[`, 1L, 'complete_count'), c(5L, 5L))
  expect_false(any(c('total_count','patient_risk','aggregate_score') %in% names(x)))
  expect_true(phase_c_schema('result')(x))
  expect_setequal(vapply(x$contexts, `[[`, '', 'source_subject_id'), c('donor-1','donor-2'))
  expect_setequal(vapply(x$contexts, `[[`, '', 'self_subject_id'), c('patient-1','patient-2'))
  ea <- a$public_explanation$contexts[[1L]]
  eb <- Neoantimon2:::indirect_public_explanation(b$result, b$foundation)$contexts[[1L]]
  forward <- Neoantimon2:::indirect_explanation_contexts(list(ea,eb))
  backward <- Neoantimon2:::indirect_explanation_contexts(list(eb,ea))
  expect_identical(forward, backward); expect_length(forward, 2L)
  multi <- a$public_explanation; multi$contexts <- forward
  expect_true(phase_c_schema('explanation')(multi))
  expect_error(Neoantimon2:::indirect_result_collection(list(first, first)), 'Duplicate')
  # The frozen runtime still rejects duplicate scenario contexts; reporting does not weaken it.
  input <- indirect_input(); input$a$contexts <- rep(input$a$contexts, 2L)
  expect_error(indirect_test(input), 'unique scenario_id')
})
test_that('no-allo results and explain delegate with exact historical structure', {
  p <- phase_c_runtime_project()
  expect_identical(Neoantimon2:::allo_results_append(neo_results(p), p), neo_results(p))
  ids <- 'unmatched:synthetic-snapshot'
  expect_identical(neo_explain(p, candidate_id = ids), Neoantimon2:::neo_explain_core(p, candidate_id = ids))
  expect_error(neo_explain(p), 'At least one explanation query')
  expect_null(attr(neo_results(p), 'indirect_hla_result'))
})
test_that('Phase A only retains exact foundation summary and explanation', {
  p <- allo_test_enable(phase_c_runtime_project())
  expect_identical(attr(neo_results(p), 'allo_foundation'), Neoantimon2:::allo_summary(p$allo_foundation))
  expect_identical(neo_explain(p)$allo_foundation, p$allo_foundation)
  expect_null(attr(neo_results(p), 'indirect_hla_result')); expect_null(neo_explain(p)$indirect_hla_explanation)
  expect_match(neo_explain(p, format = 'text'), 'ALLO PHASE A')
})
test_that('Phase C leaves the Phase B internal object and old compact metadata intact', {
  p <- phase_c_project(); before <- serialize(p$indirect_hla, NULL)
  r <- neo_results(p); e <- neo_explain(p)
  expect_identical(serialize(p$indirect_hla, NULL), before)
  expect_identical(e$indirect_hla_baseline, p$indirect_hla)
  expect_identical(attr(r, 'indirect_hla_baseline'), Neoantimon2:::indirect_summary(p$indirect_hla))
  expect_identical(Neoantimon2:::indirect_read(p), p$indirect_hla)
})
test_that('somatic ordering fields identities and grades do not change under reporting', {
  p <- phase_c_runtime_project(); before <- neo_results(p); q <- indirect_enable(p)
  after <- neo_results(q); attr(after, 'allo_foundation') <- NULL
  attr(after, 'indirect_hla_baseline') <- NULL; attr(after, 'indirect_hla_result') <- NULL
  expect_identical(after, before)
  expect_identical(p$analysis_id, q$analysis_id); expect_identical(p$config_fingerprint, q$config_fingerprint)
  expect_identical(p$input_fingerprint, q$input_fingerprint)
})
test_that('S4 v4.2 records do not infer an equality comparator or exact domains', {
  x <- Neoantimon2:::indirect_compatibility()
  t <- Neoantimon2:::indirect_report_one(x$targets, 'target_label', 'neo-indirect-hla-compatibility/pirche-4.2-described-affinity')
  expect_identical(t$application_version, '4.2'); expect_false(t$official_reproduction)
  d <- function(n) Neoantimon2:::indirect_report_one(t$dimensions, 'dimension', n)
  expect_identical(d('threshold_operator')$status, 'unresolved')
  expect_identical(d('reference_domains')$status, 'differs')
  expect_identical(d('full_peptide_vs_core')$status, 'matched')
  expect_identical(d('aggregation')$status, 'unresolved')
})
test_that('S9 v4.2 Frost records expose provider metric differences and missing artifacts', {
  x <- Neoantimon2:::indirect_compatibility()
  t <- Neoantimon2:::indirect_report_one(x$targets, 'target_label', 'neo-indirect-hla-compatibility/pirche-4.2-s9-frost-described-concepts')
  d <- function(n) Neoantimon2:::indirect_report_one(t$dimensions, 'dimension', n)
  expect_identical(d('provider_family')$status, 'differs'); expect_identical(d('threshold_units')$status, 'differs')
  expect_identical(d('provider_version_pinning')$status, 'unresolved')
  expect_match(d('binding_threshold')$documented_behavior, '300 permille')
})
test_that('S10 build-unpinned evidence never inherits a v4.2 release label', {
  t <- Neoantimon2:::indirect_report_one(Neoantimon2:::indirect_compatibility()$targets, 'target_label', 'neo-indirect-hla-compatibility/pirche-t2-s10-frost-build-unpinned-concepts')
  expect_null(t$application_version); expect_identical(t$matrix_profile, 'G7')
  expect_match(t$configuration, 'build unpinned')
  d <- Neoantimon2:::indirect_report_one(t$dimensions, 'dimension', 'provider_version_pinning')
  expect_identical(d$status, 'unresolved'); expect_false(t$official_reproduction)
})
test_that('later v4x v4.5 and current TxPredictor remain outside frozen audit scope', {
  t <- Filter(function(t) t$audit_scope == 'outside_audit_scope', Neoantimon2:::indirect_compatibility()$targets)
  expect_length(t, 1L); expect_match(t[[1L]]$configuration, 'v4.5')
  expect_true(all(vapply(t[[1L]]$dimensions, function(d) d$status == 'outside_audit_scope', TRUE)))
  expect_null(t[[1L]]$application_version); expect_false(t[[1L]]$official_reproduction)
})
test_that('legacy S7 inclusive boundary differs from the native strict policy', {
  t <- Filter(function(t) t$matrix_profile == 'G2', Neoantimon2:::indirect_compatibility()$targets)[[1L]]
  d <- Neoantimon2:::indirect_report_one(t$dimensions, 'dimension', 'threshold_operator')
  expect_identical(d$status, 'differs'); expect_match(d$documented_behavior, '<=')
})
test_that('every compatibility target is claim-safe and U or assumption cells remain unresolved', {
  x <- Neoantimon2:::indirect_compatibility()
  expect_false(x$official_reproduction); expect_length(x$targets, 14L)
  for (t in x$targets) {
    expect_false(t$official_reproduction); expect_length(t$dimensions, 19L)
    expect_false(any(c('compatible','compatibility_score','percentage_compatibility') %in% names(t)))
    for (d in t$dimensions) if (d$classification %in% c('A','U','A/U','R/A')) expect_identical(d$status, 'unresolved')
  }
  expect_identical(x$audit_source_sha256, phase_c_json_fixture('compatibility-expectations.json')$matrix_sha256)
})
test_that('compatibility fixture independently pins audited configuration distinctions', {
  expected <- phase_c_json_fixture('compatibility-expectations.json')
  x <- Neoantimon2:::indirect_compatibility()
  for (c in expected$checks) {
    t <- Neoantimon2:::indirect_report_one(x$targets, 'target_label', c$target_label)
    d <- Neoantimon2:::indirect_report_one(t$dimensions, 'dimension', c$dimension)
    expect_identical(d$status, c$status); expect_false(t$official_reproduction)
  }
})
test_that('reference_release is exactly one independent dimension on all 14 targets', {
  x <- neo_explain(phase_c_project())$indirect_hla_explanation$compatibility
  expected <- phase_c_json_fixture('compatibility-expectations.json')$reference_release_checks
  expect_length(x$targets, 14L); expect_length(expected, 14L)
  expect_identical(sum(vapply(x$targets, function(t) length(t$dimensions), 0L)), 266L)
  for (c in expected) {
    t <- Neoantimon2:::indirect_report_one(x$targets, 'target_label', c$target_label)
    names <- vapply(t$dimensions, `[[`, '', 'dimension')
    expect_length(names, 19L); expect_identical(anyDuplicated(names), 0L)
    expect_identical(sum(names == 'reference_release'), 1L)
    d <- Neoantimon2:::indirect_report_one(t$dimensions, 'dimension', 'reference_release')
    expect_identical(d$status, c$status); expect_identical(d$classification, c$classification)
    expect_match(d$neoantimon2_behavior, 'Explicitly pinned analysis reference')
    expect_match(d$locator, 'pirche_reconstruction_matrix.md', fixed = TRUE)
    if (!is.null(c$release)) expect_match(d$documented_behavior, c$release, fixed = TRUE)
    expect_false(t$official_reproduction)
  }
  expect_false(x$official_reproduction)
})
test_that('reference_release preserves distinct S9 versions and the S4 same-reference condition', {
  x <- Neoantimon2:::indirect_compatibility()
  ref <- function(label) {
    t <- Neoantimon2:::indirect_report_one(x$targets, 'target_label', paste0('neo-indirect-hla-compatibility/', label))
    Neoantimon2:::indirect_report_one(t$dimensions, 'dimension', 'reference_release')
  }
  v3 <- ref('pirche-v3-s9-described-affinity')
  v42 <- ref('pirche-4.2-s9-frost-described-concepts')
  expect_match(v3$documented_behavior, 'IMGT 3.47', fixed = TRUE)
  expect_match(v42$documented_behavior, 'IMGT 3.54', fixed = TRUE)
  expect_false(grepl('3.54', v3$documented_behavior, fixed = TRUE))
  expect_false(grepl('3.47', v42$documented_behavior, fixed = TRUE))
  expect_false(identical(v3$documented_behavior, v42$documented_behavior))
  expect_identical(v3$status, 'unresolved'); expect_identical(v42$status, 'unresolved')
  s4 <- Filter(function(t) t$matrix_profile == 'G5', x$targets)
  expect_length(s4, 4L)
  for (t in s4) {
    d <- Neoantimon2:::indirect_report_one(t$dimensions, 'dimension', 'reference_release')
    expect_match(d$documented_behavior, 'Same IPD-IMGT/HLA 3.47', fixed = TRUE)
    expect_match(d$documented_behavior, 'controlled application 3.3 / 4.2 SOT and SCT', fixed = TRUE)
    expect_identical(d$status, 'unresolved')
  }
})
test_that('reference_release unknowns and later-version audit boundaries cannot become matches', {
  x <- Neoantimon2:::indirect_compatibility()
  unknown <- Filter(function(t) t$matrix_profile %in% c('G0','G1','G2','G3','G4','G7'), x$targets)
  expect_length(unknown, 7L)
  for (t in unknown) {
    d <- Neoantimon2:::indirect_report_one(t$dimensions, 'dimension', 'reference_release')
    expect_identical(d$status, 'unresolved'); expect_identical(d$classification, 'U')
    expect_false(grepl('3.47|3.54', d$documented_behavior))
    expect_match(d$documented_behavior, 'not pinned|not provided')
    bad <- x; i <- which(vapply(bad$targets, function(q) q$target_label == t$target_label, TRUE))
    j <- which(vapply(t$dimensions, function(d) d$dimension == 'reference_release', TRUE))
    bad$targets[[i]]$dimensions[[j]]$status <- 'matched'
    expect_false(phase_c_schema('compatibility')(bad))
  }
  s10 <- Filter(function(t) t$matrix_profile == 'G7', unknown)
  expect_length(s10, 2L)
  for (t in s10) {
    d <- Neoantimon2:::indirect_report_one(t$dimensions, 'dimension', 'reference_release')
    expect_match(d$documented_behavior, 'completed sequence/model snapshot not provided', fixed = TRUE)
  }
  outside <- Filter(function(t) t$audit_scope == 'outside_audit_scope', x$targets)
  expect_length(outside, 1L)
  d <- Neoantimon2:::indirect_report_one(outside[[1L]]$dimensions, 'dimension', 'reference_release')
  expect_identical(d$status, 'outside_audit_scope')
  expect_false(grepl('3.47|3.54', d$documented_behavior))
  expect_match(d$documented_behavior, 'no inherited v4.2 reference-release judgment', fixed = TRUE)
  expect_false(any(vapply(x$targets, function(t)
    Neoantimon2:::indirect_report_one(t$dimensions, 'dimension', 'reference_release')$status == 'matched', TRUE)))
})
test_that('reference_release schema rejects missing and semantically duplicate dimensions', {
  x <- Neoantimon2:::indirect_compatibility(); valid <- phase_c_schema('compatibility')
  schema <- function(name) jsonlite::read_json(system.file('schema', paste0('indirect-hla-', name, '-1.0.0.json'), package = 'Neoantimon2'))
  expect_identical(schema('explanation')$properties$compatibility, schema('compatibility'))
  expect_true(valid(x))
  # Check every target, including the outside-scope marker, and every known name.
  for (i in seq_along(x$targets)) {
    j <- which(vapply(x$targets[[i]]$dimensions, function(d) d$dimension == 'reference_release', TRUE))
    bad <- x; bad$targets[[i]]$dimensions <- bad$targets[[i]]$dimensions[-j]
    expect_false(valid(bad))
    bad <- x; k <- if (j == 1L) 2L else 1L
    bad$targets[[i]]$dimensions[[k]] <- bad$targets[[i]]$dimensions[[j]]
    bad$targets[[i]]$dimensions[[k]]$documented_behavior <- 'Different note, same reference_release dimension'
    expect_false(valid(bad))
    bad <- x; bad$targets[[i]]$dimensions[[j]] <- bad$targets[[i]]$dimensions[[k]]
    expect_false(valid(bad))
  }
  for (j in seq_along(x$targets[[1L]]$dimensions)) {
    bad <- x; k <- if (j == 1L) 2L else 1L
    bad$targets[[1L]]$dimensions[[j]] <- bad$targets[[1L]]$dimensions[[k]]
    expect_false(valid(bad))
  }
})
test_that('reference_release affects only audit IDs and keeps target dimension order deterministic', {
  x <- Neoantimon2:::indirect_compatibility()
  old <- x
  for (i in seq_along(old$targets)) old$targets[[i]]$dimensions <-
    Filter(function(d) d$dimension != 'reference_release', old$targets[[i]]$dimensions)
  old <- Neoantimon2:::indirect_audit_projection(old)
  expect_true(all(vapply(seq_along(x$targets), function(i)
    x$targets[[i]]$audit_record_id != old$targets[[i]]$audit_record_id, TRUE)))
  changed <- x
  for (i in seq_along(changed$targets)) {
    j <- which(vapply(changed$targets[[i]]$dimensions, function(d) d$dimension == 'reference_release', TRUE))
    changed$targets[[i]]$dimensions[[j]]$documented_behavior <- 'Test-only alternative reference-release annotation'
  }
  changed <- Neoantimon2:::indirect_audit_projection(changed)
  expect_true(all(vapply(seq_along(x$targets), function(i)
    x$targets[[i]]$audit_record_id != changed$targets[[i]]$audit_record_id, TRUE)))
  shuffled <- x; shuffled$targets <- rev(shuffled$targets)
  for (i in seq_along(shuffled$targets)) shuffled$targets[[i]]$dimensions <- rev(shuffled$targets[[i]]$dimensions)
  expect_identical(Neoantimon2:::indirect_audit_projection(shuffled), x)
  for (t in x$targets) expect_identical(vapply(t$dimensions, `[[`, '', 'dimension'),
    sort(vapply(t$dimensions, `[[`, '', 'dimension'), method = 'radix'))
  expect_identical(vapply(x$targets, `[[`, '', 'target_label'), vapply(old$targets, `[[`, '', 'target_label'))
  for (name in c('complete', 'incomplete_self', 'missing_source', 'union')) {
    p <- phase_c_project(name); before <- serialize(p$indirect_hla, NULL)
    current <- neo_explain(p)$indirect_hla_explanation; result <- neo_results(p)
    prior <- local({
      testthat::local_mocked_bindings(indirect_compatibility = function() old, .package = 'Neoantimon2')
      list(explanation = neo_explain(p)$indirect_hla_explanation, result = neo_results(p))
    })
    expect_identical(prior$explanation$compatibility, old)
    expect_identical(current$compatibility, x)
    current$compatibility <- NULL; prior$explanation$compatibility <- NULL
    # Includes every biological graph ID, unit, count, decision, reconciliation
    # record and physical provenance field, across complete and incomplete cases.
    expect_identical(current, prior$explanation)
    expect_identical(result, prior$result)
    expect_identical(serialize(p$indirect_hla, NULL), before)
  }
  expect_identical(phase_c_summary()$complete_count, 5L)
  expect_null(phase_c_summary('incomplete_self')$complete_count)
})
test_that('compatibility overlay cannot change biological identities or baseline counts', {
  p <- phase_c_project(); before <- serialize(p$indirect_hla, NULL)
  e <- neo_explain(p)$indirect_hla_explanation
  e$compatibility$targets[[1L]]$limitations <- list('test-only descriptive annotation')
  expect_identical(serialize(p$indirect_hla, NULL), before)
  expect_identical(e$contexts[[1L]]$summary$allo_analysis_id, p$indirect_hla$allo_analysis_id)
  expect_identical(e$contexts[[1L]]$summary$allo_policy_id, p$indirect_hla$allo_policy_id)
  expect_identical(e$contexts[[1L]]$summary$complete_count, 5L)
  expect_true(all(vapply(e$compatibility$targets, function(t) startsWith(t$audit_record_id, 'audit-record-sha256:'), TRUE)))
})
test_that('graph references and embedded IDs resolve for complete and incomplete evidence', {
  for (name in c('complete','incomplete_self','incomplete_scope','accessory','provider_failure','missing_provider','missing_source','missing_universe','union')) {
    e <- phase_c_detail(name)
    expect_silent(Neoantimon2:::indirect_graph_integrity(e$graph))
    nodes <- e$graph$nodes
    for (table in nodes) expect_false(anyDuplicated(vapply(table, `[[`, '', 'id')) > 0L)
  }
})
test_that('orphan graph and embedded references are rejected', {
  g <- phase_c_detail()$graph; g$edges[[1L]]$to_id <- 'missing'
  expect_error(Neoantimon2:::indirect_graph_integrity(g), 'Orphan')
  g <- phase_c_detail()$graph; g$nodes$pmhc[[1L]]$data$occurrence_ids <- list(paste0('indirect-peptide-occurrence-id-v1:', strrep('0',64)))
  expect_error(Neoantimon2:::indirect_graph_integrity(g), 'Orphan embedded')
})
test_that('all counted units are reachable from the context through evidence links', {
  e <- phase_c_detail(); edges <- e$graph$edges
  reachable <- paste('context', e$allo_context_id)
  repeat {
    old <- reachable
    for (edge in edges) if (paste(edge$from_type, edge$from_id) %in% reachable)
      reachable <- unique(c(reachable, paste(edge$to_type, edge$to_id)))
    if (identical(old, reachable)) break
  }
  expect_true(all(paste('aggregate_unit', vapply(e$included_units, `[[`, '', 'counted_unit_key')) %in% reachable))
  expect_true(all(c('source_origin','generates','molecular_ligand','presented_as','core_assignment','compared_by','supported_by','informs_decision','aggregation_membership') %in% vapply(edges, `[[`, '', 'relation')))
})
test_that('detailed reconciliation follows frozen units without a conservation equation', {
  for (name in c('complete','incomplete_self','incomplete_scope','accessory','provider_failure','nonbinder','union','two_presenters')) {
    e <- phase_c_detail(name); r <- e$reconciliation
    expect_identical(r$status, 'reconciled')
    expect_identical(r$conditional_count, length(e$included_units))
    expect_identical(r$self_excluded_unique_core_presenter_units, e$summary$self_excluded_count)
    expect_identical(r$binding_policy_excluded_pmhc_records, e$summary$binding_policy_excluded_count)
    expect_identical(r$unevaluable_pmhc_records, e$summary$unevaluable_count)
    if (e$summary$completeness_status == 'complete') expect_identical(r$complete_count, length(e$included_units)) else expect_null(r$complete_count)
  }
  bad <- phase_c_case()$result$contexts[[1L]]; bad$summary$complete_count <- 99L
  expect_error(Neoantimon2:::indirect_reconciliation(bad), 'scientific contract must be reviewed')
})
test_that('supported nonbinding at exactly 1000 nM stays separate from missingness', {
  e <- phase_c_detail('nonbinder'); d <- Filter(function(d) d$state == 'excluded_nonbinder', e$decisions)
  expect_length(d, 1L); expect_identical(d[[1L]]$affinity_nM, 1000)
  expect_false(d[[1L]]$affinity_under_frozen_policy); expect_identical(d[[1L]]$outcome_kind, 'provider_not_binding_under_policy')
  expect_identical(e$summary$complete_count, 4L); expect_identical(e$summary$binding_policy_excluded_count, 1L)
})
test_that('peptide occurrences binding origins included and excluded rows order by semantic IDs', {
  x <- phase_c_case('union'); y <- x$result
  for (k in c('pmhcs','comparisons','baseline_units','presenters','self_searches')) y$contexts[[1L]][[k]] <- rev(y$contexts[[1L]][[k]])
  y$contexts[[1L]]$generation$occurrences <- rev(y$contexts[[1L]]$generation$occurrences)
  y$contexts[[1L]]$generation$ledger <- rev(y$contexts[[1L]]$generation$ledger)
  y$provider$rows <- rev(y$provider$rows)
  a <- x$foundation; a$reference$records <- rev(a$reference$records); a$source_occurrences <- rev(a$source_occurrences); a$typing <- rev(a$typing)
  expect_identical(Neoantimon2:::indirect_public_result(y, a), x$public_result)
  expect_identical(Neoantimon2:::indirect_public_explanation(y, a), x$public_explanation)
  for (table in x$public_explanation$contexts[[1L]]$graph$nodes) {
    ids <- vapply(table, `[[`, '', 'id'); expect_identical(ids, sort(ids, method = 'radix'))
  }
})
test_that('actual input row order and physical provenance do not change public semantic content', {
  x <- phase_c_case(); input <- indirect_input()
  m <- indirect_fixture('binding.json'); m$rows <- rev(m$rows)
  for (i in seq_along(m$rows)) {
    m$rows[[i]]$raw_output$path <- '/physical/alias'; m$rows[[i]]$provenance$timestamp <- '2099-01-01'
  }
  m$provenance$metadata <- list(path = '/physical/provider/alias')
  r <- indirect_fixture('reference.json'); r$records <- rev(r$records); r$provenance$retrieved_at <- '2099-01-01'
  input$a$typing <- rev(input$a$typing); input$a$subjects <- rev(input$a$subjects)
  y <- indirect_test(indirect_with_binding(m, indirect_with_reference(r, input)))
  expect_identical(Neoantimon2:::indirect_public_result(y$result, y$foundation), x$public_result)
  e <- Neoantimon2:::indirect_public_explanation(y$result, y$foundation)
  expect_identical(phase_c_semantic(e), phase_c_semantic(x$public_explanation))
  expect_false(identical(e$provenance, x$public_explanation$provenance))
})
test_that('audit target and dimension row order cannot change non-biological audit IDs', {
  x <- Neoantimon2:::indirect_compatibility(); y <- x; y$targets <- rev(y$targets)
  for (i in seq_along(y$targets)) y$targets[[i]]$dimensions <- rev(y$targets[[i]]$dimensions)
  expect_identical(Neoantimon2:::indirect_audit_projection(y), x)
})
test_that('versioned public result schema accepts complete incomplete and conditional fixtures', {
  valid <- phase_c_schema('result')
  for (name in c('complete','incomplete_self','incomplete_scope','accessory','provider_failure','union','missing_source')) expect_true(valid(phase_c_case(name)$public_result))
  for (name in c('complete','incomplete_self','incomplete_scope','accessory','provider_failure')) {
    expected <- phase_c_json_fixture(paste0(name, '-result.json'))
    expect_identical(jsonlite::toJSON(phase_c_case(name)$public_result, auto_unbox = TRUE, null = 'null', digits = NA), jsonlite::toJSON(expected, auto_unbox = TRUE, null = 'null', digits = NA))
  }
})
test_that('public schema rejects false zero incomplete counts unknown states and claims', {
  valid <- phase_c_schema('result'); x <- phase_c_case('incomplete_self')$public_result
  x$contexts[[1L]]$complete_count <- 0L; expect_false(valid(x))
  x <- phase_c_case()$public_result; x$contexts[[1L]]['complete_count'] <- list(NULL); expect_false(valid(x))
  x <- phase_c_case()$public_result; x$contexts[[1L]]$reason_codes <- NULL; expect_false(valid(x))
  x <- phase_c_case()$public_result; x$contexts[[1L]]$analysis_status <- 'negative'; expect_false(valid(x))
  x <- phase_c_case()$public_result; x$official_reproduction <- TRUE; expect_false(valid(x))
})
test_that('explain schema validates all typed graphs and rejects malformed decisions', {
  valid <- phase_c_schema('explanation')
  for (name in c('complete','incomplete_self','incomplete_scope','accessory','provider_failure','union','missing_universe','reference_failure')) expect_true(valid(phase_c_case(name)$public_explanation))
  x <- phase_c_case()$public_explanation; x$contexts[[1L]]$decisions[[1L]]$state <- 'negative'; expect_false(valid(x))
  x <- phase_c_case()$public_explanation; x$contexts[[1L]]$graph$nodes$pmhc[[1L]]$id <- 'orphan'; expect_false(valid(x))
})
test_that('compatibility schema requires explicit non-official claim and open gaps', {
  valid <- phase_c_schema('compatibility'); x <- Neoantimon2:::indirect_compatibility()
  expect_true(valid(x))
  y <- x; y$targets[[1L]]$official_reproduction <- TRUE; expect_false(valid(y))
  y <- x; y$targets[[1L]]$compatible <- TRUE; expect_false(valid(y))
  y <- x; y$targets[[1L]]$dimensions[[1L]]$status <- 'mostly_compatible'; expect_false(valid(y))
  y <- x; i <- which(vapply(y$targets[[1L]]$dimensions, function(d) d$classification == 'U', TRUE))[[1L]]
  y$targets[[1L]]$dimensions[[i]]$status <- 'matched'; expect_false(valid(y))
})
test_that('text explanation reports null honestly with profile and reasons', {
  t <- neo_explain(phase_c_project('incomplete_self'), format = 'text')
  expect_match(t, 'complete_count=unavailable \\(NULL\\)'); expect_match(t, 'inspectable conditional result only')
  expect_match(t, 'self_search_incomplete'); expect_match(t, 'not official PIRCHE reproduction')
  expect_match(t, 'neo-literature-derived-indirect-hla-ii-dr-strict-v1')
})
test_that('exact approved public exports remain and no clinical interpretation fields exist', {
  expect_identical(sort(getNamespaceExports('Neoantimon2')), sort(current_public_api))
  x <- phase_c_case()$public_result
  expect_false(any(c('PIRCHE_score','pirche_score','clinical_risk','rejection_probability','donor_accepted','transplant_suitable') %in% names(x$contexts[[1L]])))
})
