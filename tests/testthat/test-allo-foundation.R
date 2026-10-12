test_that("Phase A validates SOT and both distinct frozen HSCT orientations", {
  roles <- allo_fixture("roles.json")
  ids <- character()
  for (k in c("sot", "hsct_gvh", "hsct_hvg")) {
    x <- allo_test_input(FALSE); x$contexts <- list(roles[[k]]); x$self_search <- list()
    a <- allo_test_build(x); c <- a$contexts[[1L]]
    expect_identical(c$semantics, roles[[k]])
    expect_identical(c$frozen_projection$patient_id, "patient-1")
    expect_identical(c$semantics$presenter_subject_id, c$semantics$self_subject_id)
    ids <- c(ids, c$allo_context_id)
  }
  expect_length(unique(ids), 3L)
  expect_match(ids[[1L]], "^allo-context-id-v1:")
})
test_that("contradictions missing roles and ambiguous orientation are rejected", {
  roles <- allo_fixture("roles.json")
  for (k in c("contradictory", "missing_role")) {
    x <- allo_test_input(); x$contexts <- list(roles[[k]])
    expect_error(allo_test_build(x), "Contradictory|missing required")
  }
  x <- allo_test_input(); x$contexts[[1L]]$direction <- "ambiguous"
  expect_error(allo_test_build(x), "Contradictory")
  x <- allo_test_input(); x$contexts[[1L]]$scenario_id <- "mixed_chimerism"
  expect_error(allo_test_build(x), "deferred")
  x <- allo_test_input(); x$contexts[[1L]]$source_subject_role <- "recipient"
  expect_error(allo_test_build(x), "Contradictory transplant role")
  x <- allo_test_input(); x$contexts[[1L]]$transplant_recipient_subject_id <- "donor-1"
  expect_error(allo_test_build(x), "Invalid allo subject")
  x <- allo_test_input(); x$contexts[[1L]]$claim_level <- "observed_immune_response"
  expect_error(allo_test_build(x), "hypothetical")
})
test_that("exact two-field and sequence-equivalent inputs never impute genotype", {
  exact <- allo_one_resolution(" a*02:01 ")
  expect_identical(exact$evaluability$status, "evaluated")
  expect_identical(exact$resolved_allele, "HLA-A*02:01")
  equivalent <- allo_one_resolution("A*24:02")
  expect_identical(equivalent$evaluability$status, "evaluated")
  expect_null(equivalent$resolved_allele)
  expect_length(equivalent$compatible_alleles, 2L)
  expect_length(equivalent$sequence_ids, 1L)
  expect_identical(equivalent$resolution_operation, "pinned_sequence_equivalence_without_genotype_imputation")
  expect_identical(allo_one_resolution("A*24:02:01")$resolved_allele, "HLA-A*24:02:01")
})
test_that("ambiguity insufficient typing and group labels remain fail closed", {
  ambiguous <- allo_one_resolution("B*07:02")
  expect_identical(ambiguous$evaluability$status, "ambiguous")
  expect_identical(ambiguous$evaluability$reason_codes, list("hla_typing_ambiguous"))
  expect_length(ambiguous$sequence_ids, 2L)
  expect_null(ambiguous$resolved_allele)
  for (allele in c("A*02", "A*02:01G", "A*02:01P")) {
    t <- allo_one_resolution(allele)
    expect_identical(t$evaluability$status, "insufficient")
    expect_null(t$resolved_allele)
    if (grepl("[GP]$", allele)) expect_length(t$compatible_alleles, 0L)
  }
  expect_identical(allo_one_resolution("A*02:01", complete = FALSE)$evaluability$status, "insufficient")
  expect_error(allo_one_resolution("A*02:01~A*24:02"), "single-chain")
})
test_that("unsupported missing partial and provider failure are never negative", {
  expect_identical(allo_one_resolution("B*99:99")$evaluability$status, "unsupported")
  expect_identical(allo_one_resolution("C*04:01")$evaluability$status, "unsupported")
  unavailable <- allo_one_resolution("C*03:04")
  expect_identical(unavailable$evaluability$status, "unavailable")
  expect_null(unavailable$evaluability$value)
  m <- allo_fixture("reference.json")
  m$records[[1L]]$sequence_status <- "partial"
  expect_identical(allo_one_resolution("A*02:01", m)$evaluability$status, "insufficient")
  m$status <- "failed"
  t <- allo_one_resolution("A*02:01", m)
  expect_identical(t$evaluability$status, "failed")
  expect_true("reference_provider_failure" %in% unlist(t$evaluability$reason_codes))
})
test_that("null and low-expression annotations do not assert antigen absence", {
  for (suffix in c("N", "L")) {
    t <- allo_one_resolution(paste0("A*01:01", suffix))
    expect_identical(t$evaluability$status, "evaluated")
    expect_identical(t$expression_annotations, list(suffix))
    expect_length(t$sequence_ids, 1L)
    expect_null(t$evaluability$value)
  }
  expect_identical(allo_one_resolution("A*01:01")$evaluability$status, "ambiguous")
})
test_that("literal residue coverage and coordinate contracts forbid completion", {
  m <- allo_fixture("reference.json"); m$records[[1L]]$protein_sequence <- "MAXLKDPEFGHIRSTWQNYC"
  expect_error(allo_test_build(allo_test_manifest(m)), "Unknown residues|Complete sequence")
  m$records[[1L]]$sequence_status <- "partial"
  m$records[[1L]]$covered_regions <- list(list(name = "known_1", start = 1L, end = 2L), list(name = "known_2", start = 4L, end = 20L))
  a <- allo_test_build(allo_test_manifest(m))
  expect_identical(a$reference$records[[which(vapply(a$reference$records, function(r) r$allele == "HLA-A*02:01", TRUE))]]$unknown_positions, list(3L))
  expect_identical(a$self_universes[[1L]]$evaluability$status, "insufficient_search")
  m$completion_policy <- "homology"
  expect_error(allo_test_build(allo_test_manifest(m)), "completion/imputation forbidden")
  m <- allo_fixture("reference.json"); m$parser$name <- "alignment"
  expect_error(allo_test_build(allo_test_manifest(m)), "alignment")
})
test_that("self envelopes separate coverage readiness from unsupported absence", {
  complete <- allo_test_build(allo_test_input(FALSE))$self_universes[[1L]]
  expect_true(complete$coverage_complete_for_requested_scope)
  expect_false(complete$biological_self_completeness_established)
  expect_identical(complete$evaluability$status, "evaluated")
  expect_identical(complete$comparison_status, "not_evaluated")
  expect_false(complete$search_performed)
  expect_null(complete$evaluability$value)
  x <- allo_test_input(FALSE); x$self_search[[1L]]$catalog_complete <- FALSE
  incomplete <- allo_test_build(x)$self_universes[[1L]]
  expect_false(incomplete$coverage_complete_for_requested_scope)
  expect_identical(incomplete$evaluability$status, "insufficient_search")
  expect_identical(incomplete$comparison_status, "not_evaluated")
  x$self_search[[1L]]$catalog_complete <- TRUE; x$self_search[[1L]]$required_loci <- list("A", "C")
  expect_identical(allo_test_build(x)$self_universes[[1L]]$evaluability$status, "insufficient_search")
  x$self_search[[1L]]$self_subject_id <- "donor-1"
  expect_error(allo_test_build(x), "binding conflict")
  expect_setequal(Neoantimon2:::allo_comparison_states(), c("supported_nonself_within_scope", "self_match", "insufficient_search", "not_evaluated", "ambiguous", "not_applicable"))
})
test_that("semantic identities ignore physical paths checksums clocks row and key order", {
  x <- allo_test_input(); a <- allo_test_build(x)
  m <- allo_fixture("reference.json"); m$records <- rev(m$records)
  m$provenance$retrieved_at <- "2099-01-01T00:00:00Z"
  y <- allo_test_manifest(m, x); y$typing <- rev(y$typing); y$subjects <- rev(y$subjects)
  y$contexts[[1L]] <- rev(y$contexts[[1L]])
  y$contexts[[1L]]$provenance <- list(path = "/other/physical", sha256 = paste(rep("a", 64), collapse = ""), sampling_time = "2099")
  b <- allo_test_build(y)
  expect_identical(a$allo_analysis_id, b$allo_analysis_id)
  expect_identical(lapply(a$contexts, `[[`, "allo_context_id"), lapply(b$contexts, `[[`, "allo_context_id"))
  expect_identical(a$reference$reference_policy_id, b$reference$reference_policy_id)
  expect_identical(lapply(a$self_universes, `[[`, "self_universe_id"), lapply(b$self_universes, `[[`, "self_universe_id"))
  expect_identical(lapply(a$source_occurrences, `[[`, "hla_source_occurrence_id"), lapply(b$source_occurrences, `[[`, "hla_source_occurrence_id"))
  expect_false(identical(a$reference$provenance$artifact, b$reference$provenance$artifact))
  y$typing <- c(y$typing, y$typing[1L])
  expect_identical(b$allo_analysis_id, allo_test_build(y)$allo_analysis_id)
})
test_that("reference release pins policy identity while content identity is reused", {
  x <- allo_test_input(FALSE); a <- allo_test_build(x)
  m <- allo_fixture("reference.json"); m$reference_release <- "fixture-2"
  b <- allo_test_build(allo_test_manifest(m, x))
  expect_false(identical(a$reference$reference_policy_id, b$reference$reference_policy_id))
  expect_false(identical(a$allo_analysis_id, b$allo_analysis_id))
  expect_identical(lapply(a$reference$records, `[[`, "hla_source_sequence_id"), lapply(b$reference$records, `[[`, "hla_source_sequence_id"))
  expect_identical(a$self_universes[[1L]]$self_universe_id, b$self_universes[[1L]]$self_universe_id)
})
test_that("bound references reject mutation and derived overlay tampering", {
  x <- allo_test_manifest(allo_fixture("reference.json")); cfg <- allo_test_config(x)
  p <- neo_project(config = neo_config(env = character(), allo = cfg))
  expect_silent(Neoantimon2:::allo_read(p))
  p$allo_foundation$analysis_status <- "evaluated_negative"
  expect_error(Neoantimon2:::allo_read(p), "derived contract")
  writeLines("changed", cfg$reference$path)
  expect_error(Neoantimon2:::allo_config(cfg), "integrity mismatch")
})
test_that("opt-in config project resume and accessors preserve no-allo contracts", {
  cfg <- neo_config(env = character(), paths = list(run_root = tempdir())); p <- neo_project(config = cfg, fixture_kind = "synthetic")
  expect_null(p$config$allo)
  expect_false("allo_foundation" %in% names(p))
  after <- neo_project(config = neo_config(env = character(), paths = cfg$paths, allo = allo_test_input(FALSE)), fixture_kind = "synthetic")
  expect_identical(p$analysis_id, after$analysis_id)
  expect_identical(p$config_fingerprint, after$config_fingerprint)
  expect_identical(p$input_fingerprint, after$input_fingerprint)
  after <- neo_run(after)
  expect_identical(neo_project(path = after$path)$allo_foundation, after$allo_foundation)
  expect_identical(neo_explain(after)$allo_foundation, after$allo_foundation)
  expect_match(neo_explain(after, format = "text"), "analysis=not_evaluated")
  expect_identical(attr(neo_results(after), "allo_foundation")$analysis_status, "not_evaluated")
  p <- neo_run(p)
  expect_error(neo_explain(p), "At least one")
  expect_null(attr(neo_results(p), "allo_foundation"))
  expect_setequal(getNamespaceExports("Neoantimon2"), current_public_api)
})
test_that("Phase A E2E preserves full frozen somatic source longitudinal and vaccine state", {
  e <- vaccine_real_e2e(); before <- allo_frozen_surface(e$after)
  p <- allo_test_enable(e$after)
  after <- allo_frozen_surface(p)
  expect_true(nrow(before$results) > 0L)
  expect_true(any(before$results$final_grade == "NOT_GRADED"))
  attr(after$results, "allo_foundation") <- NULL
  after$explanation$allo_foundation <- NULL
  expect_identical(after, before)
  expect_identical(p$analysis_id, e$after$analysis_id)
  expect_identical(p$config_fingerprint, e$after$config_fingerprint)
  expect_identical(p$allo_foundation$readiness$status, "ambiguous")
  expect_identical(p$allo_foundation$analysis_status, "not_evaluated")
  expect_identical(p$allo_foundation$self_universes[[1L]]$comparison_status, "not_evaluated")
  expect_false(any(grepl("indirect-peptide|indirect-pmhc", unlist(p$allo_foundation))))
})
test_that("explicit alternatives and Class II chains preserve scoped hypotheses", {
  x <- allo_test_input(FALSE)
  t <- allo_test_typing("B*07:02:01"); t$alleles <- list("B*07:02:01", "B*07:02:02")
  x$typing <- list(t); x$self_search <- list()
  a <- allo_test_build(x)$typing[[1L]]
  expect_identical(a$evaluability$status, "ambiguous")
  expect_length(a$normalized_alleles, 2L)
  expect_null(a$resolved_allele)
  t <- allo_test_typing("A*02:01"); t$alleles <- list("A*02:01", "A*99:99N")
  x$typing <- list(t)
  expect_identical(allo_test_build(x)$typing[[1L]]$evaluability$status, "unsupported")
  for (chain in c("DRA*01:01", "DRB1*04:01")) {
    a <- allo_one_resolution(chain)
    expect_identical(a$hla_class, "II")
    expect_identical(a$evaluability$status, "evaluated")
    expect_false("mhc_key" %in% names(a))
  }
})
test_that("missing self envelopes and identical HLA never establish zero or foreignness", {
  x <- allo_test_input(FALSE); a <- allo_test_build(x)
  expect_identical(a$readiness$status, "evaluated")
  expect_identical(a$analysis_status, "not_evaluated")
  expect_null(a$readiness$value)
  expect_false(any(c("count", "score", "risk", "foreign") %in% names(a)))
  expect_identical(a$contexts[[1L]]$semantics$source_subject_id, "donor-1")
  expect_identical(a$contexts[[1L]]$semantics$self_subject_id, "patient-1")
  x$self_search <- list()
  expect_identical(allo_test_build(x)$readiness$status, "not_evaluated")
  x$typing <- list()
  expect_length(allo_test_build(x)$source_occurrences, 0L)
  expect_identical(allo_test_build(x)$analysis_status, "not_evaluated")
})
test_that("Phase A is reconstructed by the public constructor with frozen identities", {
  p <- reconstruction_test_project()
  cfg <- p$config; cfg$allo <- allo_test_input(FALSE)
  q <- neo_project(samples = p$samples, hla = if (is.null(p$hla)) NULL else lapply(p$hla, unlist), config = cfg,
    analysis_context = p$analysis_context, inputs = recognition_clone_inputs(p), fixture_kind = p$fixture_kind)
  expect_identical(p$analysis_id, q$analysis_id)
  expect_identical(p$config_fingerprint, q$config_fingerprint)
  expect_identical(p$input_fingerprint, q$input_fingerprint)
  expect_identical(p$analysis_context, q$analysis_context)
  expect_identical(p$hla, q$hla)
  expect_match(q$allo_foundation$allo_analysis_id, "^allo-analysis-id-v1:")
})
test_that("equivalent reference aliases share one source occurrence with preserved provenance", {
  x <- allo_test_input(FALSE); m <- allo_fixture("reference.json")
  a <- allo_test_build(x)
  alias <- m$records[[1L]]; alias$sequence_identifier <- "alternate_sequence_alias"
  alias$resolution_provenance <- "independent_alias_mapping"
  m$records <- c(m$records, list(alias))
  b <- allo_test_build(allo_test_manifest(m, x))
  expect_length(b$source_occurrences, 1L)
  expect_identical(a$source_occurrences[[1L]]$hla_source_occurrence_id, b$source_occurrences[[1L]]$hla_source_occurrence_id)
  expect_length(b$source_occurrences[[1L]]$resolution_provenance, 2L)
  expect_identical(allo_one_resolution("A*02:01", m)$evaluability$status, "evaluated")
})
test_that("conflicting completeness for the same source hypothesis is rejected", {
  x <- allo_test_input(FALSE); alternative <- x$typing[[1L]]; alternative$complete <- FALSE
  x$typing <- c(x$typing, list(alternative))
  expect_error(allo_test_build(x), "Conflicting source occurrence evaluability")
})
test_that("subject bindings normalize NFC before validation and semantic hashing", {
  make <- function(subject) {
    x <- allo_test_input(FALSE); x$subjects[[1L]]$subject_id <- subject
    x$contexts[[1L]]$transplant_donor_subject_id <- subject
    x$contexts[[1L]]$source_subject_id <- subject; x$typing[[1L]]$subject_id <- subject
    allo_test_build(x)
  }
  a <- make("donor-\u00e9"); b <- make("donor-e\u0301")
  expect_identical(a$allo_analysis_id, b$allo_analysis_id)
  expect_identical(a$contexts[[1L]]$allo_context_id, b$contexts[[1L]]$allo_context_id)
  expect_identical(a$source_occurrences[[1L]]$hla_source_occurrence_id, b$source_occurrences[[1L]]$hla_source_occurrence_id)
})
test_that("per-record provenance cannot reorder the semantic reference pin", {
  x <- allo_test_input(FALSE); m <- allo_fixture("reference.json")
  alias <- m$records[[1L]]; alias$sequence_identifier <- "another_alias"
  alias$resolution_provenance <- "zzzz_provenance"
  m$records[[1L]]$resolution_provenance <- "aaaa_provenance"
  m$records <- c(m$records, list(alias)); a <- allo_test_build(allo_test_manifest(m, x))
  m$records[[1L]]$resolution_provenance <- "zzzz_provenance"
  m$records[[length(m$records)]]$resolution_provenance <- "aaaa_provenance"
  b <- allo_test_build(allo_test_manifest(m, x))
  expect_identical(a$reference$reference_policy_id, b$reference$reference_policy_id)
  expect_identical(a$allo_analysis_id, b$allo_analysis_id)
})
test_that("source occurrence identity excludes alternative-list and copy/phase evidence", {
  x <- allo_test_input(FALSE)
  x$typing[[1L]]$alleles <- list("A*02:01", "A*24:02")
  a <- allo_test_build(x)
  occurrence <- function(a) Filter(function(o) o$semantics$allele_hypothesis == "HLA-A*02:01", a$source_occurrences)[[1L]]
  id <- occurrence(a)$hla_source_occurrence_id
  changes <- list(list("A*24:02", "A*02:01"), list("A*02:01", "A*24:02", "A*02:01"), list("A*02:01"))
  for (alleles in changes) {
    y <- x; y$typing[[1L]]$alleles <- alleles; b <- allo_test_build(y)
    expect_identical(occurrence(b)$hla_source_occurrence_id, id)
    expect_identical(occurrence(b)$typing_evidence[[1L]]$raw_alleles, alleles)
  }
  for (field in c("copy_number", "phase")) {
    y <- x; y$typing[[1L]][[field]] <- if (field == "copy_number") 2L else "haplotype:2"
    b <- allo_test_build(y)
    expect_identical(occurrence(b)$hla_source_occurrence_id, id)
    expect_false(identical(a$allo_analysis_id, b$allo_analysis_id))
    expect_identical(occurrence(b)$typing_evidence[[1L]][[field]], y$typing[[1L]][[field]])
  }
  expect_false(any(c("normalized_typing", "copy_number", "phase") %in% names(occurrence(a)$semantics)))
})
test_that("self universe identity depends on content not typing resolution evidence", {
  x <- allo_test_input(FALSE)
  own <- which(vapply(x$typing, function(t) t$subject_id == "patient-1", TRUE))
  x$typing[[own]]$alleles <- list("A*24:02:01")
  a <- allo_test_build(x)
  x$typing[[own]]$alleles <- list("A*24:02")
  b <- allo_test_build(x)
  expect_identical(a$self_universes[[1L]]$self_universe_id, b$self_universes[[1L]]$self_universe_id)
  expect_identical(a$self_universes[[1L]]$semantics, b$self_universes[[1L]]$semantics)
  expect_false(identical(a$self_universes[[1L]]$typing_evidence, b$self_universes[[1L]]$typing_evidence))
  expect_false(identical(a$allo_analysis_id, b$allo_analysis_id))
  x$self_search[[1L]]$catalog_complete <- FALSE; c <- allo_test_build(x)
  expect_identical(b$self_universes[[1L]]$self_universe_id, c$self_universes[[1L]]$self_universe_id)
  expect_identical(c$self_universes[[1L]]$evaluability$status, "insufficient_search")
  expect_false(identical(b$allo_analysis_id, c$allo_analysis_id))
})
test_that("duplicate ordered raw self sets and typing alternatives preserve normalized identities", {
  valid <- allo_test_validator()
  x <- allo_test_input(); a <- allo_test_build(x)
  x$self_search[[1L]]$required_loci <- c(rev(x$self_search[[1L]]$required_loci), list("A"))
  x$self_search[[1L]]$required_domains <- rep(x$self_search[[1L]]$required_domains, 2L)
  for (i in seq_along(x$typing)) x$typing[[i]]$alleles <- c(rev(x$typing[[i]]$alleles), x$typing[[i]]$alleles)
  bound <- allo_test_config(x); b <- Neoantimon2:::allo_foundation(bound)
  expect_true(valid(x)); expect_true(valid(bound))
  expect_identical(bound$self_search, x$self_search)
  expect_identical(bound$typing, x$typing)
  expect_identical(a$allo_analysis_id, b$allo_analysis_id)
  expect_identical(a$self_universes[[1L]]$semantics, b$self_universes[[1L]]$semantics)
  expect_identical(lapply(a$source_occurrences, `[[`, "hla_source_occurrence_id"), lapply(b$source_occurrences, `[[`, "hla_source_occurrence_id"))
})
test_that("present optional null and invalid provenance/reason types agree with raw schemas", {
  input_valid <- allo_test_validator(); manifest_valid <- allo_test_validator("hla-reference", "manifest")
  x <- allo_test_input(FALSE); x$contexts[[1L]]["provenance"] <- list(NULL)
  expect_false(input_valid(x)); expect_error(allo_test_build(x), "null field")
  for (field in c("sha256", "manifest")) {
    x <- allo_test_input(FALSE); x$reference[field] <- list(NULL)
    expect_false(input_valid(x)); expect_error(allo_test_build(x), "null field")
  }
  for (value in list(NULL, "invalid", list(), list(path = NULL), list(path = 3L))) {
    m <- allo_fixture("reference.json"); m["provenance"] <- list(value)
    expect_false(manifest_valid(m)); expect_error(allo_test_build(allo_test_manifest(m)))
  }
  for (value in list(NULL, "invalid", list(3L))) {
    m <- allo_fixture("reference.json"); m["reason_codes"] <- list(value)
    expect_false(manifest_valid(m)); expect_error(allo_test_build(allo_test_manifest(m)))
  }
  for (status in c("evaluated", "failed")) {
    m <- allo_fixture("reference.json"); m$status <- status; m$reason_codes <- list("declared_reason")
    expect_true(manifest_valid(m)); expect_silent(allo_test_build(allo_test_manifest(m)))
  }
  m <- allo_fixture("reference.json"); m$provenance <- stats::setNames(list(), character())
  expect_true(manifest_valid(m)); expect_silent(allo_test_build(allo_test_manifest(m)))
  x <- allo_test_input(FALSE); x$contexts[[1L]]$provenance <- list(path = "declared")
  expect_true(input_valid(x)); expect_silent(allo_test_build(x))
})
test_that("declared nullable typing and coordinate fields remain valid", {
  valid <- allo_test_validator(); mv <- allo_test_validator("hla-reference", "manifest")
  x <- allo_test_input(FALSE)
  for (i in seq_along(x$typing)) x$typing[[i]][c("copy_number", "phase")] <- list(NULL, NULL)
  expect_true(valid(x)); expect_silent(allo_test_build(x))
  m <- allo_fixture("reference.json")
  expect_true(mv(m)); expect_silent(allo_test_build(allo_test_manifest(m)))
  m$sequence_policy <- "mature"
  for (i in seq_along(m$records)) m$records[[i]]$coordinate_map <- list(sequence_coordinates = "mature", precursor_start = NULL, mature_start = 1L)
  x$self_search[[1L]]$sequence_policy <- "mature"
  expect_true(mv(m)); expect_silent(allo_test_build(allo_test_manifest(m, x)))
})
test_that("both coordinate starts reject zero negative fractional and wrong origin values", {
  valid <- allo_test_validator("hla-reference", "manifest")
  for (field in c("precursor_start", "mature_start")) for (value in c(0L, -1L, 1.5)) {
    m <- allo_fixture("reference.json"); m$records[[1L]]$coordinate_map[[field]] <- value
    expect_false(valid(m)); expect_error(allo_test_build(allo_test_manifest(m)), "Invalid coordinate map")
  }
  m <- allo_fixture("reference.json"); m$records[[1L]]$coordinate_map$precursor_start <- 2L
  expect_false(valid(m)); expect_error(allo_test_build(allo_test_manifest(m)), "coordinate origin")
  m$records[[1L]]$coordinate_map$precursor_start <- 1L; m$records[[1L]]$coordinate_map$mature_start <- 2L
  expect_true(valid(m)); expect_silent(allo_test_build(allo_test_manifest(m)))
})
test_that("source typing scenarios fail closed without selecting or unioning", {
  x <- allo_test_input(); a <- allo_test_build(x)
  expect_true(all(c("A", "B") %in% vapply(Filter(function(t) t$subject_id == "donor-1", a$typing), `[[`, "", "locus")))
  expect_true(all(vapply(a$source_occurrences, function(o) o$semantics$typing_scenario_id == "typing:1", TRUE)))
  other <- x$typing[[1L]]; other$typing_scenario_id <- "typing:alternative"
  x$typing <- c(x$typing, list(other))
  expect_error(allo_test_build(x), "Multiple source typing scenarios unsupported")
  x$typing <- rev(x$typing)
  expect_error(allo_test_build(x), "Multiple source typing scenarios unsupported")
  x <- allo_test_input(FALSE); other <- x$typing[[2L]]; other$typing_scenario_id <- "typing:alternative"
  x$typing <- c(x$typing, list(other)); b <- allo_test_build(x)
  expect_identical(b$self_universes[[1L]]$typing_scenario_id, "typing:1")
  expect_true(all(vapply(b$self_universes[[1L]]$typing_evidence, function(t) t$typing_scenario_id == "typing:1", TRUE)))
})
test_that("context identity uses the exact frozen payload and retains claim metadata", {
  x <- allo_test_input(FALSE); a <- allo_test_build(x); context <- a$contexts[[1L]]
  fields <- c("transplant_donor_subject_id", "transplant_recipient_subject_id",
    "source_subject_id", "self_subject_id", "presenter_subject_id", "responding_subject_id",
    "source_subject_role", "self_subject_role", "presenter_subject_role", "responding_subject_role",
    "transplant_type", "scenario_id", "direction", "tissue", "apc_compartment", "occasion_id")
  expected <- x$contexts[[1L]][fields]
  namespace <- "allo-context-id-v1"
  id <- paste0(namespace, ":", Neoantimon2:::infra_hash(list(namespace = namespace, payload = expected)))
  payload <- Neoantimon2:::allo_context_identity_payload(context$semantics)
  expect_identical(payload, expected)
  expect_identical(context$allo_context_id, id)
  expect_identical(context$semantics$claim_level, "hypothetical_molecular_scenario")
  expect_false("claim_level" %in% names(payload))
  future <- context$semantics; future$future_record_metadata <- "inspectable_only"
  expect_identical(Neoantimon2:::allo_context_identity_payload(future), expected)
  x$contexts[[1L]]$claim_level <- "observed_immune_response"
  expect_error(allo_test_build(x), "hypothetical molecular scenarios")
})
test_that("complete and partial status share molecular identity but preserve distinct evidence", {
  x <- allo_test_input(FALSE); m <- allo_fixture("reference.json")
  a <- allo_test_build(allo_test_manifest(m, x))
  m$records[[1L]]$sequence_status <- "partial"
  b <- allo_test_build(allo_test_manifest(m, x))
  sequence <- function(a) Filter(function(r) r$allele == "HLA-A*02:01", a$reference$records)[[1L]]
  r <- sequence(a); s <- sequence(b)
  payload <- Neoantimon2:::allo_sequence_identity_payload(a$reference$semantics$species, r)
  expect_identical(payload, Neoantimon2:::allo_sequence_identity_payload(b$reference$semantics$species, s))
  expect_false("sequence_status" %in% names(payload))
  expect_identical(r$hla_source_sequence_id, s$hla_source_sequence_id)
  expect_identical(r$sequence_status, "complete")
  expect_identical(s$sequence_status, "partial")
  expect_false(identical(a$reference$reference_policy_id, b$reference$reference_policy_id))
  expect_false(identical(a$allo_analysis_id, b$allo_analysis_id))
  expect_true(all(vapply(a$typing, function(t) t$evaluability$status == "evaluated", TRUE)))
  expect_true(all(vapply(b$typing, function(t) t$evaluability$status == "insufficient", TRUE)))
  expect_true(all(vapply(b$typing, function(t) "source_sequence_partial" %in% unlist(t$evaluability$reason_codes), TRUE)))
  expect_identical(a$readiness$status, "evaluated")
  expect_identical(b$readiness$status, "insufficient")
  expect_identical(a$source_occurrences[[1L]]$hla_source_occurrence_id, b$source_occurrences[[1L]]$hla_source_occurrence_id)
  valid <- allo_test_validator("hla-reference", "manifest")
  complete <- allo_fixture("reference.json")
  expect_true(valid(complete)); expect_true(valid(m))
  future <- r; future$future_record_metadata <- "inspectable_only"
  expect_identical(Neoantimon2:::allo_sequence_identity_payload(a$reference$semantics$species, future), payload)
})
test_that("four frozen namespace payload contracts cannot silently widen", {
  x <- allo_test_input(FALSE); a <- allo_test_build(x)
  context <- a$contexts[[1L]]; occurrence <- a$source_occurrences[[1L]]; universe <- a$self_universes[[1L]]
  record <- Filter(function(r) r$allele == "HLA-A*02:01", a$reference$records)[[1L]]
  frozen_id <- function(namespace, payload) paste0(namespace, ":",
    Neoantimon2:::infra_hash(list(namespace = namespace, payload = payload)))
  context_payload <- list(transplant_donor_subject_id = "donor-1", transplant_recipient_subject_id = "patient-1",
    source_subject_id = "donor-1", self_subject_id = "patient-1", presenter_subject_id = "patient-1",
    responding_subject_id = "patient-1", source_subject_role = "donor", self_subject_role = "recipient",
    presenter_subject_role = "recipient", responding_subject_role = "recipient", transplant_type = "SOT",
    scenario_id = "sot_recipient_response", direction = "donor_to_recipient", tissue = "declared_test_tissue",
    apc_compartment = "declared_test_apc", occasion_id = "occasion:phase-a-1")
  sequence_payload <- list(species = "NCBITaxon:9606", amino_acids = record$protein_sequence,
    unknown_positions = record$unknown_positions, protein_definition = record$protein_definition,
    isoform = record$isoform, coordinate_map = record$coordinate_map, domains = record$domains,
    covered_regions = record$covered_regions)
  occurrence_payload <- list(source_subject_id = "donor-1", allele_hypothesis = "HLA-A*02:01",
    sequence_id = record$hla_source_sequence_id, domains = record$domains, coordinate_map = record$coordinate_map,
    typing_scenario_id = "typing:1")
  universe_payload <- list(self_subject_id = "patient-1", scenario = context$allo_context_id,
    scope = "all_self_hla", required_loci = list("A"), required_domains = list("full_precursor"),
    isoform_policy = "explicit_mapped_isoforms", sequence_policy = "precursor",
    sequences = list(list(sequence_id = record$hla_source_sequence_id, domains = record$domains,
      covered_regions = record$covered_regions)))
  expect_identical(Neoantimon2:::allo_context_identity_payload(context$semantics), context_payload)
  expect_identical(context$allo_context_id, frozen_id("allo-context-id-v1", context_payload))
  expect_identical(Neoantimon2:::allo_sequence_identity_payload(a$reference$semantics$species, record), sequence_payload)
  expect_identical(record$hla_source_sequence_id, frozen_id("hla-source-sequence-id-v1", sequence_payload))
  expect_identical(occurrence$semantics, occurrence_payload)
  expect_identical(occurrence$hla_source_occurrence_id, frozen_id("hla-source-occurrence-id-v1", occurrence_payload))
  expect_identical(universe$semantics, universe_payload)
  expect_identical(universe$self_universe_id, frozen_id("self-universe-id-v1", universe_payload))
})
