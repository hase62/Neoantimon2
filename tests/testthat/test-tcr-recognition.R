test_that("AIRR bulk chain counts retain distinct nullable measures", {
  r <- recognition_fixture_repertoire()
  expect_length(r$chains, 4L)
  expect_length(r$profiles, 3L)
  expect_length(r$observations, 4L)
  by_id <- setNames(r$chains, vapply(r$chains, `[[`, "", "sequence_id"))
  expect_identical(by_id$b1$duplicate_count, 10)
  expect_identical(by_id$b1$umi_count, 3)
  expect_identical(by_id$b1$consensus_count, 2)
  expect_null(by_id$missing$duplicate_count)
  expect_identical(by_id$zero$duplicate_count, 0)
  expect_true(by_id$b1$productive)
  expect_false(by_id$b1$stop_codon)
  expect_true(by_id$b1$vj_in_frame)
  expect_null(by_id$b1$airr_schema_version)
  expect_true(all(vapply(r$observations, function(o) is.null(o$receptor_hash), TRUE)))
})
test_that("bulk alpha and beta do not manufacture pairs", {
  r <- Neoantimon2:::tcr_import(recognition_fixture_project(list(recognition_fixture_input("bulk-unpaired.tsv"))))
  expect_length(r$pairing_hypotheses, 0L)
  expect_false(any(vapply(r$profiles, function(p) p$completeness == "paired_cdr3_profile", TRUE)))
  expect_true(all(vapply(r$observations, function(o) o$pairing_status == "unpaired_chain", TRUE)))
})
test_that("observed virtual and multiple chains preserve every hypothesis", {
  r <- recognition_fixture_repertoire(TRUE)
  observed <- Filter(function(o) identical(o$cell_id, "observed"), r$pairing_hypotheses)
  virtual <- Filter(function(o) identical(o$cell_id, "virtual"), r$pairing_hypotheses)
  multi <- Filter(function(o) identical(o$cell_id, "multiple"), r$pairing_hypotheses)
  expect_length(observed, 1L); expect_identical(observed[[1L]]$pairing_status, "observed_pair")
  expect_length(virtual, 1L); expect_identical(virtual[[1L]]$pairing_status, "virtual_or_inferred_pair")
  expect_length(multi, 2L)
  expect_true(all(vapply(multi, function(o) o$pairing_status == "ambiguous_pairing", TRUE)))
  expect_length(unique(vapply(multi, `[[`, "", "receptor_feature_profile_id")), 2L)
  expect_length(r$repertoires, 1L)
  expect_identical(r$repertoires[[1L]]$airr_schema_version, "2.0")
  expect_identical(observed[[1L]]$t_cell_type, "CD8")
  expect_null(virtual[[1L]]$t_cell_type)
})
test_that("full Receptor provider hash remains independent of CDR3 profiles", {
  r <- recognition_fixture_repertoire(TRUE)
  domain <- Filter(function(p) p$completeness == "paired_variable_domain_receptor", r$profiles)
  expect_length(domain, 1L)
  expect_setequal(unlist(domain[[1L]]$available_features), c("alpha_variable_aa", "beta_variable_aa"))
  o <- Filter(function(o) o$receptor_feature_profile_id == domain[[1L]]$receptor_feature_profile_id, r$observations)[[1L]]
  expect_identical(o$receptor_hash, r$receptors[[1L]]$raw$receptor_hash)
  expect_match(o$receptor_hash, "^[a-f0-9]{64}$")
  expect_null(o$receptor_record$raw$junction_aa)
  expect_false("beta_cdr3_aa" %in% unlist(domain[[1L]]$available_features))
})
test_that("receptor profile identity ignores abundance source and row ordering", {
  original <- recognition_fixture_repertoire()
  input <- recognition_fixture_input("bulk.tsv")
  lines <- readLines(input$path)
  lines[-1L] <- sub("\t10\t3\t2\t", "\t500\t600\t700\t", lines[-1L])
  path <- tempfile(fileext = ".tsv"); writeLines(c(lines[[1L]], rev(lines[-1L])), path)
  input$path <- path
  changed <- Neoantimon2:::tcr_import(recognition_fixture_project(list(input)))
  expect_identical(original$profiles, changed$profiles)
  expect_false(identical(original$observations, changed$observations))
  cap <- recognition_test_cap()
  a <- Neoantimon2:::recognition_prepare(list(recognition_test_candidate()), original, cap, recognition_test_settings())
  b <- Neoantimon2:::recognition_prepare(list(recognition_test_candidate()), changed, cap, recognition_test_settings())
  expect_identical(a$queries, b$queries)
})
test_that("paired alpha changes distinguish receptors sharing beta", {
  p <- recognition_paired_profile(); q <- p$features; q$alpha_cdr3_aa <- "CAVSDNYQLIW"
  new <- Neoantimon2:::recognition_profile(q)
  expect_false(identical(p$receptor_feature_profile_id, new$receptor_feature_profile_id))
  expect_identical(p$features$beta_cdr3_aa, new$features$beta_cdr3_aa)
  beta <- Neoantimon2:::recognition_profile(list(beta_cdr3_aa = p$features$beta_cdr3_aa))
  expect_identical(beta$completeness, "beta_cdr3_profile")
  expect_false(identical(beta$receptor_feature_profile_id, p$receptor_feature_profile_id))
})
test_that("same receptor in tumor and normal retains separate observations", {
  p <- recognition_fixture_project(list(recognition_fixture_input("bulk.tsv"),
    recognition_fixture_input("bulk.tsv", id = "normal-bulk", sample = "normal")))
  r <- Neoantimon2:::tcr_import(p)
  expect_length(r$profiles, 3L); expect_length(r$observations, 8L)
  expect_setequal(vapply(r$observations, `[[`, "", "sample_id"), c("sample-1", "normal"))
  expect_length(unique(vapply(r$observations, `[[`, "", "repertoire_observation_id")), 8L)
})
test_that("AIRR named header reordering and unknown fields are safe", {
  input <- recognition_fixture_input("bulk.tsv")
  fields <- lapply(readLines(input$path), function(s) stringi::stri_split_fixed(s, "\t", omit_empty = FALSE)[[1L]])
  path <- tempfile(fileext = ".tsv")
  writeLines(vapply(fields, function(x) paste(rev(x), collapse = "\t"), ""), path)
  input$path <- path
  changed <- Neoantimon2:::tcr_import(recognition_fixture_project(list(input)))
  expect_identical(changed$profiles, recognition_fixture_repertoire()$profiles)
  expect_true(all(vapply(changed$chains, function(c) identical(c$raw_metadata$source, "synthetic"), TRUE)))
})
test_that("malformed required sequences stay unsupported without replacement", {
  cap <- recognition_test_cap("pmtnet"); c <- recognition_test_candidate()
  for (seq in c("CASSX", "CASSU", "CASS*", "cassy", "", "CASS-")) {
    p <- Neoantimon2:::recognition_profile(list(beta_cdr3_aa = seq))
    expect_true(Neoantimon2:::recognition_input_status(c, p, cap) %in% c("unsupported_tcr_character", "incomplete_features"))
    expect_identical(p$features$beta_cdr3_aa, seq)
  }
  expect_identical(Neoantimon2:::recognition_input_status(c, NULL, cap), "incomplete_features")
  expect_error(Neoantimon2:::recognition_number("1.5", TRUE), "Invalid")
  expect_error(Neoantimon2:::recognition_boolean("maybe"), "Invalid")
})
test_that("MHC-conditioned query identity preserves canonical presenter", {
  p <- recognition_paired_profile(); cap <- recognition_test_cap("pmtnet")
  a <- Neoantimon2:::recognition_query(recognition_test_candidate(), p, cap)
  b <- Neoantimon2:::recognition_query(recognition_test_candidate(hla = "HLA-A*01:01"), p, cap)
  expect_false(identical(a$recognition_query_id, b$recognition_query_id))
  expect_identical(a$provider_hla, "HLA-A02:01")
  expect_identical(a$presenter, "HLA-A*02:01")
  ii <- Neoantimon2:::recognition_query(recognition_test_candidate(hla = "HLA-DRA*01:01~HLA-DRB1*04:01"),
    p, recognition_test_cap())
  expect_identical(ii$input_status, "supported")
  expect_identical(ii$presenter, "HLA-DRA*01:01~HLA-DRB1*04:01")
})
test_that("query batch order and size do not alter prediction identities", {
  c <- list(recognition_test_candidate(), recognition_test_candidate(hla = "HLA-A*01:01"))
  r <- recognition_fixture_repertoire(TRUE); cap <- recognition_test_cap()
  a <- Neoantimon2:::recognition_prepare(c, r, cap, recognition_test_settings(batch_size = 1))
  r$profiles <- rev(r$profiles); r$observations <- rev(r$observations)
  b <- Neoantimon2:::recognition_prepare(rev(c), r, cap, recognition_test_settings(batch_size = 100))
  expect_identical(a$queries, b$queries)
  expect_identical(Neoantimon2:::recognition_predict(a, cap)$evidence,
    Neoantimon2:::recognition_predict(b, cap)$evidence)
})
test_that("query and pairing limits refuse without truncating", {
  cap <- recognition_test_cap(); r <- recognition_fixture_repertoire(TRUE)
  a <- Neoantimon2:::recognition_prepare(list(recognition_test_candidate()), r, cap, recognition_test_settings(max_queries = 1))
  expect_identical(a$status, "not_run"); expect_identical(a$reason, "query_limit")
  expect_length(a$queries, 0L); expect_gt(a$requested, 1)
  p <- recognition_fixture_project(list(recognition_fixture_input("paired.tsv"), recognition_fixture_input("airr.json", "tcr_airr")))
  limited <- Neoantimon2:::tcr_import(p, max_hypotheses = 1)
  expect_true(any(vapply(limited$pairing_statuses, function(x) x$status == "not_run/query_limit", TRUE)))
  expect_lte(length(limited$pairing_hypotheses), 1L)
  expect_length(limited$chains, 7L)
})
test_that("PRIME retains per-allele raw metrics without receptor duplication", {
  result <- recognition_test_predict("prime", "prime.tsv")
  expect_length(result$evidence, 3L)
  expect_true(all(vapply(result$evidence, function(e) is.null(e$receptor_feature_profile_id), TRUE)))
  expect_true(all(vapply(result$evidence, function(e) e$scope == "peptide_intrinsic", TRUE)))
  expect_setequal(vapply(result$evidence, `[[`, "", "metric"), c("PRIME.rank", "PRIME.score", "MixMHCpred.rank_in_PRIME_output"))
  expect_identical(result$evidence[[1L]]$metadata$raw_fields$BestAllele, "A0201")
  c <- recognition_test_candidate(hla = "HLA-DRA*01:01~HLA-DRB1*04:01")
  expect_identical(Neoantimon2:::recognition_input_status(c, NULL, recognition_test_cap("prime")), "unsupported_class")
  cap <- recognition_test_cap("prime")
  expect_false(cap$execution_available); expect_true(cap$parser_available)
  expect_setequal(unlist(cap$dependencies), c("PRIME 2.1", "MixMHCpred >=3.0"))
})
test_that("PRIME best-allele-only data cannot leak into another presenter", {
  path <- tempfile(fileext = ".tsv")
  writeLines(c("Peptide\t%Rank_bestAllele\tBestAllele", "GILGFVFTL\t1\tA0201"), path)
  expect_error(Neoantimon2:::recognition_parse(path, recognition_test_cap("prime")), "per-allele")
})
test_that("pMTnet rank remains rank and beta-only consumed features are explicit", {
  result <- recognition_test_predict("pmtnet", "pmtnet.csv")$evidence[[1L]]
  expect_identical(result$rank, 0.01)
  expect_null(result$raw_score); expect_null(result$probability)
  expect_identical(result$score_direction, "lower")
  expect_identical(names(result$metadata$features_used), "beta_cdr3_aa")
  expect_match(result$metadata$ranking_background, "no probability")
  expect_identical(recognition_test_cap("pmtnet")$required_features, list("beta_cdr3_aa"))
  expect_error(recognition_test_cap("pmtnet", "paired"), "beta CDR3 only")
})
test_that("ERGO-II configurations never fall back silently", {
  p <- recognition_paired_profile(); c <- recognition_test_candidate()
  beta <- Neoantimon2:::recognition_profile(list(beta_cdr3_aa = p$features$beta_cdr3_aa))
  expect_identical(Neoantimon2:::recognition_input_status(c, beta, recognition_test_cap("ergo_ii", "beta_only")), "supported")
  expect_identical(Neoantimon2:::recognition_input_status(c, beta, recognition_test_cap("ergo_ii", "paired")), "incomplete_features")
  expect_identical(Neoantimon2:::recognition_input_status(c, p, recognition_test_cap("ergo_ii", "paired")), "supported")
  expect_identical(Neoantimon2:::recognition_input_status(c, p, recognition_test_cap("ergo_ii", "paired_vj")), "supported")
  for (model in c("beta_only", "paired", "paired_vj")) {
    e <- recognition_test_predict("ergo_ii", "ergo.csv", model, use_mhc = FALSE)$evidence[[1L]]
    expect_identical(e$mhc_conditioning, "absent")
    expect_null(e$metadata$presenter)
    expect_identical(e$raw_score, 0.05)
  }
  expect_identical(recognition_test_predict("ergo_ii", "ergo.csv", "paired")$evidence[[1L]]$mhc_conditioning, "present")
})
test_that("NetTCR-2.0 has exact A0201 nine-mer release gates", {
  p <- recognition_paired_profile(); cap <- recognition_test_cap("nettcr_2_0")
  expect_identical(Neoantimon2:::recognition_input_status(recognition_test_candidate(), p, cap), "supported")
  expect_identical(Neoantimon2:::recognition_input_status(recognition_test_candidate("GILGFVFTLA"), p, cap), "unsupported_length")
  expect_identical(Neoantimon2:::recognition_input_status(recognition_test_candidate(hla = "HLA-A*01:01"), p, cap), "unsupported_presenter")
  expect_identical(Neoantimon2:::recognition_input_status(recognition_test_candidate(hla = "HLA-DRA*01:01~HLA-DRB1*04:01"), p, cap), "unsupported_class")
  expect_length(recognition_test_predict("nettcr_2_0", "nettcr.csv", "paired")$evidence, 1L)
  expect_identical(recognition_test_predict("nettcr_2_0", "nettcr.csv")$evidence[[1L]]$mhc_conditioning, "required_context")
})
test_that("NetTCR score does not transfer to another HLA allele", {
  cap <- recognition_test_cap("nettcr_2_0"); r <- recognition_fixture_repertoire(TRUE)
  a <- Neoantimon2:::recognition_prepare(list(recognition_test_candidate(hla = "HLA-A*01:01")), r, cap,
    recognition_test_settings(profile_ids = list(recognition_paired_profile(r)$receptor_feature_profile_id)))
  p <- recognition_fixture_project(list(recognition_fixture_input("nettcr.csv", "recognition_output")))
  result <- Neoantimon2:::recognition_predict(a, cap, p$inputs[[1L]])
  expect_length(result$evidence, 0L)
  expect_identical(result$statuses[[1L]]$reason, "unsupported_presenter")
})
test_that("provider model versions options and consumed sequence fingerprint output", {
  p <- recognition_paired_profile(); c <- recognition_test_candidate()
  a <- Neoantimon2:::recognition_query(c, p, recognition_test_cap())
  for (cap in list(recognition_test_cap(provider_version = "different"), recognition_test_cap(model_version = "different"),
    recognition_test_cap(options = list(background = "different")))) {
    b <- Neoantimon2:::recognition_query(c, p, cap)
    expect_false(identical(a$fingerprint, b$fingerprint))
  }
  file <- tempfile(); writeLines("model-a", file)
  cap_a <- recognition_test_cap(model_files = list(file))
  writeLines("model-b", file)
  cap_b <- recognition_test_cap(model_files = list(file))
  expect_false(identical(Neoantimon2:::recognition_query(c, p, cap_a)$fingerprint,
    Neoantimon2:::recognition_query(c, p, cap_b)$fingerprint))
  expect_identical(recognition_test_cap()$version_status, "unknown")
})
test_that("cross-model conflicting output stays independent", {
  records <- c(recognition_test_predict("pmtnet", "pmtnet.csv")$evidence,
    recognition_test_predict("ergo_ii", "ergo.csv")$evidence,
    recognition_test_predict("nettcr_2_0", "nettcr.csv")$evidence)
  expect_length(records, 3L)
  expect_setequal(vapply(records, `[[`, "", "provider"), c("pmtnet", "ergo_ii", "nettcr_2_0"))
  expect_identical(records[[1L]]$rank, 0.01)
  expect_identical(records[[2L]]$raw_score, 0.05)
  expect_identical(records[[3L]]$raw_score, 0.95)
  expect_true(all(vapply(records, function(e) is.null(e$probability), TRUE)))
  expect_false(any(c("consensus", "winner", "immunogenicity_score") %in% names(records[[1L]])))
})
test_that("experimental positive negative and continuous measurements coexist", {
  p <- recognition_fixture_project(list(recognition_fixture_input("paired.tsv"), recognition_fixture_input("airr.json", "tcr_airr"),
    recognition_fixture_input("reactivity.json", "tcr_reactivity")))
  r <- Neoantimon2:::tcr_import(p)
  e <- Neoantimon2:::recognition_experiments(list(recognition_test_candidate()), r)
  expect_length(e$measurements, 5L); expect_length(e$evidence, 4L)
  measurements <- setNames(e$measurements, vapply(e$measurements, `[[`, "", "reactivity_id"))
  expect_identical(measurements$positive$qualitative_annotation, "positive")
  expect_identical(measurements$negative$qualitative_annotation, "negative")
  expect_identical(measurements$positive$unit, "SFU")
  expect_identical(measurements$continuous$value, 38)
  expect_null(measurements$continuous$threshold)
  expect_identical(measurements$continuous$target_level, "cell")
  expect_identical(measurements$continuous$receptor_assignment, "ambiguous_pairing")
  expect_length(measurements$continuous$possible_receptor_profile_ids, 2L)
  ambiguous <- Filter(function(e) e$metadata$reactivity_id == "continuous", e$evidence)[[1L]]
  expect_null(ambiguous$receptor_feature_profile_id)
  expect_match(measurements$positive$receptor_assignment, "causality_not_inferred")
  expect_identical(measurements$unknown$presenter$status, "unknown")
  expect_identical(measurements$positive$threshold, 10L)
  expect_identical(measurements$positive$external_references, list("synthetic:study"))
})
test_that("experimental and predicted disagreement retains both records", {
  p <- recognition_fixture_project(list(recognition_fixture_input("paired.tsv"), recognition_fixture_input("airr.json", "tcr_airr"),
    recognition_fixture_input("reactivity.json", "tcr_reactivity")))
  e <- Neoantimon2:::recognition_experiments(list(recognition_test_candidate()), Neoantimon2:::tcr_import(p))$evidence
  predictions <- recognition_test_predict("ergo_ii", "ergo.csv")$evidence
  all <- c(e, predictions)
  expect_length(all, 5L)
  expect_true(any(vapply(all, function(e) identical(e$metadata$qualitative_annotation, "positive"), TRUE)))
  expect_true(any(vapply(all, function(e) identical(e$metadata$qualitative_annotation, "negative"), TRUE)))
  expect_true(any(vapply(all, function(e) e$provider == "ergo_ii" && e$raw_score == 0.05, TRUE)))
})
test_that("Class II experiment preserves full heterodimer and unknown phenotype", {
  p <- recognition_fixture_project(list(recognition_fixture_input("paired.tsv"), recognition_fixture_input("airr.json", "tcr_airr"),
    recognition_fixture_input("reactivity.json", "tcr_reactivity")))
  r <- Neoantimon2:::tcr_import(p)
  e <- Neoantimon2:::recognition_experiments(list(recognition_test_candidate(hla = "HLA-DRA*01:01~HLA-DRB1*04:01")), r)
  expect_length(e$evidence, 2L) # Class II assay plus explicitly unknown restriction.
  known <- Filter(function(e) e$mhc_conditioning == "present", e$evidence)[[1L]]
  expect_identical(known$metadata$presenter$mhc_key, "HLA-DRA*01:01~HLA-DRB1*04:01")
  expect_false("inferred_CD4" %in% names(known$metadata))
})
test_that("VDJdb exact paired beta and unknown restriction preserve annotation semantics", {
  p <- recognition_fixture_project(list(recognition_fixture_input("vdjdb.tsv", "tcr_vdjdb")))
  annotations <- Neoantimon2:::recognition_vdjdb(p$inputs[[1L]])
  r <- recognition_fixture_repertoire(TRUE)
  matches <- Neoantimon2:::recognition_curated_matches(list(recognition_test_candidate()), r, annotations)
  expect_gt(length(matches$evidence), 0L)
  expect_true(all(vapply(matches$evidence, function(e) e$scope == "external_annotation" && !e$metadata$patient_assay, TRUE)))
  refs <- vapply(matches$evidence, function(e) e$metadata$annotation$raw_fields$reference.id, "")
  expect_true("synthetic:paired" %in% refs)
  expect_true("synthetic:beta" %in% refs)
  expect_true("synthetic:unknown" %in% refs)
  expect_false("synthetic:conflict" %in% refs)
  expect_true(all(vapply(matches$evidence, function(e) e$provider_version == "synthetic-2.0", TRUE)))
  expect_true(all(vapply(matches$observations, function(o) o$status == "receptor_observed", TRUE)))
  unknown <- Filter(function(e) e$mhc_conditioning == "unknown", matches$evidence)
  expect_gt(length(unknown), 0L)
  expect_match(unknown[[1L]]$metadata$match_status, "restriction_unknown")
  expect_identical(annotations[[1L]]$raw_fields$extra, "synthetic")
  expect_identical(annotations[[1L]]$raw_fields$vdjdb.score, "3")
})
test_that("VDJdb feature matching is exact and row-order invariant", {
  input <- recognition_fixture_input("vdjdb.tsv", "tcr_vdjdb")
  p <- recognition_fixture_project(list(input)); a <- Neoantimon2:::recognition_vdjdb(p$inputs[[1L]])
  lines <- readLines(input$path); path <- tempfile(fileext = ".tsv")
  writeLines(c(lines[[1L]], rev(lines[-1L])), path); input$path <- path
  q <- recognition_fixture_project(list(input)); b <- Neoantimon2:::recognition_vdjdb(q$inputs[[1L]])
  expect_setequal(vapply(a, `[[`, "", "annotation_id"), vapply(b, `[[`, "", "annotation_id"))
  r <- recognition_fixture_repertoire(TRUE)
  x <- Neoantimon2:::recognition_curated_matches(list(recognition_test_candidate()), r, a)
  y <- Neoantimon2:::recognition_curated_matches(list(recognition_test_candidate()), r, b)
  expect_identical(vapply(x$evidence, `[[`, "", "recognition_evidence_id"), vapply(y$evidence, `[[`, "", "recognition_evidence_id"))
  for (i in seq_along(r$profiles)) if (!is.null(r$profiles[[i]]$features$beta_cdr3_aa)) {
    f <- r$profiles[[i]]$features; f$beta_cdr3_aa <- paste0(f$beta_cdr3_aa, "A")
    r$profiles[[i]] <- Neoantimon2:::recognition_profile(f)
  }
  absent <- Neoantimon2:::recognition_curated_matches(list(recognition_test_candidate()), r, a)
  expect_length(absent$evidence, 0L); expect_length(absent$observations, 0L)
})
test_that("repertoire absence and no database match create no negative evidence", {
  p <- recognition_fixture_project(list(recognition_fixture_input("bulk.tsv"), recognition_fixture_input("vdjdb.tsv", "tcr_vdjdb")))
  c <- recognition_test_candidate("AAAAAAAAA")
  x <- Neoantimon2:::recognition_construct(p, list(c))
  expect_length(x$evidence, 0L); expect_length(x$repertoire_matches, 0L)
  expect_identical(Neoantimon2:::recognition_candidate_view(x, c$candidate_id)$summary$recognition_evidence_status, "no_evidence_supplied")
  expect_identical(c$final_grade, "NOT_GRADED")
})
test_that("unconfigured execution and licenses are passive capabilities", {
  for (provider in c("prime", "pmtnet", "ergo_ii", "nettcr_2_0")) {
    cap <- recognition_test_cap(provider)
    expect_true(cap$parser_available); expect_false(cap$execution_available)
    expect_false(cap$local_dependencies_detected)
    expect_identical(cap$license_status, "not_verified")
  }
  p <- recognition_test_cap("prime")$config
  command <- Neoantimon2:::recognition_prime_command(p, "in.txt", "out.txt", "A0201")
  expect_identical(command$status, "unavailable")
  file <- tempfile(); writeLines("synthetic executable placeholder", file)
  p$executable <- file; p$model_files <- list(file); p$dependencies <- list(file); p$license_status <- "user_confirmed"
  command <- Neoantimon2:::recognition_prime_command(p, "in.txt", "out.txt", "A0201")
  expect_false(command$executed)
  expect_identical(command$status, "configured_command_only")
})
test_that("frozen somatic IDs grades ranks eligibility and priority origins are exact invariants", {
  p <- recognition_e2e_pair()
  before <- neo_results(p$before); after <- neo_results(p$after)
  expect_gt(nrow(before), 0L)
  expect_identical(recognition_invariance_fields(before), recognition_invariance_fields(after))
  expect_true(all(c("A", "B", "C", "E") %in% before$final_grade))
  expect_true(all(after$policy_version == "somatic-priority-v1"))
  expect_true(all(after$priority_origin_selection_policy == "best-supported-eligible-origin-v1"))
  expect_true(any(after$tcr_specific_prediction_count > 0))
  expect_true(any(after$peptide_intrinsic_evidence_count > 0))
  before_predictions <- Neoantimon2:::somatic_prediction_runs(p$before)
  after_predictions <- Neoantimon2:::somatic_prediction_runs(p$after)
  project_metrics <- function(runs) {
    rows <- unlist(lapply(runs, `[[`, "evidence"), recursive = FALSE)
    x <- lapply(rows, function(r) r[c("candidate_id", "task", "metric", "raw_score", "percentile_rank", "predictor", "predictor_version")])
    Neoantimon2:::sorted_records(x)
  }
  expect_identical(project_metrics(before_predictions), project_metrics(after_predictions))
})
test_that("TCR prediction E2E reaches compact results and separate explanation families", {
  p <- recognition_e2e_pair()$after
  rows <- neo_results(p)
  expect_identical(nrow(rows), length(Neoantimon2:::recognition_read(p)$candidate_ids))
  expect_true(all(rows$tcr_specific_prediction_count == 1L))
  expect_true(all(rows$recognition_evidence_status == "evidence_present"))
  explain <- neo_explain(p, candidate_id = rows$candidate_id[[1L]])$candidates[[1L]]
  expect_length(explain$recognition$tcr_specific, 1L)
  expect_length(explain$recognition$receptor_feature_profiles, 3L)
  expect_length(explain$recognition$repertoire_observations, 3L)
  predicted_profile <- explain$recognition$tcr_specific[[1L]]$receptor_feature_profile_id
  paired_observations <- Filter(function(o) o$receptor_feature_profile_id == predicted_profile,
    explain$recognition$repertoire_observations)
  expect_length(paired_observations, 1L)
  observed <- paired_observations[[1L]]
  expect_identical(observed$status, "receptor_observed")
  expect_identical(observed$observation$cell_id, "observed")
  expect_length(observed$observation$support, 2L)
  expect_true(all(vapply(observed$observation$support, function(s) is.null(s$duplicate_count), TRUE)))
  expect_true(explain$recognition$tcr_specific[[1L]]$input_feature_completeness == "paired_cdr3_profile")
  text <- neo_explain(p, candidate_id = rows$candidate_id[[1L]], format = "text")
  for (heading in c("PRESENTATION", "PEPTIDE-INTRINSIC RECOGNITION EVIDENCE", "TCR-SPECIFIC PREDICTION",
    "EXPERIMENTAL REACTIVITY", "CURATED SPECIFICITY MATCHES", "REPERTOIRE OBSERVATION")) expect_match(text, heading, fixed = TRUE)
  expect_false(grepl("candidate is immunogenic", text, fixed = TRUE))
  expect_identical(neo_run(p)$run_id, p$run_id)
})
test_that("experimental reactivity E2E retains units and documented threshold", {
  p <- recognition_e2e_pair()
  row <- neo_results(p$after, candidate_id = p$measured_candidate)
  expect_identical(row$experimental_reactivity_count, 1L)
  e <- neo_explain(p$after, candidate_id = p$measured_candidate)$candidates[[1L]]$recognition$experimental_reactivity[[1L]]
  expect_identical(e$raw_score, 120)
  expect_identical(e$metadata$unit, "SFU")
  expect_identical(e$metadata$threshold, 10L)
  expect_identical(e$metadata$qualitative_annotation, "positive")
  expect_identical(e$metadata$target_level, "cell")
  expect_true(startsWith(e$metadata$receptor_assignment, "unique_cell_profile_link"))
})
test_that("recognition sidecar verifies artifacts and output hashes", {
  p <- recognition_e2e_pair()$after
  x <- Neoantimon2:::recognition_read(p); expect_silent(Neoantimon2:::recognition_validate(x))
  snapshot <- Neoantimon2:::read_runtime_snapshot(p)
  fp <- Neoantimon2:::recognition_fingerprints(p, Neoantimon2:::recognition_candidates(p, snapshot))
  root <- Neoantimon2:::recognition_root(p, fp); file <- file.path(root, "evidence.json")
  saved <- readBin(file, "raw", n = file.info(file)$size)
  on.exit(writeBin(saved, file), add = TRUE)
  writeLines("[]", file)
  expect_error(neo_results(p), "integrity")
  writeBin(saved, file)
  expect_identical(neo_results(p)$candidate_id, neo_results(recognition_e2e_pair()$before)$candidate_id)
})
test_that("fusion splice long-read and noncanonical source grades stay NOT_GRADED", {
  base <- long_read_fixture_project(list(long_read_fixture_input("novel-exon.gtf"),
    long_read_fixture_input("settings-denovo.json", "long_read_settings", "long-read-settings"),
    list(input_id = "fusion", input_role = "star_fusion", sample_id = "sample-1",
      path = normalizePath(test_path("..", "fixtures", "source", "star-fusion.tsv")), version = "synthetic-v1"),
    list(input_id = "splice", input_role = "splice_structure", sample_id = "sample-1",
      path = normalizePath(test_path("..", "fixtures", "source", "splice-structure.tsv")), version = "synthetic-v1")),
    hla = list("patient-1" = "HLA-A*02:01"))
  cfg <- base$config; cfg$peptides <- Neoantimon2:::peptide_policy(list(class_i_lengths = as.list(6:11)))
  base <- neo_project(samples = base$samples, config = cfg, hla = lapply(base$hla, unlist),
    inputs = recognition_clone_inputs(base))
  before <- neo_run(base); prior <- neo_results(before)
  after <- neo_run(recognition_clone_with(before, list(recognition_fixture_input("unique-paired.tsv"),
    recognition_fixture_input("unique-cell.json", "tcr_airr"), recognition_temp_settings())))
  post <- neo_results(after)
  expect_gt(nrow(prior), 0L)
  expect_identical(recognition_invariance_fields(prior), recognition_invariance_fields(post))
  expect_true(all(post$final_grade == "NOT_GRADED"))
  expect_true(all(post$tcr_specific_prediction_count == 1L))
  lr <- Neoantimon2:::long_read_read(after)
  expect_true(any(vapply(lr$hypotheses, function(h) h$orf$evidence_type == "de_novo_aug_orf", TRUE)))
  expect_true(any(grepl("long_read_isoform", post$source_type)))
  expect_true(any(grepl("fusion", post$source_type)))
  expect_true(any(grepl("novel_junction", post$source_type)))
})
test_that("public exports and recognition contract remain independent", {
  expect_setequal(getNamespaceExports("Neoantimon2"), current_public_api)
  expect_identical(Neoantimon2:::recognition_version(), "tcr-recognition-1.0.0")
  expect_identical(Neoantimon2:::recognition_id_versions()$profile, "tcr-feature-profile-id-v1")
  expect_identical(Neoantimon2:::recognition_id_versions()$query, "recognition-query-id-v1")
  expect_identical(Neoantimon2:::somatic_priority_version(), "somatic-priority-v1")
})
test_that("local provider input rendering uses only declared features and stable order", {
  profile <- recognition_paired_profile()
  for (provider in c("pmtnet", "ergo_ii", "nettcr_2_0")) {
    cap <- recognition_test_cap(provider)
    queries <- lapply(list(recognition_test_candidate(), recognition_test_candidate("AAAAAAAAA")),
      function(c) Neoantimon2:::recognition_query(c, profile, cap))
    a <- Neoantimon2:::recognition_provider_input(queries, cap)
    b <- Neoantimon2:::recognition_provider_input(rev(queries), cap)
    expect_identical(a, b)
    expect_length(a$rows, 2L)
    if (provider == "pmtnet") {
      expect_identical(names(a$rows[[1L]]), c("CDR3", "Antigen", "HLA"))
      expect_identical(a$rows[[1L]]$HLA, "HLA-A02:01")
    } else expect_null(a$rows[[1L]][[if (provider == "ergo_ii") "TRA" else "CDR3a"]])
  }
  cap <- recognition_test_cap("ergo_ii", "paired", use_mhc = FALSE)
  q <- Neoantimon2:::recognition_query(recognition_test_candidate(), profile, cap)
  rendered <- Neoantimon2:::recognition_provider_input(list(q), cap)$rows[[1L]]
  expect_identical(rendered$TRA, profile$features$alpha_cdr3_aa)
  expect_null(rendered$MHC)
})
test_that("ERGO-II explicit phenotype is consumed without class-based imputation", {
  r <- recognition_fixture_repertoire(TRUE)
  p <- Filter(function(p) identical(p$features$t_cell_type, "CD8") && p$completeness == "paired_cdr3_profile", r$profiles)[[1L]]
  cap <- recognition_test_cap("ergo_ii", "custom", features = list("beta_cdr3_aa", "t_cell_type"))
  q <- Neoantimon2:::recognition_query(recognition_test_candidate(), p, cap)
  expect_identical(q$features_used$t_cell_type, "CD8")
  expect_identical(q$input_status, "supported")
  expect_error(recognition_test_cap("ergo_ii", "custom",
    features = list("beta_cdr3_aa", "beta_variable_aa")), "Unknown receptor feature")
  unknown <- recognition_paired_profile(r)
  expect_identical(Neoantimon2:::recognition_query(recognition_test_candidate(), unknown, cap)$input_status, "incomplete_features")
})
test_that("missing Cell pairing declaration never implies observed pairing", {
  airr <- jsonlite::read_json(recognition_fixture_path("unique-cell.json"))
  airr$Cell[[1L]]$virtual_pairing <- NULL
  path <- tempfile(fileext = ".json"); jsonlite::write_json(airr, path, auto_unbox = TRUE)
  input <- recognition_fixture_input("unique-cell.json", "tcr_airr"); input$path <- path
  r <- Neoantimon2:::tcr_import(recognition_fixture_project(list(recognition_fixture_input("unique-paired.tsv"), input)))
  expect_length(r$pairing_hypotheses, 0L)
  expect_identical(r$pairing_statuses[[1L]]$status, "unpaired_chain")
})
test_that("AIRR inferred and annotated methods do not become patient assays", {
  airr <- jsonlite::read_json(recognition_fixture_path("reactivity.json"))
  airr$Reactivity <- airr$Reactivity[1:2]
  airr$Reactivity[[1L]]$reactivity_method <- "inferred"
  airr$Reactivity[[2L]]$reactivity_method <- "annotated"
  path <- tempfile(fileext = ".json"); jsonlite::write_json(airr, path, auto_unbox = TRUE)
  input <- recognition_fixture_input("reactivity.json", "tcr_reactivity"); input$path <- path
  r <- Neoantimon2:::tcr_import(recognition_fixture_project(list(recognition_fixture_input("unique-paired.tsv"),
    recognition_fixture_input("unique-cell.json", "tcr_airr"), input)))
  e <- Neoantimon2:::recognition_experiments(list(recognition_test_candidate()), r)
  expect_length(e$measurements, 2L)
  expect_length(e$evidence, 1L)
  expect_identical(e$evidence[[1L]]$scope, "external_annotation")
  expect_identical(e$evidence[[1L]]$metadata$assay_type, "annotated")
})
test_that("unknown database HLA does not become a pMHC-specific restriction", {
  input <- recognition_fixture_project(list(recognition_fixture_input("vdjdb.tsv", "tcr_vdjdb")))$inputs[[1L]]
  a <- Neoantimon2:::recognition_vdjdb(input)
  a <- Filter(function(a) a$presenter$status == "unknown", a)
  c <- recognition_test_candidate(hla = "HLA-A*01:01")
  result <- Neoantimon2:::recognition_curated_matches(list(c), recognition_fixture_repertoire(TRUE), a)
  expect_gt(length(result$evidence), 0L)
  expect_true(all(vapply(result$evidence, function(e) e$mhc_conditioning == "unknown", TRUE)))
})
test_that("AIRR direct receptor reference retains declared level and official hash", {
  p <- recognition_fixture_project(list(recognition_fixture_input("paired.tsv"), recognition_fixture_input("airr.json", "tcr_airr")))
  r <- Neoantimon2:::tcr_import(p)
  measurement <- list(reactivity_id = "domain-assay", receptor_id = "full", peptide_sequence_aa = "GILGFVFTL",
    mhc_allele_1 = "HLA-A*02:01", reactivity_method = "affinity_assay", reactivity_readout = "dissociation",
    reactivity_value = 3, reactivity_unit = "micromolar")
  r$reactivities <- list(list(sample_id = "sample-1", raw = measurement, artifact = r$receptors[[1L]]$artifact))
  e <- Neoantimon2:::recognition_experiments(list(recognition_test_candidate()), r)$evidence[[1L]]
  expect_identical(e$metadata$target_level, "receptor_declared")
  expect_identical(e$input_feature_completeness, "declared_receptor_profile_link")
  expect_identical(e$receptor_feature_profile_id, Filter(function(p) p$completeness == "paired_variable_domain_receptor", r$profiles)[[1L]]$receptor_feature_profile_id)
})
test_that("AIRR unknown or conflicting HLA never receives a known restriction", {
  expect_identical(Neoantimon2:::recognition_reactivity_presenter(list(mhc_allele_1 = "HLA-A*02"))$status, "unresolved")
  expect_identical(Neoantimon2:::recognition_reactivity_presenter(list())$status, "unknown")
  expect_identical(Neoantimon2:::recognition_reactivity_presenter(list(mhc_class = "II", mhc_allele_1 = "HLA-DQA1*05:01",
    mhc_allele_2 = "HLA-DQB1*02:01"))$mhc_key, "HLA-DQA1*05:01~HLA-DQB1*02:01")
})
test_that("AIRR Reactivity named TSV parser retains raw extra columns", {
  path <- tempfile(fileext = ".tsv")
  writeLines(c("reactivity_unit\treactivity_id\treactivity_method\treactivity_readout\treactivity_value\textra",
    "AU\tr1\tfluorescence\tintensity\t0\tsynthetic"), path)
  input <- recognition_fixture_input("reactivity.json", "tcr_reactivity"); input$path <- path
  r <- Neoantimon2:::tcr_import(recognition_fixture_project(list(input)))
  expect_length(r$reactivities, 1L)
  expect_identical(r$reactivities[[1L]]$raw$extra, "synthetic")
  expect_null(r$reactivities[[1L]]$airr_schema_version)
})
test_that("declared version families retain the frozen adapter scope", {
  c <- recognition_test_candidate(); p <- recognition_paired_profile()
  for (version in c("2.1", "2.1.1", "unknown"))
    expect_identical(Neoantimon2:::recognition_input_status(c, NULL,
      recognition_test_cap("prime", provider_version = version)), "supported")
  expect_identical(Neoantimon2:::recognition_input_status(c, NULL,
    recognition_test_cap("prime", provider_version = "3.0")), "unsupported_provider_version")
  expect_identical(Neoantimon2:::recognition_input_status(c, p,
    recognition_test_cap("nettcr_2_0", provider_version = "2.2")), "unsupported_provider_version")
})
test_that("intrinsic fixture scope requires no receptor feature and is not replicated", {
  cap <- recognition_test_cap(scope = "peptide_intrinsic")
  prepared <- Neoantimon2:::recognition_prepare(list(recognition_test_candidate()), recognition_fixture_repertoire(TRUE),
    cap, recognition_test_settings())
  expect_length(prepared$queries, 1L)
  expect_null(prepared$queries[[1L]]$receptor_feature_profile_id)
  expect_identical(prepared$queries[[1L]]$features_used, list())
  expect_length(Neoantimon2:::recognition_predict(prepared, cap)$evidence, 1L)
})

test_that("AIRR functional evidence gates every TCR-specific provider without deleting observations", {
  for (case in recognition_functionality_cases()) {
    r <- recognition_functionality_repertoire(case$flags)
    beta <- Filter(function(c) c$locus == "TRB", r$chains)[[1L]]
    expect_identical(beta$functional_status, case$expected, info = case$name)
    for (key in c("productive", "stop_codon", "vj_in_frame"))
      expect_identical(beta[[key]], case$flags[[key]], info = paste(case$name, key))
    expect_length(r$chains, 2L); expect_length(r$pairing_hypotheses, 1L)
    expect_identical(r$pairing_hypotheses[[1L]]$functional_status, case$expected)
    p <- recognition_paired_profile(r)
    for (provider in c("fixture", "pmtnet", "ergo_ii", "nettcr_2_0")) {
      cap <- recognition_test_cap(provider)
      prepared <- Neoantimon2:::recognition_prepare(list(recognition_test_candidate()), r, cap,
        recognition_test_settings(profile_ids = list(p$receptor_feature_profile_id)))
      q <- prepared$queries[[1L]]
      expect_identical(q$input_status, case$input_status, info = paste(case$name, provider))
      expect_identical(q$receptor_functionality$status, case$expected)
      output <- if (provider == "fixture") NULL else recognition_fixture_project(list(recognition_fixture_input(
        switch(provider, pmtnet = "pmtnet.csv", ergo_ii = "ergo.csv", nettcr_2_0 = "nettcr.csv"), "recognition_output")))$inputs[[1L]]
      result <- Neoantimon2:::recognition_predict(prepared, cap, output)
      expect_length(result$evidence, if (case$input_status == "supported") 1L else 0L)
      expect_identical(result$statuses[[1L]]$status, if (case$input_status == "supported") "available" else case$input_status)
      rendered <- Neoantimon2:::recognition_provider_input(prepared$queries, cap)
      expect_length(rendered$rows, if (case$input_status == "supported") 1L else 0L)
    }
  }
})

test_that("functionality changes observation evidence but never sequence profile or query identity", {
  original <- recognition_functionality_repertoire(list(productive = TRUE, stop_codon = FALSE, vj_in_frame = TRUE))
  base <- recognition_paired_profile(original)
  cap <- recognition_test_cap("fixture", "paired")
  prepare <- function(r) Neoantimon2:::recognition_prepare(list(recognition_test_candidate()), r, cap,
    recognition_test_settings(profile_ids = list(base$receptor_feature_profile_id)))$queries[[1L]]
  q <- prepare(original)
  for (flags in list(list(productive = FALSE), list(stop_codon = TRUE), list(vj_in_frame = FALSE), list())) {
    changed <- recognition_functionality_repertoire(flags)
    expect_identical(original$profiles, changed$profiles)
    expect_false(identical(original$chains, changed$chains))
    expect_false(identical(original$observations, changed$observations))
    expect_identical(q$recognition_query_id, prepare(changed)$recognition_query_id)
    expect_identical(q$fingerprint, prepare(changed)$fingerprint)
    expect_identical(q$candidate_id, prepare(changed)$candidate_id)
  }
})

test_that("one nonfunctional required chain blocks paired hypotheses independently", {
  for (locus in c("TRA", "TRB")) {
    r <- recognition_functionality_repertoire(list(productive = FALSE), locus)
    expect_length(r$pairing_hypotheses, 1L)
    expect_identical(r$pairing_hypotheses[[1L]]$pairing_status, "observed_pair")
    cap <- recognition_test_cap("fixture", "paired")
    result <- Neoantimon2:::recognition_predict(Neoantimon2:::recognition_prepare(
      list(recognition_test_candidate()), r, cap, recognition_test_settings()), cap)
    paired <- Filter(function(s) s$receptor_feature_profile_id == recognition_paired_profile(r)$receptor_feature_profile_id, result$statuses)
    expect_identical(paired[[1L]]$status, "nonfunctional_receptor")
    expect_length(result$evidence, 0L)
  }
  input <- recognition_fixture_input("paired.tsv")
  rows <- Neoantimon2:::source_table(input$path, c("sequence_id", "locus", "junction_aa"), "AIRR")
  alpha <- Filter(function(c) c$cell_id == "multiple" && c$locus == "TRA", rows)[[2L]]$sequence_id
  for (i in seq_along(rows)) if (rows[[i]]$sequence_id == alpha) rows[[i]]$productive <- FALSE
  x <- jsonlite::read_json(recognition_fixture_path("airr.json")); x$Rearrangement <- rows
  path <- tempfile(fileext = ".json"); jsonlite::write_json(x, path, auto_unbox = TRUE)
  input$path <- path; input$input_role <- "tcr_airr"
  r <- Neoantimon2:::tcr_import(recognition_fixture_project(list(input)))
  multi <- Filter(function(o) identical(o$cell_id, "multiple"), r$pairing_hypotheses)
  expect_length(multi, 2L)
  expect_setequal(vapply(multi, `[[`, "", "functional_status"), c("functional", "explicitly_nonfunctional"))
  cap <- recognition_test_cap("fixture", "paired")
  prepared <- Neoantimon2:::recognition_prepare(list(recognition_test_candidate()), r, cap,
    recognition_test_settings(profile_ids = as.list(vapply(multi, `[[`, "", "receptor_feature_profile_id"))))
  expect_setequal(vapply(prepared$queries, `[[`, "", "input_status"), c("supported", "nonfunctional_receptor"))
})

test_that("shared sequence profile retains distinct functional observations in its sample", {
  r <- recognition_functionality_repertoire(list(productive = FALSE))
  functional <- recognition_functionality_repertoire(list(productive = TRUE))
  r$observations <- c(r$observations, functional$observations)
  expect_identical(r$profiles, functional$profiles)
  p <- recognition_paired_profile(r)
  own <- Filter(function(o) o$receptor_feature_profile_id == p$receptor_feature_profile_id, r$observations)
  expect_length(own, 2L)
  expect_length(unique(vapply(own, `[[`, "", "repertoire_observation_id")), 2L)
  cap <- recognition_test_cap("fixture", "paired")
  settings <- recognition_test_settings(profile_ids = list(p$receptor_feature_profile_id))
  q <- Neoantimon2:::recognition_prepare(list(recognition_test_candidate()), r, cap, settings)$queries[[1L]]
  expect_identical(q$input_status, "supported")
  expect_setequal(unlist(q$receptor_functionality$observation_statuses), c("functional", "explicitly_nonfunctional"))
  for (i in seq_along(r$observations)) if (r$observations[[i]]$functional_status == "functional")
    r$observations[[i]]$sample_id <- "normal"
  q <- Neoantimon2:::recognition_prepare(list(recognition_test_candidate()), r, cap, settings)$queries[[1L]]
  expect_identical(q$input_status, "nonfunctional_receptor")
})

test_that("PRIME queries and evidence are identical for every receptor functional status", {
  cap <- recognition_test_cap("prime")
  output <- recognition_fixture_project(list(recognition_fixture_input("prime.tsv", "recognition_output")))$inputs[[1L]]
  results <- lapply(recognition_functionality_cases(), function(case) {
    prepared <- Neoantimon2:::recognition_prepare(list(recognition_test_candidate()),
      recognition_functionality_repertoire(case$flags), cap, recognition_test_settings())
    list(queries = prepared$queries, result = Neoantimon2:::recognition_predict(prepared, cap, output))
  })
  for (result in results) expect_identical(result, results[[1L]])
  expect_null(results[[1L]]$queries[[1L]]$receptor_functionality)
})

test_that("experimental and curated evidence survive nonfunctional and conflicting rearrangements", {
  annotations <- Neoantimon2:::recognition_vdjdb(recognition_fixture_project(list(
    recognition_fixture_input("vdjdb.tsv", "tcr_vdjdb")))$inputs[[1L]])
  for (flags in list(list(productive = FALSE), list(productive = TRUE, stop_codon = TRUE), list())) {
    project <- recognition_fixture_project(list(recognition_functionality_input(flags),
      recognition_fixture_input("reactivity.json", "tcr_reactivity")))
    r <- Neoantimon2:::tcr_import(project)
    experiments <- Neoantimon2:::recognition_experiments(list(recognition_test_candidate()), r)
    expect_length(experiments$measurements, 5L); expect_length(experiments$evidence, 4L)
    positive <- Filter(function(e) e$metadata$reactivity_id == "positive", experiments$evidence)[[1L]]
    expect_identical(positive$raw_score, 120)
    expect_identical(positive$metadata$qualitative_annotation, "positive")
    expect_identical(positive$metadata$receptor_functionality_context[[1L]]$functional_status,
      Neoantimon2:::tcr_functional_status(flags))
    curated <- Neoantimon2:::recognition_curated_matches(list(recognition_test_candidate()), r, annotations)
    expect_gt(length(curated$evidence), 0L)
    expect_true(all(vapply(curated$evidence, function(e) length(e$metadata$observations) > 0L, TRUE)))
    unknown <- Filter(function(e) e$mhc_conditioning == "unknown", curated$evidence)
    expect_gt(length(unknown), 0L)
    expect_match(unknown[[1L]]$metadata$match_status, "restriction_unknown")
  }
})

test_that("AIRR MHC class normalization preserves complete canonical presenters and conflicts", {
  cases <- jsonlite::read_json(recognition_fixture_path("mhc-classes.json"))
  for (case in cases) {
    p <- Neoantimon2:::recognition_reactivity_presenter(case$row)
    expect_identical(p$status, case$status, info = case$name)
    expect_identical(p$mhc_key, case$key, info = case$name)
    if (!is.null(case$reason)) expect_identical(p$reason, case$reason, info = case$name)
    expect_identical(p$raw_alleles$mhc_allele_2, case$row$mhc_allele_2)
  }
  for (raw in c("I", "MHC-I", "MHC I", "mhc-i")) expect_identical(Neoantimon2:::recognition_mhc_class(raw), "I")
  for (raw in c("II", "MHC-II", "MHC II", "mhc ii")) expect_identical(Neoantimon2:::recognition_mhc_class(raw), "II")
  expect_null(Neoantimon2:::recognition_mhc_class("CD4"))
  expect_null(Neoantimon2:::recognition_mhc_class(NULL))
  for (pair in list(c("HLA-DRA*01:01", "HLA-DRB1*04:01"), c("HLA-DPA1*01:03", "HLA-DPB1*04:01"))) {
    p <- Neoantimon2:::recognition_reactivity_presenter(list(mhc_class = "MHC-II", mhc_allele_1 = pair[[2L]], mhc_allele_2 = pair[[1L]]))
    expect_identical(p$mhc_key, paste(pair, collapse = "~"))
  }
  expect_identical(Neoantimon2:::recognition_reactivity_presenter(list(mhc_class = "unspecified", mhc_allele_1 = "HLA-A*02:01", t_cell_type = "CD8"))$status, "unresolved")
})

test_that("partial conflicting and unknown AIRR restrictions preserve raw measurements safely", {
  cases <- jsonlite::read_json(recognition_fixture_path("mhc-classes.json"))
  for (case in cases) {
    row <- c(list(reactivity_id = case$name, reactivity_method = "ELISpot", reactivity_readout = "spots",
      reactivity_value = 12, reactivity_unit = "SFU", peptide_sequence_aa = "GILGFVFTL"), case$row)
    r <- list(reactivities = list(list(sample_id = "sample-1", raw = row, artifact = list(version = "synthetic"))), observations = list())
    candidates <- list(recognition_test_candidate(), recognition_test_candidate(hla = "HLA-DQA1*05:01~HLA-DQB1*02:01"))
    result <- Neoantimon2:::recognition_experiments(candidates, r)
    expect_length(result$measurements, 1L)
    expect_identical(result$measurements[[1L]]$raw_fields, row)
    expect_length(result$evidence, if (case$status == "known") 1L else 0L)
  }
  r <- list(reactivities = list(list(sample_id = "sample-1", raw = list(reactivity_id = "unknown",
    reactivity_method = "ELISpot", reactivity_readout = "spots", reactivity_value = 12, reactivity_unit = "SFU",
    peptide_sequence_aa = "GILGFVFTL", mhc_class = "unspecified"), artifact = list())), observations = list())
  result <- Neoantimon2:::recognition_experiments(list(recognition_test_candidate()), r)
  expect_length(result$evidence, 1L)
  expect_identical(result$evidence[[1L]]$mhc_conditioning, "unknown")
})
