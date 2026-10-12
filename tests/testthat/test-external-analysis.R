test_that("derived direct imports rank explicitly and preserve native files and results", {
  p <- peptide_test_run()
  b <- neo_export_bundle(p)
  q <- b$graph$candidates[[1]]
  x <- external_test_row(q$peptide_sequence, q$mhc_key)
  path <- external_test_write(x)
  native <- neo_results(p)
  explanation <- neo_explain(p, candidate_id = q$candidate_id)
  originals <- c(Neoantimon2:::manifest_file(p$path), Neoantimon2:::peptide_stage_path(p),
    Neoantimon2:::runtime_snapshot_path(p))
  hashes <- vapply(originals, Neoantimon2:::external_sha, "")
  policy <- list(mode = "external_rank", provider = "new-model", run_id = "external-run-1",
    task = "external:example", metric = "score", score_direction = "higher_is_better")
  a <- neo_analyze_external(p, path, policy)
  r <- neo_external_results(a)
  expect_equal(nrow(r), nrow(native))
  expect_true(any(r$external_status == "NOT_PREDICTED"))
  expect_equal(r$external_rank[r$candidate_id == q$candidate_id], 1)
  expect_true(all(is.na(r$external_rank[r$external_status == "NOT_PREDICTED"])))
  expect_identical(r$final_grade, vapply(b$graph$candidates, `[[`, "", "final_grade"))
  expect_identical(neo_external_results(a$analysis_id, neo_project(path = p$path)), r)
  e <- neo_external_explain(a, q$candidate_id)
  expect_length(e$candidates[[1]]$origins, length(q$origin_ids))
  expect_identical(e$candidates[[1]]$evidence[[1]]$source_row, "X:1")
  expect_identical(vapply(originals, Neoantimon2:::external_sha, ""), hashes)
  expect_identical(neo_results(p), native)
  expect_identical(neo_explain(p, candidate_id = q$candidate_id), explanation)
  expect_false(dir.exists(Neoantimon2:::predictor_artifact_root(p)))
  only <- neo_analyze_external(p, path, list(mode = "evidence_only"))
  expect_true(all(is.na(neo_external_results(only)$external_rank)))
  x$provider <- "second-model"
  second <- neo_analyze_external(p, external_test_write(x), list(mode = "evidence_only"))
  expect_false(identical(second$analysis_id, a$analysis_id))
  expect_identical(neo_external_results(a), r)
  writeLines("tampered", file.path(a$path, "results.tsv"))
  expect_error(neo_external_results(a), "INTEGRITY_FAILURE")
})

test_that("directional ranks reuse sample/class ordering with stable ties and missingness", {
  g <- external_test_graph()
  x <- rbind(external_test_row(), external_test_row(mhc = "HLA-B*07:02", score = "0.9", source_row = "X:2"))
  select <- function(x) Neoantimon2:::external_policy(list(mode = "external_rank", provider = "new-model",
    run_id = "external-run-1", task = "external:example", metric = "score", score_direction = "higher_is_better"), x)
  x <- Neoantimon2:::external_validate_predictions(x)
  m <- Neoantimon2:::external_match_predictions(g, x)
  r <- Neoantimon2:::external_derive(NULL, g, m, select(x))
  expect_equal(vapply(r, `[[`, 0L, "external_rank"), c(1L, 1L, 2L))
  x$score[x$mhc_key == "HLA-B*07:02"] <- "1.0"
  r <- Neoantimon2:::external_derive(NULL, g, Neoantimon2:::external_match_predictions(g, x), select(x))
  expect_equal(vapply(r, `[[`, 0L, "external_rank"), c(2L, 1L, 1L))
  expect_true(all(vapply(r, function(c) c$final_grade == "NOT_GRADED", TRUE)))
  expect_error(Neoantimon2:::external_policy(list(mode = "external_rank"), x), "provider")
  y <- x; y$model_version <- "2"
  expect_error(select(rbind(x, y)), "SELECTION_AMBIGUOUS")
  expect_error(Neoantimon2:::external_policy(list(mode = "somatic_grade", provider = "new-model",
    run_id = "external-run-1", task = "external:example", metric = "score", hla_class = "I",
    score_direction = "higher_is_better"), x), "GRADING_PREREQUISITES")
})

test_that("typed non-evaluated and missing predictions never fabricate negative values", {
  g <- external_test_graph()
  for (status in c("unsupported", "not_evaluated", "failed")) {
    x <- Neoantimon2:::external_validate_predictions(external_test_row(status = status, score = "", reason = "test_reason"))
    m <- Neoantimon2:::external_match_predictions(g, x)
    s <- Neoantimon2:::external_policy(list(mode = "external_rank", provider = "new-model", run_id = "external-run-1",
      task = "external:example", metric = "score", score_direction = "higher_is_better"), x)
    rows <- Neoantimon2:::external_derive(NULL, g, m, s)
    expect_identical(rows[[1]]$external_status, status)
    expect_identical(rows[[3]]$external_status, "NOT_PREDICTED")
    expect_true(all(vapply(rows, function(c) is.null(c$external_score) && is.null(c$external_rank), TRUE)))
    expect_true(all(vapply(rows, function(c) c$external_rank_status == "UNRANKED", TRUE)))
  }
  p <- peptide_test_run()
  a <- neo_analyze_external(p, external_test_write(external_test_row()[FALSE, ]), list(mode = "evidence_only"))
  expect_true(all(neo_external_results(a)$external_status == "NOT_PREDICTED"))
  expect_true(is.numeric(neo_external_results(a)$external_score))
  expect_true(is.numeric(neo_external_results(a)$external_rank))
  expect_equal(neo_external_explain(a, neo_external_results(a)$candidate_id[[1]])$import_report$input_rows, 0L)
})

test_that("declared external percent rank reuses somatic grading without native policy changes", {
  pre <- peptide_test_run()
  cfg <- pre$config
  event <- neo_events(pre)$event_id[[1]]
  dna <- tempfile(fileext = ".tsv")
  writeLines(c("event_id\ttumor_ref_reads\ttumor_alt_reads\ttumor_depth\tnormal_ref_reads\tnormal_alt_reads\tnormal_depth",
    paste(event, 20, 5, 25, 19, 1, 20, sep = "\t")), dna)
  cfg$somatic_evidence <- list(sources = list(list(kind = "dna_reads", sample_id = "sample-1", path = dna)))
  p <- neo_run(neo_project(samples = pre$samples, config = cfg,
    hla = lapply(pre$hla, unlist), inputs = lapply(pre$inputs, function(i)
      list(input_id = i$input_id, input_role = i$input_role, sample_id = i$sample_id, path = i$resolved_path))))
  q <- Filter(function(c) c$hla_class == "I", neo_export_bundle(p)$graph$candidates)[[1]]
  x <- external_test_row(q$peptide_sequence, q$mhc_key, task = "mhc1_presentation",
    metric = "declared_rank", metric_kind = "percentile_rank", metric_semantics = "mhc_presentation_percentile",
    score = "0.5", score_direction = "lower_is_better", unit = "percent")
  policy <- list(mode = "somatic_grade", provider = "new-model", run_id = "external-run-1",
    task = "mhc1_presentation", metric = "declared_rank", hla_class = "I", score_direction = "lower_is_better")
  originals <- c(Neoantimon2:::manifest_file(p$path), Neoantimon2:::peptide_stage_path(p),
    Neoantimon2:::runtime_snapshot_path(p), Neoantimon2:::somatic_priority_path(p),
    file.path(p$path, "results", "candidates.tsv"))
  hashes <- vapply(originals, Neoantimon2:::external_sha, "")
  native <- neo_results(p)
  native_explain <- neo_explain(p, candidate_id = q$candidate_id)
  a <- neo_analyze_external(p, external_test_write(x), policy)
  graph <- neo_export_bundle(p)$graph
  normalized <- Neoantimon2:::external_validate_predictions(x)
  matched <- Neoantimon2:::external_match_predictions(graph, normalized)
  selection <- Neoantimon2:::external_policy(policy, normalized)
  expect_identical(Neoantimon2:::external_analysis_read(a)$payload$rows,
    external_scaling_legacy()$legacy_external_derive(p, graph, matched, selection))
  row <- neo_external_explain(a, q$candidate_id)$candidates[[1]]$result
  expect_identical(row$external_final_grade, "B")
  expect_identical(row$final_grade, "NOT_GRADED")
  expect_identical(row$grading_status, "NOT_GRADED")
  expect_identical(row$external_priority$policy_version, "external-percentile-grade-v1")
  expect_identical(row$external_priority$primary_predictor, "new-model")
  expect_equal(row$external_priority$primary_rank, 0.5)
  expect_identical(row$external_priority$priority_origin_selection_policy, "best-supported-eligible-origin-v1")
  expect_identical(neo_results(p), native)
  expect_identical(neo_explain(p, candidate_id = q$candidate_id), native_explain)
  expect_identical(vapply(originals, Neoantimon2:::external_sha, ""), hashes)
  expect_error(Neoantimon2:::somatic_priority_policy(list(primary = list(I = list(provider = "new-model")))), "Unknown primary provider")
  x$metric_kind <- "raw_score"; x$metric_semantics <- "raw"
  expect_error(neo_analyze_external(p, external_test_write(x), policy), "GRADING_PREREQUISITES")
  x$metric_kind <- "percentile_rank"; x$metric_semantics <- "mhc_presentation_percentile"
  expect_error(neo_analyze_external(pre, external_test_write(x), policy), "GRADING_PREREQUISITES")
})
