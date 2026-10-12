test_that("hashed grouping preserves first records, complete graphs and validation", {
  old <- external_scaling_legacy()
  f <- external_scaling_workload(24L)
  f$records <- rev(f$records)
  f$origins <- c(rev(f$origins), f$origins[1L])
  f$events <- list(list(event_id = "event-z"), list(event_id = "event-a"))
  f$annotations <- list(core_peptides = list(list(provenance_occurrences = list(list(source = "one"), list(source = "two")))))
  compare <- function(f) {
    args <- f[c("records", "origins", "events", "annotations")]; names(args)[[1L]] <- "candidates"
    expect_identical(do.call(Neoantimon2:::external_merge_graph, args), do.call(old$legacy_external_merge_graph, args))
  }
  compare(f)
  for (field in c("sample_id", "peptide_sequence", "mhc_key", "hla_class")) {
    bad <- f; bad$records[[length(bad$records) + 1L]] <- bad$records[[1L]]
    bad$records[[length(bad$records)]][[field]] <- "conflict"
    expect_error(Neoantimon2:::external_merge_graph(bad$records, bad$origins), "ANNOTATION_CONFLICT")
    expect_error(old$legacy_external_merge_graph(bad$records, bad$origins), "ANNOTATION_CONFLICT")
  }
  bad <- f; bad$origins[[length(bad$origins)]]$sample_id <- "conflict"
  expect_error(Neoantimon2:::external_merge_graph(bad$records, bad$origins), "ANNOTATION_CONFLICT")
  expect_error(old$legacy_external_merge_graph(bad$records, bad$origins), "ANNOTATION_CONFLICT")
  bad <- f; bad$records[[1]]$origin_ids <- list("missing")
  expect_error(Neoantimon2:::external_merge_graph(bad$records, bad$origins), "MISSING_ORIGIN")
  expect_error(old$legacy_external_merge_graph(bad$records, bad$origins), "MISSING_ORIGIN")
})

test_that("indexed exact fan-out and all interpretation rows equal the legacy algorithms", {
  old <- external_scaling_legacy()
  g <- external_test_graph()
  pair <- "HLA-DQA1*05:01~HLA-DQB1*02:01"
  for (id in c("c4", "c5")) {
    c <- g$candidates[[1]]; c$candidate_id <- id; c$mhc_key <- pair; c$hla_class <- "II"
    c$sample_id <- if (id == "c4") "s1" else "s2"; c$origin_ids <- list(paste0(id, "a"), paste0(id, "b"))
    g$candidates[[length(g$candidates) + 1L]] <- c
  }
  missing <- g$candidates[[1]]; missing$candidate_id <- "c6"; missing$peptide_sequence <- "UUU"
  g$candidates[[6L]] <- missing
  g$candidates[[7L]] <- missing; g$candidates[[7L]]$candidate_id <- "c7"; g$candidates[[7L]]$peptide_sequence <- "LMNPQRSTV"
  x <- rbind(external_test_row(), external_test_row(mhc = "HLA-B*07:02", score = "0.1", source_row = "X:2"),
    external_test_row(mhc = pair, score = "", status = "failed", reason = "failure", source_row = "X:3"),
    external_test_row(scope = "peptide", mhc_key = "", metric = "intrinsic", source_row = "X:4"),
    external_test_row(scope = "peptide", mhc_key = "", task = "processing.class1", source_row = "X:5"),
    external_test_row(mhc = pair, task = "mhc1_binding", source_row = "X:6"),
    external_test_row(task = "mhc2_binding", source_row = "X:7"),
    external_test_row("UUU", score = "", status = "unsupported", reason = "residue", source_row = "X:8"),
    external_test_row("LMNPQRSTV", score = "", status = "not_evaluated", reason = "unrun", source_row = "X:9"),
    external_test_row("YYYY", source_row = "X:10"))
  x <- Neoantimon2:::external_validate_predictions(x)
  matched <- Neoantimon2:::external_match_predictions(g, x)
  expect_identical(matched, old$legacy_external_match_predictions(g, x))
  expect_equal(matched$report$unmatched_input_keys, 3L)
  for (direction in c("higher_is_better", "lower_is_better")) {
    y <- x; y$score_direction <- direction
    s <- Neoantimon2:::external_policy(list(mode = "external_rank", provider = "new-model",
      run_id = "external-run-1", task = "external:example", metric = "score", score_direction = direction), y)
    m <- Neoantimon2:::external_match_predictions(g, y)
    legacy_m <- old$legacy_external_match_predictions(g, y)
    rows <- Neoantimon2:::external_derive(NULL, g, m, s)
    legacy_rows <- old$legacy_external_derive(NULL, g, legacy_m, s)
    expect_identical(rows, legacy_rows)
    shuffled <- Neoantimon2:::external_validate_predictions(y[rev(seq_len(nrow(y))), ])
    expect_identical(Neoantimon2:::external_match_predictions(g, shuffled), m)
    tables <- list(rows, m$evidence, m$unmatched)
    legacy_tables <- list(legacy_rows, legacy_m$evidence, legacy_m$unmatched)
    for (i in seq_along(tables)) {
      table <- tables[[i]]
      fields <- Neoantimon2:::external_fields(table, "candidate_id")
      a <- tempfile(); b <- tempfile()
      Neoantimon2:::external_write_tsv(table, fields, a)
      Neoantimon2:::external_write_tsv(legacy_tables[[i]], fields, b)
      expect_identical(readBin(a, "raw", file.info(a)$size), readBin(b, "raw", file.info(b)$size))
      expect_identical(jsonlite::toJSON(table, auto_unbox = TRUE, null = "null", digits = NA),
        jsonlite::toJSON(legacy_tables[[i]], auto_unbox = TRUE, null = "null", digits = NA))
    }
  }
  s <- Neoantimon2:::external_policy(list(mode = "evidence_only"), x)
  expect_identical(Neoantimon2:::external_derive(NULL, g, matched, s), old$legacy_external_derive(NULL, g, matched, s))
  expect_error(Neoantimon2:::external_validate_predictions(rbind(x, x[1L, ])), "DUPLICATE_PREDICTION")
})

test_that("zero and singleton indexed workloads preserve exact output identities", {
  old <- external_scaling_legacy()
  for (n in c(0L, 1L, 37L, 1000L)) {
    f <- external_scaling_workload(n)
    g <- Neoantimon2:::external_merge_graph(f$records, f$origins)
    expect_identical(g, old$legacy_external_merge_graph(f$records, f$origins))
    m <- Neoantimon2:::external_match_predictions(g, f$predictions)
    expect_identical(m, old$legacy_external_match_predictions(g, f$predictions))
    expect_identical(Neoantimon2:::external_derive(NULL, g, m, f$selection),
      old$legacy_external_derive(NULL, g, m, f$selection))
  }
})

test_that("indexed explanation keeps origin order and the complete public schema", {
  f <- external_scaling_workload(24L)
  graph <- Neoantimon2:::external_merge_graph(f$records, f$origins,
    list(list(event_id = "global-event")), list(source = list(observations = list(list(source_id = "global-observation")))))
  m <- Neoantimon2:::external_match_predictions(graph, f$predictions)
  rows <- Neoantimon2:::external_derive(NULL, graph, m, f$selection)
  payload <- list(manifest = list(analysis_id = "test-id", fingerprint = list(bytes = "original")),
    payload = list(rows = rows, graph = graph, policy = f$selection$policy, import_report = m$report, conversion = NULL))
  old <- external_scaling_legacy()
  testthat::with_mocked_bindings({
    for (ids in list(c("c00000003", "c00000001"), "absent", vapply(rows, `[[`, "", "candidate_id"))) {
      expect_identical(neo_external_explain("test", ids), old$legacy_neo_external_explain("test", ids))
    }
  }, external_analysis_read = function(...) payload, .package = "Neoantimon2")
})

test_that("automatic bundle bytes and wide metric validation equal the baseline", {
  old <- external_scaling_legacy()
  p <- peptide_test_run(sec = TRUE)
  new <- neo_export_bundle(p)
  baseline <- old$legacy_neo_export_bundle(p, tempfile("legacy-export-"))
  expect_identical(new$manifest, baseline$manifest)
  expect_identical(new$graph, baseline$graph)
  expect_identical(vapply(new$manifest$files, `[[`, "", "sha256"),
    vapply(baseline$manifest$files, `[[`, "", "sha256"))
  x <- external_scaling_workload(37L)$predictions
  x$model_version <- as.character(seq_len(nrow(x)))
  expect_identical(Neoantimon2:::external_validate_predictions(x), old$legacy_external_validate_predictions(x))
  row <- which(x$scope == "peptide_hla")[[1L]]
  bad <- rbind(x, x[row, ]); bad$mhc_key[nrow(bad)] <- if (x$mhc_key[[row]] == "HLA-B*07:02") "HLA-A*02:01" else "HLA-B*07:02"
  bad$score_direction[nrow(bad)] <- "lower_is_better"
  expect_error(Neoantimon2:::external_validate_predictions(bad), "AMBIGUOUS_METRIC")
  expect_error(old$legacy_external_validate_predictions(bad), "AMBIGUOUS_METRIC")
})
