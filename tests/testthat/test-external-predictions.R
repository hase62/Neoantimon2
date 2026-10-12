test_that("exact query matching retains every origin and distinct samples and HLA", {
  x <- rbind(external_test_row(), external_test_row(mhc = "HLA-B*07:02", score = "0.2", source_row = "X:2"))
  x <- Neoantimon2:::external_validate_predictions(x)
  m <- Neoantimon2:::external_match_predictions(external_test_graph(), x)
  expect_length(m$evidence, 3)
  expect_identical(m$evidence[[1]]$origin_ids, list("o1", "o2"))
  expect_identical(m$evidence[[1]]$evidence_id, m$evidence[[2]]$evidence_id)
  expect_equal(m$evidence[[3]]$numeric_score, 0.2)
  expect_equal(m$report$fan_out_copies, 1)
  expect_equal(m$report$origin_associations, 4)
  expect_identical(Neoantimon2:::external_match_predictions(external_test_graph(),
    Neoantimon2:::external_validate_predictions(x[2:1, ])), m)
  peptide <- Neoantimon2:::external_validate_predictions(external_test_row(scope = "peptide", mhc_key = ""))
  expect_length(Neoantimon2:::external_match_predictions(external_test_graph(), peptide)$evidence, 3)
  unknown <- Neoantimon2:::external_match_predictions(external_test_graph(),
    Neoantimon2:::external_validate_predictions(external_test_row("UUU")))
  expect_equal(unknown$report$unmatched_input_keys, 1)
  expect_equal(unknown$report$missing_candidate_predictions, 3)
  expect_error(Neoantimon2:::external_validate_predictions(rbind(x, x[1, ])), "DUPLICATE")
})

test_that("class II requires the exact full presenter and task applicability", {
  key <- "HLA-DQA1*05:01~HLA-DQB1*02:01"
  x <- Neoantimon2:::external_validate_predictions(external_test_row(mhc = key, task = "mhc2_binding"))
  g <- external_test_graph(); g$candidates[[1]]$mhc_key <- key; g$candidates[[1]]$hla_class <- "II"
  expect_length(Neoantimon2:::external_match_predictions(g, x)$evidence, 1)
  for (bad in c("HLA-DQB1*02:01", "HLA-DQB1*02:01~HLA-DQA1*05:01", "hla-a*02:01", "A0201", ""))
    expect_error(Neoantimon2:::external_validate_predictions(external_test_row(mhc = bad)), "INVALID_PRESENTER")
  expect_error(Neoantimon2:::external_validate_predictions(external_test_row(scope = "peptide", mhc_key = "", task = "mhc1_binding")), "UNSUPPORTED_SCOPE")
})

test_that("strict values, unsupported states, scopes and declared percentile semantics", {
  validate <- Neoantimon2:::external_validate_predictions
  for (bad in c("", "NA", "NaN", "Inf", "-Inf", "abc"))
    expect_error(validate(external_test_row(score = bad)), "INVALID_SCORE")
  for (bad in c("acdef", "AX", "AC*", " AC", ""))
    expect_error(validate(external_test_row(bad)))
  expect_silent(validate(external_test_row("U", status = "unsupported", reason = "unsupported_residue", score = "")))
  expect_error(validate(external_test_row(status = "failed", reason = "error")), "INVALID_SCORE")
  expect_error(validate(external_test_row(scope = "sample")), "UNSUPPORTED_SCOPE")
  expect_error(validate(external_test_row(task = "recognition.tcr_specific")), "UNSUPPORTED_SCOPE")
  expect_error(validate(external_test_row(task = "made_up")), "UNSUPPORTED_TASK")
  expect_error(validate(external_test_row(score_direction = "unknown")), "INVALID_DIRECTION")
  expect_error(validate(external_test_row(metric_kind = "probability")), "INVALID_METRIC_KIND")
  percent <- external_test_row(task = "mhc1_presentation", metric_kind = "percentile_rank",
    metric_semantics = "mhc_presentation_percentile", unit = "percent", score_direction = "lower_is_better")
  expect_silent(validate(percent))
  for (bad in c("-1", "101")) { percent$score <- bad; expect_error(validate(percent), "INVALID_PERCENTILE") }
  percent$score <- "0.5"; percent$score_direction <- "higher_is_better"
  expect_error(validate(percent), "INVALID_PERCENTILE")
  expect_error(validate(external_test_row(metric_semantics = "mhc_presentation_percentile")), "INVALID_SEMANTICS")
  expect_equal(nrow(Neoantimon2:::external_read_predictions(external_test_write(external_test_row()))), 1)
  expect_equal(nrow(validate(external_test_row()[FALSE, ])), 0)
})

test_that("TSV quoting preserves literal values and cannot sanitize malformed peptides", {
  x <- external_test_row(source_row = 'X:quoted "value"')
  expect_identical(Neoantimon2:::external_read_predictions(external_test_write(x))$source_row, x$source_row)
  expect_error(Neoantimon2:::external_read_predictions(external_test_write(external_test_row('AC"DE"'))), "INVALID_PEPTIDE")
  p <- external_test_write(external_test_row())
  lines <- readLines(p)
  lines[[2]] <- sub('"ACDEFGHIK"', 'AC"DEF"GHIK', lines[[2]], fixed = TRUE)
  writeLines(lines, p)
  expect_error(Neoantimon2:::external_read_predictions(p), "INVALID_TABLE")
  lines[[2]] <- paste0(lines[[2]], "\textra")
  writeLines(lines, p)
  expect_error(Neoantimon2:::external_read_predictions(p), "INVALID_TABLE")
  out <- tempfile(fileext = ".tsv")
  Neoantimon2:::external_write_tsv(list(list(origin_ids = list("o1", "o2"))), "origin_ids", out)
  cells <- read.delim(out, colClasses = "character")
  expect_identical(jsonlite::fromJSON(cells$origin_ids[[1]], simplifyVector = FALSE), list("o1", "o2"))
})
