test_that("source and long-read shared candidates export and fan out every origin", {
  model <- tempfile(fileext = ".gtf")
  writeLines(c('chr1\tx\texon\t1\t7\t.\t+\t.\ttranscript_id "shared";',
    'chr1\tx\texon\t21\t33\t.\t+\t.\ttranscript_id "shared";'), model)
  lr <- long_read_fixture_input(); lr$path <- model
  src <- list(input_id = "splice", input_role = "splice_structure", sample_id = "sample-1",
    path = normalizePath(test_path("..", "fixtures", "source", "splice-structure.tsv")), version = "v1")
  base <- long_read_fixture_project(list(lr, src))
  cfg <- base$config; cfg$peptides <- list(class_i_lengths = list(6L), class_ii_lengths = list(13L))
  p <- neo_run(neo_project(samples = base$samples, config = cfg, inputs = list(lr, src),
    hla = list("patient-1" = "HLA-A*02:01")))
  b <- neo_export_bundle(p)
  q <- Filter(function(c) length(c$origin_ids) > 1L, b$graph$candidates)[[1]]
  expect_gt(length(b$graph$annotations$source$observations), 0)
  expect_gt(length(b$graph$annotations$long_read$events), 0)
  expect_setequal(vapply(b$graph$candidates, `[[`, "", "candidate_id"), neo_results(p)$candidate_id)
  x <- external_test_row(q$peptide_sequence, q$mhc_key)
  a <- neo_analyze_external(p, external_test_write(x), list(mode = "evidence_only"))
  e <- neo_external_explain(a, q$candidate_id)$candidates[[1]]
  expect_length(e$origins, length(q$origin_ids))
  expect_gte(length(e$origins), 2)
  expect_length(e$evidence, 1)
  expect_setequal(unlist(e$evidence[[1]]$origin_ids), unlist(q$origin_ids))
})

test_that("public repository quickstart executes from preserved inputs to disposable runs", {
  repo <- normalizePath(test_path("..", ".."))
  script <- file.path(repo, "scripts", "development", "external_prediction_quickstart.R")
  skip_if_not(file.exists(script), "Repository-only quickstart excluded from installed package")
  old <- setwd(repo); on.exit(setwd(old), add = TRUE)
  data <- tempfile("quickstart-inputs-"); runs <- tempfile("quickstart-runs-")
  out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
    c(shQuote(script), "--data-root", shQuote(data), "--run-root", shQuote(runs)), stdout = TRUE, stderr = TRUE,
    timeout = 120))
  expect_null(attr(out, "status"), info = paste(out, collapse = "\n"))
  dirs <- list.dirs(runs, recursive = FALSE)
  expect_length(dirs, 1)
  summary <- jsonlite::read_json(file.path(dirs[[1]], "quickstart-summary.json"))
  expect_identical(summary$status, "PASS")
  expect_gte(summary$origin_count, 2L)
  expect_true(summary$native_hashes_unchanged)
})
