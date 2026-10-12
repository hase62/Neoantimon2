t004_config <- function() {
  root <- tempfile("t004-")
  neo_config(env = character(), paths = list(run_root = file.path(root, "runs"),
    data_root = file.path(root, "data")))
}

t004_fixture <- function() jsonlite::read_json(
  testthat::test_path("..", "fixtures", "truth", "simple_somatic.json"),
  simplifyVector = FALSE)

t004_project <- function(initialize = FALSE) {
  fixture <- t004_fixture()
  neo_project(samples = data.frame(sample_id = "sample-1", patient_id = "patient-1"),
    analysis_context = fixture$run_manifest$analysis_contexts,
    config = t004_config(), initialize = initialize, fixture_kind = "synthetic")
}

t004_adapter <- function(status = "NOT_RUN") {
  force(status)
  function(project) {
    x <- t004_fixture()
    x$run_manifest <- Neoantimon2:::read_manifest(project$path)$run_manifest
    for (table in c("events", "peptide_origins", "neoantigen_results"))
      x[[table]] <- lapply(x[[table]], function(row) {
        row$analysis_id <- project$analysis_id
        row
      })
    result <- x$neoantigen_results[[1L]]
    result$prediction_status <- status
    if (status %in% c("FAILED", "PARTIAL"))
      result$not_evaluated_codes <- list("PREDICTION_FAILED")
    x$neoantigen_results[[1L]] <- result
    x
  }
}

test_that("minimal project-run-accessor-resume path preserves ungraded candidates", {
  p <- Neoantimon2:::run_runtime(t004_project(), adapter = t004_adapter())
  results <- neo_results(p)
  expect_s3_class(results, "data.frame")
  expect_equal(nrow(results), 1L)
  expect_identical(results$prediction_status, "NOT_RUN")
  expect_identical(results$grading_status, "NOT_GRADED")
  expect_identical(results$final_grade, "NOT_GRADED")
  expect_equal(nrow(neo_events(p)), 1L)
  expect_equal(nrow(neo_origins(p)), 1L)
  expect_equal(nrow(neo_multisample(p)), 0L)

  manifest <- Neoantimon2:::read_manifest(p$path)
  expect_identical(manifest$status, "COMPLETED")
  expect_identical(manifest$checkpoints[[1L]]$status, "COMPLETED")
  resumed <- neo_project(path = p$path)
  expect_identical(neo_results(resumed)$candidate_id, results$candidate_id)
  expect_identical(Neoantimon2:::run_runtime(resumed)$run_id, p$run_id)
})

test_that("result filters are non-destructive views", {
  p <- Neoantimon2:::run_runtime(t004_project(), adapter = t004_adapter())
  all <- neo_results(p)
  expect_equal(nrow(neo_results(p, grade = "NOT_GRADED")), 1L)
  expect_equal(nrow(neo_results(p, prediction_status = "SUCCESS")), 0L)
  expect_equal(nrow(neo_results(p, reason = "NOT_RUN")), 1L)
  expect_equal(nrow(neo_results(p)), nrow(all))
  filtered <- neo_results(p, grade = "NOT_GRADED")
  filtered$final_grade[[1L]] <- "A"
  expect_identical(neo_results(p)$final_grade, "NOT_GRADED")
})

test_that("unexecuted projects have no results; executed empty output has a schema", {
  for (p in list(t004_project(), t004_project(initialize = TRUE))) {
    for (accessor in list(neo_results, neo_events, neo_origins, neo_multisample))
      expect_error(accessor(p), "Results are not available")
    expect_error(neo_explain(p, peptide = "MISSING"), "Results are not available")
    expect_error(neo_explain(p), "Results are not available")
  }
  p <- neo_run(t004_project())
  for (table in c(neoantigen_results = "neo_results", events = "neo_events",
                  peptide_origins = "neo_origins", multisample_summary = "neo_multisample")) {
    result <- get(table)(p)
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 0L)
  }
  expect_identical(names(neo_results(p)), vapply(
    Neoantimon2:::schema_contract()$tables$neoantigen_results$fields, `[[`, "", "name"))
  expect_length(neo_explain(p, peptide = "MISSING")$candidates, 0L)
  expect_true(file.exists(Neoantimon2:::runtime_snapshot_path(p)))
})

test_that("invalid input, optional omission, and requested tool failure stay distinct", {
  missing <- tempfile("missing-input-")
  expect_error(neo_project(samples = data.frame(sample_id = "sample-1",
    patient_id = "patient-1"), config = t004_config(),
    inputs = list(list(input_id = "bad", input_role = "somatic_variants",
      sample_id = "sample-1", path = missing))), "Required input file missing")

  not_run <- Neoantimon2:::run_runtime(t004_project(), adapter = t004_adapter("NOT_RUN"))
  expect_identical(neo_results(not_run)$prediction_status, "NOT_RUN")
  failed <- Neoantimon2:::run_runtime(t004_project(), adapter = t004_adapter("FAILED"))
  expect_identical(neo_results(failed)$prediction_status, "FAILED")
  expect_identical(neo_results(failed)$final_grade, "NOT_GRADED")
  expect_identical(neo_results(failed)$not_evaluated_codes[[1L]],
    list("PREDICTION_FAILED"))

  partial <- Neoantimon2:::run_runtime(t004_project(), adapter = t004_adapter("PARTIAL"))
  expect_identical(neo_results(partial)$prediction_status, "PARTIAL")
  for (p in list(failed, partial))
    expect_identical(Neoantimon2:::read_manifest(p$path)$status, "COMPLETED")
})

test_that("invalid mock output fails without becoming an empty success", {
  p <- t004_project(initialize = TRUE)
  bad <- function(project) {
    x <- t004_adapter()(project)
    x$neoantigen_results[[1L]]$final_grade <- "E"
    x
  }
  expect_error(Neoantimon2:::run_runtime(p, adapter = bad), "Grade/status mismatch")
  manifest <- Neoantimon2:::read_manifest(p$path)
  expect_identical(manifest$status, "FAILED")
  expect_identical(manifest$checkpoints[[1L]]$status, "FAILED")
  expect_false(file.exists(Neoantimon2:::runtime_snapshot_path(p)))

  real <- neo_project(config = t004_config(), fixture_kind = "real")
  expect_error(Neoantimon2:::run_runtime(real, adapter = function(project) NULL),
    "restricted to synthetic/mock")

  cfg <- t004_config()
  input <- file.path(cfg$paths$data_root, "variants.vcf")
  dir.create(dirname(input), recursive = TRUE)
  writeLines("synthetic placeholder", input)
  production <- neo_project(samples = data.frame(sample_id = "sample-1",
    patient_id = "patient-1"), config = cfg,
    inputs = list(list(input_id = "variants", input_role = "somatic_variants",
      sample_id = "sample-1", path = input)), initialize = TRUE)
  expect_error(Neoantimon2:::run_runtime(production), "not implemented in T004")
  expect_identical(Neoantimon2:::read_manifest(production$path)$status, "FAILED")
})

test_that("neo_explain searches reasons and traces every origin and event", {
  p <- Neoantimon2:::run_runtime(t004_project(), adapter = t004_adapter())
  id <- neo_results(p)$candidate_id[[1L]]
  explanation <- neo_explain(p, candidate_id = id)
  expect_s3_class(explanation, "neo_explanation")
  expect_length(explanation$candidates, 1L)
  expect_length(explanation$candidates[[1L]]$origins, 1L)
  expect_length(explanation$candidates[[1L]]$events, 1L)
  expect_length(neo_explain(p, reason = "NOT_RUN")$candidates, 1L)
  expect_match(neo_explain(p, peptide = "ACDEFGHIK", format = "text"), id,
    fixed = TRUE)
  expect_identical(neo_explain(p, peptide = "MISSING", format = "text"),
    "No matching candidate.")
  expect_error(neo_explain(p), "query")
})

test_that("exports exactly the approved current milestone API", {
  expect_identical(names(formals(neo_run)), "project")
  exports <- getNamespaceExports("Neoantimon2")
  expect_setequal(exports, current_public_api)
  expect_length(setdiff(exports, approved_public_api), 0L)
})

test_that("real no-input execution fails explicitly", {
  p <- neo_project(config = t004_config(), fixture_kind = "real", initialize = TRUE)
  expect_error(neo_run(p), "not implemented in T004")
  expect_identical(Neoantimon2:::read_manifest(p$path)$status, "FAILED")
  expect_error(neo_results(p), "Results are not available")
})

test_that("the authoritative run_manifest cannot be substituted", {
  p <- t004_project(initialize = TRUE)
  divergent <- function(project) {
    x <- t004_adapter()(project)
    x$run_manifest$package_version <- "99.0.0"
    Neoantimon2:::validate_snapshot(x)
    x
  }
  expect_error(Neoantimon2:::run_runtime(p, adapter = divergent),
    "run_manifest differs")
  expect_identical(Neoantimon2:::read_manifest(p$path)$status, "FAILED")
  expect_false(file.exists(Neoantimon2:::runtime_snapshot_path(p)))

  p <- neo_run(t004_project())
  path <- Neoantimon2:::runtime_snapshot_path(p)
  x <- jsonlite::read_json(path, simplifyVector = FALSE)
  x$run_manifest$package_version <- "99.0.0"
  Neoantimon2:::atomic_json(x, path)
  expect_error(neo_results(p), "run_manifest differs")
  expect_error(neo_run(p), "run_manifest differs")
  expect_identical(Neoantimon2:::read_manifest(p$path)$status, "COMPLETED")
})

test_that("missing, corrupt, and valid-but-altered completed output is detected", {
  for (damage in c("missing", "corrupt", "altered")) {
    p <- Neoantimon2:::run_runtime(t004_project(), adapter = t004_adapter())
    path <- Neoantimon2:::runtime_snapshot_path(p)
    if (damage == "missing") unlink(path)
    if (damage == "corrupt") writeLines("{broken", path)
    if (damage == "altered") {
      x <- jsonlite::read_json(path, simplifyVector = FALSE)
      x$neoantigen_results <- list()
      Neoantimon2:::validate_snapshot(x)
      Neoantimon2:::atomic_json(x, path)
    }
    expect_error(neo_results(p))
    expect_error(neo_explain(p, peptide = "MISSING"))
    expect_error(neo_run(p))
    expect_identical(Neoantimon2:::read_manifest(p$path)$status, "COMPLETED")
  }
})

test_that("reuse integrity errors after RUNNING transition mark execution FAILED", {
  for (damage in c("missing", "corrupt", "manifest")) {
    p <- t004_project(initialize = TRUE)
    Neoantimon2:::update_run_status(p, "RUNNING")
    fp <- Neoantimon2:::runtime_fingerprints(p)
    Neoantimon2:::checkpoint_update(p, "pipeline", fingerprints = fp)
    Neoantimon2:::checkpoint_update(p, "pipeline", "RUNNING", fp)
    path <- Neoantimon2:::runtime_snapshot_path(p)
    x <- Neoantimon2:::empty_runtime_snapshot(p)
    if (damage == "manifest") x$run_manifest$package_version <- "99.0.0"
    Neoantimon2:::atomic_json(x, path)
    Neoantimon2:::checkpoint_update(p, "pipeline", "COMPLETED", fp,
      output_records = list(Neoantimon2:::runtime_output_record(p)))
    Neoantimon2:::update_run_status(p, "INTERRUPTED")
    if (damage == "missing") unlink(path)
    if (damage == "corrupt") writeLines("{broken", path)
    called <- FALSE
    expect_error(Neoantimon2:::run_runtime(p, adapter = function(project) {
      called <<- TRUE
      stop("must not recompute")
    }))
    expect_false(called)
    expect_identical(Neoantimon2:::read_manifest(p$path)$status, "FAILED")
  }
})

test_that("intact completed output reopens without recomputation or manifest mutation", {
  p <- Neoantimon2:::run_runtime(t004_project(), adapter = t004_adapter())
  before <- Neoantimon2:::read_manifest(p$path)
  reopened <- neo_project(path = p$path)
  expect_identical(neo_results(reopened), neo_results(p))
  resumed <- Neoantimon2:::run_runtime(reopened, adapter = function(project)
    stop("must not recompute"))
  expect_identical(resumed$run_id, p$run_id)
  expect_identical(Neoantimon2:::read_manifest(p$path), before)
})
