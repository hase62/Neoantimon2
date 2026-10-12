direction_config <- function() neo_config(env = character(), paths = list(
  run_root = file.path(tempdir(), "direction-runs"),
  data_root = file.path(tempdir(), "direction-data")))

test_that("somatic direction is explicit and preserves scientific set identity", {
  samples <- data.frame(sample_id = c("s3", "s2", "s1"),
    patient_id = c("p2", "p1", "p1"))
  cfg <- direction_config()
  p <- neo_project(samples = samples, config = cfg)
  contexts <- p$analysis_context
  expect_identical(vapply(contexts, `[[`, "", "direction"),
    rep("within_patient", 2))
  explicit <- neo_project(samples = samples, config = cfg,
    analysis_context = contexts)
  expect_identical(explicit$analysis_id, p$analysis_id)
  expect_identical(explicit$config_fingerprint, p$config_fingerprint)
  reordered <- lapply(rev(contexts), function(x) {
    x$sample_bindings <- rev(x$sample_bindings)
    x
  })
  expect_identical(neo_project(samples = samples[3:1, ], config = cfg,
    analysis_context = reordered)$analysis_id, p$analysis_id)
  for (direction in c("tumor_to_patient", "donor_to_recipient",
                      "recipient_to_donor", "unknown")) {
    invalid <- contexts
    invalid[[1]]$direction <- direction
    expect_error(neo_project(samples = samples, config = cfg,
      analysis_context = invalid), "Invalid somatic direction")
    expect_identical(invalid[[1]]$direction, direction)
  }
  # Self may be patient rather than normal, without changing subject identity.
  contexts[[1]]$self_subject_role <- "patient"
  expect_s3_class(neo_project(samples = samples, config = cfg,
    analysis_context = contexts), "neo_project")
})

test_that("default project contexts pass the frozen T002 snapshot validator", {
  p <- neo_project(samples = data.frame(sample_id = "s1", patient_id = "p1"),
    config = direction_config(), initialize = TRUE, fixture_kind = "synthetic")
  m <- Neoantimon2:::read_manifest(p$path)
  snapshot <- list(run_manifest = m$run_manifest, events = list(),
    peptide_origins = list(), neoantigen_results = list(), multisample_summary = list())
  expect_silent(Neoantimon2:::validate_snapshot(snapshot))
  expect_identical(m$run_manifest$analysis_contexts, p$analysis_context)
  expect_identical(neo_project(path = p$path)$analysis_id, p$analysis_id)
})

test_that("stored legacy somatic direction rejects resume without rewriting artifacts", {
  p <- neo_project(samples = data.frame(sample_id = "s1", patient_id = "p1"),
    config = direction_config(), initialize = TRUE)
  # A sealed development envelope carrying the old value must fail at context
  # validation, before resume identity comparisons or any metadata mutation.
  m <- Neoantimon2:::read_manifest(p$path)
  m$run_manifest$analysis_contexts[[1]]$direction <- "tumor_to_patient"
  Neoantimon2:::write_manifest(m, p$path)
  path <- Neoantimon2:::manifest_file(p$path)
  before <- readBin(path, "raw", n = file.info(path)$size)
  expect_error(neo_project(path = p$path), "Invalid somatic direction")
  expect_identical(readBin(path, "raw", n = file.info(path)$size), before)
  expect_true(dir.exists(p$path))
  expect_false(dir.exists(file.path(p$path, ".metadata-lock")))
})
