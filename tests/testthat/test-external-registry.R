external_test_registry <- function() {
  root <- tempfile("registry-"); dir.create(root)
  script <- system.file("extdata", "external-example-converter.R", package = "Neoantimon2")
  file.copy(script, file.path(root, "converter.R"))
  entry <- list(id = "synthetic-example@1.0.0", runtime = "R", script = "converter.R",
    sha256 = Neoantimon2:::external_sha(script),
    metadata = list(output_schema = "external-prediction-1.0.0", description = "Synthetic example",
      usage = "Rscript converter.R --input X --output Y --bundle B",
      arguments = list("--input", "--output", "--bundle"), license = "MIT"),
    supported_predictor_versions = list("synthetic-model@1"), source_revision = "local-fixture-v1",
    issue_url = NULL, ok = 0L, ng = 0L)
  index <- file.path(root, "index.json")
  Neoantimon2:::atomic_json(list(schema_version = "1.0.0", converters = list(entry)), index)
  list(index = index, root = root, entry = entry)
}

test_that("registry IDs require trust, verify pins, and reuse offline", {
  p <- peptide_test_run(); q <- neo_export_bundle(p)$graph$candidates[[1]]
  reg <- external_test_registry()
  raw <- external_test_write(data.frame(peptide = q$peptide_sequence, hla = q$mhc_key, value = "0.8"))
  call <- function(trust = FALSE) neo_analyze_external(p, raw, list(mode = "evidence_only"),
    converter_id = "synthetic-example@1.0.0", registry = reg$index, trust = trust)
  expect_error(call(), "TRUST_REQUIRED")
  expect_false(dir.exists(file.path(p$path, "cache", "external_conversions")))
  a <- call(TRUE)
  expect_identical(Neoantimon2:::external_analysis_read(a)$manifest$fingerprint$converter$registry_pin$id,
    "synthetic-example@1.0.0")
  # Loss of the index and source script cannot block an already approved pin.
  unlink(reg$root, recursive = TRUE)
  expect_identical(call()$analysis_id, a$analysis_id)
  pins <- list.dirs(file.path(p$path, "cache", "external_registry"), recursive = FALSE)
  writeLines("tampered", file.path(pins[[1]], "converter.R"))
  expect_error(call(), "CHECKSUM_MISMATCH")
})

test_that("registry parser rejects bad integrity, schema, ID and path without execution", {
  cfg <- neo_config(env = character(), paths = list(run_root = tempfile("registry-runs-")))
  p <- neo_run(neo_project(config = cfg, fixture_kind = "synthetic"))
  reg <- external_test_registry()
  resolve <- function(id = "synthetic-example@1.0.0") Neoantimon2:::external_registry_resolve(p, id, reg$index)
  expect_error(resolve("missing@1.0.0"), "REGISTRY_ID_MISSING")
  index <- list(schema_version = "1.0.0", converters = list(reg$entry))
  index$converters[[1]]$sha256 <- paste(rep("0", 64), collapse = "")
  Neoantimon2:::atomic_json(index, reg$index)
  expect_error(resolve(), "CHECKSUM_MISMATCH")
  index$converters[[1]]$metadata$output_schema <- "2.0.0"
  expect_error(Neoantimon2:::external_registry_validate(index), "REGISTRY_SCHEMA")
  index$converters[[1]] <- reg$entry; index$converters[[1]]$script <- "../converter.R"
  expect_error(Neoantimon2:::external_registry_validate(index), "UNSAFE_PATH")
  index$converters[[1]] <- reg$entry; index$converters[[2]] <- reg$entry
  expect_error(Neoantimon2:::external_registry_validate(index), "REGISTRY_ID")
  expect_error(Neoantimon2:::external_registry_resolve(p, "bad", reg$index), "REGISTRY_ID")
  expect_error(Neoantimon2:::external_registry_resolve(p, "synthetic-example@1.0.0", reg$index, NA), "TRUST_REQUIRED")
})

test_that("mutable feedback remains compatible and cannot authorize or replace a converter pin", {
  cfg <- neo_config(env = character(), paths = list(run_root = tempfile("feedback-runs-")))
  p <- neo_run(neo_project(config = cfg, fixture_kind = "synthetic"))
  reg <- external_test_registry()
  resolve <- function(trust = FALSE) Neoantimon2:::external_registry_resolve(p, reg$entry$id, reg$index, trust)
  expect_error(resolve(), "TRUST_REQUIRED")
  entry <- reg$entry; entry$issue_url <- "https://github.com/example/converters/issues/42"
  entry$ok <- 12L; entry$ng <- 3L
  index <- list(schema_version = "1.0.0", converters = list(entry))
  expect_identical(Neoantimon2:::external_registry_validate(index), index)
  Neoantimon2:::atomic_json(index, reg$index)
  expect_error(resolve(), "TRUST_REQUIRED")
  pin <- resolve(TRUE)
  expect_identical(Neoantimon2:::external_sha(pin$path), reg$entry$sha256)
  expect_null(pin$provenance$entry$issue_url)
  expect_equal(pin$provenance$entry$ok, 0L)
  # A second run sees the newer discovery snapshot but still needs explicit trust.
  other <- neo_run(neo_project(config = cfg, fixture_kind = "synthetic"))
  expect_error(Neoantimon2:::external_registry_resolve(other, reg$entry$id, reg$index), "TRUST_REQUIRED")
  newer <- Neoantimon2:::external_registry_resolve(other, reg$entry$id, reg$index, TRUE)
  expect_identical(newer$provenance$entry$issue_url, entry$issue_url)
  expect_equal(newer$provenance$entry$ok, 12L)
  expect_identical(Neoantimon2:::external_sha(newer$path), reg$entry$sha256)
})
