test_that("native no-predictor export is automatic, deterministic and verified", {
  p <- peptide_test_run()
  b <- neo_export_bundle(p)
  expect_identical(b$manifest$contract, "external-export-1.0.0")
  expect_identical(neo_export_bundle(p)$manifest, b$manifest)
  expect_identical(neo_export_bundle(neo_project(path = p$path))$manifest, b$manifest)
  copy <- neo_export_bundle(p, tempfile("export-copy-"))
  expect_identical(copy$manifest, b$manifest)
  expect_setequal(vapply(b$graph$candidates, `[[`, "", "candidate_id"), neo_results(p)$candidate_id)
  expect_equal(length(b$graph$links), sum(vapply(b$graph$candidates,
    function(c) length(c$origin_ids), 0L)))
  expect_gt(length(b$graph$annotations$core_peptides), 0)
  expect_true(all(vapply(b$graph$annotations$core_peptides,
    function(p) length(p$provenance_occurrences) > 0, TRUE)))
  writeLines("tampered", file.path(copy$path, "peptides.tsv"))
  expect_error(neo_export_bundle(p, copy$path), "INTEGRITY_FAILURE")
  expect_error(neo_export_bundle(p, file.path(p$path, "intermediate", "export")), "UNSAFE_PATH")
})

test_that("legacy unspecified presenters have explicit query coverage without inferred HLA", {
  cfg <- neo_config(env = character(), paths = list(run_root = tempfile("legacy-runs-")))
  p <- neo_project(config = cfg, fixture_kind = "synthetic",
    samples = data.frame(sample_id = "sample-1", patient_id = "patient-1"),
    analysis_context = fixture("simple_somatic")$run_manifest$analysis_contexts)
  p <- Neoantimon2:::run_runtime(p, function(p) {
    x <- fixture("simple_somatic"); x$run_manifest <- Neoantimon2:::read_manifest(p$path)$run_manifest
    for (table in c("events", "peptide_origins", "neoantigen_results"))
      for (i in seq_along(x[[table]])) x[[table]][[i]]$analysis_id <- p$analysis_id
    x
  })
  b <- neo_export_bundle(p)
  expect_equal(b$manifest$unsupported_presenter_queries, 1L)
  expect_equal(nrow(read.delim(file.path(b$path, "pmhc_queries.tsv"))), 0L)
  expect_identical(read.delim(file.path(b$path, "query_coverage.tsv"))$pmhc_query_status,
    "UNSUPPORTED_PRESENTER")
  expect_length(b$graph$candidates, 1L)
})

test_that("empty bundles have headers and preserve literal U", {
  cfg <- neo_config(env = character(), paths = list(run_root = tempfile("empty-runs-")))
  p <- neo_run(neo_project(config = cfg, fixture_kind = "synthetic"))
  b <- neo_export_bundle(p)
  expect_length(b$graph$candidates, 0)
  for (f in c("peptides", "pmhc_queries", "candidates", "origins", "events", "candidate_origin_links"))
    expect_equal(nrow(utils::read.delim(file.path(b$path, paste0(f, ".tsv")))), 0)
  sec <- neo_export_bundle(peptide_test_run(sec = TRUE))
  expect_true(any(grepl("U", vapply(sec$graph$candidates, `[[`, "", "peptide_sequence"), fixed = TRUE)))
})
