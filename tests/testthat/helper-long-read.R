long_read_fixture_input <- function(file = "models.gtf", role = "long_read_gtf",
                                    id = "model", sample = "sample-1") list(input_id = id,
  input_role = role, sample_id = sample, path = normalizePath(testthat::test_path(
    "..", "fixtures", "long_read", file)), checksum = NULL, version = "synthetic-v1")
long_read_fixture_project <- function(inputs = list(long_read_fixture_input()),
                                       hla = list("patient-1" = c("HLA-A*02:01",
                                         "HLA-DPA1*01:03~HLA-DPB1*04:01")),
                                       samples = data.frame(sample_id = "sample-1", patient_id = "patient-1")) {
  base <- reconstruction_test_project()
  values <- unclass(base$config)
  root <- normalizePath(testthat::test_path("..", "fixtures", "long_read", "reference"))
  for (n in names(values$resources)) values$resources[[n]]$path <- file.path(root,
    basename(values$resources[[n]]$path))
  cfg <- do.call(neo_config, c(values, list(env = character())))
  neo_project(samples = samples, inputs = inputs, hla = hla, config = cfg,
    path = file.path(cfg$paths$run_root, "long-read-test"))
}
long_read_fixture_construct <- function(inputs = list(long_read_fixture_input())) {
  p <- long_read_fixture_project(inputs)
  x <- Neoantimon2:::long_read_extension_construct(p)
  Neoantimon2:::long_read_validate(x, p)
  x
}
long_read_test_event <- function(x, tx) {
  o <- Filter(function(o) o$provider_transcript_id == tx, x$observations)[[1L]]
  Filter(function(e) e$long_read_source_event_id == o$long_read_source_event_id, x$events)[[1L]]
}
long_read_test_orfs <- function(x, tx, type = "reference_anchored_orf") {
  e <- long_read_test_event(x, tx)
  lapply(Filter(function(h) h$long_read_source_event_id == e$long_read_source_event_id &&
    h$orf$evidence_type == type, x$hypotheses), `[[`, "orf")
}
