allo_fixture <- function(file = "e2e.json") jsonlite::read_json(testthat::test_path("..", "fixtures", "allo_foundation", file))
allo_test_input <- function(ambiguous = TRUE) {
  x <- allo_fixture()
  x$reference$path <- normalizePath(testthat::test_path("..", "fixtures", "allo_foundation", "reference.json"))
  if (!ambiguous) {
    x$typing <- Filter(function(t) t$locus == "A", x$typing)
    x$self_search[[1L]]$required_loci <- list("A")
  }
  x
}
allo_test_config <- function(x = allo_test_input(), base = NULL) Neoantimon2:::allo_config(x, base)
allo_test_build <- function(x = allo_test_input()) Neoantimon2:::allo_foundation(allo_test_config(x))
allo_test_manifest <- function(m, x = allo_test_input(FALSE)) {
  f <- tempfile(fileext = ".json")
  jsonlite::write_json(m, f, pretty = TRUE, auto_unbox = TRUE, null = "null", digits = NA)
  x$reference <- list(path = f)
  x
}
allo_test_typing <- function(allele, complete = TRUE) list(subject_id = "patient-1", typing_scenario_id = "typing:1",
  locus = sub("\\*.*", "", sub("^HLA-", "", toupper(trimws(allele)))), alleles = list(allele), complete = complete, copy_number = 1L, phase = NULL)
allo_one_resolution <- function(allele, manifest = allo_fixture("reference.json"), complete = TRUE) {
  x <- allo_test_input(FALSE); x$typing <- list(allo_test_typing(allele, complete))
  x$self_search <- list()
  allo_test_build(allo_test_manifest(manifest, x))$typing[[1L]]
}
allo_test_enable <- function(p, input = allo_test_input()) {
  p$config$allo <- allo_test_config(input)
  p$allo_foundation <- Neoantimon2:::allo_foundation(p$config$allo)
  p
}
allo_frozen_surface <- function(p) list(results = neo_results(p), events = neo_events(p), origins = neo_origins(p),
  multisample = neo_multisample(p), explanation = neo_explain(p, candidate_id = neo_results(p)$candidate_id),
  longitudinal = Neoantimon2:::longitudinal_read(p),
  vaccine = if (Neoantimon2:::vaccine_enabled(p)) Neoantimon2:::vaccine_read(p) else NULL)
allo_test_validator <- function(schema = "allo-foundation", definition = "foundation_input") {
  path <- system.file("schema", paste0(schema, "-1.0.0.json"), package = "Neoantimon2")
  s <- jsonlite::read_json(path); s$`$ref` <- paste0("#/$defs/", definition)
  validator <- jsonvalidate::json_validator(jsonlite::toJSON(s, auto_unbox = TRUE, null = "null"), engine = "ajv")
  function(x) as.logical(validator(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", digits = NA)))
}
