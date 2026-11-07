autorun_setup <- function(mode = "good") {
  data <- indirect_input(); root <- tempfile("indirect autorun ; literal "); dir.create(root)
  tool <- file.path(root, "mock tool ; literal")
  source <- system.file("extdata", "indirect-autorun", "mock-netmhciipan.py", package = "Neoantimon2")
  file.copy(source, tool); Sys.chmod(tool, "0755")
  model <- file.path(root, "model ; literal"); dir.create(file.path(model, "data"), recursive = TRUE)
  jsonlite::write_json(indirect_fixture("binding.json")$rows, file.path(model, "data/binding-map.json"),
    auto_unbox = TRUE, null = "null", digits = NA)
  writeLines(mode, file.path(model, "mode.txt"))
  list(predictor_version = "3.2", allo = data$a, presenters = vapply(data$b$presenters, `[[`, "", "presenter"),
    reference = list(path = data$a$reference$path, sha256 = digest::digest(file = data$a$reference$path, algo = "sha256", serialize = FALSE)),
    executable = tool, model = list(path = model, identity = "controlled mock data; not vendor models"),
    output_dir = file.path(root, "run ; literal"),
    authorization = list(status = "mock_only", terms = "Synthetic test only, not real model authorization"))
}
autorun_call <- function(x = autorun_setup(), ...) do.call(neo_indirect_autorun, c(x, list(...)))
autorun_ctx <- function(x) x$result$contexts[[1L]]
