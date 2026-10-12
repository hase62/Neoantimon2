external_converter_runtime <- function(path) {
  if (grepl("[.]R$", path, ignore.case = TRUE)) "R" else
    if (grepl("[.]py$", path, ignore.case = TRUE)) "Python" else
      external_error("CONVERTER_RUNTIME", "Converter must be an R or Python script")
}

external_converter_interpreter <- function(runtime) {
  path <- if (runtime == "R") file.path(R.home("bin"), "Rscript") else Sys.which("python3")
  if (!length(path) || !nzchar(path) || !file.exists(path))
    external_error("INTERPRETER_MISSING", paste(runtime, "interpreter is unavailable"))
  path
}

external_convert <- function(project, input, converter, bundle, timeout = 60, registry_pin = NULL) {
  if (!is.numeric(timeout) || length(timeout) != 1L || !is.finite(timeout) || timeout < 1 || timeout != floor(timeout))
    external_error("INVALID_TIMEOUT", "Require a finite positive integer timeout in seconds")
  for (p in c(input, converter)) if (!file.exists(p) || dir.exists(p))
    external_error("CONVERTER_INPUT_MISSING", "Raw input or converter script not found")
  runtime <- external_converter_runtime(converter)
  interpreter <- external_converter_interpreter(runtime)
  version <- paste(suppressWarnings(system2(interpreter, "--version", stdout = TRUE, stderr = TRUE, timeout = 10)), collapse = "\n")
  identity <- list(contract = "external-converter-1.0.0", schema_version = external_version(),
    raw_input_sha256 = external_sha(input), converter_sha256 = external_sha(converter),
    runtime = runtime, interpreter_version = version, interpreter_sha256 = external_sha(interpreter),
    bundle_id = bundle$manifest$bundle_id, timeout = as.numeric(timeout), registry_pin = registry_pin)
  id <- infra_hash(identity)
  parent <- file.path(project$path, "cache", "external_conversions")
  if (!within_path(resolve_path(parent), resolve_path(project$path)))
    external_error("UNSAFE_PATH", "Conversion cache escapes run")
  target <- file.path(parent, id)
  if (file.exists(file.path(target, "manifest.json"))) {
    old <- jsonlite::read_json(file.path(target, "manifest.json"))
    if (!identical(infra_hash(old$identity), id)) external_error("INTEGRITY_FAILURE", "Conversion identity mismatch")
    envelope <- old; envelope$metadata_digest <- NULL
    if (!identical(old$metadata_digest, infra_hash(envelope))) external_error("INTEGRITY_FAILURE", "Conversion envelope mismatch")
    external_verify_files(target, old$files)
    external_read_predictions(file.path(target, "standard.tsv"))
    return(list(path = file.path(target, "standard.tsv"), metadata = old))
  }
  if (file.exists(target)) external_error("INTEGRITY_FAILURE", "Incomplete conversion cache")
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(".conversion-", tmpdir = parent); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  script_name <- if (runtime == "R") "converter.R" else "converter.py"
  raw_path <- file.path(tmp, "input-X")
  script <- file.path(tmp, script_name); output <- file.path(tmp, "standard.tsv")
  if (!file.copy(input, raw_path) || !file.copy(converter, script) ||
      external_sha(raw_path) != identity$raw_input_sha256 || external_sha(script) != identity$converter_sha256)
    external_error("INTEGRITY_FAILURE", "Input/script changed while staging conversion")
  args <- c(shQuote(script), "--input", shQuote(raw_path), "--output", shQuote(output), "--bundle", shQuote(bundle$path))
  stdout <- file.path(tmp, "stdout.txt"); stderr <- file.path(tmp, "stderr.txt")
  status <- tryCatch(suppressWarnings(system2(interpreter, args, stdout = stdout, stderr = stderr,
    timeout = timeout)), error = function(e) {
      writeLines(conditionMessage(e), stderr); 1L
    })
  if (!file.exists(stdout)) file.create(stdout)
  if (!file.exists(stderr)) file.create(stderr)
  problem <- if (status == 124L) "CONVERTER_TIMEOUT" else if (status != 0L) "CONVERTER_FAILED" else
    if (!file.exists(output) || dir.exists(output)) "CONVERTER_OUTPUT_MISSING" else NULL
  if (is.null(problem)) problem <- tryCatch({ external_read_predictions(output); NULL },
    error = function(e) paste0("CONVERTER_INVALID_OUTPUT: ", conditionMessage(e)))
  meta <- list(identity = identity, execution = list(exit_status = as.integer(status),
    command = as.list(c(runtime, script_name, "--input", "INPUT_X", "--output", "STANDARD_TSV", "--bundle", "EXPORT_BUNDLE")),
    stdout = paste(readLines(stdout, warn = FALSE), collapse = "\n"),
    stderr = paste(readLines(stderr, warn = FALSE), collapse = "\n"), problem = problem),
    files = external_records(tmp, c(script_name, "input-X", "stdout.txt", "stderr.txt",
      if (file.exists(output) && !dir.exists(output)) "standard.tsv" else character())))
  meta$metadata_digest <- infra_hash(meta)
  atomic_json(meta, file.path(tmp, "manifest.json"))
  if (!is.null(problem)) {
    failed <- file.path(parent, paste0("failed-", basename(tmp)))
    if (!file.rename(tmp, failed)) stop("Could not store conversion diagnostics")
    external_error(sub(":.*", "", problem), paste("See local diagnostics at", failed))
  }
  if (!file.rename(tmp, target)) stop("Could not atomically store conversion cache")
  list(path = file.path(target, "standard.tsv"), metadata = meta)
}
