external_registry_validate <- function(index) {
  reject_duplicate_keys(index)
  strict_keys(index, c("schema_version", "converters"), "converter registry")
  if (!identical(index$schema_version, external_version()) || !is.list(index$converters) ||
      !is.null(names(index$converters))) external_error("REGISTRY_SCHEMA", "Invalid registry schema/array")
  ids <- character()
  for (entry in index$converters) {
    fields <- c("id", "runtime", "script", "sha256", "metadata", "supported_predictor_versions",
      "source_revision", "issue_url", "ok", "ng")
    strict_keys(entry, fields, "converter entry")
    if (!setequal(names(entry), fields)) external_error("REGISTRY_SCHEMA", "Incomplete registry entry")
    for (key in c("id", "runtime", "script", "sha256", "source_revision")) scalar_text(entry[[key]], key)
    if (!grepl("^[a-z0-9][a-z0-9-]*@[0-9]+[.][0-9]+[.][0-9]+$", entry$id))
      external_error("REGISTRY_ID", "Require an immutable name@major.minor.patch ID")
    if (!entry$runtime %in% c("R", "Python") || !grepl("^[a-f0-9]{64}$", entry$sha256))
      external_error("REGISTRY_SCHEMA", "Invalid runtime or SHA256")
    relative <- grepl("^[A-Za-z0-9_./@-]+$", entry$script) && !startsWith(entry$script, "/") &&
      !any(strsplit(entry$script, "/", fixed = TRUE)[[1L]] %in% c("", ".", ".."))
    if (!relative && !grepl("^https://[^[:space:]]+$", entry$script))
      external_error("UNSAFE_PATH", "Script must be HTTPS or a safe relative path")
    if (!identical(external_converter_runtime(entry$script), entry$runtime))
      external_error("REGISTRY_SCHEMA", "Script extension/runtime mismatch")
    strict_keys(entry$metadata, c("output_schema", "description", "usage", "arguments", "license"), "converter metadata")
    if (!identical(entry$metadata$output_schema, "external-prediction-1.0.0"))
      external_error("REGISTRY_SCHEMA", "Incompatible converter output schema")
    for (key in c("description", "usage", "license")) scalar_text(entry$metadata[[key]], key)
    if (!identical(unlist(entry$metadata$arguments, use.names = FALSE), c("--input", "--output", "--bundle")))
      external_error("REGISTRY_SCHEMA", "Unsupported converter arguments")
    if (!is.list(entry$supported_predictor_versions) || !is.null(names(entry$supported_predictor_versions)) ||
        any(!vapply(entry$supported_predictor_versions, function(v) is.character(v) && length(v) == 1L && !is.na(v) && nzchar(v), TRUE)))
      external_error("REGISTRY_SCHEMA", "Predictor versions must be a string array")
    if (!is.null(entry$issue_url) && (!is.character(entry$issue_url) || length(entry$issue_url) != 1L ||
        !grepl("^https://github[.]com/[^/]+/[^/]+/issues/[0-9]+$", entry$issue_url)))
      external_error("REGISTRY_SCHEMA", "Issue URL must name a GitHub issue, or be null before publication")
    for (key in c("ok", "ng")) if (!is.null(entry[[key]]) && (!is.numeric(entry[[key]]) ||
      length(entry[[key]]) != 1L || !is.finite(entry[[key]]) || entry[[key]] < 0 || entry[[key]] != floor(entry[[key]])))
      external_error("REGISTRY_SCHEMA", "Invalid feedback count")
    ids <- c(ids, entry$id)
  }
  if (anyDuplicated(ids)) external_error("REGISTRY_ID", "Duplicate immutable IDs")
  index
}

external_fetch <- function(source, destination) {
  if (grepl("^https://", source)) {
    old <- getOption("timeout"); options(timeout = 30); on.exit(options(timeout = old), add = TRUE)
    status <- tryCatch(suppressWarnings(utils::download.file(source, destination, mode = "wb", quiet = TRUE)),
      error = function(e) 1L)
    if (status != 0L || !file.exists(destination)) external_error("REGISTRY_UNAVAILABLE", "Could not fetch registry resource")
  } else {
    if (grepl("^[A-Za-z][A-Za-z0-9+.-]*://", source)) external_error("UNSAFE_URL", "Only HTTPS remote registry resources are allowed")
    if (!file.exists(source) || dir.exists(source) || !file.copy(source, destination))
      external_error("REGISTRY_UNAVAILABLE", "Local registry resource missing")
  }
  invisible(destination)
}

external_registry_resolve <- function(project, id, registry, trust = FALSE) {
  scalar_text(id, "converter_id"); scalar_text(registry, "registry")
  if (!grepl("^[a-z0-9][a-z0-9-]*@[0-9]+[.][0-9]+[.][0-9]+$", id)) external_error("REGISTRY_ID", "Invalid converter ID")
  if (!is.logical(trust) || length(trust) != 1L || is.na(trust)) external_error("TRUST_REQUIRED", "trust must be explicit TRUE or FALSE")
  if (!grepl("^https://", registry)) registry <- resolve_path(registry)
  # Pin by configured index and immutable ID. Reuse needs no network or index.
  key <- infra_hash(list(registry = registry, id = id))
  parent <- file.path(project$path, "cache", "external_registry")
  if (!within_path(resolve_path(parent), resolve_path(project$path))) external_error("UNSAFE_PATH", "Registry cache escapes run")
  root <- file.path(parent, key)
  if (file.exists(file.path(root, "pin.json"))) {
    pin <- jsonlite::read_json(file.path(root, "pin.json"))
    envelope <- pin; envelope$metadata_digest <- NULL
    if (!identical(pin$metadata_digest, infra_hash(envelope)) || pin$id != id || pin$registry != registry)
      external_error("INTEGRITY_FAILURE", "Registry pin envelope mismatch")
    external_registry_validate(list(schema_version = external_version(), converters = list(pin$entry)))
    script <- file.path(root, if (pin$entry$runtime == "R") "converter.R" else "converter.py")
    if (!file.exists(script) || !identical(external_sha(script), pin$entry$sha256))
      external_error("CHECKSUM_MISMATCH", "Cached community converter bytes changed")
  } else {
    if (file.exists(root)) external_error("INTEGRITY_FAILURE", "Incomplete registry cache")
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(".registry-", tmpdir = parent); dir.create(tmp)
    on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
    index_path <- file.path(tmp, "index.json")
    external_fetch(registry, index_path)
    index <- external_registry_validate(jsonlite::read_json(index_path))
    entry <- Filter(function(e) e$id == id, index$converters)
    if (length(entry) != 1L) external_error("REGISTRY_ID_MISSING", "Converter ID is absent from index")
    entry <- entry[[1L]]
    location <- if (grepl("^https://", entry$script)) entry$script else if (grepl("^https://", registry))
      paste0(sub("[^/]*$", "", registry), entry$script) else resolve_path(entry$script, dirname(registry))
    if (!grepl("^https://", location) && !within_path(location, dirname(registry)))
      external_error("UNSAFE_PATH", "Registry script escapes index directory")
    filename <- if (entry$runtime == "R") "converter.R" else "converter.py"
    external_fetch(location, file.path(tmp, filename))
    if (!identical(external_sha(file.path(tmp, filename)), entry$sha256))
      external_error("CHECKSUM_MISMATCH", "Community converter bytes do not match index SHA256")
    pin <- list(id = id, registry = registry, entry = entry, approved = FALSE)
    pin$metadata_digest <- infra_hash(pin)
    atomic_json(pin, file.path(tmp, "pin.json"))
    if (!file.rename(tmp, root)) stop("Could not store registry pin")
    script <- file.path(root, filename)
  }
  if (!isTRUE(pin$approved)) {
    if (!isTRUE(trust)) external_error("TRUST_REQUIRED",
      paste("Review pinned converter", id, "SHA256", pin$entry$sha256,
        "at", script, "; trust=TRUE explicitly authorizes arbitrary local code. Community OK votes are not certification."))
    pin$approved <- TRUE; pin$metadata_digest <- NULL; pin$metadata_digest <- infra_hash(pin)
    atomic_json(pin, file.path(root, "pin.json"))
  }
  list(path = script, provenance = list(id = pin$id, entry = pin$entry,
    registry_identity = infra_hash(list(registry = registry)), approved = TRUE))
}
