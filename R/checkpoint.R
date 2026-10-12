# Checkpoints live in the same atomic envelope as their summary, avoiding a
# two-file transaction. checkpoints/ is reserved for future step payloads.
validate_checkpoint <- function(cp) {
  fields <- c("step_id", "status", "input_fingerprint", "config_fingerprint", "dependency_fingerprint",
    "output_records", "started_at", "completed_at", "error")
  strict_keys(cp, fields, "checkpoint")
  if (length(setdiff(fields, names(cp)))) stop("Incomplete checkpoint")
  for (n in c("step_id", "input_fingerprint", "config_fingerprint", "dependency_fingerprint")) scalar_text(cp[[n]], n)
  if (!grepl("^[A-Za-z0-9_.-]+$", cp$step_id)) stop("Invalid step_id")
  if (!cp$status %in% c("PENDING", "RUNNING", "COMPLETED", "FAILED", "INVALIDATED")) stop("Invalid checkpoint status")
  if (!is.list(cp$output_records) || !is.null(names(cp$output_records))) stop("Invalid output records")
  for (o in cp$output_records) {
    strict_keys(o, c("path", "checksum"), "checkpoint output")
    scalar_text(o$path, "output path")
    if (grepl("^[/~]|^[A-Za-z]:|(^|/)\\.\\.(/|$)", o$path)) stop("Output must be run-relative")
    parts <- strsplit(o$path, "/", fixed = TRUE)[[1]]
    if (!grepl("^(intermediate|results|checkpoints)/.+", o$path) ||
        grepl("\\", o$path, fixed = TRUE) || any(parts %in% c("", ".", "..")) ||
        endsWith(o$path, "/")) stop("Invalid checkpoint output namespace")
    resource_spec(list(checksum = o$checksum))
    if (is.null(o$checksum)) stop("Output checksum required")
  }
  for (n in c("started_at", "completed_at")) if (!is.null(cp[[n]])) scalar_text(cp[[n]], n)
  if (cp$status == "COMPLETED" && (is.null(cp$started_at) || is.null(cp$completed_at))) stop("Missing checkpoint timestamps")
  if (!is.null(cp$error)) {
    strict_keys(cp$error, c("code"), "checkpoint error")
    scalar_text(cp$error$code, "error code")
    if (!grepl("^[A-Z][A-Z0-9_]*$", cp$error$code)) stop("Use a non-secret error code")
  }
  if (cp$status == "FAILED" && is.null(cp$error)) stop("Failed checkpoint requires error metadata")
  invisible(TRUE)
}

checkpoint_fingerprints <- function(inputs, config, dependencies) {
  # Callers supply only the scientific slice used by this step.
  list(input_fingerprint = infra_hash(inputs), config_fingerprint = infra_hash(config),
       dependency_fingerprint = infra_hash(dependencies))
}

checkpoint_update <- function(project, step_id, status = "PENDING", fingerprints,
                              output_records = list(), error = NULL) {
  strict_keys(fingerprints, c("input_fingerprint", "config_fingerprint", "dependency_fingerprint"), "fingerprints")
  if (length(fingerprints) != 3L) stop("Three fingerprints required")
  mutate_manifest(project, function(m) {
    if (m$status == "COMPLETED") stop("Completed run is immutable")
    ids <- vapply(m$checkpoints, `[[`, "", "step_id")
    index <- match(step_id, ids)
    old <- if (is.na(index)) NULL else m$checkpoints[[index]]
    transitions <- list(PENDING = c("RUNNING", "INVALIDATED"), RUNNING = c("COMPLETED", "FAILED", "INVALIDATED"),
      FAILED = c("RUNNING", "INVALIDATED"), COMPLETED = "INVALIDATED", INVALIDATED = "PENDING")
    if (is.null(old)) {
      if (status != "PENDING") stop("New checkpoint must be PENDING")
    } else {
      if (!status %in% transitions[[old$status]]) stop("Invalid checkpoint transition")
      matches <- all(vapply(names(fingerprints), function(k) identical(fingerprints[[k]], old[[k]]), TRUE))
      if (!matches && old$status != "INVALIDATED") stop("Checkpoint fingerprint mismatch; invalidate before replacement")
    }
    cp <- c(list(step_id = step_id, status = status), fingerprints,
      list(output_records = output_records,
        started_at = if (status == "RUNNING") utc_now() else old$started_at,
        completed_at = if (status == "COMPLETED") utc_now() else NULL, error = error))
    validate_checkpoint(cp)
    if (status == "COMPLETED") verify_checkpoint_outputs(cp, project$path)
    if (is.na(index)) m$checkpoints[[length(m$checkpoints) + 1L]] <- cp else m$checkpoints[[index]] <- cp
    m
  })
}
verify_checkpoint_outputs <- function(cp, path) {
  for (o in cp$output_records) {
    p <- resolve_path(o$path, path)
    namespace <- file.path(path, strsplit(o$path, "/", fixed = TRUE)[[1]][1])
    if (!within_path(p, namespace) || !file.exists(p) || dir.exists(p)) stop("Checkpoint output missing/unsafe")
    integrity_record(p, "sha256", o$checksum)
  }
  invisible(TRUE)
}
checkpoint_reusable <- function(project, step_id, fingerprints) {
  if (is.null(project$path)) stop("Run is not initialized")
  check_run_location(project$path, project$config)
  with_run_lock(project$path, function() {
    m <- read_manifest(project$path)
    if (m$run_id != project$run_id || m$run_manifest$analysis_id != project$analysis_id) stop("Stale project identity")
    cp <- Filter(function(x) identical(x$step_id, step_id), m$checkpoints)
    if (!length(cp)) return(FALSE)
    strict_keys(fingerprints, c("input_fingerprint", "config_fingerprint", "dependency_fingerprint"), "fingerprints")
    if (length(fingerprints) != 3L || !all(vapply(names(fingerprints), function(k)
      identical(fingerprints[[k]], cp[[1]][[k]]), TRUE))) stop("Checkpoint fingerprint mismatch; cannot reuse")
    if (cp[[1]]$status != "COMPLETED") return(FALSE)
    verify_checkpoint_outputs(cp[[1]], project$path)
    TRUE
  })
}
