utc_now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC")

atomic_json <- function(value, path, before_rename = NULL) {
  if (!dir.exists(dirname(path))) stop("Metadata parent directory missing")
  tmp <- tempfile(".metadata-", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  # Complete serialization precedes opening any file; destination is never truncated.
  json <- jsonlite::toJSON(value, auto_unbox = TRUE, null = "null", digits = NA, pretty = TRUE)
  con <- file(tmp, open = "wb")
  tryCatch({writeBin(charToRaw(paste0(enc2utf8(json), "\n")), con); flush(con)}, finally = close(con))
  jsonlite::read_json(tmp) # reject an incomplete serialization before commit
  if (!is.null(before_rename)) before_rename(tmp, path)
  if (!file.rename(tmp, path)) stop("Atomic metadata rename failed; old destination retained")
  invisible(path)
}

check_run_location <- function(path, config) {
  path <- resolve_path(path)
  root <- resolve_path(config$paths$run_root)
  if (!within_path(path, root) || identical(path, root)) stop("Run path must be strictly inside run_root")
  forbidden <- c(config$paths$data_root, config$paths$repository)
  for (p in forbidden) if (!is.null(p)) {
    p <- resolve_path(p)
    if (within_path(path, p) || within_path(p, path)) stop("Run path overlaps data/repository")
  }
  # Detect repositories even if repository configuration is omitted.
  ancestor <- path
  repeat {
    if (dir.exists(file.path(ancestor, ".git")) || file.exists(file.path(ancestor, ".git")))
      stop("Runtime output inside repository is forbidden")
    parent <- dirname(ancestor)
    if (parent == ancestor) break
    ancestor <- parent
  }
  path
}

with_run_lock <- function(path, action) {
  path <- resolve_path(path)
  lock <- file.path(path, ".metadata-lock")
  if (!dir.create(lock, showWarnings = FALSE)) stop("Run is locked; inspect interrupted writer before manual recovery")
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  action()
}

seal_manifest <- function(m) {
  m$metadata_digest <- NULL
  m$metadata_digest <- infra_hash(m)
  m
}
manifest_file <- function(path) file.path(path, "manifest", "run.json")
write_manifest <- function(m, path) {
  validate_manifest(m)
  atomic_json(seal_manifest(m), manifest_file(path))
}
validate_manifest <- function(m) {
  required <- c("manifest_version", "checkpoint_version", "config_contract_version", "run_manifest", "run_id",
    "project_id", "created_at", "updated_at", "status", "config", "config_fingerprint", "input_fingerprint",
    "samples", "hla", "input_records", "resource_records", "checksum_policy", "platform", "git_commit", "checkpoints")
  strict_keys(m, c(required, "metadata_digest"), "execution manifest")
  if (length(setdiff(required, names(m)))) stop("Incomplete execution manifest")
  for (n in c("manifest_version", "checkpoint_version", "config_contract_version"))
    if (!identical(m[[n]], infra_version())) stop("Unsupported manifest/schema version")
  validate_table(list(m$run_manifest), "run_manifest")
  for (n in c("run_id", "project_id", "created_at", "updated_at", "config_fingerprint", "input_fingerprint")) scalar_text(m[[n]], n)
  if (!grepl("^run-[a-f0-9]{32}$", m$run_id)) stop("Invalid run_id")
  if (!m$status %in% c("INITIALIZED", "RUNNING", "COMPLETED", "FAILED", "INTERRUPTED")) stop("Invalid run status")
  if (!is.list(m$checkpoints) || !is.null(names(m$checkpoints))) stop("Invalid checkpoints")
  for (cp in m$checkpoints) validate_checkpoint(cp)
  if (anyDuplicated(vapply(m$checkpoints, `[[`, "", "step_id"))) stop("Duplicate checkpoint")
  invisible(TRUE)
}
read_manifest <- function(path) {
  file <- manifest_file(path)
  if (!file.exists(file)) stop("Run manifest missing")
  m <- jsonlite::read_json(file)
  validate_manifest(m)
  if (is.null(m$metadata_digest) || !identical(m$metadata_digest, seal_manifest(m)$metadata_digest))
    stop("Manifest integrity mismatch")
  m
}

initialize_run <- function(p, path = NULL) {
  # tempfile's process-unique component plus clock/PID produces execution identity,
  # without touching scientific identity or the user's RNG state.
  id <- paste0("run-", substr(infra_hash(list(utc_now(), as.numeric(Sys.getpid()), tempfile())), 1, 32))
  if (is.null(path)) path <- file.path(p$config$paths$run_root, id)
  path <- check_run_location(resolve_path(path, p$config$paths$run_root), p$config)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!dir.create(path, showWarnings = FALSE)) stop("Run directory already exists; refusing overwrite")
  p$path <- resolve_path(path); p$run_id <- id
  with_run_lock(p$path, function() {
    for (d in c("manifest", "logs", "checkpoints", "intermediate", "cache", "results"))
      if (!dir.create(file.path(p$path, d))) stop("Cannot initialize run directory")
    c <- schema_contract()
    ready <- Filter(function(r) isTRUE(r$integrity$exists) && !is.null(r$version) && !is.null(r$checksum), p$resources)
    biological <- list(schema_version = c$schema_version, id_schema_version = c$id_schema_version,
      analysis_id = p$analysis_id, package_version = as.character(utils::packageVersion("Neoantimon2")),
      fixture_kind = p$fixture_kind, analysis_contexts = p$analysis_context, predictions = list(),
      resources = lapply(ready, function(r) list(resource_id = if (is.null(r$resource_id)) r$logical_name else r$resource_id,
        release = r$version, checksum = paste0("sha256:", r$checksum$value))), allo_policy = NULL,
      output_schema_versions = list(events = c$schema_version, peptide_origins = c$schema_version,
        neoantigen_results = c$schema_version, multisample_summary = c$schema_version), uncertainty = list())
    git <- NULL
    repo <- p$config$paths$repository
    if (!is.null(repo) && dir.exists(file.path(repo, ".git")) && nzchar(Sys.which("git"))) {
      result <- suppressWarnings(system2("git", c("-C", shQuote(repo), "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE))
      if (length(result) == 1L && grepl("^[a-f0-9]{40,64}$", result)) git <- result
    }
    now <- utc_now()
    m <- list(manifest_version = infra_version(), checkpoint_version = infra_version(),
      config_contract_version = infra_version(), run_manifest = biological, run_id = id,
      project_id = p$project_id, created_at = now, updated_at = now, status = "INITIALIZED",
      config = unclass(p$config), config_fingerprint = p$config_fingerprint, input_fingerprint = p$input_fingerprint,
      samples = p$samples, hla = p$hla, input_records = p$inputs, resource_records = p$resources,
      checksum_policy = list(mode = p$config$checksum, declared = "sha256/file; verified only in sha256 mode",
        directory = "metadata is entry-only; recursive checksum unsupported"),
      platform = list(R = R.version.string, platform = R.version$platform, os = .Platform$OS.type),
      git_commit = git, checkpoints = list())
    write_manifest(m, p$path)
  })
  p
}

resume_project <- function(path, samples = NULL, hla = NULL, reference = NULL,
                           config = NULL, analysis_context = NULL, inputs = NULL) {
  initial <- read_manifest(path)
  initial_config <- if (is.null(config)) do.call(neo_config, c(native_allo_restore_pinned_config(initial$config), list(env = character()))) else config
  check_run_location(path, initial_config)
  with_run_lock(path, function() {
    m <- read_manifest(path)
    stored_config <- do.call(neo_config, c(native_allo_restore_pinned_config(m$config), list(env = character())))
    if (is.null(config)) config <- stored_config
    check_run_location(path, config)
    for (d in c("manifest", "logs", "checkpoints", "intermediate", "cache", "results")) {
      target <- resolve_path(file.path(path, d))
      if (!dir.exists(target) || !within_path(target, path)) stop("Unsafe/missing run subdirectory")
    }
    if (is.null(samples)) samples <- m$samples
    if (is.null(hla) && !is.null(m$hla)) hla <- lapply(m$hla, unlist)
    if (is.null(analysis_context)) analysis_context <- m$run_manifest$analysis_contexts
    if (is.null(inputs)) inputs <- lapply(m$input_records, function(r) list(input_id = r$input_id,
      input_role = r$input_role, sample_id = r$sample_id, path = r$declared_path, checksum = r$integrity$declared, version = r$version))
    p <- neo_project(samples = samples, hla = hla, reference = reference, config = config,
      project_id = m$project_id, analysis_context = analysis_context, inputs = inputs,
      fixture_kind = m$run_manifest$fixture_kind)
    if (!identical(infra_hash(p$config$allo), infra_hash(stored_config$allo)))
      stop("Resume allo configuration mismatch")
    if (!identical(infra_hash(p$config$indirect_hla), infra_hash(stored_config$indirect_hla)))
      stop("Resume indirect HLA configuration mismatch")
    if (!identical(infra_hash(p$config$allo_native_evidence), infra_hash(stored_config$allo_native_evidence)))
      stop("Resume native allo evidence configuration mismatch")
    if (!identical(p$analysis_id, m$run_manifest$analysis_id) ||
        !identical(p$config_fingerprint, m$config_fingerprint) ||
        !identical(p$input_fingerprint, m$input_fingerprint)) stop("Resume identity/fingerprint mismatch")
    # Even portable content declarations do not permit silent external path substitution.
    if (!identical(infra_hash(p$inputs), infra_hash(m$input_records)) ||
        !identical(infra_hash(p$resources), infra_hash(m$resource_records))) stop("Resume path/integrity mismatch")
    p$path <- path; p$run_id <- m$run_id
    p
  })
}

mutate_manifest <- function(project, mutate) {
  if (is.null(project$path)) stop("Run is not initialized")
  check_run_location(project$path, project$config)
  for (d in c("manifest", "checkpoints"))
    if (!within_path(resolve_path(file.path(project$path, d)), resolve_path(project$path))) stop("Unsafe metadata directory")
  with_run_lock(project$path, function() {
    m <- read_manifest(project$path)
    if (!identical(m$run_id, project$run_id) || !identical(m$run_manifest$analysis_id, project$analysis_id)) stop("Stale project identity")
    m <- mutate(m)
    m$updated_at <- utc_now()
    write_manifest(m, project$path)
    invisible(m)
  })
}
update_run_status <- function(project, status) {
  scalar_text(status, "status")
  transitions <- list(INITIALIZED = c("RUNNING", "FAILED", "INTERRUPTED"),
    RUNNING = c("COMPLETED", "FAILED", "INTERRUPTED"), FAILED = "RUNNING", INTERRUPTED = "RUNNING", COMPLETED = character())
  mutate_manifest(project, function(m) {
    if (!status %in% transitions[[m$status]]) stop("Invalid run status transition")
    if (status == "COMPLETED" && any(vapply(m$checkpoints, function(cp) cp$status != "COMPLETED", TRUE)))
      stop("Incomplete checkpoints")
    m$status <- status
    m
  })
}
