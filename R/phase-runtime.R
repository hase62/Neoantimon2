phase_stage_path <- function(project) file.path(project$path, "intermediate", "phase_resolution.json")

phase_stage_fingerprints <- function(project) checkpoint_fingerprints(
  inputs = list(input_fingerprint = project$input_fingerprint,
    event_output = event_stage_output(project)),
  config = list(config_fingerprint = project$config_fingerprint,
    interpretation = "t007-relative-phase-v1"),
  dependencies = list(inputs = sorted_records(identity_records(project$inputs))))

phase_stage_output <- function(project) list(path = "intermediate/phase_resolution.json",
  checksum = list(algorithm = "sha256", scope = "file",
    value = digest::digest(file = phase_stage_path(project), algo = "sha256", serialize = FALSE)))

phase_stage_run <- function(project, event_rows) {
  fp <- phase_stage_fingerprints(project)
  cp <- Filter(function(x) identical(x$step_id, "phase_resolution"),
    read_manifest(project$path)$checkpoints)
  if (length(cp) && cp[[1L]]$status == "COMPLETED" &&
      checkpoint_reusable(project, "phase_resolution", fp)) {
    rows <- jsonlite::read_json(phase_stage_path(project), simplifyVector = FALSE)
    validate_table(rows, "events")
    return(rows)
  }
  if (!length(cp)) checkpoint_update(project, "phase_resolution", fingerprints = fp)
  checkpoint_update(project, "phase_resolution", "RUNNING", fp)
  tryCatch({
    rows <- phase_resolve_events(event_rows)
    atomic_json(rows, phase_stage_path(project))
    checkpoint_update(project, "phase_resolution", "COMPLETED", fp,
      output_records = list(phase_stage_output(project)))
    rows
  }, error = function(e) {
    checkpoint_update(project, "phase_resolution", "FAILED", fp,
      error = list(code = runtime_error_code(e)))
    stop(e)
  })
}
