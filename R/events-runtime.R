event_stage_path <- function(project) file.path(project$path, "intermediate", "event_import.json")

event_stage_output <- function(project) list(path = "intermediate/event_import.json",
  checksum = list(algorithm = "sha256", scope = "file",
    value = digest::digest(file = event_stage_path(project), algo = "sha256", serialize = FALSE)))

event_stage_fingerprints <- function(project) checkpoint_fingerprints(
  inputs = list(input_fingerprint = project$input_fingerprint),
  config = list(config_fingerprint = project$config_fingerprint,
                assembly = project$config$scientific$genome_assembly,
                normalization = "t006-primary-small-variant-v1"),
  dependencies = list(reference = sorted_records(identity_records(project$resources))))

event_stage_run <- function(project) {
  fp <- event_stage_fingerprints(project)
  cp <- Filter(function(x) identical(x$step_id, "event_import"),
               read_manifest(project$path)$checkpoints)
  if (length(cp) && cp[[1L]]$status == "COMPLETED" &&
      checkpoint_reusable(project, "event_import", fp)) {
    rows <- jsonlite::read_json(event_stage_path(project), simplifyVector = FALSE)
    validate_table(rows, "events")
    return(rows)
  }
  if (!length(cp)) checkpoint_update(project, "event_import", fingerprints = fp)
  checkpoint_update(project, "event_import", "RUNNING", fp)
  tryCatch({
    rows <- event_import_project(project)
    validate_table(rows, "events")
    atomic_json(rows, event_stage_path(project))
    checkpoint_update(project, "event_import", "COMPLETED", fp,
                      output_records = list(event_stage_output(project)))
    rows
  }, error = function(e) {
    checkpoint_update(project, "event_import", "FAILED", fp,
                      error = list(code = runtime_error_code(e)))
    stop(e)
  })
}
