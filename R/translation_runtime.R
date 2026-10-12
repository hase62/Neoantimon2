translation_stage_path <- function(project)
  file.path(project$path, "intermediate", "translation_resolution.rds")

translation_stage_fingerprints <- function(project) checkpoint_fingerprints(
  inputs = list(reconstruction_output = reconstruction_stage_output(project),
                phase_output = phase_stage_output(project)),
  config = list(policy = translation_policy(), analysis_id = project$analysis_id),
  dependencies = list(reference = sorted_records(identity_records(project$resources))))

translation_stage_output <- function(project) list(
  path = "intermediate/translation_resolution.rds",
  checksum = list(algorithm = "sha256", scope = "file",
    value = digest::digest(file = translation_stage_path(project), algo = "sha256",
                           serialize = FALSE)))

translation_stage_run <- function(project, phase_rows, contexts) {
  fp <- translation_stage_fingerprints(project)
  cp <- Filter(function(x) identical(x$step_id, "translation_resolution"),
               read_manifest(project$path)$checkpoints)
  if (length(cp) && cp[[1L]]$status == "COMPLETED" &&
      checkpoint_reusable(project, "translation_resolution", fp))
    return(readRDS(translation_stage_path(project)))
  if (!length(cp)) checkpoint_update(project, "translation_resolution", fingerprints = fp)
  checkpoint_update(project, "translation_resolution", "RUNNING", fp)
  tryCatch({
    resolved <- if (!length(contexts)) list() else
      translation_contexts(project, reference_bundle(project), phase_rows, contexts)
    path <- translation_stage_path(project)
    tmp <- tempfile("translation-resolution-", tmpdir = dirname(path))
    on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
    saveRDS(resolved, tmp, version = 3L)
    if (!file.rename(tmp, path)) stop("Could not atomically store translation artifact")
    checkpoint_update(project, "translation_resolution", "COMPLETED", fp,
                      output_records = list(translation_stage_output(project)))
    resolved
  }, error = function(e) {
    checkpoint_update(project, "translation_resolution", "FAILED", fp,
                      error = list(code = runtime_error_code(e)))
    stop(e)
  })
}
