reconstruction_stage_path <- function(project)
  file.path(project$path, "intermediate", "sequence_reconstruction.rds")

reconstruction_stage_fingerprints <- function(project) checkpoint_fingerprints(
  inputs = list(phase_output = phase_stage_output(project)),
  config = list(policy = reconstruction_policy(), analysis_id = project$analysis_id),
  dependencies = list(reference = sorted_records(identity_records(project$resources))))

reconstruction_stage_output <- function(project) list(
  path = "intermediate/sequence_reconstruction.rds",
  checksum = list(algorithm = "sha256", scope = "file",
    value = digest::digest(file = reconstruction_stage_path(project), algo = "sha256",
                           serialize = FALSE)))

reconstruction_stage_run <- function(project, phase_rows) {
  fp <- reconstruction_stage_fingerprints(project)
  cp <- Filter(function(x) identical(x$step_id, "sequence_reconstruction"),
               read_manifest(project$path)$checkpoints)
  if (length(cp) && cp[[1L]]$status == "COMPLETED" &&
      checkpoint_reusable(project, "sequence_reconstruction", fp))
    return(readRDS(reconstruction_stage_path(project)))
  if (!length(cp)) checkpoint_update(project, "sequence_reconstruction", fingerprints = fp)
  checkpoint_update(project, "sequence_reconstruction", "RUNNING", fp)
  tryCatch({
    # A valid zero-event run does not need to reopen the large reference.
    small <- any(vapply(phase_rows, function(x)
      !is.null(x$contig) && x$event_type %in% c("SNV", "MNV", "INS", "DEL", "INDEL"), TRUE))
    contexts <- if (!small) list() else
      reconstruction_contexts(project, reference_bundle(project), phase_rows)
    path <- reconstruction_stage_path(project)
    tmp <- tempfile("sequence-reconstruction-", tmpdir = dirname(path))
    on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
    saveRDS(contexts, tmp, version = 3L)
    if (!file.rename(tmp, path)) stop("Could not atomically store reconstruction artifact")
    checkpoint_update(project, "sequence_reconstruction", "COMPLETED", fp,
                      output_records = list(reconstruction_stage_output(project)))
    contexts
  }, error = function(e) {
    checkpoint_update(project, "sequence_reconstruction", "FAILED", fp,
                      error = list(code = runtime_error_code(e)))
    stop(e)
  })
}
