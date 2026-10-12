# Temporary T004 bootstrap transport, not a public/frozen result-file contract.
runtime_snapshot_path <- function(project) {
  if (!inherits(project, "neo_project")) stop("project must be a neo_project", call. = FALSE)
  if (is.null(project$path)) return(NULL)
  file.path(project$path, "intermediate", "t004_runtime_snapshot.json")
}

empty_runtime_snapshot <- function(project) {
  manifest <- read_manifest(project$path)$run_manifest
  list(run_manifest = manifest, events = list(), peptide_origins = list(),
       neoantigen_results = list(), multisample_summary = list())
}

default_runtime_adapter <- function(project) {
  snapshot <- empty_runtime_snapshot(project)
  if (longitudinal_aggregate(project)) return(snapshot)
  if (project$fixture_kind == "real") {
    imported <- event_stage_run(project)
    snapshot$events <- phase_stage_run(project, imported)
    contexts <- reconstruction_stage_run(project, snapshot$events)
    translated <- translation_stage_run(project, snapshot$events, contexts)
    peptide <- peptide_stage_run(project, translated, snapshot$events)
    extension <- source_stage_run(project)
    long_read <- long_read_stage_run(project)
    prediction_input <- peptide
    if (!is.null(extension)) {
      prediction_input$candidates <- c(peptide$candidates, extension$candidates)
      if (length(prediction_input$candidates))
        prediction_input$candidates <- prediction_input$candidates[
          !duplicated(vapply(prediction_input$candidates, `[[`, "", "candidate_id"))]
    }
    if (!is.null(long_read)) {
      prediction_input$candidates <- c(prediction_input$candidates, long_read$candidates)
      if (length(prediction_input$candidates))
        prediction_input$candidates <- prediction_input$candidates[
          !duplicated(vapply(prediction_input$candidates, `[[`, "", "candidate_id"))]
    }
    prediction_stage_run(project, prediction_input)
    somatic <- somatic_stage_run(project, snapshot$events, peptide)
    if (!is.null(somatic)) somatic_priority_stage_run(project, snapshot$events, peptide, somatic)
    mutate_manifest(project, function(m) {
      m$run_manifest$output_schema_versions <- peptide_contract()$output_schema_versions
      m
    })
    snapshot$run_manifest <- read_manifest(project$path)$run_manifest
    snapshot$peptide_origins <- peptide$origins
    snapshot$neoantigen_results <- peptide$candidates
  }
  else if (length(project$inputs)) stop("Fixture inputs require an injected internal adapter", call. = FALSE)
  snapshot
}

runtime_fingerprints <- function(project) checkpoint_fingerprints(
  inputs = list(input_fingerprint = project$input_fingerprint),
  config = list(config_fingerprint = project$config_fingerprint,
                schema_version = project$schema_version,
                peptide_policy = if (project$fixture_kind == "real")
                  peptide_policy(project$config$peptides) else NULL,
                somatic_evidence = if (project$fixture_kind == "real")
                  project$config$somatic_evidence else NULL,
                somatic_priority = if (project$fixture_kind == "real")
                  project$config$somatic_priority else NULL),
  dependencies = list())

runtime_output_record <- function(project) {
  path <- runtime_snapshot_path(project)
  list(path = "intermediate/t004_runtime_snapshot.json",
       checksum = list(algorithm = "sha256", scope = "file",
                       value = digest::digest(file = path, algo = "sha256", serialize = FALSE)))
}

validate_runtime_manifest <- function(snapshot, project) {
  manifest <- read_manifest(project$path)
  if (!identical(manifest$run_id, project$run_id) ||
      !identical(manifest$run_manifest$analysis_id, project$analysis_id))
    stop("Run output analysis identity mismatch", call. = FALSE)
  if (!identical(snapshot$run_manifest, manifest$run_manifest))
    stop("Runtime snapshot run_manifest differs from authoritative manifest", call. = FALSE)
  invisible(TRUE)
}

read_runtime_snapshot <- function(project) {
  path <- runtime_snapshot_path(project)
  if (is.null(path) || !file.exists(path))
    stop("Results are not available: runtime output is missing", call. = FALSE)
  if (!within_path(resolve_path(path), file.path(project$path, "intermediate")))
    stop("Unsafe run output path", call. = FALSE)
  snapshot <- jsonlite::read_json(path, simplifyVector = FALSE)
  validate_snapshot(snapshot, runtime_contract(snapshot))
  validate_runtime_manifest(snapshot, project)
  cp <- Filter(function(x) identical(x$step_id, "pipeline"),
    read_manifest(project$path)$checkpoints)
  if (!length(cp) || cp[[1L]]$status != "COMPLETED" ||
      !identical(cp[[1L]]$output_records, list(runtime_output_record(project))))
    stop("Runtime snapshot checkpoint integrity failure", call. = FALSE)
  verify_checkpoint_outputs(cp[[1L]], project$path)
  event_cp <- Filter(function(x) identical(x$step_id, "event_import"),
    read_manifest(project$path)$checkpoints)
  if (length(event_cp)) {
    if (length(event_cp) != 1L || event_cp[[1L]]$status != "COMPLETED")
      stop("Event import checkpoint integrity failure", call. = FALSE)
    verify_checkpoint_outputs(event_cp[[1L]], project$path)
  }
  phase_cp <- Filter(function(x) identical(x$step_id, "phase_resolution"),
    read_manifest(project$path)$checkpoints)
  if (length(phase_cp)) {
    if (length(phase_cp) != 1L || phase_cp[[1L]]$status != "COMPLETED")
      stop("Phase resolution checkpoint integrity failure", call. = FALSE)
    verify_checkpoint_outputs(phase_cp[[1L]], project$path)
  }
  reconstruction_cp <- Filter(function(x) identical(x$step_id, "sequence_reconstruction"),
    read_manifest(project$path)$checkpoints)
  if (project$fixture_kind == "real" && length(event_cp) && !length(reconstruction_cp))
    stop("Sequence reconstruction checkpoint missing from completed real run", call. = FALSE)
  if (length(reconstruction_cp)) {
    if (length(reconstruction_cp) != 1L || reconstruction_cp[[1L]]$status != "COMPLETED")
      stop("Sequence reconstruction checkpoint integrity failure", call. = FALSE)
    verify_checkpoint_outputs(reconstruction_cp[[1L]], project$path)
  }
  translation_cp <- Filter(function(x) identical(x$step_id, "translation_resolution"),
    read_manifest(project$path)$checkpoints)
  if (project$fixture_kind == "real" && length(event_cp) && !length(translation_cp))
    stop("Translation resolution checkpoint missing from completed real run", call. = FALSE)
  if (length(translation_cp)) {
    if (length(translation_cp) != 1L || translation_cp[[1L]]$status != "COMPLETED")
      stop("Translation resolution checkpoint integrity failure", call. = FALSE)
    verify_checkpoint_outputs(translation_cp[[1L]], project$path)
  }
  peptide_cp <- Filter(function(x) x$step_id == "peptide_candidates",
    read_manifest(project$path)$checkpoints)
  if (project$fixture_kind == "real" && !is.null(runtime_contract(snapshot)$peptide_version) &&
      (!length(peptide_cp) || peptide_cp[[1L]]$status != "COMPLETED"))
    stop("Peptide candidate checkpoint integrity failure", call. = FALSE)
  if (length(peptide_cp)) verify_checkpoint_outputs(peptide_cp[[1L]], project$path)
  if (length(source_inputs(project))) source_read(project)
  if (length(long_read_inputs(project))) long_read_read(project)
  prediction_cps <- Filter(function(x) startsWith(x$step_id, "prediction_"),
    read_manifest(project$path)$checkpoints)
  for (prediction_cp in prediction_cps) {
    if (prediction_cp$status != "COMPLETED")
      stop("Prediction checkpoint integrity failure", call. = FALSE)
    verify_checkpoint_outputs(prediction_cp, project$path)
  }
  if (project$fixture_kind=="real" && !is.null(project$config$somatic_evidence))
    somatic_read_sidecar(project,snapshot)
  if (project$fixture_kind=="real" && !is.null(project$config$somatic_evidence))
    somatic_priority_read(project,snapshot,allow_original=TRUE)
  snapshot
}

runtime_error_code <- function(condition) {
  code <- attr(condition, "neo_error_code", exact = TRUE)
  if (is.character(code) && length(code) == 1L && grepl("^[A-Z][A-Z0-9_]*$", code)) code
  else "PIPELINE_FAILED"
}

neo_run <- function(project) run_runtime(project)

# Adapter injection is internal and restricted to synthetic/mock fixtures.
run_runtime <- function(project, adapter = NULL) {
  if (!inherits(project, "neo_project")) stop("project must be a neo_project", call. = FALSE)
  if (!is.null(adapter)) {
    if (!is.function(adapter)) stop("adapter must be a function", call. = FALSE)
    if (!project$fixture_kind %in% c("synthetic", "mock"))
      stop("Injected adapters are restricted to synthetic/mock projects", call. = FALSE)
  } else adapter <- default_runtime_adapter
  if (is.null(project$path)) project <- initialize_run(project)

  manifest <- read_manifest(project$path)
  if (manifest$status == "COMPLETED") {
    snapshot <- read_runtime_snapshot(project)
    somatic_priority_postrun(project,snapshot)
    hla_loh_postrun(project,snapshot)
    recognition_postrun(project,snapshot)
    presentation_postrun(project,snapshot)
    longitudinal_postrun(project,snapshot)
    vaccine_postrun(project)
    neo_export_bundle(project)
    return(project)
  }
  update_run_status(project, "RUNNING")
  result <- tryCatch({
    fp <- runtime_fingerprints(project)
    cp <- Filter(function(x) identical(x$step_id, "pipeline"),
      read_manifest(project$path)$checkpoints)
    reusable <- length(cp) && cp[[1L]]$status == "COMPLETED" &&
      checkpoint_reusable(project, "pipeline", fp)
    if (!reusable) {
      if (!length(cp)) checkpoint_update(project, "pipeline", fingerprints = fp)
      checkpoint_update(project, "pipeline", "RUNNING", fp)
      snapshot <- adapter(project)
      if (!is.list(snapshot) || !identical(sort(names(snapshot)),
        sort(c("run_manifest", "events", "peptide_origins", "neoantigen_results",
               "multisample_summary"))))
        stop("Adapter must return one complete runtime snapshot")
      validate_snapshot(snapshot, runtime_contract(snapshot))
      validate_runtime_manifest(snapshot, project)
      atomic_json(snapshot, runtime_snapshot_path(project))
      checkpoint_update(project, "pipeline", "COMPLETED", fp,
        output_records = list(runtime_output_record(project)))
    }
    read_runtime_snapshot(project)
    update_run_status(project, "COMPLETED")
    project
  }, error = function(e) {
    # Mark execution failed first, even if checkpoint cleanup itself fails.
    update_run_status(project, "FAILED")
    cp <- Filter(function(x) identical(x$step_id, "pipeline"),
      read_manifest(project$path)$checkpoints)
    if (length(cp) && cp[[1L]]$status == "RUNNING")
      checkpoint_update(project, "pipeline", "FAILED",
        cp[[1L]][c("input_fingerprint", "config_fingerprint", "dependency_fingerprint")],
        error = list(code = runtime_error_code(e)))
    stop(conditionMessage(e), call. = FALSE)
  })
  hla_loh_postrun(result,read_runtime_snapshot(result))
  recognition_postrun(result,read_runtime_snapshot(result))
  presentation_postrun(result,read_runtime_snapshot(result))
  longitudinal_postrun(result,read_runtime_snapshot(result))
  vaccine_postrun(result)
  neo_export_bundle(result)
  result
}
