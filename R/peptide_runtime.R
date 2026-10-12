peptide_stage_path <- function(project)
  file.path(project$path, "intermediate", "peptide_candidates.rds")

peptide_stage_output <- function(project) list(path = "intermediate/peptide_candidates.rds",
  checksum = list(algorithm = "sha256", scope = "file",
    value = digest::digest(file = peptide_stage_path(project), algo = "sha256", serialize = FALSE)))

peptide_reference_checksum <- function(project) {
  found <- Filter(function(x) x$logical_name == "reference.protein_fasta", project$resources)
  if (!length(found) || is.null(found[[1L]]$resolved_path) ||
      !file.exists(found[[1L]]$resolved_path)) return(NULL)
  list(algorithm = "sha256", value = digest::digest(file = found[[1L]]$resolved_path,
    algo = "sha256", serialize = FALSE))
}

peptide_stage_fingerprints <- function(project) checkpoint_fingerprints(
  inputs = list(translation_output = translation_stage_output(project)),
  config = list(implementation = "peptide-candidate-v1", policy = peptide_policy(project$config$peptides),
    contract = peptide_contract()$peptide_version, hla = project$hla),
  dependencies = list(reference = sorted_records(identity_records(project$resources)),
    self_proteome_content = peptide_reference_checksum(project)))

peptide_stage_run <- function(project, contexts, events) {
  # The complete T009 artifact checksum covers proteins, mapping and provenance.
  fp <- peptide_stage_fingerprints(project)
  cp <- Filter(function(x) x$step_id == "peptide_candidates", read_manifest(project$path)$checkpoints)
  if (length(cp) && cp[[1L]]$status == "COMPLETED" &&
      checkpoint_reusable(project, "peptide_candidates", fp)) return(readRDS(peptide_stage_path(project)))
  if (!length(cp)) checkpoint_update(project, "peptide_candidates", fingerprints = fp)
  checkpoint_update(project, "peptide_candidates", "RUNNING", fp)
  tryCatch({
    resources <- Filter(function(x) x$logical_name == "reference.protein_fasta", project$resources)
    resource <- if (length(resources)) resources[[1L]] else NULL
    if (!is.null(resource) && !identical(resource$integrity,
        integrity_record(resource$resolved_path, resource$integrity$mode, resource$integrity$declared)))
      stop("Reference proteome changed since project binding")
    # T009 already validated the bundle. Its bound catalog sequence and Sec
    # annotation establish the literal-U representation without a new GTF parse.
    sec_specs <- Filter(function(x) length(x$translation$spec$sec_aa_positions) > 0L, contexts)
    literal_u <- length(sec_specs) > 0L && all(vapply(sec_specs, function(x) {
      s <- x$translation$spec
      !is.null(s$reference_protein_catalog) && all(substring(s$reference_protein_catalog,
        s$sec_aa_positions, s$sec_aa_positions) == "U")
    }, TRUE))
    result <- peptide_construct(project, contexts, events, resource$resolved_path,
      identity_records(resources), if (literal_u) "LITERAL_U" else "UNDECLARED")
    result$translation_artifact <- translation_stage_output(project)
    result$self_proteome_checksum <- peptide_reference_checksum(project)
    target <- peptide_stage_path(project)
    tmp <- tempfile("peptide-", tmpdir = dirname(target))
    on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
    saveRDS(result, tmp, version = 3L)
    if (!file.rename(tmp, target)) stop("Could not atomically store peptide output")
    checkpoint_update(project, "peptide_candidates", "COMPLETED", fp,
      output_records = list(peptide_stage_output(project)))
    result
  }, error = function(e) {
    checkpoint_update(project, "peptide_candidates", "FAILED", fp,
      error = list(code = runtime_error_code(e)))
    stop(e)
  })
}
