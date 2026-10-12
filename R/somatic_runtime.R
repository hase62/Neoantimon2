somatic_raw_path <- function(project)
  file.path(project$path,"intermediate","somatic_evidence_raw.rds")

somatic_eligibility_path <- function(project)
  file.path(project$path,"intermediate","somatic_eligibility.rds")

somatic_stage_output <- function(path,relative) list(path=relative,
  checksum=list(algorithm="sha256",scope="file",value=digest::digest(
    file=path,algo="sha256",serialize=FALSE)))

somatic_raw_fingerprints <- function(project,sources) checkpoint_fingerprints(
  inputs=list(events=phase_stage_output(project),peptides=peptide_stage_output(project),
    sources=lapply(sources,function(s) list(kind=s$kind,sample_id=s$sample_id,
      source=s$source,version=s$source_version,reference_build=s$reference_build,
      sha256=s$sha256))),
  config=list(schema_version=somatic_evidence_version(),
    ccf_model_version=somatic_ccf_model_version(),facets_adapter="facets-cncf-v1"),
  dependencies=list())

somatic_eligibility_fingerprints <- function(project,policy) checkpoint_fingerprints(
  inputs=list(raw=somatic_stage_output(somatic_raw_path(project),
    "intermediate/somatic_evidence_raw.rds")),
  config=list(policy=policy,policy_version=somatic_eligibility_version()),
  dependencies=list())

somatic_stage_run <- function(project,events,peptide) {
  policy <- project$config$somatic_evidence
  if (is.null(policy)) return(NULL)
  sources <- somatic_source_files(project)
  raw_fp <- somatic_raw_fingerprints(project,sources)
  raw_cp <- Filter(function(x) x$step_id=="somatic_evidence_raw",
    read_manifest(project$path)$checkpoints)
  if (length(raw_cp) && raw_cp[[1L]]$status=="COMPLETED" &&
      checkpoint_reusable(project,"somatic_evidence_raw",raw_fp)) {
    raw <- readRDS(somatic_raw_path(project))
  } else {
    if (!length(raw_cp)) checkpoint_update(project,"somatic_evidence_raw",fingerprints=raw_fp)
    checkpoint_update(project,"somatic_evidence_raw","RUNNING",raw_fp)
    tryCatch({
      raw <- somatic_raw_construct(events,peptide,somatic_source_rows(sources),project$run_id)
      somatic_validate_raw(raw,events,peptide)
      tmp <- tempfile("somatic-raw-",tmpdir=dirname(somatic_raw_path(project)))
      on.exit(if (file.exists(tmp)) unlink(tmp),add=TRUE)
      saveRDS(raw,tmp,version=3L)
      if (!file.rename(tmp,somatic_raw_path(project))) stop("Could not store somatic raw evidence")
      checkpoint_update(project,"somatic_evidence_raw","COMPLETED",raw_fp,
        output_records=list(somatic_stage_output(somatic_raw_path(project),
          "intermediate/somatic_evidence_raw.rds")))
    },error=function(e) {
      checkpoint_update(project,"somatic_evidence_raw","FAILED",raw_fp,
        error=list(code="SOMATIC_EVIDENCE_FAILED"))
      stop(e)
    })
  }
  somatic_validate_raw(raw,events,peptide)
  eligibility_fp <- somatic_eligibility_fingerprints(project,policy$eligibility)
  eligibility_cp <- Filter(function(x) x$step_id=="somatic_eligibility",
    read_manifest(project$path)$checkpoints)
  if (length(eligibility_cp) && eligibility_cp[[1L]]$status=="COMPLETED" &&
      checkpoint_reusable(project,"somatic_eligibility",eligibility_fp)) {
    eligibility <- readRDS(somatic_eligibility_path(project))
  } else {
    if (!length(eligibility_cp)) checkpoint_update(project,"somatic_eligibility",
      fingerprints=eligibility_fp)
    checkpoint_update(project,"somatic_eligibility","RUNNING",eligibility_fp)
    tryCatch({
      eligibility <- somatic_eligibility_construct(raw,policy$eligibility)
      somatic_validate_eligibility(eligibility,raw)
      tmp <- tempfile("somatic-eligibility-",tmpdir=dirname(somatic_eligibility_path(project)))
      on.exit(if (file.exists(tmp)) unlink(tmp),add=TRUE)
      saveRDS(eligibility,tmp,version=3L)
      if (!file.rename(tmp,somatic_eligibility_path(project)))
        stop("Could not store somatic eligibility")
      checkpoint_update(project,"somatic_eligibility","COMPLETED",eligibility_fp,
        output_records=list(somatic_stage_output(somatic_eligibility_path(project),
          "intermediate/somatic_eligibility.rds")))
    },error=function(e) {
      checkpoint_update(project,"somatic_eligibility","FAILED",eligibility_fp,
        error=list(code="SOMATIC_ELIGIBILITY_FAILED"))
      stop(e)
    })
  }
  somatic_validate_eligibility(eligibility,raw)
  list(raw=raw,eligibility=eligibility)
}

somatic_read_sidecar <- function(project,snapshot=NULL) {
  if (is.null(project$config$somatic_evidence)) return(NULL)
  manifest <- read_manifest(project$path)
  for (step in c("somatic_evidence_raw","somatic_eligibility")) {
    cp <- Filter(function(x) x$step_id==step,manifest$checkpoints)
    if (length(cp)!=1L || cp[[1L]]$status!="COMPLETED")
      stop("Somatic evidence checkpoint missing or incomplete")
    verify_checkpoint_outputs(cp[[1L]],project$path)
  }
  raw <- readRDS(somatic_raw_path(project))
  eligibility <- readRDS(somatic_eligibility_path(project))
  if (!identical(raw$schema_version,somatic_evidence_version()) ||
      !identical(raw$run_id,project$run_id)) stop("Somatic evidence identity mismatch")
  if (!is.null(snapshot) && (!setequal(
      vapply(raw$event_evidence,`[[`,"","event_id"),
      vapply(snapshot$events,`[[`,"","event_id")) ||
      !setequal(vapply(raw$candidate_links,`[[`,"","candidate_id"),
        vapply(snapshot$neoantigen_results,`[[`,"","candidate_id"))))
    stop("Somatic evidence snapshot identity mismatch")
  somatic_validate_eligibility(eligibility,raw)
  list(raw=raw,eligibility=eligibility)
}
