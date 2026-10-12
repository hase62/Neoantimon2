predictor_model_identity <- function(path) {
  if (is.null(path) || !dir.exists(path)) return(NULL)
  files <- sort(list.files(path, recursive = TRUE, full.names = TRUE,
    all.files = TRUE, no.. = TRUE), method = "radix")
  files <- files[!dir.exists(files)]
  list(path = normalizePath(path, mustWork = TRUE), files = lapply(files, function(f)
    list(relative_path = substring(f,nchar(path)+2L), sha256 = digest::digest(
      file=f,algo="sha256",serialize=FALSE))))
}

predictor_run_identity <- function(prepared, options = list()) {
  predictor_options_validate(prepared$provider,prepared$task,options)
  cap <- prepared$capability
  inputs <- lapply(prepared$candidates, function(x) list(candidate_id=x$candidate_id,
    peptide_id=x$sample_peptide_id, sequence=x$ligand_sequence, presenter=x$mhc_key,
    analysis_context=x$analysis_context, sample_id=x$sample_id))
  paste0("predictor-run-", infra_hash(list(candidates=inputs, provider=prepared$provider,
    requested_task=prepared$task, predictor_version=cap$predictor_version,
    model_version=cap$model_version, model_identity=cap$model_identity,
    executable_checksum=cap$executable_checksum, availability=cap$availability_status,
    adapter_version=predictor_adapter_version(),
    evidence_schema_version=predictor_evidence_schema_version(), options=options)))
}

predictor_options_validate <- function(provider, task, options) {
  if (!is.list(options)) stop("predictor options must be a record")
  strict_keys(options,c("affinity_only","timeout_seconds","retry_token"),"predictor options")
  if (!is.null(options$affinity_only)) {
    if (!is.logical(options$affinity_only) || length(options$affinity_only)!=1L ||
        is.na(options$affinity_only)) stop("Invalid affinity_only")
    if (provider!="mhcflurry") stop("affinity_only is supported only for MHCflurry")
    if (isTRUE(options$affinity_only) && task!="mhc1_binding")
      stop("affinity_only requires an MHCflurry binding task")
  }
  if (!is.null(options$timeout_seconds) && (!is.numeric(options$timeout_seconds) ||
      length(options$timeout_seconds)!=1L || !is.finite(options$timeout_seconds) ||
      options$timeout_seconds < 1)) stop("Invalid timeout_seconds")
  if (!is.null(options$retry_token)) scalar_text(options$retry_token,"retry_token")
  invisible(TRUE)
}

predictor_policy <- function(x) {
  if (is.null(x)) return(NULL)
  if (!is.list(x) || !is.null(names(x))) stop("predictors must be a record array")
  for (i in seq_along(x)) {
    p <- x[[i]]
    strict_keys(p,c("provider","task","executable","model_path","options"),"predictor request")
    if (!p$provider %in% names(predictor_registry())) stop("Unknown predictor")
    if (!p$task %in% peptide_prediction_tasks()) stop("Unknown predictor task")
    if (p$provider == "fixture") stop("Fixture predictor cannot be configured for production")
    for (k in c("executable","model_path")) if (!is.null(p[[k]])) scalar_text(p[[k]],k)
    if (is.null(p$options)) p$options <- list()
    predictor_options_validate(p$provider,p$task,p$options)
    x[[i]] <- p
  }
  x
}

predictor_artifact_root <- function(project) file.path(project$path,"intermediate","prediction_runs")

predictor_result_path <- function(project, run_id)
  file.path(predictor_artifact_root(project),run_id,"result.rds")

predictor_output_record <- function(project, run_id) {
  path <- predictor_result_path(project,run_id)
  list(path=file.path("intermediate","prediction_runs",run_id,"result.rds"),
    checksum=list(algorithm="sha256",scope="file",value=digest::digest(
      file=path,algo="sha256",serialize=FALSE)))
}

predictor_fingerprints <- function(project, prepared, options) checkpoint_fingerprints(
  inputs=list(peptide_artifact=peptide_stage_output(project),
    candidates=lapply(prepared$candidates,function(x) x[c("candidate_id","ligand_sequence","mhc_key","sample_peptide_id")])),
  config=list(provider=prepared$provider, task=prepared$task, options=options,
    adapter_version=predictor_adapter_version(),
    evidence_schema_version=predictor_evidence_schema_version()),
  dependencies=list(capability=prepared$capability[c("availability_status","predictor_version",
    "model_version","model_identity","executable_checksum")]))

predictor_read_run <- function(artifact_root, run_id) {
  target <- file.path(artifact_root,run_id)
  stored <- file.path(target,"result.rds")
  checksum_path <- file.path(target,"result.sha256")
  if (!file.exists(stored) || !file.exists(checksum_path) ||
      !identical(readLines(checksum_path,warn=FALSE),
        digest::digest(file=stored,algo="sha256",serialize=FALSE)))
    stop("Predictor result cache integrity failure")
  old <- readRDS(stored)
  if (!identical(old$run_id,run_id)) stop("Predictor cache identity mismatch")
  predictor_validate_result(old)
  for (record in old$raw_artifacts) {
    path <- file.path(target,record$file)
    if (!file.exists(path) || !identical(digest::digest(file=path,
        algo="sha256",serialize=FALSE),record$sha256))
      stop("Predictor raw cache integrity failure")
  }
  for (command in old$commands) {
    path <- file.path(target,command$input_file)
    if (!file.exists(path) || !identical(digest::digest(file=path,
        algo="sha256",serialize=FALSE),command$input_sha256))
      stop("Predictor input cache integrity failure")
  }
  old
}

predictor_run <- function(prepared, artifact_root, options = list(), fixture_raw = NULL) {
  predictor_options_validate(prepared$provider,prepared$task,options)
  run_id <- predictor_run_identity(prepared,options)
  target <- file.path(artifact_root,run_id)
  stored <- file.path(target,"result.rds")
  if (file.exists(stored)) return(predictor_read_run(artifact_root,run_id))
  dir.create(artifact_root,recursive=TRUE,showWarnings=FALSE)
  temp <- tempfile("predictor-",tmpdir=artifact_root)
  dir.create(temp)
  on.exit(if (dir.exists(temp)) unlink(temp,recursive=TRUE),add=TRUE)
  raw <- list(); raw_refs <- list(); artifacts <- list(); commands <- list()
  for (allele in names(prepared$batches)) {
    tag <- digest::digest(allele,algo="sha256",serialize=FALSE)
    input_name <- paste0("input-",tag,if (prepared$provider=="mhcflurry") ".csv" else ".txt")
    raw_name <- paste0("raw-",tag,if (prepared$provider=="mhcflurry") ".csv" else ".txt")
    input_path <- file.path(temp,input_name); output_path <- file.path(temp,raw_name)
    writeLines(predictor_input(prepared,allele),input_path,useBytes=TRUE)
    command <- predictor_command(prepared,allele,input_path,output_path)
    # Only relative artifact names and a digest are stored; temp paths are
    # excluded from deterministic run identity and persisted command evidence.
    commands[[allele]] <- list(executable=command$executable,
      arguments=if (prepared$provider=="fixture") character() else
        if (prepared$provider=="mhcflurry") c("INPUT.csv","--out","OUTPUT.csv",
        if (isTRUE(options$affinity_only)) "--affinity-only" else character()) else
        c(if (prepared$provider=="netmhcpan") "-p" else "-f", "INPUT.txt",
          if (prepared$provider=="netmhciipan") c("-inptype","1") else character(),
          "-a",allele,"-BA"),
      input_file=input_name,input_sha256=digest::digest(file=input_path,algo="sha256",serialize=FALSE))
    output <- NULL
    if (!is.null(fixture_raw)) output <- fixture_raw[[allele]]
    else if (prepared$provider=="fixture") output <- NULL
    else if (prepared$capability$available) {
      args <- command$args
      if (prepared$provider=="mhcflurry") {
        if (isTRUE(options$affinity_only)) args <- c(args,"--affinity-only")
        if (!is.null(prepared$capability$model_identity))
          args <- c(args,"--models",prepared$capability$model_identity$path)
      }
      attempt <- tryCatch({
        stdout <- suppressWarnings(system2(command$executable,args,stdout=TRUE,stderr=TRUE,
          timeout=if (is.null(options$timeout_seconds)) 600 else options$timeout_seconds))
        list(log=paste(stdout,collapse="\n"),
          exit_status=if (is.null(attr(stdout,"status"))) 0L else as.integer(attr(stdout,"status")))
      },error=function(e) list(log=conditionMessage(e),exit_status=1L))
      log_name <- paste0("execution-",tag,".txt")
      writeLines(attempt$log,file.path(temp,log_name),useBytes=TRUE)
      artifacts[[length(artifacts)+1L]] <- list(allele=allele,file=log_name,kind="execution_log",
        sha256=digest::digest(file=file.path(temp,log_name),algo="sha256",serialize=FALSE))
      commands[[allele]]$exit_status <- attempt$exit_status
      if (attempt$exit_status!=0L && file.exists(output_path))
        artifacts[[length(artifacts)+1L]] <- list(allele=allele,file=raw_name,
          kind="partial_raw_output",sha256=digest::digest(file=output_path,
            algo="sha256",serialize=FALSE))
      if (attempt$exit_status==0L) output <- if (prepared$provider=="mhcflurry") {
        if (file.exists(output_path)) paste(readLines(output_path,warn=FALSE),collapse="\n") else NULL
      } else attempt$log
    }
    if (!is.null(output)) {
      writeLines(output,output_path,useBytes=TRUE)
      raw[[allele]] <- output
      artifacts[[length(artifacts)+1L]] <- list(allele=allele,file=raw_name,kind="raw_output",
        sha256=digest::digest(file=output_path,algo="sha256",serialize=FALSE))
      raw_refs[[allele]] <- file.path("intermediate","prediction_runs",run_id,raw_name)
    }
  }
  normalized <- predictor_normalize(prepared,raw,run_id,raw_refs,
    prepared$capability$model_version)
  result <- c(list(run_id=run_id, provider=prepared$provider, task=prepared$task,
    predictor_version=prepared$capability$predictor_version,
    model_version=prepared$capability$model_version,
    capability=prepared$capability, candidate_ids=as.list(vapply(prepared$candidates,`[[`,"","candidate_id")),
    commands=commands, raw_artifacts=artifacts, adapter_version=predictor_adapter_version(),
    evidence_schema_version=predictor_evidence_schema_version(),
    options=options),normalized)
  predictor_validate_result(result)
  saveRDS(result,file.path(temp,"result.rds"),version=3L)
  writeLines(digest::digest(file=file.path(temp,"result.rds"),algo="sha256",
    serialize=FALSE),file.path(temp,"result.sha256"))
  if (!file.rename(temp,target)) stop("Could not atomically store predictor run")
  result
}

prediction_stage_run <- function(project, peptide_result) {
  policy <- predictor_policy(project$config$predictors)
  if (is.null(policy) || !length(policy)) return(list())
  out <- list()
  for (request in policy) {
    cap <- predictor_capability(request$provider,request$executable,request$model_path)
    prepared <- predictor_prepare(peptide_result$candidates,request$provider,request$task,cap)
    options <- request$options
    run_id <- predictor_run_identity(prepared,options)
    step <- paste0("prediction_",sub("^predictor-run-","",run_id))
    fp <- predictor_fingerprints(project,prepared,options)
    cp <- Filter(function(x) x$step_id==step,read_manifest(project$path)$checkpoints)
    if (length(cp) && cp[[1L]]$status=="COMPLETED" &&
        checkpoint_reusable(project,step,fp)) {
      out[[run_id]] <- predictor_read_run(predictor_artifact_root(project),run_id); next
    }
    if (!length(cp)) checkpoint_update(project,step,fingerprints=fp)
    checkpoint_update(project,step,"RUNNING",fp)
    tryCatch({
      result <- predictor_run(prepared,predictor_artifact_root(project),options)
      checkpoint_update(project,step,"COMPLETED",fp,
        output_records=list(predictor_output_record(project,run_id)))
      out[[run_id]] <- result
    },error=function(e) {
      checkpoint_update(project,step,"FAILED",fp,error=list(code="PREDICTOR_ADAPTER_FAILED"))
      stop(e)
    })
  }
  out
}

# Completed biological runs are immutable. A later installation or model
# update may create a new, separately fingerprinted prediction sidecar without
# changing the frozen peptide artifact or reopening the run manifest.
prediction_postrun <- function(project, requests) {
  if (!inherits(project,"neo_project") || is.null(project$path))
    stop("Initialized project required")
  read_runtime_snapshot(project)
  policy <- predictor_policy(requests)
  peptide <- readRDS(peptide_stage_path(project))
  out <- list()
  for (request in policy) {
    cap <- predictor_capability(request$provider,request$executable,request$model_path)
    prepared <- predictor_prepare(peptide$candidates,request$provider,request$task,cap)
    result <- predictor_run(prepared,predictor_artifact_root(project),request$options)
    out[[result$run_id]] <- result
  }
  out
}
