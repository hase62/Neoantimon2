somatic_priority_path <- function(project)
  file.path(project$path,"intermediate","somatic_priority.rds")

somatic_priority_fingerprints <- function(project, policy, runs) checkpoint_fingerprints(
  inputs=list(peptides=peptide_stage_output(project),
    raw=somatic_stage_output(somatic_raw_path(project),"intermediate/somatic_evidence_raw.rds"),
    eligibility=somatic_stage_output(somatic_eligibility_path(project),
      "intermediate/somatic_eligibility.rds"),
    prediction_runs=lapply(runs,function(x) list(run_id=x$run_id,
      checksum=digest::digest(file=predictor_result_path(project,x$run_id),
        algo="sha256",serialize=FALSE)))),
  config=policy,dependencies=list())

somatic_priority_files <- function() c("candidates.tsv","candidates.json","events.tsv",
  "evidence.tsv","predictions.tsv","run_summary.json","report.md")

somatic_flat <- function(rows, fields) {
  out <- lapply(fields,function(field) vapply(rows,function(row) {
    value <- row[[field]]
    if (is.null(value) || !length(value)) return(NA_character_)
    if (length(value)!=1L || is.list(value))
      return(as.character(jsonlite::toJSON(value,auto_unbox=TRUE,null="null",digits=NA)))
    as.character(value)
  },""))
  names(out) <- fields
  as.data.frame(out,check.names=FALSE,stringsAsFactors=FALSE)
}

somatic_write_tsv <- function(rows,fields,path) {
  utils::write.table(somatic_flat(rows,fields),path,sep="\t",row.names=FALSE,
    col.names=TRUE,quote=TRUE,na="")
}

somatic_priority_summary <- function(project,events,peptide,rows,runs) {
  grades <- c("A","B","C","D","E","INVALID","NOT_GRADED")
  counts <- stats::setNames(as.list(vapply(grades,function(g)
    sum(vapply(rows,function(r) identical(r$final_grade,g),TRUE)),0L)),grades)
  eligibility <- c("eligible","ineligible","not_evaluated")
  eligible_counts <- stats::setNames(as.list(vapply(eligibility,function(g)
    sum(vapply(rows,function(r) identical(r$technical_eligibility,g),TRUE)),0L)),eligibility)
  list(run_id=project$run_id,analysis_id=project$analysis_id,
    samples=as.list(sort(unique(vapply(project$samples,`[[`,"","sample_id")),method="radix")),
    events=length(events),peptide_origins=length(peptide$origins),
    candidates=length(rows),grade_counts=counts,eligibility_counts=eligible_counts,
    not_graded_reasons=as.list(sort(unique(as.character(unlist(lapply(Filter(function(r)
      r$final_grade=="NOT_GRADED",rows),`[[`,"grade_reason_codes")))),method="radix")),
    predictor_runs=lapply(runs,function(x) list(run_id=x$run_id,provider=x$provider,
      task=x$task,predictor_version=x$predictor_version,model_version=x$model_version,
      availability_status=x$capability$availability_status,
      evaluated_rows=length(x$evidence),status_rows=length(x$statuses),
      status_reasons=as.list(table(vapply(x$statuses,`[[`,"","reason"))))),
    rna_status_counts=as.list(table(vapply(rows,`[[`,"","rna_status"))),
    ccf_status_counts=as.list(table(vapply(rows,`[[`,"","ccf_status"))),
    cn_status_counts=as.list(table(vapply(somatic_priority_cn_records(rows),`[[`,"","status"))),
    unresolved_candidates=sum(vapply(rows,function(r) r$final_grade=="NOT_GRADED" ||
      r$ccf_status %in% c("ambiguous","not_evaluable","out_of_range") ||
      r$primary_prediction_status=="ambiguous",TRUE)),
    policy_version=somatic_priority_version())
}

somatic_priority_cn_records <- function(rows) lapply(rows,function(r)
  if (is.null(r$local_cn)) list(status="unavailable") else r$local_cn)

somatic_priority_materialize <- function(project,events,peptide,somatic,runs,rows,policy,
                                        root=file.path(project$path,"results")) {
  if (!dir.exists(root)) dir.create(root,recursive=TRUE,showWarnings=FALSE)
  fields <- c("candidate_id","sample_id","hla_class","mhc_key","ligand_sequence",
    "origin_ids","source_event_ids","priority_origin_id",
    "priority_origin_selection_policy","priority_source_event_ids",
    "final_grade","grade_status","grade_reason_codes",
    "technical_eligibility","primary_predictor","primary_metric","primary_rank",
    "pmhc_band","rna_status","rna_support","rna_vaf","transcript_tpm","gene_tpm",
    "ccf_status","ccf_value","tumor_vaf","tumor_alt_reads","exact_self_match",
    "self_match_status","predictor_disagreement","rank_within_sample_class","rank_status",
    "policy_version")
  somatic_write_tsv(rows,fields,file.path(root,"candidates.tsv"))
  atomic_json(rows,file.path(root,"candidates.json"))
  somatic_write_tsv(events,c("event_id","sample_id","contig","pos","end","event_role"),
    file.path(root,"events.tsv"))
  somatic_write_tsv(somatic$raw$evidence_rows,c("evidence_id","target_type","target_id",
    "evidence_domain","metric","value_number","value_text","unit","status","source"),
    file.path(root,"evidence.tsv"))
  predictions <- unlist(lapply(runs,function(run) c(run$evidence,
    lapply(run$statuses,function(s) c(list(run_id=run$run_id,
      predictor=run$provider),s)))),recursive=FALSE)
  somatic_write_tsv(predictions,c("run_id","candidate_id","predictor","task","metric",
    "percentile_rank","raw_score","predictor_version","model_version",
    "raw_output_ref","raw_row_index","status","reason"),file.path(root,"predictions.tsv"))
  summary <- somatic_priority_summary(project,events,peptide,rows,runs)
  atomic_json(summary,file.path(root,"run_summary.json"))
  top <- function(cl) {
    selected <- Filter(function(r) r$hla_class==cl && r$final_grade!="INVALID",rows)
    if (!length(selected)) return("- None")
    paste(vapply(utils::head(selected,10L),function(r) paste0("- ",r$candidate_id,
      " (",r$final_grade,", ",r$ligand_sequence,", rank ",
      if (is.null(r$rank_within_sample_class)) "unranked" else r$rank_within_sample_class,
      ")"),""),collapse="\n")
  }
  report <- c("# Somatic candidate results", "",
    paste0("Run: `",project$run_id,"`; analysis: `",project$analysis_id,"`."),"",
    paste0("Events: ",length(events),"; peptide origins: ",length(peptide$origins),
      "; candidates: ",length(rows),"."),"",
    paste0("Grades: ",paste(paste(names(summary$grade_counts),
      unlist(summary$grade_counts),sep="="),collapse=", "),"."),"",
    paste0("Eligibility: ",paste(paste(names(summary$eligibility_counts),
      unlist(summary$eligibility_counts),sep="="),collapse=", "),"."),"",
    paste0("Predictor runs: ",length(runs)," [",
      paste(vapply(summary$predictor_runs,function(x) paste0(x$provider,":",
        x$availability_status," (",x$evaluated_rows," rows)"),""),collapse=", "),
      "]; RNA states: ",
      paste(names(summary$rna_status_counts),unlist(summary$rna_status_counts),
        sep="=",collapse=", "),"."),"",
    paste0("CCF states: ",paste(names(summary$ccf_status_counts),
      unlist(summary$ccf_status_counts),sep="=",collapse=", "),
      "; CN states: ",paste(names(summary$cn_status_counts),
      unlist(summary$cn_status_counts),sep="=",collapse=", "),
      "; unresolved candidates: ",summary$unresolved_candidates,"."),"",
    paste0("NOT_GRADED reasons: ",paste(unlist(summary$not_graded_reasons),collapse=", "),"."),"",
    "## Class I", "",top("I"),"","## Class II","",top("II"),"",
    "A-E are evidence tiers under somatic-priority-v1, not calibrated clinical or immunogenicity probabilities.")
  writeLines(report,file.path(root,"report.md"),useBytes=TRUE)
}

somatic_priority_stage_run <- function(project,events,peptide,somatic) {
  policy <- project$config$somatic_priority
  runs <- somatic_prediction_runs(project)
  fp <- somatic_priority_fingerprints(project,policy,runs)
  cp <- Filter(function(x)x$step_id=="somatic_priority",
    read_manifest(project$path)$checkpoints)
  if (length(cp) && cp[[1L]]$status=="COMPLETED" &&
      checkpoint_reusable(project,"somatic_priority",fp))
    return(readRDS(somatic_priority_path(project)))
  if (!length(cp)) checkpoint_update(project,"somatic_priority",fingerprints=fp)
  checkpoint_update(project,"somatic_priority","RUNNING",fp)
  tryCatch({
    rows <- somatic_priority_construct(peptide$candidates,somatic,runs,policy)
    if (!setequal(vapply(rows,`[[`,"","candidate_id"),
      vapply(peptide$candidates,`[[`,"","candidate_id")))
      stop("Priority candidate identity mismatch")
    result <- list(schema_version=somatic_priority_schema_version(),
      policy=policy,run_id=project$run_id,rows=rows)
    tmp <- tempfile("priority-",tmpdir=dirname(somatic_priority_path(project)))
    on.exit(if (file.exists(tmp)) unlink(tmp),add=TRUE)
    saveRDS(result,tmp,version=3L)
    if (!file.rename(tmp,somatic_priority_path(project))) stop("Could not store priority result")
    somatic_priority_materialize(project,events,peptide,somatic,runs,rows,policy)
    files <- c("intermediate/somatic_priority.rds",paste0("results/",somatic_priority_files()))
    outputs <- lapply(files,function(relative) somatic_stage_output(
      file.path(project$path,relative),relative))
    checkpoint_update(project,"somatic_priority","COMPLETED",fp,output_records=outputs)
    result
  },error=function(e) {
    checkpoint_update(project,"somatic_priority","FAILED",fp,
      error=list(code="SOMATIC_PRIORITY_FAILED"))
    stop(e)
  })
}

somatic_priority_original_read <- function(project,snapshot=NULL) {
  if (is.null(project$config$somatic_evidence)) return(NULL)
  cp <- Filter(function(x)x$step_id=="somatic_priority",
    read_manifest(project$path)$checkpoints)
  if (length(cp)!=1L || cp[[1L]]$status!="COMPLETED")
    stop("Priority checkpoint missing or incomplete")
  verify_checkpoint_outputs(cp[[1L]],project$path)
  x <- readRDS(somatic_priority_path(project))
  if (!identical(x$schema_version,somatic_priority_schema_version()) ||
      !identical(x$run_id,project$run_id))
    stop("Priority run identity mismatch")
  if (!is.null(snapshot) && !setequal(vapply(x$rows,`[[`,"","candidate_id"),
    vapply(snapshot$neoantigen_results,`[[`,"","candidate_id")))
    stop("Priority snapshot candidate identity mismatch")
  x
}

somatic_priority_variant_dir <- function(project,fp) file.path(project$path,"results",
  "policies",infra_hash(fp))

somatic_priority_variant_verify <- function(project,fp,snapshot=NULL) {
  root <- somatic_priority_variant_dir(project,fp)
  checkpoint <- file.path(root,"checkpoint.json")
  if (!file.exists(checkpoint)) stop("Priority policy result is not materialized")
  meta <- jsonlite::read_json(checkpoint)
  if (!identical(meta$fingerprints,fp) || !identical(meta$policy_version,
    somatic_priority_version())) stop("Priority variant fingerprint mismatch")
  for (record in meta$outputs) {
    path <- file.path(root,record$file)
    if (!file.exists(path) || !identical(digest::digest(file=path,algo="sha256",
      serialize=FALSE),record$sha256)) stop("Priority variant output integrity failure")
  }
  result <- readRDS(file.path(root,"priority.rds"))
  if (!identical(result$policy,project$config$somatic_priority) ||
      !identical(result$run_id,project$run_id))
    stop("Priority variant policy identity mismatch")
  if (!is.null(snapshot) && !setequal(vapply(result$rows,`[[`,"","candidate_id"),
    vapply(snapshot$neoantigen_results,`[[`,"","candidate_id")))
    stop("Priority variant candidate identity mismatch")
  result
}

somatic_priority_postrun <- function(project,snapshot) {
  if (is.null(project$config$somatic_evidence)) return(invisible(NULL))
  original <- somatic_priority_original_read(project,snapshot)
  if (identical(original$policy,project$config$somatic_priority)) return(invisible(original))
  stored <- do.call(neo_config,c(read_manifest(project$path)$config,list(env=character())))
  if (!identical(infra_hash(stored$somatic_evidence),
    infra_hash(project$config$somatic_evidence)))
    stop("Completed run somatic evidence configuration changed")
  runs <- somatic_prediction_runs(project)
  fp <- somatic_priority_fingerprints(project,project$config$somatic_priority,runs)
  root <- somatic_priority_variant_dir(project,fp)
  if (file.exists(file.path(root,"checkpoint.json")))
    return(invisible(somatic_priority_variant_verify(project,fp,snapshot)))
  if (dir.exists(root)) stop("Incomplete priority variant directory")
  if (!dir.create(root,recursive=TRUE,showWarnings=FALSE))
    stop("Could not create priority variant directory")
  somatic <- somatic_read_sidecar(project,snapshot)
  peptide <- readRDS(peptide_stage_path(project))
  rows <- somatic_priority_construct(peptide$candidates,somatic,runs,
    project$config$somatic_priority)
  result <- list(schema_version=somatic_priority_schema_version(),
    policy=project$config$somatic_priority,run_id=project$run_id,rows=rows)
  saveRDS(result,file.path(root,"priority.rds"),version=3L)
  somatic_priority_materialize(project,snapshot$events,peptide,somatic,runs,rows,
    project$config$somatic_priority,root=root)
  files <- c("priority.rds",somatic_priority_files())
  outputs <- lapply(files,function(file) list(file=file,sha256=digest::digest(
    file=file.path(root,file),algo="sha256",serialize=FALSE)))
  atomic_json(list(policy_version=somatic_priority_version(),fingerprints=fp,
    outputs=outputs),file.path(root,"checkpoint.json"))
  invisible(somatic_priority_variant_verify(project,fp,snapshot))
}

somatic_priority_read <- function(project,snapshot=NULL,allow_original=FALSE) {
  original <- somatic_priority_original_read(project,snapshot)
  if (is.null(original) || allow_original ||
      identical(original$policy,project$config$somatic_priority)) return(original)
  runs <- somatic_prediction_runs(project)
  fp <- somatic_priority_fingerprints(project,project$config$somatic_priority,runs)
  somatic_priority_variant_verify(project,fp,snapshot)
}
