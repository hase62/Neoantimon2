# Explicit opt-in paired BA queries using an already configured native provider.
# This version-pinned example uses the package's tested internal adapters; it is
# not an additional public API. neo_priorities() itself remains read-only.
library(Neoantimon2)
run_paper_binding_queries <- function(project, output, request=1L,
    metric=c("BA_nM","BA_rank"), execute=FALSE) {
  metric <- match.arg(metric)
  if (!is.logical(execute) || length(execute)!=1L || is.na(execute)) stop("execute must be TRUE or FALSE")
  internal <- function(name) getFromNamespace(name,"Neoantimon2")
  graph <- internal("external_graph")(project)
  configured <- project$config$predictors
  if (!length(configured) || length(request)!=1L || !is.numeric(request) ||
      !is.finite(request) || request!=floor(request) || request<1 || request>length(configured))
    stop("Select an existing configured native predictor request; no software is installed automatically")
  cfg <- configured[[request]]
  if (!cfg$provider %in% c("netmhcpan","netmhciipan","mhcflurry")) stop("Unsupported paired BA provider")
  task <- if(startsWith(cfg$task,"mhc1_")) "mhc1_binding" else "mhc2_binding"
  class <- if(task=="mhc1_binding") "I" else "II"
  cap <- internal("predictor_capability")(cfg$provider,cfg$executable,cfg$model_path)
  queries <- list();mapping <- list();excluded <- list()
  for (c in graph$candidates) {
    if (c$hla_class!=class) next
    for (oid in unlist(c$origin_ids)) {
      origins <- Filter(function(o)o$origin_id==oid,graph$origins)
      if (length(origins)!=1L) stop("Ambiguous Origin identity")
      wt <- origins[[1L]]$wt_peptide_sequence
      if (is.null(wt) || is.na(wt) || !nzchar(wt) || nchar(wt)!=nchar(c$peptide_sequence)) {
        excluded[[length(excluded)+1L]] <- list(candidate_id=c$candidate_id,origin_id=oid,reason="missing_equal_length_WT")
        next
      }
      for (kind in c("mutant","wild_type")) {
        sequence <- if(kind=="mutant") c$peptide_sequence else wt
        id <- paste0("paper-query:",digest::digest(list(sequence,c$mhc_key,c$candidate_id,oid,kind),algo="sha256"))
        q <- c;q$candidate_id <- id;q$sample_peptide_id <- id
        q$ligand_sequence <- sequence;q$peptide_sequence <- sequence;q$origin_ids <- list(oid)
        queries[[id]] <- q
        mapping[[length(mapping)+1L]] <- list(kind=kind,query_id=id,candidate_id=c$candidate_id,
          origin_id=oid,sample_id=c$sample_id,peptide_sequence=sequence,mhc_key=c$mhc_key)
      }
    }
  }
  if (!length(queries)) stop("No equal-length mutant/WT candidate pairs for this provider class")
  prepared <- internal("predictor_prepare")(unname(queries),cfg$provider,task,cap)
  plan <- list(provider=cfg$provider,task=task,metric=metric,query_pairs=length(mapping)/2L,unique_peptide_hla_queries=length(prepared$queries),
    available=isTRUE(cap$available),unavailable_queries=prepared$statuses,
    excluded_origins=excluded,execution_authorized=execute,
    scope="new independent BA sidecar; original native predictions, Grades and ranks unchanged")
  if (!execute) return(plan) # Inspection/version detection only; no predictions or output writes.
  if (!plan$available || length(prepared$statuses))
    stop("Paired BA prerequisites unavailable; inspect execute=FALSE plan. No predictions were invented")
  output <- internal("resolve_path")(output)
  if (file.exists(output)) stop("Use a new independent output directory")
  roots <- c(project$path,project$config$paths$data_root,project$config$paths$repository)
  for (root in roots) if (!is.null(root) &&
      (internal("within_path")(output,internal("resolve_path")(root)) ||
       internal("within_path")(internal("resolve_path")(root),output)))
    stop("Paired BA output must not overlap the completed project, data or repository")
  ancestor <- dirname(output)
  repeat {
    if (file.exists(file.path(ancestor,".git"))) stop("Do not write prediction output inside a repository")
    parent <- dirname(ancestor);if(parent==ancestor) break;ancestor <- parent
  }
  dir.create(output,recursive=TRUE)
  if (!dir.exists(output)) stop("Cannot create paired BA output")
  jsonlite::write_json(plan,file.path(output,"plan.json"),auto_unbox=TRUE,pretty=TRUE,null="null")
  saveRDS(mapping,file.path(output,"query-origin-mapping.rds"))
  result <- internal("predictor_run")(prepared,file.path(output,"intermediate","prediction_runs"),cfg$options)
  evidence <- Filter(function(e)e$status=="evaluated" && e$task==task && e$metric==metric,result$evidence)
  ids <- vapply(evidence,`[[`,"","candidate_id")
  if (length(result$statuses) || anyDuplicated(ids) || !setequal(ids,names(queries)))
    stop("Provider did not return one unambiguous requested BA metric per query; raw sidecar retained")
  if (any(vapply(evidence,function(e)e$raw_score<0,TRUE)))
    stop("Invalid negative BA measurement; raw sidecar retained, no complete TSV pair")
  rows <- lapply(evidence,function(e) {
    q <- queries[[e$candidate_id]]
    list(schema_version="1.0.0",peptide_sequence=q$ligand_sequence,scope="peptide_hla",mhc_key=q$mhc_key,
      task=e$task,provider=e$predictor,predictor_version=if(is.null(e$predictor_version)) "" else e$predictor_version,
      model_version=if(is.null(e$model_version)) "" else e$model_version,run_id=e$run_id,metric=e$metric,
      metric_kind=if(e$unit=="percent") "percentile_rank" else "raw_score",
      metric_semantics=if(e$unit=="percent") "mhc_binding_percentile" else "raw",
      score=format(e$raw_score,digits=17,scientific=FALSE,trim=TRUE),score_direction=e$score_direction,
      unit=e$unit,status="evaluated",reason="",source_row=paste0(e$run_id,":",e$raw_output_ref,":",e$raw_row_index))
  })
  for (kind in c("mutant","wild_type")) {
    selected <- unique(vapply(Filter(function(x)x$kind==kind,mapping),`[[`,"","query_id"))
    # The existing adapter evaluates each distinct peptide/HLA once and fans
    # out to source-specific queries. Deduplicate only identical measurements;
    # every source relationship remains in query-origin-mapping.rds.
    table <- unique(do.call(rbind,lapply(rows[ids %in% selected],as.data.frame,stringsAsFactors=FALSE)))
    table <- internal("external_validate_predictions")(table)
    utils::write.table(table,file.path(output,paste0(kind,".tsv")),sep="\t",quote=TRUE,row.names=FALSE,na="")
  }
  paths <- c(mutant=file.path(output,"mutant.tsv"),wild_type=file.path(output,"wild_type.tsv"))
  jsonlite::write_json(list(status="PAIRED_BA_OUTPUT_COMPLETE",plan=plan,run_id=result$run_id,
    output_sha256=lapply(c(paths,query_origin_mapping=file.path(output,"query-origin-mapping.rds"),
      plan=file.path(output,"plan.json")),function(p)digest::digest(file=p,algo="sha256",serialize=FALSE)),
    capability=result$capability),file.path(output,"completion.json"),auto_unbox=TRUE,pretty=TRUE,null="null")
  list(plan=plan,run_id=result$run_id,mutant=unname(paths[["mutant"]]),
    wild_type=unname(paths[["wild_type"]]),output=output)
}
