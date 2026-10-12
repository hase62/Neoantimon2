somatic_priority_version <- function() "somatic-priority-v1"
somatic_priority_schema_version <- function() "1.0.0"
somatic_origin_selection_version <- function() "best-supported-eligible-origin-v1"

somatic_priority_policy <- function(x = NULL) {
  if (is.null(x)) x <- list()
  strict_keys(x, c("policy_version", "schema_version", "thresholds", "primary",
    "rna_min_alt_reads", "ranking_fields", "origin_selection_policy"),
    "somatic priority")
  if (!is.null(x$policy_version) && !identical(x$policy_version,somatic_priority_version()))
    stop("Unsupported somatic priority policy version")
  if (!is.null(x$schema_version) && !identical(x$schema_version,somatic_priority_schema_version()))
    stop("Unsupported somatic priority schema version")
  if (!is.null(x$origin_selection_policy) &&
      !identical(x$origin_selection_policy,somatic_origin_selection_version()))
    stop("Unsupported origin selection policy")
  defaults <- list(I=list(strong=0.5, binder=2), II=list(strong=1, binder=5))
  if (!is.null(x$thresholds)) {
    strict_keys(x$thresholds, c("I", "II"), "priority thresholds")
    for (cl in names(x$thresholds)) {
      strict_keys(x$thresholds[[cl]], c("strong", "binder"), "class thresholds")
      for (key in names(x$thresholds[[cl]])) defaults[[cl]][[key]] <- x$thresholds[[cl]][[key]]
    }
  }
  for (cl in names(defaults)) {
    z <- defaults[[cl]]
    if (any(!vapply(z, function(v) is.numeric(v) && length(v)==1L &&
      is.finite(v) && v>=0 && v<=100, TRUE)) || z$strong>z$binder)
      stop("Invalid priority rank thresholds")
  }
  primary <- list(I=list(provider="netmhcpan", task="mhc1_presentation",
    metric="EL_rank", run_id=NULL, model_version=NULL),
    II=list(provider="netmhciipan", task="mhc2_presentation",
      metric="EL_rank", run_id=NULL, model_version=NULL))
  if (!is.null(x$primary)) {
    strict_keys(x$primary,c("I","II"),"primary prediction")
    for (cl in names(x$primary)) {
      strict_keys(x$primary[[cl]],names(primary[[cl]]),"primary class prediction")
      for (key in names(x$primary[[cl]])) primary[[cl]][key] <-
        list(x$primary[[cl]][[key]])
    }
  }
  for (cl in names(primary)) {
    p <- primary[[cl]]
    for (key in c("provider","task","metric")) scalar_text(p[[key]],key)
    if (!p$provider %in% names(predictor_registry())) stop("Unknown primary provider")
    if (!p$task %in% peptide_prediction_tasks()) stop("Unknown primary task")
    if (!p$metric %in% c("EL_rank","BA_rank","presentation_percentile"))
      stop("Primary metric must be a normalized percentile rank")
    for (key in c("run_id","model_version")) if (!is.null(p[[key]])) scalar_text(p[[key]],key)
  }
  rna <- if (is.null(x$rna_min_alt_reads)) 1L else x$rna_min_alt_reads
  if (!is.numeric(rna) || length(rna)!=1L || !is.finite(rna) ||
      rna!=floor(rna) || rna<1) stop("Invalid RNA support threshold")
  fields <- c("final_grade","primary_rank","rna_status","transcript_tpm",
    "gene_tpm","ccf_value","tumor_alt_reads","exact_self_match","candidate_id")
  ranking <- if (is.null(x$ranking_fields)) fields else unlist(x$ranking_fields,
    use.names=FALSE)
  if (!identical(sort(ranking),sort(fields)))
    stop("Ranking fields must contain each named component exactly once")
  list(policy_version=somatic_priority_version(),schema_version=somatic_priority_schema_version(),
    thresholds=defaults,primary=primary,rna_min_alt_reads=as.integer(rna),
    ranking_fields=as.list(ranking),origin_selection_policy=somatic_origin_selection_version())
}

somatic_pmhc_band <- function(rank, hla_class, policy) {
  z <- policy$thresholds[[hla_class]]
  if (is.null(z)) stop("Unknown HLA class")
  if (is.null(rank) || length(rank)!=1L || !is.finite(rank) || rank<0 || rank>100)
    return(list(band=NULL,comparisons=list(),status="unavailable"))
  comparisons <- lapply(c("strong","binder"), function(name)
    list(metric="percentile_rank",observed=rank,threshold=z[[name]],
      operator="<=",satisfied=rank<=z[[name]],policy_version=policy$policy_version))
  names(comparisons) <- c("strong","binder")
  list(band=if (comparisons$strong$satisfied) "STRONG" else
    if (comparisons$binder$satisfied) "BINDER" else "OUTSIDE",
    comparisons=comparisons,status="evaluated")
}

somatic_prediction_runs <- function(project) {
  root <- predictor_artifact_root(project)
  if (!dir.exists(root)) return(list())
  ids <- sort(list.dirs(root, recursive=FALSE, full.names=FALSE),method="radix")
  out <- lapply(ids,function(id) predictor_read_run(root,id))
  names(out) <- ids
  out
}

somatic_select_primary <- function(candidate, runs, policy) {
  p <- policy$primary[[candidate$hla_class]]
  matching <- Filter(function(run) identical(run$provider,p$provider) &&
    (is.null(p$run_id) || identical(run$run_id,p$run_id)) &&
    (is.null(p$model_version) || identical(run$model_version,p$model_version)),runs)
  rows <- unlist(lapply(matching,function(run) Filter(function(e)
    identical(e$candidate_id,candidate$candidate_id) && identical(e$task,p$task) &&
      identical(e$metric,p$metric) && identical(e$status,"evaluated"),run$evidence)),
    recursive=FALSE)
  statuses <- unlist(lapply(matching,function(run) Filter(function(s)
    (is.null(s$candidate_id) || identical(s$candidate_id,candidate$candidate_id)) &&
      (is.null(s$task) || identical(s$task,p$task)),run$statuses)),recursive=FALSE)
  duplicate_output <- any(vapply(statuses,function(s)
    identical(s$reason,"duplicate_output_row"),TRUE))
  state <- if (length(rows)>1L || duplicate_output) "ambiguous" else if (
    length(rows)==1L && !is.null(rows[[1L]]$percentile_rank)) "evaluated" else if (
      length(rows)==1L) "unavailable" else if (length(statuses)) {
      if (any(vapply(statuses,function(s) identical(s$status,"unsupported"),TRUE)))
        "unsupported" else if (any(vapply(statuses,function(s) identical(s$status,"failed"),TRUE)))
        "failed" else "unavailable"
    } else "unavailable"
  list(status=state, selected=if (state=="evaluated") rows[[1L]] else NULL,
    applicable_rows=rows, statuses=statuses, selection=p,
    reason=if (state=="ambiguous") "PRIMARY_PREDICTION_AMBIGUOUS" else
      if (state!="evaluated") paste0("PRIMARY_PREDICTION_",toupper(state)) else NULL)
}

somatic_secondary_predictions <- function(candidate,runs,primary,policy) {
  all <- unlist(lapply(runs,function(run) Filter(function(e)
    identical(e$candidate_id,candidate$candidate_id) &&
      e$metric %in% c("EL_rank","BA_rank","presentation_percentile") &&
      !is.null(e$percentile_rank),run$evidence)),recursive=FALSE)
  chosen <- primary$selected
  secondary <- Filter(function(e) is.null(chosen) ||
    !identical(c(e$run_id,e$task,e$metric,e$raw_row_index),
      c(chosen$run_id,chosen$task,chosen$metric,chosen$raw_row_index)),all)
  bands <- vapply(all,function(e) somatic_pmhc_band(e$percentile_rank,
    candidate$hla_class,policy)$band,"")
  providers <- vapply(all,`[[`,"","predictor")
  pairs <- if (length(all)>=2L) utils::combn(seq_along(all),2L) else
    matrix(integer(),nrow=2L)
  list(rows=secondary,evaluable_providers=length(unique(vapply(all,`[[`,"","predictor"))),
    disagreement=if (ncol(pairs)) any(apply(pairs,2L,function(i)
      providers[[i[[1L]]]]!=providers[[i[[2L]]]] &&
        bands[[i[[1L]]]]!=bands[[i[[2L]]]])) else FALSE)
}

somatic_priority_scalar <- function(x) if (is.null(x) || !length(x)) NULL else x[[1L]]

somatic_origin_priority_summary <- function(path,technical_status,event_index,
                                            expression_index,policy) {
  ids <- unlist(path$somatic_event_ids,use.names=FALSE)
  events <- event_index[ids]
  if (any(vapply(events,is.null,TRUE))) stop("Priority origin event evidence missing")
  event_summaries <- lapply(events,function(e) list(event_id=e$event_id,
    evidence_id=e$evidence_id,rna_status=e$rna$status,
    rna_alt_reads=e$rna$alt_reads,
    rna_support=identical(e$rna$status,"observed_support") &&
      !is.null(e$rna$alt_reads) && e$rna$alt_reads>=policy$rna_min_alt_reads,
    tumor_status=e$tumor$status,tumor_alt_reads=e$tumor$alt_reads))
  support <- vapply(event_summaries,`[[`,TRUE,"rna_support")
  rna_states <- vapply(event_summaries,`[[`,"","rna_status")
  rna_status <- if (!length(events)) "not_evaluated" else if (all(support))
    "observed_support" else if (any(support)) "partial_support" else if (
      all(rna_states %in% c("observed_support","observed_no_alt")))
      "observed_no_alt" else if ("not_evaluated" %in% rna_states)
        "not_evaluated" else "unavailable"
  expression <- expression_index[unlist(path$expression_evidence_ids,use.names=FALSE)]
  best <- function(field) {
    values <- unlist(lapply(expression,function(e)e[[field]]),use.names=FALSE)
    if (length(values)) max(values) else NULL
  }
  ccf <- if (length(events)==1L &&
    events[[1L]]$ccf$status %in% c("conditional","external"))
      events[[1L]]$ccf$bounded_value else NULL
  tumor_counts <- lapply(event_summaries,`[[`,"tumor_alt_reads")
  tumor_evaluable <- length(events)>0L && all(vapply(event_summaries,function(e)
    !is.null(e$tumor_alt_reads) && e$tumor_status %in%
      c("observed_support","observed_no_alt"),TRUE))
  list(origin_id=path$origin_id,technical_eligibility=technical_status,
    source_event_ids=as.list(ids),event_evidence_ids=path$event_evidence_ids,
    expression_evidence_ids=path$expression_evidence_ids,
    event_summaries=event_summaries,rna_status=rna_status,
    rna_support=length(events)>0L && all(support),
    rna_alt_reads=lapply(events,function(e)e$rna$alt_reads),
    rna_vaf=if (length(events)==1L) events[[1L]]$rna$vaf else NULL,
    transcript_tpm=best("transcript_tpm"),gene_tpm=best("gene_tpm"),
    tumor_vaf=if (length(events)==1L) events[[1L]]$tumor$vaf else NULL,
    tumor_alt_reads=if (tumor_evaluable) min(unlist(tumor_counts)) else NULL,
    tumor_alt_reads_summary="minimum_required_event_tumor_alt_reads",
    ccf_status=if (length(events)==1L) events[[1L]]$ccf$status else
      if (length(events)>1L) "multiple_events" else "not_evaluated",
    ccf_value=ccf,ccf_solutions=if (length(events)==1L) events[[1L]]$ccf$solutions else list(),
    purity=if (length(events)==1L) events[[1L]]$purity$value else NULL,
    local_cn=if (length(events)==1L) events[[1L]]$cn else NULL,
    generic_loh=if (length(events)==1L) events[[1L]]$loh else NULL)
}

somatic_select_priority_origin <- function(summaries) {
  eligible <- Filter(function(x) identical(x$technical_eligibility,"eligible"),summaries)
  if (!length(eligible)) return(NULL)
  rna_order <- c(observed_support=1L,partial_support=2L,unavailable=3L,
    not_evaluated=3L,observed_no_alt=4L)
  descending <- function(field) vapply(eligible,function(x)
    if (is.null(x[[field]])) Inf else -x[[field]],0)
  ix <- order(vapply(eligible,function(x)rna_order[[x$rna_status]],0L),
    descending("transcript_tpm"),descending("gene_tpm"),
    descending("ccf_value"),descending("tumor_alt_reads"),
    vapply(eligible,`[[`,"","origin_id"),method="radix")
  eligible[[ix[[1L]]]]
}

somatic_candidate_modifiers <- function(candidate,link,eligibility,raw,policy) {
  event_index <- stats::setNames(raw$event_evidence,vapply(raw$event_evidence,`[[`,"","event_id"))
  expression_index <- stats::setNames(raw$expression_evidence,
    vapply(raw$expression_evidence,`[[`,"","evidence_id"))
  eligibility_index <- stats::setNames(eligibility$origin_results,
    vapply(eligibility$origin_results,`[[`,"","origin_id"))
  summaries <- lapply(link$origin_paths,function(path) {
    technical <- eligibility_index[[path$origin_id]]
    if (is.null(technical)) stop("Priority origin eligibility missing")
    somatic_origin_priority_summary(path,technical$status,event_index,
      expression_index,policy)
  })
  summaries <- summaries[order(vapply(summaries,`[[`,"","origin_id"),method="radix")]
  selected <- somatic_select_priority_origin(summaries)
  field <- function(name,default=NULL) if (is.null(selected)) default else selected[[name]]
  list(origin_id=field("origin_id"),priority_origin_id=field("origin_id"),
    priority_origin_selection_policy=policy$origin_selection_policy,
    priority_source_event_ids=field("source_event_ids",list()),
    origin_priority_summaries=summaries,
    source_event_ids=as.list(sort(unique(as.character(unlist(lapply(link$origin_paths,
      `[[`,"somatic_event_ids")))),method="radix")),
    event_evidence_ids=as.list(sort(unique(as.character(unlist(lapply(link$origin_paths,
      `[[`,"event_evidence_ids")))),method="radix")),
    expression_evidence_ids=field("expression_evidence_ids",list()),
    rna_status=field("rna_status","unavailable"),rna_support=field("rna_support",FALSE),
    rna_alt_reads=field("rna_alt_reads",list()),rna_vaf=field("rna_vaf"),
    transcript_tpm=field("transcript_tpm"),gene_tpm=field("gene_tpm"),
    tumor_vaf=field("tumor_vaf"),tumor_alt_reads=field("tumor_alt_reads"),
    tumor_alt_reads_summary="minimum_required_event_tumor_alt_reads",
    ccf_status=field("ccf_status","not_evaluated"),ccf_value=field("ccf_value"),
    ccf_solutions=field("ccf_solutions",list()),purity=field("purity"),
    local_cn=field("local_cn"),generic_loh=field("generic_loh"),
    exact_self_match=candidate$self_exact_match,self_match_status=candidate$self_status,
    sec_status=if (grepl("U",candidate$ligand_sequence,fixed=TRUE)) "contains_U" else "none")
}

somatic_grade_candidate <- function(candidate,eligibility,primary,band,modifiers,policy) {
  if (identical(eligibility$eligibility_status,"ineligible"))
    return(list(grade="INVALID",status="INVALID",reasons=list("TECHNICAL_INELIGIBLE")))
  if (!identical(eligibility$eligibility_status,"eligible"))
    return(list(grade="NOT_GRADED",status="NOT_GRADED",
      reasons=list("TECHNICAL_NOT_EVALUATED")))
  if (primary$status!="evaluated" || is.null(band$band))
    return(list(grade="NOT_GRADED",status="NOT_GRADED",reasons=list(primary$reason)))
  grade <- if (band$band=="STRONG") {
    if (modifiers$rna_support) "A" else "B"
  } else if (band$band=="BINDER") {
    if (modifiers$rna_support) "C" else "D"
  } else "E"
  list(grade=grade,status="GRADED",reasons=list(paste0("PMHC_",band$band),
    paste0("RNA_",toupper(modifiers$rna_status))))
}

somatic_priority_construct <- function(candidates,somatic,runs,policy) {
  raw <- somatic$raw
  links <- stats::setNames(raw$candidate_links,vapply(raw$candidate_links,`[[`,"","candidate_id"))
  elig <- stats::setNames(somatic$eligibility$candidate_eligibility,
    vapply(somatic$eligibility$candidate_eligibility,`[[`,"","candidate_id"))
  rows <- lapply(candidates,function(candidate) {
    id <- candidate$candidate_id
    primary <- somatic_select_primary(candidate,runs,policy)
    rank <- if (is.null(primary$selected)) NULL else primary$selected$percentile_rank
    band <- somatic_pmhc_band(rank,candidate$hla_class,policy)
    if (length(band$comparisons)) for (name in names(band$comparisons)) {
      band$comparisons[[name]]$metric <- primary$selected$metric
      band$comparisons[[name]]$predictor <- primary$selected$predictor
      band$comparisons[[name]]$predictor_version <- primary$selected$predictor_version
      band$comparisons[[name]]$model_version <- primary$selected$model_version
      band$comparisons[[name]]$run_id <- primary$selected$run_id
    }
    modifiers <- somatic_candidate_modifiers(candidate,links[[id]],elig[[id]],raw,policy)
    grade <- somatic_grade_candidate(candidate,elig[[id]],primary,band,modifiers,policy)
    secondary <- somatic_secondary_predictions(candidate,runs,primary,policy)
    result <- candidate
    result[names(modifiers)] <- modifiers
    updates <- list(final_grade=grade$grade,grade_status=grade$status,
      grade_reason_codes=grade$reasons,technical_eligibility=elig[[id]]$eligibility_status,
      technical_reason_records=elig[[id]]$reason_records,
      primary_prediction_status=primary$status,primary_selection=primary$selection,
      primary_prediction=primary$selected,primary_applicable_rows=primary$applicable_rows,
      primary_status_records=primary$statuses,primary_predictor=primary$selection$provider,
      primary_metric=primary$selection$metric,primary_rank=rank,pmhc_band=band$band,
      pmhc_comparisons=band$comparisons,secondary_predictions=secondary$rows,
      evaluable_provider_count=secondary$evaluable_providers,
      predictor_disagreement=secondary$disagreement,
      policy_version=policy$policy_version,rank_within_sample_class=NULL,
      rank_status=if (grade$grade=="INVALID") "INVALID" else "RANKED")
    result[names(updates)] <- updates
    result$grading_status <- grade$status
    result$grade_before_warnings <- grade$grade
    result$prediction_status <- if (primary$status=="evaluated") "SUCCESS" else "NOT_RUN"
    result$not_evaluated_codes <- if (grade$grade=="NOT_GRADED") grade$reasons else list()
    result
  })
  somatic_priority_rank(rows,policy)
}

somatic_priority_rank <- function(rows,policy) {
  if (!length(rows)) return(rows)
  grade_order <- c(A=1L,B=2L,C=3L,D=4L,E=5L,NOT_GRADED=6L,INVALID=7L)
  rna_order <- c(observed_support=1L,partial_support=2L,unavailable=3L,
    not_evaluated=3L,observed_no_alt=4L)
  component <- function(row,field) {
    v <- row[[field]]
    if (field=="final_grade") return(grade_order[[v]])
    if (field=="rna_status") return(if (v %in% names(rna_order)) rna_order[[v]] else 2L)
    if (field=="candidate_id") return(v)
    if (field=="exact_self_match") return(if (is.null(v)) 2L else if (isTRUE(v)) 1L else 0L)
    if (is.null(v) || length(v)!=1L || is.na(v)) return(Inf)
    if (field %in% c("transcript_tpm","gene_tpm","ccf_value","tumor_alt_reads"))
      return(-v)
    v
  }
  scopes <- unique(vapply(rows,function(r)paste(r$sample_id,r$hla_class,sep="\037"),""))
  for (scope in scopes) {
    indexes <- which(vapply(rows,function(r)
      identical(paste(r$sample_id,r$hla_class,sep="\037"),scope),TRUE))
    sortable <- Filter(function(i) rows[[i]]$final_grade!="INVALID",as.list(indexes))
    if (!length(sortable)) next
    ix <- unlist(sortable)
    columns <- lapply(unlist(policy$ranking_fields),function(field)
      vapply(rows[ix],component,if (field=="candidate_id") "" else 0,field=field))
    order_ix <- do.call(order,c(columns,list(method="radix")))
    for (j in seq_along(order_ix)) rows[[ix[[order_ix[[j]]]]]]$rank_within_sample_class <-
      as.integer(j)
  }
  rows[order(vapply(rows,`[[`,"","sample_id"),
    vapply(rows,`[[`,"","hla_class"),
    vapply(rows,function(r)if (is.null(r$rank_within_sample_class)) Inf else
      r$rank_within_sample_class,0),vapply(rows,`[[`,"","candidate_id"),method="radix")]
}
