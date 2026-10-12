neo_explain_core <- function(project, candidate_id = NULL, reason = NULL,
                        peptide = NULL, sample_id = NULL,
                        format = c("structured", "text")) {
  snapshot <- read_runtime_snapshot(project)
  if (!isTRUE(project$longitudinal_accessor_guard) && longitudinal_aggregate(project))
    return(longitudinal_aggregate_explain(project, candidate_id, reason, peptide, sample_id, match.arg(format)))
  somatic <- somatic_read_sidecar(project)
  priority <- somatic_priority_read(project,snapshot)
  hla <- hla_loh_read(project)
  core_evidence <- if (is.null(priority)) list() else
    unlist(lapply(somatic_prediction_runs(project),`[[`,"core_evidence"),recursive=FALSE)
  format <- match.arg(format)
  for (x in c("candidate_id", "reason", "peptide", "sample_id"))
    validate_filter(get(x), x)
  if (all(vapply(list(candidate_id, reason, peptide, sample_id), is.null, logical(1))))
    stop("At least one explanation query is required", call. = FALSE)
  results <- filter_records(if (is.null(priority)) snapshot$neoantigen_results else priority$rows,
    list(candidate_id = candidate_id, peptide_sequence = peptide, sample_id = sample_id), reason)
  candidates <- lapply(results, function(result) {
    origins <- Filter(function(x) x$origin_id %in% unlist(result$origin_ids),
                      snapshot$peptide_origins)
    event_ids <- unique(unlist(lapply(origins, function(x)
      unlist(x[c("source_event_ids", "context_event_ids", "germline_context_ids")]))))
    events <- Filter(function(x) x$event_id %in% event_ids, snapshot$events)
    biological <- NULL
    if (!is.null(somatic)) {
      link <- Filter(function(x) x$candidate_id==result$candidate_id,
        somatic$raw$candidate_links)
      eligible <- Filter(function(x) x$candidate_id==result$candidate_id,
        somatic$eligibility$candidate_eligibility)
      biological <- list(candidate_link=link[[1L]],eligibility=eligible[[1L]],
        event_evidence=Filter(function(x) x$event_id %in% event_ids,
          somatic$raw$event_evidence),
        expression_evidence=Filter(function(x) x$evidence_id %in%
          unlist(lapply(link[[1L]]$origin_paths,`[[`,"expression_evidence_ids")),
          somatic$raw$expression_evidence))
    }
    item <- list(result = result, origins = origins, events = events,
         prediction = if (is.null(result$primary_prediction_ref)) NULL else
           Filter(function(x) identical(x$prediction_ref, result$primary_prediction_ref),
                  snapshot$run_manifest$predictions))
    if (!is.null(biological)) item$somatic_evidence <- biological
    if (!is.null(priority)) {
      item$priority <- result[c("final_grade","grade_status","grade_reason_codes",
        "policy_version","technical_eligibility","technical_reason_records",
        "primary_prediction_status","primary_selection","primary_prediction",
        "primary_applicable_rows","primary_status_records","pmhc_band",
        "pmhc_comparisons","rna_status","rna_support","rna_vaf",
        "transcript_tpm","gene_tpm","ccf_status","ccf_value","ccf_solutions",
        "purity","local_cn","generic_loh","exact_self_match",
        "self_match_status","secondary_predictions","evaluable_provider_count",
        "predictor_disagreement","source_event_ids","priority_origin_id",
        "priority_origin_selection_policy","priority_source_event_ids",
        "origin_priority_summaries","rank_within_sample_class")]
      item$protein_provenance <- lapply(origins,function(x) x[c("origin_id",
        "protein_id","sequence_id","mutant_sequence_ref","normal_sequence_ref",
        "sequence_context_id","source_event_ids")])
      item$predictor_core_evidence <- Filter(function(x)
        identical(x$candidate_id,result$candidate_id),core_evidence)
    }
    if (!is.null(hla)) {
      link <- Filter(function(x) identical(x$candidate_id,result$candidate_id),hla$links)
      if (length(link)!=1L) stop("HLA LOH candidate linkage mismatch")
      evidence <- Filter(function(x) identical(x$evidence_id,
        link[[1L]]$presenter_hla_loh_evidence_id),hla$evidence)
      item$hla_loh <- list(presenter_allele=result$hla_allele,
        candidate_link=link[[1L]],allele_evidence=if (length(evidence)) evidence[[1L]] else NULL,
        generic_loh=if (is.null(priority)) NULL else result$generic_loh,
        frozen_priority_note="somatic-priority-v1 grade and rank unchanged",
        summary=hla$summary)
    }
    item
  })
  extension <- source_read(project)
  if (!is.null(extension)) {
    source_rows <- filter_records(source_candidates_with_prediction(project, extension),
      list(candidate_id = candidate_id, peptide_sequence = peptide,
        sample_id = sample_id), reason)
    for (result in source_rows) {
      origins <- Filter(function(x) x$origin_id %in% unlist(result$origin_ids),
        extension$origins)
      meta <- source_result_metadata(extension, result)
      meta$prediction_evidence <- source_prediction_evidence(project,
        result$candidate_id)
      existing <- which(vapply(candidates, function(x)
        identical(x$result$candidate_id, result$candidate_id), TRUE))
      if (length(existing)) {
        i <- existing[[1L]]
        candidates[[i]]$origins <- c(candidates[[i]]$origins, origins)
        candidates[[i]]$source_extension <- meta
      } else candidates[[length(candidates) + 1L]] <- list(
        result = result, origins = origins,
        events = meta$provider_observations, source_extension = meta,
        priority = NULL, prediction = NULL)
    }
  }
  candidates <- long_read_explain_append(candidates, project,
    list(candidate_id = candidate_id, peptide_sequence = peptide, sample_id = sample_id), reason)
  candidates <- recognition_explain_append(candidates, project, snapshot)
  candidates <- presentation_explain_append(candidates, project, snapshot)
  candidates <- longitudinal_explain_append(candidates, project)
  candidates <- vaccine_explain_append(candidates, project)
  explanation <- structure(list(query = list(candidate_id = candidate_id,
    reason = reason, peptide = peptide, sample_id = sample_id), candidates = candidates),
    class = "neo_explanation")
  if (vaccine_enabled(project)) explanation$vaccine_set <- vaccine_read(project)
  if (format == "structured") return(explanation)
  if (!length(candidates)) return("No matching candidate.")
  paste(vapply(candidates, function(x) paste0(x$result$candidate_id, ": ",
    x$result$peptide_sequence, "; prediction=", x$result$prediction_status,
    "; grade=", x$result$final_grade,
    if (is.null(x$priority)) "" else paste0(" (",x$priority$policy_version,
      "; primary=",x$priority$primary_prediction_status,
      if (is.null(x$result$primary_rank)) "" else paste0(" rank ",x$result$primary_rank),
      "; RNA=",x$priority$rna_status,
      "; technical=",x$priority$technical_eligibility,")"),
    if (is.null(x$hla_loh)) "" else paste0("; HLA LOH=",
      x$hla_loh$candidate_link$presenter_hla_loh_status,
      "; escape=",x$hla_loh$candidate_link$presenter_escape_status),
    "; origins=", length(x$origins),
    "; events=", length(x$events), "; reasons=",
    paste(reason_values(x$result), collapse = ","), recognition_explain_text(x), presentation_explain_text(x), longitudinal_explain_text(x), vaccine_explain_text(x)), ""), collapse = "\n")
}
