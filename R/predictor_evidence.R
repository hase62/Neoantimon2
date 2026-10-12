# Internal versioned sidecar contract. Frozen T002 snapshot tables are not
# broadened or repurposed for multi-provider raw prediction rows.
predictor_validate_result <- function(result) {
  scalar_text(result$run_id,"prediction run_id")
  scalar_text(result$provider,"prediction provider")
  scalar_text(result$task,"prediction task")
  if (!identical(result$evidence_schema_version,predictor_evidence_schema_version()))
    stop("Unsupported prediction evidence schema version")
  if (!identical(result$adapter_version,predictor_adapter_version()))
    stop("Unsupported prediction adapter version")
  if (!is.list(result$evidence) || !is.list(result$core_evidence) ||
      !is.list(result$statuses) || !is.list(result$raw_artifacts))
    stop("Prediction result arrays required")
  ids <- unlist(result$candidate_ids)
  for (row in result$evidence) {
    required <- c("run_id","analysis_context","sample_id","peptide_id","candidate_id",
      "presenter","external_allele","task","predictor","predictor_version",
      "model_version","model_version_status","metric","raw_score","percentile_rank",
      "score_direction","unit","status","message","raw_output_ref","raw_row_index",
      "predictor_column","adapter_version","evidence_schema_version")
    if (!all(required %in% names(row))) stop("Incomplete prediction evidence row")
    if (!row$candidate_id %in% ids || !identical(row$run_id,result$run_id) ||
        !identical(row$predictor,result$provider) || row$status!="evaluated")
      stop("Prediction evidence identity mismatch")
    if (!is.numeric(row$raw_score) || length(row$raw_score)!=1L ||
        !is.finite(row$raw_score)) stop("Prediction score must be finite")
    if (!row$score_direction %in% c("lower_is_better","higher_is_better"))
      stop("Prediction score direction required")
    if (!row$task %in% peptide_prediction_tasks()) stop("Unknown prediction task")
    if (!is.null(row$percentile_rank) && (!is.numeric(row$percentile_rank) ||
        row$percentile_rank<0 || row$percentile_rank>100)) stop("Invalid prediction percentile")
    scalar_text(row$raw_output_ref,"raw output reference")
    if (!is.integer(row$raw_row_index) || row$raw_row_index<1L)
      stop("Invalid raw output row index")
    if (!identical(row$adapter_version,result$adapter_version) ||
        !identical(row$evidence_schema_version,result$evidence_schema_version))
      stop("Prediction row version mismatch")
  }
  for (core in result$core_evidence) {
    if (!core$candidate_id %in% ids || !identical(core$run_id,result$run_id))
      stop("Core evidence identity mismatch")
    if (!core$status %in% c("evaluated","unresolved_mapping"))
      stop("Invalid core evidence status")
    if (!core$core_mapping_status %in% c("mapped","unresolved") ||
        !core$epitope_key_status %in% c("derived","deferred_orientation","unavailable") ||
        !core$binding_orientation %in% c("forward","reverse","unknown"))
      stop("Invalid core mapping or identity status")
    if (core$core_mapping_status=="mapped") {
      if (core$status!="evaluated" || is.null(core$binding_core_start) ||
          !is.integer(core$binding_core_start) || core$binding_core_start<1L ||
          core$binding_orientation=="unknown") stop("Mapped core coordinate required")
      if (core$binding_orientation=="forward" &&
          (core$epitope_key_status!="derived" || is.null(core$epitope_key)))
        stop("Forward core identity required")
      if (core$binding_orientation=="reverse" &&
          (core$epitope_key_status!="deferred_orientation" || !is.null(core$epitope_key)))
        stop("Reverse core identity must be deferred")
    } else if (core$status!="unresolved_mapping" ||
        !is.null(core$binding_core_start) || !is.null(core$epitope_key) ||
        core$epitope_key_status!="unavailable")
      stop("Unresolved core cannot have coordinate or epitope key")
    scalar_text(core$raw_output_ref,"core raw output reference")
  }
  for (status in result$statuses) {
    if (!status$status %in% c("unsupported","unavailable","failed"))
      stop("Invalid prediction status")
    scalar_text(status$reason,"prediction status reason")
    if (!is.null(status$candidate_id) && !status$candidate_id %in% ids)
      stop("Status candidate identity mismatch")
  }
  invisible(TRUE)
}
