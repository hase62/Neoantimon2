schema_table_frame <- function(rows, table, contract = schema_contract()) {
  fields <- contract$tables[[table]]$fields
  if (is.null(fields)) stop("Unknown result table", call. = FALSE)
  value <- function(row, field) {
    x <- row[[field$name]]
    if (field$type %in% c("array", "object", "json")) return(list(x))
    if (!is_missing(x)) return(x)
    switch(field$type, string = NA_character_, boolean = NA,
      integer = NA_integer_, number = NA_real_, scalar = NA, NA)
  }
  columns <- lapply(fields, function(field) {
    if (!length(rows)) {
      if (field$type %in% c("array", "object", "json")) return(I(list()))
      return(switch(field$type, string = character(), boolean = logical(),
        integer = integer(), number = numeric(), scalar = logical(), logical()))
    }
    values <- lapply(rows, value, field = field)
    if (field$type %in% c("array", "object", "json")) return(I(lapply(values, `[[`, 1L)))
    unlist(values, recursive = FALSE, use.names = FALSE)
  })
  names(columns) <- vapply(fields, `[[`, "", "name")
  as.data.frame(columns, check.names = FALSE, stringsAsFactors = FALSE,
                optional = TRUE)
}

validate_filter <- function(x, name) {
  if (is.null(x)) return(invisible(NULL))
  if (!is.character(x) || !length(x) || anyNA(x) || any(!nzchar(x)))
    stop(name, " must contain nonempty strings", call. = FALSE)
}

reason_values <- function(row) unlist(row[c("warning_codes", "downgrade_reasons",
  "exclusion_reasons", "not_evaluated_codes", "grade_reason_codes")],
  recursive = TRUE, use.names = FALSE)

filter_records <- function(rows, filters = list(), reason = NULL) {
  for (name in names(filters)) {
    requested <- filters[[name]]
    validate_filter(requested, name)
    if (!is.null(requested)) rows <- Filter(function(row) !is.null(row[[name]]) &&
      length(row[[name]]) == 1L && row[[name]] %in% requested, rows)
  }
  validate_filter(reason, "reason")
  if (!is.null(reason)) rows <- Filter(function(row) any(reason_values(row) %in% reason), rows)
  rows
}

neo_results <- function(project, grade = NULL, prediction_status = NULL,
                        grading_status = NULL, reason = NULL, sample_id = NULL,
                        candidate_id = NULL) {
  snapshot <- read_runtime_snapshot(project)
  priority <- somatic_priority_read(project,snapshot)
  rows <- if (is.null(priority)) snapshot$neoantigen_results else priority$rows
  rows <- filter_records(rows, list(final_grade = grade,
    prediction_status = prediction_status, grading_status = grading_status,
    sample_id = sample_id, candidate_id = candidate_id), reason)
  out <- schema_table_frame(rows, "neoantigen_results")
  if (!is.null(priority)) {
    fields <- c("grade_status","grade_reason_codes","technical_eligibility",
      "primary_prediction_status","primary_predictor","primary_metric","primary_rank",
      "pmhc_band","pmhc_comparisons","rna_status","rna_support","rna_vaf",
      "transcript_tpm","gene_tpm","tumor_vaf","tumor_alt_reads","ccf_status",
      "ccf_value","purity","local_cn","generic_loh","exact_self_match",
      "self_match_status","evaluable_provider_count","predictor_disagreement",
      "rank_within_sample_class","rank_status","policy_version","source_event_ids",
      "priority_origin_id","priority_origin_selection_policy",
      "priority_source_event_ids")
    for (field in fields) {
      values <- lapply(rows,function(r) r[[field]])
      if (field %in% c("grade_reason_codes","pmhc_comparisons","local_cn",
        "generic_loh","source_event_ids","priority_source_event_ids")) out[[field]] <- I(values)
      else if (field %in% c("primary_rank","rna_vaf","transcript_tpm","gene_tpm",
        "tumor_vaf","tumor_alt_reads","ccf_value","purity",
        "rank_within_sample_class","evaluable_provider_count"))
        out[[field]] <- vapply(values,function(x)if (is.null(x)) NA_real_ else as.numeric(x),0)
      else if (field %in% c("rna_support","exact_self_match","predictor_disagreement"))
        out[[field]] <- vapply(values,function(x)if (is.null(x)) NA else isTRUE(x),TRUE)
      else out[[field]] <- vapply(values,function(x)if (is.null(x)) NA_character_ else x,"")
    }
  }
  hla <- hla_loh_read(project,rows)
  if (!is.null(hla)) {
    links <- stats::setNames(hla$links,vapply(hla$links,`[[`,"","candidate_id"))
    fields <- c("presenter_hla_loh_status","presenter_hla_loh_call",
      "presenter_hla_copy_number","presenter_hla_loh_pval_unique",
      "presenter_hla_loh_evidence_id","presenter_escape_status",
      "paired_hla_loss_status","generic_hla_loh_consistency")
    for (field in fields) out[[field]] <- vapply(rows,function(row) {
      value <- links[[row$candidate_id]][[field]]
      if (is.null(value)) return(if (field %in% c("presenter_hla_copy_number",
        "presenter_hla_loh_pval_unique")) NA_real_ else NA_character_)
      value
    },if (field %in% c("presenter_hla_copy_number",
      "presenter_hla_loh_pval_unique")) 0 else "")
  }
  filters <- list(values = list(final_grade = grade,
    prediction_status = prediction_status, grading_status = grading_status,
    sample_id = sample_id, candidate_id = candidate_id), reason = reason)
  out <- presentation_results_append(recognition_results_append(
    long_read_results_append(source_results_append(out, project, filters), project, filters),
    project, snapshot), project, snapshot)
  allo_results_append(vaccine_results_append(longitudinal_results_append(out, project, filters), project), project)
}

neo_events <- function(project, sample_id = NULL, event_id = NULL) {
  rows <- filter_records(read_runtime_snapshot(project)$events,
    list(sample_id = sample_id, event_id = event_id))
  schema_table_frame(rows, "events")
}

neo_origins <- function(project, sample_id = NULL, origin_id = NULL,
                        peptide = NULL) {
  rows <- read_runtime_snapshot(project)$peptide_origins
  if (length(long_read_inputs(project))) {
    source <- source_read(project)
    rows <- c(rows, source$origins, long_read_read(project)$origins)
  }
  rows <- filter_records(rows,
    list(sample_id = sample_id, origin_id = origin_id, peptide_sequence = peptide))
  schema_table_frame(rows, "peptide_origins")
}

neo_multisample <- function(project, patient_id = NULL, peptide = NULL) {
  overlay <- longitudinal_read(project)
  if (!is.null(overlay) && longitudinal_multiple_tumors(overlay))
    return(longitudinal_multisample(overlay, patient_id, peptide))
  rows <- filter_records(read_runtime_snapshot(project)$multisample_summary,
    list(patient_id = patient_id, peptide_sequence = peptide))
  schema_table_frame(rows, "multisample_summary")
}
