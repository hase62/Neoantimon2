long_read_result_metadata <- function(extension, candidate) {
  records <- Filter(function(p) p$origin_id %in% unlist(candidate$origin_ids), extension$peptides)
  ids <- unique(vapply(records, `[[`, "", "long_read_source_event_id"))
  events <- Filter(function(e) e$long_read_source_event_id %in% ids, extension$events)
  obs_ids <- unique(unlist(lapply(events, `[[`, "observation_ids")))
  obs <- Filter(function(o) o$long_read_observation_id %in% obs_ids, extension$observations)
  orf_ids <- unique(vapply(records, `[[`, "", "orf_hypothesis_id"))
  orfs <- lapply(Filter(function(h) h$orf$orf_hypothesis_id %in% orf_ids &&
    h$long_read_source_event_id %in% ids, extension$hypotheses), `[[`, "orf")
  join <- function(x) paste(sort(unique(as.character(unlist(x)))), collapse = ";")
  classes <- unlist(lapply(events, function(e) lapply(e$structural_comparison$comparisons, `[[`, "classes")))
  list(source_type = "long_read_isoform", isoform_structure_id = join(lapply(events, `[[`, "isoform_structure_id")),
    long_read_source_event_id = join(ids), provider_observation_count = as.integer(length(obs)),
    structural_novelty_class = join(classes), orf_hypothesis_id = join(orf_ids),
    orf_evidence_type = join(lapply(orfs, `[[`, "evidence_type")),
    orf_start_type = "AUG", orf_start_nt = join(lapply(orfs, `[[`, "start_nt")),
    orf_stop_state = join(lapply(orfs, `[[`, "termination")),
    matched_normal_status = join(lapply(events, `[[`, "matched_normal_status")),
    long_read_provider = join(lapply(obs, `[[`, "provider")),
    sequence_comparison_status = join(lapply(obs, function(o) o$fasta_comparison$status)),
    long_read_classification = join(lapply(obs, function(o) c(
      lapply(o$metadata, function(m) c(m$raw_fields$structural_category, m$raw_fields$classification)),
      o$model$attributes$transcript_type))),
    long_read_support = lapply(obs, function(o) list(observation_id = o$long_read_observation_id,
      provider = o$provider, fields = lapply(o$metadata, function(m) m$raw_fields))),
    provider_observations = obs, source_events = events, orf_hypotheses = orfs,
    normal_observations = Filter(function(o) o$long_read_observation_id %in%
      unlist(lapply(events, `[[`, "normal_observation_ids")), extension$observations),
    peptide_origins = records, grading_status = "NOT_GRADED",
    grade_reason = "LONG_READ_PRIORITY_POLICY_NOT_IMPLEMENTED")
}

long_read_results_append <- function(out, project, filters) {
  extension <- long_read_read(project)
  if (is.null(extension)) return(out)
  all <- filter_records(source_candidates_with_prediction(project, extension), filters$values, filters$reason)
  rows <- Filter(function(r) !r$candidate_id %in% out$candidate_id, all)
  added <- schema_table_frame(rows, "neoantigen_results", peptide_contract())
  for (name in setdiff(names(out), names(added))) added[[name]] <- out[[name]][rep(NA_integer_, nrow(added))]
  for (name in setdiff(names(added), names(out))) out[[name]] <- added[[name]][rep(NA_integer_, nrow(out))]
  out <- rbind(out, added[names(out)])
  fields <- c("source_type", "isoform_structure_id", "long_read_source_event_id",
    "provider_observation_count", "structural_novelty_class", "orf_hypothesis_id",
    "orf_evidence_type", "orf_start_type", "orf_start_nt", "orf_stop_state",
    "matched_normal_status", "long_read_provider", "long_read_classification", "sequence_comparison_status", "grade_reason")
  metas <- lapply(out$candidate_id, function(id) {
    hit <- Filter(function(r) r$candidate_id == id, all)
    if (length(hit)) long_read_result_metadata(extension, hit[[1L]]) else NULL
  })
  for (field in fields) {
    if (!field %in% names(out)) out[[field]] <- if (field == "provider_observation_count")
      rep(NA_integer_, nrow(out)) else rep(NA_character_, nrow(out))
    out$long_read_support <- I(lapply(metas, function(m) if (is.null(m)) NULL else m$long_read_support))
  for (i in seq_along(metas)) if (!is.null(metas[[i]])) {
      existing <- out[[field]][[i]]; value <- metas[[i]][[field]]
      if (field %in% c("source_type", "orf_hypothesis_id") && !is.na(existing) && nzchar(existing))
        value <- paste(sort(unique(c(strsplit(existing, ";", fixed = TRUE)[[1L]], value))), collapse = ";")
      # Shared core/source candidates keep their frozen grade applicability.
      if (field == "grade_reason" && out$final_grade[[i]] != "NOT_GRADED") next
      if (field == "grade_reason" && !is.na(existing) && nzchar(existing)) next
      out[[field]][[i]] <- value
    }
  }
  out$long_read_support <- I(lapply(metas, function(m) if (is.null(m)) NULL else m$long_read_support))
  for (i in seq_along(metas)) if (!is.null(metas[[i]])) {
    hit <- Filter(function(r) r$candidate_id == out$candidate_id[[i]], all)[[1L]]
    out$origin_ids[[i]] <- as.list(unique(c(unlist(out$origin_ids[[i]]), unlist(hit$origin_ids))))
  }
  out
}
long_read_explain_append <- function(candidates, project, filters, reason) {
  extension <- long_read_read(project)
  if (is.null(extension)) return(candidates)
  rows <- filter_records(source_candidates_with_prediction(project, extension), filters, reason)
  for (result in rows) {
    meta <- long_read_result_metadata(extension, result)
    meta$prediction_evidence <- source_prediction_evidence(project, result$candidate_id)
    meta$self_comparison <- extension$self_matches[[result$ligand_sequence]]
    origins <- Filter(function(o) o$origin_id %in% unlist(result$origin_ids), extension$origins)
    hit <- which(vapply(candidates, function(c) c$result$candidate_id == result$candidate_id, TRUE))
    if (length(hit)) {
      i <- hit[[1L]]; candidates[[i]]$origins <- c(candidates[[i]]$origins, origins)
      candidates[[i]]$long_read_extension <- meta
    } else candidates[[length(candidates) + 1L]] <- list(result = result, origins = origins,
      events = meta$provider_observations, long_read_extension = meta, priority = NULL, prediction = NULL)
  }
  candidates
}
