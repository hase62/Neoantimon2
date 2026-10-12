longitudinal_validate <- function(x) {
  if (!identical(x$schema_version, longitudinal_version()) || !identical(x$external_clonal_runs, 0L)) stop("Invalid longitudinal overlay boundary")
  contexts <- x$context_graph$contexts
  own <- function(subject, sample) any(vapply(contexts, function(c) c$subject_id == subject && c$sample_id == sample, TRUE))
  for (l in c(x$event_links, x$candidate_links)) if (!own(l$subject_id, l$sample_id)) stop("Cross-subject feature linkage")
  for (l in x$event_links) if (!identical(l$feature_key, longitudinal_event_key(l$evidence))) stop("Cross-sample event identity mismatch")
  for (l in x$candidate_links) {
    if (!identical(l$sample_id, l$sample_result$sample_id) || !identical(l$candidate_id, l$sample_result$candidate_id)) stop("Sample candidate evidence mismatch")
    if (!is.null(l$evidence$result$candidate_key) && !identical(l$feature_key, l$evidence$result$candidate_key)) stop("Frozen candidate key mismatch")
  }
  for (e in x$clonal_evidence) {
    if (!own(e$subject_id, e$sample_id)) stop("Cross-subject clonal evidence")
    if (e$mapping_status == "mapped") for (id in e$event_ids) {
      l <- Filter(function(l) l$event_id == id && l$sample_id == e$sample_id && l$subject_id == e$subject_id, x$event_links)
      if (!length(l)) stop("Clonal mapping crosses sample/subject boundary")
    }
  }
  for (t in c(x$event_trajectories, x$candidate_trajectories, x$receptor_trajectories)) {
    ct <- Filter(function(c) c$subject_id == t$subject_id && !identical(c$role, "normal"), contexts)
    order <- longitudinal_order(ct)
    if (!identical(t$ordering, order$ordering) || !identical(t$sample_ids, as.list(vapply(order$contexts, `[[`, "", "sample_id")))) stop("Trajectory sample ordering mismatch")
    for (cell in t$cells) {
      if (!own(t$subject_id, cell$sample_id)) stop("Cross-subject trajectory")
      if (!cell$state %in% c("observed", "not_observed_evaluable", "not_observed_insufficient", "not_evaluated", "not_applicable", "ambiguous")) stop("Invalid typed observation state")
      if (cell$state == "not_observed_evaluable" && (is.null(cell$evaluability) ||
        !identical(cell$evaluability$state, "evaluated_negative") || !isTRUE(cell$evaluability$adequate))) stop("Unsupported absence claim")
      for (o in cell$observations) if (!identical(o$subject_id, t$subject_id) || !identical(o$sample_id, cell$sample_id) ||
          !identical(o$feature_key, t$feature_key)) stop("Trajectory evidence bleeding")
    }
    summary <- longitudinal_summary(t$cells, order$contexts, x$context_graph$relations)
    if (t$domain == "repertoire") { summary$transitions <- list(); summary$trajectory_state <- "sample_specific_receptor_observations" }
    if (!identical(summary, t$summary)) stop("Trajectory summary derivation mismatch")
    if (!identical(t$trajectory_evidence_id, longitudinal_trajectory_id(t, x$context_graph))) stop("Trajectory identity mismatch")
  }
  invisible(TRUE)
}
