longitudinal_multiple_tumors <- function(x) {
  length(Filter(function(c) !identical(c$role, "normal"), x$context_graph$contexts)) > 1L
}
longitudinal_frame <- function(rows) {
  if (!length(rows)) return(data.frame())
  fields <- sort(unique(unlist(lapply(rows, names))))
  columns <- stats::setNames(lapply(fields, function(f) {
    values <- lapply(rows, `[[`, f)
    scalar <- all(vapply(values, function(v) is.atomic(v) && length(v) == 1L, TRUE))
    if (scalar) unlist(values, use.names = FALSE) else I(values)
  }), fields)
  as.data.frame(columns, stringsAsFactors = FALSE, optional = TRUE)
}
longitudinal_results_append <- function(out, project, filters) {
  x <- longitudinal_read(project)
  if (is.null(x)) return(out)
  if (longitudinal_aggregate(project)) {
    contexts <- unlist(lapply(x$context_graph$subjects, function(s) vapply(longitudinal_order(Filter(function(c)
      c$subject_id == s$subject_id, x$context_graph$contexts))$contexts, `[[`, "", "sample_id")))
    links <- x$candidate_links[order(vapply(x$candidate_links, function(l) match(l$sample_id, contexts), 0L),
      vapply(x$candidate_links, `[[`, 0L, "sample_result_order"))]
    all_rows <- lapply(links, `[[`, "sample_result")
    # Apply the same frozen filters to the referenced sample rows.
    rows <- filter_records(all_rows, filters$values, filters$reason)
    if (length(all_rows)) out <- if (length(rows)) longitudinal_frame(rows) else longitudinal_frame(all_rows)[FALSE, , drop = FALSE]
  }
  # Single-sample summaries do not change the historical result shape.
  if (!longitudinal_multiple_tumors(x)) return(out)
  links <- lapply(out$candidate_id, function(id) Filter(function(l) l$candidate_id == id, x$candidate_links)[[1L]])
  trajectories <- lapply(links, function(l) Filter(function(t) t$subject_id == l$subject_id && t$feature_key == l$feature_key,
    x$candidate_trajectories)[[1L]])
  out$subject_id <- vapply(links, `[[`, "", "subject_id")
  out$cross_sample_candidate_key <- vapply(links, `[[`, "", "feature_key")
  for (f in c("sample_count_observed", "evaluable_sample_count", "trajectory_state"))
    out[[f]] <- vapply(trajectories, function(t) t$summary[[f]], if (f == "trajectory_state") "" else 0L)
  out$clonal_evidence_status <- vapply(links, function(l) if (any(vapply(x$clonal_evidence, function(e)
    e$sample_id == l$sample_id && e$subject_id == l$subject_id && e$mapping_status == "mapped" &&
      any(unlist(e$event_ids) %in% vapply(Filter(function(event) event$sample_id == l$sample_id && event$subject_id == l$subject_id && !is.null(event$feature_key) &&
        event$feature_key %in% unlist(lapply(l$origin_equivalence, `[[`, "required_event_keys")), x$event_links), `[[`, "", "event_id")), TRUE)))
    "mapped_provider_evidence" else "no_mapped_evidence", "")
  out$longitudinal_escape_status <- vapply(trajectories, function(t) t$escape$status, "")
  out
}
longitudinal_multisample <- function(x, patient_id, peptide) {
  validate_filter(patient_id, "patient_id"); validate_filter(peptide, "peptide")
  rows <- lapply(x$candidate_trajectories, function(t) {
    observations <- Filter(function(l) l$subject_id == t$subject_id && l$feature_key == t$feature_key, x$candidate_links)
    c(list(subject_id = t$subject_id, patient_id = t$subject_id, cross_sample_candidate_key = t$feature_key,
      peptide_sequence = observations[[1L]]$ligand_sequence, mhc_key = observations[[1L]]$mhc_key,
      trajectory_evidence_id = t$trajectory_evidence_id, ordering = t$ordering,
      sample_observations = t$cells), t$summary)
  })
  longitudinal_frame(filter_records(rows, list(patient_id = patient_id, peptide_sequence = peptide)))
}
longitudinal_explain_append <- function(items, project) {
  x <- longitudinal_read(project)
  if (is.null(x) || !longitudinal_multiple_tumors(x)) return(items)
  lapply(items, function(item) {
    l <- Filter(function(l) l$candidate_id == item$result$candidate_id, x$candidate_links)
    if (length(l) != 1L) stop("Longitudinal candidate link mismatch")
    l <- l[[1L]]
    t <- Filter(function(t) t$subject_id == l$subject_id && t$feature_key == l$feature_key, x$candidate_trajectories)[[1L]]
    source_keys <- unique(unlist(lapply(l$origin_equivalence, `[[`, "required_event_keys")))
    item$longitudinal <- list(biological_identity = l[c("subject_id", "feature_key", "ligand_sequence", "mhc_key")],
      candidate_trajectory = t, source_event_trajectories = Filter(function(e) e$subject_id == l$subject_id && e$feature_key %in% source_keys, x$event_trajectories),
      clonal_evidence = unlist(lapply(t$cells, `[[`, "clonal_evidence"), recursive = FALSE),
      receptor_trajectories = Filter(function(r) r$subject_id == l$subject_id, x$receptor_trajectories),
      contexts = Filter(function(s) s$subject_id == l$subject_id, x$context_graph$contexts),
      relations = Filter(function(r) r$from %in% unlist(t$sample_ids), x$context_graph$relations),
      treatments = Filter(function(s) s$subject_id == l$subject_id, x$context_graph$treatments),
      hla_genotype = Filter(function(s) s$subject_id == l$subject_id, x$hla_genotype),
      uncertainty = "Missing observations are not negative; ordered descriptive differences do not establish treatment causation.")
    item
  })
}
longitudinal_explain_text <- function(item) {
  x <- item$longitudinal
  if (is.null(x)) return("")
  t <- x$candidate_trajectory
  paste0("\nBIOLOGICAL CANDIDATE IDENTITY: ", x$biological_identity$feature_key,
    "\nSAMPLE OBSERVATIONS: ", paste(vapply(t$cells, function(c) paste(c$sample_id, c$state), ""), collapse = "; "),
    "\nSOURCE EVENT TRAJECTORY: ", length(x$source_event_trajectories), " linked groups",
    "\nCLONAL / CELLULAR-PREVALENCE EVIDENCE: ", length(x$clonal_evidence), " provider records; separate from event CCF",
    "\nRNA / EXPRESSION TRAJECTORY: sample-specific values and uncertainty retained",
    "\nHLA LOH / HLA EXPRESSION: ", t$escape$status,
    "\nAPM / PRESENTATION ESCAPE: component observations retained by sample and class",
    "\nIMMUNOPEPTIDOMICS: observed ligands remain sample-specific; missing MS is not supported absence",
    "\nTCR / REPERTOIRE: ", length(x$receptor_trajectories), " exact profile groups; no disappearance or statistical expansion inference",
    "\nTREATMENT/TIME CONTEXT: ", t$ordering, "; ", length(x$treatments), " exposure records",
    "\nEVALUABILITY / UNCERTAINTY: ", t$summary$trajectory_state, "; ", x$uncertainty)
}
longitudinal_aggregate_explain <- function(project, candidate_id, reason, peptide, sample_id, format) {
  for (n in c("candidate_id", "reason", "peptide", "sample_id")) validate_filter(get(n), n)
  if (all(vapply(list(candidate_id, reason, peptide, sample_id), is.null, TRUE))) stop("At least one explanation query is required")
  x <- longitudinal_read(project)
  links <- Filter(function(l) (is.null(candidate_id) || l$candidate_id %in% candidate_id) &&
    (is.null(sample_id) || l$sample_id %in% sample_id) && (is.null(peptide) || l$ligand_sequence %in% peptide) &&
    (is.null(reason) || any(reason_values(l$sample_result) %in% reason)), x$candidate_links)
  items <- longitudinal_explain_append(lapply(links, `[[`, "evidence"), project)
  items <- vaccine_explain_append(items, project)
  if (format == "text") return(if (!length(items)) "No matching candidate." else paste(vapply(items, function(i)
    paste0(i$result$candidate_id, "; grade=", i$result$final_grade, longitudinal_explain_text(i), vaccine_explain_text(i)), ""), collapse = "\n"))
  out <- structure(list(query = list(candidate_id = candidate_id, reason = reason, peptide = peptide, sample_id = sample_id), candidates = items), class = "neo_explanation")
  if (vaccine_enabled(project)) out$vaccine_set <- vaccine_read(project)
  out
}
