# Human decisions are separate downloadable records; never part of scientific output.
review_normalize <- function(x) {
  if (is.data.frame(x)) x <- report_rows(x)
  if (!is.list(x)) return(x)
  if (!is.null(names(x))) {
    # Reporting row IDs can depend on projection order; biological/evidence IDs
    # and every displayed scientific value remain part of the comparison.
    x <- x[setdiff(names(x), "prediction_record_id")]
    return(lapply(x[sort(names(x), method = "radix")], review_normalize))
  }
  x <- lapply(x, review_normalize)
  # Collections of records/identifiers are unordered for evidence comparison.
  # Numeric positional arrays retain order.
  if (length(x) && all(vapply(x, function(v) is.list(v) || is.character(v), TRUE)))
    x <- x[order(vapply(x, report_json, ""), method = "radix")]
  x
}

review_payload <- function(tables, manifest) {
  rows <- lapply(tables, report_rows)
  candidates <- rows$pmhc_candidates
  records <- lapply(candidates, function(c) {
    id <- c$candidate_id
    links <- Filter(function(l) l$candidate_id == id, rows$candidate_origin_links)
    origins <- Filter(function(o) o$origin_id %in% vapply(links, `[[`, "", "origin_id"), rows$origins)
    event_links <- Filter(function(l) l$origin_id %in% vapply(origins, `[[`, "", "origin_id"), rows$origin_event_links)
    events <- Filter(function(e) e$event_ref %in% vapply(event_links, `[[`, "", "event_ref"), rows$events)
    predictions <- Filter(function(p) p$candidate_id == id, rows$pmhc_predictions)
    reference <- Filter(function(r) id %in% unlist(r$candidate_ids), rows$external_reference_percentiles)
    analyses <- sort(unique(vapply(Filter(function(p) p$source == "External", predictions), `[[`, "", "analysis_id")), method = "radix")
    evidence <- list(version = "candidate-review-evidence-1.0.0", candidate = c,
      origins = origins, candidate_origin_links = links, origin_event_links = event_links,
      events = events, predictions = predictions, reference_percentiles = reference)
    normalized <- review_normalize(evidence)
    list(candidate_id = id, evidence_sha256 = report_hash(normalized), evidence_json = report_json(normalized),
      displayed_external_analysis_ids = as.list(analyses))
  })
  list(schema_version = "candidate-review-1.0.0",
    project_identity = list(project_run_id = manifest$project_run_id,
      project_analysis_id = manifest$project_analysis_id,
      source_manifest_sha256 = manifest$source_manifest$sha256),
    evidence = records)
}
