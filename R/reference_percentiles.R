# Optional local reference calibration. Never changes provider evidence or ranks.
reference_percentile_method <- function() "strict-better-plus-half-ties-v1"
reference_percentile_keys <- function() c("provider", "predictor_version", "model_version",
  "task", "metric", "unit", "score_direction", "mhc_key", "peptide_length")
reference_model_unknown <- function(x) is.null(x) || !length(x) ||
  tolower(trimws(x)) %in% c("", "unknown", "na", "n/a", "not_provided", "unavailable", "unspecified", "unversioned")

reference_percentile_load <- function(path) {
  bad <- function(message) external_error("INVALID_REFERENCE_PERCENTILES", message)
  if (!is.character(path) || length(path) != 1L || !file.exists(path) || dir.exists(path)) bad("Require a local manifest file")
  manifest_sha256 <- external_sha(path)
  manifest <- jsonlite::read_json(path, simplifyVector = FALSE)
  if (!identical(manifest$schema_version, "reference-percentiles-1.0.0") ||
      !is.list(manifest$distributions) || !is.null(names(manifest$distributions)) ||
      !length(manifest$distributions)) bad("Unknown schema or empty distribution list")
  ids <- keys <- character(); distributions <- list()
  for (d in manifest$distributions) {
    for (name in c("reference_id", setdiff(reference_percentile_keys(), "peptide_length"),
        "reference_set", "source", "source_version", "original_sha256", "scores_file", "scores_sha256", "license"))
      if (!is.character(d[[name]]) || length(d[[name]]) != 1L || is.na(d[[name]]) ||
          !nzchar(trimws(d[[name]]))) bad(paste("Missing or unknown", name))
    if (!identical(d$species, "NCBITaxon:9606") || !identical(d$hla_class, "I")) bad("Only human Class I references are supported")
    if (reference_model_unknown(d$model_version)) bad("Unknown reference model version")
    if (!d$score_direction %in% c("higher_is_better", "lower_is_better")) bad("Invalid direction")
    for (name in c("peptide_length", "count")) if (!is.numeric(d[[name]]) ||
        length(d[[name]]) != 1L || !is.finite(d[[name]]) || d[[name]] < 1 ||
        d[[name]] != floor(d[[name]])) bad(paste("Invalid", name))
    presenters <- peptide_presenters(d$mhc_key)
    if (length(presenters) != 1L || presenters[[1L]]$status != "RESOLVED" ||
        presenters[[1L]]$hla_class != "I" || presenters[[1L]]$mhc_key != d$mhc_key) bad("Require canonical Class I HLA")
    if (any(!grepl("^[a-f0-9]{64}$", c(d$original_sha256, d$scores_sha256)))) bad("Invalid SHA256")
    key <- infra_hash(d[reference_percentile_keys()])
    if (d$reference_id %in% ids || key %in% keys) bad("Duplicate reference ID or matching key")
    file <- resolve_path(d$scores_file, dirname(resolve_path(path)))
    if (!file.exists(file) || dir.exists(file) || external_sha(file) != d$scores_sha256) bad("Reference score hash mismatch")
    lines <- readLines(file, warn = FALSE)
    if (!length(lines) || lines[[1L]] != "score") bad("Score file requires a single score column")
    values <- suppressWarnings(as.numeric(lines[-1L]))
    if (length(values) != d$count) bad("Reference count mismatch")
    if (any(!is.finite(values))) bad("Non-finite reference score")
    if (!identical(external_sha(file), d$scores_sha256)) bad("Reference changed while loading")
    ids <- c(ids, d$reference_id); keys <- c(keys, key)
    distributions[[key]] <- list(metadata = d, scores = sort(values, method = "radix"))
  }
  if (!identical(external_sha(path), manifest_sha256)) bad("Reference manifest changed while loading")
  list(identity = list(manifest_sha256 = manifest_sha256, method = reference_percentile_method(),
    distributions = lapply(distributions[order(names(distributions))], `[[`, "metadata")), distributions = distributions)
}

reference_percentile_value <- function(sorted, value, direction) {
  # Binary insertion bounds keep ties (including valid zero) as frequencies.
  bound <- function(right) {
    lo <- 0L; hi <- length(sorted)
    while (lo < hi) {
      mid <- lo + (hi - lo) %/% 2L + 1L
      if (sorted[[mid]] < value || (right && sorted[[mid]] == value)) lo <- mid else hi <- mid - 1L
    }
    lo
  }
  left <- bound(FALSE); right <- bound(TRUE)
  better <- if (direction == "lower_is_better") left else length(sorted) - right
  100 * (better + 0.5 * (right - left)) / length(sorted)
}

reference_percentile_derive <- function(reference, normalized, matched, graph) {
  evidence_index <- external_group_index(vapply(matched$evidence, `[[`, "", "evidence_id"))
  candidate_index <- external_value_index(vapply(graph$candidates, `[[`, "", "candidate_id"), graph$candidates)
  measurements <- lapply(seq_len(nrow(normalized)), function(i) {
    e <- as.list(normalized[i, , drop = FALSE]); key_values <- e[setdiff(reference_percentile_keys(), "peptide_length")]
    key_values$peptide_length <- nchar(e$peptide_sequence)
    d <- reference$distributions[[infra_hash(key_values)]]
    evidence_id <- paste0("external-evidence:", infra_hash(e))
    hits <- matched$evidence[external_index_get(evidence_index, evidence_id)]
    supported <- e$scope == "peptide_hla" && e$metric_kind == "raw_score" &&
      all(vapply(lapply(vapply(hits, `[[`, "", "candidate_id"), function(id) candidate_index[[id]]),
        function(c) c$hla_class == "I" && c$species == "NCBITaxon:9606", TRUE)) &&
      length(peptide_presenters(e$mhc_key)) == 1L && identical(peptide_presenters(e$mhc_key)[[1L]]$hla_class, "I")
    reason <- if (!supported) "UNSUPPORTED_SCOPE_OR_METRIC" else
      if (reference_model_unknown(e$model_version)) "MODEL_VERSION_UNKNOWN" else
      if (e$status != "evaluated") paste0("SOURCE_", toupper(e$status)) else
      if (is.null(d)) "EXACT_REFERENCE_KEY_NOT_FOUND" else NULL
    list(evidence_id = paste0("external-evidence:", infra_hash(e)), source_row = e$source_row,
      candidate_ids = as.list(sort(unique(vapply(hits, `[[`, "", "candidate_id")))),
      peptide_sequence = e$peptide_sequence, mhc_key = e$mhc_key, matching_key = key_values,
      original_score_text = e$score, raw_score = if (e$status == "evaluated") as.numeric(e$score) else NULL,
      status = if (is.null(reason)) "CALCULATED" else "NOT_CALCULATED", reason = reason,
      reference_percentile = if (is.null(reason)) reference_percentile_value(d$scores, as.numeric(e$score), e$score_direction) else NULL,
      reference_id = if (is.null(d)) NULL else d$metadata$reference_id,
      reference_sha256 = if (is.null(d)) NULL else d$metadata$scores_sha256,
      reference_count = if (is.null(d)) NULL else d$metadata$count,
      method = reference_percentile_method(), percentile_direction = "lower_is_better")
  })
  measurements
}
