external_prediction_columns <- function() c("schema_version", "peptide_sequence",
  "scope", "mhc_key", "task", "provider", "predictor_version", "model_version",
  "run_id", "metric", "metric_kind", "metric_semantics", "score", "score_direction",
  "unit", "status", "reason", "source_row")

external_read_predictions <- function(path) {
  if (!file.exists(path) || dir.exists(path)) external_error("INPUT_MISSING", "Prediction table not found")
  # One physical line per record: control characters in values are forbidden.
  # Decode only complete quoted fields, never permissively strip peptide quotes.
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (!length(lines)) external_error("INVALID_TABLE", "A standard header is required")
  decode <- function(value) {
    if (startsWith(value, "\"")) {
      if (!grepl('^"([^"]|"")*"$', value)) external_error("INVALID_TABLE", "Malformed quoted TSV field; double internal quotes")
      return(gsub('""', '"', substring(value, 2L, nchar(value) - 1L), fixed = TRUE))
    }
    if (grepl('"', value, fixed = TRUE)) external_error("INVALID_TABLE", "Quotes inside unquoted TSV fields are forbidden")
    value
  }
  fields <- lapply(lines, function(line) {
    cells <- strsplit(paste0(line, "\tSENTINEL"), "\t", fixed = TRUE)[[1L]]
    vapply(utils::head(cells, -1L), decode, "")
  })
  header <- fields[[1L]]
  if (anyDuplicated(header) || !setequal(header, external_prediction_columns()))
    external_error("INVALID_COLUMNS", "Use exactly the external-prediction-1.0.0 columns")
  rows <- utils::tail(fields, -1L)
  if (any(vapply(rows, length, 0L) != length(header))) external_error("INVALID_TABLE", "Every record must have the header field count")
  columns <- lapply(seq_along(header), function(i) vapply(rows, `[[`, "", i))
  names(columns) <- header
  x <- as.data.frame(columns, check.names = FALSE, stringsAsFactors = FALSE)
  external_validate_predictions(x)
}

external_logical_fields <- function() c("scope", "peptide_sequence", "mhc_key", "provider",
  "predictor_version", "model_version", "run_id", "task", "metric")

external_logical_keys <- function(x) {
  fields <- external_logical_fields()
  if (!nrow(x)) return(character())
  apply(x[fields], 1L, paste, collapse = "\037")
}

external_validate_predictions <- function(x) {
  columns <- external_prediction_columns()
  if (!is.data.frame(x) || anyDuplicated(names(x)) || !setequal(names(x), columns))
    external_error("INVALID_COLUMNS", "Use exactly the external-prediction-1.0.0 columns")
  x <- x[columns]
  if (anyNA(x)) external_error("INVALID_VALUE", "Use empty cells for absent values, not NA")
  x[] <- lapply(x, as.character)
  if (any(vapply(x, function(c) any(grepl("[\t\r\n\037]", c)), TRUE)))
    external_error("INVALID_VALUE", "Control characters are not allowed")
  required <- setdiff(columns, c("mhc_key", "predictor_version", "model_version", "score", "unit", "reason"))
  if (any(vapply(x[required], function(c) any(!nzchar(c) | trimws(c) != c), TRUE)))
    external_error("INVALID_VALUE", "Required fields must be nonempty and unpadded")
  if (any(x$schema_version != external_version())) external_error("SCHEMA_VERSION", "Unsupported input schema")
  if (any(!grepl("^[ACDEFGHIKLMNPQRSTVWYU]+$", x$peptide_sequence)))
    external_error("INVALID_PEPTIDE", "Require literal nonempty uppercase AA20+U")
  if (any(!x$scope %in% c("peptide", "peptide_hla")))
    external_error("UNSUPPORTED_SCOPE", "Only peptide and peptide_hla are supported; context/TCR scopes cannot be broadcast")
  if (any(x$scope == "peptide" & nzchar(x$mhc_key)))
    external_error("INVALID_PRESENTER", "Peptide-only rows must have an empty mhc_key")
  for (key in unique(x$mhc_key[x$scope == "peptide_hla"])) {
    p <- peptide_presenters(key)
    if (length(p) != 1L || p[[1L]]$status != "RESOLVED" || !identical(p[[1L]]$mhc_key, key))
      external_error("INVALID_PRESENTER", "Require full canonical mhc_key, including alpha~beta for class II")
  }
  known <- x$task %in% peptide_prediction_tasks()
  if (any(!known & !grepl("^external:[A-Za-z0-9][A-Za-z0-9_.:-]*$", x$task)))
    external_error("UNSUPPORTED_TASK", "New tasks must use the external: namespace")
  if (any(x$task == "recognition.tcr_specific"))
    external_error("UNSUPPORTED_SCOPE", "TCR-specific evidence requires additional keys")
  intrinsic <- c("processing.class1", "transport.tap", "recognition.peptide_intrinsic")
  if (any(known & !x$task %in% intrinsic & x$scope != "peptide_hla"))
    external_error("UNSUPPORTED_SCOPE", "MHC-dependent tasks require peptide_hla scope")
  if (any(!x$metric_kind %in% c("raw_score", "percentile_rank")))
    external_error("INVALID_METRIC_KIND", "Require raw_score or percentile_rank")
  if (any(!x$score_direction %in% c("higher_is_better", "lower_is_better")))
    external_error("INVALID_DIRECTION", "Explicit score direction required")
  if (any(!x$status %in% c("evaluated", "unsupported", "not_evaluated", "failed")))
    external_error("INVALID_STATUS", "Unknown evaluation status")
  evaluated <- x$status == "evaluated"
  if (any(!evaluated & (nzchar(x$score) | !nzchar(x$reason))))
    external_error("INVALID_SCORE", "Non-evaluated rows need an empty score and explicit reason")
  scores <- suppressWarnings(as.numeric(x$score))
  if (any(evaluated & !is.finite(scores)))
    external_error("INVALID_SCORE", "Evaluated scores must be finite numeric values")
  percentile <- x$metric_kind == "percentile_rank"
  if (any(percentile & (x$score_direction != "lower_is_better" | x$unit != "percent")) ||
      any(percentile & evaluated & (scores < 0 | scores > 100)))
    external_error("INVALID_PERCENTILE", "Percent ranks require [0,100], percent units and lower_is_better")
  semantics <- c("mhc_binding_percentile", "mhc_presentation_percentile", "external_percentile")
  if (any(percentile & !x$metric_semantics %in% semantics) ||
      any(!percentile & x$metric_semantics != "raw"))
    external_error("INVALID_SEMANTICS", "Metric kind and declared meaning disagree")
  for (i in which(percentile & x$metric_semantics != "external_percentile")) {
    expected <- if (x$metric_semantics[[i]] == "mhc_binding_percentile")
      c("mhc1_binding", "mhc2_binding") else c("mhc1_presentation", "mhc2_presentation")
    if (!x$task[[i]] %in% expected) external_error("INVALID_SEMANTICS", "MHC percentile task mismatch")
  }
  keys <- external_logical_keys(x)
  if (anyDuplicated(keys)) external_error("DUPLICATE_PREDICTION", "Duplicate/conflicting logical provider/run/task/metric query")
  # Meaning is constant within one provider/model/run/task/metric channel.
  channels <- if (nrow(x)) apply(x[c("provider", "predictor_version", "model_version", "run_id", "task", "metric")],
    1L, paste, collapse = "\037") else character()
  channel_index <- external_group_index(channels)
  for (ch in unique(channels)) {
    meanings <- unique(x[external_index_get(channel_index, ch),
      c("scope", "metric_kind", "metric_semantics", "score_direction", "unit")])
    if (nrow(meanings) != 1L) external_error("AMBIGUOUS_METRIC", "A metric channel has inconsistent meaning")
  }
  x[order(keys, method = "radix"), , drop = FALSE]
}

external_match_predictions <- function(graph, x) {
  candidates <- graph$candidates
  sequences <- vapply(candidates, `[[`, "", "peptide_sequence")
  mhc <- vapply(candidates, `[[`, "", "mhc_key")
  classes <- vapply(candidates, `[[`, "", "hla_class")
  peptide_index <- external_group_index(sequences)
  pmhc_index <- external_group_index(paste(sequences, mhc, sep = "\037"))
  chunks <- unmatched <- vector("list", nrow(x))
  matched_rows <- logical(nrow(x)); query_keys <- character(nrow(x))
  for (i in seq_len(nrow(x))) {
    row <- as.list(x[i, , drop = FALSE])
    indices <- if (row$scope == "peptide") external_index_get(peptide_index, row$peptide_sequence) else
      external_index_get(pmhc_index, paste(row$peptide_sequence, row$mhc_key, sep = "\037"))
    if (startsWith(row$task, "mhc1_") || row$task %in% c("processing.class1", "transport.tap", "stability.pmhc1"))
      indices <- indices[classes[indices] == "I"]
    if (startsWith(row$task, "mhc2_")) indices <- indices[classes[indices] == "II"]
    matches <- candidates[indices]
    if (!length(matches)) {
      unmatched[[i]] <- c(row, list(diagnostic = "UNMATCHED_QUERY"))
      next
    }
    matched_rows[[i]] <- TRUE
    query_keys[[i]] <- paste(row$scope, row$peptide_sequence, row$mhc_key, sep = "\037")
    # A shared evidence ID denotes ONE measurement, even after biological fan-out.
    eid <- paste0("external-evidence:", infra_hash(row))
    chunks[[i]] <- lapply(matches, function(c) c(row, list(
      evidence_id = eid, candidate_id = c$candidate_id, sample_id = c$sample_id,
      origin_ids = c$origin_ids, matched_mhc_key = c$mhc_key,
      numeric_score = if (row$status == "evaluated") as.numeric(row$score) else NULL)))
  }
  evidence <- unlist(chunks[matched_rows], recursive = FALSE, use.names = FALSE)
  if (!length(evidence)) evidence <- list()
  unmatched <- unmatched[!matched_rows]
  matched_ids <- unique(vapply(evidence, `[[`, "", "candidate_id"))
  report <- list(schema_version = external_version(), scopes = as.list(sort(unique(x$scope))),
    input_rows = nrow(x), evaluated_rows = sum(x$status == "evaluated"),
    failed_rows = sum(x$status == "failed"), unsupported_rows = sum(x$status == "unsupported"),
    not_evaluated_rows = sum(x$status == "not_evaluated"), matched_input_rows = sum(matched_rows),
    matched_unique_query_keys = length(unique(query_keys[matched_rows])), matched_candidates = length(matched_ids),
    fan_out_copies = length(evidence) - sum(matched_rows),
    origin_associations = sum(vapply(evidence, function(e) length(e$origin_ids), 0L)),
    unmatched_input_keys = length(unmatched), missing_candidate_predictions = length(graph$candidates) - length(matched_ids),
    duplicate_ambiguous_keys = 0L)
  list(evidence = evidence, unmatched = unmatched, report = report)
}
