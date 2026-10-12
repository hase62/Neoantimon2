external_policy <- function(policy, predictions) {
  strict_keys(policy, c("mode", "provider", "run_id", "task", "metric", "score_direction",
    "predictor_version", "model_version", "hla_class", "policy_version"), "external policy")
  if (is.null(policy$mode) || !policy$mode %in% c("evidence_only", "external_rank", "somatic_grade"))
    external_error("INVALID_POLICY", "Choose evidence_only, external_rank or somatic_grade")
  version <- if (policy$mode == "somatic_grade") "external-percentile-grade-v1" else "external-analysis-v1"
  if (!is.null(policy$policy_version) && policy$policy_version != version)
    external_error("INVALID_POLICY", "Unsupported policy version")
  policy$policy_version <- version
  if (policy$mode == "evidence_only") {
    if (length(setdiff(names(policy), c("mode", "policy_version"))))
      external_error("INVALID_POLICY", "evidence_only retains all channels without selection")
    return(list(policy = policy, selected = predictions[FALSE, , drop = FALSE]))
  }
  required <- c("provider", "run_id", "task", "metric", "score_direction")
  for (key in required) scalar_text(policy[[key]], paste("External policy", key))
  if (!policy$score_direction %in% c("higher_is_better", "lower_is_better"))
    external_error("INVALID_DIRECTION", "Invalid policy direction")
  if (policy$mode == "external_rank" && !is.null(policy$hla_class))
    external_error("INVALID_POLICY", "External ranks are independently computed for every sample/class")
  selected <- predictions
  for (field in c(required, "predictor_version", "model_version")) {
    if (!is.null(policy[[field]])) {
      if (!is.character(policy[[field]]) || length(policy[[field]]) != 1L || is.na(policy[[field]]))
        external_error("INVALID_POLICY", paste("Invalid", field))
      selected <- selected[selected[[field]] == policy[[field]], , drop = FALSE]
    }
  }
  if (!nrow(selected)) external_error("SELECTION_MISSING", "No channel matches the explicit selection")
  channels <- unique(selected[c("provider", "run_id", "task", "metric", "predictor_version", "model_version")])
  if (nrow(channels) != 1L) external_error("SELECTION_AMBIGUOUS", "Select predictor_version and model_version explicitly")
  policy$predictor_version <- channels$predictor_version[[1L]]
  policy$model_version <- channels$model_version[[1L]]
  if (policy$mode == "somatic_grade") {
    scalar_text(policy$hla_class, "External grading hla_class")
    if (!policy$hla_class %in% c("I", "II")) external_error("INVALID_POLICY", "Choose I or II")
    tasks <- if (policy$hla_class == "I") c("mhc1_binding", "mhc1_presentation") else c("mhc2_binding", "mhc2_presentation")
    expected <- if (endsWith(policy$task, "_binding")) "mhc_binding_percentile" else "mhc_presentation_percentile"
    if (!policy$task %in% tasks || any(selected$scope != "peptide_hla") ||
        any(selected$metric_kind != "percentile_rank") || any(selected$metric_semantics != expected) ||
        policy$score_direction != "lower_is_better")
      external_error("GRADING_PREREQUISITES", "Somatic grading requires declared MHC percent-rank semantics for the selected class/task")
  }
  list(policy = policy, selected = selected)
}

external_evidence_adapter <- function(e) {
  # This adapter is a new interpretation boundary, not native evidence v1.0.1.
  list(run_id = e$run_id, candidate_id = e$candidate_id, task = e$task,
    metric = e$metric, predictor = e$provider, predictor_version = e$predictor_version,
    model_version = e$model_version, raw_score = e$numeric_score,
    percentile_rank = if (e$metric_kind == "percentile_rank") e$numeric_score else NULL,
    status = e$status, score_direction = e$score_direction, unit = e$unit,
    raw_row_index = e$source_row, raw_output_ref = "normalized.tsv",
    evidence_id = e$evidence_id, external_contract = "external-prediction-1.0.0")
}

external_derive <- function(project, graph, matched, selection) {
  policy <- selection$policy
  selected_keys <- external_logical_keys(selection$selected)
  keys <- vapply(matched$evidence, function(e)
    paste(unlist(e[external_logical_fields()], use.names = FALSE), collapse = "\037"), "")
  selected <- matched$evidence[keys %in% selected_keys]
  evidence_index <- external_group_index(vapply(matched$evidence, `[[`, "", "candidate_id"))
  selected_index <- external_group_index(vapply(selected, `[[`, "", "candidate_id"))
  rows <- lapply(graph$candidates, function(c) {
    own <- matched$evidence[external_index_get(evidence_index, c$candidate_id)]
    chosen <- selected[external_index_get(selected_index, c$candidate_id)]
    if (length(chosen) > 1L) external_error("SELECTION_AMBIGUOUS", "More than one selected value per candidate")
    c$external_status <- if (policy$mode == "evidence_only") {
      if (!length(own)) "NOT_PREDICTED" else if (any(vapply(own, function(e) e$status == "evaluated", TRUE)))
        "evaluated" else "NON_EVALUATED_EVIDENCE"
    } else if (!length(chosen)) "NOT_PREDICTED" else chosen[[1L]]$status
    c["external_score"] <- list(if (length(chosen)) chosen[[1L]]$numeric_score else NULL)
    c["external_rank"] <- list(NULL)
    c$external_rank_status <- "UNRANKED"
    c$external_policy <- policy$policy_version
    c$external_evidence <- own
    c$external_selected_evidence_id <- if (length(chosen)) chosen[[1L]]$evidence_id else NULL
    c
  })
  if (policy$mode == "external_rank") {
    ready <- Filter(function(c) c$external_status == "evaluated", rows)
    ready <- lapply(ready, function(c) {
      c$final_grade <- "NOT_GRADED"
      c$external_sort_score <- if (policy$score_direction == "higher_is_better") -c$external_score else c$external_score
      c
    })
    # Reuse the existing deterministic sample/class ordering engine.
    ranked <- somatic_priority_rank(ready, list(ranking_fields = list("external_sort_score", "candidate_id")))
    rank_index <- external_value_index(vapply(ranked, `[[`, "", "candidate_id"),
      lapply(ranked, `[[`, "rank_within_sample_class"))
    rows <- lapply(rows, function(c) {
      c["external_rank"] <- list(rank_index[[c$candidate_id]])
      c$external_rank_status <- if (is.null(c$external_rank)) "UNRANKED" else "USER_SELECTED_EXTERNAL_RANK"
      c
    })
  }
  if (policy$mode == "somatic_grade") {
    snapshot <- external_completed(project)
    somatic <- somatic_read_sidecar(project, snapshot)
    if (is.null(somatic) || !file.exists(peptide_stage_path(project)))
      external_error("GRADING_PREREQUISITES", "Original somatic technical/evidence sidecars are required")
    candidates <- Filter(function(c) c$hla_class == policy$hla_class,
      readRDS(peptide_stage_path(project))$candidates)
    # A separately versioned selection policy reuses the frozen constructor.
    # The native validator and its registry/primary metric constraints are intact.
    grading <- somatic_priority_original_read(project, snapshot)$policy
    grading$policy_version <- policy$policy_version
    grading$primary[[policy$hla_class]] <- list(provider = policy$provider, task = policy$task,
      metric = policy$metric, run_id = policy$run_id, model_version = policy$model_version)
    evidence <- lapply(Filter(function(e) e$status == "evaluated", selected), external_evidence_adapter)
    statuses <- lapply(Filter(function(e) e$status != "evaluated", selected), function(e)
      list(candidate_id = e$candidate_id, task = e$task, status = if (e$status == "not_evaluated") "unavailable" else e$status,
        reason = e$reason))
    runs <- list(list(provider = policy$provider, run_id = policy$run_id,
      model_version = policy$model_version, evidence = evidence, statuses = statuses))
    graded <- somatic_priority_construct(candidates, somatic, runs, grading)
    index <- external_value_index(vapply(graded, `[[`, "", "candidate_id"), graded)
    rows <- lapply(rows, function(c) {
      g <- index[[c$candidate_id]]
      c$external_grade_status <- if (is.null(g)) "NOT_APPLICABLE" else g$grade_status
      c$external_grade_origin_coverage <- if (is.null(g)) "NO_SOMATIC_PREREQUISITES" else "CORE_SOMATIC_ORIGINS"
      c["external_final_grade"] <- list(if (is.null(g)) NULL else g$final_grade)
      if (!is.null(g)) {
        # Preserve the full export origin graph, including extension paths.
        c$external_priority <- g
        c["external_rank"] <- list(g$rank_within_sample_class)
        c$external_rank_status <- g$rank_status
      }
      c
    })
  }
  external_plain(rows)
}

external_analysis_root <- function(project, id) {
  if (!is.character(id) || length(id) != 1L || !grepl("^[a-f0-9]{64}$", id))
    external_error("INVALID_ANALYSIS_ID", "Require the explicit 64-character analysis fingerprint")
  root <- resolve_path(file.path(project$path, "results", "external_predictions", id))
  if (!within_path(root, resolve_path(file.path(project$path, "results"))))
    external_error("UNSAFE_PATH", "Analysis path escapes results")
  root
}

external_analysis_read <- function(analysis, project = NULL) {
  if (inherits(analysis, "neo_external_analysis")) {
    project <- analysis$project; id <- analysis$analysis_id
  } else id <- analysis
  external_completed(project)
  root <- external_analysis_root(project, id)
  meta <- jsonlite::read_json(file.path(root, "manifest.json"))
  identity <- meta; identity$metadata_digest <- NULL
  if (!identical(meta$metadata_digest, infra_hash(identity)) || !identical(meta$analysis_id, id) ||
      !identical(infra_hash(meta$fingerprint), id) || !identical(meta$schema_version, external_version()) ||
      !identical(meta$fingerprint$source_manifest_sha256, external_sha(manifest_file(project$path))))
    external_error("INTEGRITY_FAILURE", "External analysis source/envelope mismatch")
  external_verify_files(root, meta$files)
  payload <- jsonlite::read_json(file.path(root, "analysis.json"))
  if (!identical(payload$policy, meta$fingerprint$policy))
    external_error("INTEGRITY_FAILURE", "Analysis policy mismatch")
  list(path = root, manifest = meta, payload = payload)
}

external_frame <- function(rows, required = c("candidate_id", "external_status", "external_score", "external_rank")) {
  fields <- external_fields(rows, required)
  columns <- lapply(fields, function(field) {
    values <- lapply(rows, `[[`, field)
    present <- Filter(function(v) !is.null(v), values)
    if (any(vapply(present, is.list, TRUE))) return(I(values))
    numeric <- field %in% c("external_score", "external_rank") ||
      (length(present) && all(vapply(present, is.numeric, TRUE)))
    logical <- length(present) && all(vapply(present, is.logical, TRUE))
    vapply(values, function(v) if (is.null(v)) {
      if (numeric) NA_real_ else if (logical) NA else NA_character_
    } else v, if (numeric) 0 else if (logical) TRUE else "")
  })
  names(columns) <- fields
  as.data.frame(columns, stringsAsFactors = FALSE, check.names = FALSE, optional = TRUE)
}

neo_analyze_external <- function(project, predictions, policy, converter = NULL,
                                 converter_id = NULL, registry = NULL, trust = FALSE, timeout = 60) {
  bundle <- neo_export_bundle(project)
  if (!is.null(converter) && !is.null(converter_id)) external_error("INVALID_CONVERTER", "Choose converter or converter_id")
  registry_pin <- NULL
  if (!is.null(converter_id)) {
    if (is.null(registry)) registry <- getOption("Neoantimon2.converter_registry")
    resolved <- external_registry_resolve(project, converter_id, registry, trust)
    converter <- resolved$path; registry_pin <- resolved$provenance
  }
  raw_input_sha256 <- external_sha(predictions)
  conversion <- NULL
  if (!is.null(converter)) {
    converted <- external_convert(project, predictions, converter, bundle, timeout, registry_pin)
    predictions <- converted$path; conversion <- converted$metadata
  }
  normalized <- external_read_predictions(predictions)
  selection <- external_policy(policy, normalized)
  matched <- external_match_predictions(bundle$graph, normalized)
  rows <- external_derive(project, bundle$graph, matched, selection)
  fingerprint <- list(schema_version = external_version(), bundle_id = bundle$manifest$bundle_id,
    source_manifest_sha256 = bundle$manifest$source_manifest_sha256,
    raw_input_sha256 = raw_input_sha256, normalized_sha256 = external_sha(predictions),
    converter = if (is.null(conversion)) NULL else conversion$identity, policy = selection$policy)
  id <- infra_hash(fingerprint)
  root <- external_analysis_root(project, id)
  handle <- structure(list(analysis_id = id, project = project, path = root), class = "neo_external_analysis")
  if (file.exists(file.path(root, "manifest.json"))) {
    external_analysis_read(handle)
    return(handle)
  }
  if (file.exists(root)) external_error("ANALYSIS_EXISTS", "Incomplete derived analysis target")
  dir.create(dirname(root), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(".external-", tmpdir = dirname(root)); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  if (!file.copy(predictions, file.path(tmp, "normalized.tsv"))) stop("Could not copy normalized evidence")
  if (!identical(external_sha(file.path(tmp, "normalized.tsv")), fingerprint$normalized_sha256))
    external_error("INTEGRITY_FAILURE", "Normalized input changed while storing analysis")
  payload <- list(schema_version = external_version(), policy = selection$policy,
    graph = bundle$graph, predictions = normalized, evidence = matched$evidence,
    unmatched = matched$unmatched, import_report = matched$report, rows = rows, conversion = conversion)
  atomic_json(payload, file.path(tmp, "analysis.json"))
  external_write_tsv(rows, external_fields(rows, c("candidate_id", "external_status", "external_score", "external_rank")),
    file.path(tmp, "results.tsv"))
  external_write_tsv(matched$evidence, c(external_prediction_columns(), "evidence_id", "candidate_id", "sample_id", "origin_ids"),
    file.path(tmp, "matched_evidence.tsv"))
  external_write_tsv(matched$unmatched, c(external_prediction_columns(), "diagnostic"), file.path(tmp, "unmatched.tsv"))
  atomic_json(matched$report, file.path(tmp, "import_report.json"))
  writeLines(c("# External prediction analysis", "", paste("Policy:", selection$policy$policy_version),
    paste("Mode:", selection$policy$mode), paste("Source bundle:", bundle$manifest$bundle_id),
    paste("Selected channel:", paste(unlist(selection$policy[c("provider", "run_id", "task", "metric", "score_direction")]), collapse = "; ")),
    "Missing scores remain visible and unranked. Fan-out is one measurement attached to multiple origins.",
    "External ordering is user selected. A-E, where explicitly selected, reuse somatic rules and are not clinical probabilities.",
    "Vaccine, longitudinal, receptor, presentation and allo overlays continue to select native results."),
    file.path(tmp, "report.md"))
  files <- c("normalized.tsv", "analysis.json", "results.tsv", "matched_evidence.tsv", "unmatched.tsv", "import_report.json", "report.md")
  meta <- list(schema_version = external_version(), analysis_id = id,
    fingerprint = fingerprint, files = external_records(tmp, files))
  meta$metadata_digest <- infra_hash(meta)
  atomic_json(meta, file.path(tmp, "manifest.json"))
  if (!file.rename(tmp, root)) stop("Could not atomically store derived analysis")
  external_analysis_read(handle)
  handle
}

neo_external_results <- function(analysis, project = NULL) {
  external_frame(external_analysis_read(analysis, project)$payload$rows)
}

neo_external_explain <- function(analysis, candidate_id, project = NULL) {
  validate_filter(candidate_id, "candidate_id")
  if (is.null(candidate_id)) stop("candidate_id is required")
  x <- external_analysis_read(analysis, project)
  rows <- x$payload$rows[vapply(x$payload$rows, `[[`, "", "candidate_id") %in% candidate_id]
  origins <- x$payload$graph$origins; links <- x$payload$graph$links
  origin_index <- external_group_index(vapply(origins, `[[`, "", "origin_id"))
  link_index <- external_group_index(vapply(links, `[[`, "", "candidate_id"))
  list(analysis_id = x$manifest$analysis_id, policy = x$payload$policy,
    provenance = x$manifest$fingerprint, import_report = x$payload$import_report,
    conversion = x$payload$conversion,
    candidates = lapply(rows, function(c) list(result = c,
      origins = origins[sort(unique(unlist(lapply(c$origin_ids, function(id)
        external_index_get(origin_index, id)), use.names = FALSE)))],
      candidate_origin_links = links[external_index_get(link_index, c$candidate_id)],
      evidence = c$external_evidence,
      events = x$payload$graph$events,
      annotation_provenance = x$payload$graph$annotations)))
}
