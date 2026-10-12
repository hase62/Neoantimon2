# T064 equivalence oracle: unchanged algorithms from b3a078b9d8e090a02a1fe16d9db0e2fefc3abd10.
# Authored reference code; never used by package production paths.
legacy_external_merge_graph <- function(candidates, origins, events = list(), annotations = list()) {
  # An ID shared by extensions denotes the same candidate with more origins.
  ids <- vapply(candidates, `[[`, "", "candidate_id")
  merged <- lapply(sort(unique(ids), method = "radix"), function(id) {
    group <- candidates[ids == id]; first <- group[[1L]]
    for (c in group) for (field in c("sample_id", "peptide_sequence", "mhc_key", "hla_class"))
      if (!identical(c[[field]], first[[field]]))
        external_error("ANNOTATION_CONFLICT", paste(id, field))
    first$origin_ids <- as.list(sort(unique(as.character(unlist(lapply(group,
      `[[`, "origin_ids")))), method = "radix"))
    first$n_origins <- length(first$origin_ids)
    first
  })
  oids <- vapply(origins, `[[`, "", "origin_id")
  for (id in unique(oids[duplicated(oids)])) {
    group <- origins[oids == id]
    if (length(unique(vapply(group, function(o) infra_hash(external_plain(o)), ""))) != 1L)
      external_error("ANNOTATION_CONFLICT", paste("Origin", id))
  }
  origins <- external_sort(origins[!duplicated(oids)], "origin_id")
  links <- unlist(lapply(merged, function(c) lapply(c$origin_ids, function(id) {
    if (!id %in% oids) external_error("MISSING_ORIGIN", paste(c$candidate_id, id))
    list(candidate_id = c$candidate_id, origin_id = id, sample_id = c$sample_id)
  })), recursive = FALSE)
  external_plain(list(candidates = merged, origins = origins, links = links,
    events = external_sort(events, "event_id"), annotations = annotations))
}
legacy_external_match_predictions <-
function (graph, x)
{
    evidence <- list()
    unmatched <- list()
    matched_rows <- integer()
    query_keys <- character()
    for (i in seq_len(nrow(x))) {
        row <- as.list(x[i, , drop = FALSE])
        matches <- Filter(function(c) identical(c$peptide_sequence, row$peptide_sequence) && (row$scope ==
            "peptide" || identical(c$mhc_key, row$mhc_key)) && (!startsWith(row$task, "mhc1_") || c$hla_class ==
            "I") && (!startsWith(row$task, "mhc2_") || c$hla_class == "II") && (!row$task %in% c("processing.class1",
            "transport.tap", "stability.pmhc1") || c$hla_class == "I"), graph$candidates)
        if (!length(matches)) {
            unmatched[[length(unmatched) + 1L]] <- c(row, list(diagnostic = "UNMATCHED_QUERY"))
            next
        }
        matched_rows <- c(matched_rows, i)
        query_keys <- c(query_keys, paste(row$scope, row$peptide_sequence, row$mhc_key, sep = "\037"))
        eid <- paste0("external-evidence:", infra_hash(row))
        for (c in matches) evidence[[length(evidence) + 1L]] <- c(row, list(evidence_id = eid, candidate_id = c$candidate_id,
            sample_id = c$sample_id, origin_ids = c$origin_ids, matched_mhc_key = c$mhc_key, numeric_score = if (row$status ==
                "evaluated") as.numeric(row$score) else NULL))
    }
    matched_ids <- unique(vapply(evidence, `[[`, "", "candidate_id"))
    report <- list(schema_version = external_version(), scopes = as.list(sort(unique(x$scope))), input_rows = nrow(x),
        evaluated_rows = sum(x$status == "evaluated"), failed_rows = sum(x$status == "failed"), unsupported_rows = sum(x$status ==
            "unsupported"), not_evaluated_rows = sum(x$status == "not_evaluated"), matched_input_rows = length(matched_rows),
        matched_unique_query_keys = length(unique(query_keys)), matched_candidates = length(matched_ids),
        fan_out_copies = length(evidence) - length(matched_rows), origin_associations = sum(vapply(evidence,
            function(e) length(e$origin_ids), 0L)), unmatched_input_keys = length(unmatched), missing_candidate_predictions = length(graph$candidates) -
            length(matched_ids), duplicate_ambiguous_keys = 0L)
    list(evidence = evidence, unmatched = unmatched, report = report)
}
legacy_external_derive <-
function (project, graph, matched, selection)
{
    policy <- selection$policy
    selected_keys <- external_logical_keys(selection$selected)
    selected <- Filter(function(e) {
        row <- as.data.frame(e[external_prediction_columns()], stringsAsFactors = FALSE)
        external_logical_keys(row) %in% selected_keys
    }, matched$evidence)
    rows <- lapply(graph$candidates, function(c) {
        own <- Filter(function(e) e$candidate_id == c$candidate_id, matched$evidence)
        chosen <- Filter(function(e) e$candidate_id == c$candidate_id, selected)
        if (length(chosen) > 1L)
            external_error("SELECTION_AMBIGUOUS", "More than one selected value per candidate")
        c$external_status <- if (policy$mode == "evidence_only") {
            if (!length(own))
                "NOT_PREDICTED"
            else if (any(vapply(own, function(e) e$status == "evaluated", TRUE)))
                "evaluated"
            else "NON_EVALUATED_EVIDENCE"
        }
        else if (!length(chosen))
            "NOT_PREDICTED"
        else chosen[[1L]]$status
        c["external_score"] <- list(if (length(chosen)) chosen[[1L]]$numeric_score else NULL)
        c["external_rank"] <- list(NULL)
        c$external_rank_status <- "UNRANKED"
        c$external_policy <- policy$policy_version
        c$external_evidence <- own
        c$external_selected_evidence_id <- if (length(chosen))
            chosen[[1L]]$evidence_id
        else NULL
        c
    })
    if (policy$mode == "external_rank") {
        ready <- Filter(function(c) c$external_status == "evaluated", rows)
        ready <- lapply(ready, function(c) {
            c$final_grade <- "NOT_GRADED"
            c$external_sort_score <- if (policy$score_direction == "higher_is_better")
                -c$external_score
            else c$external_score
            c
        })
        ranked <- somatic_priority_rank(ready, list(ranking_fields = list("external_sort_score", "candidate_id")))
        rank_index <- stats::setNames(lapply(ranked, `[[`, "rank_within_sample_class"), vapply(ranked,
            `[[`, "", "candidate_id"))
        rows <- lapply(rows, function(c) {
            c["external_rank"] <- list(rank_index[[c$candidate_id]])
            c$external_rank_status <- if (is.null(c$external_rank))
                "UNRANKED"
            else "USER_SELECTED_EXTERNAL_RANK"
            c
        })
    }
    if (policy$mode == "somatic_grade") {
        snapshot <- external_completed(project)
        somatic <- somatic_read_sidecar(project, snapshot)
        if (is.null(somatic) || !file.exists(peptide_stage_path(project)))
            external_error("GRADING_PREREQUISITES", "Original somatic technical/evidence sidecars are required")
        candidates <- Filter(function(c) c$hla_class == policy$hla_class, readRDS(peptide_stage_path(project))$candidates)
        grading <- somatic_priority_original_read(project, snapshot)$policy
        grading$policy_version <- policy$policy_version
        grading$primary[[policy$hla_class]] <- list(provider = policy$provider, task = policy$task, metric = policy$metric,
            run_id = policy$run_id, model_version = policy$model_version)
        evidence <- lapply(Filter(function(e) e$status == "evaluated", selected), external_evidence_adapter)
        statuses <- lapply(Filter(function(e) e$status != "evaluated", selected), function(e) list(candidate_id = e$candidate_id,
            task = e$task, status = if (e$status == "not_evaluated") "unavailable" else e$status, reason = e$reason))
        runs <- list(list(provider = policy$provider, run_id = policy$run_id, model_version = policy$model_version,
            evidence = evidence, statuses = statuses))
        graded <- somatic_priority_construct(candidates, somatic, runs, grading)
        index <- stats::setNames(graded, vapply(graded, `[[`, "", "candidate_id"))
        rows <- lapply(rows, function(c) {
            g <- index[[c$candidate_id]]
            c$external_grade_status <- if (is.null(g))
                "NOT_APPLICABLE"
            else g$grade_status
            c$external_grade_origin_coverage <- if (is.null(g))
                "NO_SOMATIC_PREREQUISITES"
            else "CORE_SOMATIC_ORIGINS"
            c["external_final_grade"] <- list(if (is.null(g)) NULL else g$final_grade)
            if (!is.null(g)) {
                c$external_priority <- g
                c["external_rank"] <- list(g$rank_within_sample_class)
                c$external_rank_status <- g$rank_status
            }
            c
        })
    }
    external_plain(rows)
}
legacy_neo_external_explain <-
function (analysis, candidate_id, project = NULL)
{
    validate_filter(candidate_id, "candidate_id")
    if (is.null(candidate_id))
        stop("candidate_id is required")
    x <- external_analysis_read(analysis, project)
    rows <- Filter(function(c) c$candidate_id %in% candidate_id, x$payload$rows)
    list(analysis_id = x$manifest$analysis_id, policy = x$payload$policy, provenance = x$manifest$fingerprint,
        import_report = x$payload$import_report, conversion = x$payload$conversion, candidates = lapply(rows,
            function(c) list(result = c, origins = Filter(function(o) o$origin_id %in% unlist(c$origin_ids),
                x$payload$graph$origins), candidate_origin_links = Filter(function(l) l$candidate_id ==
                c$candidate_id, x$payload$graph$links), evidence = c$external_evidence, events = x$payload$graph$events,
                annotation_provenance = x$payload$graph$annotations)))
}
legacy_neo_export_bundle <-
function (project, dir = NULL)
{
    graph <- external_graph(project)
    if (is.null(dir))
        dir <- file.path(project$path, "results", "export_bundle")
    dir <- resolve_path(dir)
    forbidden <- c(file.path(project$path, "intermediate"), file.path(project$path, "manifest"), project$config$paths$repository,
        vapply(project$inputs, `[[`, "", "resolved_path"), vapply(Filter(function(r) !is.null(r$resolved_path),
            project$resources), `[[`, "", "resolved_path"))
    if (any(vapply(forbidden, function(p) within_path(dir, resolve_path(p)) || within_path(resolve_path(p),
        dir), TRUE)))
        external_error("UNSAFE_PATH", "Export overlaps source inputs")
    if (file.exists(file.path(dir, "manifest.json"))) {
        out <- external_bundle_read(project, dir)
        if (!identical(infra_hash(out$graph), infra_hash(graph)))
            external_error("BUNDLE_MISMATCH", "Persisted annotation graph changed")
        return(structure(out, class = "neo_export_bundle"))
    }
    if (file.exists(dir))
        external_error("EXPORT_EXISTS", "Refusing an existing incomplete export target")
    dir.create(dirname(dir), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(".export-", tmpdir = dirname(dir))
    if (!dir.create(tmp))
        stop("Could not create export staging directory")
    on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
    sequences <- sort(unique(vapply(graph$candidates, `[[`, "", "peptide_sequence")), method = "radix")
    peptides <- lapply(sequences, function(s) list(peptide_sequence = s))
    query_status <- lapply(graph$candidates, function(c) {
        p <- peptide_presenters(c$mhc_key)
        ready <- length(p) == 1L && p[[1L]]$status == "RESOLVED" && identical(p[[1L]]$mhc_key, c$mhc_key)
        list(candidate_id = c$candidate_id, pmhc_query_status = if (ready) "ELIGIBLE" else "UNSUPPORTED_PRESENTER")
    })
    ready_ids <- vapply(Filter(function(s) s$pmhc_query_status == "ELIGIBLE", query_status), `[[`, "",
        "candidate_id")
    queries <- lapply(Filter(function(c) c$candidate_id %in% ready_ids, graph$candidates), function(c) c[c("peptide_sequence",
        "mhc_key", "hla_class", "hla_alpha", "hla_beta", "hla_allele")])
    if (length(queries)) {
        keys <- vapply(queries, function(q) paste(q$peptide_sequence, q$mhc_key, sep = "\037"), "")
        queries <- queries[!duplicated(keys)]
        queries <- queries[order(keys[!duplicated(keys)], method = "radix")]
    }
    external_write_tsv(peptides, "peptide_sequence", file.path(tmp, "peptides.tsv"))
    external_write_tsv(queries, c("peptide_sequence", "mhc_key", "hla_class", "hla_alpha", "hla_beta",
        "hla_allele"), file.path(tmp, "pmhc_queries.tsv"))
    for (table in c("candidates", "origins", "events")) {
        required <- switch(table, candidates = c("candidate_id", "sample_id", "peptide_sequence", "mhc_key",
            "hla_class", "origin_ids"), origins = c("origin_id", "sample_id", "peptide_sequence"), events = c("event_id",
            "sample_id"))
        external_write_tsv(graph[[table]], external_fields(graph[[table]], required), file.path(tmp,
            paste0(table, ".tsv")))
    }
    external_write_tsv(graph$links, c("candidate_id", "origin_id", "sample_id"), file.path(tmp, "candidate_origin_links.tsv"))
    external_write_tsv(query_status, c("candidate_id", "pmhc_query_status"), file.path(tmp, "query_coverage.tsv"))
    source_events <- external_annotation_records(graph, "source", "canonical_events")
    long_read_events <- external_annotation_records(graph, "long_read", "events")
    external_write_tsv(source_events, external_fields(source_events, "source_event_id"), file.path(tmp,
        "source_events.tsv"))
    external_write_tsv(long_read_events, external_fields(long_read_events, "long_read_source_event_id"),
        file.path(tmp, "long_read_events.tsv"))
    atomic_json(graph, file.path(tmp, "annotations.json"))
    files <- c("peptides.tsv", "pmhc_queries.tsv", "candidates.tsv", "origins.tsv", "events.tsv", "source_events.tsv",
        "long_read_events.tsv", "query_coverage.tsv", "candidate_origin_links.tsv", "annotations.json")
    meta <- list(schema_version = external_version(), contract = "external-export-1.0.0", package_version = as.character(utils::packageVersion("Neoantimon2")),
        run_id = project$run_id, analysis_id = project$analysis_id, source_manifest_sha256 = external_sha(manifest_file(project$path)),
        query_semantics = "literal-AA20+U; exact canonical mhc_key; no normalization", coverage = "core/source/long-read; extension events and provenance in annotations.json",
        unsupported_presenter_queries = length(query_status) - length(ready_ids), files = external_records(tmp,
            files))
    meta$bundle_id <- infra_hash(meta)
    atomic_json(meta, file.path(tmp, "manifest.json"))
    if (!file.rename(tmp, dir))
        stop("Could not atomically store export")
    structure(external_bundle_read(project, dir), class = "neo_export_bundle")
}
legacy_external_validate_predictions <-
function (x)
{
    columns <- external_prediction_columns()
    if (!is.data.frame(x) || anyDuplicated(names(x)) || !setequal(names(x), columns))
        external_error("INVALID_COLUMNS", "Use exactly the external-prediction-1.0.0 columns")
    x <- x[columns]
    if (anyNA(x))
        external_error("INVALID_VALUE", "Use empty cells for absent values, not NA")
    x[] <- lapply(x, as.character)
    if (any(vapply(x, function(c) any(grepl("[\t\r\n\037]", c)), TRUE)))
        external_error("INVALID_VALUE", "Control characters are not allowed")
    required <- setdiff(columns, c("mhc_key", "predictor_version", "model_version", "score", "unit",
        "reason"))
    if (any(vapply(x[required], function(c) any(!nzchar(c) | trimws(c) != c), TRUE)))
        external_error("INVALID_VALUE", "Required fields must be nonempty and unpadded")
    if (any(x$schema_version != external_version()))
        external_error("SCHEMA_VERSION", "Unsupported input schema")
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
    if (any(percentile & (x$score_direction != "lower_is_better" | x$unit != "percent")) || any(percentile &
        evaluated & (scores < 0 | scores > 100)))
        external_error("INVALID_PERCENTILE", "Percent ranks require [0,100], percent units and lower_is_better")
    semantics <- c("mhc_binding_percentile", "mhc_presentation_percentile", "external_percentile")
    if (any(percentile & !x$metric_semantics %in% semantics) || any(!percentile & x$metric_semantics !=
        "raw"))
        external_error("INVALID_SEMANTICS", "Metric kind and declared meaning disagree")
    for (i in which(percentile & x$metric_semantics != "external_percentile")) {
        expected <- if (x$metric_semantics[[i]] == "mhc_binding_percentile")
            c("mhc1_binding", "mhc2_binding")
        else c("mhc1_presentation", "mhc2_presentation")
        if (!x$task[[i]] %in% expected)
            external_error("INVALID_SEMANTICS", "MHC percentile task mismatch")
    }
    keys <- external_logical_keys(x)
    if (anyDuplicated(keys))
        external_error("DUPLICATE_PREDICTION", "Duplicate/conflicting logical provider/run/task/metric query")
    channels <- if (nrow(x))
        apply(x[c("provider", "predictor_version", "model_version", "run_id", "task", "metric")], 1L,
            paste, collapse = "\037")
    else character()
    for (ch in unique(channels)) {
        meanings <- unique(x[channels == ch, c("scope", "metric_kind", "metric_semantics", "score_direction",
            "unit")])
        if (nrow(meanings) != 1L)
            external_error("AMBIGUOUS_METRIC", "A metric channel has inconsistent meaning")
    }
    x[order(keys, method = "radix"), , drop = FALSE]
}
