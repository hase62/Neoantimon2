# Independent interoperability envelopes; never part of biological identity.
external_version <- function() "1.0.0"
external_plain <- function(x) jsonlite::fromJSON(jsonlite::toJSON(x,
  auto_unbox = TRUE, null = "null", digits = NA), simplifyVector = FALSE)
external_sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
external_error <- function(code, message) {
  e <- simpleError(paste0(code, ": ", message))
  attr(e, "neo_error_code") <- code
  stop(e)
}
external_sort <- function(rows, field) {
  if (!length(rows)) return(list())
  rows[order(vapply(rows, `[[`, "", field), method = "radix")]
}
external_fields <- function(rows, required) sort(unique(c(required,
  unlist(lapply(rows, names), use.names = FALSE))), method = "radix")
external_write_tsv <- function(rows, fields, path) utils::write.table(
  somatic_flat(rows, fields), path, sep = "\t", row.names = FALSE,
  col.names = TRUE, quote = TRUE, qmethod = "double", na = "")

external_completed <- function(project) {
  if (!inherits(project, "neo_project") || is.null(project$path))
    external_error("INCOMPLETE_RUN", "Require an initialized neo_project")
  check_run_location(project$path, project$config)
  for (name in c("results", "cache")) if (!within_path(resolve_path(file.path(project$path, name)),
      resolve_path(project$path))) external_error("UNSAFE_PATH", "Derived artifact namespace escapes run")
  snapshot <- read_runtime_snapshot(project)
  if (read_manifest(project$path)$status != "COMPLETED")
    external_error("INCOMPLETE_RUN", "External operations require a completed run")
  snapshot
}

external_graph <- function(project) {
  snapshot <- external_completed(project)
  if (longitudinal_aggregate(project)) {
    graphs <- longitudinal_sample_graphs(project, longitudinal_settings(project))
    parts <- lapply(graphs, function(g) {
      x <- external_graph(g$project)
      x$candidates <- Filter(function(c) c$sample_id == g$sample_id, x$candidates)
      x$origins <- Filter(function(o) o$sample_id == g$sample_id, x$origins)
      x$events <- Filter(function(e) e$sample_id == g$sample_id, x$events)
      x
    })
    candidates <- unlist(lapply(parts, `[[`, "candidates"), recursive = FALSE)
    origins <- unlist(lapply(parts, `[[`, "origins"), recursive = FALSE)
    events <- unlist(lapply(parts, `[[`, "events"), recursive = FALSE)
    annotations <- list(sample_graphs = parts)
  } else {
    source <- source_read(project); lr <- long_read_read(project)
    priority <- somatic_priority_original_read(project, snapshot)
    candidates <- c(if (is.null(priority)) snapshot$neoantigen_results else priority$rows,
      source$candidates, lr$candidates)
    origins <- c(snapshot$peptide_origins, source$origins, lr$origins)
    events <- snapshot$events
    peptide <- if (file.exists(peptide_stage_path(project))) readRDS(peptide_stage_path(project)) else NULL
    annotations <- list(core_peptides = peptide$peptides, source = source, long_read = lr)
  }
  external_merge_graph(candidates, origins, events, annotations)
}

external_merge_graph <- function(candidates, origins, events = list(), annotations = list()) {
  # An ID shared by extensions denotes the same candidate with more origins.
  ids <- vapply(candidates, `[[`, "", "candidate_id")
  candidate_index <- external_group_index(ids)
  merged <- lapply(sort(unique(ids), method = "radix"), function(id) {
    group <- candidates[external_index_get(candidate_index, id)]; first <- group[[1L]]
    for (c in group) for (field in c("sample_id", "peptide_sequence", "mhc_key", "hla_class"))
      if (!identical(c[[field]], first[[field]]))
        external_error("ANNOTATION_CONFLICT", paste(id, field))
    first$origin_ids <- as.list(sort(unique(as.character(unlist(lapply(group,
      `[[`, "origin_ids")))), method = "radix"))
    first$n_origins <- length(first$origin_ids)
    first
  })
  oids <- vapply(origins, `[[`, "", "origin_id")
  origin_index <- external_group_index(oids)
  for (id in unique(oids[duplicated(oids)])) {
    group <- origins[external_index_get(origin_index, id)]
    if (length(unique(vapply(group, function(o) infra_hash(external_plain(o)), ""))) != 1L)
      external_error("ANNOTATION_CONFLICT", paste("Origin", id))
  }
  origins <- external_sort(origins[!duplicated(oids)], "origin_id")
  links <- unlist(lapply(merged, function(c) lapply(c$origin_ids, function(id) {
    if (!exists(id, envir = origin_index, inherits = FALSE))
      external_error("MISSING_ORIGIN", paste(c$candidate_id, id))
    list(candidate_id = c$candidate_id, origin_id = id, sample_id = c$sample_id)
  })), recursive = FALSE)
  external_plain(list(candidates = merged, origins = origins, links = links,
    events = external_sort(events, "event_id"), annotations = annotations))
}

external_records <- function(root, files) lapply(files, function(f)
  list(file = f, sha256 = external_sha(file.path(root, f))))

external_annotation_records <- function(graph, family, table) {
  rows <- graph$annotations[[family]][[table]]
  for (part in graph$annotations$sample_graphs)
    rows <- c(rows, external_annotation_records(part, family, table))
  if (length(rows)) rows <- rows[!duplicated(vapply(rows, function(r) infra_hash(external_plain(r)), ""))]
  rows
}

external_verify_files <- function(root, records) {
  files <- vapply(records, `[[`, "", "file")
  if (anyDuplicated(files)) external_error("INTEGRITY_FAILURE", "Duplicate artifact names")
  for (r in records) {
    # Only flat, explicitly named artifacts can be referenced by an envelope.
    if (!grepl("^[A-Za-z0-9_.-]+$", r$file) || r$file %in% c(".", ".."))
      external_error("UNSAFE_PATH", "Invalid artifact name")
    path <- resolve_path(file.path(root, r$file))
    if (!within_path(path, resolve_path(root)) || !file.exists(path) || dir.exists(path) ||
        !identical(external_sha(path), r$sha256))
      external_error("INTEGRITY_FAILURE", paste("Artifact", r$file))
  }
  invisible(TRUE)
}

external_bundle_read <- function(project, root) {
  external_completed(project)
  meta <- jsonlite::read_json(file.path(root, "manifest.json"))
  if (!identical(meta$schema_version, external_version()) ||
      !identical(meta$run_id, project$run_id) ||
      !identical(meta$analysis_id, project$analysis_id) ||
      !identical(meta$source_manifest_sha256, external_sha(manifest_file(project$path))))
    external_error("BUNDLE_MISMATCH", "Export source/version mismatch")
  identity <- meta; identity$bundle_id <- NULL
  if (!identical(meta$bundle_id, infra_hash(identity)))
    external_error("INTEGRITY_FAILURE", "Bundle envelope mismatch")
  external_verify_files(root, meta$files)
  graph <- jsonlite::read_json(file.path(root, "annotations.json"))
  list(path = resolve_path(root), manifest = meta, graph = graph)
}

neo_export_bundle <- function(project, dir = NULL) {
  graph <- external_graph(project)
  if (is.null(dir)) dir <- file.path(project$path, "results", "export_bundle")
  dir <- resolve_path(dir)
  # Export may be copied outside runs, but never overwrite or overlap originals.
  forbidden <- c(file.path(project$path, "intermediate"), file.path(project$path, "manifest"),
    project$config$paths$repository, vapply(project$inputs, `[[`, "", "resolved_path"),
    vapply(Filter(function(r) !is.null(r$resolved_path), project$resources), `[[`, "", "resolved_path"))
  if (any(vapply(forbidden, function(p) within_path(dir, resolve_path(p)) ||
      within_path(resolve_path(p), dir), TRUE))) external_error("UNSAFE_PATH", "Export overlaps source inputs")
  if (file.exists(file.path(dir, "manifest.json"))) {
    out <- external_bundle_read(project, dir)
    if (!identical(infra_hash(out$graph), infra_hash(graph)))
      external_error("BUNDLE_MISMATCH", "Persisted annotation graph changed")
    return(structure(out, class = "neo_export_bundle"))
  }
  if (file.exists(dir)) external_error("EXPORT_EXISTS", "Refusing an existing incomplete export target")
  dir.create(dirname(dir), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(".export-", tmpdir = dirname(dir))
  if (!dir.create(tmp)) stop("Could not create export staging directory")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  sequences <- sort(unique(vapply(graph$candidates, `[[`, "", "peptide_sequence")), method = "radix")
  peptides <- lapply(sequences, function(s) list(peptide_sequence = s))
  query_status <- lapply(graph$candidates, function(c) {
    p <- peptide_presenters(c$mhc_key)
    ready <- length(p) == 1L && p[[1L]]$status == "RESOLVED" && identical(p[[1L]]$mhc_key, c$mhc_key)
    list(candidate_id = c$candidate_id, pmhc_query_status = if (ready) "ELIGIBLE" else "UNSUPPORTED_PRESENTER")
  })
  eligible <- vapply(query_status, function(s) s$pmhc_query_status == "ELIGIBLE", TRUE)
  queries <- lapply(graph$candidates[eligible], function(c) c[c("peptide_sequence", "mhc_key",
    "hla_class", "hla_alpha", "hla_beta", "hla_allele")])
  if (length(queries)) {
    keys <- vapply(queries, function(q) paste(q$peptide_sequence, q$mhc_key, sep = "\037"), "")
    queries <- queries[!duplicated(keys)]
    queries <- queries[order(keys[!duplicated(keys)], method = "radix")]
  }
  external_write_tsv(peptides, "peptide_sequence", file.path(tmp, "peptides.tsv"))
  external_write_tsv(queries, c("peptide_sequence", "mhc_key", "hla_class", "hla_alpha",
    "hla_beta", "hla_allele"), file.path(tmp, "pmhc_queries.tsv"))
  for (table in c("candidates", "origins", "events")) {
    required <- switch(table, candidates = c("candidate_id", "sample_id", "peptide_sequence",
      "mhc_key", "hla_class", "origin_ids"), origins = c("origin_id", "sample_id", "peptide_sequence"),
      events = c("event_id", "sample_id"))
    external_write_tsv(graph[[table]], external_fields(graph[[table]], required),
      file.path(tmp, paste0(table, ".tsv")))
  }
  external_write_tsv(graph$links, c("candidate_id", "origin_id", "sample_id"),
    file.path(tmp, "candidate_origin_links.tsv"))
  external_write_tsv(query_status, c("candidate_id", "pmhc_query_status"), file.path(tmp, "query_coverage.tsv"))
  # Extension events have their own canonical ID namespaces, preserved verbatim.
  source_events <- external_annotation_records(graph, "source", "canonical_events")
  long_read_events <- external_annotation_records(graph, "long_read", "events")
  external_write_tsv(source_events,
    external_fields(source_events, "source_event_id"), file.path(tmp, "source_events.tsv"))
  external_write_tsv(long_read_events,
    external_fields(long_read_events, "long_read_source_event_id"), file.path(tmp, "long_read_events.tsv"))
  atomic_json(graph, file.path(tmp, "annotations.json"))
  files <- c("peptides.tsv", "pmhc_queries.tsv", "candidates.tsv", "origins.tsv",
    "events.tsv", "source_events.tsv", "long_read_events.tsv", "query_coverage.tsv", "candidate_origin_links.tsv", "annotations.json")
  meta <- list(schema_version = external_version(), contract = "external-export-1.0.0",
    package_version = as.character(utils::packageVersion("Neoantimon2")),
    run_id = project$run_id, analysis_id = project$analysis_id,
    source_manifest_sha256 = external_sha(manifest_file(project$path)),
    query_semantics = "literal-AA20+U; exact canonical mhc_key; no normalization",
    coverage = "core/source/long-read; extension events and provenance in annotations.json",
    unsupported_presenter_queries = length(query_status) - sum(eligible),
    files = external_records(tmp, files))
  meta$bundle_id <- infra_hash(meta)
  atomic_json(meta, file.path(tmp, "manifest.json"))
  if (!file.rename(tmp, dir)) stop("Could not atomically store export")
  structure(external_bundle_read(project, dir), class = "neo_export_bundle")
}
