# T037-T039 overlay. Biological keys never contain observation or subject data.
longitudinal_version <- function() "longitudinal-integration-1.0.0"
longitudinal_namespaces <- function() list(context = "subject-context-id-v1",
  event = "cross-sample-event-key-v1", candidate = "cross-sample-candidate-key-v1",
  clonal = "clonal-evidence-id-v1", trajectory = "trajectory-evidence-id-v1")
longitudinal_id <- function(kind, payload) paste0(kind, ":", infra_hash(list(
  namespace = longitudinal_namespaces()[[kind]], payload = longitudinal_plain(payload))))
longitudinal_roles <- function() c("longitudinal_settings", "clonal_evidence", "pyclone_vi", "clonal_mapping")
longitudinal_plain <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) return(lapply(seq_len(nrow(x)), function(i) lapply(x, function(c) longitudinal_plain(c[[i]]))))
  if (is.list(x)) {
    n <- names(x)
    if (!is.null(n) && (anyNA(n) || anyDuplicated(n) || any(!nzchar(n))))
      return(lapply(seq_along(x), function(i) list(key = n[[i]], value = longitudinal_plain(x[[i]]))))
    return(lapply(x, longitudinal_plain))
  }
  if (length(x) != 1L) return(lapply(as.list(x), longitudinal_plain))
  x
}
longitudinal_inputs <- function(project) Filter(function(i) i$input_role %in% longitudinal_roles(), project$inputs)
longitudinal_sort <- function(x) {
  if (!length(x)) return(list())
  x <- unique(x)
  x[order(vapply(x, function(v) canonical_bytes(longitudinal_plain(v)), ""), method = "radix")]
}
longitudinal_settings <- function(project) {
  inputs <- Filter(function(i) i$input_role == "longitudinal_settings", project$inputs)
  if (!length(inputs)) return(NULL)
  if (length(inputs) != 1L) stop("Exactly one longitudinal settings artifact required")
  x <- jsonlite::read_json(inputs[[1L]]$resolved_path)
  reject_duplicate_keys(x)
  strict_keys(x, c("version", "subjects", "contexts", "relations", "treatments",
    "evaluability", "providers", "sample_runs"), "longitudinal settings")
  if (!identical(x$version, longitudinal_version())) stop("Unsupported longitudinal schema")
  for (k in setdiff(c("subjects", "contexts", "relations", "treatments", "evaluability", "providers", "sample_runs"), names(x))) x[[k]] <- list()
  for (k in setdiff(names(x), "version")) {
    if (!is.list(x[[k]]) || !is.null(names(x[[k]]))) stop("Longitudinal records must be arrays")
    x[[k]] <- longitudinal_sort(x[[k]])
  }
  longitudinal_graph_validate(x, project)
  x
}
longitudinal_graph_validate <- function(x, project) {
  subjects <- vapply(x$subjects, function(s) { scalar_text(s$subject_id, "subject_id"); s$subject_id }, "")
  samples <- vapply(x$contexts, function(s) {
    scalar_text(s$sample_id, "sample_id"); scalar_text(s$subject_id, "subject_id")
    if (!s$subject_id %in% subjects) stop("Unknown explicit subject assignment")
    bound <- Filter(function(b) b$sample_id == s$sample_id, project$samples)
    if (length(bound) != 1L || bound[[1L]]$patient_id != s$subject_id) stop("Subject/sample binding conflict")
    longitudinal_time(s)
    if (!is.null(s$relative_time)) {
      recognition_number(s$relative_time); scalar_text(s$time_unit, "relative time_unit")
      scalar_text(s$time_origin, "relative time_origin")
    }
    if (!is.null(s$ordinal_time)) recognition_number(s$ordinal_time)
    if (!is.null(s$collection_time)) longitudinal_time(s)
    if (!is.null(s$hla_genotype)) peptide_presenters(s$hla_genotype)
    s$sample_id
  }, "")
  if (anyDuplicated(subjects) || anyDuplicated(samples)) stop("Duplicate subject/context")
  if (!setequal(samples, vapply(project$samples, `[[`, "", "sample_id"))) stop("Every sample requires explicit longitudinal context")
  subject <- function(id) x$contexts[[match(id, samples)]]$subject_id
  allowed <- c("same_subject", "same_lesion", "different_region", "primary_to_metastasis",
    "pre_to_post", "serial_sample", "matched_normal", "technical_replicate")
  for (r in x$relations) {
    if (!r$type %in% allowed || !all(c(r$from, r$to) %in% samples) || r$from == r$to) stop("Invalid context relationship")
    if (subject(r$from) != subject(r$to)) stop("Cross-subject relationship prohibited")
  }
  for (r in x$sample_runs) {
    if (!r$sample_id %in% samples) stop("Unknown sample run binding")
    scalar_text(r$path, "completed sample run path")
  }
  if (anyDuplicated(vapply(x$sample_runs, `[[`, "", "sample_id"))) stop("Duplicate sample run")
  if (length(x$sample_runs) && !setequal(vapply(x$sample_runs, `[[`, "", "sample_id"), samples))
    stop("Aggregate projects require a completed run binding for every sample")
  for (t in x$treatments) {
    scalar_text(t$treatment_id, "treatment_id")
    if (!t$subject_id %in% subjects) stop("Unknown treatment subject")
    for (k in c("before_samples", "after_samples")) for (id in t[[k]])
      if (!id %in% samples || subject(id) != t$subject_id) stop("Treatment subject boundary conflict")
  }
  if (anyDuplicated(vapply(x$treatments, `[[`, "", "treatment_id"))) stop("Duplicate treatment ID")
  for (a in x$evaluability) {
    if (!a$sample_id %in% samples || !a$domain %in% c("event", "candidate", "rna", "ms", "repertoire", "hla_loh", "apm")) stop("Invalid evaluability binding/domain")
    scalar_text(a$feature_key, "evaluability feature_key")
    if (!a$state %in% c("evaluated_positive", "evaluated_negative", "not_evaluated", "insufficient_evidence", "not_applicable", "ambiguous")) stop("Invalid evaluability state")
    if (a$state %in% c("evaluated_positive", "evaluated_negative")) {
      if (!isTRUE(a$adequate)) stop("Positive/negative assertion requires explicit adequate opportunity")
      scalar_text(a$assay, "assay"); scalar_text(a$evidence_ref, "evaluability evidence_ref")
      scalar_text(a$policy_version, "evaluability policy_version")
    }
  }
  keys <- vapply(x$evaluability, function(a) paste(a$sample_id, a$domain, a$feature_key), "")
  if (anyDuplicated(keys)) stop("Conflicting duplicate evaluability declarations")
  invisible(TRUE)
}
longitudinal_time <- function(s) {
  axes <- c("collection_time", "relative_time", "ordinal_time")
  supplied <- axes[vapply(axes, function(k) !is.null(s[[k]]), TRUE)]
  if (length(supplied) > 1L) stop("Conflicting temporal axes: at most one of collection_time, relative_time, ordinal_time per context")
  if (!is.null(s$collection_time)) {
    v <- s$collection_time
    if (!grepl("^\\d{4}-\\d{2}-\\d{2}$|^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}Z$", v)) stop("Collection time requires ISO date or UTC timestamp")
    date <- as.Date(substr(v, 1, 10), format = "%Y-%m-%d")
    if (is.na(date) || format(date, "%Y-%m-%d") != substr(v, 1, 10)) stop("Invalid collection date")
    stamp <- as.POSIXct(v, format = if (nchar(v) == 10L) "%Y-%m-%d" else "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    if (is.na(stamp)) stop("Invalid collection timestamp")
    return(list(axis = "absolute_utc", value = as.numeric(stamp)))
  }
  if (!is.null(s$relative_time)) return(list(axis = paste("relative", s$time_origin, s$time_unit), value = as.numeric(s$relative_time)))
  if (!is.null(s$ordinal_time)) return(list(axis = "ordinal", value = as.numeric(s$ordinal_time)))
  list(axis = NULL, value = NULL)
}
longitudinal_order <- function(contexts) {
  t <- lapply(contexts, longitudinal_time)
  known <- length(t) > 1L && all(vapply(t, function(t) !is.null(t$value), TRUE)) &&
    length(unique(vapply(t, `[[`, "", "axis"))) == 1L
  i <- if (known) order(vapply(t, `[[`, 0, "value"), vapply(contexts, `[[`, "", "sample_id"), method = "radix") else
    order(vapply(contexts, `[[`, "", "sample_id"), method = "radix")
  list(contexts = contexts[i], ordering = if (known) "explicit_time" else "unordered")
}
longitudinal_event_key <- function(e) {
  if (!is.null(e$isoform_structure_id)) return(e$isoform_structure_id)
  if (!is.null(e$structure)) {
    s <- e$structure; s$sample_id <- NULL
    if (identical(e$source_type, "fusion")) {
      gene <- function(g) if (!is.null(g$id) && nzchar(g$id) && g$id != ".") g$id else g$name
      s$gene_5prime <- gene(s$genes[[1L]]); s$gene_3prime <- gene(s$genes[[2L]]); s$genes <- NULL
    }
    return(longitudinal_id("event", s))
  }
  # Canonical small variants are already normalized by T006.
  if (is.null(e$canonical_event) || is.null(e$contig) || is.null(e$ref) || is.null(e$alt)) return(NULL)
  if (is.null(e$genome_assembly)) return(NULL)
  longitudinal_id("event", list(reference_build = e$genome_assembly,
    canonical_event = e$canonical_event))
}
longitudinal_candidate_key <- function(c) {
  scalar_text(c$candidate_key, "frozen candidate_key")
  c$candidate_key
}
longitudinal_origin_basis <- function(o, event_links) {
  ids <- unique(c(unlist(o[c("source_event_ids", "context_event_ids", "germline_context_ids")]),
    unlist(o$origin_identity$source_event_keys), o$source_antigen_key))
  keys <- lapply(ids, function(id) {
    hit <- Filter(function(e) (e$event_id == id || identical(e$evidence$event_key, id)) && e$sample_id == o$sample_id, event_links)
    keys <- unique(unlist(lapply(hit, `[[`, "feature_key")))
    if (length(keys) != 1L) return(NULL)
    keys[[1L]]
  })
  list(status = if (length(ids) && any(vapply(keys, is.null, TRUE))) "unresolved" else "mapped",
    required_event_keys = as.list(sort(unique(unlist(keys)))), transcript_id = o$transcript_id,
    ligand_sequence = o$ligand_sequence, peptide_start_aa = o$peptide_start_aa,
    peptide_end_aa = o$peptide_end_aa, phase_status = o$phase_status,
    reconstructed_sequence_checksum = o$origin_identity$reconstructed_sequence_checksum)
}
