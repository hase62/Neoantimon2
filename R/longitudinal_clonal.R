longitudinal_clonal_import <- function(project, settings, event_links) {
  out <- list()
  for (input in Filter(function(i) i$input_role %in% c("pyclone_vi", "clonal_evidence"), project$inputs)) {
    cfg <- Filter(function(p) p$input_id == input$input_id, settings$providers)
    if (length(cfg) != 1L) stop("Explicit clonal provider binding required")
    cfg <- cfg[[1L]]
    for (k in c("provider", "provider_version", "run_id", "subject_id")) scalar_text(cfg[[k]], k)
    if (!cfg$subject_id %in% vapply(settings$subjects, `[[`, "", "subject_id")) stop("Unknown clonal subject")
    if (input$input_role == "pyclone_vi" && cfg$provider != "PyClone-VI") stop("PyClone-VI provider mismatch")
    rows <- presentation_table(input$resolved_path, c("mutation_id", "sample_id",
      if (input$input_role == "pyclone_vi") "cluster_id"))
    mapping <- list()
    mapping_artifact <- NULL
    mapping_semantics <- NULL
    if (!is.null(cfg$mapping_input_id)) {
      m <- Filter(function(i) i$input_id == cfg$mapping_input_id && i$input_role == "clonal_mapping", project$inputs)
      if (length(m) != 1L) stop("Unresolved clonal mapping artifact")
      mapping <- presentation_table(m[[1L]]$resolved_path, c("mutation_id", "sample_id", "event_id"))
      mapping_artifact <- recognition_artifact(m[[1L]])
      mapping_semantics <- c(mapping_artifact[c("input_id", "sample_id", "role", "version")],
        list(rows = longitudinal_sort(mapping)))
    }
    for (r in rows) {
      raw <- r
      if (!is.null(r$subject_id) && !identical(r$subject_id, cfg$subject_id)) stop("Clonal row subject conflict")
      provider_sample <- r$sample_id
      if (!is.null(cfg$sample_mapping)) {
        mapped_sample <- cfg$sample_mapping[[provider_sample]]
        if (is.null(mapped_sample)) stop("Unmapped provider sample")
        r$sample_id <- mapped_sample
      }
      ct <- Filter(function(s) s$sample_id == r$sample_id && s$subject_id == cfg$subject_id, settings$contexts)
      if (length(ct) != 1L) stop("Clonal evidence subject/sample conflict")
      mutations <- Filter(function(e) e$sample_id == r$sample_id && e$subject_id == cfg$subject_id, event_links)
      ids <- if (is.null(cfg$mapping_input_id)) if (is.null(r$event_id)) r$mutation_id else r$event_id else unique(vapply(Filter(function(m)
        m$sample_id == r$sample_id && m$mutation_id == r$mutation_id, mapping), `[[`, "", "event_id"))
      hits <- Filter(function(e) e$event_id %in% ids, mutations)
      hit_ids <- sort(unique(vapply(hits, `[[`, "", "event_id")))
      status <- if (length(ids) > 1L || length(hit_ids) > 1L) "ambiguous" else
        if (length(hit_ids) == 1L) "mapped" else "unmapped"
      for (f in c("cellular_prevalence", "cellular_prevalence_std", "cluster_assignment_prob")) r[f] <- list(recognition_number(r[[f]]))
      for (f in c("cellular_prevalence", "cluster_assignment_prob"))
        if (!is.null(r[[f]]) && (r[[f]] < 0 || r[[f]] > 1)) stop("Provider proportion out of range")
      if (!is.null(r$cellular_prevalence_std) && r$cellular_prevalence_std < 0) stop("Negative provider uncertainty")
      scalar_text(r$mutation_id, "provider mutation_id")
      if (!is.null(r$cluster_id)) scalar_text(r$cluster_id, "provider cluster_id")
      record <- list(subject_id = cfg$subject_id, sample_id = r$sample_id,
        provider = cfg$provider, provider_version = cfg$provider_version, provider_run_id = cfg$run_id,
        configuration = cfg$options, reference_build = cfg$reference_build,
        sample_mapping = cfg$sample_mapping, mapping_policy = if (is.null(cfg$mapping_input_id)) "direct-event-id-v1" else "explicit-table-v1",
        mapping_input_id = cfg$mapping_input_id, mapping_artifact = mapping_artifact,
        mapping_semantics = mapping_semantics, artifact = recognition_artifact(input),
        provider_sample_id = provider_sample, requested_event_ids = as.list(sort(ids)),
        provider_mutation_id = r$mutation_id, provider_cluster_id = r$cluster_id,
        cluster_scope = list(provider = cfg$provider, provider_version = cfg$provider_version,
          run_id = cfg$run_id, subject_id = cfg$subject_id, options = cfg$options),
        event_ids = as.list(hit_ids), mapping_status = status,
        cellular_prevalence = r$cellular_prevalence, cellular_prevalence_std = r$cellular_prevalence_std,
        cluster_assignment_prob = r$cluster_assignment_prob, provider_state = r$provider_state,
        raw_fields = raw)
      record$clonal_evidence_id <- longitudinal_id("clonal", longitudinal_clonal_semantics(record))
      out[[length(out) + 1L]] <- record
    }
  }
  out <- unique(out)
  out[order(vapply(out, `[[`, "", "clonal_evidence_id"), method = "radix")]
}

# recognition_artifact exposes role/version. Paths and byte checksums are audit
# provenance; normalized mapping rows and declared metadata identify the import.
longitudinal_clonal_semantics <- function(record) {
  record$artifact <- record$artifact[c("input_id", "sample_id", "role", "version")]
  record$mapping_artifact <- NULL
  record
}
