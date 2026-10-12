longitudinal_aggregate <- function(project) {
  s <- longitudinal_settings(project)
  !is.null(s) && length(s$sample_runs) > 0L
}
longitudinal_sample_graphs <- function(project, settings) {
  graphs <- list()
  for (ct in settings$contexts) {
    binding <- Filter(function(r) r$sample_id == ct$sample_id, settings$sample_runs)
    p <- if (length(binding)) neo_project(path = resolve_path(binding[[1L]]$path, project$config$paths$input_root)) else project
    if (length(binding) && (identical(p$path, project$path) || length(longitudinal_inputs(p)))) stop("Nested/cyclic longitudinal sample graph prohibited")
    found <- Filter(function(s) s$sample_id == ct$sample_id && s$patient_id == ct$subject_id, p$samples)
    if (length(found) != 1L) stop("Sample graph subject binding conflict")
    subject <- Filter(function(s) s$subject_id == ct$subject_id, settings$subjects)[[1L]]
    if (!is.null(subject$species) && !identical(subject$species, p$config$scientific$species)) stop("Subject/sample species conflict")
    if (read_manifest(p$path)$status != "COMPLETED") stop("Sample graph must be completed")
    snap <- read_runtime_snapshot(p)
    if (length(binding) && !is.null(p$config$somatic_evidence)) {
      original <- somatic_priority_original_read(p, snap)
      # JSON chooses integer/double differently. Use the frozen policy object
      # only when its scientific fingerprint exactly equals the stored config.
      if (identical(infra_hash(original$policy), infra_hash(p$config$somatic_priority)))
        p$config$somatic_priority <- original$policy
    }
    results <- if (length(binding)) neo_results(p, sample_id = ct$sample_id) else longitudinal_base_results(p, snap)
    candidates <- recognition_candidates(p, snap)
    candidates <- Filter(function(c) c$sample_id == ct$sample_id, candidates)
    events <- Filter(function(e) e$sample_id == ct$sample_id, snap$events)
    somatic <- somatic_read_sidecar(p, snap)
    source <- source_read(p); lr <- long_read_read(p)
    for (o in source$observations) if (identical(o$structure$sample_id, ct$sample_id)) {
      o$event_id <- o$source_event_id; o$sample_id <- ct$sample_id
      events[[length(events) + 1L]] <- o
    }
    for (o in lr$observations) if (identical(o$context$sample_id, ct$sample_id) && !is.null(o$isoform_structure_id)) {
      o$event_id <- o$long_read_source_event_id; o$sample_id <- ct$sample_id
      events[[length(events) + 1L]] <- o
    }
    origins <- Filter(function(o) o$sample_id == ct$sample_id, c(snap$peptide_origins, source$origins, lr$origins))
    recognition <- recognition_read(p, snap)
    graphs[[length(graphs) + 1L]] <- list(sample_id = ct$sample_id, subject_id = ct$subject_id,
      project = p, snapshot = snap, results = results, candidates = candidates, events = events,
      origins = origins, repertoire = recognition$repertoire, quantitative_evidence = somatic$raw$event_evidence,
      binding = list(analysis_id = p$analysis_id, run_id = p$run_id, path = p$path))
  }
  graphs
}
longitudinal_base_results <- function(project, snapshot) {
  # The guard prevents accessor recursion; it does not change sample evidence.
  project$longitudinal_accessor_guard <- TRUE
  neo_results(project)
}
longitudinal_base_explain <- function(project, id) {
  project$longitudinal_accessor_guard <- TRUE
  neo_explain(project, candidate_id = id)
}
longitudinal_construct <- function(project, settings, graphs) {
  event_links <- list(); candidate_links <- list(); receptor_links <- list()
  for (g in graphs) {
    for (e in g$events) event_links[[length(event_links) + 1L]] <- list(subject_id = g$subject_id,
      sample_id = g$sample_id, event_id = e$event_id, feature_key = longitudinal_event_key(e), evidence = e,
      quantitative_evidence = Filter(function(q) q$event_id == e$event_id, g$quantitative_evidence))
    explanations <- if (length(g$candidates)) longitudinal_base_explain(g$project,
      vapply(g$candidates, `[[`, "", "candidate_id"))$candidates else list()
    for (c in g$candidates) {
      item <- Filter(function(i) i$result$candidate_id == c$candidate_id, explanations)
      if (length(item) != 1L) stop("Ambiguous frozen candidate explanation")
      row <- g$results[g$results$candidate_id == c$candidate_id, , drop = FALSE]
      candidate_links[[length(candidate_links) + 1L]] <- list(subject_id = g$subject_id,
        sample_id = g$sample_id, candidate_id = c$candidate_id,
        feature_key = longitudinal_candidate_key(c), ligand_sequence = c$ligand_sequence,
        mhc_key = c$mhc_key, hla_class = c$hla_class, origin_ids = c$origin_ids,
        sample_result_order = match(c$candidate_id, g$results$candidate_id),
        sample_result = lapply(row, function(col) if (length(col)) col[[1L]] else NULL), evidence = item[[1L]])
    }
    for (o in g$repertoire$observations) if (o$sample_id == g$sample_id)
      receptor_links[[length(receptor_links) + 1L]] <- list(subject_id = g$subject_id,
        sample_id = g$sample_id, feature_key = o$receptor_feature_profile_id, evidence = o,
        chain_observations = Filter(function(c) c$chain_observation_id %in% unlist(o$chain_observation_ids), g$repertoire$chains))
  }
  event_links <- longitudinal_sort(event_links); candidate_links <- longitudinal_sort(candidate_links)
  for (i in seq_along(candidate_links)) {
    o <- candidate_links[[i]]$evidence$origins
    candidate_links[[i]]$origin_equivalence <- lapply(o, longitudinal_origin_basis, event_links = event_links)
  }
  clonal <- longitudinal_clonal_import(project, settings, event_links)
  hla_contexts <- settings$contexts
  for (i in seq_along(hla_contexts)) {
    g <- Filter(function(g) g$sample_id == hla_contexts[[i]]$sample_id, graphs)[[1L]]
    hla_contexts[[i]]$observed_hla_genotype <- g$project$hla[[g$subject_id]]
  }
  trajectories <- function(links, domain) {
    out <- list()
    for (s in settings$subjects) {
      contexts <- Filter(function(c) c$subject_id == s$subject_id && !identical(c$role, "normal"), settings$contexts)
      order <- longitudinal_order(contexts); contexts <- order$contexts
      own <- Filter(function(l) l$subject_id == s$subject_id && !is.null(l$feature_key), links)
      keys <- sort(unique(vapply(own, `[[`, "", "feature_key")))
      for (key in keys) {
        cells <- lapply(contexts, function(c) {
          observed <- Filter(function(l) l$feature_key == key && l$sample_id == c$sample_id, own)
          cell <- longitudinal_cell(settings, c$sample_id, domain, key, observed)
          if (domain == "candidate") {
            cell$candidate_generated <- length(observed) > 0L
            source_keys <- unique(unlist(lapply(observed, function(l) lapply(l$origin_equivalence, `[[`, "required_event_keys"))))
            cell$source_event_states <- lapply(sort(source_keys), function(k) longitudinal_cell(settings,
              c$sample_id, "event", k, Filter(function(e) e$subject_id == s$subject_id && e$sample_id == c$sample_id &&
                identical(e$feature_key, k), event_links)))
          }
          event_ids <- if (domain == "event") vapply(observed, `[[`, "", "event_id") else if (domain == "candidate")
            vapply(Filter(function(e) e$sample_id == c$sample_id && e$subject_id == s$subject_id && !is.null(e$feature_key) &&
              e$feature_key %in% source_keys, event_links), `[[`, "", "event_id") else character()
          cell$clonal_evidence <- Filter(function(e) e$sample_id == c$sample_id && e$subject_id == s$subject_id &&
            e$mapping_status == "mapped" && any(unlist(e$event_ids) %in% event_ids), clonal)
          # Absence declarations apply only to their explicit evidence domain.
          cell$domain_evaluability <- lapply(c("rna", "ms", "repertoire", "hla_loh", "apm"), function(d) {
            evidence <- if (domain != "candidate") list() else unlist(lapply(observed, function(o) {
              view <- o$evidence$presentation_evidence
              switch(d, ms = Filter(function(l) l$status == "exact_unmodified_ligand_match", view$immunopeptidomics$links),
                rna = if (isTRUE(o$sample_result$rna_support)) list(o$evidence$somatic_evidence$event_evidence) else list(),
                hla_loh = if (is.null(o$evidence$hla_loh)) list() else list(o$evidence$hla_loh),
                apm = Filter(function(l) l$applicability == "applicable", view$antigen_presentation_machinery$links), list())
            }), recursive = FALSE)
            list(domain = d, cell = longitudinal_cell(settings, c$sample_id, d, key, evidence))
          })
          cell
        })
        summary <- longitudinal_summary(cells, contexts, settings$relations)
        if (domain == "repertoire") {
          summary$transitions <- list()
          summary$trajectory_state <- "sample_specific_receptor_observations"
        }
        record <- list(subject_id = s$subject_id, feature_key = key, domain = domain,
          ordering = order$ordering, sample_ids = as.list(vapply(contexts, `[[`, "", "sample_id")),
          cells = cells, summary = summary,
          escape = if (domain == "candidate") longitudinal_escape(cells, contexts, settings$relations) else NULL)
        record$trajectory_evidence_id <- longitudinal_trajectory_id(record, settings)
        out[[length(out) + 1L]] <- record
      }
    }
    out[order(vapply(out, `[[`, "", "trajectory_evidence_id"), method = "radix")]
  }
  x <- list(schema_version = longitudinal_version(), identity_namespaces = longitudinal_namespaces(),
    candidate_identity_basis = "frozen candidate_key", context_graph = settings,
    subject_context_ids = lapply(settings$subjects, function(s) list(subject_id = s$subject_id,
      subject_context_id = longitudinal_id("context", s))),
    sample_graph_refs = lapply(graphs, function(g) list(sample_id = g$sample_id, binding = g$binding)),
    event_links = event_links, candidate_links = candidate_links, clonal_evidence = clonal,
    event_trajectories = trajectories(event_links, "event"), candidate_trajectories = trajectories(candidate_links, "candidate"),
    receptor_trajectories = trajectories(longitudinal_sort(receptor_links), "repertoire"),
    hla_genotype = longitudinal_hla_consistency(hla_contexts, settings$subjects), external_clonal_runs = 0L)
  x$sample_matrix <- lapply(c(x$event_trajectories, x$candidate_trajectories), function(t)
    t[c("subject_id", "domain", "feature_key", "sample_ids", "cells")])
  x
}
longitudinal_fingerprint <- function(project, settings, graphs) {
  artifacts <- lapply(longitudinal_inputs(project), function(i) {
    if (!identical(integrity_record(i$resolved_path, project$config$checksum, i$integrity$declared), i$integrity)) stop("Longitudinal input changed since binding")
    recognition_artifact(i)
  })
  infra_hash(longitudinal_plain(list(schema = longitudinal_version(), namespaces = longitudinal_namespaces(), settings = settings,
    artifacts = longitudinal_sort(artifacts), graphs = lapply(graphs, function(g)
      g[setdiff(names(g), c("project", "binding"))]))))
}
longitudinal_path <- function(project) file.path(project$path, "results", "longitudinal")
longitudinal_postrun <- function(project, snapshot) {
  settings <- longitudinal_settings(project)
  if (is.null(settings)) {
    if (length(longitudinal_inputs(project))) stop("Clonal imports require explicit longitudinal settings")
    return(invisible(NULL))
  }
  graphs <- longitudinal_sample_graphs(project, settings)
  fp <- longitudinal_fingerprint(project, settings, graphs)
  root <- file.path(longitudinal_path(project), fp)
  if (file.exists(file.path(root, "checkpoint.json"))) return(invisible(longitudinal_read(project)))
  parent <- dirname(root); dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  if (!within_path(resolve_path(parent), project$path)) stop("Unsafe longitudinal output path")
  tmp <- tempfile(".longitudinal-", tmpdir = parent); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  x <- longitudinal_construct(project, settings, graphs)
  longitudinal_validate(x)
  saveRDS(x, file.path(tmp, "longitudinal.rds"), version = 3L)
  # Pair-record representation preserves frozen payloads with unnamed slots.
  atomic_json(longitudinal_plain(x), file.path(tmp, "longitudinal.json"))
  files <- list.files(tmp)
  atomic_json(list(schema_version = longitudinal_version(), fingerprint = fp,
    outputs = lapply(files, function(f) list(file = f, sha256 = digest::digest(file = file.path(tmp, f), algo = "sha256", serialize = FALSE)))), file.path(tmp, "checkpoint.json"))
  if (dir.exists(root) || !file.rename(tmp, root)) stop("Could not publish longitudinal overlay")
  invisible(longitudinal_read(project))
}
longitudinal_read <- function(project) {
  if (isTRUE(project$longitudinal_accessor_guard)) return(NULL)
  s <- longitudinal_settings(project)
  if (is.null(s)) return(NULL)
  g <- longitudinal_sample_graphs(project, s)
  fp <- longitudinal_fingerprint(project, s, g); root <- file.path(longitudinal_path(project), fp)
  path <- file.path(root, "checkpoint.json")
  if (!file.exists(path)) stop("Longitudinal sidecar missing; run neo_run")
  cp <- jsonlite::read_json(path)
  if (!identical(cp$fingerprint, fp) || !identical(cp$schema_version, longitudinal_version())) stop("Longitudinal fingerprint mismatch")
  for (o in cp$outputs) {
    p <- resolve_path(file.path(root, o$file))
    if (!within_path(p, root) || !file.exists(p) || !identical(o$sha256, digest::digest(file = p, algo = "sha256", serialize = FALSE))) stop("Longitudinal sidecar integrity failure")
  }
  x <- readRDS(file.path(root, "longitudinal.rds"))
  longitudinal_validate(x)
  x
}

longitudinal_trajectory_id <- function(record, settings) {
  semantic <- record
  semantic$trajectory_evidence_id <- NULL
  semantic$cells <- lapply(record$cells, function(c) {
    c$observations <- lapply(c$observations, function(o) {
      if (record$domain == "candidate") {
        row <- o$sample_result
        row[c("analysis_id", "analysis_context", "primary_prediction_ref", "origin_ids")] <- NULL
        return(list(candidate_id = o$candidate_id, sample_result = row, origin_equivalence = o$origin_equivalence))
      }
      if (record$domain == "event") return(list(event_id = o$event_id,
        canonical = o$feature_key, quantities = o$evidence[c("expression_json", "cn_by_breakpoint_json", "ccf_candidates_json")],
        quantitative_evidence = o$quantitative_evidence))
      list(observation = o$evidence[c("receptor_feature_profile_id", "sample_id", "support", "functional_status", "pairing_status")],
        abundance_fields = lapply(o$chain_observations, function(c) c$raw[c("duplicate_count", "umi_count", "consensus_count", "frequency", "relative_abundance")]))
    })
    c$clonal_evidence <- lapply(c$clonal_evidence, longitudinal_clonal_semantics)
    c
  })
  # Escape values are captured in sample_result; detailed provenance stays in the record.
  semantic$escape <- NULL
  longitudinal_id("trajectory", list(record = semantic,
    settings = settings[setdiff(names(settings), "sample_runs")], schema = longitudinal_version()))
}
