# Frozen records are copied, never regraded or rewritten.
vaccine_projection <- function(project) {
  guarded <- project; guarded$vaccine_accessor_guard <- TRUE
  overlay <- longitudinal_read(guarded)
  if (!is.null(overlay)) {
    records <- lapply(overlay$candidate_links, function(l) list(result = l$sample_result,
      evidence = l$evidence, subject_id = l$subject_id, candidate_key = l$feature_key,
      trajectory = Filter(function(t) t$subject_id == l$subject_id && t$feature_key == l$feature_key,
        overlay$candidate_trajectories)[[1L]]))
    contexts <- overlay$context_graph$contexts; relations <- overlay$context_graph$relations
    run_paths <- unique(vapply(overlay$sample_graph_refs, function(r) r$binding$path, ""))
  } else {
    rows <- neo_results(guarded)
    items <- if (nrow(rows)) neo_explain(guarded, candidate_id = rows$candidate_id)$candidates else list()
    records <- lapply(items, function(i) list(result = i$result, evidence = i,
      subject_id = Filter(function(s) s$sample_id == i$result$sample_id, project$samples)[[1L]]$patient_id,
      candidate_key = i$result$candidate_key, trajectory = NULL))
    contexts <- lapply(project$samples, function(s) list(sample_id = s$sample_id, subject_id = s$patient_id, role = s$role))
    relations <- list()
    run_paths <- project$path
  }
  # Core evidence can exist for source candidates even when their frozen
  # NOT_GRADED explanation does not include a somatic priority block.
  cores <- unlist(lapply(run_paths, function(path) unlist(lapply(somatic_prediction_runs(list(path = path)),
    `[[`, "core_evidence"), recursive = FALSE)), recursive = FALSE)
  records <- lapply(records, function(r) {
    r$evidence$predictor_core_evidence <- unique(c(r$evidence$predictor_core_evidence,
      Filter(function(c) identical(c$candidate_id, r$result$candidate_id), cores)))
    r
  })
  vaccine_project_records(records, contexts, relations, project$config$vaccine_set$normalized)
}
vaccine_project_records <- function(records, contexts, relations, policy) {
  projections <- lapply(records, function(rec) {
    r <- rec$result; e <- rec$evidence; t <- rec$trajectory
    scalar_text(rec$candidate_key, "frozen candidate_key")
    p <- list(subject_id = rec$subject_id, candidate_id = r$candidate_id, candidate_key = rec$candidate_key,
      sample_id = r$sample_id, ligand_sequence = r$ligand_sequence, presenter = r$mhc_key,
      hla_class = r$hla_class, epitope_key = if (is_missing(r$epitope_key)) NULL else r$epitope_key,
      frozen_result = r, origins = e$origins, evidence = e, trajectory = t,
      clonal_evidence = if (is.null(t)) list() else Filter(function(c) c$mapping_status == "mapped" && c$sample_id == r$sample_id,
        unlist(lapply(t$cells, `[[`, "clonal_evidence"), recursive = FALSE)))
    cores <- Filter(function(c) !is_missing(c$epitope_key) && identical(c$epitope_key_status, "derived") && identical(c$status, "evaluated") &&
      identical(c$presenter, p$presenter) && identical(c$candidate_id, p$candidate_id), e$predictor_core_evidence)
    p$epitope_keys <- vaccine_strings(c(if (is_missing(r$epitope_key)) list() else list(r$epitope_key),
      lapply(cores, `[[`, "epitope_key")), "frozen epitope keys")
    p$eligibility <- vaccine_eligibility(p, policy)
    p$escape <- vaccine_escape(p)
    p
  })
  ids <- vapply(projections, `[[`, "", "candidate_id")
  for (id in unique(ids[duplicated(ids)])) {
    duplicates <- projections[ids == id]
    if (length(unique(vapply(duplicates, function(p) canonical_bytes(longitudinal_plain(p)), ""))) != 1L)
      stop("Conflicting duplicate frozen candidate identity")
  }
  projections <- projections[!duplicated(ids)]
  projections <- projections[order(vapply(projections, `[[`, "", "candidate_id"), method = "radix")]
  tumor <- Filter(function(c) !identical(c$role, "normal"), contexts)
  list(candidates = projections, contexts = longitudinal_sort(tumor), relations = longitudinal_sort(relations),
    subject_ids = as.list(sort(unique(vapply(contexts, `[[`, "", "subject_id")), method = "radix")))
}
vaccine_escape <- function(p) {
  e <- p$evidence; h <- e$hla_loh; out <- list(hla_loh = list(), apm_loss = list(), apm_abnormality = list())
  link <- h$candidate_link
  # Consume the frozen allele-specific link; paired homolog loss is not this allele's loss.
  if (p$hla_class == "I" && identical(link$presenter_hla_loh_call, "loss_supported") &&
      identical(h$presenter_allele, p$presenter)) out$hla_loh <- list(list(sample_id = p$sample_id,
        presenter = p$presenter, evidence_id = link$presenter_hla_loh_evidence_id, evidence = h))
  apm <- e$presentation_evidence$antigen_presentation_machinery
  for (l in apm$links) if (identical(l$applicability, "applicable")) {
    obs <- Filter(function(o) identical(o$apm_observation_id, l$apm_observation_id) &&
      identical(o$hla_class, p$hla_class) && identical(o$context$sample_id, p$sample_id), apm$observations)
    if (!length(obs)) next
    o <- obs[[1L]]
    chains <- strsplit(p$presenter, "~", fixed = TRUE)[[1L]]
    if ((!is.null(o$allele) && !o$allele %in% chains) ||
        (grepl("_hla$", o$role) && !o$gene %in% sub("\\*.*$", "", chains))) next
    item <- list(sample_id = p$sample_id, hla_class = p$hla_class, gene = o$gene,
      evidence_id = o$apm_observation_id, evidence = o)
    if (identical(l$state, "direct_loss_supported")) out$apm_loss <- c(out$apm_loss, list(item))
    if (identical(l$state, "abnormality_observed")) out$apm_abnormality <- c(out$apm_abnormality, list(item))
  }
  out
}
# Collapse only explicitly declared biological contexts/technical replicates.
# This creates coverage units, never transitions or competing trajectory states.
vaccine_context_units <- function(contexts, relations) {
  ids <- sort(vapply(contexts, `[[`, "", "sample_id"), method = "radix")
  edges <- lapply(Filter(function(r) r$type == "technical_replicate", relations), function(r) c(r$from, r$to))
  nodes <- sort(unique(c(ids, unlist(edges))), method = "radix")
  unit <- stats::setNames(nodes, nodes)
  labels <- unique(unlist(lapply(contexts, `[[`, "biological_context_id")))
  for (label in labels) edges <- c(edges, list(vapply(Filter(function(c) identical(c$biological_context_id, label), contexts), `[[`, "", "sample_id")))
  repeat {
    before <- unit
    for (edge in edges) {
      edge <- intersect(edge, nodes)
      if (!length(edge)) next
      old <- unique(unit[edge]); unit[unit %in% old] <- sort(old, method = "radix")[[1L]]
    }
    if (identical(before, unit)) break
  }
  # Non-tumor bridge nodes can establish transitive replicate membership but
  # never become a tumor coverage representative or add an observation.
  result <- unit[ids]
  for (component in unique(result)) result[result == component] <- sort(ids[unit[ids] == component], method = "radix")[[1L]]
  result
}
vaccine_atom_fields <- function() c("presenters", "contexts", "regions", "shared", "persistent", "clonal")
vaccine_atoms_empty <- function() stats::setNames(rep(list(list()), length(vaccine_atom_fields())), vaccine_atom_fields())
vaccine_pair_atom <- function(pair) as.character(jsonlite::toJSON(as.list(pair), auto_unbox = TRUE))
# Vaccine selection support only: frozen validity and T039 states remain untouched.
vaccine_observation_supported <- function(p) {
  if (is.null(p$trajectory)) return(TRUE)
  cells <- Filter(function(cell) identical(cell$sample_id, p$sample_id), p$trajectory$cells)
  length(cells) == 1L && identical(cells[[1L]]$state, "observed")
}
vaccine_atoms <- function(observations, context_units, contexts) {
  atoms <- vaccine_atoms_empty()
  supported <- Filter(vaccine_observation_supported, observations)
  if (!length(supported)) return(atoms)
  atoms$presenters <- as.list(unique(vapply(supported, `[[`, "", "presenter")))
  observed_samples <- vapply(supported, `[[`, "", "sample_id")
  # Only eligible frozen observations contribute; no unobserved cells are coverage.
  atoms$contexts <- as.list(unique(unname(context_units[observed_samples])))
  observed_contexts <- Filter(function(c) c$sample_id %in% observed_samples, contexts)
  # Replicate groups provide at most one region atom; conflicting raw region
  # labels remain evidence and never inflate supported spatial coverage.
  atoms$regions <- unlist(lapply(sort(unique(unname(context_units[observed_samples])), method = "radix"), function(unit) {
    regions <- sort(unique(as.character(unlist(lapply(Filter(function(c) context_units[[c$sample_id]] == unit, observed_contexts), `[[`, "region_id")))), method = "radix")
    if (length(regions)) list(regions[[1L]]) else list()
  }), recursive = FALSE)
  t <- supported[[1L]]$trajectory
  if (!is.null(t) && t$summary$trajectory_state != "insufficient_to_classify") {
    if (t$summary$spatial_state %in% c("observed_all_regions", "observed_subset_regions") && length(atoms$regions) > 1L)
      atoms$shared <- lapply(utils::combn(sort(unlist(atoms$regions)), 2L, simplify = FALSE), vaccine_pair_atom)
    if (identical(t$ordering, "explicit_time")) {
      for (tr in t$summary$transitions) if (tr$state == "observed_at_both_timepoints" &&
          all(c(tr$from, tr$to) %in% observed_samples) && context_units[[tr$from]] != context_units[[tr$to]])
        atoms$persistent <- c(atoms$persistent, list(vaccine_pair_atom(c(context_units[[tr$from]], context_units[[tr$to]]))))
    }
  }
  for (o in supported) for (e in o$clonal_evidence) {
    # Provider-run-local links are coverage units; numeric quantities are never objectives.
    atoms$clonal <- c(atoms$clonal, list(infra_hash(longitudinal_plain(list(
      subject = o$subject_id, context = context_units[[o$sample_id]], provider = e$cluster_scope,
      cluster = e$provider_cluster_id, mutation = if (is.null(e$provider_cluster_id)) e$provider_mutation_id else NULL)))))
  }
  lapply(atoms, vaccine_strings, label = "coverage atoms")
}
vaccine_units <- function(projection, subject) {
  contexts <- Filter(function(c) c$subject_id == subject, projection$contexts)
  context_units <- vaccine_context_units(contexts, projection$relations)
  own <- Filter(function(p) p$subject_id == subject && p$sample_id %in% names(context_units) &&
    p$eligibility$status == "eligible" && vaccine_observation_supported(p), projection$candidates)
  keys <- sort(unique(vapply(own, `[[`, "", "candidate_key")), method = "radix")
  lapply(keys, function(key) {
    obs <- Filter(function(p) p$candidate_key == key, own); first <- obs[[1L]]
    if (length(unique(vapply(obs, function(p) paste(p$ligand_sequence, p$presenter, p$hla_class), ""))) != 1L)
      stop("Frozen candidate key content conflict")
    ep <- sort(unique(as.character(unlist(lapply(obs, `[[`, "epitope_keys")))), method = "radix")
    groups <- lapply(ep, function(ep) vaccine_id("redundancy", list(type = "same_epitope_key", epitope_key = ep, presenter = first$presenter)))
    list(candidate_key = key, candidate_ids = as.list(sort(vapply(obs, `[[`, "", "candidate_id"), method = "radix")),
      ligand_sequence = first$ligand_sequence, presenter = first$presenter, hla_class = first$hla_class,
      observations = obs, strong_groups = groups, atoms = vaccine_atoms(obs, context_units, contexts),
      escape = lapply(names(first$escape), function(k) longitudinal_sort(unlist(lapply(obs, function(o) o$escape[[k]]), recursive = FALSE))))
  }) |> lapply(function(u) { names(u$escape) <- c("hla_loh", "apm_loss", "apm_abnormality"); u })
}
vaccine_source_ids <- function(o) sort(unique(as.character(unlist(o[c("source_event_ids", "context_event_ids", "source_antigen_key")]))), method = "radix")
vaccine_redundancy <- function(projection, subject) {
  own <- Filter(function(p) p$subject_id == subject, projection$candidates)
  out <- list()
  if (length(own) < 2L) return(out)
  region_basis <- lapply(own, function(p) unique(lapply(p$origins, function(o) list(
    source_event_ids = as.list(vaccine_source_ids(o)), transcript_id = o$transcript_id,
    peptide_start_aa = o$peptide_start_aa, peptide_end_aa = o$peptide_end_aa,
    coordinate_context = if (!is.null(o$sequence_context_id)) o$sequence_context_id else
      if (!is.null(o$source_protein_id)) o$source_protein_id else o$origin_identity$reconstructed_sequence_checksum))))
  seen <- new.env(parent = emptyenv()); group_cache <- new.env(parent = emptyenv())
  add <- function(a, b, type, strong, detail = NULL) {
    group <- if (type == "same_epitope_key" && strong) list(type = type, epitope_key = detail$epitope_key, presenter = a$presenter) else
      list(type = type, keys = as.list(sort(unique(c(a$candidate_key, b$candidate_key)), method = "radix")), detail = detail)
    candidates <- as.list(sort(c(a$candidate_id, b$candidate_id), method = "radix"))
    # Internal memoization only; persisted IDs still use the versioned canonical contract.
    cache_key <- digest::digest(group, algo = "sha256", serialize = TRUE)
    record_key <- paste(cache_key, paste(unlist(candidates), collapse = "~"), sep = "~")
    if (exists(record_key, seen, inherits = FALSE)) return(invisible(NULL))
    assign(record_key, TRUE, seen)
    if (!exists(cache_key, group_cache, inherits = FALSE)) assign(cache_key, vaccine_id("redundancy", group), group_cache)
    out[[length(out) + 1L]] <<- list(group_id = get(cache_key, group_cache, inherits = FALSE),
      type = type, strong = strong, candidate_ids = candidates,
      candidate_keys = as.list(sort(unique(c(a$candidate_key, b$candidate_key)), method = "radix")), detail = detail)
  }
  for (i in seq_len(length(own) - 1L)) for (j in seq.int(i + 1L, length(own))) {
    a <- own[[i]]; b <- own[[j]]
    if (a$candidate_key == b$candidate_key) add(a, b, "same_frozen_candidate_key", TRUE)
    if (a$ligand_sequence == b$ligand_sequence && a$presenter != b$presenter) add(a, b, "same_ligand_different_presenter", FALSE)
    for (ep in intersect(unlist(a$epitope_keys), unlist(b$epitope_keys))) add(a, b, "same_epitope_key", a$presenter == b$presenter,
      list(epitope_key = ep, distinct_presenters = a$presenter != b$presenter))
    for (oa in region_basis[[i]]) for (ob in region_basis[[j]]) {
      source <- intersect(vaccine_source_ids(oa), vaccine_source_ids(ob))
      if (!length(source)) next
      add(a, b, "same_source_event", FALSE, list(source_event_ids = as.list(source)))
      if (identical(oa$transcript_id, ob$transcript_id) && !is.null(oa$transcript_id) &&
          identical(oa$coordinate_context, ob$coordinate_context) &&
          all(vapply(list(oa$peptide_start_aa, oa$peptide_end_aa, ob$peptide_start_aa, ob$peptide_end_aa),
            function(v) is.numeric(v) && length(v) == 1L && !is.na(v), TRUE)) &&
          max(oa$peptide_start_aa, ob$peptide_start_aa) <= min(oa$peptide_end_aa, ob$peptide_end_aa))
        add(a, b, "source_region_overlap", FALSE, list(source_event_ids = as.list(source), transcript_id = oa$transcript_id,
          coordinate_context = oa$coordinate_context,
          intervals = longitudinal_sort(list(list(start = oa$peptide_start_aa, end = oa$peptide_end_aa),
            list(start = ob$peptide_start_aa, end = ob$peptide_end_aa)))))
    }
  }
  out[order(vapply(out, function(r) paste(r$group_id, r$type, paste(unlist(r$candidate_ids), collapse = "~")), ""), method = "radix")]
}
