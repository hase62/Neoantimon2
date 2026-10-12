presentation_candidates <- function(project, snapshot, catalog = presentation_source_catalog(project, snapshot)) {
  candidates <- recognition_candidates(project, snapshot)
  lapply(candidates, function(c) {
    # A shared frozen candidate may be present in several source sidecars.
    own <- Filter(function(o) identical(o$sample_id, c$sample_id) &&
      identical(o$sample_peptide_id, c$sample_peptide_id) && identical(o$ligand_sequence, c$ligand_sequence), catalog$origins)
    sample <- Filter(function(s) identical(s$sample_id, c$sample_id), project$samples)
    if (length(sample) == 1L) for (field in c("timepoint", "tissue")) c[field] <- sample[[1L]][field]
    c$origin_ids <- as.list(sort(unique(c(unlist(c$origin_ids), vapply(own, `[[`, "", "origin_id")))))
    c
  })
}
# Read authoritative translated sources; never reconstruct flanks from another origin.
presentation_source_catalog <- function(project, snapshot) {
  origins <- snapshot$peptide_origins; proteins <- list(); boundaries <- list()
  path <- translation_stage_path(project)
  if (file.exists(path)) {
    cp <- Filter(function(x) x$step_id == "translation_resolution", read_manifest(project$path)$checkpoints)
    if (length(cp) != 1L || cp[[1L]]$status != "COMPLETED") stop("Translation checkpoint missing")
    verify_checkpoint_outputs(cp[[1L]], project$path)
    for (context in readRDS(path)) {
      t <- context$translation
      proteins[[t$translation_id]] <- t$mutant_protein
      boundaries[[t$translation_id]] <- list(protein_origin_id = if (is.null(t$mutant_protein)) NULL else
        paste0("protein-sha256:", peptide_digest(t$mutant_protein)), source_record_id = t$translation_id,
        n_terminus_known = identical(t$start_effect, "intact"),
        c_terminus_known = !is.null(t$mutant_stop_aa), translation_boundary = t$mutant_boundary)
    }
  }
  for (extension in list(source_read(project), long_read_read(project))) if (!is.null(extension)) {
    origins <- c(origins, extension$origins)
    for (h in extension$hypotheses) if (!is.null(h$orf$protein_origin_id)) {
      orf <- h$orf
      proteins[[orf$orf_hypothesis_id]] <- orf$protein
      boundaries[[orf$orf_hypothesis_id]] <- list(protein_origin_id = orf$protein_origin_id,
        n_terminus_known = TRUE, c_terminus_known = !is.null(orf$stop_nt), translation_boundary = orf$termination)
    }
  }
  list(origins = recognition_unique(origins, "origin_id"), proteins = proteins, boundaries = boundaries)
}
presentation_fingerprints <- function(project, snapshot) {
  catalog <- presentation_source_catalog(project, snapshot)
  checkpoint_fingerprints(inputs = list(artifacts = sorted_records(presentation_artifacts(project)),
      candidates = presentation_candidates(project, snapshot, catalog), origins = catalog$origins, proteins = catalog$proteins, boundaries = catalog$boundaries),
    config = list(schema = presentation_version(), identities = presentation_id_versions(),
      settings = presentation_settings(project), component_map = apm_component_map(), parser_version = "presentation-parser-v1"),
    dependencies = list(capabilities = lapply(presentation_settings(project)$providers, processing_capability,
      base = project$config$paths$input_root)))
}
presentation_construct <- function(project, candidates, catalog) {
  settings <- presentation_settings(project); inputs <- presentation_inputs(project)
  artifacts <- presentation_artifacts(project)
  queries <- list(); links <- list(); evidence <- list(); statuses <- list(); runs <- list()
  for (cfg in settings$providers) {
    cap <- processing_capability(cfg, project$config$paths$input_root)
    preparation <- processing_prepare(candidates, catalog, cap, settings$max_queries)
    prepared <- preparation$queries
    requested <- preparation$requested_queries; limited <- preparation$limited
    links <- c(links, preparation$links)
    output <- NULL
    if (!is.null(cfg$output_input_id)) {
      found <- Filter(function(i) i$input_id == cfg$output_input_id && i$input_role == "processing_output", inputs)
      if (length(found) != 1L) stop("Unresolved processing output binding")
      output <- found[[1L]]
    }
    result <- processing_predict(prepared, cap, output)
    queries <- c(queries, prepared); evidence <- c(evidence, result$evidence); statuses <- c(statuses, result$statuses)
    runs[[length(runs) + 1L]] <- list(capability = cap, requested_queries = requested,
      limit = settings$max_queries, status = if (limited) "not_run/query_limit" else "prepared",
      external_runs = 0L, local_execution = FALSE)
  }
  psms <- list(); apm <- list()
  for (i in inputs) {
    if (i$input_role %in% c("immunopeptidomics_psm", "immunopeptidomics_mzidentml")) {
      exp <- ms_experiment(i, settings)
      psms <- c(psms, if (i$input_role == "immunopeptidomics_psm") ms_generic_import(i, exp) else ms_mzidentml_import(i, exp))
    }
    if (i$input_role %in% c("apm_evidence", "hla_expression", "arcas_hla_quant")) apm <- c(apm, apm_import(i))
  }
  ms <- ms_observations(psms)
  evidence <- recognition_unique(evidence, "processing_evidence_id")
  for (i in seq_along(evidence)) evidence[[i]]$run_id <- project$run_id
  apm <- recognition_unique(apm, "apm_observation_id")
  list(schema_version = presentation_version(), identity_versions = presentation_id_versions(),
    immunopeptidomics_schema = "immunopeptidomics-evidence-1.0.0", component_map_version = apm_component_version(),
    run_id = project$run_id, artifacts = artifacts, queries = recognition_unique(queries, "processing_query_id"),
    processing_links = sorted_records(unique(links)), processing_evidence = evidence, psms = ms$psms, observations = ms$observations,
    ms_links = ms_candidate_links(candidates, ms$observations), apm_observations = apm,
    apm_links = apm_candidate_links(candidates, apm), statuses = sorted_records(statuses), provider_runs = runs,
    candidate_contexts = lapply(candidates, function(c) {
      fields <- c("candidate_id", "peptide_key", "sample_peptide_id", "sample_id", "ligand_sequence", "mhc_key", "hla_class", "origin_ids", "timepoint", "tissue")
      stats::setNames(lapply(fields, function(f) c[[f]]), fields)
    }),
    candidate_ids = as.list(vapply(candidates, `[[`, "", "candidate_id")))
}
presentation_validate <- function(x) {
  contract <- presentation_contract()
  if (!identical(contract$extends, peptide_contract()$peptide_version)) stop("Presentation overlay base mismatch")
  if (!identical(x$schema_version, presentation_version()) || !identical(x$identity_versions, presentation_id_versions()) ||
      !identical(x$component_map_version, apm_component_version()) || !all(unlist(contract$required) %in% names(x)))
    stop("Presentation schema mismatch")
  for (table in names(contract$records)) for (row in x[[table]])
    if (!all(unlist(contract$records[[table]]) %in% names(row))) stop("Incomplete presentation record: ", table)
  check_ids <- function(rows, field, kind, exclude = character()) {
    ids <- vapply(rows, `[[`, "", field)
    if (anyDuplicated(ids)) stop("Duplicate presentation ID")
    for (r in rows) {
      payload <- r[setdiff(names(r), c(field, exclude))]
      if (!identical(r$fingerprint, infra_hash(payload[setdiff(names(payload), "fingerprint")])))
        stop("Presentation record fingerprint mismatch")
      if (!identical(r[[field]], presentation_id(kind, payload))) stop("Presentation identity mismatch: ", field)
    }
  }
  for (q in x$queries) {
    payload <- q[setdiff(names(q), c("processing_query_id", "fingerprint", "input_status"))]
    if (!identical(q$processing_query_id, presentation_id("query", processing_query_identity(payload))) ||
        !identical(q$fingerprint, infra_hash(payload))) stop("Processing query identity mismatch")
  }
  if (anyDuplicated(vapply(x$queries, `[[`, "", "processing_query_id"))) stop("Duplicate processing query")
  link_keys <- vapply(x$processing_links, function(l) paste(l$processing_query_id, l$candidate_id, sep = "|"), "")
  if (anyDuplicated(link_keys) || !identical(x$processing_links, sorted_records(x$processing_links)))
    stop("Duplicate or unordered processing candidate links")
  for (l in x$processing_links) {
    q <- Filter(function(q) q$processing_query_id == l$processing_query_id, x$queries)
    candidate <- Filter(function(c) c$candidate_id == l$candidate_id, x$candidate_contexts)
    if (length(q) != 1L || length(candidate) != 1L) stop("Unresolved processing linkage")
    q <- q[[1L]]; candidate <- candidate[[1L]]
    applicable <- candidate$hla_class == "I"
    if (!identical(l$applicability, if (applicable) "applicable" else "not_applicable") ||
        !identical(l$status, if (applicable) "applicable" else "not_applicable_class_ii") ||
        !identical(q$peptide, candidate$ligand_sequence) ||
        (!is.null(q$presenter) && !identical(q$presenter, candidate$mhc_key)))
      stop("Unsafe processing candidate linkage")
    for (context in l$origin_contexts) {
      if (!context$origin_id %in% unlist(candidate$origin_ids)) stop("Unknown linked processing origin")
      if (q$task != "pmhc_stability") {
        payload <- q; payload$context <- context
        payload$n_flank <- if (q$options$use_n_flank) context$n_flank else NULL
        payload$c_flank <- if (q$options$use_c_flank) context$c_flank else NULL
        if (!identical(processing_query_identity(payload), processing_query_identity(q)))
          stop("Processing origin linkage mismatch")
      }
    }
  }
  if (anyDuplicated(vapply(x$processing_evidence, `[[`, "", "processing_evidence_id")))
    stop("Duplicate processing evidence ID")
  for (e in x$processing_evidence) {
    payload <- e[setdiff(names(e), c("processing_evidence_id", "fingerprint", "run_id"))]
    if (!identical(e$fingerprint, infra_hash(payload))) stop("Presentation record fingerprint mismatch")
    if (!identical(e$processing_evidence_id, presentation_id("evidence", processing_evidence_identity(e))))
      stop("Processing evidence identity mismatch")
  }
  check_ids(x$psms, "psm_id", "psm", "assignment_status")
  check_ids(x$observations, "immunopeptidomics_observation_id", "observation")
  check_ids(x$apm_observations, "apm_observation_id", "apm")
  expected_ms <- ms_observations(x$psms)
  if (!identical(expected_ms$psms, x$psms) || !identical(expected_ms$observations, x$observations))
    stop("MS observation derivation mismatch")
  if (!identical(ms_candidate_links(x$candidate_contexts, x$observations), x$ms_links) ||
      !identical(apm_candidate_links(x$candidate_contexts, x$apm_observations), x$apm_links))
    stop("Presentation candidate linkage mismatch")
  for (e in x$processing_evidence) {
    q <- Filter(function(q) q$processing_query_id == e$processing_query_id, x$queries)
    if (length(q) != 1L || q[[1L]]$input_status != "supported")
      stop("Ineligible processing evidence")
    presentation_number(e$raw_value)
  }
  for (o in x$observations) if (!all(unlist(o$psm_ids) %in% vapply(x$psms, `[[`, "", "psm_id"))) stop("Unresolved observation PSM")
  for (l in x$ms_links) {
    o <- Filter(function(o) o$immunopeptidomics_observation_id == l$observation_id, x$observations)
    if (length(o) != 1L || !l$candidate_id %in% unlist(x$candidate_ids)) stop("Unresolved MS linkage")
    if (l$status == "exact_unmodified_ligand_match" &&
        (o[[1L]]$modification_status != "unmodified" || o[[1L]]$assignment_status != "unique_assignment")) stop("Unsafe direct MS claim")
    if (l$direct_hla_assignment == "direct_allele_resolved" && o[[1L]]$restriction$status != "direct_allele_resolved")
      stop("Unsafe direct HLA claim")
  }
  for (l in x$apm_links) if (!l$apm_observation_id %in% vapply(x$apm_observations, `[[`, "", "apm_observation_id") ||
      !l$candidate_id %in% unlist(x$candidate_ids)) stop("Unresolved APM linkage")
  invisible(TRUE)
}
presentation_root <- function(project, fp) file.path(project$path, "results", "presentation", infra_hash(fp))
presentation_read <- function(project, snapshot = NULL) {
  if (!length(presentation_inputs(project))) return(NULL)
  if (is.null(snapshot)) snapshot <- read_runtime_snapshot(project)
  fp <- presentation_fingerprints(project, snapshot); root <- presentation_root(project, fp)
  path <- file.path(root, "checkpoint.json")
  if (!file.exists(path)) stop("Presentation sidecar missing; run neo_run")
  cp <- jsonlite::read_json(path)
  if (!identical(cp$fingerprints, fp) || !identical(cp$schema_version, presentation_version())) stop("Presentation fingerprint mismatch")
  for (f in cp$outputs) {
    file <- resolve_path(file.path(root, f$file))
    if (!within_path(file, root) || !file.exists(file) || !identical(f$sha256,
        digest::digest(file = file, algo = "sha256", serialize = FALSE))) stop("Presentation sidecar integrity failure")
  }
  x <- readRDS(file.path(root, "presentation.rds"))
  if (!identical(x$run_id, project$run_id)) stop("Presentation run mismatch")
  presentation_validate(x); x
}
presentation_postrun <- function(project, snapshot) {
  if (!length(presentation_inputs(project))) return(invisible(NULL))
  fp <- presentation_fingerprints(project, snapshot); root <- presentation_root(project, fp)
  if (file.exists(file.path(root, "checkpoint.json"))) return(invisible(presentation_read(project, snapshot)))
  parent <- dirname(root); dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  if (!within_path(resolve_path(parent), project$path)) stop("Unsafe presentation path")
  tmp <- tempfile(".presentation-", tmpdir = parent); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  x <- presentation_construct(project, presentation_candidates(project, snapshot), presentation_source_catalog(project, snapshot))
  presentation_validate(x)
  saveRDS(x, file.path(tmp, "presentation.rds"), version = 3L)
  for (field in c("queries", "processing_links", "processing_evidence", "psms", "observations", "ms_links", "apm_observations", "apm_links", "statuses", "provider_runs"))
    atomic_json(x[[field]], file.path(tmp, paste0(field, ".json")))
  files <- list.files(tmp)
  outputs <- lapply(files, function(f) list(file = f, sha256 = digest::digest(file = file.path(tmp, f), algo = "sha256", serialize = FALSE)))
  atomic_json(list(schema_version = presentation_version(), fingerprints = fp, outputs = outputs), file.path(tmp, "checkpoint.json"))
  if (dir.exists(root) || !file.rename(tmp, root)) stop("Could not publish presentation sidecar")
  invisible(presentation_read(project, snapshot))
}
