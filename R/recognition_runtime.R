recognition_candidates <- function(project, snapshot) {
  priority <- somatic_priority_read(project, snapshot)
  candidates <- if (is.null(priority)) snapshot$neoantigen_results else priority$rows
  for (extension in list(source_read(project), long_read_read(project)))
    if (!is.null(extension)) candidates <- c(candidates, extension$candidates)
  if (length(candidates)) candidates <- candidates[!duplicated(vapply(candidates, `[[`, "", "candidate_id"))]
  recognition_sort(candidates, "candidate_id")
}
recognition_fingerprints <- function(project, candidates) checkpoint_fingerprints(
  inputs = list(artifacts = sorted_records(recognition_artifacts(project)),
    candidates = lapply(candidates, function(c) c[c("candidate_id", "sample_id", "ligand_sequence", "mhc_key", "hla_class")])),
  config = list(contract = recognition_version(), identities = recognition_id_versions(),
    settings = recognition_settings(project), parser = "recognition-provider-parser-v1"),
  dependencies = list(providers = lapply(recognition_settings(project)$providers, recognition_capability,
    base = project$config$paths$input_root)))
recognition_construct <- function(project, candidates) {
  settings <- recognition_settings(project)
  repertoire <- tcr_import(project, settings$max_pairing_hypotheses)
  queries <- list(); predictions <- list(); statuses <- list(); runs <- list()
  for (config in settings$providers) {
    cap <- recognition_capability(config, project$config$paths$input_root)
    prepared <- recognition_prepare(candidates, repertoire, cap, settings)
    output <- NULL
    if (!is.null(config$output_input_id)) {
      found <- Filter(function(x) x$input_id == config$output_input_id && x$input_role == "recognition_output",
        recognition_inputs(project))
      if (length(found) != 1L) stop("Unresolved recognition output binding")
      output <- found[[1L]]
    }
    prediction <- recognition_predict(prepared, cap, output)
    queries <- c(queries, prepared$queries); predictions <- c(predictions, prediction$evidence)
    statuses <- c(statuses, prediction$statuses)
    runs[[length(runs) + 1L]] <- list(capability = cap, status = prepared$status,
      reason = prepared$reason, requested_queries = prepared$requested, limit = settings$max_queries,
      candidate_ids = if (is.null(prepared$candidate_ids)) as.list(sort(vapply(candidates, `[[`, "", "candidate_id"))) else prepared$candidate_ids,
      query_ids = as.list(vapply(prepared$queries, `[[`, "", "recognition_query_id")), executed_external = FALSE)
  }
  experiments <- recognition_experiments(candidates, repertoire)
  annotations <- unlist(lapply(Filter(function(x) x$input_role == "tcr_vdjdb", recognition_inputs(project)),
    recognition_vdjdb), recursive = FALSE)
  curated <- recognition_curated_matches(candidates, repertoire, annotations)
  evidence <- recognition_unique(c(predictions, experiments$evidence, curated$evidence), "recognition_evidence_id")
  for (i in seq_along(evidence)) evidence[[i]]$run_id <- project$run_id
  list(schema_version = recognition_version(), identity_versions = recognition_id_versions(),
    run_id = project$run_id, artifacts = recognition_artifacts(project), repertoire = repertoire,
    queries = recognition_unique(queries, "recognition_query_id"), evidence = evidence,
    measurements = experiments$measurements, annotations = recognition_sort(annotations, "annotation_id"),
    repertoire_matches = curated$observations, statuses = sorted_records(statuses), provider_runs = runs,
    candidate_ids = as.list(vapply(candidates, `[[`, "", "candidate_id")))
}
recognition_validate <- function(x) {
  if (!identical(x$schema_version, recognition_version()) || !identical(x$identity_versions, recognition_id_versions()))
    stop("Recognition schema/identity mismatch")
  contract <- jsonlite::read_json(system.file("schema", "tcr-recognition-1.0.0.json", package = "Neoantimon2"))
  if (!identical(contract$extends, peptide_contract()$peptide_version)) stop("Recognition overlay base mismatch")
  if (!all(unlist(contract$required) %in% names(x))) stop("Incomplete recognition sidecar")
  for (table in names(contract$records)) for (row in x[[table]])
    if (!all(unlist(contract$records[[table]]) %in% names(row))) stop("Incomplete recognition record: ", table)
  ids <- vapply(x$repertoire$profiles, `[[`, "", "receptor_feature_profile_id")
  if (anyDuplicated(ids)) stop("Duplicate receptor profile")
  for (p in x$repertoire$profiles) if (!identical(p$receptor_feature_profile_id,
    recognition_profile(p$features, p$completeness)$receptor_feature_profile_id)) stop("Profile identity mismatch")
  for (o in x$repertoire$observations) {
    if (!o$functional_status %in% c("functional", "explicitly_nonfunctional",
        "functionality_unknown", "functionality_conflicting")) stop("Invalid receptor functionality")
    if (!o$receptor_feature_profile_id %in% ids ||
      !all(unlist(o$chain_observation_ids) %in% vapply(x$repertoire$chains, `[[`, "", "chain_observation_id")))
      stop("Unresolved repertoire observation")
  }
  for (q in x$queries) {
    if (q$scope == "tcr_specific") {
      own <- Filter(function(o) o$sample_id == q$sample_id &&
        o$receptor_feature_profile_id == q$receptor_feature_profile_id, x$repertoire$observations)
      functionality <- recognition_functionality(own)
      if (!identical(q$receptor_functionality, functionality)) stop("Receptor functionality context mismatch")
      if (functionality$status %in% c("explicitly_nonfunctional", "functionality_conflicting") &&
          q$input_status == "supported") stop("Nonfunctional receptor prediction input")
    }
    fields <- c("candidate_id", "peptide", "receptor_feature_profile_id", "presenter", "provider",
      "provider_version", "model_version", "model", "scope", "features_used", "required_features",
      "mhc_conditioning", "options", "model_identity", "implementation_identity", "background_identity")
    payload <- q[fields]
    if (!identical(q$recognition_query_id, recognition_id("query", payload)) ||
      !identical(q$fingerprint, infra_hash(list(contract_version = recognition_version(), payload = payload))))
      stop("Recognition query identity mismatch")
    if (!q$candidate_id %in% unlist(x$candidate_ids) || (!is.null(q$receptor_feature_profile_id) &&
      !q$receptor_feature_profile_id %in% ids)) stop("Unresolved recognition query")
    if (q$mhc_conditioning == "absent" && !is.null(q$presenter)) stop("Unconditioned query presenter claim")
    if (q$mhc_conditioning != "absent" && q$input_status == "supported" &&
      recognition_presenter(q$presenter)$status != "known") stop("Invalid conditioned recognition presenter")
  }
  if (anyDuplicated(vapply(x$queries, `[[`, "", "recognition_query_id"))) stop("Duplicate recognition query")
  if (anyDuplicated(vapply(x$evidence, `[[`, "", "recognition_evidence_id"))) stop("Duplicate recognition evidence")
  for (e in x$evidence) {
    if (!e$scope %in% unlist(contract$scopes) || !e$score_direction %in% c("higher", "lower", "not_defined"))
      stop("Invalid recognition evidence semantics")
    scalar_text(e$metric, "recognition metric")
    for (field in c("raw_score", "rank", "probability")) if (!is.null(e[[field]])) recognition_number(e[[field]])
    if (!is.null(e$query_id) && !e$query_id %in% vapply(x$queries, `[[`, "", "recognition_query_id"))
      stop("Unresolved recognition evidence query")
    if (e$scope == "tcr_specific" && !is.null(e$query_id)) {
      query <- Filter(function(q) q$recognition_query_id == e$query_id, x$queries)[[1L]]
      if (query$input_status != "supported") stop("Prediction evidence for ineligible receptor query")
    }
    if (!e$candidate_id %in% unlist(x$candidate_ids)) stop("Unknown recognition candidate")
    if (!is.null(e$receptor_feature_profile_id) && !e$receptor_feature_profile_id %in% ids) stop("Unknown recognition profile")
    if (e$scope == "peptide_intrinsic" && !is.null(e$receptor_feature_profile_id)) stop("Intrinsic receptor mismatch")
    if (e$scope == "tcr_specific" && is.null(e$receptor_feature_profile_id)) stop("Missing TCR profile")
  }
  invisible(TRUE)
}
recognition_root <- function(project, fp) file.path(project$path, "results", "recognition", infra_hash(fp))
recognition_read <- function(project, snapshot = NULL) {
  if (!length(recognition_inputs(project))) return(NULL)
  if (is.null(snapshot)) snapshot <- read_runtime_snapshot(project)
  fp <- recognition_fingerprints(project, recognition_candidates(project, snapshot))
  root <- recognition_root(project, fp)
  path <- file.path(root, "checkpoint.json")
  if (!file.exists(path)) stop("Recognition sidecar missing; run neo_run")
  cp <- jsonlite::read_json(path)
  if (!identical(cp$fingerprints, fp) || !identical(cp$schema_version, recognition_version())) stop("Recognition fingerprint mismatch")
  for (f in cp$outputs) {
    file <- resolve_path(file.path(root, f$file))
    if (!within_path(file, root) || !file.exists(file) || !identical(f$sha256,
        digest::digest(file = file, algo = "sha256", serialize = FALSE))) stop("Recognition sidecar integrity failure")
  }
  x <- readRDS(file.path(root, "recognition.rds"))
  if (!identical(x$run_id, project$run_id)) stop("Recognition run identity mismatch")
  recognition_validate(x)
  x
}
recognition_postrun <- function(project, snapshot) {
  if (!length(recognition_inputs(project))) return(invisible(NULL))
  candidates <- recognition_candidates(project, snapshot)
  fp <- recognition_fingerprints(project, candidates)
  root <- recognition_root(project, fp)
  if (file.exists(file.path(root, "checkpoint.json"))) return(invisible(recognition_read(project, snapshot)))
  # Stage all files before publishing the sidecar. Interrupted writes can retry.
  parent <- dirname(root); dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  if (!within_path(resolve_path(parent), project$path)) stop("Unsafe recognition result path")
  tmp <- tempfile(".recognition-", tmpdir = parent); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  x <- recognition_construct(project, candidates); recognition_validate(x)
  saveRDS(x, file.path(tmp, "recognition.rds"), version = 3L)
  atomic_json(x$evidence, file.path(tmp, "evidence.json"))
  atomic_json(x$repertoire, file.path(tmp, "repertoire.json"))
  atomic_json(x$queries, file.path(tmp, "queries.json"))
  atomic_json(list(statuses = x$statuses, provider_runs = x$provider_runs), file.path(tmp, "status.json"))
  files <- c("recognition.rds", "evidence.json", "repertoire.json", "queries.json", "status.json")
  outputs <- lapply(files, function(f) list(file = f, sha256 = digest::digest(
    file = file.path(tmp, f), algo = "sha256", serialize = FALSE)))
  atomic_json(list(schema_version = recognition_version(), fingerprints = fp, outputs = outputs), file.path(tmp, "checkpoint.json"))
  if (dir.exists(root) || !file.rename(tmp, root)) stop("Could not publish recognition sidecar")
  invisible(recognition_read(project, snapshot))
}
recognition_candidate_view <- function(x, candidate_id) {
  evidence <- Filter(function(e) e$candidate_id == candidate_id, x$evidence)
  statuses <- Filter(function(s) s$candidate_id == candidate_id, x$statuses)
  families <- stats::setNames(lapply(c("peptide_intrinsic", "tcr_specific", "experimental_reactivity", "external_annotation"),
    function(scope) Filter(function(e) e$scope == scope, evidence)),
    c("peptide_intrinsic", "tcr_specific", "experimental_reactivity", "external_annotation"))
  matches <- Filter(function(m) m$candidate_id == candidate_id, x$repertoire_matches)
  queries <- Filter(function(q) q$candidate_id == candidate_id, x$queries)
  observed <- Filter(function(o) any(vapply(queries, function(q)
    identical(o$sample_id, q$sample_id) && identical(o$receptor_feature_profile_id,
      q$receptor_feature_profile_id), TRUE)), x$repertoire$observations)
  # Exact profile membership supplies observation context, never specificity.
  observations <- c(lapply(observed, function(o) list(candidate_id = candidate_id,
    receptor_feature_profile_id = o$receptor_feature_profile_id,
    repertoire_observation_id = o$repertoire_observation_id, status = "receptor_observed", observation = o)), matches)
  if (length(observations)) observations <- observations[!duplicated(vapply(observations, `[[`, "", "repertoire_observation_id"))]
  linked_ids <- unique(c(vapply(observations, `[[`, "", "receptor_feature_profile_id"),
    unlist(lapply(evidence, `[[`, "receptor_feature_profile_id"))))
  profiles <- Filter(function(p) p$receptor_feature_profile_id %in% linked_ids, x$repertoire$profiles)
  limits <- Filter(function(r) identical(r$reason, "query_limit") && candidate_id %in% unlist(r$candidate_ids), x$provider_runs)
  state <- if (length(evidence)) "evidence_present" else if (length(limits) || any(vapply(statuses,
    function(s) s$status == "not_run", TRUE))) "not_run" else if (length(statuses)) "unsupported_provider" else "no_evidence_supplied"
  summary <- list(recognition_evidence_status = state,
    peptide_intrinsic_evidence_count = length(families$peptide_intrinsic),
    tcr_specific_prediction_count = length(unique(vapply(families$tcr_specific, function(e)
      if (is.null(e$query_id)) e$recognition_evidence_id else e$query_id, ""))),
    experimental_reactivity_count = length(families$experimental_reactivity),
    curated_specificity_match_count = sum(vapply(families$external_annotation, function(e) e$provider == "vdjdb", TRUE)),
    repertoire_receptor_match_count = length(unique(vapply(observations, `[[`, "", "repertoire_observation_id"))))
  c(list(summary = summary), families, list(receptor_feature_profiles = profiles, repertoire_observations = observations,
    exact_curated_repertoire_matches = matches, query_statuses = statuses,
    query_limits = limits, grade_policy_note = "Recognition evidence does not modify frozen grade, rank, eligibility or priority origin."))
}
recognition_results_append <- function(out, project, snapshot) {
  x <- recognition_read(project, snapshot)
  if (is.null(x)) return(out)
  views <- lapply(out$candidate_id, function(id) recognition_candidate_view(x, id)$summary)
  fields <- names(recognition_candidate_view(x, "")$summary)
  for (f in fields) out[[f]] <- vapply(views, `[[`, if (f == "recognition_evidence_status") "" else 0L, f)
  out
}
recognition_explain_append <- function(candidates, project, snapshot) {
  x <- recognition_read(project, snapshot)
  if (is.null(x)) return(candidates)
  lapply(candidates, function(c) {
    c$recognition <- recognition_candidate_view(x, c$result$candidate_id)
    c
  })
}
recognition_explain_text <- function(item) {
  x <- item$recognition
  if (is.null(x)) return("")
  family <- function(rows) if (!length(rows)) "no evidence supplied" else paste(vapply(rows, function(e) paste0(
    e$provider, " ", e$metric, if (is.null(e$rank)) "" else paste0(" rank=", e$rank),
    if (is.null(e$raw_score)) "" else paste0(" value=", e$raw_score),
    " [", e$status, "; MHC=", e$mhc_conditioning, "]"), ""), collapse = "; ")
  paste0("\nPRESENTATION: ", item$result$prediction_status,
    "\nPEPTIDE-INTRINSIC RECOGNITION EVIDENCE: ", family(x$peptide_intrinsic),
    "\nTCR-SPECIFIC PREDICTION: ", family(x$tcr_specific),
    "\nEXPERIMENTAL REACTIVITY: ", family(x$experimental_reactivity),
    "\nCURATED SPECIFICITY MATCHES: ", family(x$external_annotation),
    "\nREPERTOIRE OBSERVATION: ", length(x$repertoire_observations), " exact linked observations; presence is not recognition.",
    "\nRECOGNITION STATUS: ", x$summary$recognition_evidence_status)
}
