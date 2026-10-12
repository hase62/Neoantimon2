peptide_blank <- function(table, contract = peptide_contract()) {
  fields <- contract$tables[[table]]$fields
  stats::setNames(lapply(fields, function(f) if (f$type == "array") list() else
    if (f$type == "object") structure(list(), names = character()) else NULL),
    vapply(fields, `[[`, "", "name"))
}

peptide_binding <- function(project, sample_id) {
  ct <- Filter(function(x) sample_id %in% unlist(x$sample_bindings), project$analysis_context)
  if (length(ct) != 1L) stop("Ambiguous peptide analysis context")
  ct[[1L]]
}

peptide_origin <- function(p, context, project, events, contract = peptide_contract()) {
  ct <- peptide_binding(project, context$sample_id)
  by_id <- stats::setNames(events, vapply(events, `[[`, "", "event_id"))
  keys <- function(ids) {
    if (any(!ids %in% names(by_id))) stop("Peptide contributing event is missing")
    as.list(sort(unique(vapply(by_id[ids], `[[`, "", "event_key")), method = "radix"))
  }
  somatic <- keys(context$somatic_event_ids); germline <- keys(context$germline_event_ids)
  phase <- lapply(context$phase_context, function(x) list(event_key = by_id[[x$event_id]]$event_key,
    relation = x$relation))
  haplotype <- paste0("haplotype:", infra_hash(list(somatic = somatic, germline = germline,
    phase = sorted_records(phase), phase_state = context$phase_state,
    normal_protein = peptide_digest(context$translation$patient_normal_protein))))
  tx <- context$transcript_id
  version <- if (grepl("\\.[0-9]+$", tx)) sub(".*\\.", "", tx) else NULL
  identity <- list(transcript_id = sub("\\.[0-9]+$", "", tx), transcript_version = version,
    reconstructed_sequence_checksum = peptide_digest(context$tumor_transcript),
    peptide_start_aa = p$peptide_start_aa, peptide_end_aa = p$peptide_end_aa,
    source_event_keys = somatic, phase_namespace = NULL, phase_set = NULL,
    haplotype_id = haplotype, source_subject_id = ct$source_subject_id,
    source_antigen_key = NULL, scenario_id = ct$scenario_id)
  key <- biological_id("peptide_key", list(ligand_sequence = p$peptide_sequence), contract)
  sample_key <- biological_id("sample_peptide_id", list(context_id = ct$context_id,
    patient_id = ct$self_subject_id, sample_id = context$sample_id, peptide_key = key), contract)
  id <- biological_id("origin_id", list(sample_peptide_id = sample_key, origin_identity = identity), contract)
  source <- context$translation$mutant_protein
  wt <- p$normal_correspondence
  sequence_ref <- paste0("intermediate/translation_resolution.rds#", context$translation$translation_id)
  values <- list(schema_version = contract$peptide_version, analysis_id = project$analysis_id,
    patient_id = ct$self_subject_id, sample_id = context$sample_id, origin_id = id,
    sample_peptide_id = sample_key, peptide_key = key, peptide_sequence = p$peptide_sequence,
    ligand_sequence = p$peptide_sequence, origin_type = ct$mode,
    sequence_id = paste0("protein-sha256:", p$protein_checksum),
    sequence_context_id = context$translation$translation_id,
    peptide_length = p$peptide_length, peptide_start_aa = p$peptide_start_aa,
    peptide_end_aa = p$peptide_end_aa, transcript_id = identity$transcript_id,
    transcript_version = version, protein_id = context$protein_id,
    mane_status = if (length(context$mane_status))
      paste(sort(unique(context$mane_status), method = "radix"), collapse = ";") else NULL,
    consequence = context$translation$protein_effect,
    source_event_ids = as.list(sort(context$somatic_event_ids, method = "radix")),
    germline_context_ids = as.list(sort(context$germline_event_ids, method = "radix")),
    mutation_positions_json = p$altered_intervals, junction_positions_json = p$junction_positions,
    wt_peptide_sequence = wt$sequence, wt_definition = "PATIENT_NORMAL_T009_MAPPING",
    wt_comparison_status = wt$status, normal_comparison_status = wt$sequence_comparison,
    upstream_flank = if (p$peptide_start_aa > 1L) substring(source,
      max(1L, p$peptide_start_aa - 5L), p$peptide_start_aa - 1L) else NULL,
    downstream_flank = if (p$peptide_end_aa < nchar(source)) substring(source,
      p$peptide_end_aa + 1L, min(nchar(source), p$peptide_end_aa + 5L)) else NULL,
    mutant_sequence_ref = sequence_ref, normal_sequence_ref = sequence_ref,
    haplotype_id = haplotype, phase_status = context$phase_state,
    sequence_certainty = "T009_DETERMINISTIC", translation_status = context$translation$status,
    isoform_source = "BOUND_REFERENCE", origin_validity = "VALID",
    origin_evidence_status = "NOT_EVALUATED",
    origin_evaluation_ref = paste0("intermediate/peptide_candidates.rds#", id),
    analysis_context = ct$context_id, source_protein_id = context$protein_id,
    origin_identity = identity)
  for (role in c("source", "self", "presenter")) for (suffix in c("subject_id", "subject_role")) {
    name <- paste(role, suffix, sep = "_"); values[name] <- ct[name]
  }
  out <- peptide_blank("peptide_origins", contract); out[names(values)] <- values
  out
}

peptide_candidate <- function(origins, presenter, project, self, contract = peptide_contract()) {
  first <- origins[[1L]]
  key <- biological_id("candidate_key", list(species = project$config$scientific$species,
    hla_class = presenter$hla_class, ligand_sequence = first$ligand_sequence,
    mhc_key = presenter$mhc_key), contract)
  id <- biological_id("candidate_id", list(context_id = first$analysis_context,
    sample_id = first$sample_id, candidate_key = key), contract)
  values <- list(schema_version = contract$peptide_version, analysis_id = project$analysis_id,
    patient_id = first$patient_id, sample_id = first$sample_id, candidate_id = id,
    candidate_key = key, sample_peptide_id = first$sample_peptide_id,
    peptide_key = first$peptide_key, peptide_sequence = first$peptide_sequence,
    ligand_sequence = first$ligand_sequence, peptide_length = first$peptide_length,
    n_origins = as.integer(length(origins)),
    origin_ids = as.list(sort(vapply(origins, `[[`, "", "origin_id"), method = "radix")),
    mhc_key = presenter$mhc_key, hla_class = presenter$hla_class,
    hla_allele = presenter$hla_allele, hla_alpha = presenter$hla_alpha,
    hla_beta = presenter$hla_beta, hla_raw_input = presenter$raw_input,
    hla_pairing_status = presenter$pairing_status, species = project$config$scientific$species,
    prediction_status = "NOT_RUN", grading_status = "NOT_GRADED", final_grade = "NOT_GRADED",
    grade_before_warnings = "NOT_GRADED", not_evaluated_codes = list("NOT_RUN"),
    default_include = FALSE, analysis_context = first$analysis_context,
    self_status = if (self$status == "EVALUATED") "EVALUATED" else "NOT_EVALUABLE",
    self_exact_match = self$exact_match,
    self_match_protein_count = if (self$status == "EVALUATED")
      as.integer(length(unique(vapply(self$matches, `[[`, "", "protein_id")))) else NULL,
    self_match_details_ref = paste0("intermediate/peptide_candidates.rds#", first$peptide_key))
  for (role in c("source", "self", "presenter")) for (suffix in c("subject_id", "subject_role")) {
    name <- paste(role, suffix, sep = "_"); values[name] <- first[name]
  }
  out <- peptide_blank("neoantigen_results", contract); out[names(values)] <- values
  out$evidence_records_json <- list(peptide_self_evidence(out, project, self))
  out
}

peptide_self_evidence <- function(candidate, project, self) {
  resource <- Filter(function(x) x$logical_name == "reference.protein_fasta", project$resources)
  evaluated <- self$status == "EVALUATED"
  detail <- if (evaluated) list(self_scope = "HUMAN_REFERENCE_PROTEOME",
    self_subject_id = candidate$self_subject_id, query_identity = candidate$peptide_key,
    query_sequence = candidate$ligand_sequence, self_policy_id = "literal-exact-all-isoforms-v1",
    reference_release = if (length(resource) && !is.null(resource[[1L]]$version))
      resource[[1L]]$version else "DECLARED_REFERENCE_SEE_ARTIFACT",
    reference_checksum = self$reference_checksum, presenter_scope = "INDEPENDENT_OF_PRESENTER",
    binding_filter = list(applied = FALSE), search_coverage = list(all_supplied_isoforms = TRUE),
    matches = self$matches) else NULL
  list(evidence_id = paste0("self-evidence:", infra_hash(list(candidate = candidate$candidate_id,
      checksum = self$reference_checksum, status = self$status))),
    analysis_context = candidate$analysis_context, target_type = "candidate",
    target_id = candidate$candidate_id, evidence_type = "whole_proteome_self_match",
    assessment_status = if (evaluated) "EVALUATED" else "NOT_EVALUABLE",
    finding = if (!evaluated) NULL else if (self$exact_match) "POSITIVE" else "NEGATIVE",
    value = self$exact_match, unit = if (evaluated) "boolean" else NULL,
    metric_id = if (evaluated) "literal_exact_sequence_match" else NULL,
    reason = if (evaluated) NULL else self$status, method = "LITERAL_EXACT_ALL_ISOFORMS",
    source_ref = candidate$self_match_details_ref, source_checksum = self$reference_checksum,
    model_version = NULL, assay = NULL, tissue = NULL, timepoint = NULL,
    applicability = "APPLICABLE", uncertainty = list("REFERENCE_NOT_PERSONALIZED"),
    details = list(scope = "bound human reference protein catalog; no tolerance inference",
      self_match = detail))
}

peptide_construct <- function(project, contexts, events, self_path = NULL,
                              self_identity = NULL, sec_encoding = "UNDECLARED") {
  contract <- peptide_contract(); policy <- peptide_policy(project$config$peptides)
  origins <- list(); peptides <- list(); statuses <- list()
  if (length(contexts)) contexts <- contexts[order(vapply(contexts, `[[`, "", "reconstruction_id"), method = "radix")]
  for (context in contexts) {
    generated <- peptide_enumerate(context, policy)
    statuses[[length(statuses) + 1L]] <- generated[setdiff(names(generated), "peptides")]
    for (p in generated$peptides) {
      origin <- peptide_origin(p, context, project, events, contract)
      p$peptide_id <- origin$origin_id; p$sequence_key <- origin$peptide_key
      p$provenance_occurrences <- list(p$provenance)
      p$sample_peptide_id <- origin$sample_peptide_id
      p$sample_id <- origin$sample_id; p$analysis_context <- origin$analysis_context
      p$context_type <- peptide_binding(project, origin$sample_id)$mode
      p$transcript_id <- context$transcript_id; p$source_protein_id <- context$protein_id
      p$source_event_ids <- origin$source_event_ids; p$germline_context_ids <- origin$germline_context_ids
      for (role in c("source", "self", "presenter")) for (suffix in c("subject_id", "subject_role")) {
        name <- paste(role, suffix, sep = "_"); p[name] <- origin[name]
      }
      origins[[length(origins) + 1L]] <- origin; peptides[[length(peptides) + 1L]] <- p
    }
  }
  if (length(origins)) {
    ids <- vapply(origins, `[[`, "", "origin_id")
    for (i in which(duplicated(ids))) {
      first <- match(ids[[i]], ids)
      peptides[[first]]$provenance_occurrences <- c(peptides[[first]]$provenance_occurrences,
        peptides[[i]]$provenance_occurrences)
    }
    keep <- !duplicated(ids)
    origins <- origins[keep]; peptides <- peptides[keep]
    ordering <- order(ids[keep], method = "radix")
    origins <- origins[ordering]; peptides <- peptides[ordering]
  }
  self <- peptide_self_search(vapply(origins, `[[`, "", "ligand_sequence"),
    self_path, self_identity, sec_encoding)
  registry <- lapply(project$hla, peptide_presenters)
  candidates <- list(); requests <- list()
  groups <- split(origins, vapply(origins, `[[`, "", "sample_peptide_id"))
  for (group in groups) {
    ct <- peptide_binding(project, group[[1L]]$sample_id)
    presenters <- Filter(function(x) x$status == "RESOLVED", registry[[ct$presenter_subject_id]])
    if (length(presenters)) presenters <- presenters[!duplicated(vapply(presenters, `[[`, "", "mhc_key"))]
    for (presenter in presenters) {
      lengths <- policy[[if (presenter$hla_class == "I") "class_i_lengths" else "class_ii_lengths"]]
      if (!group[[1L]]$peptide_length %in% unlist(lengths)) next
      candidate <- peptide_candidate(group, presenter, project, self[[group[[1L]]$ligand_sequence]], contract)
      candidates[[length(candidates) + 1L]] <- candidate
      requests[[length(requests) + 1L]] <- list(candidate_id = candidate$candidate_id,
        sample_peptide_id = candidate$sample_peptide_id, biological_peptide_ids = candidate$origin_ids,
        peptide_sequence = candidate$ligand_sequence, presenter = candidate$mhc_key,
        hla_class = candidate$hla_class, supported_residue_check = "ADAPTER_PENDING",
        task = NULL, predictor = NULL, binding_core = NULL, binding_core_start = NULL,
        epitope_key = NULL, schema_version = contract$peptide_version)
    }
  }
  if (length(candidates)) {
    ordering <- order(vapply(candidates, `[[`, "", "candidate_id"), method = "radix")
    candidates <- candidates[ordering]; requests <- requests[ordering]
  }
  list(schema_version = contract$peptide_version, policy = policy,
    statuses = statuses, peptides = peptides, origins = origins, candidates = candidates,
    self_matches = self, presenter_registry = registry, predictor_inputs = requests,
    task_vocabulary = as.list(peptide_prediction_tasks()), prediction_evidence = peptide_prediction_evidence())
}
