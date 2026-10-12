long_read_stage_path <- function(project)
  file.path(project$path, "intermediate", "long_read_extension.rds")
long_read_observations <- function(project, bundle) {
  out <- list()
  for (b in long_read_bindings(project)) {
    models <- long_read_models(b$input$resolved_path, bundle$profile$assembly_release)
    metadata <- unlist(lapply(b$metadata, long_read_metadata), recursive = FALSE)
    fasta <- if (length(b$fasta)) long_read_fasta(b$fasta[[1L]]$resolved_path) else list()
    checksum <- digest::digest(file = b$input$resolved_path, algo = "sha256", serialize = FALSE)
    for (model in models) {
      id <- model$provider_transcript_id
      meta <- Filter(function(m) m$provider_transcript_id == id, metadata)
      s <- model$structure
      for (m in meta) if (b$provider == "sqanti3" && !is.null(s)) {
        r <- m$raw_fields
        numeric <- vapply(c(r$length, r$exons), long_read_integer_evidence, 0L)
        if (anyNA(numeric) || r$chrom != s$contig || r$strand != s$strand ||
            numeric[[1L]] != length(long_read_coordinates(s)) || numeric[[2L]] != length(s$exons))
          stop("SQANTI3 metadata/structure conflict")
      }
      ct <- peptide_binding(project, b$input$sample_id)
      structural_id <- if (is.null(s)) NULL else long_read_structure_id(s)
      event <- if (is.null(s)) NULL else long_read_id("long_read_source_event_id", list(
        structure = structural_id, context = b$context, analysis_context = ct$context_id))
      obs <- long_read_id("long_read_observation_id", list(input_id = b$input$input_id,
        checksum = checksum, provider_transcript_id = id, context = b$context,
        auxiliary_artifacts = lapply(c(b$metadata, b$fasta), function(x) list(input_id = x$input_id,
          sha256 = digest::digest(file = x$resolved_path, algo = "sha256", serialize = FALSE)))))
      out[[length(out) + 1L]] <- list(long_read_observation_id = obs,
        isoform_structure_id = structural_id, long_read_source_event_id = event,
        provider = b$provider, provider_version = b$input$version,
        provider_transcript_id = id, provider_gene_ids = model$provider_gene_ids,
        context = b$context, patient_id = ct$self_subject_id,
        sequencing_technology = b$sequencing_technology,
        matched_normal_status = b$matched_normal_status, de_novo_aug = b$de_novo_aug,
        reference_build = bundle$profile$assembly_release, reference_identity = bundle$resource_identity,
        parser_version = "long-read-provider-parser-v1", model = model, metadata = meta,
        artifact = list(input_id = b$input$input_id, path = b$input$declared_path, sha256 = checksum),
        provider_fasta_artifact = if (!length(b$fasta)) NULL else list(
          input_id = b$fasta[[1L]]$input_id, path = b$fasta[[1L]]$declared_path,
          sha256 = digest::digest(file = b$fasta[[1L]]$resolved_path, algo = "sha256", serialize = FALSE)),
        provider_sequence = fasta[[id]], fasta_supplied = length(b$fasta) > 0L)
    }
    unknown <- setdiff(vapply(metadata, `[[`, "", "provider_transcript_id"),
      vapply(models, `[[`, "", "provider_transcript_id"))
    if (length(unknown)) stop("Long-read metadata references absent structural model: ", unknown[[1L]])
  }
  if (length(out)) out <- out[order(vapply(out, `[[`, "", "long_read_observation_id"))]
  out
}

long_read_extension_construct <- function(project, bundle = reference_bundle(project, retain_indexes = TRUE)) {
  observations <- long_read_observations(project, bundle)
  parsed <- Filter(function(o) !is.null(o$model$structure), observations)
  structures <- lapply(split(parsed, vapply(parsed, `[[`, "", "isoform_structure_id")),
    function(g) g[[1L]]$model$structure)
  reference <- event_reference_open(bundle); on.exit(event_reference_close(reference))
  refs <- long_read_reference_models(bundle, structures, reference)
  refs <- Filter(Negate(is.null), refs)
  reconstructed <- lapply(structures, long_read_reconstruct, reference = reference)
  comparisons <- lapply(names(structures), function(id)
    long_read_compare(structures[[id]], reconstructed[[id]], refs))
  names(comparisons) <- names(structures)
  events <- list(); hypotheses <- list(); origins <- list(); peptides <- list()
  for (group in split(parsed, vapply(parsed, `[[`, "", "long_read_source_event_id"))) {
    o <- group[[1L]]; id <- o$isoform_structure_id; rec <- reconstructed[[id]]
    cmp <- comparisons[[id]]
    normal <- Filter(function(n) identical(n$isoform_structure_id, id) &&
      n$patient_id == o$patient_id && n$context$observation_role == "normal", observations)
    normal_states <- unique(vapply(group, `[[`, "", "matched_normal_status"))
    state <- if (length(normal)) "normal_observed" else if (length(normal_states) == 1L)
      normal_states[[1L]] else "unresolved"
    orfs <- long_read_orfs(id, rec, cmp, group)
    ids <- as.list(vapply(group, `[[`, "", "long_read_observation_id"))
    for (n in seq_along(group)) {
      group[[n]]$fasta_comparison <- long_read_fasta_comparison(rec$sequence,
        group[[n]]$provider_sequence, group[[n]]$fasta_supplied)
      hit <- which(vapply(observations, `[[`, "", "long_read_observation_id") == group[[n]]$long_read_observation_id)
      observations[[hit]] <- group[[n]]
    }
    record <- list(long_read_source_event_id = o$long_read_source_event_id,
      isoform_structure_id = id, context = o$context, observation_ids = ids,
      normal_observation_ids = as.list(vapply(normal, `[[`, "", "long_read_observation_id")),
      matched_normal_status = state, structural_comparison = cmp,
      altered_transcript = c(rec, list(long_read_source_event_id = o$long_read_source_event_id,
        isoform_structure_id = id, sample_context = o$context, reference_identity = bundle$resource_identity,
        junctions = cmp$junctions, reference_comparison = cmp$status,
        novel_nt_intervals = cmp$novel_nt_intervals,
        structural_novelty = lapply(cmp$comparisons, function(c) c[c("transcript_id", "classes")]))),
      orf_unavailable = orfs$unavailable)
    events[[length(events) + 1L]] <- record
    for (orf in orfs$hypotheses) {
      hypotheses[[length(hypotheses) + 1L]] <- list(long_read_source_event_id = o$long_read_source_event_id,
        isoform_structure_id = id, orf = orf, observation_ids = ids)
      if (o$context$observation_role != "tumor" || !isTRUE(orf$peptide_eligible)) next
      for (p in source_peptide_windows(orf, peptide_policy(project$config$peptides))) {
        # Frozen origin identity construction is reused with an isoform and ORF
        # identity in the transcript/haplotype slots; no new peptide namespace.
        event <- list(source_type = "long_read_isoform", source_event_id = o$long_read_source_event_id,
          structure = list(sample_id = o$context$sample_id))
        altered <- list(sequence = rec$sequence,
          source_models = list(reference = list(model = list(transcript_id = id))))
        origin <- source_origin(project, event, orf$orf_hypothesis_id, altered, orf, p)
        for (field in c("mutant_sequence_ref", "origin_evaluation_ref", "source_sequence_ref"))
          origin[[field]] <- sub("source_event_extension.rds", "long_read_extension.rds", origin[[field]], fixed = TRUE)
        origins[[length(origins) + 1L]] <- origin
        peptides[[length(peptides) + 1L]] <- list(origin_id = origin$origin_id,
          peptide_key = origin$peptide_key, long_read_source_event_id = o$long_read_source_event_id,
          isoform_structure_id = id, orf_hypothesis_id = orf$orf_hypothesis_id,
          evidence_type = orf$evidence_type, peptide = p)
      }
    }
  }
  self_resource <- Filter(function(x) x$logical_name == "reference.protein_fasta", project$resources)
  self <- peptide_self_search(vapply(origins, `[[`, "", "ligand_sequence"),
    if (length(self_resource)) self_resource[[1L]]$resolved_path else NULL,
    identity_records(self_resource), sec_encoding = "LITERAL_U")
  candidates <- list()
  for (g in split(origins, vapply(origins, `[[`, "", "sample_peptide_id"))) {
    ct <- peptide_binding(project, g[[1L]]$sample_id)
    for (presenter in Filter(function(p) p$status == "RESOLVED",
        peptide_presenters(project$hla[[ct$presenter_subject_id]]))) {
      lengths <- peptide_policy(project$config$peptides)[[if (presenter$hla_class == "I")
        "class_i_lengths" else "class_ii_lengths"]]
      if (!g[[1L]]$peptide_length %in% unlist(lengths)) next
      row <- peptide_candidate(g, presenter, project, self[[g[[1L]]$ligand_sequence]])
      row$not_evaluated_codes <- list("NOT_APPLICABLE")
      row$self_match_details_ref <- sub("peptide_candidates.rds", "long_read_extension.rds",
        row$self_match_details_ref, fixed = TRUE)
      candidates[[length(candidates) + 1L]] <- row
    }
  }
  list(schema_version = long_read_version(), identity_schema_version = long_read_identity_version(),
    artifacts = long_read_artifacts(project), observations = observations,
    structures = lapply(names(structures), function(id) list(isoform_structure_id = id, structure = structures[[id]])),
    events = events, hypotheses = hypotheses, origins = origins, peptides = peptides,
    candidates = candidates, self_matches = self)
}

long_read_reference_artifacts <- function(project) lapply(Filter(function(r)
  startsWith(r$logical_name, "reference.") && !is.null(r$resolved_path), project$resources), function(r) {
  list(role = r$logical_name, version = r$version, sha256 = digest::digest(
    file = r$resolved_path, algo = "sha256", serialize = FALSE))
})
long_read_canonical_input_ids <- function(project) sort(unique(unlist(lapply(
  long_read_bindings(project), function(b) vapply(Filter(function(m) !is.null(m$structure),
    long_read_models(b$input$resolved_path, project$config$scientific$genome_assembly)),
    function(m) long_read_structure_id(m$structure), "")))))
long_read_fingerprints <- function(project) checkpoint_fingerprints(
  inputs = list(artifacts = sorted_records(long_read_artifacts(project)),
    canonical_structure_ids = as.list(long_read_canonical_input_ids(project))),
  config = list(schema = long_read_version(), parser = "long-read-provider-parser-v1",
    comparison = "long-read-reference-comparison-v1", reconstruction = "long-read-reconstruction-v1",
    orf = list("long-read-reference-orf-v1", "long-read-provider-orf-v1", long_read_orf_policy()),
    settings = long_read_bindings(project), sec_context = "exact-reference-context-v1",
    peptides = peptide_policy(project$config$peptides), hla = project$hla,
    contexts = project$analysis_context, samples = project$samples),
  dependencies = list(reference = sorted_records(identity_records(project$resources)),
    reference_content = sorted_records(long_read_reference_artifacts(project))))
long_read_validate <- function(x, project) {
  if (!identical(x$schema_version, long_read_version()) ||
      !identical(x$identity_schema_version, long_read_identity_version()) ||
      !identical(x$artifacts, long_read_artifacts(project))) stop("Long-read extension identity/artifact mismatch")
  contract <- jsonlite::read_json(system.file("schema", "long-read-source-1.0.0.json", package = "Neoantimon2"))
  if (!identical(contract$extends, peptide_contract()$peptide_version) ||
      !all(unlist(contract$required) %in% names(x))) stop("Invalid long-read schema overlay")
  for (table in names(contract$records)) for (row in x[[table]])
    if (!all(unlist(contract$records[[table]]) %in% names(row))) stop("Missing long-read record field")
  for (pair in list(c("structures", "isoform_structure_id"), c("observations", "long_read_observation_id"),
      c("events", "long_read_source_event_id")))
    if (anyDuplicated(vapply(x[[pair[[1L]]]], `[[`, "", pair[[2L]]))) stop("Duplicate long-read record identity")
  for (e in x$events) {
    if (!all(unlist(e$observation_ids) %in% vapply(x$observations, `[[`, "", "long_read_observation_id")))
      stop("Unresolved long-read observation")
    if (!identical(e$altered_transcript$sequence_authority, "REFERENCE_GENOME") ||
        nchar(e$altered_transcript$sequence) != length(e$altered_transcript$coordinates))
      stop("Invalid long-read sequence authority/mapping")
  }
  for (s in x$structures) if (!identical(s$isoform_structure_id, long_read_structure_id(s$structure)))
    stop("Long-read structural identity mismatch")
  validate_table(x$origins, "peptide_origins", peptide_contract())
  validate_table(x$candidates, "neoantigen_results", peptide_contract())
  origins <- vapply(x$origins, `[[`, "", "origin_id")
  if (anyDuplicated(origins)) stop("Duplicate long-read origin")
  for (c in x$candidates) if (c$final_grade != "NOT_GRADED" ||
    !all(unlist(c$origin_ids) %in% origins)) stop("Long-read candidate grade/provenance mismatch")
  invisible(TRUE)
}
long_read_stage_run <- function(project) {
  if (!length(long_read_inputs(project))) return(NULL)
  fp <- long_read_fingerprints(project)
  step <- "long_read_extension"
  cp <- Filter(function(x) x$step_id == step, read_manifest(project$path)$checkpoints)
  if (length(cp) && cp[[1L]]$status == "COMPLETED" && checkpoint_reusable(project, step, fp))
    return(readRDS(long_read_stage_path(project)))
  if (!length(cp)) checkpoint_update(project, step, fingerprints = fp)
  checkpoint_update(project, step, "RUNNING", fp)
  tryCatch({
    result <- long_read_extension_construct(project); long_read_validate(result, project)
    result$fingerprints <- fp
    tmp <- tempfile("long-read-", tmpdir = dirname(long_read_stage_path(project)))
    on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
    saveRDS(result, tmp, version = 3L)
    if (!file.rename(tmp, long_read_stage_path(project))) stop("Could not store long-read extension")
    checksum <- digest::digest(file = long_read_stage_path(project), algo = "sha256", serialize = FALSE)
    checkpoint_update(project, step, "COMPLETED", fp, output_records = list(list(
      path = "intermediate/long_read_extension.rds", checksum = list(
        algorithm = "sha256", scope = "file", value = checksum))))
    result
  }, error = function(e) {
    checkpoint_update(project, step, "FAILED", fp, error = list(code = "LONG_READ_EXTENSION_FAILED")); stop(e)
  })
}
long_read_read <- function(project) {
  if (!length(long_read_inputs(project))) return(NULL)
  cp <- Filter(function(x) x$step_id == "long_read_extension", read_manifest(project$path)$checkpoints)
  if (length(cp) != 1L || cp[[1L]]$status != "COMPLETED") stop("Long-read checkpoint missing")
  verify_checkpoint_outputs(cp[[1L]], project$path)
  x <- readRDS(long_read_stage_path(project)); long_read_validate(x, project)
  if (!identical(x$fingerprints, long_read_fingerprints(project))) stop("Long-read fingerprint changed")
  x
}
