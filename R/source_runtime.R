source_input_roles <- function() c("star_fusion", "arriba", "splice_structure",
  "spliceai_vcf", "star_junction")
source_inputs <- function(project) Filter(function(x) x$input_role %in% source_input_roles(),
  project$inputs)
source_artifacts <- function(project) lapply(source_inputs(project), function(input) {
  if (is.null(input$version))
    stop("Source provider version declaration required: ", input$input_id)
  observed <- integrity_record(input$resolved_path, project$config$checksum,
    input$integrity$declared)
  if (!identical(observed, input$integrity))
    stop("Source input changed since project binding: ", input$input_id)
  list(input_id = input$input_id, sample_id = input$sample_id,
    input_role = input$input_role, version = input$version,
    sha256 = digest::digest(file = input$resolved_path, algo = "sha256",
      serialize = FALSE))
})
source_stage_path <- function(project)
  file.path(project$path, "intermediate", "source_event_extension.rds")
source_stage_output <- function(project) list(path = "intermediate/source_event_extension.rds",
  checksum = list(algorithm = "sha256", scope = "file", value = digest::digest(
    file = source_stage_path(project), algo = "sha256", serialize = FALSE)))

source_origin <- function(project, event, hypothesis, altered, orf, peptide) {
  ct <- peptide_binding(project, event$structure$sample_id)
  contract <- peptide_contract()
  sequence <- peptide$sequence
  key <- biological_id("peptide_key", list(ligand_sequence = sequence), contract)
  sample_key <- biological_id("sample_peptide_id", list(context_id = ct$context_id,
    patient_id = ct$self_subject_id, sample_id = event$structure$sample_id,
    peptide_key = key), contract)
  tx <- if (event$source_type == "fusion")
    paste0(altered$source_models$left$model$transcript_id, "~",
      altered$source_models$right$model$transcript_id) else
    altered$source_models$reference$model$transcript_id
  identity <- list(transcript_id = tx, transcript_version = NULL,
    reconstructed_sequence_checksum = peptide_digest(altered$sequence),
    peptide_start_aa = peptide$start_aa, peptide_end_aa = peptide$end_aa,
    source_event_keys = list(event$source_event_id), phase_namespace = NULL,
    phase_set = NULL,
    haplotype_id = hypothesis, source_subject_id = ct$source_subject_id,
    source_antigen_key = event$source_event_id, scenario_id = ct$scenario_id)
  id <- biological_id("origin_id", list(sample_peptide_id = sample_key,
    origin_identity = identity), contract)
  ref <- paste0("intermediate/source_event_extension.rds#", orf$protein_origin_id)
  values <- list(schema_version = contract$peptide_version,
    analysis_id = project$analysis_id, patient_id = ct$self_subject_id,
    sample_id = event$structure$sample_id, origin_id = id,
    sample_peptide_id = sample_key, peptide_key = key,
    peptide_sequence = sequence, ligand_sequence = sequence,
    origin_type = if (event$source_type == "fusion") "fusion" else "sequence",
    sequence_id = paste0("protein-sha256:", peptide_digest(orf$protein)),
    sequence_context_id = orf$orf_hypothesis_id,
    peptide_length = as.integer(nchar(sequence)),
    peptide_start_aa = peptide$start_aa, peptide_end_aa = peptide$end_aa,
    transcript_id = tx, transcript_version = NULL,
    protein_id = orf$protein_origin_id, consequence = orf$frame_state,
    source_event_ids = list(), germline_context_ids = list(),
    mutation_positions_json = as.list(orf$novel_aa),
    junction_positions_json = as.list(orf$junction_aa),
    wt_peptide_sequence = NULL,
    wt_definition = orf$reference_comparison,
    wt_comparison_status = peptide$reference_comparison,
    normal_comparison_status = "NO_SINGLE_COUNTERPART",
    upstream_flank = if (peptide$start_aa > 1L) substring(orf$protein,
      max(1L, peptide$start_aa - 5L), peptide$start_aa - 1L) else NULL,
    downstream_flank = if (peptide$end_aa < nchar(orf$protein))
      substring(orf$protein, peptide$end_aa + 1L,
        min(nchar(orf$protein), peptide$end_aa + 5L)) else NULL,
    mutant_sequence_ref = ref, normal_sequence_ref = NULL,
    haplotype_id = hypothesis, phase_status = "NOT_APPLICABLE",
    sequence_certainty = "SOURCE_STRUCTURE_DETERMINISTIC",
    translation_status = orf$termination, isoform_source = "BOUND_REFERENCE",
    origin_validity = "VALID", origin_evidence_status = "NOT_EVALUATED",
    origin_evaluation_ref = paste0("intermediate/source_event_extension.rds#", id),
    analysis_context = ct$context_id,
    source_protein_id = orf$protein_origin_id,
    source_antigen_key = event$source_event_id, source_sequence_ref = ref,
    origin_identity = identity)
  for (role in c("source", "self", "presenter"))
    for (suffix in c("subject_id", "subject_role")) {
      name <- paste(role, suffix, sep = "_")
      values[name] <- ct[name]
    }
  out <- peptide_blank("peptide_origins", contract)
  out[names(values)] <- values
  out
}

source_extension_construct <- function(project, bundle = reference_bundle(project,
  retain_indexes = TRUE)) {
  inputs <- source_inputs(project)
  events <- unlist(lapply(inputs, function(input) switch(input$input_role,
    star_fusion = source_fusion_rows(input, bundle$profile$assembly_release),
    arriba = source_fusion_rows(input, bundle$profile$assembly_release),
    splice_structure = source_splice_rows(input, bundle$profile$assembly_release),
    spliceai_vcf = list(), star_junction = list())), recursive = FALSE)
  predictions <- unlist(lapply(Filter(function(x) x$input_role == "spliceai_vcf", inputs),
    source_spliceai_rows), recursive = FALSE)
  junctions <- unlist(lapply(Filter(function(x) x$input_role == "star_junction", inputs),
    source_star_junction_rows), recursive = FALSE)
  groups <- split(events, vapply(events, `[[`, "", "source_event_id"))
  models <- source_model_catalog(bundle, events)
  reference <- event_reference_open(bundle)
  on.exit(event_reference_close(reference))
  transcripts <- reconstruction_fasta_select(
    bundle$resources$transcript_fasta$resolved_path, names(models))
  model_seq <- lapply(models, source_model_sequence, bundle = bundle,
    reference = reference, transcript_catalog = transcripts)
  model_seq <- Filter(Negate(is.null), model_seq)
  hypotheses <- list(); origins <- list(); peptide_records <- list()
  for (group in groups) {
    event <- group[[1L]]
    structures <- list()
    if (event$source_type == "fusion") {
      for (left in model_seq) for (right in model_seq) {
        altered <- source_fusion_structure(event, left, right, bundle)
        if (!is.null(altered)) structures[[length(structures) + 1L]] <- altered
      }
    } else for (model in model_seq) {
      altered <- source_splice_structure(event, model, reference)
      if (!is.null(altered)) structures[[length(structures) + 1L]] <- altered
    }
    if (!length(structures)) {
      unsupported_mito <- if (event$source_type == "fusion")
        any(c(event$structure$left$contig,
          event$structure$right$contig) == "chrM") else
        event$structure$contig == "chrM"
      hypotheses[[length(hypotheses) + 1L]] <- list(source_event_id = event$source_event_id,
        status = if (unsupported_mito) "UNSUPPORTED_MITOCHONDRIAL" else
          "NO_COMPATIBLE_TRANSCRIPT", observation_ids = as.list(vapply(group,
          `[[`, "", "observation_id")))
      next
    }
    for (altered in structures) {
      ids <- if (event$source_type == "fusion") c(
        altered$source_models$left$model$transcript_id,
        altered$source_models$right$model$transcript_id) else
        altered$source_models$reference$model$transcript_id
      hypothesis <- source_id("transcript_structure_id", list(
        event = event$source_event_id, transcripts = as.list(ids),
        reference = bundle$resource_identity,
        policy = "source-reconstruction-v1"))
      orf <- source_orf(event, altered, hypothesis)
      record <- list(source_event_id = event$source_event_id,
        transcript_structure_id = hypothesis,
        altered_transcript_id = if (is.null(altered$sequence)) NULL else
          source_id("altered_transcript_id", list(
            hypothesis = hypothesis, sequence = altered$sequence)),
        transcript_ids = as.list(ids), status = orf$status,
        observation_ids = as.list(sort(vapply(group, `[[`, "", "observation_id"))),
        altered_transcript = c(altered[c("sequence", "coordinates", "contigs",
          "coding_start", "junction_nt", "source_resume_nt", "blocks",
          "novel_nt_intervals", "skipped_reference_coordinates", "retained_start", "status")],
          list(reference_identity = bundle$resource_identity,
            reconstruction_policy = "source-reconstruction-v1",
            translation_policy = "fusion-splice-orf-v1")),
        orf = orf)
      hypotheses[[length(hypotheses) + 1L]] <- record
      for (p in source_peptide_windows(orf, peptide_policy(project$config$peptides))) {
        origin <- source_origin(project, event, hypothesis, altered, orf, p)
        origins[[length(origins) + 1L]] <- origin
        peptide_records[[length(peptide_records) + 1L]] <- list(
          origin_id = origin$origin_id, peptide_key = origin$peptide_key,
          source_event_id = event$source_event_id,
          transcript_structure_id = hypothesis,
          orf_hypothesis_id = orf$orf_hypothesis_id,
          source_type = event$source_type,
          peptide = p)
      }
    }
  }
  if (length(origins)) {
    keep <- !duplicated(vapply(origins, `[[`, "", "origin_id"))
    origins <- origins[keep]
    peptide_records <- peptide_records[keep]
  }
  self_resource <- Filter(function(x) x$logical_name == "reference.protein_fasta",
    project$resources)
  self_path <- if (length(self_resource)) self_resource[[1L]]$resolved_path else NULL
  self <- peptide_self_search(vapply(origins, `[[`, "", "ligand_sequence"),
    self_path, identity_records(self_resource), sec_encoding = "LITERAL_U")
  registry <- lapply(project$hla, peptide_presenters)
  candidates <- list()
  for (group in split(origins, vapply(origins, `[[`, "", "sample_peptide_id"))) {
    ct <- peptide_binding(project, group[[1L]]$sample_id)
    presenters <- Filter(function(x) x$status == "RESOLVED", registry[[ct$presenter_subject_id]])
    for (presenter in presenters) {
      lengths <- peptide_policy(project$config$peptides)[[
        if (presenter$hla_class == "I") "class_i_lengths" else "class_ii_lengths"]]
      if (!group[[1L]]$peptide_length %in% unlist(lengths)) next
      candidate <- peptide_candidate(group, presenter, project,
        self[[group[[1L]]$ligand_sequence]])
      candidate$not_evaluated_codes <- list("NOT_APPLICABLE")
      candidates[[length(candidates) + 1L]] <- candidate
    }
  }
  list(schema_version = source_extension_version(),
    identity_schema_version = source_identity_schema_version(),
    source_artifacts = source_artifacts(project),
    observations = events, splice_predictions = predictions,
    junction_observations = junctions,
    canonical_events = lapply(groups, function(g) {
      id <- g[[1L]]$source_event_id
      mapped <- Filter(function(x) x$source_event_id == id, hypotheses)
      list(source_event_id = id, structure = g[[1L]]$structure,
        mapping_status = if (length(mapped) && mapped[[1L]]$status %in%
          c("NO_COMPATIBLE_TRANSCRIPT", "UNSUPPORTED_MITOCHONDRIAL"))
          mapped[[1L]]$status else "MAPPED",
        transcript_hypothesis_status = if (length(mapped) > 1L) "AMBIGUOUS" else
          if (length(mapped) == 1L && !mapped[[1L]]$status %in%
              c("NO_COMPATIBLE_TRANSCRIPT", "UNSUPPORTED_MITOCHONDRIAL"))
            "RESOLVED" else "ABSENT",
        observation_ids = as.list(sort(vapply(g, `[[`, "", "observation_id"))))
    }),
    hypotheses = hypotheses, origins = origins, peptides = peptide_records,
    candidates = candidates, self_matches = self)
}

source_validate_extension <- function(x, project) {
  if (!identical(x$schema_version, source_extension_version()))
    stop("Source extension version mismatch")
  if (!identical(x$source_artifacts, source_artifacts(project)))
    stop("Source extension input artifact changed")
  contract <- jsonlite::read_json(system.file("schema", "source-event-2.0.0.json",
    package = "Neoantimon2"))
  if (!identical(contract$version, x$schema_version) ||
      !identical(contract$extends, peptide_contract()$peptide_version))
    stop("Source extension overlay mismatch")
  if (!identical(x$identity_schema_version, source_identity_schema_version()) ||
      !identical(contract$identity_schema_version, x$identity_schema_version))
    stop("Source biological identity schema mismatch")
  unique_field <- function(rows, field) {
    ids <- vapply(rows, `[[`, "", field)
    if (anyDuplicated(ids)) stop("Duplicate source extension identity: ", field)
    ids
  }
  observations <- unique_field(x$observations, "observation_id")
  events <- unique_field(x$canonical_events, "source_event_id")
  hypotheses <- vapply(Filter(function(h) !is.null(h$transcript_structure_id),
    x$hypotheses), `[[`, "", "transcript_structure_id")
  if (anyDuplicated(hypotheses)) stop("Duplicate transcript hypothesis")
  for (e in x$canonical_events) {
    if (!all(unlist(e$observation_ids) %in% observations))
      stop("Unresolved source observation")
  }
  if (any(!vapply(x$observations, `[[`, "", "source_event_id") %in% events))
    stop("Unresolved canonical source event")
  validate_table(x$origins, "peptide_origins", peptide_contract())
  validate_table(x$candidates, "neoantigen_results", peptide_contract())
  ids <- vapply(x$origins, `[[`, "", "origin_id")
  for (candidate in x$candidates)
    if (!all(unlist(candidate$origin_ids) %in% ids) ||
        candidate$final_grade != "NOT_GRADED")
      stop("Invalid source candidate origin or grade")
  invisible(TRUE)
}

source_stage_run <- function(project) {
  if (!length(source_inputs(project))) return(NULL)
  fp <- checkpoint_fingerprints(inputs = list(sources = sorted_records(
    source_artifacts(project))),
    config = list(version = source_extension_version(),
      parser = "source-provider-adapters-v1",
      reconstruction = "source-reconstruction-v1", orf = "fusion-splice-orf-v1",
      peptides = peptide_policy(project$config$peptides), hla = project$hla),
    dependencies = list(reference = sorted_records(identity_records(project$resources))))
  cp <- Filter(function(x) x$step_id == "source_event_extension",
    read_manifest(project$path)$checkpoints)
  if (length(cp) && cp[[1L]]$status == "COMPLETED" &&
      checkpoint_reusable(project, "source_event_extension", fp))
    return(readRDS(source_stage_path(project)))
  if (!length(cp)) checkpoint_update(project, "source_event_extension", fingerprints = fp)
  checkpoint_update(project, "source_event_extension", "RUNNING", fp)
  tryCatch({
    result <- source_extension_construct(project)
    source_validate_extension(result, project)
    tmp <- tempfile("source-extension-", tmpdir = dirname(source_stage_path(project)))
    on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
    saveRDS(result, tmp, version = 3L)
    if (!file.rename(tmp, source_stage_path(project))) stop("Could not store source extension")
    checkpoint_update(project, "source_event_extension", "COMPLETED", fp,
      output_records = list(source_stage_output(project)))
    result
  }, error = function(e) {
    checkpoint_update(project, "source_event_extension", "FAILED", fp,
      error = list(code = "SOURCE_EXTENSION_FAILED"))
    stop(e)
  })
}

source_read <- function(project) {
  if (!length(source_inputs(project))) return(NULL)
  cp <- Filter(function(x) x$step_id == "source_event_extension",
    read_manifest(project$path)$checkpoints)
  if (length(cp) != 1L || cp[[1L]]$status != "COMPLETED")
    stop("Source extension checkpoint missing")
  verify_checkpoint_outputs(cp[[1L]], project$path)
  x <- readRDS(source_stage_path(project))
  if (!identical(x$schema_version, source_extension_version()))
    stop("Source extension schema mismatch")
  source_validate_extension(x, project)
  x
}

source_result_metadata <- function(extension, candidate) {
  origin_ids <- unlist(candidate$origin_ids)
  records <- Filter(function(x) x$origin_id %in% origin_ids, extension$peptides)
  ids <- unique(vapply(records, `[[`, "", "source_event_id"))
  events <- Filter(function(x) x$source_event_id %in% ids, extension$observations)
  canonical <- Filter(function(x) x$source_event_id %in% ids,
    extension$canonical_events)
  structures <- Filter(function(x) x$source_event_id %in% ids,
    extension$hypotheses)
  fusion <- Filter(function(x) x$structure$source_type == "fusion", canonical)
  splice <- Filter(function(x) x$structure$source_type != "fusion", canonical)
  bp <- function(x) paste(x$contig, x$position, x$strand, sep = ":")
  list(source_type = paste(sort(unique(vapply(records, `[[`, "", "source_type"))),
    collapse = ";"), source_event_id = paste(sort(ids), collapse = ";"),
    transcript_structure_id = paste(sort(unique(vapply(records, `[[`, "",
      "transcript_structure_id"))), collapse = ";"),
    orf_hypothesis_id = paste(sort(unique(vapply(records, `[[`, "",
      "orf_hypothesis_id"))), collapse = ";"),
    fusion_partners = paste(vapply(fusion, function(x) paste(vapply(
      x$structure$genes, `[[`, "", "name"), collapse = "--"), ""),
      collapse = ";"),
    fusion_breakpoints = paste(vapply(fusion, function(x) paste(
      bp(x$structure$left), bp(x$structure$right), sep = "--"), ""),
      collapse = ";"),
    splice_event_type = paste(vapply(splice, function(x)
      x$structure$source_type, ""), collapse = ";"),
    splice_donor = paste(vapply(splice, function(x)
      as.character(x$structure$donor), ""), collapse = ";"),
    splice_acceptor = paste(vapply(splice, function(x)
      as.character(x$structure$acceptor), ""), collapse = ";"),
    junction_support = paste(vapply(events, function(x) if (is.null(
      x$support$junction_reads)) "" else as.character(x$support$junction_reads),
      ""), collapse = ";"),
    provider_observations = events, transcript_hypotheses = structures,
    canonical_events = canonical,
    grade_reason = "SOURCE_PRIORITY_POLICY_NOT_IMPLEMENTED")
}

source_prediction_evidence <- function(project, candidate_id, runs = NULL) {
  if (is.null(runs)) runs <- somatic_prediction_runs(project)
  list(evidence = unlist(lapply(runs, function(run) Filter(function(x)
    identical(x$candidate_id, candidate_id), run$evidence)), recursive = FALSE),
    statuses = unlist(lapply(runs, function(run) Filter(function(x)
      identical(x$candidate_id, candidate_id), run$statuses)), recursive = FALSE))
}

source_candidates_with_prediction <- function(project, extension) {
  runs <- somatic_prediction_runs(project)
  lapply(extension$candidates, function(candidate) {
    prediction <- source_prediction_evidence(project, candidate$candidate_id, runs)
    if (length(prediction$evidence)) candidate$prediction_status <- "SUCCESS"
    else if (length(prediction$statuses) && any(vapply(prediction$statuses,
        function(x) identical(x$status, "failed"), TRUE)))
      candidate$prediction_status <- "FAILED"
    candidate
  })
}

source_results_append <- function(out, project, filters) {
  extension <- source_read(project)
  if (is.null(extension)) return(out)
  all_source <- filter_records(source_candidates_with_prediction(project, extension),
    filters$values, filters$reason)
  rows <- Filter(function(x) !x$candidate_id %in% out$candidate_id, all_source)
  added <- schema_table_frame(rows, "neoantigen_results", peptide_contract())
  for (name in setdiff(names(out), names(added)))
    added[[name]] <- out[[name]][rep(NA_integer_, nrow(added))]
  for (name in setdiff(names(added), names(out)))
    out[[name]] <- added[[name]][rep(NA_integer_, nrow(out))]
  added <- added[names(out)]
  merged <- rbind(out, added)
  all_rows <- c(lapply(out$candidate_id, function(id) {
    hit <- Filter(function(x) x$candidate_id == id, all_source)
    if (length(hit)) source_result_metadata(extension, hit[[1L]]) else NULL
  }), lapply(rows,
    function(row) source_result_metadata(extension, row)))
  for (field in c("source_type", "source_event_id", "transcript_structure_id",
    "orf_hypothesis_id", "fusion_partners", "fusion_breakpoints",
    "splice_event_type", "splice_donor", "splice_acceptor",
    "junction_support", "grade_reason"))
    merged[[field]] <- vapply(all_rows, function(x) if (is.null(x)) NA_character_
      else x[[field]], "")
  merged
}
