# Supplied peptides are input occurrences, never inferred molecular Events.
peptide_fasta_inputs <- function(project) Filter(function(x)
  identical(x$input_role, "peptide_fasta"), project$inputs)

peptide_fasta_verify <- function(project) {
  for (input in peptide_fasta_inputs(project))
    if (!identical(input$integrity, integrity_record(input$resolved_path,
        "sha256", input$integrity$declared)))
      stop("Peptide FASTA changed since project binding: ", input$input_id)
  invisible(TRUE)
}

peptide_fasta_read <- function(input) {
  fail <- function(line, reason) stop(input$resolved_path, ":", line, ": ", reason, call. = FALSE)
  con <- reference_connection(input$resolved_path); on.exit(close(con))
  records <- list(); current <- NULL; seen <- character(); line <- 0L
  finish <- function() {
    if (is.null(current)) return()
    if (!nzchar(current$sequence)) fail(current$header_line, paste("Empty FASTA record", current$record_id))
    current$sequence <- toupper(current$sequence)
    records[[length(records) + 1L]] <<- current
  }
  repeat {
    chunk <- readLines(con, n = 10000L, warn = FALSE)
    if (!length(chunk)) break
    for (value in chunk) {
      line <- line + 1L
      if (startsWith(value, ">")) {
        finish()
        header <- substring(value, 2L)
        id <- strsplit(trimws(header), "[[:space:]]+")[[1L]][1L]
        if (is.na(id) || !nzchar(id)) fail(line, "Empty FASTA record ID")
        if (id %in% names(seen)) fail(line, paste0("Duplicate FASTA record ID '", id,
          "' (first header at line ", seen[[id]], ")"))
        seen[[id]] <- as.character(line)
        current <- list(record_id = id, full_header = header, header_line = line,
          record_order = length(records) + 1L, sequence = "")
      } else {
        if (is.null(current)) fail(line, "Sequence before first FASTA header")
        if (!nzchar(value)) fail(line, "Empty sequence line")
        if (grepl("[[:space:]]", value)) fail(line, "Internal whitespace in sequence line")
        tryCatch(peptide_sequence_validate(toupper(value)), error = function(e)
          fail(line, paste("Invalid peptide residue:", conditionMessage(e))))
        current$sequence <- paste0(current$sequence, value)
      }
    }
  }
  finish()
  if (!length(records)) fail(1L, "Empty FASTA file")
  records
}

peptide_fasta_construct <- function(project, contract = peptide_contract()) {
  peptide_fasta_verify(project)
  origins <- peptides <- statuses <- list()
  policy <- peptide_policy(project$config$peptides)
  for (input in peptide_fasta_inputs(project)) {
    records <- peptide_fasta_read(input)
    file_sha256 <- external_sha(input$resolved_path)
    for (record in records) {
    ct <- peptide_binding(project, input$sample_id); sequence <- record$sequence
    provenance <- c(list(input_id = input$input_id, input_sha256 = file_sha256,
      input_role = "peptide_fasta", sample_id = input$sample_id), record)
    record_key <- paste0("fasta-record:", infra_hash(provenance))
    key <- biological_id("peptide_key", list(ligand_sequence = sequence), contract)
    sample_key <- biological_id("sample_peptide_id", list(context_id = ct$context_id,
      patient_id = ct$self_subject_id, sample_id = input$sample_id, peptide_key = key), contract)
    identity <- list(transcript_id = NULL, transcript_version = NULL,
      reconstructed_sequence_checksum = peptide_digest(sequence), peptide_start_aa = 1L,
      peptide_end_aa = nchar(sequence), source_event_keys = list(), phase_namespace = NULL,
      phase_set = NULL, haplotype_id = NULL, source_subject_id = ct$source_subject_id,
      source_antigen_key = record_key, scenario_id = ct$scenario_id)
    id <- biological_id("origin_id", list(sample_peptide_id = sample_key, origin_identity = identity), contract)
    origin <- peptide_blank("peptide_origins", contract)
    values <- list(schema_version = contract$peptide_version, analysis_id = project$analysis_id,
      patient_id = ct$self_subject_id, sample_id = input$sample_id, analysis_context = ct$context_id,
      origin_id = id, sample_peptide_id = sample_key, peptide_key = key,
      peptide_sequence = sequence, ligand_sequence = sequence, peptide_length = nchar(sequence),
      peptide_start_aa = 1L, peptide_end_aa = nchar(sequence), origin_type = "sequence",
      sequence_id = paste0("protein-sha256:", peptide_digest(sequence)),
      sequence_context_id = record_key, sequence_certainty = "SUPPLIED_PEPTIDE_SEQUENCE",
      translation_status = "NOT_APPLICABLE", phase_status = "NOT_EVALUATED",
      wt_comparison_status = "NOT_EVALUATED", normal_comparison_status = "NOT_EVALUATED",
      origin_validity = "VALID", origin_evidence_status = "NOT_EVALUATED",
      origin_warning_codes = list(),
      source_domain = "provided_peptide", source_antigen_key = record_key,
      source_annotation_ref = paste0("intermediate/peptide_candidates.rds#", id), origin_identity = identity)
    for (role in c("source", "self", "presenter")) for (suffix in c("subject_id", "subject_role")) {
      name <- paste(role, suffix, sep = "_"); values[name] <- ct[name]
    }
    origin[names(values)] <- values
    origins[[length(origins) + 1L]] <- origin
    peptides[[length(peptides) + 1L]] <- list(peptide_id = id, sample_peptide_id = sample_key,
      sequence_key = key, sample_id = input$sample_id, peptide_sequence = sequence,
      analysis_context = ct$context_id, source_event_ids = list(), germline_context_ids = list(),
      provenance = list(phase_condition = list(status = "NOT_APPLICABLE"), supplied_peptide = provenance),
      provenance_occurrences = list())
    presenters <- Filter(function(p) p$status == "RESOLVED",
      peptide_presenters(project$hla[[ct$presenter_subject_id]]))
    eligible <- vapply(presenters, function(p) nchar(sequence) %in%
      unlist(policy[[if (p$hla_class == "I") "class_i_lengths" else "class_ii_lengths"]]), TRUE)
    statuses[[length(statuses) + 1L]] <- list(origin_id = id, input_id = input$input_id,
      record_id = record$record_id, peptide_length = nchar(sequence),
      status = if (any(eligible)) "CANDIDATE_LENGTH_ELIGIBLE" else "NOT_CANDIDATIZED",
      reason = if (!length(presenters)) "NO_RESOLVED_PRESENTER" else
        if (!any(eligible)) "LENGTH_OUTSIDE_CONFIGURED_CLASS_LENGTHS" else NULL)
    }
  }
  peptide_fasta_verify(project)
  list(origins = origins, peptides = peptides, statuses = statuses)
}
