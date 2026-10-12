# Explicit downstream validation overlay. The frozen authority and ID inputs
# are unchanged; no version is inserted into biological identity payloads.
peptide_contract <- local({
  cached <- NULL
  function() {
  if (!is.null(cached)) return(cached)
  c <- schema_contract()
  d <- jsonlite::read_json(system.file("schema", "peptide-1.1.0.json",
    package = "Neoantimon2"))
  stopifnot(identical(c$schema_version, d$extends),
    identical(c$id_schema_version, d$id_schema_version))
  pattern <- paste0("^[", d$alphabet, "]+$")
  for (section in c("tables", "records")) for (name in names(d[[section]])) {
    for (i in seq_along(c[[section]][[name]]$fields)) {
      f <- c[[section]][[name]]$fields[[i]]
      if (f$name %in% unlist(d[[section]][[name]])) f$pattern <- pattern
      if (section == "tables" && f$name == "schema_version") f$enum <- list(d$version)
      c[[section]][[name]]$fields[[i]] <- f
    }
  }
  c$output_schema_versions <- list(events = "1.0.0", peptide_origins = d$version,
    neoantigen_results = d$version, multisample_summary = d$version)
  c$peptide_version <- d$version
  c$peptide_alphabet <- d$alphabet
  cached <<- c
  c
}})

runtime_contract <- function(snapshot) {
  versions <- snapshot$run_manifest$output_schema_versions
  new <- peptide_contract()
  if (identical(versions, new$output_schema_versions)) return(new)
  # Unknown, mixed, or absent versions are never implicitly treated as T010+.
  schema_contract()
}

peptide_sequence_validate <- function(sequence, nullable = FALSE) {
  validate_value(sequence, list(type = "string", nullable = nullable,
    pattern = paste0("^[", peptide_contract()$peptide_alphabet, "]+$")),
    peptide_contract(), "translated_peptide")
  invisible(TRUE)
}

peptide_translation_provenance <- function(context) {
  t <- context$translation
  list(translation_id = t$translation_id, reconstruction_id = context$reconstruction_id,
    translation_system = t$translation_system,
    translation_exception_ids = as.list(t$translation_exception_ids),
    annotated_translation_exception_ids = as.list(t$spec$translation_exception_ids),
    translation_spec_id = t$spec$spec_id,
    sec_sites_used = as.list(t$sec_sites_used), reference_identity = t$reference_identity,
    mutant_boundary = t$mutant_boundary, start_effect = t$start_effect,
    frame_effect = t$frame_effect, termination_effect = t$termination_effect,
    protein_effect = t$protein_effect)
}
