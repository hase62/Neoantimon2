# Additive T034-T036 contracts, independent of frozen biological/priority IDs.
presentation_version <- function() "presentation-evidence-1.0.0"
presentation_id_versions <- function() list(query = "processing-query-id-v1",
  evidence = "processing-evidence-id-v1", psm = "immunopeptidomics-psm-id-v1",
  observation = "immunopeptidomics-observation-id-v1", apm = "apm-observation-id-v1")
presentation_id <- function(kind, payload) paste0(sub("-id-v1$", "", presentation_id_versions()[[kind]]),
  ":", infra_hash(list(namespace = presentation_id_versions()[[kind]], payload = payload)))
presentation_roles <- function() c("presentation_settings", "processing_output", "immunopeptidomics_psm",
  "immunopeptidomics_mzidentml", "apm_evidence", "hla_expression", "arcas_hla_quant")
presentation_inputs <- function(project) Filter(function(i) i$input_role %in% presentation_roles(), project$inputs)
presentation_artifacts <- function(project) lapply(presentation_inputs(project), function(i) {
  if (!identical(integrity_record(i$resolved_path, project$config$checksum, i$integrity$declared), i$integrity))
    stop("Presentation input changed since binding")
  recognition_artifact(i)
})
presentation_table <- function(path, required = character()) {
  x <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character",
    quote = "", comment.char = "", na.strings = c("", "NA"))
  if (anyDuplicated(names(x)) || !all(required %in% names(x))) stop("Missing/duplicate named presentation columns")
  lapply(seq_len(nrow(x)), function(i) lapply(as.list(x[i, , drop = FALSE]), recognition_nullable))
}
presentation_settings <- function(project) {
  found <- Filter(function(i) i$input_role == "presentation_settings", presentation_inputs(project))
  if (length(found) > 1L) stop("Multiple presentation settings")
  s <- if (length(found)) jsonlite::read_json(found[[1L]]$resolved_path) else list()
  reject_duplicate_keys(s)
  strict_keys(s, c("version", "providers", "max_queries", "experiments"), "presentation settings")
  if (length(found) && !identical(s$version, "presentation-run-settings-v1")) stop("Unsupported presentation settings")
  if (is.null(s$providers)) s$providers <- list()
  if (is.null(s$experiments)) s$experiments <- list()
  if (is.null(s$max_queries)) s$max_queries <- 10000L
  if (!is.numeric(s$max_queries) || length(s$max_queries) != 1L || !is.finite(s$max_queries) ||
      s$max_queries < 1 || s$max_queries != floor(s$max_queries)) stop("Invalid processing query limit")
  if (!is.list(s$providers) || !is.null(names(s$providers))) stop("Providers must be array")
  s$providers <- lapply(s$providers, processing_config)
  if (!is.list(s$experiments) || !is.null(names(s$experiments))) stop("Experiments must be array")
  if (anyDuplicated(vapply(s$experiments, `[[`, "", "input_id"))) stop("Duplicate experiment binding")
  s
}
presentation_context <- function(row, sample) list(sample_id = sample,
  region = recognition_nullable(row$region), timepoint = recognition_nullable(row$timepoint),
  tissue = recognition_nullable(row$tissue), cell_context = recognition_nullable(row$cell_context))
presentation_context_matches <- function(context, candidate) {
  if (!identical(context$sample_id, candidate$sample_id)) return(FALSE)
  # Additional experiment strata require an explicit matching candidate stratum.
  all(vapply(c("region", "timepoint", "tissue", "cell_context"), function(k)
    is.null(context[[k]]) || identical(context[[k]], candidate[[k]]), TRUE))
}
presentation_number <- function(x, nonnegative = FALSE) {
  n <- recognition_number(x)
  if (!is.null(n) && nonnegative && n < 0) stop("Negative presentation measurement")
  n
}
presentation_bool <- function(x) {
  if (identical(x, "1")) return(TRUE)
  if (identical(x, "0")) return(FALSE)
  recognition_boolean(x)
}
presentation_sequence <- function(x) {
  scalar_text(x, "peptide backbone")
  if (!grepl("^[ACDEFGHIKLMNPQRSTVWYU]+$", x)) stop("Invalid literal peptide backbone")
  x
}
presentation_seal <- function(row, kind, field) {
  row$fingerprint <- infra_hash(row)
  row[[field]] <- presentation_id(kind, row)
  row
}
presentation_contract <- function() jsonlite::read_json(system.file("schema",
  "presentation-evidence-1.0.0.json", package = "Neoantimon2"))
