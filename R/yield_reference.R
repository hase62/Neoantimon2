# Optional Experimental Yield v0.1 contracts; no native scientific writes.
yield_text <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(trimws(x)))
    stop(name, " must be one nonempty string", call. = FALSE)
  x
}
yield_number <- function(x, name, lower = -Inf, upper = Inf, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < lower || x > upper || (integer && x != floor(x)))
    stop("Invalid ", name, call. = FALSE)
  x
}
yield_object <- function(x, schema) {
  artifact <- NULL
  if (is.character(x)) {
    yield_text(x, "JSON path")
    if (!file.exists(x) || dir.exists(x)) stop("Yield JSON file unavailable: ", x)
    artifact <- list(path = normalizePath(x), sha256 = digest::digest(file = x,
      algo = "sha256", serialize = FALSE))
    x <- jsonlite::read_json(x, simplifyVector = FALSE)
  }
  if (!is.list(x) || !identical(x$schema_version, schema))
    stop("Unsupported yield schema: require ", schema, call. = FALSE)
  for (f in c("id", "version", "provenance")) {
    if (f == "provenance") {
      if (is.null(x[[f]]) || !length(x[[f]])) stop("Yield provenance is required")
    } else yield_text(x[[f]], f)
  }
  list(data = x, artifact = if (is.null(artifact)) list(
    sha256 = digest::digest(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null",
      digits = NA), algo = "sha256", serialize = FALSE), input = "explicit R object") else artifact)
}
yield_peptide <- function(x) {
  if (!is.character(x) || anyNA(x)) stop("Peptide sequences must be strings")
  out <- toupper(trimws(x))
  if (any(!grepl("^[ACDEFGHIKLMNPQRSTVWYU]+$", out)))
    stop("Yield peptide requires literal amino acids, no wildcard or modification")
  out
}
yield_hla <- function(x) {
  if (!is.character(x) || anyNA(x)) stop("Exact HLA strings are required")
  vapply(x, function(value) {
    parts <- strsplit(toupper(trimws(value)), "~", fixed = TRUE)[[1L]]
    if (!length(parts) || length(parts) > 2L) stop("Invalid exact HLA")
    parts <- vapply(parts, function(p) {
      p <- sub("^HLA-", "", p)
      # Only punctuation normalization. No field truncation or fuzzy matching.
      p <- sub("^([A-Z]+[0-9]*)([0-9]{2}:)", "\\1*\\2", p)
      if (!grepl("^(A|B|C|E|F|G|DRA|DRB[1345]|DPA1|DPB1|DQA1|DQB1)\\*[0-9]{2,3}:[0-9]{2,3}(:[0-9]{2,3}){0,2}[NLSCAQ]?$", p))
        stop("Unsupported exact HLA notation: ", value, call. = FALSE)
      paste0("HLA-", p)
    }, "", USE.NAMES = FALSE)
    if (length(parts) == 2L) {
      alpha <- grepl("^HLA-D(RA|PA1|QA1)\\*", parts)
      if (sum(alpha) != 1L) stop("Require explicit Class II alpha and beta")
      parts <- c(parts[alpha], parts[!alpha])
      loci <- sub("\\*.*", "", parts)
      valid <- (loci[1] == "HLA-DRA" && grepl("^HLA-DRB", loci[2])) ||
        identical(loci, c("HLA-DPA1", "HLA-DPB1")) ||
        identical(loci, c("HLA-DQA1", "HLA-DQB1"))
      if (!valid) stop("Incompatible explicit Class II chains")
    }
    paste(parts, collapse = "~")
  }, "", USE.NAMES = FALSE)
}
yield_reference <- function(reference) {
  if (is.null(reference) || (is.character(reference) && length(reference) == 1L &&
      !file.exists(reference))) return(list(available = FALSE, status = "REFERENCE_UNAVAILABLE",
        declared = reference))
  obj <- yield_object(reference, "experimental-reference-0.1.0"); db <- obj$data
  if (is.null(db$patient_identity) || !db$patient_identity %in% c("study_local", "global_verified"))
    stop("Reference requires explicit patient_identity: study_local or global_verified")
  records <- db$records
  if (is.data.frame(records)) records <- lapply(seq_len(nrow(records)), function(i) as.list(records[i, , drop = FALSE]))
  if (!is.list(records)) stop("Reference records must be an array")
  records <- lapply(records, function(r) {
    for (f in c("record_id", "peptide_sequence", "mhc_key", "endpoint", "assay",
                "result", "study_id", "citation")) yield_text(r[[f]], paste("Reference", f))
    if (!r$result %in% c("positive", "negative")) stop("Only explicitly tested positive/negative records are accepted")
    r$source_peptide_sequence <- r$peptide_sequence; r$source_mhc_key <- r$mhc_key
    r$peptide_sequence <- yield_peptide(r$peptide_sequence); r$mhc_key <- yield_hla(r$mhc_key)
    if (is.null(r$patient_id) || (length(r$patient_id) == 1L && is.na(r$patient_id))) r$patient_id <- ""
    if (!is.character(r$patient_id) || length(r$patient_id) != 1L) stop("Invalid reference patient_id")
    r
  })
  ids <- vapply(records, `[[`, "", "record_id")
  if (anyDuplicated(ids)) stop("Duplicate experimental record_id")
  db$records <- records
  list(available = TRUE, status = "AVAILABLE", data = db, artifact = obj$artifact)
}
yield_pair_input <- function(x) {
  if (!is.data.frame(x) || !all(c("peptide_sequence", "mhc_key") %in% names(x)))
    stop("Require a data frame with peptide_sequence and mhc_key")
  x$source_peptide_sequence <- x$peptide_sequence; x$source_mhc_key <- x$mhc_key
  x$peptide_sequence <- yield_peptide(x$peptide_sequence); x$mhc_key <- yield_hla(x$mhc_key)
  x
}

#' Describe exact and related experimental evidence without assigning probability
#' @param x Data frame with peptide_sequence and mhc_key. Extra columns are retained.
#' @param reference Explicit versioned experimental-reference-0.1.0 JSON path or list.
#' @param endpoint One exact experimental endpoint. Required with an available reference.
#' @param assay One exact assay type. Incompatible assays are never pooled.
#' @param min_records Minimum independent tested peptide/patient/HLA units for C.
#' @param min_patients Minimum verified independent patients for C.
#' @return Data frame, one row per input, with A-D category, both B subtypes,
#'   positive/negative/conflicting statuses and complete matched records in list columns.
#' @export
neo_yield_evidence <- function(x, reference = NULL, endpoint = NULL, assay = NULL,
    min_records = 100L, min_patients = 3L) {
  x <- yield_pair_input(x)
  yield_number(min_records, "min_records", 1, integer = TRUE)
  yield_number(min_patients, "min_patients", 1, integer = TRUE)
  ref <- yield_reference(reference); records <- list(); available <- ref$available
  if (available) {
    yield_text(endpoint, "Evidence endpoint"); yield_text(assay, "Evidence assay")
    records <- Filter(function(r) r$endpoint == endpoint && r$assay == assay, ref$data$records)
    if (!length(records)) { available <- FALSE; ref$status <- "REFERENCE_ENDPOINT_ASSAY_UNAVAILABLE" }
  }
  ref_peptides <- vapply(records, `[[`, "", "peptide_sequence")
  ref_hlas <- vapply(records, `[[`, "", "mhc_key")
  ref_lengths <- nchar(ref_peptides)
  status <- function(rows) {
    values <- unique(vapply(rows, `[[`, "", "result"))
    if (!length(values)) "none" else if (length(values) > 1L) "conflicting" else values
  }
  rows <- lapply(seq_len(nrow(x)), function(i) {
    pep <- x$peptide_sequence[i]; hla <- x$mhc_key[i]
    matches <- list(); bytype <- list(A = list(), B1 = list(), B2 = list())
    relevant <- which(ref_peptides == pep | (ref_hlas == hla & ref_lengths == nchar(pep)))
    for (r in records[relevant]) {
      exact_peptide <- identical(r$peptide_sequence, pep); exact_hla <- identical(r$mhc_key, hla)
      type <- NULL; substitutions <- list()
      if (exact_peptide && exact_hla) type <- "A"
      else if (exact_peptide && !exact_hla) type <- "B1"
      else if (exact_hla && nchar(r$peptide_sequence) == nchar(pep)) {
        a <- strsplit(r$peptide_sequence, "", fixed = TRUE)[[1L]]
        b <- strsplit(pep, "", fixed = TRUE)[[1L]]; at <- which(a != b)
        if (length(at) == 1L) {
          type <- "B2"; substitutions <- list(position_1based = at,
            reference_residue = a[at], candidate_residue = b[at])
        }
      }
      if (!is.null(type)) {
        m <- list(match_type = type, substitutions = substitutions, record = r)
        matches[[length(matches) + 1L]] <- m
        bytype[[type]][[length(bytype[[type]]) + 1L]] <- r
      }
    }
    own <- records[which(ref_hlas == hla)]
    # Do not inflate tested units with technical duplicates or conflicting replicates.
    keys <- vapply(own, function(r) paste(r$study_id, r$patient_id, r$peptide_sequence,
      r$mhc_key, sep = "\037"), "")
    ntested <- length(unique(keys)); patients <- Filter(function(r) nzchar(r$patient_id), own)
    npatients <- if (!length(patients)) 0L else if (ref$data$patient_identity == "global_verified")
      length(unique(vapply(patients, `[[`, "", "patient_id"))) else {
        studies <- split(patients, vapply(patients, `[[`, "", "study_id"))
        max(vapply(studies, function(z) length(unique(vapply(z, `[[`, "", "patient_id"))), 0L))
      }
    subtypes <- names(bytype)[names(bytype) != "A" & lengths(bytype) > 0L]
    category <- if (!available) NA_character_ else if (length(bytype$A)) "A" else
      if (length(subtypes)) "B" else if (ntested >= min_records && npatients >= min_patients) "C" else "D"
    list(evidence_category = category, evidence_status = if (available) "AVAILABLE" else ref$status,
      b_subtypes = subtypes, A_status = status(bytype$A), B1_status = status(bytype$B1),
      B2_status = status(bytype$B2), hla_tested_units = ntested, hla_verified_patients = npatients,
      matches = matches)
  })
  for (f in c("evidence_category", "evidence_status", "A_status", "B1_status", "B2_status"))
    x[[f]] <- vapply(rows, `[[`, "", f)
  for (f in c("hla_tested_units", "hla_verified_patients")) x[[f]] <- vapply(rows, `[[`, 0L, f)
  for (f in c("b_subtypes", "matches")) x[[f]] <- I(lapply(rows, `[[`, f))
  attr(x, "reference") <- ref
  attr(x, "evidence_policy") <- list(version = "yield-evidence-0.1.0", endpoint = endpoint,
    assay = assay, min_records = min_records, min_patients = min_patients,
    patient_count = if (ref$available) ref$data$patient_identity else "unavailable",
    categories_do_not_modify_probability = TRUE)
  x
}
