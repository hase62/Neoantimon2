# Additive recognition contracts; never enter frozen candidate/priority identity.
recognition_version <- function() "tcr-recognition-1.0.0"
recognition_id_versions <- function() list(profile = "tcr-feature-profile-id-v1",
  query = "recognition-query-id-v1", evidence = "recognition-evidence-id-v1",
  observation = "tcr-observation-id-v1")
recognition_id <- function(kind, payload) paste0(switch(kind, profile = "tcr-feature-profile", query = "recognition-query",
  evidence = "recognition-evidence", observation = "tcr-observation"), ":", infra_hash(list(
  namespace = recognition_id_versions()[[kind]], kind = kind, payload = payload)))
recognition_nullable <- function(x) {
  if (is.null(x) || !length(x) || (length(x) == 1L && (is.na(x) || identical(x, "")))) NULL else x
}
recognition_number <- function(x, count = FALSE) {
  x <- recognition_nullable(x)
  if (is.null(x)) return(NULL)
  n <- suppressWarnings(as.numeric(x))
  if (length(n) != 1L || !is.finite(n) || (count && (n < 0 || n != floor(n))))
    stop("Invalid recognition numeric evidence")
  n
}
recognition_boolean <- function(x) {
  x <- recognition_nullable(x)
  if (is.null(x)) return(NULL)
  if (is.logical(x) && length(x) == 1L) return(x)
  if (!as.character(x) %in% c("T", "F", "true", "false", "TRUE", "FALSE"))
    stop("Invalid AIRR boolean")
  as.character(x) %in% c("T", "true", "TRUE")
}
recognition_sequence_status <- function(x) {
  if (is.null(recognition_nullable(x))) return("missing_sequence")
  if (!is.character(x) || length(x) != 1L || !grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", x))
    return("unsupported_character")
  "valid"
}
recognition_sort <- function(rows, field) {
  if (!length(rows)) return(list())
  rows[order(vapply(rows, `[[`, "", field), method = "radix")]
}
recognition_unique <- function(rows, field) {
  rows <- recognition_sort(rows, field)
  rows[!duplicated(vapply(rows, `[[`, "", field))]
}
recognition_roles <- function() c("tcr_rearrangement", "tcr_airr", "tcr_repertoire",
  "tcr_cell", "tcr_receptor", "tcr_reactivity", "tcr_vdjdb", "recognition_settings",
  "recognition_output")
recognition_inputs <- function(project) Filter(function(x)
  x$input_role %in% recognition_roles(), project$inputs)
recognition_artifact <- function(input) list(input_id = input$input_id,
  sample_id = input$sample_id, role = input$input_role, version = input$version,
  path = input$declared_path, sha256 = digest::digest(file = input$resolved_path,
    algo = "sha256", serialize = FALSE))
recognition_artifacts <- function(project) lapply(recognition_inputs(project), function(x) {
  if (!identical(integrity_record(x$resolved_path, project$config$checksum,
    x$integrity$declared), x$integrity)) stop("Recognition input changed since project binding")
  recognition_artifact(x)
})
recognition_settings <- function(project) {
  found <- Filter(function(x) x$input_role == "recognition_settings", recognition_inputs(project))
  if (length(found) > 1L) stop("Multiple recognition settings artifacts")
  s <- if (length(found)) jsonlite::read_json(found[[1L]]$resolved_path) else list()
  reject_duplicate_keys(s)
  strict_keys(s, c("version", "providers", "candidate_ids", "profile_ids", "max_queries",
    "max_pairing_hypotheses", "batch_size"), "recognition settings")
  if (length(found) && !identical(s$version, "recognition-run-settings-v1"))
    stop("Unsupported recognition settings")
  for (key in c("max_queries", "max_pairing_hypotheses", "batch_size")) {
    value <- s[[key]]
    if (is.null(value)) value <- switch(key, max_queries = 10000, max_pairing_hypotheses = 1000, batch_size = 256)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value < 1 || value != floor(value))
      stop("Invalid recognition limit: ", key)
    s[[key]] <- value
  }
  for (key in c("candidate_ids", "profile_ids")) if (!is.null(s[[key]])) {
    validate_filter(unlist(s[[key]]), key)
    s[[key]] <- as.list(sort(unique(unlist(s[[key]])), method = "radix"))
  }
  if (is.null(s$providers)) s$providers <- list()
  if (!is.list(s$providers) || !is.null(names(s$providers))) stop("Providers must be an array")
  s$providers <- lapply(s$providers, recognition_provider_config)
  s
}
recognition_presenter <- function(raw) {
  raw <- recognition_nullable(raw)
  if (is.null(raw)) return(list(status = "unknown", mhc_key = NULL, raw = NULL))
  p <- peptide_presenters(list(raw))[[1L]]
  list(status = if (p$status == "RESOLVED") "known" else "unresolved",
    mhc_key = if (p$status == "RESOLVED") p$mhc_key else NULL, raw = raw)
}
recognition_mhc_class <- function(raw) {
  raw <- recognition_nullable(raw)
  if (is.null(raw)) return(NULL)
  value <- toupper(trimws(raw))
  if (value %in% c("I", "MHC-I", "MHC I")) return("I")
  if (value %in% c("II", "MHC-II", "MHC II")) return("II")
  NULL
}
recognition_reactivity_presenter <- function(row) {
  key <- recognition_nullable(row$mhc_key)
  a <- recognition_nullable(row$mhc_allele_1); b <- recognition_nullable(row$mhc_allele_2)
  class <- recognition_mhc_class(row$mhc_class)
  raw <- if (!is.null(key)) key else if (identical(class, "II") && !is.null(a) && !is.null(b))
    paste(a, b, sep = "~") else a
  out <- recognition_presenter(raw)
  out$source_class <- recognition_nullable(row$mhc_class)
  out$normalized_class <- class
  out$raw_alleles <- list(mhc_allele_1 = a, mhc_allele_2 = b)
  unresolved <- function(reason) {
    out$status <- "unresolved"; out$mhc_key <- NULL; out$reason <- reason
    out
  }
  if (is.null(key) && is.null(class) && (!is.null(a) || !is.null(b))) {
    # Preserve the established single canonical Class I allele boundary when
    # class is omitted. An explicitly unrecognized class never authorizes it.
    single_i <- is.null(recognition_nullable(row$mhc_class)) && !is.null(a) && is.null(b) &&
      identical(peptide_presenters(list(a))[[1L]]$hla_class, "I")
    if (!single_i) return(unresolved("unknown_mhc_class"))
  }
  if (is.null(key) && identical(class, "II") && (is.null(a) || is.null(b)))
    return(unresolved("partial_class_ii_presenter"))
  if (!is.null(raw)) {
    canonical <- peptide_presenters(list(raw))[[1L]]
    if (!is.null(class) && canonical$status == "RESOLVED" && !identical(class, canonical$hla_class))
      return(unresolved("conflicting_mhc_class"))
    if (!is.null(key) && canonical$status == "RESOLVED") {
      # Class I's second AIRR side is secondary metadata (e.g. beta2-microglobulin).
      supplied <- if (canonical$hla_class == "I") list(a) else list(a, b)
      for (allele in Filter(Negate(is.null), supplied)) {
        chain <- peptide_presenters(list(allele))[[1L]]$canonical_input
        if (!chain %in% strsplit(canonical$mhc_key, "~", fixed = TRUE)[[1L]])
          return(unresolved("conflicting_mhc_alleles"))
      }
    }
  }
  out
}
recognition_profile <- function(features, completeness = NULL) {
  features <- Filter(Negate(is.null), features)
  if (is.null(completeness)) completeness <- if (all(c("alpha_variable_aa", "beta_variable_aa") %in% names(features)))
    "paired_variable_domain_receptor" else if (all(c("alpha_cdr3_aa", "beta_cdr3_aa") %in% names(features)))
    "paired_cdr3_profile" else if ("beta_cdr3_aa" %in% names(features)) "beta_cdr3_profile" else "chain_profile"
  identity <- list(completeness = completeness, features = features)
  c(list(receptor_feature_profile_id = recognition_id("profile", identity),
    identity_version = recognition_id_versions()$profile, completeness = completeness,
    available_features = as.list(sort(names(features), method = "radix"))), list(features = features))
}
