# T043-T045: additive Phase A overlay; no peptide or presentation objects.
allo_version <- function() "allo-foundation-1.0.0"
allo_namespaces <- function() c(context = "allo-context-id-v1",
  sequence = "hla-source-sequence-id-v1", occurrence = "hla-source-occurrence-id-v1",
  universe = "self-universe-id-v1", policy = "allo-policy-id-v1",
  analysis = "allo-analysis-id-v1")
allo_id <- function(kind, payload) {
  namespace <- allo_namespaces()[[kind]]
  paste0(namespace, ":", infra_hash(list(namespace = namespace, payload = payload)))
}
allo_set <- function(x) {
  if (!length(x)) return(list())
  as.list(sort(unique(stringi::stri_trans_nfc(unlist(x))), method = "radix"))
}
allo_records <- function(x) {
  if (is.null(x)) return(list())
  if (!is.list(x) || !is.null(names(x))) stop("Allo records must be arrays")
  if (!length(x)) return(list())
  key <- vapply(x, canonical_bytes, "")
  x[order(key, method = "radix")][!duplicated(sort(key, method = "radix"))]
}
allo_keys <- function(x, required, optional = character(), label = "allo record", nullable = character()) {
  strict_keys(x, c(required, optional), label)
  if (!length(x) && is.null(names(x))) stop(label, " must be an object")
  if (length(setdiff(required, names(x)))) stop(label, ": missing required field")
  for (k in setdiff(names(x), nullable)) if (is_missing(x[[k]])) stop(label, ": null field ", k)
  invisible(TRUE)
}
allo_text <- function(x, label) {
  scalar_text(x, label)
  stringi::stri_trans_nfc(x)
}
allo_boolean <- function(x, label) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) stop(label, " must be boolean")
  x
}
allo_strings <- function(x, label, nonempty = FALSE) {
  if (!is.list(x) || !is.null(names(x))) stop(label, " must be a string array")
  for (v in x) scalar_text(v, label)
  if (nonempty && !length(x)) stop(label, " must be nonempty")
  allo_set(x)
}
allo_state <- function(status, reasons = list(), required = list(), examined = list(),
                       scope, evidence_refs = list()) list(status = status,
  reason_codes = allo_set(reasons), required_components = allo_set(required),
  examined_components = allo_set(examined), scope = scope,
  evidence_refs = allo_set(evidence_refs), value = NULL)
allo_comparison_states <- function() c("supported_nonself_within_scope", "self_match",
  "insufficient_search", "not_evaluated", "ambiguous", "not_applicable")

allo_config <- function(x, base = NULL) {
  allo_keys(x, c("version", "subjects", "contexts", "typing", "reference", "self_search"))
  if (!identical(x$version, allo_version())) stop("Unsupported allo foundation version")
  reject_duplicate_keys(x)
  for (k in c("subjects", "contexts", "typing", "self_search")) allo_records(x[[k]])
  ref <- x$reference
  allo_keys(ref, "path", c("sha256", "manifest"), "allo reference binding")
  path <- resolve_path(ref$path, base)
  if (!file.exists(path) || dir.exists(path)) stop("Pinned allo reference artifact unavailable")
  hash <- digest::digest(file = path, algo = "sha256", serialize = FALSE)
  if ("sha256" %in% names(ref) && !identical(ref$sha256, hash)) stop("Allo reference integrity mismatch")
  manifest <- jsonlite::read_json(path)
  reject_duplicate_keys(manifest)
  if ("manifest" %in% names(ref) && !identical(canonical_bytes(ref$manifest), canonical_bytes(manifest)))
    stop("Pinned allo reference content mismatch")
  # Snapshot actual mapping/coverage as well as the release. Resume never consults
  # a live database or substitutes a changed artifact for this bound input.
  x$reference <- list(path = path, sha256 = hash, manifest = manifest)
  x
}
allo_enabled <- function(project) !is.null(project$config$allo)
