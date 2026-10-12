# T040: opt-in policy, independent of every frozen biological contract.
vaccine_policy_version <- function() "vaccine-set-policy-1.0.0"
vaccine_selection_version <- function() "vaccine-set-selection-1.0.0"
vaccine_algorithm <- function() "lexicographic-greedy-v1"
vaccine_namespaces <- function() list(policy = "vaccine-policy-id-v1",
  set = "vaccine-set-id-v1", redundancy = "vaccine-redundancy-group-v1")
vaccine_id <- function(kind, x) paste0("vaccine-", kind, ":", infra_hash(list(
  namespace = vaccine_namespaces()[[kind]], content = longitudinal_plain(x))))
vaccine_strings <- function(x, label) {
  if (!length(x)) return(list())
  if (!is.null(names(x)) || !(is.list(x) || is.character(x))) stop(label, " requires an array")
  for (s in x) scalar_text(s, label)
  as.list(sort(unique(unlist(x)), method = "radix"))
}
vaccine_policy_normalize <- function(x) {
  reject_duplicate_keys(x)
  counts <- c("max_peptides", "min_class_i", "min_class_ii", "max_class_i", "max_class_ii",
    "max_per_presenter", "max_per_redundancy_group", "min_length", "max_length")
  arrays <- c("required_presenters", "allowlist", "blocklist", "forbidden_residues", "forbidden_motifs", "require_evidence")
  flags <- c("forbid_confirmed_exact_self", "fill_to_capacity", "reduce_overlap")
  strict_keys(x, c("version", counts, arrays, flags), "vaccine policy")
  for (k in names(x)) if (is.null(x[[k]])) stop("Explicit null is invalid for vaccine policy field: ", k)
  if (!identical(x$version, vaccine_policy_version())) stop("Unsupported vaccine policy version")
  if (is.null(x$max_peptides)) stop("Explicit max_peptides is required; no default")
  for (k in counts) if (!is.null(x[[k]])) {
    v <- x[[k]]
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || !is.finite(v) || v < 0 || v != floor(v) || v > .Machine$integer.max)
      stop("Invalid vaccine count: ", k)
    x[[k]] <- as.integer(v)
  }
  for (k in c("min_class_i", "min_class_ii")) if (is.null(x[[k]])) x[[k]] <- 0L
  for (k in arrays) x[[k]] <- vaccine_strings(x[[k]], k)
  for (k in flags) {
    if (is.null(x[[k]])) x[[k]] <- FALSE
    if (!is.logical(x[[k]]) || length(x[[k]]) != 1L || is.na(x[[k]])) stop("Invalid vaccine boolean: ", k)
  }
  if (any(!grepl("^[A-Z]$", unlist(x$forbidden_residues)))) stop("Forbidden residues must be single uppercase residues")
  if (any(!grepl("^[A-Z]+$", unlist(x$forbidden_motifs)))) stop("Forbidden motifs are literal uppercase sequences")
  if (length(setdiff(unlist(x$require_evidence), c("exact_self_evaluated", "longitudinal_evaluable", "linked_clonal", "ms_observed"))))
    stop("Unsupported vaccine evidence requirement")
  # Canonical presenter spelling is required: aliases never create hidden semantics.
  for (p in x$required_presenters) {
    presenter <- peptide_presenters(list(p))[[1L]]
    if (!identical(presenter$status, "RESOLVED") || !identical(presenter$mhc_key, p)) stop("Required presenter must be a canonical mhc_key")
  }
  x[order(names(x), method = "radix")]
}
vaccine_policy_bind <- function(x, base = NULL) {
  if (is.list(x) && identical(x$binding_version, "vaccine-policy-binding-1.0.0")) {
    strict_keys(x, c("binding_version", "input", "normalized", "policy_id", "physical", "raw_text"), "vaccine policy binding")
    n <- vaccine_policy_normalize(x$input)
    if (!identical(canonical_bytes(n), canonical_bytes(vaccine_policy_normalize(x$normalized))) ||
        !identical(x$policy_id, vaccine_id("policy", n))) stop("Vaccine policy binding mismatch")
    x$normalized <- n
    return(x)
  }
  physical <- list(); raw_text <- NULL
  if (is.character(x)) {
    scalar_text(x, "vaccine policy path")
    path <- resolve_path(x, base)
    bytes <- readBin(path, "raw", n = file.info(path)$size)
    raw_text <- rawToChar(bytes)
    physical <- list(path = path, sha256 = digest::digest(bytes, algo = "sha256", serialize = FALSE))
    x <- if (grepl("\\.ya?ml$", path, ignore.case = TRUE)) yaml::yaml.load(raw_text, eval.expr = FALSE) else jsonlite::fromJSON(raw_text, simplifyVector = FALSE)
  }
  n <- vaccine_policy_normalize(x)
  list(binding_version = "vaccine-policy-binding-1.0.0", input = x,
    normalized = n, policy_id = vaccine_id("policy", n), physical = physical, raw_text = raw_text)
}
vaccine_enabled <- function(project) !isTRUE(project$vaccine_accessor_guard) &&
  !isTRUE(project$longitudinal_accessor_guard) && !is.null(project$config$vaccine_set)
vaccine_reason <- function(code, constraint, detail = NULL) list(code = code, constraint = constraint, detail = detail)
vaccine_eligibility <- function(p, policy) {
  r <- p$frozen_result; seq <- p$ligand_sequence; reasons <- list(); checks <- list()
  add <- function(name, passed, detail = NULL, missing = FALSE) {
    checks[[length(checks) + 1L]] <<- list(constraint = name,
      status = if (missing) "insufficient_evidence" else if (passed) "satisfied" else "failed", detail = detail)
    if (!passed) reasons[[length(reasons) + 1L]] <<- vaccine_reason(
      if (missing) "insufficient_evidence_for_requirement" else if (name %in% c("allowlist", "blocklist")) "blocked_by_policy" else "hard_constraint_failure", name, detail)
  }
  invalid <- any(unlist(r[c("final_grade", "grading_status", "grade_status", "technical_eligibility")]) %in% c("INVALID", "invalid")) ||
    any(vapply(p$origins, function(o) identical(o$origin_validity, "INVALID"), TRUE)) &&
    all(vapply(p$origins, function(o) identical(o$origin_validity, "INVALID"), TRUE))
  add("frozen_technical_validity", !invalid, r[c("final_grade", "technical_eligibility")])
  ids <- c(p$candidate_id, p$candidate_key)
  if (length(policy$allowlist)) add("allowlist", any(ids %in% unlist(policy$allowlist)))
  if (length(policy$blocklist)) add("blocklist", !any(ids %in% unlist(policy$blocklist)))
  if (!is.null(policy$min_length)) add("min_length", nchar(seq) >= policy$min_length)
  if (!is.null(policy$max_length)) add("max_length", nchar(seq) <= policy$max_length)
  for (res in policy$forbidden_residues) add("forbidden_residue", !grepl(res, seq, fixed = TRUE), res)
  for (motif in policy$forbidden_motifs) add("forbidden_motif", !grepl(motif, seq, fixed = TRUE), motif)
  if (policy$forbid_confirmed_exact_self) add("confirmed_exact_self", !isTRUE(r$exact_self_match) && !isTRUE(r$self_exact_match),
    list(exact_self_match = r$exact_self_match, self_exact_match = r$self_exact_match, status = r$self_status))
  for (req in policy$require_evidence) {
    met <- switch(req, exact_self_evaluated = identical(r$self_status, "EVALUATED") || identical(r$self_match_status, "evaluated"),
      longitudinal_evaluable = !is.null(p$trajectory) && p$trajectory$summary$trajectory_state != "insufficient_to_classify",
      linked_clonal = length(p$clonal_evidence) > 0L,
      ms_observed = identical(p$evidence$presentation_evidence$summary$immunopeptidomics_status, "exact_unmodified_ligand_observed"))
    add(req, isTRUE(met), "Explicit policy evidence requirement; missing is not negative", !isTRUE(met))
  }
  list(status = if (length(reasons)) "ineligible" else "eligible", checks = checks, reasons = reasons)
}
