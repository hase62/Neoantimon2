indirect_provider_arguments <- function(x) {
  allo_keys(x, "BA", label = "semantic provider arguments")
  ba <- allo_boolean(x$BA, "BA semantic argument")
  if (!ba) stop("Strict profile requires BA semantic argument = true")
  list(BA = ba)
}
indirect_provider_provenance <- function(x, p) {
  if (p$evidence_kind == "synthetic_fixture") {
    allo_keys(x, c("kind", "fixture_id", "license", "source", "execution"), "metadata", "synthetic provider provenance")
    if (!identical(x$kind, "synthetic_fixture") || !identical(x$execution, "not_executed"))
      stop("Synthetic fixture must declare not_executed provenance")
    for (k in c("fixture_id", "license", "source")) allo_text(x[[k]], paste("synthetic", k))
  } else {
    allo_keys(x, c("kind", "query_artifact", "result_artifact", "executable_artifact", "model_artifact",
      "license", "execution"), "metadata", "precomputed provider provenance")
    if (!identical(x$kind, "precomputed_predictor_output")) stop("Precomputed provider provenance kind mismatch")
    for (k in c("query_artifact", "result_artifact", "executable_artifact", "model_artifact")) {
      a <- x[[k]]
      allo_keys(a, c("identity", "reference", "sha256"), label = paste(k, "provenance"))
      for (field in c("identity", "reference")) allo_text(a[[field]], field)
      if (!is.character(a$sha256) || length(a$sha256) != 1L || is.na(a$sha256) ||
        !grepl("^[0-9a-f]{64}$", a$sha256)) stop("Provider artifact checksum must be SHA-256")
    }
    if (!identical(x$executable_artifact$identity, p$build) || !identical(x$model_artifact$identity, p$model))
      stop("Provider artifact identities must match pinned build/model")
    allo_keys(x$license, c("status", "terms"), label = "provider license provenance")
    if (!identical(x$license$status, "authorized_external_precomputed")) stop("Provider license/provision status required")
    allo_text(x$license$terms, "license terms")
    allo_keys(x$execution, c("origin", "source", "provider_name", "provider_version", "command"), label = "provider execution provenance")
    if (!identical(x$execution$origin, "imported_external_execution") ||
      !identical(x$execution$provider_name, p$name) || !identical(x$execution$provider_version, p$version))
      stop("Provider execution provenance conflicts with pinned provider")
    for (k in c("source", "command")) allo_text(x$execution[[k]], paste("execution", k))
  }
  if ("metadata" %in% names(x) && (!is.list(x$metadata) || is.null(names(x$metadata))))
    stop("Physical provider metadata must be an object")
  canonical_bytes(x)
  x
}
indirect_provider <- function(m) {
  if (is.null(m$provenance)) stop("Structured provider provenance is required")
  allo_keys(m, c("version", "provider", "rows", "provenance"), label = "indirect binding import")
  if (m$version != "indirect-binding-import-1.0.0") stop("Unsupported indirect binding import version")
  if ("provenance" %in% names(m) && (!is.list(m$provenance) || is.null(names(m$provenance))))
    stop("Provider provenance must be an object")
  p <- m$provider
  allo_keys(p, c("name", "version", "build", "model", "task", "prediction_mode", "species", "arguments",
    "supported_presenters", "supported_lengths", "calibration", "batching", "evidence_kind"), label = "pinned predictor")
  for (k in c("name", "version", "build", "model", "task", "prediction_mode", "species", "calibration", "batching", "evidence_kind"))
    p[[k]] <- allo_text(p[[k]], k)
  if (p$name != "NetMHCIIpan" || !grepl("^3\\.2($|\\.)", p$version) ||
    !p$task %in% peptide_prediction_tasks() || p$task != "mhc2_binding" ||
    p$prediction_mode != "BA" || p$species != "NCBITaxon:9606")
    stop("Strict profile requires pinned NetMHCIIpan 3.2 binding-affinity evidence; other engines require a separately authorized profile")
  if (!p$evidence_kind %in% c("precomputed_predictor_output", "synthetic_fixture")) stop("Unknown imported evidence kind")
  p$arguments <- indirect_provider_arguments(p$arguments)
  provenance <- indirect_provider_provenance(m$provenance, p)
  p$supported_presenters <- allo_strings(p$supported_presenters, "provider presenter scope", TRUE)
  for (s in p$supported_presenters) {
    v <- peptide_presenters(list(s))[[1L]]
    if (v$status != "RESOLVED" || v$hla_class != "II" || !identical(v$mhc_key, s)) stop("Provider presenter must be canonical Class II pair")
  }
  if (!is.list(p$supported_lengths) || !is.null(names(p$supported_lengths)) || !length(p$supported_lengths)) stop("Provider length scope must be an array")
  p$supported_lengths <- as.list(sort(unique(vapply(p$supported_lengths, function(n) {
    if (!is.numeric(n) || length(n) != 1L || is.na(n) || n != floor(n) || n < 1L) stop("Invalid supported peptide length")
    as.integer(n)
  }, 1L))))
  config_id <- allo_id("policy", list(profile = indirect_profile(), provider = p))
  rows <- lapply(allo_records(m$rows), function(r) {
    allo_keys(r, c("peptide_sequence", "presenter", "status", "affinity_nM", "percentile_rank", "cores", "raw_output", "reason_codes"),
      "provenance", "imported binding row", nullable = c("affinity_nM", "percentile_rank"))
    if (!is.character(r$peptide_sequence) || length(r$peptide_sequence) != 1L ||
      !grepl("^[ACDEFGHIKLMNPQRSTVWY]{15}$", r$peptide_sequence)) stop("Provider query must be a literal 15-mer")
    v <- peptide_presenters(list(r$presenter))[[1L]]
    if (v$status != "RESOLVED" || !identical(v$mhc_key, r$presenter)) stop("Provider row requires canonical presenter")
    if (!r$status %in% c("evaluated", "unsupported", "unavailable", "failed", "ambiguous", "not_applicable")) stop("Invalid binding provider status")
    r$reason_codes <- allo_strings(r$reason_codes, "provider reasons")
    for (k in c("affinity_nM", "percentile_rank")) if (!is.null(r[[k]])) {
      a <- r[[k]]
      if (!is.numeric(a) || length(a) != 1L || !is.finite(a) || a < 0 || (k == "percentile_rank" && a > 100)) stop("Invalid imported provider metric")
      r[[k]] <- as.numeric(a)
    }
    if (r$status == "evaluated" && (is.null(r$affinity_nM) || !r$presenter %in% unlist(p$supported_presenters) ||
      !15L %in% unlist(p$supported_lengths))) stop("Evaluated binding row lacks supported affinity query")
    if (r$status != "evaluated" && (!is.null(r$affinity_nM) || length(r$cores))) stop("Unevaluable row cannot fabricate affinity/core evidence")
    if (!is.list(r$raw_output) || is.null(names(r$raw_output))) stop("Raw provider output must be an object")
    canonical_bytes(r$raw_output)
    r$cores <- allo_records(lapply(allo_records(r$cores), function(c) {
      allo_keys(c, c("core_sequence", "offset", "orientation", "register", "primary", "status"), label = "core assignment")
      if (!is.character(c$core_sequence) || length(c$core_sequence) != 1L || !grepl("^[ACDEFGHIKLMNPQRSTVWY]{9}$", c$core_sequence)) stop("Class II core must be nine literal residues")
      if (!is.numeric(c$offset) || length(c$offset) != 1L || is.na(c$offset) || c$offset != floor(c$offset) || c$offset < 0L || c$offset > 6L) stop("Invalid zero-based provider core offset")
      c$offset <- as.integer(c$offset); c$register <- allo_text(c$register, "register")
      c$primary <- allo_boolean(c$primary, "primary core")
      if (!c$orientation %in% c("forward", "reverse", "unknown") || !c$status %in% c("evaluated", "ambiguous", "unresolved_mapping", "unsupported")) stop("Invalid core assignment state")
      if (c$status == "evaluated" && c$orientation == "forward" &&
        substring(r$peptide_sequence, c$offset + 1L, c$offset + 9L) != c$core_sequence) stop("Core does not map to full ligand")
      c
    }))
    if ("provenance" %in% names(r)) {
      if (!is.list(r$provenance) || is.null(names(r$provenance))) stop("Row provenance must be an object")
      canonical_bytes(r$provenance)
    }
    r
  })
  # Identical duplicate evidence is normalized. Contradictory rows cannot be
  # selected by order. Physical provenance may differ without new evidence.
  queries <- vapply(rows, function(r) canonical_bytes(r[c("peptide_sequence", "presenter")]), "")
  normalized <- lapply(sort(unique(queries), method = "radix"), function(q) {
    rr <- rows[queries == q]
    sem <- lapply(rr, function(r) r[setdiff(names(r), c("provenance", "raw_output"))])
    if (length(unique(vapply(sem, canonical_bytes, ""))) != 1L) stop("Contradictory duplicate provider query")
    r <- sem[[1L]]
    r$evidence_id <- indirect_id("indirect-evidence-id-v1", list(target = r[c("peptide_sequence", "presenter")],
      dimension = "binding_affinity_and_core", provider_configuration_id = config_id,
      assay_occasion_id = NULL, result = r))
    r$raw_output_alternatives <- allo_records(lapply(rr, `[[`, "raw_output"))
    r$raw_output <- r$raw_output_alternatives[[1L]]
    r$provenance <- allo_records(lapply(rr, function(z) if (is.null(z$provenance)) list() else z$provenance))
    r
  })
  list(schema_version = "indirect-binding-import-1.0.0", semantics = p,
    provider_configuration_id = config_id, rows = allo_records(normalized),
    provenance = provenance)
}
indirect_query <- function(provider, peptide, presenter) {
  rows <- Filter(function(r) r$peptide_sequence == peptide && r$presenter == presenter, provider$rows)
  if (length(rows)) return(rows[[1L]])
  list(peptide_sequence = peptide, presenter = presenter,
    status = if (presenter %in% unlist(provider$semantics$supported_presenters) &&
      15L %in% unlist(provider$semantics$supported_lengths)) "unavailable" else "unsupported",
    affinity_nM = NULL, percentile_rank = NULL, cores = list(), raw_output = list(),
    raw_output_alternatives = list(),
    evidence_id = NULL, reason_codes = list("provider_query_unevaluable"), provenance = list())
}
indirect_primary <- function(row) {
  primary <- Filter(function(c) c$primary, row$cores)
  if (length(primary) != 1L) return(list(status = "ambiguous", core = NULL, reason_codes = list("primary_core_ambiguous")))
  c <- primary[[1L]]
  if (c$status != "evaluated" || c$orientation != "forward") return(list(status = "ambiguous", core = NULL,
    reason_codes = list("primary_core_unresolved")))
  list(status = "evaluated", core = c, reason_codes = list())
}
indirect_presenter_semantics <- function(p) {
  out <- lapply(p, function(v) allo_text(v, "presenter declaration"))
  v <- peptide_presenters(list(out$presenter))[[1L]]
  out$presenter <- if (v$status == "RESOLVED") v$mhc_key else v$canonical_input
  out
}
indirect_presenter_scope <- function(foundation, context, presenters) {
  universe <- Filter(function(u) u$semantics$scenario == context$allo_context_id, foundation$self_universes)
  scenario <- if (length(universe) == 1L) universe[[1L]]$typing_scenario_id else NULL
  calls <- Filter(function(t) t$subject_id == context$semantics$presenter_subject_id &&
    identical(t$typing_scenario_id, scenario) && t$locus %in% c("DRA", "DRB1"), foundation$typing)
  exact <- Filter(function(t) t$complete && t$evaluability$status == "evaluated" &&
    !is.null(t$resolved_allele), calls)
  alpha <- allo_set(lapply(Filter(function(t) t$locus == "DRA", exact), `[[`, "resolved_allele"))
  beta <- allo_set(lapply(Filter(function(t) t$locus == "DRB1", exact), `[[`, "resolved_allele"))
  declared <- Filter(function(p) p$status == "evaluated", presenters)
  molecules <- lapply(declared, function(p) peptide_presenters(list(p$presenter))[[1L]])
  declared_alpha <- allo_set(lapply(molecules, `[[`, "hla_alpha"))
  declared_beta <- allo_set(lapply(molecules, `[[`, "hla_beta"))
  missing <- allo_set(as.list(setdiff(unlist(beta), unlist(declared_beta))))
  reasons <- list()
  if (is.null(scenario) || !length(calls) || length(exact) != length(calls) ||
    !any(vapply(exact, function(t) t$locus == "DRB1", TRUE)) || !length(alpha))
    reasons <- c(reasons, list("presenter_scope_typing_insufficient"))
  # Phase A has chain typing, not a chain-pairing contract. A unique exactly
  # resolved alpha can be checked against explicit pairs; multiple alpha chains
  # cannot establish a complete scope by taking their Cartesian product.
  if (length(alpha) > 1L || anyDuplicated(unlist(lapply(molecules, `[[`, "hla_beta"))))
    reasons <- c(reasons, list("presenter_scope_pairing_ambiguous"))
  if (length(alpha) != 1L || !identical(alpha, declared_alpha) || length(missing) || !length(beta))
    reasons <- c(reasons, list("presenter_scope_declaration_incomplete"))
  complete <- !length(reasons)
  if (!complete) reasons <- c(reasons, list("presenter_scope_incomplete"))
  list(typing_scenario_id = scenario, complete = complete,
    status = if (complete) "evaluated" else "insufficient", evaluated_alpha_alleles = alpha,
    evaluated_beta_alleles = beta, undeclared_beta_alleles = missing, reason_codes = allo_set(reasons))
}
indirect_presenters <- function(config, foundation, context) {
  own <- Filter(function(p) p$scenario_id == context$semantics$scenario_id, config$presenters)
  out <- lapply(own, function(p) {
    if (p$subject_id != context$semantics$presenter_subject_id) stop("Presenter subject conflicts with allo context")
    v <- peptide_presenters(list(p$presenter))[[1L]]
    supported <- v$status == "RESOLVED" && v$hla_class == "II" &&
      startsWith(v$hla_alpha, "HLA-DRA*") && startsWith(v$hla_beta, "HLA-DRB1*")
    ready <- supported
    reasons <- if (supported) list() else list("presenter_out_of_profile_scope")
    if (supported) {
      universe <- Filter(function(u) u$semantics$scenario == context$allo_context_id, foundation$self_universes)
      scenario <- if (length(universe)) universe[[1L]]$typing_scenario_id else NULL
      for (chain in c(v$hla_alpha, v$hla_beta)) {
        calls <- Filter(function(t) t$subject_id == p$subject_id && identical(t$typing_scenario_id, scenario) &&
          identical(t$resolved_allele, chain) && t$evaluability$status == "evaluated", foundation$typing)
        if (!length(calls)) { ready <- FALSE; reasons <- c(reasons, list("presenter_resolution_insufficient")) }
      }
    }
    list(presenter = if (v$status == "RESOLVED") v$mhc_key else v$canonical_input,
      subject_id = p$subject_id, resolution_policy = p$resolution_policy,
      status = if (ready) "evaluated" else if (!supported) "unsupported" else "insufficient",
      reason_codes = allo_set(reasons))
  })
  allo_records(out)
}
