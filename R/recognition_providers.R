# Independent registry: none of these metrics is a presentation/priority score.
recognition_registry <- function() list(
  prime = list(scope = "peptide_intrinsic", task = "recognition.peptide_intrinsic",
    classes = list("I"), lengths = as.list(8:14), mhc_conditioning = "present",
    dependencies = list("PRIME 2.1", "MixMHCpred >=3.0"), license = "user_license_required"),
  pmtnet = list(scope = "tcr_specific", task = "recognition.tcr_specific",
    classes = list("I"), lengths = NULL, mhc_conditioning = "present",
    dependencies = list("pMTnet", "configured Python environment", "model/HLA/background library"), license = "user_managed"),
  ergo_ii = list(scope = "tcr_specific", task = "recognition.tcr_specific",
    classes = list("I", "II"), lengths = NULL, mhc_conditioning = "configuration",
    dependencies = list("ERGO-II", "configured PyTorch environment", "selected model"), license = "user_managed"),
  nettcr_2_0 = list(scope = "tcr_specific", task = "recognition.tcr_specific",
    classes = list("I"), lengths = list(9L), presenters = list("HLA-A*02:01"),
    mhc_conditioning = "required_context", dependencies = list("NetTCR-2.0", "configured TensorFlow environment", "selected model"),
    license = "academic_noncommercial_or_separate_license"),
  fixture = list(scope = "tcr_specific", task = "recognition.tcr_specific",
    classes = list("I", "II"), lengths = NULL, mhc_conditioning = "present",
    dependencies = list(), license = "synthetic_fixture_only"))
recognition_provider_config <- function(x) {
  strict_keys(x, c("provider", "provider_version", "model_version", "model", "features", "use_mhc",
    "options", "model_files", "background_files", "supported_presenters", "output_input_id",
    "executable", "dependencies", "license_status", "fixture", "scope"), "recognition provider")
  if (!x$provider %in% names(recognition_registry())) stop("Unknown recognition provider")
  if (is.null(x$model)) x$model <- if (x$provider == "prime") "prime-2.1" else "beta_only"
  if (!x$model %in% c("prime-2.1", "beta_only", "paired", "beta_vj", "paired_vj", "custom"))
    stop("Unsupported recognition feature configuration")
  if (x$provider == "prime" && x$model != "prime-2.1") stop("PRIME model must be prime-2.1")
  if (x$provider == "pmtnet" && x$model != "beta_only") stop("pMTnet interface consumes beta CDR3 only")
  if (x$provider == "nettcr_2_0" && !x$model %in% c("beta_only", "paired")) stop("Unsupported NetTCR-2.0 variant")
  required <- switch(x$model, `prime-2.1` = character(), beta_only = "beta_cdr3_aa",
    paired = c("alpha_cdr3_aa", "beta_cdr3_aa"),
    beta_vj = c("beta_cdr3_aa", "beta_v_call", "beta_j_call"),
    paired_vj = c("alpha_cdr3_aa", "beta_cdr3_aa", "alpha_v_call", "alpha_j_call", "beta_v_call", "beta_j_call"),
    custom = unlist(x$features))
  if (!is.null(x$features) && x$model != "custom" && !setequal(unlist(x$features), required))
    stop("Feature declaration conflicts with selected model")
  if (x$model == "custom" && (x$provider != "ergo_ii" || !"beta_cdr3_aa" %in% required))
    stop("Custom ERGO-II configuration requires beta CDR3")
  allowed <- c("alpha_cdr3_aa", "beta_cdr3_aa", "alpha_v_call", "alpha_j_call",
    "beta_v_call", "beta_j_call", "t_cell_type")
  if (any(!required %in% allowed)) stop("Unknown receptor feature")
  if (x$provider == "fixture" && identical(x$scope, "peptide_intrinsic")) required <- character()
  x$features <- as.list(sort(unique(required), method = "radix"))
  if (is.null(x$use_mhc)) x$use_mhc <- x$provider != "ergo_ii"
  if (!is.logical(x$use_mhc) || length(x$use_mhc) != 1L || is.na(x$use_mhc)) stop("Invalid MHC configuration")
  if (x$provider %in% c("prime", "pmtnet", "nettcr_2_0") && !x$use_mhc) stop("Provider requires HLA context")
  if (is.null(x$options)) x$options <- list()
  if (!is.list(x$options)) stop("Provider options must be an object")
  if (is.null(x$fixture)) x$fixture <- FALSE
  if (!is.logical(x$fixture) || length(x$fixture) != 1L || is.na(x$fixture)) stop("Invalid fixture flag")
  if (x$provider == "fixture" && !x$fixture) stop("Fixture provider must be explicitly marked synthetic")
  if (!is.null(x$scope) && !(x$provider == "fixture" && x$scope %in% c("peptide_intrinsic", "tcr_specific")))
    stop("Scope is fixed for named providers")
  x
}
recognition_file_identity <- function(paths, base = NULL) lapply(paths, function(p) {
  path <- resolve_path(p, base)
  list(path = path, available = file.exists(path) && !dir.exists(path), sha256 =
    if (file.exists(path) && !dir.exists(path)) digest::digest(file = path, algo = "sha256", serialize = FALSE) else NULL)
})
recognition_capability <- function(config, base = NULL) {
  x <- recognition_provider_config(config); spec <- recognition_registry()[[x$provider]]
  if (!is.null(x$scope)) spec$scope <- x$scope
  spec$task <- paste0("recognition.", spec$scope)
  spec$mhc_conditioning <- if (!x$use_mhc) "absent" else if (x$provider == "nettcr_2_0") "required_context" else "present"
  models <- recognition_file_identity(x$model_files, base)
  background <- recognition_file_identity(x$background_files, base)
  dependencies <- recognition_file_identity(x$dependencies, base)
  executable <- recognition_file_identity(if (is.null(x$executable)) list() else list(x$executable), base)
  # Passive file detection only. No probing commands or network side effects.
  declared_version <- recognition_nullable(x$provider_version)
  version_compatible <- is.null(declared_version) || identical(declared_version, "unknown") ||
    switch(x$provider, prime = grepl("^2\\.1($|[.+-])", declared_version),
      nettcr_2_0 = grepl("^2\\.0($|[.+-])", declared_version), TRUE)
  files <- c(models, background, dependencies, executable)
  ready <- length(executable) > 0L && length(models) > 0L && length(dependencies) > 0L &&
    all(vapply(files, function(f) isTRUE(f$available), TRUE))
  c(spec, list(provider = x$provider, provider_version = x$provider_version,
    version_compatible = version_compatible,
    version_status = if (is.null(x$provider_version) || identical(x$provider_version, "unknown")) "unknown" else "declared",
    model = x$model, model_version = x$model_version, required_features = x$features,
    parser_available = TRUE, dependency_availability = files,
    local_dependencies_detected = ready, dependency_probe_status = "configured_files_only/runtime_not_probed",
    implementation_identity = c(executable, dependencies), execution_available = FALSE,
    execution_status = "not_implemented_in_this_phase", license_status = if (is.null(x$license_status))
      "not_verified" else x$license_status, model_identity = models, background_identity = background,
    applicability = if (x$provider %in% c("nettcr_2_0", "prime")) "provider_scope_limited" else "unknown_generalization",
    config = x))
}
recognition_provider_hla <- function(candidate, provider) {
  p <- recognition_presenter(candidate$mhc_key)
  if (p$status != "known" || !identical(p$mhc_key, candidate$mhc_key)) return(NULL)
  if (provider == "prime") {
    if (candidate$hla_class != "I" || !grepl("^HLA-[ABC]\\*[0-9]{2,3}:[0-9]{2,3}$", p$mhc_key)) return(NULL)
    return(gsub("[:*]", "", sub("^HLA-", "", p$mhc_key)))
  }
  if (provider == "pmtnet") {
    if (candidate$hla_class != "I" || !grepl("^HLA-[ABC]\\*[0-9]{2,3}:[0-9]{2,3}$", p$mhc_key)) return(NULL)
    return(sub("\\*", "", p$mhc_key))
  }
  p$mhc_key
}
recognition_input_status <- function(candidate, profile, cap) {
  if (!isTRUE(cap$version_compatible)) return("unsupported_provider_version")
  if (!candidate$hla_class %in% unlist(cap$classes)) return("unsupported_class")
  if (!is.null(cap$lengths) && !nchar(candidate$ligand_sequence) %in% unlist(cap$lengths)) return("unsupported_length")
  if (!is.null(cap$presenters) && !candidate$mhc_key %in% unlist(cap$presenters)) return("unsupported_presenter")
  if (!is.null(cap$config$supported_presenters) && !candidate$mhc_key %in% unlist(cap$config$supported_presenters))
    return("unsupported_presenter")
  if (recognition_sequence_status(candidate$ligand_sequence) != "valid") return("unsupported_peptide_character")
  features <- unlist(cap$required_features)
  if (length(features) && (is.null(profile) || !all(features %in% names(profile$features)))) return("incomplete_features")
  for (f in features) {
    value <- profile$features[[f]]
    if (grepl("_aa$", f) && recognition_sequence_status(value) != "valid") return("unsupported_tcr_character")
    if (is.null(recognition_nullable(value))) return("incomplete_features")
  }
  if (cap$mhc_conditioning != "absent" && is.null(recognition_provider_hla(candidate, cap$provider)))
    return("unsupported_presenter")
  "supported"
}
recognition_functionality <- function(observations) {
  states <- sort(unique(vapply(observations, function(o) if (is.null(o$functional_status))
    "functionality_unknown" else o$functional_status, "")), method = "radix")
  status <- if ("functional" %in% states) "functional" else if (!length(states) ||
    all(states == "functionality_unknown")) "functionality_unknown" else if (
    all(states == "explicitly_nonfunctional")) "explicitly_nonfunctional" else "functionality_conflicting"
  list(status = status, observation_statuses = as.list(states))
}
recognition_query <- function(candidate, profile, cap, observations = list()) {
  if (cap$scope == "peptide_intrinsic") profile <- NULL
  required <- unlist(cap$required_features)
  used <- if (is.null(profile)) list() else profile$features[intersect(required, names(profile$features))]
  payload <- list(candidate_id = candidate$candidate_id, peptide = candidate$ligand_sequence,
    receptor_feature_profile_id = if (is.null(profile)) NULL else profile$receptor_feature_profile_id,
    presenter = if (cap$mhc_conditioning == "absent") NULL else candidate$mhc_key,
    provider = cap$provider, provider_version = cap$provider_version, model_version = cap$model_version,
    model = cap$model, scope = cap$scope, features_used = used, required_features = cap$required_features,
    mhc_conditioning = cap$mhc_conditioning, options = cap$config$options,
    model_identity = cap$model_identity, implementation_identity = cap$implementation_identity,
    background_identity = cap$background_identity)
  functionality <- if (cap$scope == "peptide_intrinsic") NULL else recognition_functionality(observations)
  status <- recognition_input_status(candidate, profile, cap)
  if (!is.null(functionality) && functionality$status == "explicitly_nonfunctional") status <- "nonfunctional_receptor"
  if (!is.null(functionality) && functionality$status == "functionality_conflicting") status <- "receptor_functionality_unresolved"
  c(list(recognition_query_id = recognition_id("query", payload), fingerprint = infra_hash(list(contract_version = recognition_version(), payload = payload))), payload,
    list(sample_id = candidate$sample_id, input_status = status, receptor_functionality = functionality,
      input_completeness = if (is.null(profile)) "no_receptor_input" else profile$completeness,
      provider_hla = if (cap$mhc_conditioning == "absent") NULL else recognition_provider_hla(candidate, cap$provider)))
}
recognition_prepare <- function(candidates, repertoire, cap, settings) {
  if (!is.null(settings$candidate_ids)) {
    if (any(!unlist(settings$candidate_ids) %in% vapply(candidates, `[[`, "", "candidate_id"))) stop("Unknown candidate subset")
    candidates <- Filter(function(c) c$candidate_id %in% unlist(settings$candidate_ids), candidates)
  }
  profiles <- repertoire$profiles
  if (!is.null(settings$profile_ids)) {
    if (any(!unlist(settings$profile_ids) %in% vapply(profiles, `[[`, "", "receptor_feature_profile_id"))) stop("Unknown receptor subset")
    profiles <- Filter(function(p) p$receptor_feature_profile_id %in% unlist(settings$profile_ids), profiles)
  }
  profile_ids <- vapply(profiles, `[[`, "", "receptor_feature_profile_id")
  by_sample <- split(repertoire$observations, vapply(repertoire$observations, `[[`, "", "sample_id"))
  selected_ids <- lapply(by_sample, function(observations) intersect(profile_ids,
    unique(vapply(observations, `[[`, "", "receptor_feature_profile_id"))))
  counts <- vapply(candidates, function(c) if (cap$scope == "peptide_intrinsic") 1L else
    length(selected_ids[[c$sample_id]]), 0L)
  count <- sum(counts)
  if (count > settings$max_queries) return(list(status = "not_run", reason = "query_limit", requested = count,
    limit = settings$max_queries, candidate_ids = as.list(sort(vapply(candidates, `[[`, "", "candidate_id"))),
    queries = list(), batches = list()))
  queries <- list()
  for (c in candidates) {
    selected <- if (cap$scope == "peptide_intrinsic") list(NULL) else
      profiles[match(selected_ids[[c$sample_id]], profile_ids)]
    for (p in selected) {
      observations <- if (is.null(p)) list() else Filter(function(o)
        o$receptor_feature_profile_id == p$receptor_feature_profile_id, by_sample[[c$sample_id]])
      queries[[length(queries) + 1L]] <- recognition_query(c, p, cap, observations)
    }
  }
  queries <- recognition_unique(queries, "recognition_query_id")
  batches <- if (!length(queries)) list() else unname(split(queries,
    ceiling(seq_along(queries) / settings$batch_size)))
  list(status = if (length(queries)) "prepared" else "no_queries", reason = NULL,
    requested = count, queries = queries, batches = batches)
}
# Reviewable construction only; execution requires a future local runner.
recognition_prime_command <- function(config, input, output, alleles, base = NULL) {
  cap <- recognition_capability(config, base)
  if (cap$provider != "prime" || is.null(config$executable)) return(list(status = "unavailable", command = NULL))
  list(status = if (cap$local_dependencies_detected && identical(cap$license_status, "user_confirmed"))
    "configured_command_only" else "dependencies_or_license_unverified", command = resolve_path(config$executable, base),
    args = as.list(c("-i", input, "-o", output, "-a", paste(alleles, collapse = ","))), executed = FALSE)
}
# Local input rendering exposes the exact configured consumption contract.
recognition_provider_input <- function(queries, cap) {
  accepted <- recognition_sort(Filter(function(q) q$input_status == "supported", queries), "recognition_query_id")
  rejected <- Filter(function(q) q$input_status != "supported", queries)
  rows <- lapply(accepted, function(q) {
    f <- q$features_used
    if (cap$provider == "prime") return(list(Peptide = q$peptide, Allele = q$provider_hla))
    if (cap$provider == "pmtnet") return(list(CDR3 = f$beta_cdr3_aa, Antigen = q$peptide, HLA = q$provider_hla))
    if (cap$provider == "ergo_ii") return(list(TRA = f$alpha_cdr3_aa, TRB = f$beta_cdr3_aa,
      TRAV = f$alpha_v_call, TRAJ = f$alpha_j_call, TRBV = f$beta_v_call, TRBJ = f$beta_j_call,
      `T-Cell-Type` = f$t_cell_type, Peptide = q$peptide, MHC = q$provider_hla))
    if (cap$provider == "nettcr_2_0") return(list(peptide = q$peptide, CDR3a = f$alpha_cdr3_aa, CDR3b = f$beta_cdr3_aa))
    list(query_id = q$recognition_query_id, peptide = q$peptide)
  })
  list(rows = rows, recognition_query_ids = as.list(vapply(accepted, `[[`, "", "recognition_query_id")),
    rejected_query_ids = as.list(vapply(rejected, `[[`, "", "recognition_query_id")),
    consumed_features = cap$required_features, mhc_conditioning = cap$mhc_conditioning)
}
