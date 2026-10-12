# T013 provider declarations describe implemented adapters. Availability and
# supported local models are detected separately for each execution request.
predictor_adapter_version <- function() "1.0.1"
predictor_evidence_schema_version <- function() "1.0.1"

predictor_registry <- function() list(
  netmhcpan = list(name = "NetMHCpan", executable = "netMHCpan", hla_class = "I",
    tasks = c("mhc1_binding", "mhc1_presentation"), lengths = 8:14,
    alphabet = "ACDEFGHIKLMNPQRSTVWY", allele_types = c("A", "B", "C"),
    metrics = c("BA_nM", "BA_rank", "EL_score", "EL_rank"),
    returns_core = FALSE, returns_stability = FALSE,
    requires_models = FALSE, license = "external_DTU_license", implemented = TRUE),
  netmhciipan = list(name = "NetMHCIIpan", executable = "netMHCIIpan", hla_class = "II",
    tasks = c("mhc2_binding", "mhc2_presentation"), lengths = 9:25,
    alphabet = "ACDEFGHIKLMNPQRSTVWY", allele_types = c("DR", "DQ", "DP"),
    metrics = c("BA_nM", "BA_rank", "EL_score", "EL_rank"),
    returns_core = TRUE, returns_stability = FALSE,
    requires_models = FALSE, license = "external_DTU_license", implemented = TRUE),
  mhcflurry = list(name = "MHCflurry", executable = "mhcflurry-predict", hla_class = "I",
    tasks = c("mhc1_binding", "mhc1_presentation"), lengths = 8:15,
    alphabet = "ACDEFGHIKLMNPQRSTVWY", allele_types = c("A", "B", "C"),
    metrics = c("BA_nM", "BA_rank", "presentation_score", "presentation_percentile"),
    returns_core = FALSE, returns_stability = FALSE,
    requires_models = TRUE, license = "external_open_source_models", implemented = TRUE),
  fixture = list(name = "fixture", executable = NULL, hla_class = c("I", "II"),
    tasks = c("mhc1_binding", "mhc1_presentation", "mhc2_binding", "mhc2_presentation"),
    lengths = 1:30, alphabet = "ACDEFGHIKLMNPQRSTVWYU", allele_types = c("A", "B", "C", "DR", "DQ", "DP"),
    metrics = c("BA_nM", "BA_rank", "EL_score", "EL_rank"),
    returns_core = TRUE, returns_stability = FALSE, requires_models = FALSE,
    license = "internal_test_only", implemented = TRUE))

predictor_version_from_text <- function(text, provider) {
  pattern <- if (provider == "netmhcpan") "NetMHCpan[^0-9]*([0-9]+\\.[0-9]+(?:\\.[0-9]+)?(?:[A-Za-z][A-Za-z0-9._-]*)?)" else
    if (provider == "netmhciipan") "NetMHCIIpan[^0-9]*([0-9]+\\.[0-9]+[a-z]?)" else
      "[Mm][Hh][Cc]flurry[^0-9]*([0-9]+\\.[0-9]+(\\.[0-9]+)?[A-Za-z0-9.]*)"
  m <- regexec(pattern, paste(text, collapse = "\n"), perl = TRUE)
  hit <- regmatches(paste(text, collapse = "\n"), m)[[1L]]
  if (length(hit) < 2L) NULL else hit[[2L]]
}

predictor_capability <- function(provider, executable = NULL, model_path = NULL) {
  spec <- predictor_registry()[[provider]]
  if (is.null(spec)) stop("Unknown predictor adapter")
  if (provider == "fixture") return(c(spec, list(available = TRUE,
    availability_status = "available", executable_path = NULL,
    executable_checksum = NULL, predictor_version = "fixture-1",
      version_status = "detected", model_version = "fixture-model-1",
      model_version_status = "detected",
    model_identity = "fixture-model-1", adapter_version = predictor_adapter_version())))
  path <- if (is.null(executable)) Sys.which(spec$executable) else executable
  if (!length(path) || !nzchar(path) || !file.exists(path) || dir.exists(path)) path <- NULL
  if (!is.null(path)) path <- normalizePath(path, mustWork = TRUE)
  # Installed binaries can still lack a runnable model/resource. No model is
  # inferred from the adapter declaration or from a package version.
  model_identity <- predictor_model_identity(model_path)
  model_ok <- !spec$requires_models || (!is.null(model_identity) &&
    length(model_identity$files)>0L)
  version <- NULL
  if (!is.null(path)) {
    flag <- if (provider == "mhcflurry") "--version" else "-h"
    output <- tryCatch(suppressWarnings(system2(path, flag, stdout = TRUE, stderr = TRUE)),
      error = function(e) character())
    version <- predictor_version_from_text(output, provider)
  }
  supported_version <- is.null(version) || if (provider == "netmhcpan")
    grepl("^4\\.[12](?:\\.[0-9]+)?(?:[A-Za-z][A-Za-z0-9._-]*)?$",version) else if (provider == "netmhciipan")
      startsWith(version,"4.3") else startsWith(version,"2.")
  c(spec, list(available = !is.null(path) && model_ok && !is.null(version) && supported_version,
    availability_status = if (is.null(path)) "predictor_unavailable" else
      if (is.null(version)) "version_unresolved" else
      if (!supported_version) "unsupported_version" else
      if (!model_ok) "model_unavailable" else "available",
    executable_path = path,
    executable_checksum = if (is.null(path)) NULL else
      digest::digest(file = path, algo = "sha256", serialize = FALSE),
    predictor_version = version, version_status = if (is.null(version)) "unknown" else "detected",
    model_version = NULL, model_version_status = "unknown",
    model_identity = model_identity,
    adapter_version = predictor_adapter_version()))
}

predictor_allele <- function(candidate, provider) {
  if (provider == "fixture") return(candidate$mhc_key)
  if (provider %in% c("netmhcpan", "mhcflurry")) {
    if (candidate$hla_class != "I" || is.null(candidate$hla_allele)) return(NULL)
    allele <- candidate$hla_allele
    if (!grepl("^HLA-[ABC]\\*[0-9]{2,3}:[0-9]{2,3}$", allele)) return(NULL)
    if (provider == "netmhcpan") return(sub("\\*", "", allele))
    return(allele)
  }
  if (candidate$hla_class != "II" || is.null(candidate$hla_alpha) ||
      is.null(candidate$hla_beta)) return(NULL)
  alpha <- candidate$hla_alpha; beta <- candidate$hla_beta
  if (!grepl("^HLA-(DRA|DQA1|DPA1)\\*[0-9]{2}:[0-9]{2}$", alpha) ||
      !grepl("^HLA-(DRB[1345]|DQB1|DPB1)\\*[0-9]{2}:[0-9]{2}$", beta)) return(NULL)
  compact <- function(x) gsub("[:*]", "", sub("^HLA-", "", x))
  if (startsWith(beta, "HLA-DRB")) {
    if (alpha != "HLA-DRA*01:01") return(NULL)
    return(sub("^(DRB[1345])([0-9]{4})$", "\\1_\\2", compact(beta)))
  }
  paste0("HLA-", compact(alpha), "-", compact(beta))
}

predictor_candidate_status <- function(candidate, capability, provider, task) {
  if (!task %in% capability$tasks) return(list(status = "unsupported", reason = "unsupported_task"))
  if (!candidate$hla_class %in% capability$hla_class)
    return(list(status = "unsupported", reason = "unsupported_presenter"))
  if ((startsWith(task,"mhc1_") && candidate$hla_class!="I") ||
      (startsWith(task,"mhc2_") && candidate$hla_class!="II"))
    return(list(status = "unsupported", reason = "unsupported_task"))
  n <- nchar(candidate$ligand_sequence)
  if (!n %in% capability$lengths) return(list(status = "unsupported", reason = "unsupported_length"))
  residues <- strsplit(candidate$ligand_sequence, "", fixed = TRUE)[[1L]]
  bad <- sort(unique(residues[!residues %in% strsplit(capability$alphabet, "", fixed = TRUE)[[1L]]]),
    method = "radix")
  if (length(bad)) return(list(status = "unsupported", reason = "unsupported_residue",
    residue = paste(bad, collapse = "")))
  external <- predictor_allele(candidate, provider)
  if (is.null(external)) return(list(status = "unsupported", reason = "unsupported_allele_format"))
  if (!capability$available) return(list(status = "unavailable",
    reason = capability$availability_status))
  list(status = "ready", reason = NULL, external_allele = external)
}
