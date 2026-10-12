sample_records <- function(samples) {
  if (is.null(samples)) return(list())
  if (is.data.frame(samples)) samples <- lapply(seq_len(nrow(samples)), function(i)
    lapply(samples, function(col) col[[i]]))
  if (!is.list(samples) || !is.null(names(samples))) stop("samples must be a data.frame or record array")
  for (s in samples) {
    strict_keys(s, c("sample_id", "patient_id", "timepoint", "tissue", "role"), "sample")
    scalar_text(s$sample_id, "sample_id"); scalar_text(s$patient_id, "patient_id")
    for (n in names(s)) scalar_text(s[[n]], n)
  }
  if (anyDuplicated(vapply(samples, `[[`, "", "sample_id"))) stop("Duplicate sample_id")
  sorted_records(samples)
}

input_records <- function(inputs, config, samples) {
  if (!is.list(inputs) || !is.null(names(inputs))) stop("inputs must be a record array")
  out <- lapply(inputs, function(x) {
    strict_keys(x, c("input_id", "input_role", "sample_id", "path", "checksum", "version"), "input")
    scalar_text(x$input_role, "input_role")
    if (!grepl("^[a-z][a-z0-9_]*$", x$input_role)) stop("Invalid input_role")
    scalar_text(x$input_id, "input_id"); scalar_text(x$sample_id, "sample_id")
    if (!x$sample_id %in% vapply(samples, `[[`, "", "sample_id")) stop("Input sample binding missing")
    spec <- resource_spec(list(path = x$path, required = TRUE, checksum = x$checksum, version = x$version))
    path <- resolve_path(spec$path, config$paths$input_root)
    if (!file.exists(path) || dir.exists(path)) stop("Required input file missing")
    list(input_id = x$input_id, input_role = x$input_role, sample_id = x$sample_id, declared_path = x$path,
      resolved_path = path, version = x$version,
      integrity = integrity_record(path, config$checksum, spec$checksum))
  })
  if (anyDuplicated(vapply(out, `[[`, "", "input_id"))) stop("Duplicate input_id")
  sorted_records(out)
}

project_contexts <- function(contexts, samples, mode) {
  if (is.null(contexts)) {
    patients <- sort(unique(vapply(samples, `[[`, "", "patient_id")), method = "radix")
    contexts <- lapply(patients, function(p) list(context_id = paste0("somatic:", p), mode = mode,
      source_subject_id = p, source_subject_role = "tumor", self_subject_id = p,
      self_subject_role = "normal", presenter_subject_id = p, presenter_subject_role = "patient",
      sample_bindings = as.list(vapply(Filter(function(s) s$patient_id == p, samples), `[[`, "", "sample_id")),
      direction = "within_patient", organ = NULL, tissue = NULL, timepoint = NULL,
      scenario_id = "default", missing_reasons = list(organ = "NOT_PROVIDED", tissue = "NOT_PROVIDED", timepoint = "NOT_PROVIDED")))
  }
  if (!is.list(contexts) || !is.null(names(contexts))) stop("analysis_context must be a record array")
  c <- schema_contract()
  for (i in seq_along(contexts)) {
    x <- contexts[[i]]
    validate_record(x, c$records$analysis_context, c)
    if (x$mode != "somatic") stop("unsupported analysis mode")
    if (x$direction != "within_patient") stop("Invalid somatic direction; expected within_patient")
    if (!x$source_subject_role %in% c("tumor") || !x$self_subject_role %in% c("normal", "patient") ||
        x$presenter_subject_role != "patient" || x$source_subject_id != x$self_subject_id ||
        x$source_subject_id != x$presenter_subject_id) stop("Invalid somatic role binding")
    for (s in x$sample_bindings) {
      found <- Filter(function(a) a$sample_id == s && a$patient_id == x$source_subject_id, samples)
      if (length(found) != 1L) stop("Invalid context sample identity binding")
    }
    contexts[[i]]$sample_bindings <- as.list(sort(unlist(x$sample_bindings), method = "radix"))
  }
  if (anyDuplicated(vapply(contexts, `[[`, "", "context_id"))) stop("Duplicate context_id")
  sorted_records(contexts)
}

identity_records <- function(records) lapply(records, function(r) {
  # Paths are fallback identities only when no portable content declaration exists.
  checksum <- r$integrity$declared
  observed <- r$integrity$observed
  content <- if (!is.null(checksum)) checksum else if (observed$algorithm == "sha256") observed else
    list(path = r$resolved_path, observed = observed)
  if (!is.null(r$input_id)) return(list(sample_id = r$sample_id,
    input_role = r$input_role, version = r$version, content = content))
  list(id = r$logical_name,
       sample_id = r$sample_id, resource_id = r$resource_id, version = r$version,
       type = r$type, required = r$required, exists = r$integrity$exists, content = content)
})

neo_project <- function(samples = NULL, hla = NULL, reference = NULL, config = NULL,
                        path = NULL, mode = "somatic", project_id = "project",
                        analysis_context = NULL, inputs = list(), initialize = FALSE,
                        fixture_kind = "real") {
  if (!identical(mode, "somatic")) stop("unsupported analysis mode; allo execution is not implemented")
  if (!is.logical(initialize) || length(initialize) != 1L || is.na(initialize)) stop("Invalid initialize")
  if (!is.null(path)) {
    base <- if (is.null(config)) NULL else config$paths$run_root
    resolved <- resolve_path(path, base)
    if (dir.exists(resolved)) return(resume_project(resolved,
      samples = if (missing(samples)) NULL else samples, hla = hla, reference = reference,
      config = config, analysis_context = analysis_context,
      inputs = if (missing(inputs)) NULL else inputs))
    initialize <- TRUE
  }
  if (is.null(config)) config <- neo_config()
  if (!inherits(config, "neo_config")) stop("config must be created by neo_config")
  config <- do.call(neo_config, c(unclass(config), list(env = character())))
  if (!is.null(reference)) {
    strict_keys(reference, names(reference), "reference")
    if (length(reference)) names(reference) <- paste0("reference.", names(reference))
    reference <- validate_resources(reference)
    for (n in names(reference)) config$resources[[n]] <- reference[[n]]
  }
  scalar_text(project_id, "project_id")
  scalar_text(fixture_kind, "fixture_kind")
  if (!fixture_kind %in% c("real", "synthetic", "mock") || length(fixture_kind) != 1L) stop("Invalid fixture_kind")
  samples <- sample_records(samples)
  contexts <- project_contexts(analysis_context, samples, mode)
  if (!is.null(hla)) {
    strict_keys(hla, unique(vapply(samples, `[[`, "", "patient_id")), "hla subject")
    hla <- lapply(hla, function(x) {
      if (!is.character(x) || anyNA(x) || any(!nzchar(x))) stop("hla must contain declared allele strings")
      as.list(sort(unique(x), method = "radix"))
    })
  }
  resources <- resolve_resources(config)
  records <- input_records(inputs, config, samples)
  c <- schema_contract()
  scientific <- list(contract_version = infra_version(), schema_version = c$schema_version,
    id_schema_version = c$id_schema_version, scientific = config$scientific, hla = hla,
    contexts = contexts, resources = sorted_records(identity_records(resources)))
  if (!is.null(config$peptides)) scientific$peptides <- config$peptides
  identity <- list(config = scientific, samples = samples, inputs = sorted_records(identity_records(records)))
  p <- structure(list(project_id = project_id, analysis_context = contexts, samples = samples,
    config = config, resources = resources, inputs = records, hla = hla,
    analysis_id = paste0("analysis:", infra_hash(identity)), config_fingerprint = infra_hash(scientific),
    input_fingerprint = infra_hash(list(samples = samples, inputs = identity$inputs)),
    run_id = NULL, path = NULL, schema_version = c$schema_version,
    id_schema_version = c$id_schema_version, fixture_kind = fixture_kind), class = "neo_project")
  if (!is.null(config$allo)) p$allo_foundation <- allo_foundation(config$allo)
  if (!is.null(config$indirect_hla)) p$indirect_hla <- indirect_baseline(p$allo_foundation, config$indirect_hla)
  if (native_allo_enabled(p)) p$allo_native_evidence <- native_allo_build(p)
  if (initialize) p <- initialize_run(p, path)
  p
}

print.neo_project <- function(x, ...) {
  cat("<neo_project> ", length(x$samples), " samples; ", length(x$resources), " resources\n",
      "analysis_id: ", x$analysis_id, "\nrun_id: ", if (is.null(x$run_id)) "(not initialized)" else x$run_id,
      "\n", sep = "")
  invisible(x)
}
