# Opt-in diagnostics only: no project creation, provider process or scientific stage.
first_run_vcf_header <- function(path, sample_id) {
  con <- if (grepl("[.]gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  for (i in seq_len(10000L)) {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) break
    if (startsWith(line, "#CHROM")) {
      h <- strsplit(line, "\t", fixed = TRUE)[[1L]]
      if (length(h) < 8L || !identical(h[1:8], c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))) return("MALFORMED_HEADER")
      if (length(h) > 10L) return("MULTIPLE_SAMPLE_COLUMNS")
      if (length(h) == 8L) return("SITES_ONLY_EXPLICIT_BINDING")
      if (length(h) != 10L || h[[9L]] != "FORMAT") return("MALFORMED_HEADER")
      return(if (identical(h[[10L]], sample_id)) "SAMPLE_MATCH" else "SAMPLE_MISMATCH")
    }
    if (!startsWith(line, "#")) break
  }
  "MISSING_OR_OVERSIZED_HEADER"
}

first_run_preflight <- function(config, samples = NULL, hla = NULL, inputs = list(), output_dir = config$paths$run_root,
                                include_private = FALSE) {
  if (!inherits(config, "neo_config")) stop("config must be created by neo_config()")
  setup_flag(include_private, "include_private")
  checks <- list()
  add <- function(field, status, expected, observed, action = "") {
    checks[[length(checks) + 1L]] <<- data.frame(field = field, status = status,
      expected = expected, observed = observed, action = action, stringsAsFactors = FALSE)
  }
  # Keep setup API's fail-closed file/release/checksum checks; never expose its paths/messages by default.
  resources <- tryCatch(suppressWarnings(neo_check_resources(config = config,
    checksum = FALSE, validate_reference = FALSE)), error = function(e) NULL)
  if (is.null(resources)) add("config$resources", "FAIL", "compatible complete reference roles",
    "RESOURCE_CHECK_FAILED", "Inspect neo_check_resources(config=config) locally; replace no data automatically.") else {
    for (i in seq_len(nrow(resources$files))) {
      r <- resources$files[i, ]
      add(paste0("resources$", r$role), if (r$status == "OK") "PASS" else "FAIL",
        "readable role with compatible declared version/integrity", r$status,
        if (r$status == "OK") "" else "Inspect this role with neo_check_resources(config=config) locally; use the reference setup guide.")
    }
  }
  ss <- tryCatch(sample_records(samples), error = function(e) NULL)
  valid_samples <- !is.null(ss) && length(ss) > 0L
  add("samples", if (valid_samples) "PASS" else "FAIL", "unique sample_id and patient_id bindings",
    if (valid_samples) paste(length(ss), "declared samples") else "INVALID_OR_EMPTY",
    if (valid_samples) "" else "Inspect samples columns and duplicates; identifiers are not printed.")
  patients <- if (valid_samples) unique(vapply(ss, `[[`, "", "patient_id")) else character()
  hla_map <- is.list(hla) && !is.null(names(hla)) && !anyDuplicated(names(hla)) &&
    setequal(names(hla), patients) && length(patients) > 0L
  add("hla", if (hla_map) "PASS" else "FAIL", "one HLA declaration per declared patient",
    if (hla_map) "PATIENT_KEYS_MATCH" else "MISSING_OR_MISMATCHED_KEYS",
    if (hla_map) "" else "Key hla by the declared patient_id locally, not by sample_id; never infer typing.")
  if (hla_map) for (i in seq_along(patients)) {
    values <- hla[[patients[[i]]]]
    valid <- is.character(values) && length(values) > 0L && !anyNA(values) && all(nzchar(values))
    pp <- if (valid) peptide_presenters(list(values)) else list()
    resolved <- length(pp) > 0L && all(vapply(pp, function(p) p$status == "RESOLVED", TRUE))
    add(paste0("hla[[patient-slot-", i, "]]"), if (resolved) "PASS" else "FAIL",
      "resolved Class I allele or explicit Class II alpha~beta pair",
      if (resolved) paste(length(pp), "resolved presenters") else "UNRESOLVED_OR_EMPTY",
      if (resolved) "Syntax/chain resolution only, not proof of allele existence or model support." else
        "Example HLA-A*02:01 or HLA-DRA*01:01~HLA-DRB1*01:01. Inspect typing privately; do not guess a partner.")
  }
  valid_inputs <- is.list(inputs) && is.null(names(inputs)) && length(inputs) > 0L
  add("inputs", if (valid_inputs) "PASS" else "FAIL", "explicit input record array", if (valid_inputs) paste(length(inputs), "input records") else "INVALID_OR_EMPTY")
  if (valid_inputs) for (i in seq_along(inputs)) {
    x <- inputs[[i]];slot <- paste0("inputs[[", i, "]]")
    valid <- is.list(x) && !is.null(x$sample_id) && is.character(x$sample_id) && length(x$sample_id) == 1L &&
      !is.na(x$sample_id) && x$sample_id %in% vapply(ss, `[[`, "", "sample_id")
    add(paste0(slot, "$sample_id"), if (valid) "PASS" else "FAIL", "existing samples$sample_id", if (valid) "BINDING_FOUND" else "BINDING_MISSING",
      if (valid) "" else "Inspect this input's declared sample_id and samples privately.")
    path <- tryCatch(resolve_path(x$path, config$paths$input_root), error = function(e) NULL)
    readable <- !is.null(path) && file.exists(path) && !dir.exists(path) && file.access(path, 4L) == 0L
    add(paste0(slot, "$path"), if (readable) "PASS" else "FAIL", "existing readable regular file", if (readable) "READABLE" else "MISSING_OR_UNREADABLE",
      if (readable) "" else paste0("Inspect ", slot, "$path locally with file.exists() and file.access(...,4); paths omitted."))
    if (readable && valid && identical(x$input_role, "somatic_vcf")) {
      header <- tryCatch(suppressWarnings(first_run_vcf_header(path, x$sample_id)), error = function(e) "HEADER_UNREADABLE")
      ok <- header %in% c("SAMPLE_MATCH", "SITES_ONLY_EXPLICIT_BINDING")
      add(paste0(slot, ":VCF header"), if (!ok) "FAIL" else if (header == "SAMPLE_MATCH") "PASS" else "NOTE",
        "one exact matching sample column, or supported sites-only VCF", header,
        if (!ok) "Inspect #CHROM and the declared input binding locally; do not rename samples silently." else if (header == "SITES_ONLY_EXPLICIT_BINDING")
          "No sample genotype exists in this VCF; only the explicit input binding is checked. No GT/read counts inferred." else
            "Header only: variant bodies, REF alleles and full assembly/transcript compatibility are checked by the normal run.")
    } else if (readable && !identical(x$input_role, "somatic_vcf")) add(paste0(slot, "$input_role"), "FAIL",
      "somatic_vcf for this introductory preflight", "OUTSIDE_INTRODUCTORY_SCOPE", "Use the dedicated guide for other roles; the package's broader input contracts are unchanged.")
  }
  policy <- tryCatch(peptide_policy(config$peptides), error = function(e) NULL)
  if (is.null(policy)) add("config$peptides", "FAIL", "validated positive length vectors", "INVALID_LENGTH_POLICY", "Rebuild with neo_config(peptides=list(class_i_lengths=8:12,class_ii_lengths=15L)).") else {
    for (k in names(policy)) add(paste0("config$peptides$", k), "PASS", "effective normalized lengths",
      paste(unlist(policy[[k]]), collapse = ","), "Display only; the configuration and historical defaults are unchanged.")
  }
  requests <- config$predictors
  if (!length(requests)) add("config$predictors", "NOTE", "explicit predictor request if scoring is intended", "NOT_REQUESTED",
    "Candidate generation may run without a predictor; NOT_RUN is not zero, negative or grade E.")
  for (i in seq_along(requests)) {
    p <- requests[[i]];spec <- predictor_registry()[[p$provider]]
    path <- if (is.null(spec)) NULL else setup_executable(p$executable, spec$executable)
    add(paste0("config$predictors[[", i, "]]"), if (is.null(path)) "FAIL" else "NOTE",
      "configured executable plus compatible authorized model/runtime",
      if (is.null(path)) "EXECUTABLE_MISSING" else "EXECUTABLE_PRESENT_NOT_PROBED",
      "No provider is started. After reviewing the trusted tool, explicitly call neo_check_predictors(config=config,smoke_test=FALSE) for help/version/runtime probes; otherwise choose candidate-only config deliberately.")
    if (!is.null(spec) && isTRUE(spec$requires_models)) {
      model <- p$model_path
      exists <- is.character(model) && length(model) == 1L && !is.na(model) && dir.exists(model) && length(list.files(model, all.files = TRUE, no.. = TRUE)) > 0L
      add(paste0("config$predictors[[", i, "]]$model_path"), if (exists) "NOTE" else "FAIL",
        "nonempty explicitly supplied model directory", if (exists) "PRESENT_NOT_VALIDATED" else "MODEL_MISSING",
        "No weights are loaded or downloaded. Existence alone does not establish model usability.")
    }
  }
  out <- tryCatch(resolve_path(output_dir), error = function(e) NULL)
  parent <- out
  if (!is.null(parent)) while (!file.exists(parent) && !identical(parent, dirname(parent))) parent <- dirname(parent)
  writable <- !is.null(out) && (!file.exists(out) || dir.exists(out)) && !is.null(parent) && dir.exists(parent) && file.access(parent, 2L) == 0L
  add("output_dir", if (writable) "PASS" else "FAIL", "directory with writable existing ancestor", if (writable) "PERMISSION_CHECK_PASSED_NO_WRITE" else "NOT_WRITABLE_OR_INVALID",
    "No directory or probe file was created. Permission bits do not guarantee free space/quota or future write success; retain originals and use a new run.")
  frame <- do.call(rbind, checks)
  assembly <- config$scientific$genome_assembly
  if (is.null(assembly) || !grepl("^GRCh[0-9]+([.]p[0-9]+)?$", assembly)) assembly <- "UNRECOGNIZED_OR_UNDECLARED"
  remote <- utils::packageDescription("Neoantimon2")$RemoteSha
  if (is.null(remote) || !grepl("^[a-f0-9]{40}$", remote)) remote <- NA_character_
  list(ok = !any(frame$status == "FAIL"), checks = frame,
    software = list(R = R.version.string, package_version = as.character(utils::packageVersion("Neoantimon2")),
      exported_api_count = length(getNamespaceExports("Neoantimon2")), installed_source_commit = remote),
    effective = list(assembly = assembly, peptides = policy, lengths_source = if (is.null(config$peptides)) "package_default" else "validated_config_policy_including_filled_defaults"),
    scope = "Read-only introductory somatic checks; no inference/download/project creation; not full scientific/reference/provider validation. Declared resource checksums are honored and may read large files.",
    private = if (include_private) list(output_dir = out, samples = samples, hla = hla,
      inputs = inputs, resource_report = resources) else NULL)
}
