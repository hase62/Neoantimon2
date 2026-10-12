# Registry is separate from both frozen pMHC and TCR predictor registries.
processing_registry <- function() list(
  mhcflurry_processing = list(version = "2", task = "antigen_processing", origin_dependent = TRUE,
    hla_conditioned = FALSE, lengths = as.list(8:15), model = "class1_processing"),
  netchop = list(version = "3.1", task = "proteasomal_cleavage", origin_dependent = TRUE,
    hla_conditioned = FALSE, lengths = as.list(1:100), model = "C-term-3.0"),
  netctlpan = list(version = "1.1", task = "integrated_class_i", origin_dependent = TRUE,
    hla_conditioned = TRUE, lengths = as.list(8:11), model = "netctlpan-1.1"),
  netmhcstabpan = list(version = "1.0", task = "pmhc_stability", origin_dependent = FALSE,
    hla_conditioned = TRUE, lengths = as.list(8:14), model = "netmhcstabpan-1.0"),
  fixture = list(version = "synthetic-1", task = "antigen_processing", origin_dependent = TRUE,
    hla_conditioned = FALSE, lengths = as.list(1:100), model = "deterministic-context-v1"))
processing_config <- function(x) {
  strict_keys(x, c("provider", "provider_version", "model", "options", "executable", "model_path",
    "output_input_id", "fixture", "protein_aliases"), "processing provider")
  reg <- processing_registry()[[x$provider]]
  if (is.null(reg)) stop("Unknown processing provider")
  if (is.null(x$provider_version)) x$provider_version <- reg$version
  if (is.null(x$model)) x$model <- reg$model
  if (is.null(x$options)) x$options <- list()
  strict_keys(x$options, c("use_n_flank", "use_c_flank", "flank_length", "threshold"), "processing options")
  for (k in c("use_n_flank", "use_c_flank")) {
    if (is.null(x$options[[k]])) x$options[[k]] <- TRUE
    if (!is.logical(x$options[[k]]) || length(x$options[[k]]) != 1L || is.na(x$options[[k]])) stop("Invalid flank switch")
  }
  if (is.null(x$options$flank_length)) x$options$flank_length <- 5L
  n <- x$options$flank_length
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 0 || n != floor(n)) stop("Invalid flank length")
  if (!is.null(x$options$threshold)) presentation_number(x$options$threshold)
  if (x$provider == "netchop" && !x$model %in% c("C-term-3.0", "20S-3.0")) stop("Unsupported NetChop model")
  if (x$provider != "netchop" && !identical(x$model, reg$model)) stop("Unsupported processing model")
  x
}
processing_capability <- function(config, base = getwd()) {
  config <- processing_config(config); reg <- processing_registry()[[config$provider]]
  local <- function(p) if (is.null(p)) NULL else resolve_path(p, base)
  exe <- local(config$executable); model <- local(config$model_path)
  hash <- function(p) if (!is.null(p) && file.exists(p) && !dir.exists(p))
    digest::digest(file = p, algo = "sha256", serialize = FALSE) else NULL
  # Detection never invokes a binary, installs models, or changes environment.
  list(config = config, task = reg$task, origin_dependent = reg$origin_dependent,
    hla_conditioned = reg$hla_conditioned, supported_classes = list("I"), lengths = reg$lengths,
    residue_alphabet = "ACDEFGHIKLMNPQRSTVWY", executable = exe,
    executable_sha256 = hash(exe), model_sha256 = if (!is.null(model) && dir.exists(model)) infra_hash(predictor_model_identity(model)) else hash(model),
    execution_status = if (config$provider == "fixture") "fixture" else if (!is.null(exe) && file.access(exe, 1) == 0)
      "local_executable_detected_not_executed" else "provider_unavailable",
    scope_note = if (config$provider %in% c("netctlpan", "netmhcstabpan"))
      "Non-9-mer predictions use provider-specific approximations." else "Class I only; processing is not binding or a generation verdict.")
}
processing_command <- function(cap, input, presenter = NULL, length = NULL) {
  if (is.null(cap$executable)) stop("No explicitly configured local executable")
  p <- cap$config$provider
  args <- switch(p, netchop = c("-m", if (cap$config$model == "C-term-3.0") "0" else "1",
      if (!is.null(cap$config$options$threshold)) c("-t", as.character(cap$config$options$threshold)), input),
    netctlpan = c("-f", input, "-a", presenter, "-l", as.character(length)),
    netmhcstabpan = c("-p", input, "-a", presenter, "-l", as.character(length)),
    stop("MHCflurry processing uses local Python Class1ProcessingPredictor API; no guessed CLI"))
  list(executable = cap$executable, args = as.list(args), execute = FALSE)
}
processing_origin_context <- function(origin, protein, flank_length = 5L, source_metadata = NULL) {
  start <- origin$peptide_start_aa; end <- origin$peptide_end_aa
  if (is.null(source_metadata)) source_metadata <- list(protein_origin_id = origin$sequence_context_id,
    n_terminus_known = TRUE, c_terminus_known = TRUE, translation_boundary = "explicit_complete_test_protein")
  base <- list(origin_id = origin$origin_id, protein_origin_id = source_metadata$protein_origin_id, source_metadata = source_metadata,
    source_protein_id = origin$source_protein_id, peptide_id = origin$peptide_key,
    peptide_start_aa = start, peptide_end_aa = end, source_sequence_ref = origin$mutant_sequence_ref,
    protein_sequence = protein, protein_checksum = if (is.null(protein)) NULL else peptide_digest(protein))
  if (is.null(protein)) return(c(base, list(status = "missing_context", n_flank = NULL, c_flank = NULL,
    n_terminus = NULL, c_terminus = NULL, available_n_flank = NULL, available_c_flank = NULL)))
  if (is.null(start) || is.null(end) || start < 1 || end > nchar(protein) ||
      start != floor(start) || end != floor(end) ||
      substring(protein, start, end) != origin$ligand_sequence) stop("Origin protein/coordinate mismatch")
  n <- if (start == 1L || flank_length == 0) "" else substring(protein, max(1L, start - flank_length), start - 1L)
  c <- if (end == nchar(protein) || flank_length == 0) "" else substring(protein, end + 1L, min(nchar(protein), end + flank_length))
  c(base, list(status = "available", n_flank = n, c_flank = c, n_terminus = if (start == 1L && !isTRUE(source_metadata$n_terminus_known)) NULL else start == 1L,
    c_terminus = if (end == nchar(protein) && !isTRUE(source_metadata$c_terminus_known)) NULL else end == nchar(protein), available_n_flank = nchar(n), available_c_flank = nchar(c)))
}
processing_query <- function(candidate, origin = NULL, protein = NULL, cap, source_metadata = NULL) {
  cfg <- cap$config
  context <- if (cap$origin_dependent && !is.null(origin))
    processing_origin_context(origin, protein, cfg$options$flank_length, source_metadata) else NULL
  n <- if (!is.null(context) && cfg$options$use_n_flank) context$n_flank else NULL
  c <- if (!is.null(context) && cfg$options$use_c_flank) context$c_flank else NULL
  mode <- if (length(n) && nchar(n) && length(c) && nchar(c)) "both_flanks" else
    if (length(n) && nchar(n)) "n_flank_only" else if (length(c) && nchar(c)) "c_flank_only" else "no_flanks"
  payload <- list(peptide_id = candidate$peptide_key,
    peptide = candidate$ligand_sequence, presenter = if (cap$hla_conditioned) candidate$mhc_key else NULL,
    context = context, n_flank = n, c_flank = c, context_completeness = mode,
    provider = cfg$provider, provider_version = cfg$provider_version, model = cfg$model,
    options = cfg$options, executable_sha256 = cap$executable_sha256, model_sha256 = cap$model_sha256,
    task = cap$task, scope_note = cap$scope_note)
  status <- if (!nchar(candidate$ligand_sequence) %in% unlist(cap$lengths)) "unsupported_length" else "supported"
  input <- paste0(candidate$ligand_sequence, n, c)
  if (cfg$provider %in% c("netchop", "netctlpan")) input <- if (is.null(context)) NULL else context$protein_sequence
  if (status == "supported" && cap$origin_dependent && (is.null(context) || context$status != "available")) status <- "missing_context"
  if (status == "supported" && cfg$provider %in% c("netchop", "netctlpan") &&
      (is.null(context$n_terminus) || is.null(context$c_terminus))) status <- "missing_context"
  if (status == "supported" && !grepl(paste0("^[", cap$residue_alphabet, "]+$"), input)) status <- "unsupported_residue"
  c(list(processing_query_id = presentation_id("query", processing_query_identity(payload)), fingerprint = infra_hash(payload),
    input_status = status), payload)
}
# Header-driven DTU whitespace output and MHCflurry CSV/TSV; all raw fields kept.
processing_parse <- function(path, provider) {
  if (provider == "mhcflurry_processing") {
    x <- utils::read.table(path, header = TRUE, sep = if (grepl(",", readLines(path, n = 1L))) "," else "\t",
      check.names = FALSE, colClasses = "character", quote = '"', comment.char = "", na.strings = c("", "NA"))
    if (anyDuplicated(names(x)) || !any(c("peptide", "Peptide") %in% names(x))) stop("Invalid MHCflurry processing header")
    return(lapply(seq_len(nrow(x)), function(i) lapply(as.list(x[i, , drop = FALSE]), recognition_nullable)))
  }
  lines <- readLines(path, warn = FALSE); header <- NULL; rows <- list()
  for (line in lines) {
    words <- strsplit(trimws(sub("^\\s*#\\s*", "", line)), "[[:space:]]+")[[1L]]
    if ((provider == "netchop" && all(c("pos", "AA", "score", "Ident") %in% words)) ||
        (provider == "netctlpan" && all(c("N", "Peptide", "TAP", "Comb") %in% words)) ||
        (provider == "netmhcstabpan" && all(c("pos", "peptide", "Thalf(h)") %in% words))) {
      # NetCTLpan's native 'Sequence Name' is a single identifier column.
      if (provider == "netctlpan") words <- strsplit(gsub("Sequence Name", "protein", paste(words, collapse = " "), fixed = TRUE), " ")[[1L]]
      header <- words; next
    }
    if (is.null(header) || !length(words) || !grepl("^[0-9]+$", words[[1L]])) next
    required <- setdiff(header, "BindLevel")
    if (length(words) < length(required)) stop("Truncated processing output")
    rows[[length(rows) + 1L]] <- stats::setNames(as.list(words[seq_along(required)]), required)
  }
  if (is.null(header)) stop("Unsupported processing output header")
  rows
}
processing_metric_map <- function(provider) switch(provider,
  mhcflurry_processing = list(processing_score = c("processing_score", "score"), processing_percentile = c("processing_percentile", "processing_percentile_rank")),
  netchop = list(cleavage_score = "score"),
  netctlpan = list(mhc_binding_score = "MHC", tap_transport_score = "TAP", c_terminal_cleavage_score = "Cle",
    provider_integrated_score = "Comb", integrated_percent_rank = "%Rank"),
  netmhcstabpan = list(stability_score = "Pred", half_life = "Thalf(h)", stability_percent_rank = "%Rank_Stab",
    provider_binding_score = "1-log50k", provider_affinity = "Aff(nM)", provider_binding_rank = "%Rank_aff",
    provider_integrated_score = "Combined", provider_integrated_rank = "Combined_%rank"),
  fixture = list(processing_score = "processing_score"))
processing_row_matches <- function(row, q, cfg) {
  peptide <- if (!is.null(row[["peptide"]])) row[["peptide"]] else row[["Peptide"]]
  if (!is.null(row[["processing_query_id"]])) {
    if (!identical(row[["processing_query_id"]], q$processing_query_id) ||
        (!is.null(peptide) && !identical(peptide, q$peptide))) return(FALSE)
    for (f in c("n_flank", "c_flank")) if (f %in% names(row)) {
      normalize <- function(x) if (is.null(x)) "" else x
      if (!identical(normalize(row[[f]]), normalize(q[[f]]))) return(FALSE)
    }
    return(TRUE)
  }
  if (cfg$provider == "netchop") {
    alias <- if (is.null(cfg$protein_aliases)) NULL else cfg$protein_aliases[[row[["Ident"]]]]
    return(!is.null(alias) && identical(alias, q$context$protein_origin_id))
  }
  if (!identical(peptide, q$peptide)) return(FALSE)
  if (!is.null(q$presenter)) {
    hla <- if (!is.null(row[["HLA"]])) row[["HLA"]] else row[["Allele"]]
    if (is.null(hla) || !identical(recognition_presenter(sub("^HLA-([ABC])([0-9])", "HLA-\\1*\\2", hla))$mhc_key, q$presenter)) return(FALSE)
  }
  if (cfg$provider == "netctlpan") {
    alias <- if (is.null(cfg$protein_aliases)) NULL else cfg$protein_aliases[[row[["protein"]]]]
    return(!is.null(alias) && identical(alias, q$context$protein_origin_id) &&
      identical(as.numeric(row[["N"]]) + 1, as.numeric(q$context$peptide_start_aa)))
  }
  if (cfg$provider == "mhcflurry_processing") {
    # A peptide-only processing output cannot be assigned to biological flanks.
    if (!all(c("n_flank", "c_flank") %in% names(row))) return(FALSE)
    null_empty <- function(x) if (is.null(x)) "" else x
    return(identical(null_empty(row[["n_flank"]]), null_empty(q$n_flank)) &&
      identical(null_empty(row[["c_flank"]]), null_empty(q$c_flank)))
  }
  TRUE
}
processing_predict <- function(queries, cap, output = NULL) {
  cfg <- cap$config; evidence <- list(); statuses <- list()
  rows <- if (is.null(output)) list() else processing_parse(output$resolved_path, cfg$provider)
  for (q in recognition_unique(sorted_records(queries), "processing_query_id")) {
    status <- q$input_status
    hits <- list()
    if (status == "supported") {
      if (cfg$provider == "fixture") hits <- list(list(processing_score =
        as.character(strtoi(substr(infra_hash(list(q$peptide, q$n_flank, q$c_flank)), 1, 6), 16L) / 16777215))) else
        hits <- Filter(function(r) processing_row_matches(r, q, cfg), rows)
      if (!length(hits)) status <- if (is.null(output)) "not_run/provider_unavailable" else "not_run/missing_output"
    }
    statuses[[length(statuses) + 1L]] <- list(processing_query_id = q$processing_query_id, provider = cfg$provider, status = status)
    if (status != "supported") next
    boundaries <- if (cfg$provider == "netchop") list(c_terminal = q$context$peptide_end_aa,
      n_terminal = q$context$peptide_start_aa - 1L) else list(not_applicable = NULL)
    for (boundary in names(boundaries)) {
      selected <- hits
      if (cfg$provider == "netchop") {
        pos <- boundaries[[boundary]]
        if (pos == 0L) next # Protein N terminus has no upstream cleavage bond.
        selected <- Filter(function(r) as.numeric(r[["pos"]]) == pos, hits)
        if (length(selected) && any(vapply(selected, function(r)
          r[["AA"]] != substring(q$context$protein_sequence, pos, pos), TRUE))) stop("NetChop residue/position mismatch")
      }
      if (length(selected) > 1L) stop("Ambiguous processing output for exact query/boundary")
      for (r in selected) for (metric in names(processing_metric_map(cfg$provider))) {
        fields <- processing_metric_map(cfg$provider)[[metric]]
        found <- Filter(function(f) !is.null(r[[f]]), as.list(fields))
        if (!length(found)) next
        value <- presentation_number(r[[found[[1L]]]])
        rank <- grepl("rank|percentile", metric)
        e <- list(processing_query_id = q$processing_query_id,
          peptide_id = q$peptide_id, protein_origin_id = q$context$protein_origin_id,
          provider = cfg$provider, provider_version = cfg$provider_version, model = cfg$model, options = cfg$options,
          task = cap$task, metric = metric, raw_value = value,
          unit = if (metric == "half_life") "hours" else if (rank) "percent" else if (metric == "provider_affinity") "nM" else "provider_score",
          rank = if (rank) value else NULL, score_direction = if (rank || metric == "provider_affinity") "lower" else "higher",
          boundary = boundary, position = boundaries[[boundary]], residue = r[["AA"]],
          provider_threshold = cfg$options$threshold, provider_binary_call = r[["C"]],
          n_flank = q$n_flank, c_flank = q$c_flank, context_completeness = q$context_completeness,
          presenter = q$presenter, applicability_status = "applicable", capability_status = "supported",
          raw_artifact = if (is.null(output)) list(kind = "deterministic_fixture") else recognition_artifact(output),
          raw_fields = r, scope_note = cap$scope_note)
        e$fingerprint <- infra_hash(e)
        e$processing_evidence_id <- presentation_id("evidence", processing_evidence_identity(e))
        evidence[[length(evidence) + 1L]] <- e
      }
    }
  }
  list(evidence = recognition_unique(evidence, "processing_evidence_id"), statuses = statuses)
}

# Only actual biological/model inputs identify a computation. Candidate/sample and
# source-record provenance live on links and never enter query or evidence IDs.
processing_query_identity <- function(payload) {
  common <- payload[c("peptide", "presenter", "provider", "provider_version", "model",
    "executable_sha256", "model_sha256")]
  p <- payload$provider
  options <- payload$options
  if (p == "netmhcstabpan") return(common)
  context <- payload$context
  origin <- if (is.null(context)) NULL else context[c("protein_origin_id",
    "peptide_start_aa", "peptide_end_aa")]
  if (p %in% c("netchop", "netctlpan")) {
    common$source_context <- origin
    common$sequence_supplied <- context$protein_sequence
    common$boundary_state <- if (is.null(context)) NULL else
      context[c("status", "n_terminus", "c_terminus")]
    # Flank switches do not change the full-protein DTU inputs.
    common$options <- if (p == "netchop") options[intersect(names(options), "threshold")] else list()
  } else {
    common$source_context <- origin
    common$n_flank <- payload$n_flank
    common$c_flank <- payload$c_flank
    common$context_completeness <- payload$context_completeness
    common$boundary_state <- if (is.null(context)) NULL else
      context[c("status", "n_terminus", "c_terminus")]
    # Flank configuration has already resolved to the actual model inputs.
    common$options <- list()
  }
  common
}

# One query/candidate link keeps all origin provenance paths, including for
# origin-independent stability. Link applicability cannot promote Class II.
processing_prepare <- function(candidates, catalog, cap, max_queries = Inf) {
  queries <- list(); paths <- list()
  for (candidate in sorted_records(candidates)) {
    origins <- sorted_records(Filter(function(o) o$origin_id %in% unlist(candidate$origin_ids), catalog$origins))
    if (!length(origins)) origins <- list(NULL)
    for (origin in origins) {
      protein <- if (is.null(origin)) NULL else catalog$proteins[[origin$sequence_context_id]]
      metadata <- if (is.null(origin)) NULL else catalog$boundaries[[origin$sequence_context_id]]
      q <- processing_query(candidate, origin, protein, cap, metadata)
      queries[[length(queries) + 1L]] <- q
      provenance <- if (is.null(origin)) NULL else
        processing_origin_context(origin, protein, cap$config$options$flank_length, metadata)
      paths[[length(paths) + 1L]] <- list(processing_query_id = q$processing_query_id,
        candidate_id = candidate$candidate_id, applicability = if (candidate$hla_class == "I") "applicable" else "not_applicable",
        status = if (candidate$hla_class == "I") "applicable" else "not_applicable_class_ii", origin_context = provenance)
    }
  }
  queries <- recognition_unique(sorted_records(queries), "processing_query_id")
  links <- list()
  keys <- vapply(paths, function(l) paste(l$processing_query_id, l$candidate_id, sep = "|"), "")
  for (key in sort(unique(keys))) {
    own <- paths[keys == key]; l <- own[[1L]]
    l$origin_context <- NULL
    l$origin_contexts <- sorted_records(unique(Filter(Negate(is.null), lapply(own, `[[`, "origin_context"))))
    links[[length(links) + 1L]] <- l
  }
  requested <- sum(vapply(queries, function(q) any(vapply(links, function(l)
    l$processing_query_id == q$processing_query_id && l$applicability == "applicable", TRUE)), TRUE))
  limited <- requested > max_queries
  queries <- lapply(queries, function(q) {
    applicable <- any(vapply(links, function(l) l$processing_query_id == q$processing_query_id &&
      l$applicability == "applicable", TRUE))
    if (!applicable) q$input_status <- "not_applicable_class_ii"
    else if (limited) q$input_status <- "not_run/query_limit"
    q
  })
  list(queries = queries, links = sorted_records(links), requested_queries = requested, limited = limited)
}

# Result identity excludes artifact/source-wrapper provenance, which remains in
# the record fingerprint. One logical metric/query/boundary has one evidence ID.
processing_evidence_identity <- function(e) {
  result <- e[c("processing_query_id", "metric", "raw_value", "unit", "rank",
    "score_direction", "boundary", "position", "residue")]
  if (e$provider == "netchop") {
    result$provider_threshold <- e$provider_threshold
    result$provider_binary_call <- e$provider_binary_call
  }
  result
}
