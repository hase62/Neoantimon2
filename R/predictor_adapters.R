predictor_pair_key <- function(allele, sequence)
  paste(toupper(gsub("[:*_-]", "", allele)), sequence, sep = "\037")

predictor_prepare <- function(candidates, provider, task, capability = predictor_capability(provider)) {
  if (!is.list(candidates) || (length(candidates) && anyDuplicated(vapply(candidates, `[[`, "", "candidate_id"))))
    stop("Candidates must have unique canonical IDs")
  if (!task %in% peptide_prediction_tasks()) stop("Unknown biological prediction task")
  statuses <- list(); ready <- list()
  if (length(candidates)) candidates <- candidates[order(vapply(candidates, `[[`, "", "candidate_id"), method = "radix")]
  for (candidate in candidates) {
    peptide_sequence_validate(candidate$ligand_sequence)
    state <- predictor_candidate_status(candidate, capability, provider, task)
    if (state$status == "ready") {
      key <- predictor_pair_key(state$external_allele, candidate$ligand_sequence)
      ready[[key]] <- c(ready[[key]], list(candidate))
    } else statuses[[length(statuses) + 1L]] <- c(list(candidate_id = candidate$candidate_id,
      peptide_id = candidate$sample_peptide_id, presenter = candidate$mhc_key,
      provider = provider, task = task), state)
  }
  if (length(ready)) ready <- ready[sort(names(ready), method = "radix")]
  # One external row per distinct allele/sequence pair. Biological duplicates
  # remain distinct and are fanned out by candidate_id after parsing.
  queries <- lapply(names(ready), function(key) {
    x <- ready[[key]][[1L]]
    list(pair_key = key, allele = predictor_allele(x, provider),
      peptide = x$ligand_sequence,
      candidate_ids = as.list(vapply(ready[[key]], `[[`, "", "candidate_id")))
  })
  batches <- split(queries, vapply(queries, `[[`, "", "allele"))
  list(provider = provider, task = task, capability = capability, candidates = candidates,
    statuses = statuses, queries = queries, batches = batches, reverse_map = ready)
}

predictor_input <- function(prepared, allele) {
  queries <- prepared$batches[[allele]]
  if (is.null(queries) || !length(queries)) stop("Unknown predictor batch")
  if (prepared$provider == "mhcflurry") {
    # MHCflurry accepts an allele/peptide CSV and returns those two columns.
    return(paste(c("allele,peptide", vapply(queries, function(q)
      paste(q$allele, q$peptide, sep = ","), "")), collapse = "\n"))
  }
  paste(vapply(queries, `[[`, "", "peptide"), collapse = "\n")
}

predictor_command <- function(prepared, allele, input_path, output_path = NULL) {
  provider <- prepared$provider
  if (provider == "fixture") return(list(executable = NULL, args = character(),
    input_path = input_path, output_path = output_path, allele = allele))
  path <- prepared$capability$executable_path
  if (is.null(path)) stop("Predictor executable unavailable")
  args <- if (provider == "netmhcpan") c("-p", shQuote(input_path), "-a", shQuote(allele), "-BA") else
    if (provider == "netmhciipan") c("-f", shQuote(input_path), "-inptype", "1", "-a", shQuote(allele), "-BA") else
      c(shQuote(input_path), "--out", shQuote(output_path))
  list(executable = path, args = args, input_path = input_path,
    output_path = output_path, allele = allele)
}

predictor_parse_text <- function(text, provider) {
  if (provider == "mhcflurry") {
    parsed <- tryCatch(utils::read.csv(text = text, check.names = FALSE,
      stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(parsed) || !all(c("allele", "peptide") %in% names(parsed)))
      return(list(rows = list(), problems = list("malformed_header")))
    return(list(rows = lapply(seq_len(nrow(parsed)), function(i) as.list(parsed[i,,drop=FALSE])),
      problems = list()))
  }
  lines <- strsplit(text, "\n", fixed = TRUE)[[1L]]
  header <- NULL; rows <- list(); problems <- list()
  for (line in lines) {
    line <- trimws(line)
    if (!nzchar(line) || startsWith(line, "#") || grepl("^-{4,}$", line)) next
    words <- strsplit(line, "[[:space:]]+")[[1L]]
    if (length(words) >= 4L && words[[1L]] == "Pos" && words[[2L]] %in% c("HLA", "MHC") &&
        words[[3L]] == "Peptide") { header <- words; next }
    if (!grepl("^[0-9]+$", words[[1L]])) next
    if (!is.null(header) && utils::tail(header,1L)=="BindLevel" &&
        length(words)==length(header)+1L &&
        words[[length(words)-1L]]=="<=" && words[[length(words)]] %in% c("SB","WB"))
      words <- c(utils::head(words,-2L),paste(utils::tail(words,2L),collapse=" "))
    if (is.null(header) || length(words) < length(header) - 1L ||
        length(words) > length(header) || length(words) < 4L) {
      problems[[length(problems)+1L]] <- "malformed_row"; next
    }
    if (length(words) == length(header)-1L && utils::tail(header,1L) == "BindLevel")
      words <- c(words, "NA")
    if (length(words) != length(header)) {
      problems[[length(problems)+1L]] <- "malformed_row"; next
    }
    rows[[length(rows)+1L]] <- stats::setNames(as.list(words),header)
  }
  if (is.null(header)) problems[[length(problems)+1L]] <- "malformed_header"
  list(rows = rows, problems = problems)
}

predictor_metric_specs <- function(provider, hla_class = NULL) {
  if (provider == "mhcflurry") return(list(
    mhcflurry_affinity = list(metric = "BA_nM", task = "mhc1_binding", direction = "lower_is_better", unit = "nM"),
    mhcflurry_affinity_percentile = list(metric = "BA_rank", task = "mhc1_binding", direction = "lower_is_better", unit = "percent"),
    mhcflurry_presentation_score = list(metric = "presentation_score", task = "mhc1_presentation", direction = "higher_is_better", unit = "unitless"),
    mhcflurry_presentation_percentile = list(metric = "presentation_percentile", task = "mhc1_presentation", direction = "lower_is_better", unit = "percent")))
  class_i <- provider=="netmhcpan" || identical(hla_class,"I")
  binding <- if (class_i) "mhc1_binding" else "mhc2_binding"
  presentation <- if (class_i) "mhc1_presentation" else "mhc2_presentation"
  list(`Affinity(nM)` = list(metric="BA_nM", task=binding, direction="lower_is_better", unit="nM"),
    `Aff(nM)` = list(metric="BA_nM", task=binding, direction="lower_is_better", unit="nM"),
    `%Rank_BA` = list(metric="BA_rank", task=binding, direction="lower_is_better", unit="percent"),
    Score_BA = list(metric="BA_log_score", task=binding, direction="higher_is_better", unit="unitless"),
    Score_EL = list(metric="EL_score", task=presentation, direction="higher_is_better", unit="unitless"),
    `%Rank_EL` = list(metric="EL_rank", task=presentation, direction="lower_is_better", unit="percent"),
    Rnk_EL = list(metric="EL_rank", task="mhc1_presentation", direction="lower_is_better", unit="percent"))
}

predictor_core <- function(row, candidate, provider, predictor_version, model_version, run_id) {
  if (!provider %in% c("netmhciipan","fixture") || is.null(row$Core) || is.null(row$Of) ||
      identical(row$Core,"NA")) return(NULL)
  core <- row$Core; off <- suppressWarnings(as.integer(row$Of))
  inverted <- identical(row$Inverted, "1")
  orientation <- if (inverted) "reverse" else if (is.null(row$Inverted) ||
    identical(row$Inverted,"0")) "forward" else "unknown"
  ligand <- candidate$ligand_sequence
  core_length <- nchar(core); ligand_length <- nchar(ligand)
  start <- if (is.na(off) || off < 0L || orientation=="unknown" ||
    !grepl("^[ACDEFGHIKLMNPQRSTVWYU]+$",core)) NULL else
    if (inverted) ligand_length-off-core_length+1L else off+1L
  mapped <- !is.null(start) && start >= 1L && start+core_length-1L <= ligand_length &&
    identical(substring(ligand,start,start+core_length-1L),
      if (inverted) paste(rev(strsplit(core,"",fixed=TRUE)[[1L]]),collapse="") else core)
  if (!mapped) start <- NULL
  key <- if (mapped && !inverted) biological_id("epitope_key", list(species = candidate$species,
    hla_class = "II", binding_core = core, mhc_key = candidate$mhc_key), peptide_contract()) else NULL
  list(candidate_id = candidate$candidate_id, presenter = candidate$mhc_key,
    binding_core = core, binding_core_start = start,
    binding_orientation = orientation, raw_inverted = row$Inverted,
    status = if (mapped) "evaluated" else "unresolved_mapping",
    core_mapping_status = if (mapped) "mapped" else "unresolved",
    epitope_key_status = if (!mapped) "unavailable" else
      if (inverted) "deferred_orientation" else "derived", epitope_key = key,
    predictor = provider, predictor_version = predictor_version,
    model_version = model_version, run_id = run_id, raw_offset = row$Of)
}

predictor_normalize <- function(prepared, raw_by_allele, run_id, raw_refs, model_version = NULL) {
  provider <- prepared$provider; cap <- prepared$capability
  evidence <- list(); cores <- list(); statuses <- prepared$statuses; seen <- character()
  for (allele in names(prepared$batches)) {
    raw <- raw_by_allele[[allele]]
    if (is.null(raw)) {
      for (q in prepared$batches[[allele]]) for (c in prepared$reverse_map[[q$pair_key]])
        statuses[[length(statuses)+1L]] <- list(candidate_id=c$candidate_id, status="failed",
          reason="execution_failure", provider=provider, task=prepared$task)
      next
    }
    parsed <- predictor_parse_text(raw, provider)
    if (length(parsed$problems)) for (problem in parsed$problems)
      statuses[[length(statuses)+1L]] <- list(candidate_id=NULL, status="failed",
        reason=problem, provider=provider, task=prepared$task, batch_allele=allele)
    for (row_number in seq_along(parsed$rows)) {
      row <- parsed$rows[[row_number]]
      a <- if (provider == "mhcflurry") row$allele else if (provider == "netmhcpan")
        if (is.null(row$HLA)) row$MHC else row$HLA else row$MHC
      sequence <- if (provider == "mhcflurry") row$peptide else row$Peptide
      pair <- predictor_pair_key(a, sequence)
      if (is.null(prepared$reverse_map[[pair]]) ||
          !identical(predictor_pair_key(a,sequence),predictor_pair_key(allele,sequence)) ||
          pair %in% seen) {
        statuses[[length(statuses)+1L]] <- list(candidate_id=NULL, status="failed",
          reason=if (pair %in% seen) "duplicate_output_row" else "unmapped_output_row",
          provider=provider, task=prepared$task, batch_allele=allele)
        next
      }
      seen <- c(seen,pair)
      for (candidate in prepared$reverse_map[[pair]]) {
        metrics <- predictor_metric_specs(provider,candidate$hla_class)
        emitted <- FALSE
        for (column in intersect(names(row),names(metrics))) {
          spec <- metrics[[column]]
          # A single run may emit BA and EL; both are retained as separate
          # typed evidence, regardless of the requested primary task.
          score <- suppressWarnings(as.numeric(row[[column]]))
          if (!length(score) || !is.finite(score)) {
            statuses[[length(statuses)+1L]] <- list(candidate_id=candidate$candidate_id,
              status="failed",reason=if (row[[column]] %in% c("NA","NaN",""))
                "metric_unavailable" else "malformed_score",provider=provider,
              task=spec$task,predictor_column=column)
            next
          }
          if (spec$unit == "percent" && (score < 0 || score > 100)) {
            statuses[[length(statuses)+1L]] <- list(candidate_id=candidate$candidate_id,
              status="failed",reason="malformed_score",provider=provider,task=spec$task)
            next
          }
          emitted <- TRUE
          evidence[[length(evidence)+1L]] <- list(run_id=run_id,
            analysis_context=candidate$analysis_context, sample_id=candidate$sample_id,
            peptide_id=candidate$sample_peptide_id,
            biological_peptide_ids=candidate$origin_ids,
            candidate_id=candidate$candidate_id,
            presenter=candidate$mhc_key, external_allele=a, task=spec$task,
            predictor=provider, predictor_version=cap$predictor_version,
            model_version=model_version,
            model_version_status=cap$model_version_status,
            metric=spec$metric, raw_score=score,
            percentile_rank=if (spec$unit == "percent") score else NULL,
            score_direction=spec$direction, unit=spec$unit, status="evaluated",
            message=NULL, raw_output_ref=raw_refs[[allele]],
            raw_row_index=as.integer(row_number), predictor_column=column,
            adapter_version=predictor_adapter_version(),
            evidence_schema_version=predictor_evidence_schema_version())
        }
        if (!emitted) statuses[[length(statuses)+1L]] <- list(candidate_id=candidate$candidate_id,
          status="failed",reason="malformed_output",provider=provider,task=prepared$task)
        core <- predictor_core(row,candidate,provider,cap$predictor_version,model_version,run_id)
        if (!is.null(core)) {
          core$raw_output_ref <- raw_refs[[allele]]
          core$raw_row_index <- as.integer(row_number)
          cores[[length(cores)+1L]] <- core
        }
      }
    }
    for (q in prepared$batches[[allele]]) if (!q$pair_key %in% seen)
      for (candidate in prepared$reverse_map[[q$pair_key]])
        statuses[[length(statuses)+1L]] <- list(candidate_id=candidate$candidate_id,
          status="failed",reason="missing_output_row",provider=provider,task=prepared$task)
  }
  list(evidence=evidence, core_evidence=cores, statuses=statuses)
}
