# Header-driven AIRR import. Counts are independent nullable observations.
tcr_chain <- function(row, sample_id, artifact, airr_version = NULL) {
  get <- function(k) recognition_nullable(row[[k]])
  scalar_text(get("sequence_id"), "AIRR sequence_id")
  locus <- get("locus")
  if (is.null(locus) || !locus %in% c("TRA", "TRB", "TRG", "TRD")) stop("Unsupported TCR locus")
  chain <- list(sample_id = sample_id, repertoire_id = get("repertoire_id"), cell_id = get("cell_id"),
    sequence_id = get("sequence_id"), locus = locus, cdr3_aa = get("junction_aa"),
    junction_nt = get("junction"), v_call = get("v_call"), d_call = get("d_call"),
    j_call = get("j_call"), c_call = get("c_call"), variable_domain_aa = get("variable_domain_aa"),
    productive = recognition_boolean(get("productive")), stop_codon = recognition_boolean(get("stop_codon")),
    vj_in_frame = recognition_boolean(get("vj_in_frame")),
    duplicate_count = recognition_number(get("duplicate_count"), TRUE),
    umi_count = recognition_number(get("umi_count"), TRUE),
    consensus_count = recognition_number(get("consensus_count"), TRUE),
    t_cell_type = get("t_cell_type"), source = get("source"), provider_version = artifact$version,
    airr_schema_version = airr_version, pairing_provenance = get("pairing_status"),
    sequence_status = recognition_sequence_status(get("junction_aa")), raw_metadata = row,
    artifact = artifact)
  chain$functional_status <- tcr_functional_status(chain)
  chain$chain_observation_id <- recognition_id("observation", chain)
  chain
}
tcr_functional_status <- function(chain) {
  negative <- identical(chain$productive, FALSE) || identical(chain$stop_codon, TRUE) ||
    identical(chain$vj_in_frame, FALSE)
  if (isTRUE(chain$productive) && negative) return("functionality_conflicting")
  if (negative) return("explicitly_nonfunctional")
  if (isTRUE(chain$productive)) return("functional")
  "functionality_unknown"
}
tcr_hypothesis_functionality <- function(chains) {
  states <- vapply(chains, `[[`, "", "functional_status")
  if ("explicitly_nonfunctional" %in% states) return("explicitly_nonfunctional")
  if ("functionality_conflicting" %in% states) return("functionality_conflicting")
  if (length(states) && all(states == "functional")) return("functional")
  "functionality_unknown"
}
tcr_chain_features <- function(chain) {
  prefix <- switch(chain$locus, TRA = "alpha", TRB = "beta", TRG = "gamma", TRD = "delta")
  values <- chain[c("cdr3_aa", "v_call", "d_call", "j_call", "c_call", "variable_domain_aa")]
  names(values) <- paste0(prefix, c("_cdr3_aa", "_v_call", "_d_call", "_j_call", "_c_call", "_variable_aa"))
  Filter(Negate(is.null), values)
}
tcr_read_airr <- function(input) {
  role <- input$input_role
  if (role == "tcr_rearrangement") return(list(Rearrangement = source_table(input$resolved_path,
    c("sequence_id", "locus", "junction_aa"), "AIRR Rearrangement"), airr_version = NULL))
  if (role == "tcr_reactivity" && grepl("\\.tsv$", input$resolved_path))
    return(list(Reactivity = source_table(input$resolved_path, c("reactivity_id", "reactivity_method",
      "reactivity_readout", "reactivity_value", "reactivity_unit"), "AIRR Reactivity"), airr_version = NULL))
  x <- if (grepl("\\.ya?ml$", input$resolved_path)) yaml::read_yaml(input$resolved_path, eval.expr = FALSE) else
    jsonlite::read_json(input$resolved_path)
  reject_duplicate_keys(x)
  x$airr_version <- recognition_nullable(x$Info$version)
  x
}
tcr_import <- function(project, max_hypotheses = 1000) {
  chains <- list(); cells <- list(); receptors <- list(); repertoires <- list(); reactivities <- list()
  for (input in Filter(function(x) x$input_role %in% c("tcr_rearrangement", "tcr_airr", "tcr_repertoire",
      "tcr_cell", "tcr_receptor", "tcr_reactivity"), recognition_inputs(project))) {
    x <- tcr_read_airr(input); artifact <- recognition_artifact(input)
    for (r in x$Rearrangement) chains[[length(chains) + 1L]] <- tcr_chain(r, input$sample_id, artifact, x$airr_version)
    wrap <- function(rows) lapply(rows, function(r) list(sample_id = input$sample_id,
      airr_schema_version = x$airr_version, artifact = artifact, raw = r))
    cells <- c(cells, wrap(x$Cell)); receptors <- c(receptors, wrap(x$Receptor))
    repertoires <- c(repertoires, wrap(x$Repertoire)); reactivities <- c(reactivities, wrap(x$Reactivity))
  }
  keys <- vapply(chains, function(c) paste(c$sample_id, c$repertoire_id, c$sequence_id, sep = "\037"), "")
  if (anyDuplicated(keys)) stop("Duplicate AIRR chain observation reference")
  cell_keys <- vapply(cells, function(c) paste(c$sample_id, c$raw$repertoire_id, c$raw$cell_id, sep = "\037"), "")
  if (anyDuplicated(cell_keys)) stop("Duplicate AIRR cell reference")
  receptor_keys <- vapply(receptors, function(r) paste(r$sample_id, r$raw$receptor_id, sep = "\037"), "")
  if (anyDuplicated(receptor_keys)) stop("Duplicate AIRR Receptor reference")
  profiles <- list(); observations <- list(); hypotheses <- list(); pairing_statuses <- list()
  add <- function(selected, state, cell = NULL, extra = list(), official_hash = NULL) {
    features <- c(unlist(lapply(selected, tcr_chain_features), recursive = FALSE), extra)
    if (anyDuplicated(names(features))) stop("Ambiguous profile feature assignment")
    phenotype <- if (!is.null(cell)) recognition_nullable(cell$raw$t_cell_type) else
      if (length(selected)) recognition_nullable(selected[[1L]]$t_cell_type) else NULL
    if (!is.null(phenotype)) features$t_cell_type <- phenotype
    profile <- recognition_profile(features)
    profiles[[length(profiles) + 1L]] <<- profile
    sample <- if (length(selected)) selected[[1L]]$sample_id else cell$sample_id
    obs <- list(receptor_feature_profile_id = profile$receptor_feature_profile_id, sample_id = sample,
      repertoire_id = if (!is.null(cell)) cell$raw$repertoire_id else selected[[1L]]$repertoire_id,
      cell_id = if (!is.null(cell)) cell$raw$cell_id else selected[[1L]]$cell_id,
      chain_observation_ids = as.list(sort(vapply(selected, `[[`, "", "chain_observation_id"))),
      pairing_status = state, pairing_provenance = cell, receptor_hash = official_hash,
      functional_status = tcr_hypothesis_functionality(selected),
      chain_functionality = lapply(selected, function(c) c[c("chain_observation_id", "productive",
        "stop_codon", "vj_in_frame", "functional_status")]),
      t_cell_type = if (!is.null(cell)) cell$raw$t_cell_type else selected[[1L]]$t_cell_type,
      support = lapply(selected, function(c) c[c("sequence_id", "duplicate_count", "umi_count", "consensus_count")]))
    obs$repertoire_observation_id <- recognition_id("observation", obs)
    observations[[length(observations) + 1L]] <<- obs
    if (state != "unpaired_chain") hypotheses[[length(hypotheses) + 1L]] <<- obs
  }
  for (c in chains) add(list(c), "unpaired_chain")
  pair_count <- sum(vapply(cells, function(cell) {
    r <- cell$raw
    explicit <- !is.null(recognition_boolean(r$virtual_pairing)) ||
      any(r$pairing_status %in% c("observed_pair", "virtual_or_inferred_pair"))
    if (!explicit) return(0L)
    ids <- unlist(r$rearrangements)
    if (is.list(r$rearrangements) && length(r$rearrangements) && is.list(r$rearrangements[[1L]]))
      ids <- vapply(r$rearrangements, `[[`, "", "sequence_id")
    selected <- Filter(function(c) c$sample_id == cell$sample_id &&
      identical(c$repertoire_id, r$repertoire_id) && c$sequence_id %in% ids, chains)
    loci <- vapply(selected, `[[`, "", "locus")
    as.numeric(sum(loci == "TRA")) * sum(loci == "TRB")
  }, 0))
  for (cell in cells) {
    r <- cell$raw; scalar_text(r$cell_id, "AIRR cell_id")
    ids <- unlist(r$rearrangements)
    if (is.list(r$rearrangements) && length(r$rearrangements) && is.list(r$rearrangements[[1L]]))
      ids <- vapply(r$rearrangements, `[[`, "", "sequence_id")
    selected <- Filter(function(c) c$sample_id == cell$sample_id &&
      identical(c$repertoire_id, r$repertoire_id) && c$sequence_id %in% ids, chains)
    if (length(setdiff(ids, vapply(selected, `[[`, "", "sequence_id"))) ||
        any(vapply(selected, function(c) !is.null(c$cell_id) && c$cell_id != r$cell_id, TRUE)))
      stop("Unresolved/conflicting Cell rearrangement reference")
    virtual <- recognition_boolean(r$virtual_pairing)
    explicit <- !is.null(virtual) || any(r$pairing_status %in% c("observed_pair", "virtual_or_inferred_pair"))
    a <- Filter(function(c) c$locus == "TRA", selected); b <- Filter(function(c) c$locus == "TRB", selected)
    # Only explicit within-cell membership authorizes receptor hypotheses.
    count <- length(a) * length(b)
    state <- if (!count || !explicit) "unpaired_chain" else if (count > 1L) "ambiguous_pairing" else
      if (isTRUE(virtual) || identical(r$pairing_status, "virtual_or_inferred_pair")) "virtual_or_inferred_pair" else "observed_pair"
    if (pair_count > max_hypotheses && explicit && count > 0L) {
      state <- "not_run/query_limit"
    } else if (explicit && count) for (alpha in a) for (beta in b) add(list(alpha, beta), state, cell)
    pairing_statuses[[length(pairing_statuses) + 1L]] <- list(sample_id = cell$sample_id,
      repertoire_id = r$repertoire_id, cell_id = r$cell_id, status = state,
      possible_hypotheses = count, provenance = cell)
  }
  for (receptor in receptors) {
    r <- receptor$raw
    scalar_text(r$receptor_id, "AIRR receptor_id")
    loci <- c(r$receptor_variable_domain_1_locus, r$receptor_variable_domain_2_locus)
    domains <- c(r$receptor_variable_domain_1_aa, r$receptor_variable_domain_2_aa)
    if (length(loci) != 2L || length(domains) != 2L || !setequal(loci, c("TRA", "TRB"))) next
    features <- stats::setNames(as.list(domains), ifelse(loci == "TRA", "alpha_variable_aa", "beta_variable_aa"))
    p <- recognition_profile(features)
    profiles[[length(profiles) + 1L]] <- p
    linked <- Filter(function(c) c$sample_id == receptor$sample_id && r$receptor_id %in% unlist(c$raw$receptors), cells)
    # AIRR hash is preserved as supplied; never constructed from partial profiles.
    for (cell in linked) {
      context <- Filter(function(o) identical(o$sample_id, cell$sample_id) &&
        identical(o$repertoire_id, cell$raw$repertoire_id) && identical(o$cell_id, cell$raw$cell_id), hypotheses)
      obs <- list(receptor_feature_profile_id = p$receptor_feature_profile_id,
        sample_id = cell$sample_id, repertoire_id = cell$raw$repertoire_id, cell_id = cell$raw$cell_id,
        chain_observation_ids = list(), pairing_status = if (isTRUE(recognition_boolean(cell$raw$virtual_pairing)))
          "virtual_or_inferred_pair" else "observed_pair", pairing_provenance = cell,
        receptor_hash = recognition_nullable(r$receptor_hash), receptor_record = receptor,
        t_cell_type = cell$raw$t_cell_type, support = list(),
        functional_status = if (length(context) == 1L) context[[1L]]$functional_status else "functionality_unknown",
        functionality_context = lapply(context, function(o) o[c("repertoire_observation_id", "functional_status", "chain_functionality")]))
      obs$repertoire_observation_id <- recognition_id("observation", obs)
      observations[[length(observations) + 1L]] <- obs
    }
  }
  list(chains = recognition_sort(chains, "chain_observation_id"),
    profiles = recognition_unique(profiles, "receptor_feature_profile_id"),
    observations = recognition_unique(observations, "repertoire_observation_id"),
    pairing_hypotheses = recognition_unique(hypotheses, "repertoire_observation_id"),
    pairing_statuses = sorted_records(pairing_statuses), cells = sorted_records(cells),
    receptors = sorted_records(receptors), repertoires = sorted_records(repertoires),
    reactivities = sorted_records(reactivities))
}
