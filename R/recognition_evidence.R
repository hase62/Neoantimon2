recognition_evidence <- function(candidate_id, profile_id = NULL, scope, provider,
    metric, raw_score = NULL, rank = NULL, probability = NULL, direction = "not_defined",
    status = "available", conditioning = "unknown", completeness = "unknown",
    applicability = "supported", reason = NULL, query = NULL, provider_version = NULL,
    model_version = NULL, artifact = NULL, metadata = list(), source_id = NULL) {
  if (scope == "peptide_intrinsic" && !is.null(profile_id)) stop("Intrinsic evidence must not have a receptor")
  if (scope == "tcr_specific" && is.null(profile_id)) stop("TCR evidence requires a concrete profile")
  identity <- list(candidate_id = candidate_id, receptor_feature_profile_id = profile_id, scope = scope,
    provider = provider, provider_version = provider_version, model_version = model_version,
    query_id = query$recognition_query_id, metric = metric, raw_score = raw_score, rank = rank,
    probability = probability, status = status, source_id = source_id)
  c(list(recognition_evidence_id = recognition_id("evidence", identity)), identity,
    list(task = if (scope %in% c("peptide_intrinsic", "tcr_specific")) paste0("recognition.", scope) else scope,
      score_direction = direction, mhc_conditioning = conditioning, input_feature_completeness = completeness,
      applicability_status = applicability, capability_reason = reason,
      run_id = NULL, fingerprint = infra_hash(list(identity = identity, query_fingerprint = query$fingerprint,
        artifact = artifact, metadata = metadata)), raw_artifact = artifact, metadata = metadata))
}
recognition_output_table <- function(path) {
  lines <- reference_lines(path)
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  if (!length(lines)) stop("Empty recognition provider output")
  sep <- if (grepl("\t", lines[[1L]], fixed = TRUE)) "\t" else ","
  frame <- utils::read.table(text = paste(lines, collapse = "\n"), sep = sep, header = TRUE,
    quote = if (sep == ",") '"' else "", comment.char = "", check.names = FALSE,
    colClasses = "character", na.strings = NULL, fill = FALSE, stringsAsFactors = FALSE)
  if (anyDuplicated(names(frame)) || any(!nzchar(names(frame)))) stop("Duplicate recognition output column")
  lapply(seq_len(nrow(frame)), function(i) lapply(frame, `[[`, i))
}
recognition_parse <- function(path, cap) {
  rows <- recognition_output_table(path)
  lapply(rows, function(r) {
    require_columns <- function(fields) if (!all(fields %in% names(r))) stop("Unsupported recognition output header")
    feature <- list(); presenter <- NULL; peptide <- NULL; metrics <- list()
    value <- function(k) recognition_number(r[[k]])
    if (cap$provider == "prime") {
      require_columns("Peptide"); peptide <- r$Peptide
      alleles <- sub("^%Rank_", "", names(r)[grepl("^%Rank_[ABC][0-9]+$", names(r))])
      if (!length(alleles)) stop("PRIME output requires per-allele metrics; best allele cannot stand for other alleles")
      return(lapply(alleles, function(a) {
        require_columns(c(paste0("%Rank_", a), paste0("Score_", a), paste0("%RankBinding_", a)))
        list(peptide = peptide, provider_hla = a, features = list(), raw = r,
          metrics = list(list(metric = "PRIME.rank", rank = value(paste0("%Rank_", a)), direction = "lower"),
            list(metric = "PRIME.score", raw_score = value(paste0("Score_", a)), direction = "higher"),
            list(metric = "MixMHCpred.rank_in_PRIME_output", rank = value(paste0("%RankBinding_", a)), direction = "lower")))
      }))
    }
    if (cap$provider == "pmtnet") {
      require_columns(c("CDR3", "Antigen", "HLA", "Rank"))
      peptide <- r$Antigen; presenter <- r$HLA; feature <- list(beta_cdr3_aa = r$CDR3)
      metrics <- list(list(metric = "pMTnet.rank", rank = value("Rank"), direction = "lower"))
    } else if (cap$provider == "ergo_ii") {
      require_columns(c("TRB", "Peptide", "Score")); peptide <- r$Peptide; presenter <- recognition_nullable(r$MHC)
      mapping <- c(alpha_cdr3_aa = "TRA", beta_cdr3_aa = "TRB", alpha_v_call = "TRAV",
        alpha_j_call = "TRAJ", beta_v_call = "TRBV", beta_j_call = "TRBJ", t_cell_type = "T-Cell-Type")
      require_columns(unname(mapping[unlist(cap$required_features)]))
      feature <- stats::setNames(lapply(mapping, function(k) recognition_nullable(r[[k]])), names(mapping))
      metrics <- list(list(metric = "ERGO-II.score", raw_score = value("Score"), direction = "higher"))
    } else if (cap$provider == "nettcr_2_0") {
      require_columns(c("peptide", "CDR3b", "prediction")); peptide <- r$peptide
      if ("alpha_cdr3_aa" %in% unlist(cap$required_features)) require_columns("CDR3a")
      feature <- list(alpha_cdr3_aa = recognition_nullable(r$CDR3a), beta_cdr3_aa = r$CDR3b)
      presenter <- "HLA-A*02:01"
      metrics <- list(list(metric = "NetTCR-2.0.score", raw_score = value("prediction"), direction = "higher"))
    } else stop("Fixture provider uses deterministic in-memory output")
    list(list(peptide = peptide, provider_hla = presenter, features = feature, metrics = metrics, raw = r))
  }) |> unlist(recursive = FALSE)
}
recognition_predict <- function(prepared, cap, output = NULL) {
  records <- if (is.null(output)) list() else recognition_parse(output$resolved_path, cap)
  artifact <- if (is.null(output)) NULL else recognition_artifact(output)
  evidence <- list(); statuses <- list()
  for (batch in prepared$batches) for (q in batch) {
    status <- q$input_status; reason <- if (status == "supported") NULL else status
    metrics <- list(); raw <- NULL
    if (status == "supported" && cap$provider == "fixture") {
      score <- strtoi(substr(q$fingerprint, 1L, 6L), 16L) / 16777215
      metrics <- list(list(metric = "fixture.score", raw_score = score, direction = "higher"))
      raw <- list(synthetic_fixture = TRUE); status <- "available"
    } else if (status == "supported" && !is.null(output)) {
      match <- Filter(function(r) identical(r$peptide, q$peptide) &&
        (q$mhc_conditioning == "absent" || identical(r$provider_hla, q$provider_hla)) &&
        all(vapply(names(q$features_used), function(f) identical(r$features[[f]], q$features_used[[f]]), TRUE)), records)
      if (length(match) > 1L) stop("Ambiguous recognition provider output for query")
      if (!length(match)) { status <- "not_run"; reason <- "output_not_supplied_for_query" }
      else { status <- "available"; metrics <- match[[1L]]$metrics; raw <- match[[1L]]$raw }
    } else if (status == "supported") { status <- "not_run"; reason <- "local_execution_not_implemented/output_not_supplied" }
    statuses[[length(statuses) + 1L]] <- list(recognition_query_id = q$recognition_query_id,
      candidate_id = q$candidate_id, receptor_feature_profile_id = q$receptor_feature_profile_id,
      provider = cap$provider, status = status, reason = reason, receptor_functionality = q$receptor_functionality)
    for (m in metrics) evidence[[length(evidence) + 1L]] <- recognition_evidence(q$candidate_id,
      q$receptor_feature_profile_id, cap$scope, cap$provider, m$metric, raw_score = m$raw_score,
      rank = m$rank, probability = m$probability, direction = m$direction, query = q,
      conditioning = q$mhc_conditioning, completeness = q$input_completeness,
      applicability = cap$applicability, provider_version = cap$provider_version, model_version = cap$model_version,
      artifact = artifact, metadata = list(raw_fields = raw, model = cap$model, features_used = q$features_used,
        receptor_functionality = q$receptor_functionality,
        presenter = q$presenter, synthetic_fixture = cap$provider == "fixture" || cap$config$fixture,
        ranking_background = if (cap$provider == "pmtnet") "provider rank against background TCRs for the same pMHC; no probability conversion" else NULL))
  }
  list(evidence = recognition_sort(evidence, "recognition_evidence_id"), statuses = statuses)
}
# Cell measurements retain their target level even when one profile is linked.
recognition_experiments <- function(candidates, repertoire) {
  evidence <- list(); measurements <- list()
  for (wrapped in repertoire$reactivities) {
    r <- wrapped$raw
    for (f in c("reactivity_id", "reactivity_method", "reactivity_readout", "reactivity_unit")) scalar_text(r[[f]], f)
    measurement <- list(reactivity_id = r$reactivity_id, sample_id = wrapped$sample_id,
      repertoire_id = r$repertoire_id, cell_id = r$cell_id, receptor_id = r$receptor_id,
      peptide = r$peptide_sequence_aa, presenter = recognition_reactivity_presenter(r),
      assay_type = r$reactivity_method, readout = r$reactivity_readout,
      value = recognition_number(r$reactivity_value), unit = r$reactivity_unit,
      qualitative_annotation = recognition_nullable(r$qualitative_annotation),
      threshold = recognition_nullable(r$threshold), threshold_reference = recognition_nullable(r$threshold_reference),
      external_references = r$reactivity_refs, method_evidence_type = if (r$reactivity_method == "inferred")
        "computational_annotation_not_a_validated_prediction" else if (r$reactivity_method == "annotated")
        "external_annotation" else "experimental_reactivity", raw_fields = r, artifact = wrapped$artifact,
      target_level = if (!is.null(recognition_nullable(r$receptor_id))) "receptor_declared" else "cell")
    if (!is.null(measurement$threshold) && is.null(measurement$threshold_reference))
      stop("Reactivity threshold requires an explicit documented reference")
    declared <- recognition_nullable(r$receptor_id)
    related <- Filter(function(o) o$sample_id == wrapped$sample_id &&
      (if (!is.null(declared)) identical(o$receptor_record$raw$receptor_id, declared) else
        !is.null(recognition_nullable(r$cell_id)) && identical(o$repertoire_id, r$repertoire_id) && identical(o$cell_id, r$cell_id)) &&
      o$pairing_status != "unpaired_chain", repertoire$observations)
    ids <- sort(unique(vapply(related, `[[`, "", "receptor_feature_profile_id")))
    measurement$possible_receptor_profile_ids <- as.list(ids)
    measurement$receptor_functionality_context <- lapply(related, function(o)
      o[intersect(c("repertoire_observation_id", "functional_status", "chain_functionality", "functionality_context"), names(o))])
    measurement$receptor_assignment <- if (!is.null(declared) && length(ids) == 1L) "declared_receptor_profile_link" else if (length(ids) == 1L &&
      all(vapply(related, function(o) o$pairing_status == "observed_pair", TRUE)))
        "unique_cell_profile_link/causality_not_inferred" else if (length(ids)) "ambiguous_pairing" else "unresolved"
    measurement$measurement_id <- recognition_id("evidence", list(sample = wrapped$sample_id,
      repertoire = r$repertoire_id, record = r, source_version = wrapped$artifact$version))
    measurements[[length(measurements) + 1L]] <- measurement
    if (r$reactivity_method == "inferred") next
    matches <- Filter(function(c) c$sample_id == wrapped$sample_id &&
      identical(c$ligand_sequence, r$peptide_sequence_aa) &&
      (measurement$presenter$status == "unknown" || (measurement$presenter$status == "known" &&
        identical(c$mhc_key, measurement$presenter$mhc_key))), candidates)
    for (c in matches) evidence[[length(evidence) + 1L]] <- recognition_evidence(c$candidate_id,
      if (length(ids) == 1L && (startsWith(measurement$receptor_assignment, "unique") ||
        identical(measurement$receptor_assignment, "declared_receptor_profile_link"))) ids[[1L]] else NULL,
      if (r$reactivity_method %in% c("annotated", "inferred")) "external_annotation" else "experimental_reactivity",
      "airr_reactivity", r$reactivity_readout, raw_score = measurement$value,
      conditioning = if (measurement$presenter$status == "known") "present" else "unknown",
      completeness = measurement$receptor_assignment, provider_version = wrapped$artifact$version,
      artifact = wrapped$artifact, metadata = measurement, source_id = measurement$measurement_id)
  }
  list(measurements = recognition_sort(measurements, "measurement_id"), evidence = evidence)
}
recognition_vdjdb <- function(input) {
  rows <- source_table(input$resolved_path, c("gene", "cdr3", "antigen.epitope", "mhc.a", "species"), "VDJdb")
  artifact <- recognition_artifact(input)
  lapply(rows, function(r) {
    paired <- !is.null(recognition_nullable(r$cdr3.alpha)) && !is.null(recognition_nullable(r$cdr3.beta))
    features <- if (paired) list(alpha_cdr3_aa = r$cdr3.alpha, beta_cdr3_aa = r$cdr3.beta) else
      stats::setNames(list(r$cdr3), switch(r$gene, TRA = "alpha_cdr3_aa", TRB = "beta_cdr3_aa",
        TRG = "gamma_cdr3_aa", TRD = "delta_cdr3_aa", stop("Unsupported VDJdb chain")))
    if (paired) {
      for (prefix in c("alpha", "beta")) for (field in c("v_call", "j_call")) {
        value <- recognition_nullable(r[[paste0(substr(field, 1L, 1L), ".", prefix)]])
        if (!is.null(value)) features[[paste0(prefix, "_", field)]] <- value
      }
    } else for (field in c("v_call", "j_call")) {
      # VDJdb's distributed table uses v.segm / j.segm.
      value <- recognition_nullable(r[[paste0(substr(field, 1L, 1L), ".segm")]])
      prefix <- sub("_cdr3_aa$", "", names(features)[[1L]])
      if (!is.null(value)) features[[paste0(prefix, "_", field)]] <- value
    }
    raw_hla <- recognition_nullable(r$mhc.a)
    if (identical(r$mhc.class, "II") && !is.null(raw_hla) && !is.null(recognition_nullable(r$mhc.b)))
      raw_hla <- paste(raw_hla, r$mhc.b, sep = "~")
    out <- list(features = features, completeness = if (paired) "paired_cdr3_profile" else "single_chain_annotation",
      peptide = r$antigen.epitope, presenter = recognition_presenter(raw_hla), species = r$species,
      release = input$version, raw_fields = r, artifact = artifact,
      match_semantics = "exact_supplied_features_only")
    out$annotation_id <- recognition_id("evidence", list(release = input$version, record = r))
    out
  })
}
recognition_curated_matches <- function(candidates, repertoire, annotations) {
  evidence <- list(); observations <- list()
  for (a in annotations) {
    sequences <- a$features[grepl("_aa$", names(a$features))]
    if (any(vapply(sequences, recognition_sequence_status, "") != "valid")) next
    if (!a$species %in% c("HomoSapiens", "Homo sapiens", "human", "NCBITaxon:9606")) next
    profiles <- Filter(function(p) all(vapply(names(a$features), function(f)
      identical(p$features[[f]], a$features[[f]]), TRUE)), repertoire$profiles)
    for (p in profiles) {
      observed <- Filter(function(o) o$receptor_feature_profile_id == p$receptor_feature_profile_id,
        repertoire$observations)
      for (c in candidates) {
        own <- Filter(function(o) o$sample_id == c$sample_id, observed)
        if (!length(own) || !identical(c$ligand_sequence, a$peptide) ||
            !(a$presenter$status == "unknown" || (a$presenter$status == "known" && identical(c$mhc_key, a$presenter$mhc_key)))) next
        meta <- list(annotation = a, match_status = if (a$presenter$status == "known")
          "exact_external_pmhc_specificity_match" else "exact_external_sequence_match/restriction_unknown",
          patient_assay = FALSE, representation = a$completeness,
          observations = own, repertoire_status = "receptor_observed")
        evidence[[length(evidence) + 1L]] <- recognition_evidence(c$candidate_id, p$receptor_feature_profile_id,
          "external_annotation", "vdjdb", "curated_specificity_annotation", provider_version = a$release,
          conditioning = if (a$presenter$status == "known") "present" else "unknown",
          completeness = a$completeness, artifact = a$artifact, metadata = meta, source_id = a$annotation_id)
        for (o in own) observations[[length(observations) + 1L]] <- list(candidate_id = c$candidate_id,
          receptor_feature_profile_id = p$receptor_feature_profile_id, repertoire_observation_id = o$repertoire_observation_id,
          annotation_id = a$annotation_id, status = "receptor_observed", observation = o)
      }
    }
  }
  list(evidence = recognition_unique(evidence, "recognition_evidence_id"), observations = sorted_records(observations))
}
