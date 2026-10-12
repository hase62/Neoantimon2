# Read-only overlay over the Phase C graph. No native field enters Phase B.
native_allo_target_info <- function(graph, target) {
  node <- indirect_report_one(graph$nodes[[target$target_type]], 'id', target$target_id)
  d <- node$data
  out <- list(source_alleles = list(), source_genes = list(), peptide_keys = list(), presenters = list())
  source <- function(id) {
    s <- indirect_report_one(graph$nodes$source_occurrence, 'id', id)$data$semantics
    list(allele = s$allele_hypothesis, gene = sub('\\*.*', '', s$allele_hypothesis))
  }
  if (target$target_type == 'source_occurrence') {
    s <- source(target$target_id); out$source_alleles <- list(s$allele); out$source_genes <- list(s$gene)
  }
  if (target$target_type == 'peptide_occurrence') {
    s <- source(d$semantics$source_occurrence_id)
    out$source_alleles <- list(s$allele); out$source_genes <- list(s$gene)
    out$peptide_keys <- list(d$semantics$peptide_key)
  }
  if (target$target_type == 'peptide') out$peptide_keys <- list(target$target_id)
  if (target$target_type == 'pmhc') {
    out$presenters <- list(d$presenter); out$peptide_keys <- list(d$peptide_key)
  }
  if (target$target_type == 'core_assignment') {
    p <- indirect_report_one(graph$nodes$pmhc, 'id', d$semantics$indirect_pmhc_id)$data
    out$presenters <- list(p$presenter); out$peptide_keys <- list(p$peptide_key)
  }
  if (target$target_type == 'aggregate_unit') {
    out$presenters <- list(d$presenter); out$peptide_keys <- d$peptide_keys
  }
  out
}
native_allo_applicability <- function(r, target, context, scope, graph) {
  answer <- function(status, reason) list(target_type = target$target_type, target_id = target$target_id,
    status = status, reason_codes = allo_set(c(r$reason_codes, list(reason))),
    measurement_scope_only = TRUE, baseline_unchanged = TRUE)
  c <- context$semantics; d <- r$dimension; a <- r$analyte
  expected_type <- switch(d, source_expression = 'expression_assay', presenter_expression = 'expression_assay',
    hla_allele_loh = 'allele_loh_assay', antigen_processing_machinery = 'apm_assay',
    ms_observation = 'mass_spectrometry', tcr_recognition = 'tcr_recognition_assay')
  if (r$assay$evidence_type != expected_type) return(answer('not_applicable', 'assay_dimension_mismatch'))
  info <- native_allo_target_info(graph, target)
  allowed <- switch(d, source_expression = c('source_occurrence', 'peptide_occurrence'),
    presenter_expression = c('context', 'pmhc', 'core_assignment', 'aggregate_unit'),
    hla_allele_loh = c('source_occurrence', 'peptide_occurrence', 'context', 'pmhc', 'core_assignment', 'aggregate_unit'),
    antigen_processing_machinery = 'context', ms_observation = c('peptide', 'peptide_occurrence', 'pmhc', 'aggregate_unit'),
    tcr_recognition = c('context', 'pmhc', 'core_assignment', 'aggregate_unit'))
  if (!target$target_type %in% allowed) return(answer('not_applicable', 'dimension_target_scope_mismatch'))
  role <- r$subject_role
  expected_role <- switch(d, source_expression = 'source', presenter_expression = 'presenter',
    tcr_recognition = 'responding', hla_allele_loh = c('source', 'presenter'),
    antigen_processing_machinery = c('source', 'presenter'), ms_observation = c('source', 'presenter'))
  if (!role %in% expected_role || r$subject_id != c[[paste0(role, '_subject_id')]])
    return(answer('not_applicable', 'subject_role_mismatch'))
  compartment <- switch(role, source = scope$source_compartment, presenter = c$apc_compartment,
    responding = scope$responding_compartment)
  if (r$tissue != c$tissue || r$compartment != compartment)
    return(answer('not_applicable', 'tissue_compartment_mismatch'))
  if (r$occasion_id != c$occasion_id || (!is.null(r$sample_id) && !r$sample_id %in% unlist(scope[[paste0(role, '_sample_ids')]])))
    return(answer('insufficient_scope', 'sample_occasion_mismatch'))
  if (d == 'hla_allele_loh') {
    source_target <- target$target_type %in% c('source_occurrence', 'peptide_occurrence')
    if ((source_target && role != 'source') || (!source_target && role != 'presenter'))
      return(answer('not_applicable', 'loh_allele_role_mismatch'))
  }
  if (r$evaluability == 'failed') return(answer('failed', 'technical_failure_not_negative'))
  if (r$evaluability == 'not_applicable') return(answer('not_applicable', 'observation_not_applicable'))
  if (r$evaluability == 'not_evaluated') return(answer('not_evaluated', 'observation_not_evaluated'))
  if (r$evaluability == 'insufficient_resolution') return(answer('insufficient_scope', 'observation_resolution_insufficient'))
  if (r$evaluability == 'ambiguous') return(answer('ambiguous', 'observation_ambiguous'))
  if (d == 'source_expression' || (d == 'hla_allele_loh' && role == 'source')) {
    if (a$resolution == 'gene') {
      if (!any(unlist(a$genes) %in% unlist(info$source_genes))) return(answer('not_applicable', 'source_gene_mismatch'))
      return(answer('ambiguous', 'gene_measurement_not_allele_specific'))
    }
    if (a$resolution != 'allele') return(answer('insufficient_scope', 'source_allele_resolution_insufficient'))
    if (!identical(a$alleles, info$source_alleles)) return(answer('not_applicable', 'source_allele_mismatch'))
  }
  if (d == 'presenter_expression') {
    if (a$resolution == 'gene') {
      genes <- if (target$target_type == 'context') c('HLA-DRA', 'HLA-DRB1') else
        sub('\\*.*', '', unlist(strsplit(unlist(info$presenters), '~', fixed = TRUE)))
      if (!all(unlist(a$genes) %in% genes)) return(answer('not_applicable', 'presenter_gene_mismatch'))
      if (target$target_type != 'context') return(answer('insufficient_scope', 'gene_expression_not_molecular_surface_abundance'))
    } else {
      if (a$resolution != 'molecule' || is.null(a$presenter)) return(answer('insufficient_scope', 'presenter_molecule_resolution_insufficient'))
      if (!a$presenter %in% unlist(info$presenters)) return(answer('not_applicable', 'presenter_molecule_mismatch'))
    }
  }
  if (d == 'hla_allele_loh' && role == 'presenter') {
    if (a$resolution != 'allele') return(answer('insufficient_scope', 'loh_allele_resolution_insufficient'))
    presenters <- if (target$target_type == 'context') unlist(lapply(graph$nodes$pmhc, function(n) n$data$presenter)) else unlist(info$presenters)
    alleles <- unique(unlist(strsplit(presenters, '~', fixed = TRUE)))
    if (!a$alleles[[1L]] %in% alleles) return(answer('not_applicable', 'presenter_allele_mismatch'))
  }
  if (d == 'antigen_processing_machinery') {
    if (role != 'presenter') return(answer('not_applicable', 'apm_not_presenting_compartment'))
    if (is.null(a$hla_class)) return(answer('insufficient_scope', 'apm_class_scope_missing'))
    if (a$hla_class != 'II' || a$pathway_applicability == 'not_applicable')
      return(answer('not_applicable', 'class_ii_pathway_not_supported'))
    if (a$pathway_applicability == 'not_evaluated') return(answer('not_evaluated', 'class_ii_pathway_not_evaluated'))
    # Class I machinery cannot be relabeled as reviewed Class II processing.
    if (any(unlist(a$genes) %in% c('TAP1', 'TAP2', 'TAPBP', 'ERAP1', 'ERAP2', 'PSMB8', 'PSMB9')))
      return(answer('not_applicable', 'class_i_machinery_not_class_ii_rule'))
  }
  if (d == 'ms_observation') {
    if (!a$peptide_key %in% unlist(info$peptide_keys)) return(answer('not_applicable', 'observed_peptide_target_mismatch'))
    if (r$result_state == 'not_detected' && r$assignments$absence_scope == 'sample_only')
      return(answer('insufficient_scope', 'ms_nondetection_not_exhaustive_absence'))
    if (target$target_type != 'peptide') {
      if (r$assignments$source_status %in% c('unresolved', 'ambiguous') ||
          r$assignments$presenter_status %in% c('unresolved', 'ambiguous'))
        return(answer('ambiguous', 'ms_source_and_presenter_assignments_separate'))
      if (r$assignments$source_status != 'resolved' || r$assignments$presenter_status != 'resolved')
        return(answer('insufficient_scope', 'ms_assignment_scope_insufficient'))
      if (target$target_type %in% c('pmhc', 'aggregate_unit')) {
        node <- indirect_report_one(graph$nodes[[target$target_type]], 'id', target$target_id)$data
        occurrence_ids <- node$occurrence_ids
        origins <- vapply(occurrence_ids, function(id) indirect_report_one(graph$nodes$peptide_occurrence, 'id', id)$data$semantics$source_occurrence_id, '')
        if (!all(unlist(r$assignments$source_occurrence_ids) %in% origins))
          return(answer('not_applicable', 'ms_source_assignment_mismatch'))
      }
      if (length(info$presenters) && !all(unlist(info$presenters) %in% unlist(r$assignments$presenters)))
        return(answer('not_applicable', 'ms_presenter_assignment_mismatch'))
      if (target$target_type == 'peptide_occurrence') {
        origin <- indirect_report_one(graph$nodes$peptide_occurrence, 'id', target$target_id)$data$semantics$source_occurrence_id
        if (!origin %in% unlist(r$assignments$source_occurrence_ids)) return(answer('not_applicable', 'ms_source_assignment_mismatch'))
      }
    }
  }
  if (d == 'tcr_recognition') {
    if (target$target_type == 'context') return(answer('insufficient_scope', 'tcr_ligand_target_unresolved'))
    if (is.null(a$presenter) || !a$presenter %in% unlist(info$presenters) ||
        is.null(a$peptide_key) || !a$peptide_key %in% unlist(info$peptide_keys))
      return(answer('not_applicable', 'tcr_ligand_target_mismatch'))
  }
  negative <- r$result_state %in% c('negative', 'allele_loss_supported', 'not_detected')
  answer(if (negative) 'applicable_negative' else 'applicable_supported', 'scoped_observation_only_no_count_effect')
}
native_allo_reverse_links <- function(detail, records) {
  lapply(detail$included_units, function(u) {
    refs <- c(list(list(target_type = 'aggregate_unit', target_id = u$counted_unit_key)),
      lapply(u$supporting_pmhc_ids, function(id) list(target_type = 'pmhc', target_id = id)),
      lapply(u$peptide_keys, function(id) list(target_type = 'peptide', target_id = id)),
      lapply(u$occurrence_ids, function(id) list(target_type = 'peptide_occurrence', target_id = id)),
      lapply(u$source_occurrence_ids, function(id) list(target_type = 'source_occurrence', target_id = id)))
    cores <- Filter(function(n) n$data$semantics$indirect_pmhc_id %in% unlist(u$supporting_pmhc_ids) &&
      identical(n$data$epitope_key, u$counted_unit_key), detail$graph$nodes$core_assignment)
    refs <- c(refs, lapply(cores, function(n) list(target_type = 'core_assignment', target_id = n$id)))
    keys <- vapply(refs, canonical_bytes, '')
    associated <- Filter(function(r) any(vapply(r$targets, canonical_bytes, '') %in% keys), records)
    list(member_id = u$counted_unit_key, evidence_ids = allo_set(lapply(associated, `[[`, 'evidence_id')),
      semantics = 'association through explicit frozen relationships; no propagated support or penalty')
  })
}
native_allo_integrity <- function(overlay, explanation) {
  ids <- vapply(overlay$records, `[[`, '', 'evidence_id')
  if (anyDuplicated(ids)) stop('Duplicate native evidence identity')
  for (r in overlay$records) {
    if (r$evidence_id != native_allo_id(r$observation)) stop('Native evidence payload mismatch')
    detail <- indirect_report_one(explanation$contexts, 'allo_context_id', r$allo_context_id)
    for (ref in r$targets) native_allo_target_info(detail$graph, ref)
    if (!identical(r$targets, lapply(r$applicability, function(a) a[c('target_type', 'target_id')]))) stop('Native applicability target mismatch')
  }
  for (c in overlay$contexts) for (link in c$member_links) {
    detail <- indirect_report_one(explanation$contexts, 'allo_context_id', c$allo_context_id)
    if (!link$member_id %in% unlist(lapply(detail$included_units, `[[`, 'counted_unit_key')) ||
        !all(unlist(link$evidence_ids) %in% ids)) stop('Native reverse link reference mismatch')
  }
  invisible(TRUE)
}
native_allo_build <- function(project) {
  if (!native_allo_enabled(project)) return(NULL)
  if (!indirect_enabled(project)) stop('Native allo overlay requires frozen baseline')
  config <- native_allo_config(project$config$allo_native_evidence)
  baseline <- indirect_read(project); before <- serialize(baseline, NULL)
  foundation <- allo_read(project)
  explanation <- indirect_public_explanation(baseline, foundation)
  context_ids <- vapply(explanation$contexts, `[[`, '', 'allo_context_id')
  if (!setequal(vapply(config$scopes, `[[`, '', 'allo_context_id'), context_ids)) stop('Native scopes must match frozen contexts exactly')
  rows <- lapply(config$observations$manifest$observations, native_allo_row)
  for (r in rows) {
    if (!r$allo_context_id %in% context_ids) stop('Native context reference unresolved')
    detail <- indirect_report_one(explanation$contexts, 'allo_context_id', r$allo_context_id)
    for (ref in r$targets) native_allo_target_info(detail$graph, ref)
    for (id in r$assignments$source_occurrence_ids) indirect_report_one(detail$graph$nodes$source_occurrence, 'id', id)
    presenters <- unique(unlist(lapply(detail$graph$nodes$pmhc, function(n) n$data$presenter)))
    if (!all(unlist(r$assignments$presenters) %in% presenters)) stop('MS presenter assignment reference unresolved')
    # Supplied assignment semantics must resolve in the frozen graph. Candidate
    # membership alone is never used to create or resolve an assignment.
    for (e in r$assignments$assignment_evidence) {
      indirect_report_one(detail$graph$nodes$peptide, 'id', e$peptide_key)
      if (e$assignment_kind == 'source_origin') {
        indirect_report_one(detail$graph$nodes$source_occurrence, 'id', e$source_occurrence_id)
        matching <- Filter(function(n) n$data$semantics$peptide_key == e$peptide_key &&
          n$data$semantics$source_occurrence_id == e$source_occurrence_id, detail$graph$nodes$peptide_occurrence)
        if (!length(matching)) stop('Assignment evidence source/peptide frozen graph mismatch')
      }
      if (e$assignment_kind == 'presenter' && !any(vapply(detail$graph$nodes$pmhc,
          function(n) n$data$peptide_key == e$peptide_key && n$data$presenter == e$presenter, TRUE)))
        stop('Assignment evidence presenter/peptide frozen graph mismatch')
    }
  }
  ids <- vapply(rows, function(r) native_allo_id(r), '')
  records <- lapply(sort(unique(ids), method = 'radix'), function(id) {
    matches <- rows[ids == id]; r <- matches[[1L]]
    context <- indirect_report_one(foundation$contexts, 'allo_context_id', r$allo_context_id)
    scope <- indirect_report_one(config$scopes, 'allo_context_id', r$allo_context_id)
    detail <- indirect_report_one(explanation$contexts, 'allo_context_id', r$allo_context_id)
    r <- r[setdiff(names(r), c('provenance', 'raw_values', 'readonly_references'))]
    list(evidence_id = id, allo_context_id = r$allo_context_id, dimension = r$dimension,
      targets = r$targets, observation = r,
      applicability = lapply(r$targets, function(t) native_allo_applicability(r, t, context, scope, detail$graph)),
      raw_observations = allo_records(lapply(matches, function(z) z[c('raw_values', 'provenance', 'readonly_references')])))
  })
  contexts <- lapply(explanation$contexts, function(detail) {
    rr <- Filter(function(r) r$allo_context_id == detail$allo_context_id, records)
    list(allo_context_id = detail$allo_context_id, member_links = native_allo_reverse_links(detail, rr))
  })
  overlay <- list(schema_version = native_allo_version(), model_id = 'neo-allo-hla-evidence-v1',
    official_reproduction = FALSE, baseline_unchanged = TRUE,
    allo_policy_id = baseline$allo_policy_id, allo_analysis_id = baseline$allo_analysis_id,
    scopes = config$scopes, records = records, contexts = contexts,
    provenance = config$observations[c('path', 'sha256')])
  native_allo_integrity(overlay, explanation)
  if (!identical(before, serialize(indirect_read(project), NULL)) ||
      !identical(explanation, indirect_public_explanation(project$indirect_hla, foundation)))
    stop('Native overlay mutated frozen baseline or compatibility audit')
  overlay
}
native_allo_read <- function(project) {
  expected <- native_allo_build(project)
  if (!identical(expected, project$allo_native_evidence)) stop('Native overlay derived contract mismatch')
  expected
}
native_allo_summary <- function(x) {
  contexts <- lapply(x$contexts, function(c) {
    records <- Filter(function(r) r$allo_context_id == c$allo_context_id, x$records)
    dimensions <- lapply(native_allo_dimensions(), function(d) {
      rr <- Filter(function(r) r$dimension == d, records)
      status <- unlist(lapply(rr, function(r) vapply(r$applicability, `[[`, '', 'status')))
      list(dimension = d, record_count = length(rr),
        status_counts = stats::setNames(as.list(vapply(native_allo_states(), function(s) sum(status == s), 0L)), native_allo_states()),
        measured_record_count = sum(vapply(rr, function(r) r$observation$evaluability == 'evaluated', TRUE)),
        missing_record_count = sum(vapply(rr, function(r) r$observation$result_state == 'missing', TRUE)),
        target_count = length(unique(unlist(lapply(rr, function(r) vapply(r$targets, canonical_bytes, ''))))),
        reason_codes = allo_set(unlist(lapply(rr, function(r) unlist(lapply(r$applicability, `[[`, 'reason_codes'), recursive = FALSE)), recursive = FALSE)))
    })
    list(allo_context_id = c$allo_context_id, record_count = length(records), dimensions = dimensions,
      covered_member_count = sum(vapply(c$member_links, function(l) length(l$evidence_ids) > 0L, TRUE)),
      baseline_unchanged = TRUE)
  })
  list(schema_version = 'allo-native-evidence-summary-1.0.0', model_id = x$model_id,
    official_reproduction = FALSE, baseline_unchanged = TRUE,
    aggregation = 'separate_dimensions_and_contexts_no_score', contexts = contexts)
}
native_allo_text <- function(x) paste0('NATIVE ALLO EVIDENCE: ', x$model_id, '; records=', length(x$records),
  '; separate scoped dimensions; literature baseline unchanged; molecular evidence only')
