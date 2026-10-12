native_scope <- function(x = phase_c_case()) {
  c <- x$foundation$contexts[[1L]]
  list(allo_context_id = c$allo_context_id, source_compartment = 'declared_source_cells',
    responding_compartment = 'declared_responder_cells', source_sample_ids = list('native-source-sample'),
    presenter_sample_ids = list('native-apc-sample'), responding_sample_ids = list('native-responder-sample'))
}
native_fixture <- function(name) jsonlite::read_json(testthat::test_path('..', 'fixtures', 'native_allo_hla', name))
native_row <- function(dimension = 'source_expression', state = 'positive', x = phase_c_case()) {
  detail <- x$public_explanation$contexts[[1L]]
  u <- detail$included_units[[1L]]
  p <- Neoantimon2:::indirect_report_one(detail$graph$nodes$pmhc, 'id', u$supporting_pmhc_ids[[1L]])$data
  source <- detail$graph$nodes$source_occurrence[[1L]]
  c <- x$foundation$contexts[[1L]]$semantics
  role <- switch(dimension, source_expression = 'source', tcr_recognition = 'responding', 'presenter')
  target <- switch(dimension, source_expression = list(target_type = 'source_occurrence', target_id = source$id),
    antigen_processing_machinery = list(target_type = 'context', target_id = detail$allo_context_id),
    list(target_type = 'pmhc', target_id = p$indirect_pmhc_id))
  a <- list(resolution = 'allele', alleles = list(source$data$semantics$allele_hypothesis), genes = list(),
    presenter = NULL, peptide_key = NULL, hla_class = 'II', pathway = NULL, pathway_applicability = 'not_evaluated')
  if (dimension %in% c('presenter_expression', 'tcr_recognition')) {
    a$resolution <- 'molecule'; a$alleles <- list(); a$presenter <- p$presenter
    if (dimension == 'tcr_recognition') a$peptide_key <- p$peptide_key
  }
  if (dimension == 'hla_allele_loh') a$alleles <- list('HLA-DRB1*01:01')
  if (dimension == 'antigen_processing_machinery') {
    a$resolution <- 'gene'; a$alleles <- list(); a$genes <- list('HLA-DMA')
    a$pathway <- 'class_ii_endosomal_loading'; a$pathway_applicability <- 'supported'
  }
  if (dimension == 'ms_observation') {
    a$resolution <- 'peptide'; a$alleles <- list(); a$peptide_key <- p$peptide_key
  }
  unevaluated <- c('missing', 'failed', 'ambiguous', 'insufficient_resolution', 'not_applicable', 'not_evaluated')
  eval <- if (state == 'missing') 'not_evaluated' else if (state %in% unevaluated) state else 'evaluated'
  list(allo_context_id = detail$allo_context_id, targets = list(target), dimension = dimension,
    subject_id = c[[paste0(role, '_subject_id')]], subject_role = role, tissue = c$tissue,
    compartment = switch(role, source = 'declared_source_cells', presenter = c$apc_compartment, responding = 'declared_responder_cells'),
    sample_id = if (state == 'missing') NULL else paste0('native-', switch(role, source = 'source', presenter = 'apc', responding = 'responder'), '-sample'),
    occasion_id = c$occasion_id,
    assay = list(assay_id = paste0('synthetic-assay:', dimension), method = paste0('synthetic_', dimension),
      version = '1.0.0', evidence_type = switch(dimension, source_expression = 'expression_assay',
        presenter_expression = 'expression_assay', hla_allele_loh = 'allele_loh_assay',
        antigen_processing_machinery = 'apm_assay', ms_observation = 'mass_spectrometry', tcr_recognition = 'tcr_recognition_assay'), receptor_context_id = if (dimension == 'tcr_recognition' && state != 'missing') 'synthetic-tcr:receptor-context-1' else NULL),
    result_state = state, evaluability = eval, normalized_values = list(), analyte = a,
    assignments = list(source_status = 'not_applicable', source_occurrence_ids = list(),
      presenter_status = 'not_applicable', presenters = list(), assignment_evidence = list(), absence_scope = 'sample_only'),
    reason_codes = list('synthetic_fixture_observation'),
    provenance = list(fixture = 'synthetic; no external patient data', path = '/synthetic/native.json', sha256 = paste(rep('a', 64), collapse = '')),
    raw_values = list(assay_note = 'synthetic qualitative evidence'), readonly_references = list())
}
native_assignment_evidence <- function(row, kind) {
  e <- list(schema_version = 'allo-ms-assignment-evidence-v1', assignment_kind = kind,
    peptide_key = row$analyte$peptide_key, method = 'synthetic_independent_assignment',
    provider = 'synthetic_assignment_provider', method_version = '1.0.0', provider_version = '1.0.0',
    evaluability = 'evaluated', result_state = 'supported')
  claim <- switch(kind, source_origin = 'source_occurrence_id', presenter = 'presenter', assay_scope = 'assay_id')
  e[[claim]] <- switch(kind, source_origin = row$assignments$source_occurrence_ids[[1L]],
    presenter = row$assignments$presenters[[1L]], assay_scope = row$assay$assay_id)
  e
}
native_assigned_row <- function(kind = c('source_origin', 'presenter'), x = phase_c_case()) {
  r <- native_row('ms_observation', 'detected', x)
  r$assignments$source_status <- 'unresolved'; r$assignments$presenter_status <- 'unresolved'
  if ('source_origin' %in% kind) {
    occurrences <- Filter(function(n) n$data$semantics$peptide_key == r$analyte$peptide_key,
      x$public_explanation$contexts[[1L]]$graph$nodes$peptide_occurrence)
    r$assignments$source_status <- 'resolved'
    r$assignments$source_occurrence_ids <- list(occurrences[[1L]]$data$semantics$source_occurrence_id)
  }
  if ('presenter' %in% kind) {
    r$assignments$presenter_status <- 'resolved'
    r$assignments$presenters <- list(native_row('presenter_expression', x = x)$analyte$presenter)
  }
  r$assignments$assignment_evidence <- lapply(kind, function(k) native_assignment_evidence(r, k))
  r$provenance$assignment_evidence <- list(path = '/synthetic/assignment.tsv',
    sha256 = paste(rep('b', 64), collapse = ''), acquisition_timestamp = 'synthetic-acquisition',
    retrieval_timestamp = 'synthetic-retrieval')
  r
}
native_import <- function(rows) list(schema_version = 'neo-allo-hla-evidence-1.0.0',
  model_id = 'neo-allo-hla-evidence-v1', observations = unname(rows))
native_enable <- function(p = phase_c_project(), rows = native_fixture('complete-observations.json')$observations,
                          scopes = list(native_scope())) {
  binding <- list(path = indirect_write(native_import(rows)))
  p$config$allo_native_evidence <- Neoantimon2:::native_allo_config(list(version = 'neo-allo-hla-evidence-1.0.0', observations = binding, scopes = scopes))
  p$allo_native_evidence <- Neoantimon2:::native_allo_build(p)
  p
}
native_status <- function(row) native_enable(rows = list(row))$allo_native_evidence$records[[1L]]$applicability[[1L]]$status
native_semantic <- function(x) {
  x$provenance <- NULL
  for (i in seq_along(x$records)) x$records[[i]]$raw_observations <- NULL
  x
}
native_schema <- function(name = 'neo-allo-hla-evidence', definition = NULL) {
  s <- jsonlite::read_json(system.file('schema', paste0(name, '-1.0.0.json'), package = 'Neoantimon2'))
  if (!is.null(definition)) { s$oneOf <- NULL; s$`$ref` <- paste0('#/$defs/', definition) }
  v <- jsonvalidate::json_validator(jsonlite::toJSON(s, auto_unbox = TRUE, null = 'null'), engine = 'ajv')
  function(x) as.logical(v(jsonlite::toJSON(Neoantimon2:::longitudinal_plain(x), auto_unbox = TRUE, null = 'null', digits = NA)))
}
# Twenty-seven explicit biological-scope or missingness cases; expected states
# are independent assertions, not computed by the applicability implementation.
native_focused_rows <- function() {
  rows <- list()
  add <- function(name, row, expected) rows[[name]] <<- list(row = row, expected = expected)
  s <- native_row(); s$normalized_values <- list(list(metric = 'expression', value = 12.5, value_type = 'number', unit = 'TPM'))
  add('01 source expression positive', s, 'applicable_supported')
  s <- native_row(state = 'negative'); s$normalized_values <- list(list(metric = 'expression', value = 0.1, value_type = 'number', unit = 'TPM'))
  add('02 source expression measured low', s, 'applicable_negative')
  s <- native_row(); s$analyte$resolution <- 'gene'; s$analyte$genes <- list(sub('\\*.*', '', s$analyte$alleles[[1L]])); s$analyte$alleles <- list()
  add('03 gene expression allele ambiguity', s, 'ambiguous')
  add('04 missing source expression', native_row(state = 'missing'), 'not_evaluated')
  s <- native_row(); s$targets[[1L]]$target_id <- phase_c_detail()$graph$nodes$source_occurrence[[2L]]$id
  add('05 expression scoped to different source origin', s, 'not_applicable')
  add('06 APC presenter expression', native_row('presenter_expression'), 'applicable_supported')
  add('07 missing presenter expression', native_row('presenter_expression', 'missing'), 'not_evaluated')
  s <- native_row('presenter_expression'); s$subject_role <- 'source'; s$subject_id <- 'donor-1'
  add('08 source cannot satisfy presenter expression', s, 'not_applicable')
  s <- native_row('presenter_expression'); s$compartment <- 'unrelated_cells'
  add('09 wrong APC compartment', s, 'not_applicable')
  add('10 assessed allele retained', native_row('hla_allele_loh', 'allele_retained'), 'applicable_supported')
  add('11 assessed allele loss supported', native_row('hla_allele_loh', 'allele_loss_supported'), 'applicable_negative')
  s <- native_row('hla_allele_loh', 'ambiguous'); s$analyte$resolution <- 'unresolved'; s$analyte$alleles <- list('HLA-DRB1*01:01', 'HLA-DRB1*03:01')
  add('12 ambiguous LOH allele assignment', s, 'ambiguous')
  s <- native_row('hla_allele_loh', 'allele_loss_supported'); s$compartment <- 'tumor_cells'
  add('13 tumor LOH not global transplant APC evidence', s, 'not_applicable')
  add('14 low expression remains expression only', native_row(state = 'negative'), 'applicable_negative')
  add('15 individual Class II APM gene', native_row('antigen_processing_machinery'), 'applicable_supported')
  s <- native_row('antigen_processing_machinery'); s$analyte$genes <- list('TAP1'); s$analyte$hla_class <- 'I'
  add('16 unsupported Class I machinery for Class II', s, 'not_applicable')
  add('17 missing APM not negative', native_row('antigen_processing_machinery', 'missing'), 'not_evaluated')
  ms <- native_row('ms_observation', 'detected'); ms$assignments$source_status <- 'unresolved'; ms$assignments$presenter_status <- 'unresolved'
  s <- ms; s$assignments$presenter_status <- 'resolved'; s$assignments$presenters <- list(native_row('presenter_expression')$analyte$presenter); s$assignments$assignment_evidence <- list(native_assignment_evidence(s, 'presenter'))
  add('18 MS source unresolved', s, 'ambiguous')
  s <- ms; s$assignments$source_status <- 'resolved'; s$assignments$source_occurrence_ids <- list(native_row()$targets[[1L]]$target_id); s$assignments$assignment_evidence <- list(native_assignment_evidence(s, 'source_origin'))
  add('19 MS presenter unresolved', s, 'ambiguous')
  add('20 MS both assignments unresolved', ms, 'ambiguous')
  s <- ms; s$assignments$source_status <- 'resolved'; s$assignments$presenter_status <- 'resolved'
  s$assignments$source_occurrence_ids <- list(native_row()$targets[[1L]]$target_id)
  s$assignments$presenters <- list(native_row('presenter_expression')$analyte$presenter)
  s$assignments$assignment_evidence <- lapply(c('source_origin', 'presenter'), function(k) native_assignment_evidence(s, k))
  add('21 uniquely resolved MS supplied assignment', s, 'applicable_supported')
  s <- ms; s$result_state <- 'not_detected'
  add('22 MS nondetection not exhaustive absence', s, 'insufficient_scope')
  add('23 responder positive TCR assay', native_row('tcr_recognition'), 'applicable_supported')
  add('24 responder evaluated negative TCR assay', native_row('tcr_recognition', 'negative'), 'applicable_negative')
  add('25 missing TCR not no recognition', native_row('tcr_recognition', 'missing'), 'not_evaluated')
  s <- native_row('tcr_recognition'); s$subject_id <- 'donor-1'
  add('26 wrong responder TCR assay', s, 'not_applicable')
  add('27 MS detection remains MS only', ms, 'ambiguous')
  rows
}
# Explicit acceptance cases for the assignment-evidence review fix.
native_assignment_cases <- function() {
  cases <- list()
  add <- function(name, row, runtime_valid, schema_valid = runtime_valid, graph = 'complete')
    cases[[name]] <<- list(name = name, graph = graph, observation = row,
      runtime_valid = runtime_valid, schema_valid = schema_valid)
  for (kind in c('source_origin', 'presenter')) {
    r <- native_assigned_row(kind); r$assignments$assignment_evidence <- list()
    add(paste('no evidence', kind), r, FALSE)
    for (locator in c('/tmp/assignment.tsv', 'file:///tmp/assignment.tsv',
      paste(rep('a', 64), collapse = ''), 'run-2026-10-06/result.json')) {
      r$assignments$assignment_evidence <- list(locator)
      add(paste('string alone', kind, locator), r, FALSE)
    }
    r <- native_assigned_row(kind); r$assignments$assignment_evidence <- list(list(path = '/tmp/assignment.tsv'))
    add(paste('locator object alone', kind), r, FALSE)
    r <- native_assigned_row(kind); r$assignments$assignment_evidence <- list()
    r$assignments$assignment_evidence_refs <- list('/tmp/assignment.tsv')
    add(paste('legacy string escape rejected', kind), r, FALSE)
    r <- native_assigned_row(kind)
    r$assignments$assignment_evidence[[1L]]$result_state <- 'not_supported'
    add(paste('non-supporting result', kind), r, FALSE)
    r$assignments$assignment_evidence[[1L]]$evaluability <- 'failed'
    r$assignments$assignment_evidence[[1L]]$result_state <- 'failed'
    add(paste('failed assignment', kind), r, FALSE)
    r <- native_assigned_row(kind)
    r$assignments$assignment_evidence[[1L]]$assignment_kind <- 'assay_scope'
    claim <- if (kind == 'source_origin') 'source_occurrence_id' else 'presenter'
    r$assignments$assignment_evidence[[1L]][[claim]] <- NULL
    r$assignments$assignment_evidence[[1L]]$assay_id <- r$assay$assay_id
    add(paste('assay scope cannot resolve', kind), r, FALSE)
  }
  r <- native_assigned_row('source_origin'); ids <- lapply(phase_c_detail()$graph$nodes$source_occurrence, `[[`, 'id')
  r$assignments$source_occurrence_ids <- ids[!vapply(ids, identical, TRUE, r$assignments$source_occurrence_ids[[1L]])][1L]
  add('source A evidence observation B', r, FALSE, TRUE)
  x <- phase_c_case('two_presenters'); r <- native_assigned_row('presenter', x)
  presenters <- unique(vapply(x$public_explanation$contexts[[1L]]$graph$nodes$pmhc, function(n) n$data$presenter, ''))
  r$assignments$presenters <- as.list(setdiff(presenters, r$assignments$presenters[[1L]])[1L])
  add('presenter A evidence observation B', r, FALSE, TRUE, 'two_presenters')
  r <- native_assigned_row('source_origin'); r$assignments$assignment_evidence[[1L]]$peptide_key <- native_row()$targets[[1L]]$target_id
  add('assignment observed peptide mismatch', r, FALSE, TRUE)
  r <- native_assigned_row('source_origin'); r$assignments$source_occurrence_ids <- list('invalid:origin')
  r$assignments$assignment_evidence[[1L]]$source_occurrence_id <- 'invalid:origin'
  add('unresolved frozen source reference', r, FALSE, TRUE)
  r <- native_assigned_row('presenter'); r$assignments$presenters <- list('HLA-DRA*01:01~HLA-DRB1*99:99')
  r$assignments$assignment_evidence[[1L]]$presenter <- r$assignments$presenters[[1L]]
  add('unresolved frozen presenter reference', r, FALSE, TRUE)
  for (k in names(native_assigned_row('source_origin')$assignments$assignment_evidence[[1L]])) {
    r <- native_assigned_row('source_origin'); r$assignments$assignment_evidence[[1L]][[k]] <- NULL
    add(paste('required semantic field', k), r, FALSE)
  }
  for (k in c('path', 'sha256', 'acquisition_timestamp', 'retrieval_timestamp', 'raw_provider_output')) {
    r <- native_assigned_row('source_origin'); r$assignments$assignment_evidence[[1L]][[k]] <- '/tmp/physical-locator'
    add(paste('physical field excluded', k), r, FALSE)
  }
  for (kind in c('source_origin', 'presenter')) add(paste('independent evidence', kind), native_assigned_row(kind), TRUE)
  add('both independently evidenced', native_assigned_row(), TRUE)
  r <- native_assigned_row('source_origin'); r$assignments$presenter_status <- 'resolved'
  r$assignments$presenters <- list(native_row('presenter_expression')$analyte$presenter)
  add('source evidence cannot resolve presenter', r, FALSE)
  r <- native_assigned_row('presenter'); r$assignments$source_status <- 'resolved'
  r$assignments$source_occurrence_ids <- list(native_assigned_row('source_origin')$assignments$source_occurrence_ids[[1L]])
  add('presenter evidence cannot resolve source', r, FALSE)
  r <- native_assigned_row(); r$assignments$assignment_evidence[[1L]]$schema_version <- 'unknown-v2'
  add('unknown evidence version', r, FALSE)
  r <- native_assigned_row(); r$assignments$absence_scope <- 'supplied_assay_exhaustive_contract'
  add('assignment evidence cannot prove exhaustive assay scope', r, FALSE)
  r$assignments$assignment_evidence <- c(r$assignments$assignment_evidence, list(native_assignment_evidence(r, 'assay_scope')))
  add('explicit supplied exhaustive assay contract', r, TRUE)
  r <- native_row(); r$assignments$assignment_evidence <- native_assigned_row()$assignments$assignment_evidence
  add('assignment evidence restricted to MS', r, FALSE)
  cases
}
