# T052-T054: observation imports; frozen Phase B namespaces are never extended.
native_allo_version <- function() 'neo-allo-hla-evidence-1.0.0'
native_allo_dimensions <- function() c('source_expression', 'presenter_expression',
  'hla_allele_loh', 'antigen_processing_machinery', 'ms_observation', 'tcr_recognition')
native_allo_states <- function() c('applicable_supported', 'applicable_negative', 'ambiguous',
  'insufficient_scope', 'not_evaluated', 'not_applicable', 'failed')
native_allo_enabled <- function(project) !is.null(project$config$allo_native_evidence)
native_allo_enum <- function(x, choices, label) {
  x <- allo_text(x, label)
  if (!x %in% choices) stop('Invalid native ', label)
  x
}
native_allo_ref <- function(x) {
  allo_keys(x, c('target_type', 'target_id'), label = 'native target')
  x$target_type <- native_allo_enum(x$target_type, c('context', 'source_occurrence',
    'peptide_occurrence', 'peptide', 'pmhc', 'core_assignment', 'aggregate_unit'), 'target type')
  x$target_id <- allo_text(x$target_id, 'frozen target ID')
  x
}
# Inline semantic evidence, never a physical locator or a molecular identity.
native_allo_assignment_evidence <- function(e) {
  if (!is.list(e) || is.null(names(e))) stop('Native assignment evidence must be a structured object')
  common <- c('schema_version', 'assignment_kind', 'peptide_key', 'method', 'provider',
    'method_version', 'provider_version', 'evaluability', 'result_state')
  kind <- native_allo_enum(e$assignment_kind, c('source_origin', 'presenter', 'assay_scope'), 'assignment kind')
  claim <- switch(kind, source_origin = 'source_occurrence_id', presenter = 'presenter',
    assay_scope = 'assay_id')
  allo_keys(e, c(common, claim), label = 'native assignment evidence')
  if (e$schema_version != 'allo-ms-assignment-evidence-v1') stop('Unsupported assignment evidence version')
  for (k in c('peptide_key', 'method', 'provider', 'method_version', 'provider_version', claim))
    e[[k]] <- allo_text(e[[k]], paste('assignment evidence', k))
  e$evaluability <- native_allo_enum(e$evaluability,
    c('evaluated', 'ambiguous', 'not_evaluated', 'failed'), 'assignment evaluability')
  e$result_state <- native_allo_enum(e$result_state,
    c('supported', 'not_supported', 'ambiguous', 'not_evaluated', 'failed'), 'assignment result')
  expected <- if (e$result_state %in% c('supported', 'not_supported')) 'evaluated' else e$result_state
  if (e$evaluability != expected) stop('Assignment result/evaluability conflict')
  e
}
native_allo_row <- function(r) {
  fields <- c('allo_context_id', 'targets', 'dimension', 'subject_id', 'subject_role',
    'tissue', 'compartment', 'sample_id', 'occasion_id', 'assay', 'result_state',
    'evaluability', 'normalized_values', 'analyte', 'assignments', 'reason_codes',
    'provenance', 'raw_values', 'readonly_references')
  allo_keys(r, fields, label = 'native observation', nullable = c('sample_id', 'raw_values'))
  for (k in c('allo_context_id', 'subject_id', 'tissue', 'compartment', 'occasion_id'))
    r[[k]] <- allo_text(r[[k]], k)
  if (!is.null(r$sample_id)) r$sample_id <- allo_text(r$sample_id, 'sample_id')
  r$dimension <- native_allo_enum(r$dimension, native_allo_dimensions(), 'dimension')
  r$subject_role <- native_allo_enum(r$subject_role, c('source', 'presenter', 'responding'), 'subject role')
  r$targets <- allo_records(lapply(allo_records(r$targets), native_allo_ref))
  if (!length(r$targets)) stop('Native observation requires explicit frozen targets')
  allo_keys(r$assay, c('assay_id', 'method', 'version', 'evidence_type', 'receptor_context_id'),
    label = 'native assay', nullable = 'receptor_context_id')
  for (k in c('assay_id', 'method', 'version')) r$assay[[k]] <- allo_text(r$assay[[k]], k)
  r$assay$evidence_type <- native_allo_enum(r$assay$evidence_type, c('expression_assay', 'allele_loh_assay', 'apm_assay', 'mass_spectrometry', 'tcr_recognition_assay', 'binding_prediction', 'processing_prediction'), 'assay evidence type')
  if (!is.null(r$assay$receptor_context_id)) r$assay$receptor_context_id <- allo_text(r$assay$receptor_context_id, 'receptor context')
  r$evaluability <- native_allo_enum(r$evaluability, c('evaluated', 'ambiguous',
    'insufficient_resolution', 'not_evaluated', 'not_applicable', 'failed'), 'evaluability')
  results <- switch(r$dimension,
    hla_allele_loh = c('allele_retained', 'allele_loss_supported'),
    ms_observation = c('detected', 'not_detected'), c('positive', 'negative'))
  r$result_state <- native_allo_enum(r$result_state, c(results, 'missing', 'ambiguous',
    'insufficient_resolution', 'not_evaluated', 'not_applicable', 'failed'), 'result state')
  unresolved <- c('missing', 'not_evaluated', 'not_applicable', 'failed', 'insufficient_resolution')
  expected <- c(missing = 'not_evaluated', not_evaluated = 'not_evaluated',
    not_applicable = 'not_applicable', failed = 'failed', insufficient_resolution = 'insufficient_resolution', ambiguous = 'ambiguous')
  if (r$result_state %in% names(expected) && r$evaluability != expected[[r$result_state]])
    stop('Native result/evaluability conflict')
  if (r$result_state %in% results && !r$evaluability %in% c('evaluated', 'ambiguous', 'insufficient_resolution'))
    stop('Measured native result requires measurement evaluability')
  r$normalized_values <- allo_records(lapply(allo_records(r$normalized_values), function(v) {
    allo_keys(v, c('metric', 'value', 'value_type', 'unit'), label = 'native typed value', nullable = 'unit')
    v$metric <- allo_text(v$metric, 'metric')
    v$value_type <- native_allo_enum(v$value_type, c('number', 'string', 'boolean'), 'value type')
    ok <- switch(v$value_type, number = is.numeric(v$value) && length(v$value) == 1L && is.finite(v$value),
      string = is.character(v$value) && length(v$value) == 1L && !is.na(v$value) && nzchar(v$value),
      boolean = is.logical(v$value) && length(v$value) == 1L && !is.na(v$value))
    if (!ok) stop('Native value type mismatch')
    if (v$value_type == 'number') v$value <- as.numeric(v$value)
    if (v$value_type == 'string') v$value <- allo_text(v$value, 'value')
    if (!is.null(v$unit)) v$unit <- allo_text(v$unit, 'unit')
    if (v$value_type == 'number' && is.null(v$unit)) stop('Numeric native value requires explicit unit')
    v
  }))
  if (r$result_state %in% unresolved && length(r$normalized_values)) stop('Missing/unevaluated evidence cannot contain normalized values')
  if (r$result_state %in% results && is.null(r$sample_id)) stop('Measured evidence requires sample identity')
  if (r$dimension == 'tcr_recognition' && r$result_state %in% results && is.null(r$assay$receptor_context_id))
    stop('TCR assay requires receptor/context provenance')
  allo_keys(r$analyte, c('resolution', 'alleles', 'genes', 'presenter', 'peptide_key',
    'hla_class', 'pathway', 'pathway_applicability'), label = 'native analyte',
    nullable = c('presenter', 'peptide_key', 'hla_class', 'pathway'))
  r$analyte$resolution <- native_allo_enum(r$analyte$resolution,
    c('allele', 'gene', 'protein', 'molecule', 'peptide', 'core', 'pathway', 'unresolved'), 'analyte resolution')
  for (k in c('alleles', 'genes')) r$analyte[[k]] <- allo_strings(r$analyte[[k]], k)
  for (k in c('presenter', 'peptide_key', 'pathway')) if (!is.null(r$analyte[[k]])) r$analyte[[k]] <- allo_text(r$analyte[[k]], k)
  if (!is.null(r$analyte$hla_class)) r$analyte$hla_class <- native_allo_enum(r$analyte$hla_class, c('I', 'II'), 'HLA class')
  r$analyte$pathway_applicability <- native_allo_enum(r$analyte$pathway_applicability,
    c('supported', 'not_evaluated', 'not_applicable'), 'pathway applicability')
  if (r$analyte$resolution == 'allele' && length(r$analyte$alleles) != 1L) stop('Allele resolution requires exactly one allele')
  if (r$analyte$resolution == 'gene' && !length(r$analyte$genes)) stop('Gene resolution requires gene identity')
  if (r$dimension == 'hla_allele_loh' && r$result_state %in% results && r$analyte$resolution != 'allele')
    stop('Assessed LOH retained/loss requires allele-specific resolution')
  if (r$dimension == 'antigen_processing_machinery' && r$result_state %in% results &&
      (!length(r$analyte$genes) && is.null(r$analyte$pathway))) stop('APM requires individual gene/pathway identity')
  allo_keys(r$assignments, c('source_status', 'source_occurrence_ids', 'presenter_status',
    'presenters', 'assignment_evidence', 'absence_scope'), label = 'native assignments')
  for (k in c('source_status', 'presenter_status')) r$assignments[[k]] <- native_allo_enum(r$assignments[[k]],
    c('resolved', 'ambiguous', 'unresolved', 'not_applicable'), k)
  for (k in c('source_occurrence_ids', 'presenters'))
    r$assignments[[k]] <- allo_strings(r$assignments[[k]], k)
  evidence <- allo_records(lapply(allo_records(r$assignments$assignment_evidence), native_allo_assignment_evidence))
  r$assignments$assignment_evidence <- evidence
  r$assignments$absence_scope <- native_allo_enum(r$assignments$absence_scope,
    c('sample_only', 'supplied_assay_exhaustive_contract'), 'absence scope')
  if (r$dimension == 'ms_observation') {
    if (is.null(r$analyte$peptide_key)) stop('MS observation requires observed peptide identity')
    for (e in evidence) {
      if (e$peptide_key != r$analyte$peptide_key) stop('Assignment evidence observed peptide mismatch')
      claims <- switch(e$assignment_kind, source_origin = r$assignments$source_occurrence_ids,
        presenter = r$assignments$presenters, assay_scope = list(r$assay$assay_id))
      claim <- switch(e$assignment_kind, source_origin = e$source_occurrence_id,
        presenter = e$presenter, assay_scope = e$assay_id)
      if (!claim %in% unlist(claims)) stop('Assignment evidence ', e$assignment_kind, ' claim mismatch')
    }
    for (pair in list(c('source_status', 'source_occurrence_ids', 'source_origin'),
                      c('presenter_status', 'presenters', 'presenter'))) {
      supporting <- Filter(function(e) e$assignment_kind == pair[[3L]] &&
        e$evaluability == 'evaluated' && e$result_state == 'supported', evidence)
      if (r$assignments[[pair[[1L]]]] == 'resolved' &&
          (length(r$assignments[[pair[[2L]]]]) != 1L || !length(supporting)))
        stop('Resolved MS ', pair[[3L]], ' assignment requires structured supporting assignment evidence and one identity')
      if (r$assignments[[pair[[1L]]]] == 'unresolved' && length(r$assignments[[pair[[2L]]]]))
        stop('Unresolved MS assignment cannot claim identities')
    }
    if (r$assignments$absence_scope == 'supplied_assay_exhaustive_contract' &&
        !any(vapply(evidence, function(e) e$assignment_kind == 'assay_scope' &&
          e$evaluability == 'evaluated' && e$result_state == 'supported', TRUE)))
      stop('Exhaustive MS absence requires structured supplied assay contract evidence')
  } else if (length(evidence)) stop('Assignment evidence is only valid for MS observations')
  r$reason_codes <- allo_strings(r$reason_codes, 'native reasons')
  r$readonly_references <- allo_strings(r$readonly_references, 'read-only references')
  # Extensible physical/raw assay metadata remains separate from identity.
  for (k in c('provenance', 'raw_values')) {
    if (!is.null(r[[k]])) {
      strict_keys(r[[k]], names(r[[k]]), k)
      canonical_bytes(r[[k]])
    }
  }
  r
}
native_allo_payload <- function(r) {
  list(model_id = 'neo-allo-hla-evidence-v1', allo_context_id = r$allo_context_id,
    targets = r$targets, dimension = r$dimension,
    observation = r[setdiff(names(r), c('allo_context_id', 'targets', 'dimension',
      'provenance', 'raw_values', 'readonly_references'))])
}
# Non-biological observation-record identity; never a peptide, pMHC or member.
native_allo_id <- function(r) indirect_id('allo-native-observation-record-id-v1', native_allo_payload(r))
native_allo_config <- function(x, base = NULL) {
  allo_keys(x, c('version', 'observations', 'scopes'), label = 'native allo config')
  if (x$version != native_allo_version()) stop('Unsupported native allo version')
  allo_keys(x$observations, 'path', c('sha256', 'manifest'), 'native observation binding')
  path <- resolve_path(x$observations$path, base)
  if (!file.exists(path) || dir.exists(path)) stop('Native observation artifact unavailable')
  hash <- digest::digest(file = path, algo = 'sha256', serialize = FALSE)
  if (!is.null(x$observations$sha256) && x$observations$sha256 != hash) stop('Native observation integrity mismatch')
  m <- jsonlite::read_json(path); reject_duplicate_keys(m)
  if (!is.null(x$observations$manifest) && !identical(canonical_bytes(x$observations$manifest), canonical_bytes(m)))
    stop('Native observation content mismatch')
  allo_keys(m, c('schema_version', 'model_id', 'observations'), label = 'native observation import')
  if (m$schema_version != native_allo_version() || m$model_id != 'neo-allo-hla-evidence-v1') stop('Invalid native observation model')
  invisible(lapply(allo_records(m$observations), native_allo_row))
  x$scopes <- allo_records(lapply(allo_records(x$scopes), function(s) {
    allo_keys(s, c('allo_context_id', 'source_compartment', 'responding_compartment',
      'source_sample_ids', 'presenter_sample_ids', 'responding_sample_ids'), label = 'native context scope')
    for (k in c('allo_context_id', 'source_compartment', 'responding_compartment')) s[[k]] <- allo_text(s[[k]], k)
    for (k in c('source_sample_ids', 'presenter_sample_ids', 'responding_sample_ids')) s[[k]] <- allo_strings(s[[k]], k)
    s
  }))
  if (anyDuplicated(vapply(x$scopes, `[[`, '', 'allo_context_id'))) stop('Duplicate native context scope')
  x$observations <- list(path = path, sha256 = hash, manifest = m)
  x
}
# JSON cannot persist R's integer/binary64 distinction. Native opt-in resume
# restores typed snapshots from the *same pinned bytes*, never a live substitute.
# Legacy no-native resume behavior is left untouched.
native_allo_restore_pinned_config <- function(x) {
  if (is.null(x[['allo_native_evidence']])) return(x)
  json_numbers <- function(v) {
    if (is.list(v)) return(lapply(v, json_numbers))
    if (is.integer(v)) return(as.numeric(v))
    v
  }
  restore <- function(binding) {
    if (is.null(binding$manifest) || is.null(binding$sha256)) stop('Native resume requires pinned input snapshot')
    path <- binding$path
    if (!file.exists(path) || dir.exists(path) ||
        digest::digest(file = path, algo = 'sha256', serialize = FALSE) != binding$sha256)
      stop('Native resume pinned artifact integrity mismatch')
    manifest <- jsonlite::read_json(path); reject_duplicate_keys(manifest)
    if (!identical(canonical_bytes(json_numbers(manifest)), canonical_bytes(json_numbers(binding$manifest))))
      stop('Native resume persisted snapshot content mismatch')
    binding$manifest <- manifest
    binding
  }
  x$allo$reference <- restore(x$allo$reference)
  x$indirect_hla$provider <- restore(x$indirect_hla$provider)
  x$allo_native_evidence$observations <- restore(x$allo_native_evidence$observations)
  x
}
