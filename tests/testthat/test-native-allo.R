# T052-T054; frozen tests and fixtures are deliberately untouched.
for (case in native_fixture('focused-observations.json')$cases) local({
  x <- case
  test_that(x$name, {
    p <- native_enable(rows = list(x$observation)); r <- p$allo_native_evidence$records[[1L]]
    expect_identical(r$applicability[[1L]]$status, x$expected)
    expect_identical(r$dimension, x$observation$dimension)
    expect_identical(r$observation$result_state, x$observation$result_state)
    expect_identical(r$observation$evaluability, x$observation$evaluability)
    expect_identical(p$indirect_hla, phase_c_case()$result)
    expect_identical(p$indirect_hla$contexts[[1L]]$summary$complete_count, 5L)
    expect_true(native_schema()(p$allo_native_evidence))
    expect_true(native_schema(definition = 'import')(native_import(list(x$observation))))
  })
})
for (name in c('28 no native', '29 native added', '30 native removed', '31 native reordered', '32 native positive versus negative')) local({
  n <- name
  test_that(n, {
    before <- phase_c_project(); rows <- native_fixture('complete-observations.json')$observations
    if (startsWith(n, '28')) after <- before
    if (startsWith(n, '29')) after <- native_enable(before, rows)
    if (startsWith(n, '30')) after <- native_enable(before, list())
    if (startsWith(n, '31')) after <- native_enable(before, rev(rows))
    if (startsWith(n, '32')) { rows[[1L]]$result_state <- 'negative'; after <- native_enable(before, rows) }
    expect_identical(serialize(after$indirect_hla, NULL), serialize(before$indirect_hla, NULL))
    expect_identical(after$allo_foundation, before$allo_foundation)
    expect_identical(neo_explain(after)$indirect_hla_explanation, neo_explain(before)$indirect_hla_explanation)
    expect_identical(attr(neo_results(after), 'indirect_hla_result'), attr(neo_results(before), 'indirect_hla_result'))
    expect_identical(after$analysis_id, before$analysis_id)
    expect_identical(after$indirect_hla$allo_policy_id, before$indirect_hla$allo_policy_id)
    expect_identical(after$indirect_hla$allo_analysis_id, before$indirect_hla$allo_analysis_id)
    expect_identical(after$indirect_hla$contexts[[1L]]$summary$member_ids, before$indirect_hla$contexts[[1L]]$summary$member_ids)
  })
})
test_that('33 compact native results summary separates six dimensions', {
  p <- native_enable(); s <- attr(neo_results(p), 'allo_native_evidence')
  expect_length(s$contexts, 1L); expect_length(s$contexts[[1L]]$dimensions, 6L)
  expect_identical(s$contexts[[1L]]$record_count, 6L); expect_true(s$baseline_unchanged)
  expect_identical(s$model_id, 'neo-allo-hla-evidence-v1')
  expect_true(native_schema('allo-native-evidence-summary')(s))
  expect_false(any(c('records', 'normalized_values', 'raw_observations') %in% names(s)))
})
test_that('34 detailed native explain overlay preserves frozen explanation', {
  p <- native_enable(); e <- neo_explain(p)
  expect_identical(e$allo_native_evidence, p$allo_native_evidence)
  expect_identical(e$indirect_hla_explanation, phase_c_case()$public_explanation)
  expect_length(e$allo_native_evidence$records, 6L)
  for (r in e$allo_native_evidence$records) {
    expect_true(all(c('subject_id','subject_role','tissue','compartment','sample_id','occasion_id','assay',
      'result_state','evaluability','normalized_values','analyte','assignments','reason_codes') %in% names(r$observation)))
    expect_length(r$raw_observations, 1L)
  }
  expect_match(neo_explain(p, format = 'text'), 'literature baseline unchanged')
})
test_that('35 member reverse links retain association without spreading support', {
  p <- native_enable(); o <- p$allo_native_evidence; links <- o$contexts[[1L]]$member_links
  expect_length(links, 5L)
  source <- Filter(function(r) r$dimension == 'source_expression', o$records)[[1L]]
  apm <- Filter(function(r) r$dimension == 'antigen_processing_machinery', o$records)[[1L]]
  expect_true(all(vapply(links, function(l) source$evidence_id %in% unlist(l$evidence_ids), TRUE)))
  expect_false(any(vapply(links, function(l) apm$evidence_id %in% unlist(l$evidence_ids), TRUE)))
  expect_silent(Neoantimon2:::native_allo_integrity(o, phase_c_case()$public_explanation))
})
test_that('36 no allo project unchanged', {
  p <- phase_c_runtime_project()
  expect_null(attr(neo_results(p), 'allo_native_evidence'))
  expect_null(neo_explain(p, peptide = 'ACDEFGHIK')$allo_native_evidence)
  expect_identical(neo_results(p), Neoantimon2:::allo_results_append(neo_results(p), p))
})
test_that('37 Phase A only project unchanged', {
  p <- allo_test_enable(phase_c_runtime_project())
  expect_null(attr(neo_results(p), 'allo_native_evidence')); expect_null(neo_explain(p)$allo_native_evidence)
  expect_identical(Neoantimon2:::allo_read(p), p$allo_foundation)
})
test_that('38 Phase BC without native unchanged', {
  p <- phase_c_project()
  expect_identical(neo_explain(p)$indirect_hla_explanation, phase_c_case()$public_explanation)
  expect_null(neo_explain(p)$allo_native_evidence); expect_null(attr(neo_results(p), 'allo_native_evidence'))
})
test_that('39 multiple native allo contexts remain separate', {
  # Existing strict baseline accepts only one SOT scenario. Test the collection
  # projection with two independent valid contexts; do not invent a D1 runtime.
  one <- native_enable()$allo_native_evidence
  two <- one; two$contexts[[1L]]$allo_context_id <- 'independent-context:2'
  for (i in seq_along(two$records)) two$records[[i]]$allo_context_id <- 'independent-context:2'
  both <- one; both$records <- c(one$records, two$records); both$contexts <- c(one$contexts, two$contexts)
  s <- Neoantimon2:::native_allo_summary(both)
  expect_length(s$contexts, 2L); expect_identical(vapply(s$contexts, `[[`, 0L, 'record_count'), c(6L,6L))
  expect_identical(s$aggregation, 'separate_dimensions_and_contexts_no_score')
  bad <- native_row(); bad$allo_context_id <- 'independent-context:2'
  expect_error(native_enable(rows = list(bad)), 'context reference unresolved')
})
test_that('40 physical path changes leave evidence semantic identity unchanged', {
  row <- native_row(); p <- native_enable(rows = list(row)); row$provenance$path <- '/different/native.json'
  q <- native_enable(rows = list(row))
  expect_identical(native_semantic(p$allo_native_evidence), native_semantic(q$allo_native_evidence))
  expect_false(identical(p$allo_native_evidence$records[[1L]]$raw_observations, q$allo_native_evidence$records[[1L]]$raw_observations))
})
test_that('41 physical checksum raw metadata changes leave semantic ID unchanged', {
  r <- native_row(); p <- native_enable(rows = list(r)); r$provenance$sha256 <- paste(rep('b',64),collapse='')
  r$raw_values$extra_provider_metadata <- list(arbitrary = TRUE)
  q <- native_enable(rows = list(r))
  expect_identical(native_semantic(p$allo_native_evidence), native_semantic(q$allo_native_evidence))
})
test_that('42 semantic assay and result changes change evidence ID', {
  row <- native_row(); old <- native_enable(rows = list(row))$allo_native_evidence$records[[1L]]$evidence_id
  for (k in c('result_state','assay','normalized_values','sample_id','occasion_id','subject_id','tissue','compartment')) {
    r <- row
    if (k == 'result_state') r$result_state <- 'negative'
    if (k == 'assay') r$assay$version <- '2.0.0'
    if (k == 'normalized_values') r$normalized_values <- list(list(metric='expression', value=1.0,value_type='number',unit='TPM'))
    if (k %in% c('sample_id','occasion_id','subject_id','tissue','compartment')) r[[k]] <- paste0(r[[k]], '-other')
    expect_false(identical(native_enable(rows=list(r))$allo_native_evidence$records[[1L]]$evidence_id, old))
  }
})
test_that('43 input row order preserves semantic overlay', {
  rows <- native_fixture('complete-observations.json')$observations
  expect_identical(native_semantic(native_enable(rows=rows)$allo_native_evidence), native_semantic(native_enable(rows=rev(rows))$allo_native_evidence))
})
test_that('44 duplicate semantic observations normalize without hiding raw provenance', {
  row <- native_row(); other <- row; other$provenance$path <- '/other/duplicate.json'
  p <- native_enable(rows=list(row,other,row)); q <- native_enable(rows=list(other,row,row))
  expect_length(p$allo_native_evidence$records, 1L)
  expect_length(p$allo_native_evidence$records[[1L]]$raw_observations, 2L)
  expect_identical(native_semantic(p$allo_native_evidence), native_semantic(q$allo_native_evidence))
  expect_identical(p$allo_native_evidence$records, q$allo_native_evidence$records)
})
test_that('native namespace exact payload guards exclude physical provenance', {
  r <- Neoantimon2:::native_allo_row(native_row()); payload <- Neoantimon2:::native_allo_payload(r)
  expect_identical(names(payload), c('model_id','allo_context_id','targets','dimension','observation'))
  expect_false(any(c('provenance','raw_values','readonly_references') %in% names(payload$observation)))
  expect_identical(payload$observation$assay, r$assay)
  expect_identical(payload$observation$occasion_id, r$occasion_id)
  expect_match(native_enable(rows=list(r))$allo_native_evidence$records[[1L]]$evidence_id, '^allo-native-observation-record-id-v1:')
  expect_identical(Neoantimon2:::allo_namespaces(), c(context='allo-context-id-v1',sequence='hla-source-sequence-id-v1',
    occurrence='hla-source-occurrence-id-v1',universe='self-universe-id-v1',policy='allo-policy-id-v1',analysis='allo-analysis-id-v1'))
})
test_that('native missing negative ambiguous failed remain distinct and typed', {
  for (state in c('missing','negative','ambiguous','failed','insufficient_resolution','not_applicable')) {
    r <- native_row(state=state); p <- native_enable(rows=list(r)); out <- p$allo_native_evidence$records[[1L]]
    expect_identical(out$observation$result_state, state)
    expect_true(native_schema()(p$allo_native_evidence))
  }
  r <- native_row(state='missing'); r$normalized_values <- list(list(metric='expression',value=0,value_type='number',unit='TPM'))
  expect_error(native_enable(rows=list(r)), 'cannot contain normalized values')
  expect_false(native_schema(definition='observation')(r))
  r <- native_row(state='failed'); r$evaluability <- 'evaluated'
  expect_error(native_enable(rows=list(r)), 'result/evaluability conflict')
})
test_that('native subjects tissues occasions and gene resolution cannot leak', {
  for (k in c('subject_id','tissue','compartment','sample_id','occasion_id')) {
    r <- native_row('presenter_expression'); r[[k]] <- 'incompatible-scope'
    expect_true(native_status(r) %in% c('not_applicable','insufficient_scope'))
  }
  r <- native_row('presenter_expression'); r$analyte$resolution <- 'gene'; r$analyte$genes <- list('HLA-DRB1')
  r$analyte['presenter'] <- list(NULL)
  expect_identical(native_status(r), 'insufficient_scope')
  r$targets <- list(list(target_type='context',target_id=r$allo_context_id))
  expect_identical(native_status(r), 'applicable_supported')
  r <- native_row('antigen_processing_machinery'); r$analyte$genes <- list('TAP1')
  expect_identical(native_status(r), 'not_applicable')
})
test_that('one ambiguous MS observation links multiple candidates without duplication', {
  x <- phase_c_case('two_presenters'); r <- native_row('ms_observation','detected',x)
  pp <- Filter(function(p) p$peptide_key == r$analyte$peptide_key,x$result$contexts[[1L]]$pmhcs)
  r$targets <- lapply(pp,function(p)list(target_type='pmhc',target_id=p$indirect_pmhc_id))
  r$assignments$source_status <- 'unresolved'; r$assignments$presenter_status <- 'unresolved'
  p <- native_enable(phase_c_project('two_presenters'),list(r),list(native_scope(x)))
  expect_length(p$allo_native_evidence$records,1L); expect_length(p$allo_native_evidence$records[[1L]]$targets,2L)
  expect_true(all(vapply(p$allo_native_evidence$records[[1L]]$applicability,function(a)a$status=='ambiguous',TRUE)))
  r$targets <- rev(r$targets)
  q <- native_enable(phase_c_project('two_presenters'),list(r),list(native_scope(x)))
  expect_identical(native_semantic(p$allo_native_evidence),native_semantic(q$allo_native_evidence))
})
test_that('source origins keep independent expression including measured negatives', {
  a <- native_row(); b <- a; source <- phase_c_detail()$graph$nodes$source_occurrence[[2L]]
  b$targets[[1L]]$target_id <- source$id; b$analyte$alleles <- list(source$data$semantics$allele_hypothesis); b$result_state <- 'negative'
  p <- native_enable(rows=list(a,b)); r <- p$allo_native_evidence$records
  expect_length(r,2L); expect_setequal(vapply(r,function(x)x$applicability[[1L]]$status,''),c('applicable_supported','applicable_negative'))
  expect_identical(p$indirect_hla,phase_c_case()$result)
})
test_that('native evidence cannot complete incomplete literature counts', {
  for (name in c('incomplete_self','incomplete_scope','missing_source','provider_failure')) {
    p <- phase_c_project(name); rows <- native_fixture('complete-observations.json')$observations
    # Use a context-only APM observation valid even without source generation.
    rows <- Filter(function(r)r$dimension=='antigen_processing_machinery',rows)
    q <- native_enable(p,rows)
    expect_null(attr(neo_results(q),'indirect_hla_result')$contexts[[1L]]$complete_count)
    expect_identical(q$indirect_hla,p$indirect_hla)
  }
})
test_that('native overlay integrity rejects unresolved targets and tampering', {
  r <- native_row(); r$targets[[1L]]$target_id <- 'invalid:frozen-id'
  expect_error(native_enable(rows=list(r)), 'resolve exactly once')
  p <- native_enable(); p$allo_native_evidence$records[[1L]]$observation$result_state <- 'negative'
  expect_error(neo_results(p), 'derived contract mismatch')
  p <- native_enable(); p$indirect_hla$contexts[[1L]]$summary$complete_count <- 4L
  expect_error(neo_explain(p), 'baseline derived contract mismatch')
  r <- native_row('tcr_recognition'); r$assay['receptor_context_id'] <- list(NULL)
  expect_error(native_enable(rows=list(r)), 'receptor/context provenance')
})
test_that('native config binds local inputs and requires opt-in baseline', {
  x <- phase_c_case(); rows <- native_fixture('complete-observations.json')
  f <- indirect_write(rows); cfg <- list(version='neo-allo-hla-evidence-1.0.0', observations=list(path=f), scopes=list(native_scope()))
  expect_error(neo_config(allo_native_evidence=cfg), 'requires frozen indirect HLA')
  config <- neo_config(allo=x$input$a, indirect_hla=x$input$b, allo_native_evidence=cfg)
  p <- neo_project(config=config)
  expect_identical(p$indirect_hla,x$result); expect_length(p$allo_native_evidence$records,6L)
  cfg$observations$sha256 <- paste(rep('0',64),collapse='')
  expect_error(Neoantimon2:::native_allo_config(cfg), 'integrity mismatch')
  expect_setequal(getNamespaceExports('Neoantimon2'), current_public_api)
})
test_that('no combined native score rank or clinical output exists', {
  p <- native_enable(); before <- neo_results(phase_c_project()); after <- neo_results(p)
  attr(after,'allo_native_evidence') <- NULL
  expect_identical(after,before)
  text <- jsonlite::toJSON(p$allo_native_evidence, auto_unbox=TRUE,null='null')
  expect_false(grepl('enhanced_score|adjusted_PIRCHE|risk_score|priority_score|weighted_alloreactivity|rejection_probability|GVHD_probability|clinical_cutoff',text))
  expect_false(any(c('score','rank','grade','priority','weight') %in% names(p$allo_native_evidence)))
  for(f in c('native_allo_build','native_allo_summary','native_allo_applicability')) {
    code <- deparse(get(f, envir = asNamespace('Neoantimon2')))
    expect_false(any(grepl('weighted[.]mean|order.*(score|weight)|enhanced_score|risk_score',code)))
  }
})
test_that('binding prediction and MS cannot act as TCR recognition assay', {
  for (kind in c('binding_prediction','mass_spectrometry','expression_assay')) {
    r <- native_row('tcr_recognition'); r$assay$evidence_type <- kind
    expect_identical(native_status(r), 'not_applicable')
  }
  p <- native_enable(rows = list(native_row(), native_row('ms_observation','detected')))
  expect_false(any(vapply(p$allo_native_evidence$records,function(r)r$dimension=='tcr_recognition',TRUE)))
})
test_that('resolved MS requires independent assignment and frozen references', {
  r <- native_row('ms_observation','detected'); r$assignments$source_status <- 'resolved'
  r$assignments$source_occurrence_ids <- list(native_row()$targets[[1L]]$target_id)
  expect_error(native_enable(rows=list(r)), 'structured supporting assignment evidence')
  r$assignments$source_occurrence_ids <- list('invalid:origin')
  r$assignments$assignment_evidence <- list(native_assignment_evidence(r, 'source_origin'))
  expect_error(native_enable(rows=list(r)), 'resolve exactly once')
  r <- native_row('ms_observation','detected')
  r$targets <- list(list(target_type='peptide',target_id=r$analyte$peptide_key))
  r$assignments$source_status <- 'unresolved'; r$assignments$presenter_status <- 'unresolved'
  expect_identical(native_status(r), 'applicable_supported')
  p <- native_enable(rows=list(r))
  expect_identical(p$allo_native_evidence$records[[1L]]$observation$assignments$presenter_status,'unresolved')
})
test_that('native observation normalized units types and missingness are guarded', {
  for (v in list(list(metric='expression',value=NULL,value_type='number',unit='TPM'),
    list(metric='expression',value=0,value_type='number',unit=NULL),
    list(metric='expression',value='low',value_type='number',unit='TPM'))) {
    r <- native_row(); r$normalized_values <- list(v)
    expect_error(native_enable(rows=list(r)))
    expect_false(native_schema(definition='observation')(r))
  }
  r <- native_row('hla_allele_loh','allele_loss_supported'); r$analyte$resolution <- 'gene'; r$analyte$genes <- list('HLA-DRB1')
  expect_error(native_enable(rows=list(r)), 'allele-specific resolution')
  expect_false(native_schema(definition='observation')(r))
})
test_that('native overlay resume rejects observation configuration changes', {
  x <- phase_c_case(); f <- indirect_write(native_fixture('complete-observations.json'))
  cfg <- list(version='neo-allo-hla-evidence-1.0.0', observations=list(path=f), scopes=list(native_scope()))
  config <- neo_config(env=character(), paths=list(run_root=file.path(tempdir(),'native-allo-runs'),data_root=file.path(tempdir(),'native-allo-data')),
    allo=x$input$a, indirect_hla=x$input$b, allo_native_evidence=cfg)
  p <- neo_project(config=config,initialize=TRUE)
  q <- neo_project(path=p$path)
  expect_identical(q$allo_native_evidence,p$allo_native_evidence)
  cfg$observations$path <- indirect_write(native_import(list(native_row(state='negative'))))
  other <- neo_config(env=character(),paths=config$paths,allo=x$input$a,indirect_hla=x$input$b,allo_native_evidence=cfg)
  expect_error(neo_project(path=p$path,config=other), 'native allo evidence configuration mismatch')
})
test_that('native resume typed snapshot restoration preserves pinning and legacy opt-out', {
  p <- native_enable(); cfg <- unclass(p$config)
  json_cfg <- jsonlite::fromJSON(jsonlite::toJSON(cfg,auto_unbox=TRUE,null='null',digits=NA),simplifyVector=FALSE)
  restored <- Neoantimon2:::native_allo_restore_pinned_config(json_cfg)
  expect_identical(Neoantimon2:::canonical_bytes(restored$indirect_hla$provider$manifest),
    Neoantimon2:::canonical_bytes(cfg$indirect_hla$provider$manifest))
  legacy <- json_cfg; legacy$allo_native_evidence <- NULL
  expect_identical(Neoantimon2:::native_allo_restore_pinned_config(legacy),legacy)
  bad <- json_cfg; bad$indirect_hla$provider$sha256 <- paste(rep('0',64),collapse='')
  expect_error(Neoantimon2:::native_allo_restore_pinned_config(bad),'integrity mismatch')
  bad <- json_cfg;bad$indirect_hla$provider$manifest$rows[[1L]]$affinity_nM <- 900
  expect_error(Neoantimon2:::native_allo_restore_pinned_config(bad),'snapshot content mismatch')
})

for (case in native_assignment_cases()) local({
  c <- case
  test_that(paste('assignment evidence:', c$name), {
    validate <- function() native_enable(phase_c_project(c$graph), list(c$observation), list(native_scope(phase_c_case(c$graph))))
    if (c$runtime_valid) {
      p <- validate()
      expect_length(p$allo_native_evidence$records, 1L)
      expect_true(native_schema()(p$allo_native_evidence))
    } else expect_error(validate())
    expect_identical(native_schema(definition = 'observation')(c$observation), c$schema_valid)
    semantic <- c$observation[c(setdiff(names(c$observation), c('provenance', 'raw_values', 'readonly_references')))]
    expect_identical(native_schema(definition = 'semantic_observation')(semantic), c$schema_valid)
  })
})

test_that('assignment evidence source and presenter certainty remain independent', {
  for (kind in c('source_origin', 'presenter')) {
    r <- native_assigned_row(kind); p <- native_enable(rows = list(r))
    a <- p$allo_native_evidence$records[[1L]]$observation$assignments
    expect_identical(a$source_status, if (kind == 'source_origin') 'resolved' else 'unresolved')
    expect_identical(a$presenter_status, if (kind == 'presenter') 'resolved' else 'unresolved')
    expect_length(a$assignment_evidence, 1L)
    expect_identical(p$allo_native_evidence$records[[1L]]$applicability[[1L]]$status, 'ambiguous')
  }
})

test_that('assignment evidence physical provenance leaves semantic identity unchanged', {
  r <- native_assigned_row(); p <- native_enable(rows = list(r))
  id <- p$allo_native_evidence$records[[1L]]$evidence_id
  for (k in c('path', 'sha256', 'acquisition_timestamp', 'retrieval_timestamp', 'raw_provider_output')) {
    q <- r; q$provenance$assignment_evidence[[k]] <- paste0('changed-physical-', k)
    after <- native_enable(rows = list(q))
    expect_identical(after$allo_native_evidence$records[[1L]]$evidence_id, id)
    expect_identical(native_semantic(after$allo_native_evidence), native_semantic(p$allo_native_evidence))
    merged <- native_enable(rows = list(r, q))$allo_native_evidence$records
    expect_length(merged, 1L); expect_length(merged[[1L]]$raw_observations, 2L)
  }
  payload <- Neoantimon2:::native_allo_payload(Neoantimon2:::native_allo_row(r))
  expect_identical(payload$observation$assignments$assignment_evidence,
    Neoantimon2:::native_allo_row(r)$assignments$assignment_evidence)
  expect_false(grepl('/synthetic/assignment.tsv', Neoantimon2:::canonical_bytes(payload), fixed = TRUE))
})

test_that('assignment evidence semantic method provider version and result affect identity', {
  r <- native_assigned_row(); id <- Neoantimon2:::native_allo_id(Neoantimon2:::native_allo_row(r))
  for (k in c('method', 'provider', 'method_version', 'provider_version')) {
    q <- r; q$assignments$assignment_evidence[[1L]][[k]] <- paste0('changed-', k)
    expect_false(identical(native_enable(rows = list(q))$allo_native_evidence$records[[1L]]$evidence_id, id))
  }
  r$assignments$source_status <- 'ambiguous'
  id <- native_enable(rows = list(r))$allo_native_evidence$records[[1L]]$evidence_id
  for (state in c('not_supported', 'ambiguous', 'not_evaluated', 'failed')) {
    q <- r; q$assignments$assignment_evidence[[1L]]$result_state <- state
    q$assignments$assignment_evidence[[1L]]$evaluability <- if (state == 'not_supported') 'evaluated' else state
    expect_false(identical(native_enable(rows = list(q))$allo_native_evidence$records[[1L]]$evidence_id, id))
  }
})

test_that('assignment evidence changing source and presenter claims changes observation identity', {
  for (kind in c('source_origin', 'presenter')) {
    x <- phase_c_case('two_presenters'); r <- native_assigned_row(kind, x)
    graph <- x$public_explanation$contexts[[1L]]$graph
    if (kind == 'source_origin') {
      ids <- unique(vapply(Filter(function(n) n$data$semantics$peptide_key == r$analyte$peptide_key,
        graph$nodes$peptide_occurrence), function(n) n$data$semantics$source_occurrence_id, ''))
      field <- 'source_occurrence_ids'; claim <- 'source_occurrence_id'
    } else {
      ids <- unique(vapply(Filter(function(n) n$data$peptide_key == r$analyte$peptide_key,
        graph$nodes$pmhc), function(n) n$data$presenter, ''))
      field <- 'presenters'; claim <- 'presenter'
    }
    expect_gte(length(ids), 2L)
    q <- r; q$assignments[[field]] <- as.list(setdiff(ids, r$assignments[[field]][[1L]])[1L])
    q$assignments$assignment_evidence[[1L]][[claim]] <- q$assignments[[field]][[1L]]
    p <- native_enable(phase_c_project('two_presenters'), list(r), list(native_scope(x)))
    after <- native_enable(phase_c_project('two_presenters'), list(q), list(native_scope(x)))
    expect_false(identical(p$allo_native_evidence$records[[1L]]$evidence_id, after$allo_native_evidence$records[[1L]]$evidence_id))
    expect_identical(p$indirect_hla, after$indirect_hla)
    expect_identical(neo_explain(p)$indirect_hla_explanation$compatibility, neo_explain(after)$indirect_hla_explanation$compatibility)
  }
})
