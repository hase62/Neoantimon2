test_that('one complete source protein generates all stride-one 15-mers', {
  c <- indirect_context_test()
  expect_length(c$generation$occurrences, 30L)
  expect_length(unique(vapply(c$generation$occurrences,function(o)o$semantics$peptide_key,'')),6L)
  expect_true(c$generation$complete)
  expect_true(all(vapply(c$generation$occurrences,function(o)nchar(o$peptide_sequence)==15L,TRUE)))
})
test_that('generation includes both terminal windows in precursor coordinates', {
  c <- indirect_context_test(); o <- Filter(function(o)o$source_allele_hypothesis=='HLA-A*02:01',c$generation$occurrences)
  expect_setequal(vapply(o,function(o)o$semantics$start,1L),1:6)
  expect_setequal(vapply(o,function(o)o$semantics$end,1L),15:20)
  expect_true(all(vapply(o,function(o)o$sequence_coordinates=='precursor',TRUE)))
})
test_that('repeated ligand at different source coordinates retains occurrences', {
  m <- indirect_fixture('reference.json'); r <- m$records[[1L]]
  r$protein_sequence <- paste0('ACDEFGHIKLMNPQR','ACDEFGHIKLMNPQR'); r$domains[[1L]]$end<-30L;r$covered_regions[[1L]]$end<-30L
  g <- Neoantimon2:::indirect_sequence_windows(Neoantimon2:::allo_reference(list(manifest=within(m,records<-list(r)),path='fixture',sha256='physical'))$records[[1L]])
  matches <- Filter(function(w)w$peptide_sequence=='ACDEFGHIKLMNPQR',g$windows)
  expect_length(matches,2L); expect_setequal(vapply(matches,`[[`,1L,'start'),c(1L,16L))
  key <- Neoantimon2:::biological_id('peptide_key',list(ligand_sequence='ACDEFGHIKLMNPQR'),Neoantimon2:::peptide_contract())
  ids <- vapply(matches,function(w)Neoantimon2:::indirect_id('indirect-peptide-occurrence-id-v1',Neoantimon2:::indirect_occurrence_payload('source:1',w$start,w$end,key,'policy:1')),'')
  expect_length(unique(ids),2L)
})
test_that('same ligand from different source alleles keeps one peptide identity', {
  c<-indirect_context_test();o<-Filter(function(o)o$peptide_sequence=='ACDEFGHIKLMNPQR',c$generation$occurrences)
  expect_length(o,5L);expect_length(unique(vapply(o,function(o)o$semantics$peptide_key,'')),1L)
  expect_length(unique(vapply(o,`[[`,'','indirect_peptide_occurrence_id')),5L)
})
test_that('unknown-crossing windows are blocked and partial known windows unevaluable', {
  m<-indirect_fixture('reference.json');r<-m$records[[1L]];r$sequence_status<-'partial';r$protein_sequence<-paste0(substr(r$protein_sequence,1,18),'XY')
  r$covered_regions[[1L]]$end<-18L;m$records[[1L]]<-r
  c<-indirect_context_test(indirect_with_reference(m));expect_null(c$summary$complete_count)
  rr<-Neoantimon2:::allo_reference(list(manifest=m,path='fixture',sha256='hash'))$records
  g<-Neoantimon2:::indirect_sequence_windows(Filter(function(r)r$allele=='HLA-A*02:01',rr)[[1L]])
  expect_identical(g$status,'insufficient'); expect_length(g$blocked,2L); expect_length(g$windows,4L)
  expect_true(all(vapply(g$windows,function(w)w$status=='insufficient'&&!grepl('X',w$peptide_sequence),TRUE)))
})
test_that('unavailable source is not a zero generation result', {
  m<-indirect_fixture('reference.json');m$records[[1L]]$sequence_status<-'unavailable';m$records[[1L]]['protein_sequence']<-list(NULL)
  m$records[[1L]]$domains<-list();m$records[[1L]]$covered_regions<-list()
  c<-indirect_context_test(indirect_with_reference(m)); expect_false(c$generation$complete);expect_null(c$summary$complete_count)
  expect_true(any(vapply(c$generation$ledger,function(l)l$status=='unavailable',TRUE)))
})
test_that('reference typing and provider row order cannot change biological output', {
  a<-indirect_test(); x<-indirect_input();m<-indirect_fixture('reference.json');m$records<-rev(m$records)
  p<-indirect_fixture('binding.json');p$rows<-rev(p$rows);p$provider<-rev(p$provider)
  x$a$typing<-rev(x$a$typing);x$a$subjects<-rev(x$a$subjects);x$a$contexts[[1L]]<-rev(x$a$contexts[[1L]])
  b<-indirect_test(indirect_with_binding(p,indirect_with_reference(m,x)))
  expect_identical(a$result$allo_analysis_id,b$result$allo_analysis_id)
  expect_identical(a$result$contexts,b$result$contexts)
  expect_identical(a$foundation$allo_analysis_id,b$foundation$allo_analysis_id)
})
test_that('supported core self match is excluded despite different full ligand', {
  c<-indirect_context_test();s<-Filter(function(q)q$status=='self_match'&&q$semantics$query_type=='core_pmhc',c$comparisons)
  expect_length(s,1L); expect_gt(length(s[[1L]]$matches),0L)
  expect_identical(c$summary$excluded_self_count,1L)
  expect_true(any(vapply(c$comparisons,function(q)q$semantics$query_type=='full_peptide'&&q$status=='supported_nonself_within_scope',TRUE)))
})
test_that('supported non-self requires a completed same-presenter self search', {
  c<-indirect_context_test();s<-Filter(function(q)q$status=='supported_nonself_within_scope'&&q$semantics$query_type=='core_pmhc',c$comparisons)
  expect_length(s,5L);expect_true(all(vapply(s,`[[`,TRUE,'search_complete')))
  expect_true(all(vapply(s,function(q)!q$biological_self_completeness_established,TRUE)))
})
test_that('incomplete self catalog preserves positive hits but cannot establish absence', {
  x<-indirect_input();x$a$self_search[[1L]]$catalog_complete<-FALSE;c<-indirect_context_test(x)
  expect_null(c$summary$complete_count);expect_identical(c$summary$excluded_self_count,1L)
  expect_true(any(vapply(c$comparisons,function(q)q$status=='insufficient_search',TRUE)))
  expect_false(any(vapply(c$comparisons,function(q)q$status=='supported_nonself_within_scope',TRUE)))
})
test_that('ambiguous self typing is never selected as a genotype', {
  m<-indirect_fixture('reference.json');r<-Filter(function(r)r$allele=='HLA-A*01:01',m$records)[[1L]]
  r$allele<-'HLA-A*01:01:02';r$sequence_identifier<-'alternative';r$protein_sequence<-'YYYYYYYYYYYYYYY';m$records<-c(m$records,list(r))
  c<-indirect_context_test(indirect_with_reference(m));expect_null(c$summary$complete_count)
  expect_true(any(vapply(c$comparisons,function(q)q$status=='ambiguous',TRUE)))
})
test_that('missing self reference remains insufficient search', {
  m<-indirect_fixture('reference.json');m$records<-Filter(function(r)r$allele!='HLA-B*08:01',m$records)
  c<-indirect_context_test(indirect_with_reference(m));expect_null(c$summary$complete_count)
  expect_true(any(vapply(c$comparisons,function(q)q$status=='insufficient_search',TRUE)))
})
test_that('self comparison is molecular and merges multiple occurrence links', {
  c<-indirect_context_test();q<-Filter(function(q)q$semantics$query_type=='core_pmhc',c$comparisons)
  expect_length(q,6L);expect_true(all(vapply(q,function(q)length(q$occurrence_ids)==5L,TRUE)))
  expect_false(any(vapply(q,function(q)'occurrence_ids'%in%names(q$semantics),TRUE)))
})
test_that('physical path and hash changes cannot change comparison or pMHC identity', {
  a<-indirect_test();m<-indirect_fixture('reference.json');m$provenance$retrieved_at<-'2099-01-01'
  p<-indirect_fixture('binding.json');p$provenance$metadata<-list(path='/physical/alias');p$rows<-rev(p$rows)
  b<-indirect_test(indirect_with_binding(p,indirect_with_reference(m)))
  expect_false(identical(a$result$provenance,b$result$provenance));expect_identical(a$result$contexts,b$result$contexts)
  expect_identical(a$result$allo_analysis_id,b$result$allo_analysis_id)
  expect_identical(a$foundation$reference$reference_policy_id,b$foundation$reference$reference_policy_id)
})
test_that('same ligand same presenter has exactly one pMHC with multiple origins', {
  c<-indirect_context_test();expect_length(c$pmhcs,6L)
  expect_length(unique(vapply(c$pmhcs,`[[`,'','indirect_pmhc_id')),6L)
  expect_true(all(vapply(c$pmhcs,function(p)length(p$occurrence_ids)==5L,TRUE)))
})
test_that('same ligand on different DR presenters has distinct pMHCs', {
  c<-indirect_context_test(indirect_second_presenter());expect_identical(c$summary$complete_count,10L)
  expect_length(c$pmhcs,12L);expect_length(unique(vapply(c$pmhcs,`[[`,'','indirect_pmhc_id')),12L)
  expect_length(unique(vapply(c$pmhcs,`[[`,'','peptide_key')),6L)
})
test_that('overlapping ligands sharing a core are counted once per presenter', {
  p<-indirect_fixture('binding.json');p$rows[[1L]]$cores[[1L]]$offset<-1L;p$rows[[1L]]$cores[[1L]]$core_sequence<-'CDEFGHIKL'
  c<-indirect_context_test(indirect_with_binding(p));expect_identical(c$summary$complete_count,5L)
  expect_identical(c$summary$ligand_diagnostic_count,6L)
  expect_length(c$pmhcs,6L)
})
test_that('alternative registers are retained but not counted as primary', {
  p<-indirect_fixture('binding.json');alt<-p$rows[[2L]]$cores[[1L]];alt$offset<-1L;alt$core_sequence<-substr(p$rows[[2L]]$peptide_sequence,2L,10L);alt$primary<-FALSE;alt$register<-'alternative:1'
  p$rows[[2L]]$cores<-c(p$rows[[2L]]$cores,list(alt));c<-indirect_context_test(indirect_with_binding(p))
  expect_identical(c$summary$complete_count,5L)
  expect_true(any(vapply(c$pmhcs,function(p)length(p$core_assignments)==2L,TRUE)))
  p$rows[[2L]]$cores[[2L]]$primary<-TRUE; c<-indirect_context_test(indirect_with_binding(p))
  expect_null(c$summary$complete_count);expect_gt(c$summary$unevaluable_count,0L)
})
test_that('unsupported presenters remain explicit unevaluable strata', {
  x<-indirect_input();x$b$presenters[[1L]]$presenter<-'HLA-DPA1*01:03~HLA-DPB1*04:01';c<-indirect_context_test(x)
  expect_null(c$summary$complete_count);expect_identical(c$presenters[[1L]]$status,'unsupported')
  p<-indirect_fixture('binding.json');p$provider$supported_presenters<-list('HLA-DRA*01:01~HLA-DRB1*03:01');p$rows<-list()
  c<-indirect_context_test(indirect_with_binding(p));expect_null(c$summary$complete_count)
  expect_true(all(vapply(c$pmhcs,function(p)p$binding_evidence$status=='unsupported',TRUE)))
})
for (case in list(list(name='below',value=999.999999,expected=5L),list(name='equal',value=1000,expected=4L),list(name='above',value=1000.000001,expected=4L))) local({
  cc<-case
  test_that(paste('strict affinity boundary',cc$name,'1000 nM'), {
    p<-indirect_fixture('binding.json');p$rows[[2L]]$affinity_nM<-cc$value
    c<-indirect_context_test(indirect_with_binding(p));expect_identical(c$summary$complete_count,cc$expected)
    expect_identical(Neoantimon2:::indirect_threshold()$comparator,'<');expect_identical(Neoantimon2:::indirect_threshold()$unit,'nM')
  })
})
test_that('valid included units retain the full lower-level hierarchy', {
  c<-indirect_context_test();u<-Filter(function(u)u$state=='included',c$baseline_units)
  expect_length(u,5L);expect_true(all(vapply(u,function(u)!is.null(u$self_comparison_id)&&length(u$occurrence_ids)==5L,TRUE)))
  expect_length(c$source_presenter_summaries,5L)
})
test_that('self matches always exclude a strict baseline unit', {
  c<-indirect_context_test();u<-Filter(function(u)u$state=='excluded_self',c$baseline_units)
  expect_length(u,1L);expect_false(u[[1L]]$epitope_key%in%unlist(c$summary$member_ids))
})
test_that('provider failures are missingness not nonbinding', {
  p<-indirect_fixture('binding.json');p$rows[[2L]]$status<-'failed';p$rows[[2L]]['affinity_nM']<-list(NULL);p$rows[[2L]]$cores<-list();p$rows[[2L]]$reason_codes<-list('provider_execution_failed')
  c<-indirect_context_test(indirect_with_binding(p));expect_null(c$summary$complete_count);expect_identical(c$summary$unevaluable_count,1L)
  expect_identical(c$summary$excluded_nonbinder_count,0L)
})
test_that('failed self prediction cannot make source binders non-self', {
  p<-indirect_fixture('binding.json');p$rows<-Filter(function(r)r$peptide_sequence!='ACDEFGHIKAAAAAA',p$rows)
  c<-indirect_context_test(indirect_with_binding(p));expect_null(c$summary$complete_count)
  expect_identical(c$summary$conditional_included_count,0L);expect_identical(c$summary$unevaluable_count,6L)
})
test_that('duplicate biological occurrence and provider row inputs normalize once', {
  a<-indirect_context_test();x<-indirect_input();x$a$typing<-c(x$a$typing,x$a$typing[1L]);x$b$presenters<-rep(x$b$presenters,2L)
  p<-indirect_fixture('binding.json');p$rows<-c(p$rows,p$rows);b<-indirect_context_test(indirect_with_binding(p,x))
  expect_identical(a,b)
})
test_that('complete pair count is the core-presenter union not source sum', {
  c<-indirect_context_test();expect_identical(c$summary$complete_count,5L);expect_true(c$summary$analysis_complete)
  expect_identical(sum(vapply(c$source_presenter_summaries,`[[`,1L,'complete_count')),25L)
  expect_true(c$summary$source_attributed_counts_may_overlap)
})
test_that('incomplete pair complete count is null even when conditional members exist', {
  x<-indirect_input();x$a$typing<-Filter(function(t)!(t$subject_id=='donor-1'&&t$locus=='C'),x$a$typing)
  c<-indirect_context_test(x);expect_null(c$summary$complete_count);expect_false(c$summary$analysis_complete)
  expect_identical(c$summary$conditional_included_count,5L);expect_true(c$summary$conditional_count_is_not_a_lower_bound)
})
test_that('profile claim remains literature-derived and never official reproduction', {
  x<-indirect_test()$result;expect_identical(x$profile,'neo-literature-derived-indirect-hla-ii-dr-strict-v1')
  expect_match(x$claim,'not official PIRCHE reproduction');expect_false(any(c('PIRCHE_score','pirche_score','rejection_probability')%in%names(x)))
  expect_identical(x$provider$semantics$evidence_kind,'synthetic_fixture')
})
test_that('Phase B frozen namespace payloads are exact and exclude record metadata', {
  c<-indirect_context_test()
  for(o in c$generation$occurrences){
    expected<-list(source_occurrence_id=o$semantics$source_occurrence_id,start=as.integer(o$semantics$start),end=as.integer(o$semantics$end),peptide_key=o$semantics$peptide_key,generation_policy_id=o$semantics$generation_policy_id)
    expect_identical(o$semantics,expected)
    expect_identical(o$indirect_peptide_occurrence_id,Neoantimon2:::indirect_id('indirect-peptide-occurrence-id-v1',expected))
  }
  for(p in c$pmhcs){
    expected<-list(allo_context_id=c$allo_context_id,candidate_key=p$candidate_key,presenting_compartment='declared_test_apc',presenting_scenario='sot_recipient_response')
    expect_identical(p$semantics,expected);expect_identical(p$indirect_pmhc_id,Neoantimon2:::indirect_id('indirect-pmhc-id-v1',expected))
  }
  for(q in c$comparisons){
    expected<-list(query_identity=q$semantics$query_identity,query_type=q$semantics$query_type,self_universe_id=c$self_universe_id,
      algorithm=list(version='exact-literal-same-presenter-primary-core-1.0.0',matching='literal_exact',wildcard=FALSE,il_equivalence=FALSE,comparison_level=q$semantics$query_type,presenter=if(q$semantics$query_type=='core_pmhc')'HLA-DRA*01:01~HLA-DRB1*01:01' else NULL),
      self_presentation_policy=if(q$semantics$query_type=='core_pmhc')list(provider_configuration_id=indirect_test()$result$provider$provider_configuration_id,threshold=Neoantimon2:::indirect_threshold(),registers='designated_primary_forward_only') else NULL)
    expect_identical(q$semantics,expected);expect_identical(q$self_comparison_id,Neoantimon2:::indirect_id('self-comparison-id-v1',expected))
  }
})
test_that('all frozen Phase A namespaces and state remain exact under Phase B', {
  a<-indirect_test(); b<-Neoantimon2:::allo_foundation(Neoantimon2:::allo_config(indirect_input()$a))
  expect_identical(a$foundation,b);expect_identical(a$foundation$analysis_status,'not_evaluated')
  expect_identical(Neoantimon2:::allo_context_identity_payload(b$contexts[[1L]]$semantics),b$contexts[[1L]]$semantics[setdiff(names(b$contexts[[1L]]$semantics),'claim_level')])
  expect_setequal(getNamespaceExports('Neoantimon2'), current_public_api)
})
test_that('pinned import rejects contradictory rows and unsupported engines without inference', {
  m<-indirect_fixture('binding.json');m$provider$name<-'Frost'
  expect_error(indirect_test(indirect_with_binding(m)),'requires pinned NetMHCIIpan')
  m<-indirect_fixture('binding.json');r<-m$rows[[1L]];r$affinity_nM<-600;m$rows<-c(m$rows,list(r))
  expect_error(indirect_test(indirect_with_binding(m)),'Contradictory duplicate')
  m<-indirect_fixture('binding.json');m$rows[[1L]]$cores[[1L]]$core_sequence<-'YYYYYYYYY'
  expect_error(indirect_test(indirect_with_binding(m)),'does not map')
})
test_that('unsupported role reference and self scope cannot silently execute', {
  x<-indirect_input();x$a$self_search[[1L]]$scope<-'corresponding_hla';x$a$self_search[[1L]]$required_loci<-list('A')
  expect_error(indirect_test(x),'requires all-self-HLA')
  x<-indirect_input();x$b$presenters[[1L]]$resolution_policy<-'assume_DRA'
  expect_error(indirect_test(x),'DRA resolution must be explicit')
})
test_that('provider pinning and tampered runtime state fail closed', {
  p<-indirect_enable(reconstruction_test_project());expect_identical(Neoantimon2:::indirect_read(p),p$indirect_hla)
  p$indirect_hla$contexts[[1L]]$summary$complete_count<-0L
  expect_error(Neoantimon2:::indirect_read(p),'derived contract mismatch')
  x<-indirect_input();x$b$provider$sha256<-paste(rep('0',64),collapse='')
  expect_error(indirect_test(x),'integrity mismatch')
})
test_that('public constructor and compact accessors preserve frozen somatic identities', {
  p<-reconstruction_test_project();cfg<-p$config;x<-indirect_input();cfg$allo<-x$a;cfg$indirect_hla<-x$b
  q<-neo_project(samples=p$samples,hla=if(is.null(p$hla))NULL else lapply(p$hla,unlist),config=cfg,analysis_context=p$analysis_context,inputs=recognition_clone_inputs(p),fixture_kind=p$fixture_kind)
  expect_identical(p$analysis_id,q$analysis_id);expect_identical(p$config_fingerprint,q$config_fingerprint);expect_identical(p$input_fingerprint,q$input_fingerprint)
  expect_identical(neo_explain(q)$allo_foundation,q$allo_foundation)
  expect_identical(neo_explain(q)$indirect_hla_baseline$profile,Neoantimon2:::indirect_profile())
  expect_match(neo_explain(q,format='text'),'not official PIRCHE reproduction')
})
test_that('same core text on a different presenter is a different self-pMHC', {
  x<-indirect_second_presenter();m<-jsonlite::read_json(x$b$provider$path)
  for(i in seq_along(m$rows))if(m$rows[[i]]$peptide_sequence=='ACDEFGHIKAAAAAA'&&grepl('DRB1\\*03:01',m$rows[[i]]$presenter))m$rows[[i]]$affinity_nM<-1000
  c<-indirect_context_test(indirect_with_binding(m,x))
  expect_identical(c$summary$complete_count,11L);expect_identical(c$summary$excluded_self_count,1L)
  q<-Filter(function(q)q$semantics$query_type=='core_pmhc'&&grepl('DRB1\\*03:01',q$semantics$algorithm$presenter),c$comparisons)
  expect_true(all(vapply(q,function(q)q$status=='supported_nonself_within_scope',TRUE)))
})
test_that('provider aliases and metadata do not redefine evidence or molecular identities', {
  a<-indirect_test()$result;m<-indirect_fixture('binding.json')
  for(i in seq_along(m$rows)){
    m$rows[[i]]$raw_output$path<-paste0('/physical/row/',i)
    m$rows[[i]]$raw_output$timestamp<-'2099-01-01'
    m$rows[[i]]$provenance$sha256<-'physical-checksum'
  }
  b<-indirect_test(indirect_with_binding(m))$result
  expect_identical(a$allo_analysis_id,b$allo_analysis_id)
  expect_identical(lapply(a$provider$rows,`[[`,'evidence_id'),lapply(b$provider$rows,`[[`,'evidence_id'))
  expect_identical(a$contexts[[1L]]$summary,b$contexts[[1L]]$summary)
  expect_identical(lapply(a$contexts[[1L]]$pmhcs,`[[`,'indirect_pmhc_id'),lapply(b$contexts[[1L]]$pmhcs,`[[`,'indirect_pmhc_id'))
  m$provider$model<-'another explicit pinned model';c<-indirect_test(indirect_with_binding(m))$result
  expect_false(identical(a$allo_analysis_id,c$allo_analysis_id))
  expect_identical(lapply(a$contexts[[1L]]$pmhcs,`[[`,'indirect_pmhc_id'),lapply(c$contexts[[1L]]$pmhcs,`[[`,'indirect_pmhc_id'))
  expect_identical(a$contexts[[1L]]$generation,c$contexts[[1L]]$generation)
})
test_that('short proteins have explicit not-applicable windows and complete zero is allowed only after evaluation', {
  r<-indirect_fixture('reference.json')$records[[1L]];r$protein_sequence<-'ACDEFGHIK';r$domains[[1L]]$end<-9L;r$covered_regions[[1L]]$end<-9L
  g<-Neoantimon2:::indirect_sequence_windows(r);expect_identical(g$status,'not_applicable');expect_length(g$windows,0L)
  m<-indirect_fixture('binding.json')
  for(i in seq_along(m$rows))if(m$rows[[i]]$peptide_sequence!='ACDEFGHIKAAAAAA')m$rows[[i]]$affinity_nM<-1000
  c<-indirect_context_test(indirect_with_binding(m));expect_identical(c$summary$complete_count,0L);expect_true(c$summary$analysis_complete)
})
test_that('omitted self envelope retains explicit subject without fabricating a universe or comparison identity', {
  x<-indirect_input();x$a$self_search<-list();c<-indirect_context_test(x)
  expect_null(c$summary$complete_count);expect_length(c$comparisons,6L)
  expect_true(all(vapply(c$comparisons,function(q)q$self_subject_id=='patient-1'&&q$status=='not_evaluated',TRUE)))
  expect_true(all(vapply(c$comparisons,function(q)is.null(q$self_comparison_id)&&is.null(q$semantics$self_universe_id),TRUE)))
})
test_that('frozen presenter normalization preserves chain-order and spelling equivalence', {
  a<-indirect_test()$result;x<-indirect_input();x$b$presenters[[1L]]$presenter<-'drb1*01:01 ~ dra*01:01'
  b<-indirect_test(x)
  expect_identical(a$allo_analysis_id,b$result$allo_analysis_id)
  expect_identical(a$contexts,b$result$contexts)
  expect_identical(b$config$presenters[[1L]]$presenter,'drb1*01:01 ~ dra*01:01')
  x$b$presenters[[1L]]$scenario_id<-'unknown:scenario';expect_error(indirect_test(x),'Unknown Phase B presenter scenario')
})
test_that('HLA-identical source and self retains presentation and self hits with complete zero', {
  x<-indirect_input()
  for(i in seq_along(x$a$typing))if(x$a$typing[[i]]$subject_id=='donor-1'){
    locus<-x$a$typing[[i]]$locus
    own<-Filter(function(t)t$subject_id=='patient-1'&&t$locus==locus,x$a$typing)[[1L]]
    x$a$typing[[i]]$alleles<-own$alleles
  }
  c<-indirect_context_test(x)
  expect_identical(c$summary$complete_count,0L);expect_true(c$summary$analysis_complete)
  expect_identical(c$summary$excluded_self_count,1L);expect_length(c$pmhcs,1L)
  expect_true(all(vapply(c$comparisons,function(q)q$status=='self_match',TRUE)))
})
test_that('Phase B result metadata stays compact and identifies synthetic evidence', {
  s<-Neoantimon2:::indirect_summary(indirect_test()$result)
  expect_identical(s$evidence_kind,'synthetic_fixture')
  expect_false(any(c('coverage','member_ids')%in%names(s$summaries[[1L]])))
  expect_identical(s$summaries[[1L]]$complete_count,5L)
})
test_that('self catalog never admits untyped or out-of-scope aliases of identical content', {
  m<-indirect_fixture('reference.json');r<-Filter(function(r)r$allele=='HLA-DRA*01:01',m$records)[[1L]]
  r$allele<-'HLA-DRA*02:01';r$sequence_identifier<-'untyped-partial-alias';r$sequence_status<-'partial'
  m$records<-c(m$records,list(r));c<-indirect_context_test(indirect_with_reference(m))
  expect_identical(c$summary$complete_count,5L)
  hits<-c$self_searches[[1L]]$hits
  expect_false(any(vapply(unlist(lapply(hits,`[[`,'matches'),recursive=FALSE),function(m)startsWith(m$allele,'HLA-DRA*'),TRUE)))
  expect_length(hits[[1L]]$matches,5L)
})
test_that('import reuses frozen prediction task vocabulary and separates BA from EL', {
  p<-indirect_test()$result$provider$semantics
  expect_true(p$task%in%Neoantimon2:::peptide_prediction_tasks());expect_identical(p$task,'mhc2_binding')
  expect_identical(p$prediction_mode,'BA')
  m<-indirect_fixture('binding.json');m$provider$prediction_mode<-'EL'
  expect_error(indirect_test(indirect_with_binding(m)),'requires pinned NetMHCIIpan')
})
test_that('semantic provider arguments are explicit BA-only and fail closed', {
  m <- indirect_fixture('binding.json')
  expect_identical(Neoantimon2:::indirect_provider_arguments(m$provider$arguments),list(BA=TRUE))
  expect_identical(indirect_test()$result$provider$semantics$arguments,list(BA=TRUE))
  for (value in list(FALSE,'true',1,NULL)) {
    bad <- m; bad$provider$arguments['BA'] <- list(value)
    expect_error(Neoantimon2:::indirect_provider(bad))
    expect_false(indirect_schema_valid(bad))
  }
  bad <- m; bad$provider$arguments <- list()
  expect_error(Neoantimon2:::indirect_provider(bad),'semantic provider arguments')
  expect_false(indirect_schema_valid(bad))
  # No alternate BA value is applicable within the frozen strict profile.
  # Existing allowed model/configuration changes remain semantic identity changes.
  a <- Neoantimon2:::indirect_provider(m); m$provider$model <- 'different semantic model configuration'
  b <- Neoantimon2:::indirect_provider(m)
  expect_false(identical(a$provider_configuration_id,b$provider_configuration_id))
  expect_false(identical(lapply(a$rows,`[[`,'evidence_id'),lapply(b$rows,`[[`,'evidence_id')))
})
test_that('unknown and physical provider arguments never enter semantic identity', {
  for (key in c('input_path','output_path','executable_path','model_path','sha256','checksum',
    'timestamp','run_time','temporary_directory','query_file','result_file','unknown_option')) {
    m <- indirect_fixture('binding.json'); m$provider$arguments[[key]] <- 'physical'
    expect_error(Neoantimon2:::indirect_provider(m),'semantic provider arguments')
    expect_false(indirect_schema_valid(m))
    expect_false(indirect_schema_valid(m,'indirect-hla-baseline','import'))
  }
})
test_that('physical pinned provider artifacts cannot redefine any policy or analysis identity', {
  m <- indirect_precomputed_contract(); a <- indirect_test(indirect_with_binding(m))$result
  for (key in c('query_artifact','result_artifact','executable_artifact','model_artifact')) {
    m$provenance[[key]]$reference <- paste0('/different/physical/',key)
    m$provenance[[key]]$sha256 <- paste(rep('f',64L),collapse='')
  }
  m$provenance$execution$command <- 'different physical invocation with same BA semantics'
  m$provenance$metadata <- list(timestamp='2099-01-01',run_time=4,temporary_directory='/tmp/elsewhere',
    input_path='/physical/query',output_path='/physical/result')
  for (i in seq_along(m$rows)) {
    m$rows[[i]]$raw_output$path <- '/different/raw/result'
    m$rows[[i]]$provenance <- list(timestamp='2099-01-01',sha256=paste(rep('e',64),collapse=''))
  }
  m$rows <- rev(m$rows)
  b <- indirect_test(indirect_with_binding(m))$result
  expect_true(indirect_schema_valid(m))
  expect_identical(a$provider$provider_configuration_id,b$provider$provider_configuration_id)
  expect_identical(lapply(a$provider$rows,`[[`,'evidence_id'),lapply(b$provider$rows,`[[`,'evidence_id'))
  expect_identical(a$allo_policy_id,b$allo_policy_id)
  expect_identical(a$allo_analysis_id,b$allo_analysis_id)
  expect_identical(a$summaries,b$summaries)
  expect_false(identical(a$provider$provenance,b$provider$provenance))
  expect_false(identical(a$provenance,b$provenance))
})
test_that('precomputed evidence requires structured SHA-256 artifact and execution pinning', {
  m <- indirect_precomputed_contract()
  expect_true(indirect_schema_valid(m))
  expect_identical(Neoantimon2:::indirect_provider(m)$provenance,m$provenance)
  expect_identical(indirect_test(indirect_with_binding(m))$result$summaries[[1L]]$complete_count,5L)
  for (key in names(m$provenance)[names(m$provenance)!='metadata']) {
    bad <- m; bad$provenance[[key]] <- NULL
    expect_error(Neoantimon2:::indirect_provider(bad))
    expect_false(indirect_schema_valid(bad))
  }
  for (key in c('query_artifact','result_artifact','executable_artifact','model_artifact')) {
    for (value in c('MD5:abcd',paste(rep('g',64),collapse=''),paste(rep('a',63),collapse=''))) {
      bad <- m; bad$provenance[[key]]$sha256 <- value
      expect_error(Neoantimon2:::indirect_provider(bad),'SHA-256')
      expect_false(indirect_schema_valid(bad))
    }
    for (field in c('identity','reference','sha256')) {
      bad <- m; bad$provenance[[key]][[field]] <- NULL
      expect_error(Neoantimon2:::indirect_provider(bad))
      expect_false(indirect_schema_valid(bad))
    }
  }
  bad <- m; bad$provenance <- NULL
  expect_error(Neoantimon2:::indirect_provider(bad),'provenance'); expect_false(indirect_schema_valid(bad))
  bad <- m; bad$provenance <- list(source='unstructured claim')
  expect_error(Neoantimon2:::indirect_provider(bad)); expect_false(indirect_schema_valid(bad))
  for (key in c('executable_artifact','model_artifact')) {
    bad <- m; bad$provenance[[key]]$identity <- 'unmatched label'
    expect_error(Neoantimon2:::indirect_provider(bad),'match pinned build/model')
  }
  bad <- m; bad$provenance$execution$provider_version <- '4.2'
  expect_error(Neoantimon2:::indirect_provider(bad),'conflicts'); expect_false(indirect_schema_valid(bad))
  bad <- m; bad$provenance$license$status <- 'unknown'
  expect_error(Neoantimon2:::indirect_provider(bad),'license/provision'); expect_false(indirect_schema_valid(bad))
})
test_that('synthetic provenance explicitly lacks historical execution and model artifacts', {
  m <- indirect_fixture('binding.json')
  expect_true(indirect_schema_valid(m)); expect_identical(Neoantimon2:::indirect_provider(m)$provenance$execution,'not_executed')
  for (key in c('executable_artifact','model_artifact')) {
    bad <- m; bad$provenance[[key]] <- indirect_precomputed_contract()$provenance[[key]]
    expect_error(Neoantimon2:::indirect_provider(bad),'synthetic provider provenance')
    expect_false(indirect_schema_valid(bad))
  }
  bad <- m; bad$provenance$execution <- 'historical execution'
  expect_error(Neoantimon2:::indirect_provider(bad),'not_executed'); expect_false(indirect_schema_valid(bad))
  bad <- m; bad$provider$evidence_kind <- 'precomputed_predictor_output'
  expect_error(Neoantimon2:::indirect_provider(bad),'precomputed provider provenance'); expect_false(indirect_schema_valid(bad))
})
test_that('recipient strict DRB1 scope requires every exactly evaluated DRB1 presenter', {
  x <- indirect_second_presenter(); a <- indirect_context_test(x)
  expect_identical(a$summary$complete_count,10L); expect_true(a$presenter_scope$complete)
  expect_length(a$presenter_scope$evaluated_beta_alleles,2L)
  x$b$presenters <- x$b$presenters[1L]; b <- indirect_context_test(x)
  expect_null(b$summary$complete_count); expect_false(b$summary$analysis_complete)
  expect_identical(b$summary$conditional_included_count,5L)
  expect_true(b$summary$conditional_count_is_not_a_lower_bound)
  expect_true('presenter_scope_incomplete'%in%unlist(b$summary$reason_codes))
  expect_identical(b$presenter_scope$undeclared_beta_alleles,list('HLA-DRB1*03:01'))
  expect_true(all(vapply(b$source_presenter_summaries,function(s)is.null(s$complete_count),TRUE)))
  expect_length(b$pmhcs,6L)
  x$a$typing <- rev(x$a$typing); expect_identical(b,indirect_context_test(x))
  x <- indirect_second_presenter(); x$a$typing <- rev(x$a$typing); x$b$presenters <- rev(x$b$presenters)
  expect_identical(a,indirect_context_test(x))
})
test_that('presenter scope uses only the selected Phase A self typing scenario', {
  x <- indirect_second_presenter(); x$b$presenters <- x$b$presenters[1L]
  for (i in seq_along(x$a$typing)) if (identical(x$a$typing[[i]]$alleles,list('DRB1*03:01')))
    x$a$typing[[i]]$typing_scenario_id <- 'other:typing'
  c <- indirect_context_test(x)
  expect_true(c$presenter_scope$complete); expect_identical(c$summary$complete_count,5L)
  expect_identical(c$presenter_scope$typing_scenario_id,'typing:1')
  expect_identical(c$presenter_scope$evaluated_beta_alleles,list('HLA-DRB1*01:01'))
  x$a$self_search[[1L]]$typing_scenario_id <- 'other:typing'
  c <- indirect_context_test(x); expect_null(c$summary$complete_count)
  expect_true('presenter_scope_incomplete'%in%unlist(c$summary$reason_codes))
})
test_that('unresolved recipient DR typing cannot produce a complete zero or inferred pair', {
  for (kind in c('missing_alpha','incomplete_alpha','ambiguous_beta')) {
    x <- indirect_second_presenter(); x$b$presenters <- x$b$presenters[1L]
    if (kind=='missing_alpha') x$a$typing <- Filter(function(t)t$locus!='DRA',x$a$typing)
    if (kind=='incomplete_alpha') for (i in seq_along(x$a$typing)) if (x$a$typing[[i]]$locus=='DRA') x$a$typing[[i]]$complete <- FALSE
    if (kind=='ambiguous_beta') for (i in seq_along(x$a$typing)) if (identical(x$a$typing[[i]]$alleles,list('DRB1*03:01')))
      x$a$typing[[i]]$alleles <- list('DRB1*01:01','DRB1*03:01')
    m <- jsonlite::read_json(x$b$provider$path)
    for (i in seq_along(m$rows)) m$rows[[i]]$affinity_nM <- 1000
    c <- indirect_context_test(indirect_with_binding(m,x))
    expect_null(c$summary$complete_count); expect_false(c$presenter_scope$complete)
    expect_identical(c$summary$conditional_included_count,0L)
    expect_length(c$presenters,1L); expect_length(c$pmhcs,6L)
    expect_true('presenter_scope_incomplete'%in%unlist(c$summary$reason_codes))
  }
})
test_that('multiple resolved alpha chains fail closed without pairing inference', {
  x <- indirect_second_presenter(); m <- jsonlite::read_json(x$a$reference$path)
  r <- Filter(function(r)r$allele=='HLA-DRA*01:01',m$records)[[1L]]
  r$allele <- 'HLA-DRA*02:01'; r$sequence_identifier <- 'synthetic:second-alpha'; m$records <- c(m$records,list(r))
  t <- Filter(function(t)t$locus=='DRA',x$a$typing)[[1L]]; t$alleles <- list('DRA*02:01')
  x$a$typing <- c(x$a$typing,list(t)); x <- indirect_with_reference(m,x)
  # Both beta chains were declared with the first alpha. Do not invent pairs.
  c <- indirect_context_test(x)
  expect_null(c$summary$complete_count); expect_identical(c$summary$conditional_included_count,10L)
  expect_true('presenter_scope_pairing_ambiguous'%in%unlist(c$summary$reason_codes))
  expect_length(c$presenters,2L); expect_length(c$pmhcs,12L)
  expect_false(any(vapply(c$pmhcs,function(p)grepl('DRA\\*02:01',p$presenter),TRUE)))
  x$a$typing <- rev(x$a$typing); x$b$presenters <- rev(x$b$presenters)
  expect_identical(c,indirect_context_test(x))
})
test_that('core mismatch annotations do not widen frozen assignment identity payloads', {
  m <- indirect_fixture('binding.json'); alt <- m$rows[[1L]]$cores[[1L]]
  alt$primary <- FALSE; alt$orientation <- 'reverse'; alt$status <- 'unresolved_mapping'; alt$register <- 'reverse:fixture'
  m$rows[[1L]]$cores <- c(m$rows[[1L]]$cores,list(alt))
  c <- indirect_context_test(indirect_with_binding(m))
  expect_identical(c$summary$complete_count,5L)
  for (p in c$pmhcs) for (a in p$core_assignments) {
    expect_identical(a$mismatch_annotation,list(status='not_evaluated',reason_codes=list('alignment_not_performed')))
    expected <- list(indirect_pmhc_id=p$indirect_pmhc_id,core_sequence=a$semantics$core_sequence,
      offset=a$semantics$offset,register=a$semantics$register,orientation=a$semantics$orientation,
      provider_configuration_id=p$provider_configuration_id,
      assignment_semantics=a$semantics$assignment_semantics)
    expect_identical(a$semantics,expected)
    expect_identical(a$indirect_core_assignment_id,Neoantimon2:::indirect_id('indirect-core-assignment-id-v1',expected))
  }
  p <- c$pmhcs[[1L]]; expect_true(indirect_schema_valid(p,'indirect-pmhc-evidence'))
  bad <- p; bad$core_assignments[[1L]]$mismatch_annotation <- NULL
  expect_false(indirect_schema_valid(bad,'indirect-pmhc-evidence'))
  bad <- p; bad$core_assignments[[1L]]$mismatch_annotation$status <- 'overlap'
  expect_false(indirect_schema_valid(bad,'indirect-pmhc-evidence'))
})
test_that('nonempty provider scope is enforced by runtime and every raw/bound schema', {
  for (key in c('supported_presenters','supported_lengths')) {
    m <- indirect_fixture('binding.json'); m$provider[[key]] <- list()
    expect_error(Neoantimon2:::indirect_provider(m))
    expect_false(indirect_schema_valid(m))
    expect_false(indirect_schema_valid(m,'indirect-hla-baseline','import'))
    x <- indirect_test(); x$config$provider$manifest <- m
    expect_false(indirect_schema_valid(x$config,'indirect-hla-baseline','config'))
    x$result$provider$semantics[[key]] <- list()
    expect_false(indirect_schema_valid(x$result,'indirect-hla-baseline'))
    # Non-emptiness does not add a raw uniqueness restriction.
    m <- indirect_fixture('binding.json'); m$provider[[key]] <- rep(m$provider[[key]],2L)
    expect_true(indirect_schema_valid(m))
    expect_true(indirect_schema_valid(m,'indirect-hla-baseline','import'))
    expect_identical(Neoantimon2:::indirect_provider(m)$semantics[[key]],indirect_test()$result$provider$semantics[[key]])
  }
})
test_that('provider evidence policy and analysis namespace payloads are exact', {
  x <- indirect_test(); b <- x$result; f <- x$foundation; p <- b$provider
  provider_policy <- list(profile='neo-literature-derived-indirect-hla-ii-dr-strict-v1',provider=p$semantics)
  expect_identical(p$provider_configuration_id,Neoantimon2:::allo_id('policy',provider_policy))
  for (row in p$rows) {
    result <- row[c('peptide_sequence','presenter','status','affinity_nM','percentile_rank','cores','reason_codes')]
    payload <- list(target=result[c('peptide_sequence','presenter')],dimension='binding_affinity_and_core',
      provider_configuration_id=p$provider_configuration_id,assay_occasion_id=NULL,result=result)
    expect_identical(row$evidence_id,Neoantimon2:::indirect_id('indirect-evidence-id-v1',payload))
  }
  policy <- list(profile=Neoantimon2:::indirect_profile(),version='1.0.0',generation=Neoantimon2:::indirect_generation_policy(),
    reference_policy_id=f$reference$reference_policy_id,presenters=Neoantimon2:::allo_records(x$config$presenters),
    provider_configuration_id=p$provider_configuration_id,threshold=Neoantimon2:::indirect_threshold(),
    self='primary_forward_core_same_presenter_from_qualifying_self_15mers',
    counting='unique_core_presenter_union',uncertainty='fail_closed_no_imputation')
  expect_identical(b$policy,policy); expect_identical(b$allo_policy_id,Neoantimon2:::allo_id('policy',policy))
  payload <- list(context_ids=Neoantimon2:::allo_set(lapply(b$contexts,`[[`,'allo_context_id')),
    sequence_ids=Neoantimon2:::allo_set(lapply(f$reference$records,`[[`,'hla_source_sequence_id')),
    universe_ids=Neoantimon2:::allo_set(lapply(f$self_universes,`[[`,'self_universe_id')),
    policy_id=b$allo_policy_id,provider_configuration_ids=list(p$provider_configuration_id),
    input_evidence_ids=Neoantimon2:::allo_set(c(list(f$allo_analysis_id),lapply(p$rows,`[[`,'evidence_id'))))
  expect_identical(b$allo_analysis_id,Neoantimon2:::allo_id('analysis',payload))
})
test_that('exact accessory DR typing preserves Phase A without enlarging strict DRB1 scope', {
  for (locus in c('DRB3','DRB4','DRB5')) {
    input <- indirect_accessory_typing(locus); x <- indirect_test(input); c <- x$result$contexts[[1L]]
    base_input <- input; base_input$a$typing <- Filter(function(t)t$locus!=locus,base_input$a$typing)
    base <- indirect_context_test(base_input)
    expect_identical(c,base)
    expect_true(c$presenter_scope$complete); expect_identical(c$summary$complete_count,10L)
    expect_identical(c$presenter_scope$evaluated_beta_alleles,list('HLA-DRB1*01:01','HLA-DRB1*03:01'))
    expect_length(c$presenter_scope$undeclared_beta_alleles,0L)
    expect_true(all(vapply(c$presenters,function(p)grepl('~HLA-DRB1\\*',p$presenter),TRUE)))
    expect_true(all(vapply(c$pmhcs,function(p)grepl('~HLA-DRB1\\*',p$presenter),TRUE)))
    expect_true(all(vapply(c$baseline_units,function(u)grepl('~HLA-DRB1\\*',u$presenter),TRUE)))
    phase_a <- Neoantimon2:::allo_foundation(Neoantimon2:::allo_config(input$a))
    expect_identical(x$foundation,phase_a)
    accessory <- Filter(function(t)t$subject_id=='patient-1'&&t$locus==locus,x$foundation$typing)
    expect_length(accessory,1L); expect_identical(accessory[[1L]]$typing_scenario_id,'typing:1')
    expect_true(accessory[[1L]]$complete); expect_identical(accessory[[1L]]$evaluability$status,'evaluated')
    expect_identical(accessory[[1L]]$resolved_allele,paste0('HLA-',locus,'*01:01'))
    expect_identical(accessory[[1L]]$raw_alleles,list(paste0(locus,'*01:01')))
  }
})
test_that('explicit accessory DR presenters remain unsupported even with qualifying provider evidence', {
  for (locus in c('DRB3','DRB4','DRB5')) {
    input <- indirect_accessory_typing(locus)
    pair <- paste0('HLA-DRA*01:01~HLA-',locus,'*01:01')
    pp <- input$b$presenters[[1L]]; pp$presenter <- pair
    input$b$presenters <- c(input$b$presenters,list(pp))
    m <- jsonlite::read_json(input$b$provider$path)
    m$provider$supported_presenters <- c(m$provider$supported_presenters,list(pair))
    rows <- Filter(function(r)r$presenter==input$b$presenters[[1L]]$presenter,m$rows)
    m$rows <- c(m$rows,lapply(rows,function(r){r$presenter<-pair;r}))
    c <- indirect_context_test(indirect_with_binding(m,input))
    accessory <- Filter(function(p)p$presenter==pair,c$presenters)
    expect_length(accessory,1L); expect_identical(accessory[[1L]]$presenter,pair)
    expect_identical(accessory[[1L]]$status,'unsupported')
    expect_identical(accessory[[1L]]$reason_codes,list('presenter_out_of_profile_scope'))
    expect_true(c$presenter_scope$complete); expect_null(c$summary$complete_count)
    expect_false(c$summary$analysis_complete); expect_identical(c$summary$conditional_included_count,10L)
    expect_true('presenter_out_of_profile_scope'%in%unlist(c$summary$reason_codes))
    units <- Filter(function(u)u$presenter==pair,c$baseline_units)
    expect_length(units,6L); expect_true(all(vapply(units,function(u)u$state=='unevaluable',TRUE)))
    pmhcs <- Filter(function(p)p$presenter==pair,c$pmhcs)
    expect_length(pmhcs,6L); expect_true(all(vapply(pmhcs,function(p)p$status=='unsupported',TRUE)))
    expect_true(all(vapply(Filter(function(u)u$state=='included',c$baseline_units),
      function(u)grepl('~HLA-DRB1\\*',u$presenter),TRUE)))
    expect_identical(c$summary$member_ids,indirect_context_test(indirect_second_presenter())$summary$member_ids)
  }
})
test_that('DRA DRB1 and accessory typing order preserves strict completeness and counting', {
  for (locus in c('DRB3','DRB4','DRB5')) for (complete in c(TRUE,FALSE)) {
    x <- indirect_accessory_typing(locus)
    if (!complete) x$b$presenters <- x$b$presenters[1L]
    before <- indirect_test(x)
    x$a$typing <- rev(x$a$typing); x$b$presenters <- rev(x$b$presenters)
    after <- indirect_test(x)
    expect_identical(before$foundation,after$foundation)
    expect_identical(before$result$contexts,after$result$contexts)
    expect_identical(before$result$allo_analysis_id,after$result$allo_analysis_id)
    scope <- after$result$contexts[[1L]]$presenter_scope
    expect_identical(scope$complete,complete)
    expect_identical(scope$undeclared_beta_alleles,if(complete)list()else list('HLA-DRB1*03:01'))
    expect_identical(after$result$contexts[[1L]]$summary$complete_count,if(complete)10L else NULL)
  }
})
test_that('Phase A reference policy and analysis payload guards remain exact', {
  f <- indirect_test()$foundation
  reference <- f$reference$semantics
  reference$records <- Neoantimon2:::allo_records(lapply(f$reference$records,function(r)r[setdiff(names(r),'resolution_provenance')]))
  policy <- list(profile='neo-allo-phase-a-reference-v1',version='1.0.0',reference=reference,
    typing_imputation='none',resolution='all_compatible_pinned_sequences_v1')
  expect_identical(f$reference$reference_policy_id,Neoantimon2:::allo_id('policy',policy))
  payload <- list(phase='A',context_ids=Neoantimon2:::allo_set(lapply(f$contexts,`[[`,'allo_context_id')),
    sequence_ids=Neoantimon2:::allo_set(lapply(f$reference$records,`[[`,'hla_source_sequence_id')),
    universe_ids=Neoantimon2:::allo_set(lapply(f$self_universes,`[[`,'self_universe_id')),
    reference_policy_id=f$reference$reference_policy_id,
    self_search_evidence=Neoantimon2:::allo_records(lapply(f$self_universes,function(u)list(
      self_universe_id=u$self_universe_id,typing_scenario_id=u$typing_scenario_id,
      catalog_complete=u$catalog_complete,evaluability=u$evaluability))),
    typing=Neoantimon2:::allo_records(lapply(f$typing,function(t)t[setdiff(names(t),'raw_alleles')])),
    occurrence_ids=Neoantimon2:::allo_set(lapply(f$source_occurrences,`[[`,'hla_source_occurrence_id')))
  expect_identical(f$allo_analysis_id,Neoantimon2:::allo_id('analysis',payload))
})
