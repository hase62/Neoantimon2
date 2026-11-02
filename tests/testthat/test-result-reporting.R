test_that("all four levels and complete links survive reporting with orphan Origins", {
  f<-report_test_fixture();x<-report_test_build(f);t<-x$tables
  expect_equal(nrow(t$events),6);expect_equal(nrow(t$origins),5)
  expect_equal(nrow(t$pmhc_candidates),4);expect_equal(nrow(t$candidate_origin_links),6)
  expect_equal(nrow(t$peptide_summary),3);expect_equal(nrow(t$pmhc_predictions),12)
  expect_setequal(t$candidate_origin_links$origin_id,c("o1","o2","o3","o4"))
  expect_identical(t$origins$nested[[1]],f$graph$origins[[1]]$nested)
  expect_true(all(is.na(t$pmhc_candidates$evidence_category)))
  expect_true(all(t$pmhc_candidates$evidence_status=="REFERENCE_UNAVAILABLE"))
  expect_true(all(is.na(t$pmhc_predictions$raw_score[t$pmhc_predictions$record_kind=="coverage_status"])))
  expect_equal(t$event_summary$n_pmhc_candidates[t$event_summary$event_id=="empty"],0)
  expect_equal(t$event_summary$n_distinct_peptides[t$event_summary$event_id=="e4"],1)
})
test_that("Native and exact external scopes select independently without cross-class ranking", {
  t<-report_test_build()$tables;b<-t$event_best_by_channel
  get<-function(ev,source,cl,analysis=NULL) b[b$event_id==ev&b$source==source&b$hla_class==cl&
    (if(is.null(analysis))rep(TRUE,nrow(b)) else (!is.na(b$analysis_id)&b$analysis_id==analysis)),,drop=FALSE]
  n<-get("e1","Native","I");expect_identical(n$best_candidate_id,"c2")
  expect_identical(n$best_origin_id,"o1");expect_false(n$best_origin_is_native_priority_origin)
  expect_setequal(n$joint_source_event_refs[[1]],vapply(c("e1","e2"),function(id)Neoantimon2:::report_event_ref("core",id),""))
  expect_true(n$compound_source_context);expect_false(n$independent_causation_claimed)
  expect_identical(get("e1","External","I","high")$best_candidate_id,"c1")
  expect_identical(get("e1","External","I","low")$best_candidate_id,"c2")
  expect_identical(get("e1","Native","II")$best_candidate_id,"c4")
  expect_identical(get("e3","Native","I")$status,"UNAVAILABLE")
  expect_true(all(b$status[b$event_id=="empty"]=="UNAVAILABLE"))
  expect_equal(t$event_summary$n_pmhc_candidates[t$event_summary$event_id=="g1"],0)
  expect_identical(t$event_summary$status[t$event_summary$event_id=="g1"],"CONTEXT_ONLY_NOT_CAUSATIVE")
})
test_that("subset bounds do not drop Event nodes or non-primary HLA alternatives", {
  t<-report_test_build(candidate_id=c("c1","c2"))$tables
  expect_equal(nrow(t$events),6);expect_equal(nrow(t$origins),2);expect_equal(nrow(t$peptide_summary),1)
  expect_equal(nrow(t$candidate_origin_links),4)
  expect_setequal(t$peptide_summary$alternative_HLAs[[1]],c("HLA-A*02:01","HLA-B*07:02"))
  expect_error(report_test_build(candidate_id="unknown"),"absent")
})
test_that("link, sample, score and namespace corruption fails explicitly", {
  f<-report_test_fixture();f$graph$links<-f$graph$links[-1];expect_error(report_test_build(f),"linkage mismatch")
  f<-report_test_fixture();f$graph$origins[[1]]$sample_id<-"s2";expect_error(report_test_build(f),"sample or peptide")
  f<-report_test_fixture();f$graph$events[[1]]$sample_id<-"s2";expect_error(report_test_build(f),"Event/Origin sample")
  f<-report_test_fixture();f$external[[1]]$payload$evidence[[1]]$numeric_score<-NA_real_;expect_error(report_test_build(f),"finite")
  f<-report_test_fixture();dup<-f$external[[1]]$payload$evidence[[1]];dup$numeric_score<-.9
  f$external[[1]]$payload$evidence<-c(f$external[[1]]$payload$evidence,list(dup));expect_error(report_test_build(f),"Conflicting scores")
  expect_error(report_test_build(external_channels=list(list(provider="absent"))),"not in")
})
test_that("raw records preserve JSON null, empty array, attributes and numeric precision", {
  t<-report_test_build()$tables;path<-tempfile(fileext=".tsv")
  m<-Neoantimon2:::report_write_table(t$origins,path);d<-read.delim(path,check.names=FALSE,na.strings="")
  raw<-jsonlite::fromJSON(d$source_record_json[1],simplifyVector=FALSE)
  expect_null(raw$nested$null);expect_identical(raw$nested$empty,list());expect_true(raw$nested$truth)
  expect_identical(m$sha256,Neoantimon2:::external_sha(path));expect_equal(m$rows,5)
  z<-data.frame(raw_score=.12345678901234567);Neoantimon2:::report_write_table(z,path)
  expect_identical(read.delim(path)$raw_score,z$raw_score)
})
test_that("same literal Event IDs in source and long-read stay separate and explicitly mapped", {
  f<-report_test_fixture();f$graph$annotations<-list(source=list(canonical_events=list(list(source_event_id="same",structure=list(sample_id="s1"))),
    peptides=list(list(origin_id="o1",source_event_id="same"))),long_read=list(events=list(list(long_read_source_event_id="same",context=list(sample_id="s1"))),
    peptides=list(list(origin_id="o2",long_read_source_event_id="same"))))
  t<-report_test_build(f)$tables;l<-t$origin_event_links[t$origin_event_links$event_id=="same",]
  expect_setequal(l$event_domain,c("source","long_read"));expect_equal(length(unique(l$event_ref)),2)
  f$graph$annotations$source$peptides<-NULL;f$graph$annotations$long_read$peptides<-NULL
  f$graph$origins[[1]]$source_antigen_key<-"same"
  expect_error(report_test_build(f),"Ambiguous structural")
})
test_that("all model metrics stay long and model selection only changes representatives", {
  f<-report_test_fixture();base<-report_test_build(f)
  selected<-report_test_build(f,external_channels=list(list(analysis_id="high")))
  expect_identical(base$tables$pmhc_predictions,selected$tables$pmhc_predictions)
  expect_equal(sum(selected$tables$event_best_by_channel$source=="External"),12)
  expect_true(all(selected$tables$event_best_by_channel$analysis_id[selected$tables$event_best_by_channel$source=="External"]=="high"))
  duplicate<-f$external[[1]]$payload$evidence[[1]]
  f$external[[1]]$payload$evidence<-c(f$external[[1]]$payload$evidence,list(duplicate))
  expect_equal(nrow(report_test_build(f)$tables$pmhc_predictions),13)
})
test_that("yield plan and all pair evidence remain unchanged while 1 sequence has multiple HLAs", {
  f<-report_test_fixture();f$graph$candidates<-f$graph$candidates[1:2];f$graph$links<-f$graph$links[1:4]
  x<-yield_test_scores(c(.8,.7));x$peptide_sequence<-"ACDEFGHIK";x$mhc_key<-c("HLA-A*02:01","HLA-B*07:02")
  x$origin_id<-NULL;x$origin_ids<-I(rep(list(list("o1","o2")),2))
  y<-neo_experimental_yield(x,yield_test_profile());before<-digest::digest(y)
  t<-report_test_build(f,yield_plan=y)$tables
  expect_equal(nrow(t$pmhc_candidates),2);expect_equal(nrow(t$experimental_yield_pair_evidence),2)
  planned<-t$peptide_summary[t$peptide_summary$planning_status=="SINGLE_HLA_PROXY_FROM_UNCHANGED_YIELD_PLAN",]
  expect_equal(nrow(planned),1);expect_equal(planned$assumed_or_empirical_rate,y$units$probability)
  expect_identical(t$experimental_yield_curve,y$curve);expect_identical(t$experimental_yield_targets,y$targets)
  expect_identical(digest::digest(y),before)
  bad<-y;bad$pair_evidence$mhc_key[1]<-"HLA-C*07:01"
  expect_error(report_test_build(f,yield_plan=bad),"identity mismatch")
})
test_that("public report is read-only, deterministic, safe and keeps legacy export identity", {
  p<-peptide_test_run();legacy<-neo_export_bundle(p)$manifest;before<-neo_results(p)
  hashes<-function() {
    ff<-sort(list.files(p$path,recursive=TRUE,all.files=TRUE,full.names=TRUE));ff<-ff[!dir.exists(ff)]
    stats::setNames(vapply(ff,function(f)digest::digest(file=f,algo="sha256",serialize=FALSE),""),ff)
  }
  frozen<-hashes();dir<-tempfile("report-");r<-neo_report_bundle(p,dir)
  expect_s3_class(r,"neo_report_bundle");expect_identical(r$manifest$contract,"result-report-0.1.0")
  expect_equal(nrow(r$tables$pmhc_candidates),nrow(before))
  expect_true(all(r$tables$event_best_by_channel$status=="UNAVAILABLE"))
  for(m in r$manifest$tables) {
    expect_identical(Neoantimon2:::external_sha(file.path(dir,m$file)),m$sha256)
    expect_equal(nrow(read.delim(file.path(dir,m$file),check.names=FALSE)),m$rows)
  }
  expect_identical(jsonlite::read_json(file.path(dir,"manifest.json"))$report_id,r$manifest$report_id)
  r2<-neo_report_bundle(neo_project(path=p$path),tempfile("report-repeat-"))
  expect_identical(r$manifest$report_id,r2$manifest$report_id)
  expect_identical(neo_results(p),before);expect_identical(hashes(),frozen)
  expect_identical(neo_export_bundle(p)$manifest,legacy)
  expect_error(neo_report_bundle(p,dir),"already exists")
  expect_error(neo_report_bundle(p,file.path(p$path,"unsafe-report")),"overlap")
  expect_error(neo_report_bundle(p,analyses=list("absent")),".")
})
test_that("empty completed projects keep explicit table schemas and no invented representatives", {
  p<-neo_run(neo_project(config=neo_config(env=character(),paths=list(run_root=tempfile("empty-report-"))),fixture_kind="synthetic"))
  x<-neo_report_bundle(p,tempfile("empty-report-export-"))
  expect_equal(nrow(x$tables$events),0);expect_equal(nrow(x$tables$pmhc_candidates),0)
  expect_equal(nrow(x$tables$event_best_by_channel),0)
  expect_true(all(c("candidate_id","raw_score","prediction_status") %in% names(x$tables$pmhc_predictions)))
})
test_that("frozen native measurement records and provider statuses remain separate", {
  f<-report_test_fixture();e<-list(candidate_id="c1",predictor="synthetic-provider",predictor_version="1",model_version="1",
    task="mhc1_binding",metric="affinity",unit="nM",raw_score=42,score_direction="lower_is_better",status="evaluated")
  runs<-list(list(run_id="native-run",provider="synthetic-provider",evidence=list(e),
    statuses=list(list(candidate_id="c2",status="unsupported",reason="synthetic-unsupported"))))
  x<-Neoantimon2:::report_construct(f$graph,f$native,runs=runs)$tables
  pr<-x$pmhc_predictions
  expect_equal(pr$raw_score[pr$record_kind=="measurement"],42)
  expect_identical(pr$prediction_status[pr$record_kind=="provider_status"],"unsupported")
  expect_equal(nrow(x$native_run_statuses),1)
})
test_that("sampling and MHC class scopes prevent accidental pooled ranking", {
  f<-report_test_fixture();f$graph$events[[length(f$graph$events)+1]]<-list(event_id="other",sample_id="s2",event_role="somatic")
  o<-f$graph$origins[[3]];o$origin_id<-"s2-origin";o$sample_id<-"s2";o$source_event_ids<-list("other");o$germline_context_ids<-list()
  c<-f$graph$candidates[[3]];c$candidate_id<-"s2-candidate";c$sample_id<-"s2";c$origin_ids<-list(o$origin_id)
  f$graph$origins<-c(f$graph$origins,list(o));f$graph$candidates<-c(f$graph$candidates,list(c))
  f$graph$links<-c(f$graph$links,list(list(candidate_id=c$candidate_id,origin_id=o$origin_id,sample_id="s2")))
  n<-f$native[[1]];n$candidate_id<-c$candidate_id;n$sample_id<-"s2";n$priority_origin_id<-o$origin_id;n$rank_within_sample_class<-1
  f$native<-c(f$native,list(n));b<-report_test_build(f)$tables$event_best_by_channel
  expect_identical(b$best_candidate_id[b$event_id=="other"&b$source=="Native"],"s2-candidate")
  expect_false(any(b$best_candidate_id[b$sample_id=="s1"]=="s2-candidate",na.rm=TRUE))
})
test_that("report-only hashing and JSON preserve array-valued scientific fields", {
  x<-list(positions=c(1L,3L),weights=c(.12345678901234567,.9),empty=list(),unknown=NULL)
  encoded<-Neoantimon2:::report_json(x);y<-jsonlite::fromJSON(encoded)
  expect_identical(y$weights,x$weights);expect_equal(y$positions,x$positions)
  expect_identical(Neoantimon2:::report_hash(x),Neoantimon2:::report_hash(x))
  y<-x;y$weights[1]<-.1234567890123459
  expect_false(identical(Neoantimon2:::report_hash(x),Neoantimon2:::report_hash(y)))
})
test_that("unresolved source metadata stays unavailable and never invents an Event representative", {
  f<-report_test_fixture();f$graph$origins[[1]]$source_event_ids<-list("missing")
  t<-report_test_build(f)$tables
  expect_identical(t$events$record_status[t$events$event_id=="missing"],"SOURCE_EVENT_RECORD_UNAVAILABLE")
  expect_true(all(t$event_best_by_channel$status[t$event_best_by_channel$event_id=="missing"]=="UNAVAILABLE"))
  expect_true(any(t$origin_event_links$event_id=="missing"))
})
test_that("reported experimental matches preserve evidence and never alter the yield calculation", {
  f<-report_test_fixture();x<-yield_test_scores(c(.8,.7));x$peptide_sequence<-"ACDEFGHIK"
  x$mhc_key<-c("HLA-A*02:01","HLA-B*07:02");x$origin_id<-NULL;x$origin_ids<-I(rep(list(list("o1","o2")),2))
  f$graph$candidates<-f$graph$candidates[1:2];f$graph$links<-f$graph$links[1:4]
  ref<-yield_test_reference(list(yield_test_record("a-positive"),yield_test_record("a-negative",result="negative"),
    yield_test_record("similar",pep="ACDEYGHIK")))
  y<-neo_experimental_yield(x,yield_test_profile(),reference=ref)
  t<-report_test_build(f,yield_plan=y)$tables
  expect_identical(t$experimental_yield_curve,y$curve);expect_identical(t$experimental_yield_summary,y$summary)
  expect_identical(t$pmhc_candidates$A_experimental_status[1],"conflicting")
  expect_true("B1" %in% t$pmhc_candidates$B_subtypes[[2]])
  expect_length(t$pmhc_candidates$experimental_matches[[1]],3)
  other<-ref;other$records[[1]]$citation<-"different"
  expect_error(report_test_build(f,yield_plan=y,reference=other),"differs")
})

test_that("named empty objects and every list cell are valid strict JSON", {
  x<-list(array=list(),object=stats::setNames(list(),character()))
  encoded<-Neoantimon2:::report_json(x)
  expect_identical(encoded,'{"array":[],"object":{}}')
  expect_true(jsonlite::validate(encoded))
  f<-report_test_fixture();f$graph$origins[[1]]$source_metadata_json<-stats::setNames(list(),character())
  t<-report_test_build(f)$tables
  for(table in t)for(col in table)if(is.list(col))expect_true(all(vapply(col,function(v)jsonlite::validate(Neoantimon2:::report_json(v)),TRUE)))
  for(raw in t$origins$source_record_json)expect_true(jsonlite::validate(raw))
})
