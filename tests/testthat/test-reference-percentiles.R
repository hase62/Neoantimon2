test_that("strict better midrank uses binary bounds with ties, endpoints and valid zero", {
  x <- c(0,1,1,2)
  expect_equal(vapply(c(-1,0,1,2,3),function(q)reference_percentile_value(x,q,"lower_is_better"),0),c(0,12.5,50,87.5,100))
  expect_equal(vapply(c(-1,0,1,2,3),function(q)reference_percentile_value(x,q,"higher_is_better"),0),c(100,87.5,50,12.5,0))
})

test_that("calibration is a separate persisted measurement table with unchanged no-option fingerprints", {
  p <- review_fasta_project(sample_ids=c("s1","s2"));x <- external_test_row(score="1",score_direction="lower_is_better")
  file <- external_test_write(x);native <- neo_results(p);before <- external_sha(manifest_file(p$path))
  old <- neo_analyze_external(p,file,list(mode="evidence_only"))
  a <- neo_analyze_external(p,file,list(mode="evidence_only"),reference_percentiles=review_reference())
  z <- external_analysis_read(a);r <- z$payload$reference_percentiles
  expect_length(r,1L);expect_length(r[[1]]$candidate_ids,2L);expect_equal(r[[1]]$reference_percentile,50)
  expect_equal(r[[1]]$raw_score,1);expect_identical(r[[1]]$original_score_text,"1")
  expect_equal(z$payload$import_report$matched_input_rows,1L)
  expect_equal(nrow(neo_report_bundle(p,analyses=list(a))$tables$external_reference_percentiles),1L)
  expect_length(neo_external_explain(a,r[[1]]$candidate_ids[[1]])$candidates[[1]]$reference_percentiles,1L)
  expect_identical(neo_external_results(a$analysis_id,neo_project(path=p$path)),neo_external_results(a))
  legacy <- z$manifest$fingerprint;legacy$reference_percentiles<-NULL
  expect_identical(old$analysis_id,infra_hash(legacy))
  expect_identical(neo_analyze_external(p,file,list(mode="evidence_only"))$analysis_id,old$analysis_id)
  expect_identical(neo_results(p),native);expect_identical(external_sha(manifest_file(p$path)),before)
  expect_equal(neo_external_results(old)$external_score,neo_external_results(a)$external_score)
  writeLines("tampered",file.path(a$path,"reference_percentiles.tsv"));expect_error(neo_external_results(a),"INTEGRITY_FAILURE")
})

test_that("calibration refuses invalid distributions before committing and never falls back", {
  p<-review_fasta_project();x<-external_test_row(score="0",score_direction="lower_is_better")
  apply<-function(x,path=review_reference())neo_analyze_external(p,external_test_write(x),list(mode="evidence_only"),reference_percentiles=path)
  for (change in list(list(model_version=""),list(model_version="other"),list(mhc_key="HLA-B*07:02"),list(peptide_sequence="ACDEFGHI"),list(task="external:other"))) {
    y<-x;y[names(change)]<-change;r<-external_analysis_read(apply(y))$payload$reference_percentiles[[1]]
    expect_identical(r$status,"NOT_CALCULATED");expect_null(r$reference_percentile)
  }
  y<-x;y$status<-"not_evaluated";y$score<-"";y$reason<-"not supplied"
  expect_match(external_analysis_read(apply(y))$payload$reference_percentiles[[1]]$reason,"SOURCE_NOT_EVALUATED")
  expect_equal(external_analysis_read(apply(x))$payload$reference_percentiles[[1]]$reference_percentile,12.5)
  for (bad in list(review_reference(count=5L),review_reference(scores_sha256=paste(rep("a",64),collapse="")),review_reference(scores=c(0,Inf))))
    expect_error(apply(x,bad),"INVALID_REFERENCE_PERCENTILES")
  ref<-review_reference();m<-jsonlite::read_json(ref);m$distributions<-rep(m$distributions,2L);atomic_json(m,ref)
  expect_error(apply(x,ref),"Duplicate")
})

test_that("review evidence ignores display ordering, tracks changed analyses and is private only", {
  p<-review_fasta_project();b<-neo_report_bundle(p);x<-html_report_payload(b,privacy="private")
  changed<-b;changed$tables<-lapply(changed$tables,function(d)d[rev(seq_len(nrow(d))),,drop=FALSE])
  a<-x$review$evidence;a<-a[order(vapply(a,`[[`,"","candidate_id"))]
  z<-html_report_payload(changed,privacy="private")$review$evidence;z<-z[order(vapply(z,`[[`,"","candidate_id"))]
  expect_identical(a,z);expect_null(html_report_payload(b,privacy="share")$review)
  y<-external_test_row();analysis<-neo_analyze_external(p,external_test_write(y),list(mode="evidence_only"))
  now<-html_report_payload(neo_report_bundle(p,analyses=list(analysis)),privacy="private")$review
  expect_false(identical(x$review$evidence,now$evidence))
  one <- b$tables$pmhc_candidates$candidate_id[[1L]]
  full <- html_report_payload(neo_report_bundle(p,analyses=list(analysis)),privacy="private")$review$evidence
  partial <- html_report_payload(neo_report_bundle(p,analyses=list(analysis),candidate_id=one),privacy="private")$review$evidence
  expect_identical(Filter(function(r)r$candidate_id==one,full)[[1L]],partial[[1L]])
  expect_false(any(grepl("review",names(b$tables),ignore.case=TRUE)))
  shared_json <- html_report_json(html_report_payload(b,privacy="share"))
  expect_false(any(vapply(b$tables$pmhc_candidates$candidate_id,function(id)grepl(id,shared_json,fixed=TRUE),TRUE)))
})
