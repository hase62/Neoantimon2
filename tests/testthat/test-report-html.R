html_fixture_bundle <- function() {
  x<-report_test_build();list(tables=x$tables,manifest=c(list(contract="result-report-0.1.0",package_version="test",source_manifest=list(path="/Users/PRIVATE_PATIENT/canary",sha256=paste(rep("a",64),collapse=""))),x$metadata))
}
test_that("HTML is a view of every authoritative Origin and namespace link",{
  b<-html_fixture_bundle();before<-serialize(b,NULL);x<-html_report_payload(b,privacy="private")
  expect_identical(serialize(b,NULL),before)
  for(n in setdiff(names(b$tables),"pmhc_candidates"))expect_identical(x$tables[[n]],report_rows(b$tables[[n]]))
  for(i in seq_len(nrow(b$tables$pmhc_candidates)))expect_identical(x$tables$pmhc_candidates[[i]][names(b$tables$pmhc_candidates)],report_rows(b$tables$pmhc_candidates)[[i]])
  expect_equal(length(x$tables$origins),5);expect_equal(length(x$tables$candidate_origin_links),6)
  expect_setequal(vapply(Filter(function(l)l$candidate_id=="c1",x$tables$candidate_origin_links),`[[`,"","origin_id"),c("o1","o2"))
  expect_true(any(vapply(x$tables$origin_event_links,function(l)isTRUE(l$compound_source_context),TRUE)))
  expect_true(any(vapply(x$tables$origin_event_links,function(l)l$link_role=="germline_context",TRUE)))
  bad<-b;bad$tables$candidate_origin_links$origin_id[[1]]<-"absent";expect_error(html_report_payload(bad),"linkage")
})
test_that("share mode keeps all scientific associations but masks private content",{
  b<-html_fixture_bundle();b$tables$origins$source_metadata_json<-I(rep(list(list(secret="PRIVATE_PATIENT",path="/Users/PRIVATE_PATIENT/cache")),nrow(b$tables$origins)))
  b$tables$pmhc_candidates$patient_id<-"PRIVATE_PATIENT"
  x<-html_report_payload(b,list(subjects=c("s1","PRIVATE_PATIENT"),references=list(list(version="/Users/PRIVATE_PATIENT/secret"))),"share")
  json<-html_report_json(x)
  expect_false(grepl("PRIVATE_PATIENT|/Users/|source_metadata_json|source_record_json",json))
  expect_equal(length(x$tables$origins),5);expect_equal(length(x$tables$candidate_origin_links),6)
  map<-attr(x,"private_id_map");expect_gt(length(map),0)
  m<-function(value,family)Filter(function(r)r$original==value&&r$family==family,map)[[1]]$alias
  rows<-Filter(function(l)l$candidate_id==m("c1","Candidate"),x$tables$candidate_origin_links)
  expect_setequal(vapply(rows,`[[`,"","origin_id"),c(m("o1","Origin"),m("o2","Origin")))
  expect_true(all(vapply(x$tables$origin_event_links,function(l)l$event_ref%in%vapply(x$tables$events,`[[`,"","event_ref"),TRUE)))
  expect_true(any(unlist(x$metadata$withheld_field_counts)>0L))
  # Pseudonymization must not corrupt a scientific sequence matching a declared ID.
  y<-html_report_payload(b,list(subjects=c("ACDEFGHIK","s1")),"share")
  expect_identical(y$tables$pmhc_candidates[[1]]$peptide_sequence,"ACDEFGHIK")
  expect_identical(html_report_json(x),html_report_json(html_report_payload(b,list(subjects=c("s1","PRIVATE_PATIENT"),references=list(list(version="/Users/PRIVATE_PATIENT/secret"))),"share")))
})
test_that("Native public status overrides snapshot only in explicit HTML display fields",{
  b<-html_fixture_bundle();b$tables$pmhc_candidates$prediction_status<-"NOT_RUN"
  b$tables$pmhc_candidates$native_public_result[[1]]$prediction_status<-"SUCCESS"
  b$tables$pmhc_candidates$native_public_result[[1]]$grading_status<-"GRADED"
  x<-html_report_payload(b,privacy="share")
  expect_identical(x$tables$pmhc_candidates[[1]]$html_native_prediction_status,"SUCCESS")
  expect_identical(x$tables$pmhc_candidates[[1]]$prediction_status,"NOT_RUN")
  expect_identical(x$tables$pmhc_candidates[[1]]$html_native_grade_status,"GRADED")
})
test_that("zero, missing scores and independent channels remain distinct",{
  f<-report_test_fixture();f$external[[1]]$payload$evidence[[1]]$numeric_score<-0;f$external[[1]]$payload$evidence[[1]]$score<-"0.0000000000000000"
  x<-report_test_build(f);b<-html_fixture_bundle();b$tables<-x$tables
  p<-html_report_payload(b,privacy="private")$tables$pmhc_predictions
  raw<-report_rows(b$tables$pmhc_predictions)
  for(i in seq_along(p))expect_identical(p[[i]],raw[[i]])
  expect_true(any(vapply(p,function(r)identical(r$original_score_text,"0.0000000000000000")&&r$raw_score==0,TRUE)))
  expect_true(all(vapply(Filter(function(r)r$record_kind=="coverage_status",p),function(r)is.na(r$raw_score),TRUE)))
  expect_setequal(unique(vapply(Filter(function(r)r$source=="External",p),`[[`,"","score_direction")),c("higher_is_better","lower_is_better"))
})
test_that("malicious JSON script terminators and template markers stay data",{
  b<-html_fixture_bundle();evil<-"</script><img src=https://invalid.example/x onerror=globalThis.PWNED=1><script>globalThis.PWNED=2</script> & \" {{DATA}} {{SCRIPT}} {{NONCE}} \u2028\u2029"
  b$tables$origins$gene_symbol[[1]]<-evil
  json<-html_report_json(html_report_payload(b,privacy="private"))
  expect_false(grepl("</script>|<img|<script>",json,ignore.case=TRUE))
  expect_true(jsonlite::validate(json));expect_match(json,"\\\\u003c")
  decoded<-jsonlite::fromJSON(json,simplifyVector=FALSE);expect_identical(decoded$tables$origins[[1]]$gene_symbol,evil)
})
test_that("writer is bounded, read-only, exclusive and compatible with old reporting",{
  e<-new.env();sys.source(system.file("examples","first-real-run.R",package="Neoantimon2"),e);a<-e$first_real_run_example(tempfile());p<-neo_run(do.call(neo_project,a))
  old<-neo_report_bundle(p);saved<-ia_files(p$path);out<-tempfile(fileext=".html")
  result<-html_report_write(p,out,label="Authored synthetic software example; no real predictions")
  expect_true(file.exists(out));expect_identical(result$privacy,"share");expect_equal(result$counts$pmhc_candidates,65)
  expect_identical(neo_report_bundle(p),old);expect_identical(ia_files(p$path),saved)
  expect_equal(length(getNamespaceExports("Neoantimon2")),22)
  expect_identical(neo_results(neo_project(path=p$path)),neo_results(p))
  expect_error(html_report_write(p,out),"already exists")
  expect_error(html_report_write(p,tempfile(),privacy="private"),"confirm_private")
  expect_warning(html_report_write(p,tempfile(fileext=".html"),privacy="private",confirm_private=TRUE),"PRIVATE HTML")
  expect_error(html_report_write(p,tempfile(),max_candidates=1),"limit")
  small<-tempfile();expect_error(html_report_write(p,small,max_bytes=1024),"byte limit");expect_false(file.exists(small))
  expect_error(html_report_write(p,file.path(p$path,"overlap.html")),"overlap")
  expect_error(html_report_write(p,tempfile(),report_options=list(dir="forbidden")),"Unknown key")
  expect_error(html_report_write(p,tempfile(),privacy=NA_character_),"privacy")
  subset<-html_report_write(p,tempfile(fileext=".html"),candidate_id=old$tables$pmhc_candidates$candidate_id[[1]])
  expect_match(subset$scope,"EXPLICIT_CANDIDATE_SELECTION");expect_equal(subset$counts$pmhc_candidates,1)
  expect_equal(subset$counts$origins,sum(old$tables$candidate_origin_links$candidate_id==old$tables$pmhc_candidates$candidate_id[[1]]))
  doc<-paste(readLines(out,warn=FALSE),collapse="\n")
  expect_match(doc,"connect-src 'none'",fixed=TRUE);expect_match(doc,"Content-Security-Policy",fixed=TRUE)
  expect_false(grepl('/Users/|/private/var/|patient-1|sample-1|demo-sample|demo-subject',doc))
  expect_false(grepl('src="https?://|href="https?://',doc))
  expect_false(grepl("innerHTML|eval\\(|new Function|XMLHttpRequest|fetch\\(",doc))
})


test_that("sample units, HLA alternatives and unscored scope stay distinct",{
 f<-report_test_fixture();c<-f$graph$candidates[[1]];c$candidate_id<-"sample2-c1";c$sample_id<-"s2";c$origin_ids<-list("sample2-o1")
 o<-f$graph$origins[[1]];o$origin_id<-"sample2-o1";o$sample_id<-"s2";o$source_event_ids<-list("sample2-e1");o$germline_context_ids<-list("sample2-g1")
 e<-f$graph$events[[1]];e$event_id<-"sample2-e1";e$sample_id<-"s2"
 g<-e;g$event_id<-"sample2-g1";g$event_role<-"germline"
 f$graph$candidates<-c(f$graph$candidates,list(c));f$graph$origins<-c(f$graph$origins,list(o));f$graph$events<-c(f$graph$events,list(e,g))
 f$graph$links<-c(f$graph$links,list(list(candidate_id=c$candidate_id,origin_id=o$origin_id,sample_id="s2")))
 n<-f$native[[1]];n$candidate_id<-c$candidate_id;n$sample_id<-"s2";n$priority_origin_id<-o$origin_id;f$native<-c(f$native,list(n))
 z<-report_test_build(f);b<-html_fixture_bundle();b$tables<-z$tables;b$manifest$original_project_graph_counts<-z$metadata$original_project_graph_counts
 x<-html_report_payload(b,list(subjects=c("s1","s2")),"share")
 expect_equal(x$metadata$scoped_table_rows$pmhc_candidates,5);expect_equal(x$metadata$scoped_table_rows$peptide_summary,4)
 expect_equal(x$metadata$scoped_table_rows$unique_peptide_sequences,3)
 expect_equal(length(unique(vapply(x$tables$pmhc_candidates,`[[`,"","sample_id"))),2)
 expect_equal(sum(vapply(x$tables$pmhc_candidates,function(c)c$peptide_sequence=="ACDEFGHIK",TRUE)),4)
 expect_true(any(vapply(x$tables$pmhc_predictions,function(r)r$prediction_status=="NOT_PREDICTED"&&is.na(r$raw_score),TRUE)))
})

test_that("zero candidate output preserves empty scope without inventing measurements",{
 cfg<-neo_config(env=character(),paths=list(run_root=tempfile()))
 p<-neo_run(neo_project(config=cfg,fixture_kind="synthetic"));b<-neo_report_bundle(p)
 expect_equal(nrow(b$tables$pmhc_candidates),0)
 out<-html_report_write(p,tempfile(fileext=".html"))
 expect_equal(out$counts$pmhc_candidates,0)
 expect_match(paste(readLines(out$path,warn=FALSE),collapse="\n"),"this is not a negative prediction",fixed=TRUE)
})


test_that("identifiers coinciding with scientific values do not alter those values",{
 f<-report_test_fixture();f$external[[1]]$payload$evidence[[1]]$numeric_score<-0;f$external[[1]]$payload$evidence[[1]]$score<-"0"
 z<-report_test_build(f);b<-html_fixture_bundle();b$tables<-z$tables
 b$tables$pmhc_candidates$prediction_status<-"NOT_RUN"
 x<-html_report_payload(b,list(subjects=c("0","NOT_RUN","HLA","synthetic-model","GENE","somatic","core")),"share")
 expect_identical(x$tables$pmhc_candidates[[1]]$prediction_status,"NOT_RUN")
 expect_identical(x$tables$pmhc_candidates[[1]]$mhc_key,"HLA-A*02:01")
 expect_identical(x$tables$origins[[1]]$gene_symbol,"GENE")
 zero<-Filter(function(r)identical(r$raw_score,0),x$tables$pmhc_predictions)[[1]]
 expect_identical(zero$original_score_text,"0");expect_identical(zero$provider,"synthetic-model")
 x<-html_report_share(list(example=data.frame(gene_symbol="C:\\Users\\PRIVATE_PATIENT\\secret")),"PRIVATE_PATIENT")
 expect_false(grepl("PRIVATE_PATIENT|C:",html_report_json(x$tables)))
})
