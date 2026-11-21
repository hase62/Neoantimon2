test_that("FASTA production deduplicates candidates while preserving each full input occurrence", {
  p <- review_fasta_project(); native <- neo_results(p); graph <- external_graph(p)
  expect_equal(nrow(native), 2L); expect_length(graph$origins,3L); expect_length(graph$events,0L)
  expect_equal(sort(native$n_origins), c(1L,2L))
  expect_true(all(native$final_grade == "NOT_GRADED"))
  expect_true(all(vapply(graph$origins, function(o) o$origin_type == "sequence" &&
    o$source_domain == "provided_peptide" && is.null(o$transcript_id) && is.null(o$wt_peptide_sequence) &&
    !length(o$source_event_ids), TRUE)))
  peptides <- readRDS(peptide_stage_path(p))$peptides
  expect_setequal(vapply(peptides,function(x)x$provenance$supplied_peptide$record_id,""),paste0("record",1:3))
  expect_true(all(vapply(peptides,function(x)grepl("supplied description",x$provenance$supplied_peptide$full_header),TRUE)))
  expect_true(all(vapply(peptides,function(x)nchar(x$provenance$supplied_peptide$input_sha256)==64L,TRUE)))
  expect_identical(neo_results(neo_project(path=p$path)),native)
  expect_length(readRDS(peptide_stage_path(p))$predictor_inputs,2L)
  expect_equal(nrow(neo_report_bundle(p)$tables$origin_event_links),0L)
  two <- review_fasta_project(sample_ids=c("s1","s2")); expect_equal(nrow(neo_results(two)),4L)
})

test_that("FASTA parser accepts gzip, multiline and CRLF and diagnoses invalid source lines", {
  parse <- function(lines, gzip=FALSE) {
    file <- tempfile(fileext=if(gzip)".fa.gz" else ".fa")
    con <- if(gzip)gzfile(file,"wb") else file(file,"wb")
    if(length(lines))writeChar(paste(lines,collapse="\r\n"),con,eos=NULL,useBytes=TRUE)
    close(con)
    peptide_fasta_read(list(resolved_path=file))
  }
  expect_identical(parse(c(">id description","acdef","ghik"),TRUE)[[1]]$sequence,"ACDEFGHIK")
  cases <- list(list(character(),"Empty FASTA file"),list(c(">a"),"Empty FASTA record"),
    list("ACD","Sequence before"),list(c(">a","AC D"),"whitespace"),list(c(">a","ACX"),"Invalid peptide"),
    list(c(">a","ACD*"),"Invalid peptide"),list(c(">a","ACD",">a other","ACD"),"Duplicate.*a.*line 1"))
  for(x in cases)expect_error(parse(x[[1]]),paste0(":([1-9][0-9]*): .*",x[[2]]))
  expect_identical(parse(c(">u","ACDU"))[[1]]$sequence,"ACDU")
})

test_that("out-of-policy full sequences remain traceable and input content invalidates old runs", {
  p <- review_fasta_project(c("ACD","ACDEFGHIK"));artifact <- readRDS(peptide_stage_path(p))
  expect_equal(nrow(neo_results(p)),1L);expect_length(artifact$origins,2L)
  expect_true(any(vapply(artifact$statuses,function(s)identical(s$reason,"LENGTH_OUTSIDE_CONFIGURED_CLASS_LENGTHS"),TRUE)))
  statuses <- neo_report_bundle(p)$tables$supplied_peptide_status
  expect_equal(nrow(statuses),2L)
  expect_true(any(statuses$reason=="LENGTH_OUTSIDE_CONFIGURED_CLASS_LENGTHS",na.rm=TRUE))
  writeLines(c(">changed","LMNPQRSTV"),p$inputs[[1]]$resolved_path)
  expect_error(neo_run(p),"FASTA changed")
  expect_error(neo_project(path=p$path),"identity/fingerprint|integrity")
})

test_that("FASTA and molecular origins converge without adding DNA support", {
  old <- peptide_test_run(); graph <- external_graph(old); query <- graph$candidates[[1]]
  f <- tempfile();writeLines(c(">supplied same molecular peptide",query$peptide_sequence),f)
  inputs <- lapply(old$inputs,function(x)x[c("input_id","input_role","sample_id","version")])
  for(i in seq_along(inputs))inputs[[i]]$path <- old$inputs[[i]]$resolved_path
  inputs <- c(inputs,list(list(input_id="supplied",input_role="peptide_fasta",sample_id="sample-1",path=f)))
  mixed <- neo_run(neo_project(samples=old$samples,hla=lapply(old$hla,unlist),config=old$config,inputs=inputs))
  now <- external_graph(mixed);expect_setequal(vapply(now$candidates,`[[`,"","candidate_id"),vapply(graph$candidates,`[[`,"","candidate_id"))
  c <- Filter(function(c)c$candidate_id==query$candidate_id,now$candidates)[[1]]
  expect_equal(length(c$origin_ids),length(query$origin_ids)+1L)
  raw <- input_evidence_read(mixed);oldraw <- input_evidence_read(old)
  measurements <- function(rows)lapply(rows,function(r)r[setdiff(names(r),c("run_id","evidence_id"))])
  expect_identical(measurements(raw$event_evidence),measurements(oldraw$event_evidence))
  paths <- Filter(function(x)x$candidate_id==c$candidate_id,raw$candidate_links)[[1]]$origin_paths
  expect_equal(sum(vapply(paths,function(x)length(x$somatic_event_ids)>0L,TRUE)),
    sum(vapply(Filter(function(x)x$candidate_id==c$candidate_id,oldraw$candidate_links)[[1]]$origin_paths,function(x)length(x$somatic_event_ids)>0L,TRUE)))
})

test_that("supplied Origins never confer somatic eligibility or alter independently supported Native Grade", {
  old <- postreview_binding_project()$p;before <- neo_results(old);query <- before$peptide_sequence[[1L]]
  file <- tempfile();writeLines(c(">provided",query),file)
  inputs<-lapply(old$inputs,function(x)list(input_id=x$input_id,input_role=x$input_role,sample_id=x$sample_id,path=x$resolved_path))
  inputs<-c(inputs,list(list(input_id="fasta-control",input_role="peptide_fasta",sample_id="sample-1",path=file)))
  mixed<-neo_run(neo_project(samples=old$samples,hla=lapply(old$hla,unlist),config=old$config,inputs=inputs))
  after<-neo_results(mixed)
  for(field in c("candidate_id","final_grade","grade_status","technical_eligibility","priority_origin_id","primary_rank"))
    expect_identical(after[[field]],before[[field]])
  cfg<-neo_config(env=character(),paths=list(run_root=tempfile()),somatic_evidence=list())
  supplied<-review_fasta_project(config=cfg);r<-neo_results(supplied)
  expect_true(all(r$final_grade=="NOT_GRADED"));expect_true(all(r$technical_eligibility=="not_evaluated"))
})
