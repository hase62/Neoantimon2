first_ux_plan <- function() {
  e <- new.env(parent = globalenv())
  sys.source(system.file("examples", "first-real-run.R", package = "Neoantimon2"), e)
  list(helper = e, plan = e$first_real_run_example(tempfile("first-real-")))
}
first_ux_check <- function(a, ...) do.call(first_run_preflight, c(a, list(...)))

test_that("effective default lengths are displayed without changing configuration", {
  x <- first_ux_plan();a <- x$plan;a$config$peptides <- NULL
  before <- serialize(a, NULL);z <- first_ux_check(a)
  expect_true(z$ok)
  expect_identical(unlist(z$effective$peptides$class_i_lengths), 8:11)
  expect_identical(unlist(z$effective$peptides$class_ii_lengths), 13:25)
  expect_identical(z$effective$lengths_source, "package_default")
  expect_identical(serialize(a,NULL), before)
  expect_null(z$private);expect_equal(z$software$exported_api_count,22L)
  a$config <- neo_config(env=character(),paths=a$config$paths,resources=a$config$resources,
    scientific=a$config$scientific,peptides=list(class_i_lengths=8:12))
  z<-first_ux_check(a)
  expect_identical(unlist(z$effective$peptides$class_i_lengths),8:12)
  expect_identical(unlist(z$effective$peptides$class_ii_lengths),13:25)
})

test_that("preflight reads declarations without a process, write or input change", {
  x<-first_ux_plan();a<-x$plan;root<-tempfile();dir.create(root)
  canary<-file.path(root,"started");exe<-file.path(root,"predictor")
  writeLines(c("#!/bin/sh",paste("touch",shQuote(canary)),"echo 'NetMHCpan version 4.2'"),exe);Sys.chmod(exe,"0755")
  a$config$predictors<-list(list(provider="netmhcpan",task="mhc1_binding",executable=exe,options=list()))
  testthat::local_mocked_bindings(neo_check_predictors=function(...)stop("Unapproved probe"),
    neo_run=function(...)stop("Scientific execution"),neo_download_resources=function(...)stop("Download"),.package="Neoantimon2")
  paths<-c(a$inputs[[1]]$path,vapply(a$config$resources,`[[`,"","path"))
  bytes<-vapply(paths,function(p)digest::digest(file=p,algo="sha256",serialize=FALSE),"")
  mtimes<-file.info(paths)$mtime;before<-serialize(a,NULL);z<-first_ux_check(a)
  expect_true(z$ok);expect_false(file.exists(canary));expect_false(dir.exists(a$config$paths$run_root))
  expect_true(any(z$checks$observed=="EXECUTABLE_PRESENT_NOT_PROBED"))
  expect_identical(before,serialize(a,NULL));expect_identical(file.info(paths)$mtime,mtimes)
  expect_identical(vapply(paths,function(p)digest::digest(file=p,algo="sha256",serialize=FALSE),""),bytes)
})

test_that("bad declarations have safe field-specific actionable failures", {
  a<-first_ux_plan()$plan
  a$inputs[[1]]$path<-file.path(tempdir(),"PRIVATE-PATIENT-TOKEN","missing.vcf")
  z<-first_ux_check(a);expect_false(z$ok)
  row<-z$checks[z$checks$field=="inputs[[1]]$path",];expect_identical(row$status,"FAIL")
  expect_match(row$action,"file.exists")
  expect_error(do.call(neo_project,a),"Required input file missing.*inputs\\[\\[1\\]\\]\\$path")
  e<-tryCatch(do.call(neo_project,a),error=identity);expect_false(grepl("PRIVATE-PATIENT",conditionMessage(e)))
  a<-first_ux_plan()$plan;a$inputs[[1]]$sample_id<-"PRIVATE-OTHER"
  z<-first_ux_check(a);expect_false(z$ok);expect_true(any(z$checks$observed=="BINDING_MISSING"))
  expect_error(do.call(neo_project,a),"Input sample binding missing.*inputs")
  a<-first_ux_plan()$plan;a$hla[[1]]<-"DRB1*01:01"
  expect_false(first_ux_check(a)$ok)
  expect_identical(peptide_presenters(list(a$hla[[1]]))[[1]]$status,"UNRESOLVED")
  a<-first_ux_plan()$plan;a$hla<-NULL;expect_false(first_ux_check(a)$ok)
  a<-first_ux_plan()$plan;a$config$resources$reference.genome_fasta<-NULL;expect_false(first_ux_check(a)$ok)
  a<-first_ux_plan()$plan;a$config$predictors<-list(list(provider="netmhcpan",task="mhc1_binding",executable=file.path(tempdir(),"missing-provider")))
  expect_true(any(first_ux_check(a)$checks$observed=="EXECUTABLE_MISSING"))
  f<-tempfile();writeLines("existing file",f);a<-first_ux_plan()$plan
  expect_false(first_ux_check(a,output_dir=f)$ok)
})

test_that("VCF binding is checked without inventing GT or exposing identity", {
  a<-first_ux_plan()$plan;f<-tempfile(fileext=".vcf")
  lines<-readLines(a$inputs[[1]]$path);writeLines(sub("demo-sample","PRIVATE-OTHER",lines,fixed=TRUE),f);a$inputs[[1]]$path<-f
  z<-first_ux_check(a);expect_false(z$ok);expect_true(any(z$checks$observed=="SAMPLE_MISMATCH"))
  expect_false(grepl("PRIVATE-OTHER",jsonlite::toJSON(z,auto_unbox=TRUE)))
  # The authoritative parser accepts sites-only VCFs with explicit input binding.
  writeLines(c("##fileformat=VCFv4.2","#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"),f)
  z<-first_ux_check(a);expect_true(z$ok);expect_true(any(z$checks$observed=="SITES_ONLY_EXPLICIT_BINDING"))
  writeLines(c("##fileformat=VCFv4.2","#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tdemo-sample\tOTHER"),f)
  expect_false(first_ux_check(a)$ok)
  writeLines("not a VCF",f);expect_false(first_ux_check(a)$ok)
})

test_that("config diagnostics identify keys but not secret values or subject keys", {
  e<-tryCatch(neo_config(peptdies="PRIVATE-SECRET-VALUE"),error=identity)
  expect_s3_class(e,"error");expect_match(conditionMessage(e),"Unknown key in config: peptdies")
  expect_match(conditionMessage(e),"peptides");expect_false(grepl("PRIVATE-SECRET",conditionMessage(e)))
  expect_error(neo_config(runtime=list(threadz=1)),"Unknown key in runtime: threadz")
  e<-tryCatch(strict_keys(list("PRIVATE-PATIENT"=1),"EXPECTED-PATIENT","hla subject"),error=identity)
  expect_false(grepl("PRIVATE-PATIENT|EXPECTED-PATIENT",conditionMessage(e)))
  e<-tryCatch(strict_keys(stats::setNames(list(1),"ghp_PRIVATEVALUE"),"paths","config"),error=identity)
  expect_false(grepl("ghp_PRIVATEVALUE",conditionMessage(e)))
})

test_that("normal-route tutorial has standard lengths and honest unpredicted results", {
  x<-first_ux_plan();a<-x$plan;expect_true(first_ux_check(a)$ok)
  p<-neo_run(do.call(neo_project,a));r<-neo_results(p);o<-neo_origins(p);events<-neo_events(p)
  expect_equal(nrow(events),1L);expect_equal(nrow(o),65L);expect_equal(nrow(r),65L)
  expect_equal(as.integer(table(nchar(r$peptide_sequence))),c(8L,9L,10L,11L,12L,15L))
  expect_true(all(r$prediction_status=="NOT_RUN"));expect_true(all(r$grading_status=="NOT_GRADED"))
  expect_true(all(is.na(r$tumor_vaf)))
  before<-ia_files(p$path);reopened<-neo_project(path=p$path)
  expect_identical(neo_events(reopened),events);expect_identical(neo_origins(reopened),o);expect_identical(neo_results(reopened),r)
  report<-neo_report_bundle(reopened);expect_equal(report$manifest$original_project_graph_counts$candidates,65L)
  expect_equal(nrow(report$tables$candidate_origin_links),65L)
  expect_identical(ia_files(p$path),before)
  expect_gt(length(neo_explain(reopened,candidate_id=r$candidate_id[[1]])$candidates[[1]]$origins),0L)
})


test_that("missing declarations and privacy controls remain explicit", {
  a<-first_ux_plan()$plan
  expect_false(first_run_preflight(a$config)$ok)
  a$inputs<-list("invalid record")
  expect_false(first_ux_check(a)$ok)
  a<-first_ux_plan()$plan;a$samples$patient_id<-"PRIVATE_PATIENT_CANARY"
  names(a$hla)<-"PRIVATE_PATIENT_CANARY"
  z<-first_ux_check(a)
  expect_true(z$ok);expect_false(grepl("PRIVATE_PATIENT_CANARY",jsonlite::toJSON(z,auto_unbox=TRUE)))
  private<-first_ux_check(a,include_private=TRUE)
  expect_identical(private$private$samples,a$samples)
  expect_true(grepl("PRIVATE_PATIENT_CANARY",jsonlite::toJSON(private,auto_unbox=TRUE)))
})


test_that("unreadable VCF and actual importer mismatch keep actionable failures", {
  a<-first_ux_plan()$plan;f<-tempfile(fileext=".vcf");file.copy(a$inputs[[1]]$path,f)
  a$inputs[[1]]$path<-f;a$inputs[[1]]$checksum<-NULL
  Sys.chmod(f,"0000");denied<-file.access(f,4L)!=0L
  if(!denied){Sys.chmod(f,"0600");skip("Host does not deny owner read permission")}
  z<-tryCatch(first_ux_check(a),finally=Sys.chmod(f,"0600"))
  expect_false(z$ok);expect_true(any(z$checks$observed=="MISSING_OR_UNREADABLE"))
  lines<-readLines(f);writeLines(sub("demo-sample","PRIVATE_HEADER_ID",lines,fixed=TRUE),f)
  p<-do.call(neo_project,a)
  ref_path<-a$config$resources$reference.genome_fasta$path
  con<-file(ref_path,"rb")
  fai<-utils::read.delim(a$config$resources$reference.genome_fai$path,header=FALSE,
    col.names=c("contig","length","offset","line_bases","line_width"))
  e<-tryCatch(event_parse_vcf(p,p$inputs[[1]],list(con=con,fai=fai,assembly="GRCh38.p14")),error=identity,finally=close(con))
  expect_s3_class(e,"error");expect_match(conditionMessage(e),"VCF sample does not match")
  expect_match(conditionMessage(e),"#CHROM");expect_false(grepl("PRIVATE_HEADER_ID",conditionMessage(e)))
})
