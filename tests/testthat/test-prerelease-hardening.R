hardening_hash <- function(root) {
  files <- sort(list.files(root,recursive=TRUE,full.names=TRUE,all.files=TRUE))
  files <- files[!dir.exists(files)]
  stats::setNames(vapply(files,function(f)digest::digest(file=f,algo="sha256",serialize=FALSE),""),files)
}

test_that("exact paper pairs require same nonempty run and model semantics", {
  source(system.file("examples","paper-binding.R",package="Neoantimon2"),local=TRUE)
  f <- postreview_binding_project();p <- f$p;native <- neo_results(p)
  pair <- run_paper_binding_queries(p,tempfile("hardening-pair-"),execute=TRUE)
  a <- neo_analyze_external(p,pair$mutant,policy=list(mode="external_rank",provider=pair$plan$provider,
    run_id=pair$run_id,task=pair$plan$task,metric="BA_nM",score_direction="lower_is_better"))
  table <- Neoantimon2:::external_read_predictions(pair$wild_type)
  score <- function(t=table) {
    path <- tempfile(fileext=".tsv")
    write.table(t,path,sep="\t",quote=FALSE,row.names=FALSE,na="")
    neo_priorities(p,analysis=a,method="neoantimon-paper-2020-v1",
      wild_type_predictions=path,expression_unit="TPM")
  }
  good <- score();expect_true(all(is.finite(good$paper_with_ccf)))
  expect_true(all(vapply(good$paper_wt_provenance,function(x)
    x$pairing=="exact_run_and_model_identity" && x$paired_run_id==pair$run_id,TRUE)))
  bad <- table;bad$run_id <- paste0(bad$run_id,"-different")
  x <- score(bad);expect_true(all(is.na(x$paper_without_ccf)))
  expect_true(all(x$paper_status=="WT_run_identity_mismatch"))
  bad <- table;bad$run_id <- "";expect_error(score(bad),"Required fields")
  bad <- table;bad$run_id <- NULL;expect_error(score(bad),"columns")
  for (field in c("model_version","predictor_version","metric","unit")) {
    bad <- table;bad[[field]] <- paste0(bad[[field]],"-other")
    x <- score(bad);expect_true(all(is.na(x$paper_without_ccf)))
  }
  # Unrelated runs do not make the one exact match ambiguous.
  extra <- table;extra$run_id <- paste0(extra$run_id,"-other")
  expect_identical(score(rbind(table,extra))$paper_without_ccf,good$paper_without_ccf)
  expect_error(score(rbind(table,table)),"DUPLICATE_PREDICTION")
  extra <- table;extra$score <- "999"
  expect_error(score(rbind(table,extra)),"DUPLICATE_PREDICTION")
  expect_identical(neo_results(p),native)
})

test_that("policy identity canonicalizes only validated numeric storage", {
  p <- Neoantimon2:::somatic_priority_policy()
  q <- jsonlite::fromJSON(jsonlite::toJSON(p,auto_unbox=TRUE,null="null"),simplifyVector=FALSE)
  equal <- Neoantimon2:::somatic_priority_policy_equal
  expect_false(identical(p,q));expect_true(equal(p,q))
  q$thresholds$I$binder <- 2.0000000001;expect_false(equal(p,q))
  q <- p;q$schema_version <- "unknown";expect_error(equal(p,q),"schema")
  q <- p;q$policy_version <- "unknown";expect_error(equal(p,q),"version")
  q <- p;q$extra <- 1;expect_error(equal(p,q),"Incomplete")
  q <- p;q$thresholds$I$binder <- "2";expect_error(equal(p,q),"threshold")
  q <- p;q$primary$I$run_id <- "other";expect_false(equal(p,q))
  q <- p;q$thresholds$I$binder <- NULL;expect_error(equal(p,q),"Incomplete")
  expect_error(equal(p,NULL),"Incomplete")
})

test_that("reopening rich predicted projects reads without any filesystem changes", {
  f <- postreview_binding_project();p <- f$p
  before <- neo_results(p);events <- neo_events(p);origins <- neo_origins(p)
  priorities <- neo_priorities(p);id <- before$candidate_id[[1L]]
  explain <- neo_explain(p,candidate_id=id)
  original <- hardening_hash(p$path);calls <- readLines(f$calls)
  q <- neo_project(path=p$path)
  expect_identical(neo_results(q),before)
  expect_identical(neo_events(q),events);expect_identical(neo_origins(q),origins)
  expect_identical(neo_explain(q,candidate_id=id),explain)
  expect_identical(neo_priorities(q),priorities)
  expect_identical(hardening_hash(p$path),original);expect_identical(readLines(f$calls),calls)
  # A real policy change still needs its separately materialized result.
  q$config$somatic_priority$thresholds$I$strong <- .4
  expect_error(neo_results(q),"not materialized")
  expect_identical(hardening_hash(p$path),original)
  q <- neo_project(path=p$path)
  q$config$somatic_priority$policy_version <- "corrupt"
  expect_error(neo_results(q),"version")
  expect_identical(hardening_hash(p$path),original)
})

hardening_ascat_fixture <- function() {
  p <- reconstruction_test_project();root <- tempfile("hardening-ascat-");dir.create(root)
  segments <- data.frame(sample=c("upstream-A","upstream-B"),chr="chr1",startpos=1,endpos=33,nMajor=c(1,3),nMinor=c(1,1))
  ascat.output <- list(segments=segments,purity=c("upstream-A"=.4,"upstream-B"=.8),ploidy=c("upstream-A"=2,"upstream-B"=4))
  workspace <- file.path(root,"ASCAT_objects.Rdata");save(ascat.output,file=workspace)
  path <- file.path(root,"segments.txt");write.table(segments,path,sep="\t",quote=FALSE,row.names=FALSE)
  source <- somatic_test_source("ascat");source$path <- path;source$format <- "ascat_segments"
  source$reference_build <- "GRCh38.p14"
  source$options <- list(coordinate_system="1-based-inclusive",sample_name="upstream-A",
    metadata_path=workspace,metadata_format="ascat_rdata")
  list(p=p,source=source,object=ascat.output,root=root,workspace=workspace)
}

test_that("documented ASCAT workspaces and companion fits preserve identity and integrity", {
  f <- hardening_ascat_fixture();p <- f$p
  bind <- function(s=f$source) {
    p$config$somatic_evidence <- list(sources=list(s))
    Neoantimon2:::somatic_source_files(p)[[1L]]
  }
  s <- bind();x <- Neoantimon2:::input_ascat(s,p)
  expect_length(x$rows,1);expect_equal(x$rows[[1]]$tcn.em,"2")
  expect_equal(x$rows[[1]]$purity,"0.4");expect_equal(x$rows[[1]]$ploidy,"2")
  expect_identical(x$sample_mapping$upstream_sample,"upstream-A")
  prov <- Neoantimon2:::somatic_source_provenance(x)
  expect_identical(prov$companion$sha256,digest::digest(file=f$workspace,algo="sha256",serialize=FALSE))
  bad <- f$source;bad$options$metadata_checksum <- paste(rep("0",64),collapse="")
  expect_error(bind(bad),"checksum")
  bad <- f$source;bad$options$metadata_format <- NULL;expect_error(bind(bad),"metadata_format")
  bad <- f$source;bad$options$sample_name <- NULL
  expect_error(Neoantimon2:::input_ascat(bind(bad),p),"sample identities")
  # Direct workspace import needs no hand-built conversion file.
  direct <- f$source;direct$path <- f$workspace;direct$format <- "ascat_rdata"
  direct$options$metadata_path <- direct$options$metadata_format <- NULL
  expect_equal(Neoantimon2:::input_ascat(bind(direct),p)$rows,x$rows)
  # Missing purity/ploidy stay missing; no CN-derived estimates.
  ascat.output <- f$object;ascat.output$purity[] <- NA_real_;ascat.output$ploidy[] <- NA_real_
  save(ascat.output,file=f$workspace)
  expect_null(Neoantimon2:::input_ascat(bind(),p)$rows[[1]]$purity)
  expect_null(Neoantimon2:::input_ascat(bind(),p)$rows[[1]]$ploidy)
  # A different segmentation cannot lend its metadata to the declared fit.
  ascat.output$segments$nMajor[[1]] <- 2;save(ascat.output,file=f$workspace)
  expect_error(Neoantimon2:::input_ascat(bind(),p),"do not match")
  expect_false(identical(s$source_id,bind()$source_id))
  wrong <- 1;save(wrong,file=f$workspace)
  expect_error(Neoantimon2:::input_ascat(bind(direct),p),"exact object")
})


test_that("minimal RNA ASCAT and numerical-view examples reopen read-only", {
  source(system.file("examples","reopen-evidence.R",package="Neoantimon2"),local=TRUE)
  output <- tempfile("hardening-reopen-")
  demo <- run_reopen_evidence_demo(output)
  verified <- jsonlite::read_json(file.path(output,"reopen-verification.json"))
  expect_true(all(vapply(verified,function(x)isTRUE(x$read_only),TRUE)))
  expect_true(all(c("A","C","E") %in% names(verified)))
  expect_true(any(is.finite(demo$views$D_paper$paper_with_ccf)))
  expect_true(all(is.na(demo$views$E$dna_vaf)))
  expect_true(all(demo$views$E$evidence_adjusted==0))
  # Re-declare the same synthetic native ASCAT result as its documented workspace.
  old <- demo$C;cfg <- old$config
  ix <- which(vapply(cfg$somatic_evidence$sources,function(s)s$kind=="ascat",TRUE))
  src <- cfg$somatic_evidence$sources[[ix]]
  ascat.output <- readRDS(src$path)
  workspace <- file.path(output,"ASCAT_objects.Rdata");save(ascat.output,file=workspace)
  src$path <- workspace;src$format <- "ascat_rdata"
  cfg$somatic_evidence$sources[[ix]] <- src
  q <- neo_run(neo_project(samples=old$samples,hla=list(P1=c("HLA-A*02:01","HLA-B*07:02")),
    config=cfg,inputs=list(list(input_id="workspace-VCF",input_role="somatic_vcf",
      sample_id="sample-1",path=old$inputs[[1L]]$resolved_path)),project_id="native-ASCAT-workspace"))
  result <- neo_results(q);numeric <- neo_priorities(q);h <- hardening_hash(q$path)
  reopened <- neo_project(path=q$path)
  expect_identical(neo_results(reopened),result)
  expect_identical(neo_priorities(reopened),numeric)
  expect_true(all(numeric$total_cn==2));expect_true(all(numeric$purity==.4))
  expect_identical(hardening_hash(q$path),h)
})
