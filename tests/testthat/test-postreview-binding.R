test_that("paired configured BA queries are explicit, typed and preserve native results", {
  source(system.file("examples","paper-binding.R",package="Neoantimon2"),local=TRUE)
  fixture <- postreview_binding_project();p <- fixture$p
  before <- neo_results(p);files <- list.files(p$path,recursive=TRUE,full.names=TRUE,all.files=TRUE)
  files <- files[!dir.exists(files)];hash <- function(x)vapply(x,function(f)digest::digest(file=f,algo="sha256",serialize=FALSE),"")
  native_hash <- hash(files);inputs_hash <- hash(fixture$inputs)
  output <- tempfile("paired-output-");calls <- length(readLines(fixture$calls))
  plan <- run_paper_binding_queries(p,output)
  expect_true(plan$available);expect_gt(plan$query_pairs,0)
  expect_false(dir.exists(output));expect_equal(length(readLines(fixture$calls)),calls)
  expect_error(run_paper_binding_queries(p,file.path(p$path,"forbidden"),execute=TRUE),"overlap")
  expect_false(dir.exists(file.path(p$path,"forbidden")))
  paired <- run_paper_binding_queries(p,output,execute=TRUE)
  expect_gt(length(readLines(fixture$calls)),calls)
  wt <- Neoantimon2:::external_read_predictions(paired$wild_type)
  expect_true(all(wt$metric=="BA_nM"));expect_true(all(wt$unit=="nM"))
  expect_true(all(wt$score_direction=="lower_is_better"))
  expect_true(file.exists(file.path(output,"completion.json")))
  expect_identical(hash(files),native_hash);expect_identical(hash(fixture$inputs),inputs_hash)
  expect_identical(neo_results(p),before)
  expect_error(run_paper_binding_queries(p,output,execute=TRUE),"new independent")
  a <- neo_analyze_external(p,paired$mutant,policy=list(mode="external_rank",provider=plan$provider,
    run_id=paired$run_id,task=plan$task,metric=plan$metric,score_direction="lower_is_better"))
  score <- neo_priorities(p,analysis=a,method="neoantimon-paper-2020-v1",
    wild_type_predictions=paired$wild_type,expression_unit="TPM")
  expect_true(all(is.finite(score$paper_with_ccf)))
  expect_true(all(score$dna_vaf==.2));expect_true(all(score$gene_tpm==3))
  expect_true(all(score$paper_equation=="P_I"))
  expect_true(all(score$adjusted_priority_status=="exploratory_candidate_set_dependent"))
  expect_identical(neo_results(p),before)
  expect_gt(nrow(neo_external_results(a)),0)
  expect_true(length(neo_external_explain(a,candidate_id=before$candidate_id[[1]]))>0)
  # P_R uses a separately selected real binding-percentile output, not an EL score.
  paired_rank <- run_paper_binding_queries(p,tempfile("paired-rank-"),metric="BA_rank",execute=TRUE)
  rank <- Neoantimon2:::external_read_predictions(paired_rank$wild_type)
  expect_true(all(rank$metric_semantics=="mhc_binding_percentile"))
  expect_true(all(rank$unit=="percent"))
  ar <- neo_analyze_external(p,paired_rank$mutant,policy=list(mode="external_rank",provider=plan$provider,
    run_id=paired_rank$run_id,task=plan$task,metric="BA_rank",score_direction="lower_is_better"))
  pr <- neo_priorities(p,analysis=ar,method="neoantimon-paper-2020-v1",
    wild_type_predictions=paired_rank$wild_type,expression_unit="TPM")
  expect_true(all(is.finite(pr$paper_with_ccf)));expect_true(all(pr$paper_equation=="P_R"))
  expect_identical(neo_results(p),before)
  no_provider <- p;no_provider$config$predictors <- list()
  expect_error(run_paper_binding_queries(no_provider,tempfile()),"configured native")
  missing <- p;missing$config$predictors[[1]]$executable <- file.path(tempdir(),"absent-netMHCpan")
  expect_false(run_paper_binding_queries(missing,tempfile())$available)
  expect_error(run_paper_binding_queries(missing,tempfile(),execute=TRUE),"prerequisites unavailable")
  # An engine returning EL only cannot become fabricated BA/WT scores.
  writeLines(c("#!/bin/sh","if [ \"$1\" = '-h' ];then echo 'NetMHCpan version 4.2';exit 0;fi",
    "echo 'Pos HLA Peptide Score_EL %Rank_EL'"),fixture$tool)
  bad_output <- tempfile("paired-failed-")
  expect_error(run_paper_binding_queries(p,bad_output,execute=TRUE),"unambiguous requested BA metric")
  expect_false(file.exists(file.path(bad_output,"completion.json")))
  expect_true(dir.exists(file.path(bad_output,"intermediate","prediction_runs")))
  expect_identical(neo_results(p),before)
  # Reuse the stub's own parser but make the returned affinity impossible.
  original <- c("#!/bin/sh", "if [ \"$1\" = '-h' ]; then echo 'NetMHCpan version 4.2'; exit 0; fi",
    "input='';allele='';while [ $# -gt 0 ];do case \"$1\" in -p) input=$2;shift 2;; -a) allele=$2;shift 2;; *) shift;; esac;done",
    "echo 'Pos HLA Peptide Aff(nM) %Rank_BA';n=0;while IFS= read -r q;do n=$((n+1));printf '%s %s %s -1 0.5\\n' \"$n\" \"$allele\" \"$q\";done < \"$input\"")
  writeLines(original,fixture$tool)
  bad_negative <- tempfile("paired-negative-")
  expect_error(run_paper_binding_queries(p,bad_negative,execute=TRUE),"negative BA")
  expect_false(file.exists(file.path(bad_negative,"completion.json")))
})

test_that("public numerical view exposes conflicting expression without a maximum", {
  fixture <- postreview_binding_project(TRUE);p <- fixture$p;before <- neo_results(p)
  score <- neo_priorities(p)
  expect_true(all(is.na(score$gene_tpm)));expect_true(all(score$gene_tpm_status=="conflict"))
  expect_true(all(vapply(score$expression_evidence,function(rows)
    setequal(unlist(lapply(rows,`[[`,"gene_tpm")),c(3,9)),TRUE)))
  paper <- neo_priorities(p,method="neoantimon-paper-2020-v1",expression_unit="TPM")
  expect_true(all(is.na(paper$paper_with_ccf)))
  expect_true(all(paper$paper_status=="expression_conflict_no_paper_score"))
  expect_identical(neo_results(p),before)
})

test_that("paired query provenance remains sample and Origin specific after deduplication", {
  source(system.file("examples","paper-binding.R",package="Neoantimon2"),local=TRUE)
  fixture <- postreview_binding_project();old <- fixture$p
  vcf1 <- old$inputs[[1L]]$resolved_path;vcf2 <- tempfile(fileext=".vcf")
  writeLines(sub("\tsample-1$","\tsample-2",readLines(vcf1)),vcf2)
  cfg <- old$config;cfg$paths$run_root <- tempfile("paired-cohort-runs-")
  p <- neo_run(neo_project(samples=data.frame(sample_id=c("sample-1","sample-2"),patient_id="patient-1"),
    hla=list("patient-1"=c("HLA-A*02:01","HLA-B*07:02")),config=cfg,
    inputs=list(list(input_id="S1",input_role="somatic_vcf",sample_id="sample-1",path=vcf1),
      list(input_id="S2",input_role="somatic_vcf",sample_id="sample-2",path=vcf2))))
  before <- neo_results(p)
  paired <- run_paper_binding_queries(p,tempfile("paired-cohort-"),execute=TRUE)
  map <- readRDS(file.path(paired$output,"query-origin-mapping.rds"))
  raw <- Neoantimon2:::predictor_read_run(file.path(paired$output,"intermediate","prediction_runs"),paired$run_id)
  expect_setequal(vapply(map,`[[`,"","sample_id"),c("sample-1","sample-2"))
  for(e in raw$evidence) {
    from <- Filter(function(x)x$query_id==e$candidate_id,map)
    expect_setequal(unlist(e$biological_peptide_ids),vapply(from,`[[`,"","origin_id"))
    expect_setequal(e$sample_id,vapply(from,`[[`,"","sample_id"))
    expect_true(file.exists(file.path(paired$output,e$raw_output_ref)))
  }
  mutant <- Neoantimon2:::external_read_predictions(paired$mutant)
  expected <- unique(vapply(Filter(function(x)x$kind=="mutant",map),function(x)
    paste(x$peptide_sequence,x$mhc_key),""))
  expect_equal(nrow(mutant),length(expected))
  expect_identical(neo_results(p),before)
})
