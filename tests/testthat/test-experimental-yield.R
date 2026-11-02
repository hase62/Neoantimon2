test_that("A evidence retains positives, negatives, conflicts and all B subtypes", {
  x <- data.frame(peptide_sequence = "acdefghik", mhc_key = "A*02:01")
  records <- list(yield_test_record("a1"),yield_test_record("a2",result="negative"),
    yield_test_record("b1",hla="HLA-B*07:02",result="negative"),
    yield_test_record("b2",pep="ACDEYGHIK"))
  y <- yield_test_evidence(x,yield_test_reference(records))
  expect_identical(y$evidence_category,"A");expect_identical(y$A_status,"conflicting")
  expect_identical(y$B1_status,"negative");expect_identical(y$B2_status,"positive")
  expect_setequal(y$b_subtypes[[1]],c("B1","B2"));expect_length(y$matches[[1]],4)
  expect_identical(y$matches[[1]][[1]]$record$patient_id,"p1")
  expect_identical(y$source_peptide_sequence,"acdefghik")
})
test_that("B1 and B2 include negative/conflicting experimental results without ranking them", {
  x <- data.frame(peptide_sequence="ACDEFGHIK",mhc_key="HLA-A*02:01")
  r <- list(yield_test_record("b11",hla="HLA-B*07:02"),yield_test_record("b12",hla="HLA-B*07:02",result="negative"),
    yield_test_record("b21",pep="ACDEYGHIK"),yield_test_record("b22",pep="ACDEYGHIK",result="negative"))
  y <- yield_test_evidence(x,yield_test_reference(r));expect_identical(y$evidence_category,"B")
  expect_identical(y$B1_status,"conflicting");expect_identical(y$B2_status,"conflicting")
  expect_setequal(y$b_subtypes[[1]],c("B1","B2"));expect_length(y$matches[[1]],4)
  s <- y$matches[[1]][[3]]$substitutions
  expect_equal(s$position_1based,5);expect_identical(s$reference_residue,"Y");expect_identical(s$candidate_residue,"F")
  one <- yield_test_evidence(x,yield_test_reference(list(r[[2]])))
  expect_identical(one$evidence_category,"B");expect_identical(one$B1_status,"negative")
})
test_that("B2 requires exact HLA, equal length and exactly one substitution", {
  x <- data.frame(peptide_sequence="ACDEFGHIK",mhc_key="HLA-A*02:01")
  r <- list(yield_test_record("long",pep="ACDEFGHIKA"),yield_test_record("other",pep="ACDEYGHIK",hla="HLA-A*02:05"),
    yield_test_record("two",pep="ACDEYGHIA"))
  y <- yield_test_evidence(x,yield_test_reference(r));expect_identical(y$evidence_category,"D");expect_length(y$matches[[1]],0)
  r <- yield_test_reference(list(yield_test_record("I-to-L",pep="ACDEFGHLK")))
  expect_identical(yield_test_evidence(x,r)$B2_status,"positive")
})
test_that("conservative HLA punctuation normalization never guesses alleles or resolution", {
  h <- Neoantimon2:::yield_hla
  expect_identical(h(c("HLA-A02:01","a*02:01"," HLA-A*02:01 ")),rep("HLA-A*02:01",3))
  expect_identical(h("HLA-A*02:01:01:01"),"HLA-A*02:01:01:01")
  expect_identical(h("HLA-A*02:01N"),"HLA-A*02:01N")
  expect_error(h("A0201"),"notation");expect_error(h("A2"),"notation")
  expect_error(h("HLA-A*02:XX"),"notation")
  expect_identical(h("HLA-DPB1*04:01~HLA-DPA1*01:03"),"HLA-DPA1*01:03~HLA-DPB1*04:01")
  expect_error(h("HLA-DPA1*01:03~HLA-DQB1*03:01"),"Incompatible")
  x <- data.frame(peptide_sequence="ACDEFGHIK",mhc_key="HLA-A*02:01:01")
  y <- yield_test_evidence(x,yield_test_reference(list(yield_test_record("a"))))
  expect_false(y$evidence_category=="A")
  expect_error(neo_yield_evidence(data.frame(peptide_sequence="ACD X",mhc_key="HLA-A*02:01")),"amino")
})
test_that("C thresholds count independent tested units and patients conservatively", {
  x <- data.frame(peptide_sequence="YYYYYYYYY",mhc_key="HLA-A*02:01")
  seqs <- yield_test_scores(seq_len(100))$peptide_sequence
  rows <- lapply(seq_len(100),function(i)yield_test_record(paste0("r",i),pep=seqs[i],patient=paste0("p",i%%3)))
  rows[[101]] <- rows[[1]];rows[[101]]$record_id <- "technical-repeat";rows[[101]]$result <- "negative"
  ref <- yield_test_reference(rows);y <- yield_test_evidence(x,ref)
  expect_identical(y$evidence_category,"C");expect_equal(y$hla_tested_units,100);expect_equal(y$hla_verified_patients,3)
  expect_identical(yield_test_evidence(x,ref,min_records=101L)$evidence_category,"D")
  expect_identical(yield_test_evidence(x,ref,min_patients=4L)$evidence_category,"D")
  for(i in seq_along(rows)) {rows[[i]]$study_id <- paste0("study",i%%3);rows[[i]]$patient_id <- "local-p1"}
  y <- yield_test_evidence(x,yield_test_reference(rows));expect_equal(y$hla_verified_patients,1)
  expect_identical(y$evidence_category,"D")
  expect_error(yield_test_evidence(x,ref,min_records=0),"min_records")
})
test_that("unavailable references and incompatible assays do not create D or negative evidence", {
  x <- data.frame(peptide_sequence="ACDEFGHIK",mhc_key="HLA-A*02:01")
  y <- neo_yield_evidence(x);expect_true(is.na(y$evidence_category));expect_identical(y$evidence_status,"REFERENCE_UNAVAILABLE")
  y <- neo_yield_evidence(x,tempfile());expect_true(is.na(y$evidence_category))
  ref <- yield_test_reference(list(yield_test_record("wrong-assay",assay="ELISpot")))
  y <- yield_test_evidence(x,ref);expect_true(is.na(y$evidence_category));expect_identical(y$evidence_status,"REFERENCE_ENDPOINT_ASSAY_UNAVAILABLE")
  ref <- yield_test_reference(list(yield_test_record("right"),yield_test_record("wrong",result="negative",assay="ELISpot")))
  y <- yield_test_evidence(x,ref);expect_identical(y$A_status,"positive");expect_length(y$matches[[1]],1)
  expect_error(neo_yield_evidence(x,ref),"endpoint")
  ref$records[[1]]$result <- "untested";expect_error(yield_test_evidence(x,ref),"explicitly tested")
})
test_that("all four fixed bins and boundary ties respect actual direction", {
  x <- yield_test_scores(c(0,.3,.6,.9,1));y <- neo_experimental_yield(x,yield_test_profile(),ordering="input")
  expect_equal(y$units$score_bin,c(1,1,2,3,4));expect_equal(y$units$probability,c(.01,.01,.02,.03,.04))
  x$raw_score <- c(1000,900,600,100,0);x$score_direction <- "lower_is_better";x$unit <- "nM"
  p <- yield_test_profile(direction="lower_is_better");p$binning$thresholds <- c(900,600,100)
  y <- neo_experimental_yield(x,p,ordering="input");expect_equal(y$units$score_bin,c(1,1,2,3,4))
  p$binning$thresholds <- c(100,600,900);expect_error(neo_experimental_yield(x,p),"increasing prediction quality")
})
test_that("empirical profiles require matching model, metric, units and exact observed counts", {
  x <- yield_test_scores();p <- yield_test_profile();p$mode <- "empirical"
  p$model <- as.list(x[1,Neoantimon2:::yield_model_fields()]);p$bins <- lapply(1:4,function(i)list(bin=i,tested=10,positive=i,negative=10-i,rate=i/10))
  y <- neo_experimental_yield(x,p);expect_match(y$summary$probability_interpretation,"empirical_model_endpoint_matched")
  expect_equal(sort(y$units$probability),c(.1,.2,.3,.4))
  for(f in c("provider","model_version","predictor_version","metric","unit","metric_semantics","score_direction")) {
    z <- p;z$model[[f]] <- "wrong";expect_error(neo_experimental_yield(x,z),"mismatch")
  }
  z <- p;z$bins[[1]]$rate <- .11;expect_error(neo_experimental_yield(x,z),"observed counts")
  z <- p;z$bins[[1]]$tested <- 0;expect_error(neo_experimental_yield(x,z),"Inconsistent")
  z <- p;z$binning$basis <- "cohort_relative";expect_error(neo_experimental_yield(x,z),"frozen")
})
test_that("IMPROVE scenario is never mislabeled BigMHC calibration", {
  x <- yield_test_scores(seq(0,1,length.out=20));x$provider <- "bigmhc-im"
  y <- neo_experimental_yield(x,"improve-exploratory-v0.1")
  expect_identical(y$summary$mode,"scenario")
  expect_match(y$summary$probability_interpretation,"CROSS_MODEL_ASSUMPTION_NOT_CALIBRATED")
  expect_equal(sort(unique(y$units$probability)),c(.0139,.023,.0306,.0639))
  expect_match(paste(y$warnings,collapse=" "),"candidate-set-dependent")
})
test_that("repeated Origins never change pair bins, probabilities or synthesized counts", {
  x <- yield_test_scores();dup <- x[1,,drop=FALSE];dup$origin_id <- "additional-origin";dup$origin_ids <- I(list(list("additional-origin")))
  base <- neo_experimental_yield(x,yield_test_profile(),ordering="input")
  y <- neo_experimental_yield(rbind(x,dup),yield_test_profile(),ordering="input")
  expect_equal(y$units$probability,base$units$probability);expect_equal(y$curve,base$curve)
  expect_equal(y$summary$N_available,4);expect_equal(y$summary$distinct_synthesized_peptides,4)
  expect_setequal(y$pair_evidence$all_origin_ids[[1]],c("o1","additional-origin"))
  dup$raw_score <- .9;expect_error(neo_experimental_yield(rbind(x,dup),yield_test_profile()),"Conflicting duplicate")
})
test_that("same peptide with multiple HLA uses a labelled primary proxy and retains alternatives", {
  x <- yield_test_scores(c(.4,.8));x$peptide_sequence <- "ACDEFGHIK";x$mhc_key[2] <- "HLA-B*07:02"
  y <- neo_experimental_yield(x,yield_test_profile(c(.1,.2,.3,.4)))
  expect_equal(y$summary$N_available,1);expect_equal(y$summary$distinct_synthesized_peptides,1)
  expect_identical(y$units$mhc_key,"HLA-B*07:02");expect_equal(y$units$probability,.3)
  expect_equal(nrow(y$units$alternative_pairs[[1]]),2);expect_true(y$units$primary_HLA_proxy)
  expect_setequal(y$units$all_candidate_ids[[1]],x$candidate_id)
  expect_setequal(y$units$all_origin_ids[[1]],x$origin_id)
  y <- neo_experimental_yield(x,yield_test_profile(c(.1,.2,.3,.4)),primary_hla="lexicographic")
  expect_identical(y$units$mhc_key,"HLA-A*02:01");expect_equal(y$units$probability,.2)
  y <- neo_experimental_yield(x,yield_test_profile(c(.1,.2,.3,.4)),unit="peptide_hla")
  expect_equal(y$summary$N_available,2);expect_equal(y$summary$distinct_synthesized_peptides,1)
  expect_equal(y$summary$expected_positive_units,.5);expect_true(is.na(y$summary$probability_at_least_one))
  expect_true(all(is.na(y$targets$minimum_N)));expect_match(y$targets$reason[1],"JOINT_PROBABILITY_UNAVAILABLE")
})
test_that("A-D and B-subtype classifications cannot alter independent probabilities", {
  x <- yield_test_scores(c(.1,.4,.7,1));p <- yield_test_profile();base <- neo_experimental_yield(x,p)
  rr <- list(yield_test_record("positive",pep=x$peptide_sequence[1]),yield_test_record("negative",pep=x$peptide_sequence[2],result="negative"),yield_test_record("B1",pep=x$peptide_sequence[3],hla="HLA-B*07:02"))
  ref <- yield_test_reference(rr);y <- neo_experimental_yield(x,p,ref)
  expect_equal(y$curve,base$curve);expect_equal(y$targets,base$targets)
  expect_identical(y$units$probability,base$units$probability)
  expect_true(any(y$units$evidence_category=="A"));expect_true(any(y$units$evidence_category=="B"))
})
test_that("analytical yield, monotonicity and target minima are reproduced stably", {
  x <- yield_test_scores(rep(.5,25));y <- neo_experimental_yield(x,yield_test_profile(rep(.1,4)))
  N <- 0:25
  expect_equal(y$curve$expected_positive_units,.1*N,tolerance=1e-14)
  expect_equal(y$curve$probability_at_least_one,1-.9^N,tolerance=1e-14)
  expect_true(all(diff(y$curve$probability_at_least_one)>=0))
  expect_equal(y$targets$minimum_N,c(7,16,22,NA));expect_identical(y$targets$reachable,c(TRUE,TRUE,TRUE,FALSE))
  expect_match(y$targets$reason[4],"UNREACHABLE")
  y <- neo_experimental_yield(x[1:3,],yield_test_profile(rep(.01,4)))
  expect_true(all(is.na(y$targets$minimum_N)));expect_true(all(!y$targets$reachable))
  y <- neo_experimental_yield(x,yield_test_profile(rep(1e-16,4)))
  expect_gt(y$summary$probability_at_least_one,0)
  expect_lt(abs(y$summary$probability_at_least_one - (-expm1(25*log1p(-1e-16)))),1e-28)
})
test_that("missing predictions are unknown, genuine zero probabilities and p=1 remain valid", {
  x <- yield_test_scores(c(0,NA,1));x$prediction_status[2] <- "NOT_RUN"
  y <- neo_experimental_yield(x,yield_test_profile(c(0,.2,.3,1)),ordering="input")
  expect_equal(y$units$raw_score[1],0);expect_equal(y$units$probability[1],0)
  expect_true(is.na(y$units$probability[2]));expect_equal(y$units$probability[3],1)
  expect_true(is.na(y$summary$expected_positive_units));expect_equal(y$summary$N_unscored,1)
  expect_equal(y$summary$expected_positive_scored_subset,1);expect_equal(y$summary$probability_scored_subset,1)
  expect_true(all(is.na(y$curve$probability_at_least_one[3:4])))
  z <- neo_experimental_yield(x[c(1,3),],yield_test_profile(c(0,.2,.3,1)))
  expect_equal(z$summary$probability_at_least_one,1);expect_equal(z$targets$minimum_N,rep(1,4))
  z <- neo_experimental_yield(x[1,,drop=FALSE],yield_test_profile(rep(0,4)))
  expect_equal(z$summary$probability_at_least_one,0);expect_true(all(!z$targets$reachable))
})
test_that("invalid rates, scores, populations and identity conflicts are rejected", {
  x <- yield_test_scores();p <- yield_test_profile()
  for(r in list(c(NA,.2,.3,.4),c(-.1,.2,.3,.4),c(1.1,.2,.3,.4),c(.1,.2,.3),c(Inf,.2,.3,.4))) {
    z <- p;z$rates <- r;expect_error(neo_experimental_yield(x,z),"four valid probabilities")
  }
  y <- x;y$raw_score[1] <- NaN;expect_error(neo_experimental_yield(y,p),"raw_score")
  y <- x;y$raw_score[1] <- Inf;expect_error(neo_experimental_yield(y,p),"raw_score")
  y <- x;y$sample_id[1] <- "s2";expect_error(neo_experimental_yield(y,p),"one sample")
  y <- x;y$run_id[1] <- "other";expect_error(neo_experimental_yield(y,p),"exactly one")
  expect_equal(nrow(neo_experimental_yield(y,p,channel=list(run_id="test-run"))$units),3)
  y <- x;y$raw_score[1] <- -.1;y$unit <- "nM";expect_error(neo_experimental_yield(y,p),"Negative affinity")
})
test_that("evidence-only lookup and JSON profile/reference are reproducible", {
  x <- yield_test_scores();ref <- yield_test_reference(list(yield_test_record("a",pep=x$peptide_sequence[1])))
  path <- tempfile(fileext=".json");pp <- tempfile(fileext=".json")
  jsonlite::write_json(ref,path,auto_unbox=TRUE);jsonlite::write_json(yield_test_profile(),pp,auto_unbox=TRUE)
  y <- neo_experimental_yield(x,pp,path)
  expect_identical(y$profile$artifact$path,normalizePath(pp));expect_identical(y$reference$data$id,"synthetic-reference")
  expect_identical(y$reference$artifact$sha256,digest::digest(file=path,algo="sha256",serialize=FALSE))
  expect_identical(neo_yield_evidence(x,path,"synthetic-T-cell","synthetic-assay")$evidence_category,
    yield_test_evidence(x,ref)$evidence_category)
})
test_that("yield access never changes supplied native grades/ranks or external semantics", {
  x <- yield_test_scores();before <- serialize(x,NULL);y <- neo_experimental_yield(x,yield_test_profile())
  expect_identical(serialize(x,NULL),before);expect_identical(y$input_provenance,x)
  for(i in seq_len(nrow(y$pair_evidence))) {
    src <- y$pair_evidence$source_rows[[i]];at <- match(src$candidate_id,x$candidate_id)
    expect_identical(src$native_grade,x$native_grade[at]);expect_identical(src$native_rank,x$native_rank[at])
    expect_identical(src$raw_score,x$raw_score[at]);expect_identical(src$score_direction,x$score_direction[at])
  }
  y <- x;y$predictor_version <- "";y$model_version <- ""
  expect_s3_class(neo_experimental_yield(y,yield_test_profile()),"neo_experimental_yield")
})

test_that("direct native/external project yield is read-only and preserves Grades and provenance", {
  base <- reconstruction_test_project();cfg <- base$config
  cfg$peptides <- Neoantimon2:::peptide_policy(list(class_i_lengths=8:9,class_ii_lengths=15L))
  exe <- tempfile("explicit-synthetic-netmhcpan-")
  writeLines(c("#!/bin/sh", "if [ \"$1\" = '-h' ]; then echo '# NetMHCpan version 4.2'; exit 0; fi",
    "input=''; allele=''", "while [ $# -gt 0 ]; do case \"$1\" in -p|-f) input=$2; shift 2;; -a) allele=$2; shift 2;; *) shift;; esac; done",
    "echo 'Pos HLA Peptide Score_EL %Rank_EL'", "n=0; while IFS= read -r peptide; do n=$((n+1)); printf '%s %s %s 0.5 0.5\\n' \"$n\" \"$allele\" \"$peptide\"; done < \"$input\""),exe)
  Sys.chmod(exe,"0755")
  cfg$predictors <- list(list(provider="netmhcpan",task="mhc1_presentation",executable=exe,options=list()))
  vcf <- tempfile(fileext=".vcf")
  writeLines(c("##fileformat=VCFv4.2","##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample-1",
    "chrW\t4\t.\tA\tG\t.\tPASS\t.\tGT\t0/1"),vcf)
  p <- neo_run(neo_project(samples=base$samples,config=cfg,hla=list("patient-1"="HLA-A*02:01"),
    inputs=list(list(input_id="yield-synthetic",input_role="somatic_vcf",sample_id="sample-1",path=vcf))))
  before <- neo_results(p);z <- neo_priorities(p);expect_true(any(z$prediction_status=="evaluated"))
  hashes <- function() {
    files <- sort(list.files(p$path,recursive=TRUE,full.names=TRUE,all.files=TRUE));files <- files[!dir.exists(files)]
    stats::setNames(vapply(files,function(f)digest::digest(file=f,algo="sha256",serialize=FALSE),""),files)
  }
  frozen <- hashes();metric <- sort(unique(z$metric[!is.na(z$metric)]))[1]
  scenario <- yield_test_profile();scenario$binning <- list(basis="cohort_relative")
  plan <- neo_experimental_yield(p,scenario,channel=list(provider="netmhcpan",metric=metric))
  expect_s3_class(plan,"neo_experimental_yield");expect_identical(neo_results(p),before);expect_identical(hashes(),frozen)
  expect_setequal(unlist(plan$pair_evidence$all_origin_ids),unique(z$origin_id[z$metric==metric]))
  q <- before[1,,drop=FALSE];pred <- external_test_row(q$peptide_sequence,q$mhc_key,
    task="external:immunogenicity",provider="yield-synthetic-external")
  a <- neo_analyze_external(p,external_test_write(pred),list(mode="external_rank",
    provider="yield-synthetic-external",run_id="external-run-1",task="external:immunogenicity",metric="score",score_direction="higher_is_better"))
  ext <- neo_external_results(a);frozen <- hashes()
  plan2 <- neo_experimental_yield(p,scenario,analysis=a,channel=list(provider="yield-synthetic-external"))
  expect_identical(neo_results(p),before);expect_identical(neo_external_results(a),ext);expect_identical(hashes(),frozen)
  expect_setequal(unlist(plan2$pair_evidence$all_origin_ids),unique(neo_priorities(p,analysis=a)$origin_id))
  expect_identical(neo_results(neo_project(path=p$path)),before)
})


test_that("candidate and patient identities cannot be relabelled across observations", {
  x <- yield_test_scores();y <- x;y$candidate_id[2] <- y$candidate_id[1]
  expect_error(neo_experimental_yield(y,yield_test_profile()),"candidate_id maps to conflicting")
  y <- x;y$patient_id <- c("P1","P2","P1","P1")
  expect_error(neo_experimental_yield(y,yield_test_profile()),"patient identities")
})


test_that("experimental evidence categories never overwrite BA/EL/IM task category", {
  x <- yield_test_scores();x$category <- "IM"
  ref <- yield_test_reference(list(yield_test_record("a",pep=x$peptide_sequence[1])))
  y <- neo_experimental_yield(x,yield_test_profile(),reference=ref)
  expect_true(all(y$units$category=="IM"));expect_true(all(y$pair_evidence$category=="IM"))
  expect_true(any(y$units$evidence_category=="A"));expect_identical(y$input_provenance,x)
})
