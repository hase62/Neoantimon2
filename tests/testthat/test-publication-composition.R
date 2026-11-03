# Independent small synthetic fixtures; never model-execution evidence.
composition_script <- system.file("examples", "publication-planning.R", package="Neoantimon2")
if (!nzchar(composition_script)) composition_script <- testthat::test_path("..","..","inst","examples","publication-planning.R")
source(composition_script, local=TRUE)
composition_fixture <- function() {
  a <- yield_test_scores(c(.1,.8,.5,.8));a$hla_class <- "I"
  a$peptide_sequence <- c("ACDEFGHIK","ACDEFGHIK","ACDEFGHIL","ACDEFGHIM")
  a$mhc_key[2] <- "HLA-B*07:02";a$event_ids<-I(list(c("e1","e2"),c("e1","e2"),"e3","e4"))
  b<-a;b$provider<-"model-2";b$score_direction<-"lower_is_better";b$raw_score<-c(0,5,2,1)
  x<-rbind(a,b);x$prediction_scope<-"peptide_hla";x
}
composition_profile <- function(rates=rep(.1,4)) {
  p<-yield_test_profile(rates);p$binning<-list(basis="cohort_relative");p
}
composition_channels <- list(high=list(provider="test-model"),low=list(provider="model-2"))
composition_compare <- function(x=composition_fixture(),N=c(1,2,9),...) publication_model_selection(
  x,composition_channels,N,composition_profile(),...)

test_that("existing proxy selection handles alternative HLAs, ties and different directions",{
  y<-composition_compare()
  high<-subset(y$selected,model=="high" & requested_N==1)
  low<-subset(y$selected,model=="low" & requested_N==1)
  expect_identical(high$peptide_sequence,"ACDEFGHIK");expect_identical(high$mhc_key,"HLA-B*07:02")
  expect_identical(low$peptide_sequence,"ACDEFGHIK");expect_identical(low$mhc_key,"HLA-A*02:01")
  expect_equal(nrow(high$alternative_pairs[[1]]),2)
  expect_setequal(high$all_origin_ids[[1]],c("o1","o2"))
  expect_setequal(high$all_candidate_ids[[1]],c("c1","c2"))
  expect_equal(subset(y$summary,requested_N==9)$selected_peptides,c(3,3))
  expect_true(all(subset(y$summary,requested_N==9)$status=="FEWER_THAN_REQUESTED"))
  expect_equal(subset(y$pairwise,requested_N==1)$jaccard,1)
  # Independent handwritten sets at N=2: high K,M; low K,M.
  expect_setequal(subset(y$selected,model=="high" & requested_N==2)$peptide_sequence,c("ACDEFGHIK","ACDEFGHIM"))
  expect_equal(subset(y$pairwise,requested_N==2)$intersection,2)
})
test_that("Origins and compound Event sets remain source records without changing selection",{
  x<-composition_fixture();copy<-x[1,,drop=FALSE];copy$origin_id<-"o-extra";copy$origin_ids<-I(list(list("o-extra")))
  y<-composition_compare(rbind(x,copy));z<-composition_compare(x)
  expect_equal(y$summary,z$summary);expect_equal(y$pairwise,z$pairwise)
  expect_true("o-extra" %in% y$selected$all_origin_ids[[1]])
  rows<-y$selected$alternative_pairs[[1]]$source_rows
  expect_true(any(vapply(rows,function(r) any(vapply(r$event_ids,function(e) setequal(e,c("e1","e2")),TRUE)),TRUE)))
  expect_equal(y$input_records,rbind(x,copy))
})
test_that("different coverage, absent and missing predictions are not poor scores",{
  x<-composition_fixture();x<-x[!(x$provider=="model-2" & x$candidate_id=="c4"),]
  x$raw_score[x$provider=="model-2" & x$candidate_id=="c3"]<-NA_real_
  x$prediction_status[x$provider=="model-2" & x$candidate_id=="c3"]<-"failed"
  y<-composition_compare(x,pool="common_evaluated")
  expect_setequal(y$common_candidate_ids,c("c1","c2"));expect_true(all(y$summary$available_peptides==1))
  expect_identical(subset(y$coverage,model=="low" & candidate_id=="c4")$status,"no_record")
  expect_identical(subset(y$coverage,model=="low" & candidate_id=="c3")$status,"recorded_unavailable")
  z<-composition_compare(x,pool="model_available")
  expect_equal(subset(z$summary,requested_N==9)$available_peptides,c(3,1))
  expect_equal(subset(z$summary,model=="low")$selection_coverage,rep(.5,3))
  expect_true(all(!is.na(z$selected$raw_score)))
})
test_that("empty common pool and zero N do not invent overlap or extra peptides",{
  x<-composition_fixture();x$prediction_status[x$provider=="model-2"]<-"unsupported"
  x$raw_score[x$provider=="model-2"]<-NA_real_
  y<-composition_compare(x,N=c(0,2));expect_equal(nrow(y$selected),0)
  expect_true(all(y$summary$selected_peptides==0));expect_true(all(is.na(y$pairwise$jaccard)))
})
test_that("sample, class, scope, identity and channel ambiguity fail explicitly",{
  x<-composition_fixture();x$sample_id[1]<-"different";expect_error(composition_compare(x),"one sample")
  x<-composition_fixture();x$hla_class[1]<-"II";expect_error(composition_compare(x),"MHC class")
  x<-composition_fixture();x$prediction_scope<-NULL;expect_error(composition_compare(x),"prediction_scope")
  x<-composition_fixture();x$prediction_scope[1]<-"peptide";expect_error(composition_compare(x),"peptide-only")
  x<-composition_fixture();x$mhc_key[1]<-"HLA-A*02:05";expect_error(composition_compare(x),"conflicting identity")
  x<-composition_fixture();x$run_id[1]<-"other-run";expect_error(composition_compare(x),"exactly one")
  x<-composition_fixture();x$analysis_id<-"one";x$analysis_id[1]<-"two";expect_error(composition_compare(x),"one analysis")
  expect_error(composition_compare(N=NA_real_),"integer budgets")
  expect_error(publication_model_selection(composition_fixture(),list(a=composition_channels$high,b=composition_channels$high),1,composition_profile()),"identical channels")
})
test_that("cross-task comparison needs explicit opt-in and cannot imply accuracy",{
  x<-composition_fixture();x$task[x$provider=="model-2"]<-"external:binding"
  expect_error(composition_compare(x),"Different tasks")
  y<-composition_compare(x,allow_different_tasks=TRUE)
  expect_true(all(y$pairwise$comparison_interpretation=="MODEL_AND_TASK_CHOICE_NOT_ACCURACY"))
})
test_that("probability assumptions never change model selection or native fields",{
  x<-composition_fixture();before<-serialize(x,NULL)
  a<-publication_model_selection(x,composition_channels,2,composition_profile(rep(0,4)))
  b<-publication_model_selection(x,composition_channels,2,composition_profile(rep(1,4)))
  expect_equal(a$selected[c("peptide_sequence","mhc_key","all_origin_ids")],b$selected[c("peptide_sequence","mhc_key","all_origin_ids")])
  expect_identical(serialize(x,NULL),before);expect_equal(a$input_records$native_grade,x$native_grade)
  expect_equal(a$input_records$native_rank,x$native_rank)
})
test_that("scenario budgets look up full-pool curves without re-binning or duplicating maths",{
  x<-composition_fixture();x<-x[x$provider=="test-model",];before<-serialize(x,NULL)
  profiles<-list(zero=composition_profile(rep(0,4)),certain=composition_profile(rep(1,4)),tenth=composition_profile())
  y<-publication_scenarios(x,profiles,c(0,1,2,20))
  z<-subset(y$budgets,scenario=="tenth")
  expect_equal(z$N_tested,c(0,1,2,3));expect_equal(z$expected_positive_units,c(0,.1,.2,.3))
  expect_equal(z$probability_at_least_one,1-.9^c(0,1,2,3))
  expect_identical(z$budget_status[4],"FEWER_THAN_REQUESTED")
  expect_true(all(subset(y$targets,scenario=="zero")$reachable==FALSE))
  expect_equal(subset(y$targets,scenario=="certain")$minimum_N,rep(1,4))
  expect_true(all(diff(subset(y$curves,scenario=="tenth")$probability_at_least_one)>=0))
  expect_identical(serialize(x,NULL),before)
  expect_equal(y,publication_scenarios(x,profiles,c(0,1,2,20)))
})
test_that("incompatible bin profiles and invalid probabilities do not masquerade as scenarios",{
  x<-composition_fixture();x<-x[x$provider=="test-model",]
  p<-composition_profile();q<-p;q$binning<-yield_test_profile()$binning
  expect_error(publication_scenarios(x,list(a=p,b=q),1),"identical peptide order")
  q<-p;q$rates[1]<--.1;expect_error(publication_scenarios(x,list(a=q),1),"probabilities")
  expect_error(publication_scenarios(x,list(a=p),NA_real_),"budgets")
})
