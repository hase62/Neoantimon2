# Small internal regression controls complement the separately executed public-API chains.
bounded_contexts <- function(rows) {
  p <- reconstruction_test_project(); b <- Neoantimon2:::reference_bundle(p)
  old <- Neoantimon2:::reconstruction_contexts(p,b,rows)
  Neoantimon2:::translation_contexts(p,b,rows,old)
}
bounded_one <- function(rows,tx) Filter(function(x)x$transcript_id==tx,bounded_contexts(rows))[[1L]]

test_that("empty Sec positions have true empty annotation IDs; declared and used stay distinct", {
  x <- bounded_one(list(reconstruction_test_row("m","chrF",4L,"A","G")),"TF.1")
  expect_identical(x$translation$spec$translation_exception_ids,character())
  expect_identical(Neoantimon2:::peptide_translation_provenance(x)$annotated_translation_exception_ids,list())
  sec <- bounded_one(list(reconstruction_test_row("s","chrS2",4L,"A","G")),"TS.1")
  expect_identical(sec$translation$spec$translation_exception_ids,"GENCODE_SEC:TS.1:8")
  expect_identical(sec$translation$sec_sites_used,"GENCODE_SEC:TS.1:8")
  changed <- bounded_one(list(reconstruction_test_row("s","chrS2",22L,"T","C")),"TS.1")
  expect_identical(changed$translation$spec$translation_exception_ids,"GENCODE_SEC:TS.1:8")
  expect_identical(changed$translation$sec_sites_used,character())
})

test_that("phase background is separate from all somatic phase conditions", {
  a<-reconstruction_test_row("a","chrF",4L,"A","G")
  b<-reconstruction_test_row("b","chrF",5L,"A","C")
  x<-bounded_one(list(a,b),"TF.1")
  expect_identical(x$phase_state,"RESOLVED_OR_UNAMBIGUOUS")
  pc<-Neoantimon2:::peptide_phase_condition(x)
  expect_identical(pc$status,"CONDITIONAL")
  expect_true(any(vapply(pc$phase_context,function(y)y$relation=="UNKNOWN",TRUE)))
  expect_identical(Neoantimon2:::peptide_translation_provenance(x)$phase_condition,pc)
  a<-reconstruction_test_row("a","chrF",4L,"A","G",side="LEFT")
  for(side in c("LEFT","RIGHT")) {
    b<-reconstruction_test_row("b","chrF",5L,"A","C",side=side)
    expect_true(all(vapply(bounded_contexts(list(a,b)),function(y)
      Neoantimon2:::peptide_phase_condition(y)$status=="RESOLVED",TRUE)))
  }
})

test_that("conditional-only paths cannot rescue native eligibility; independent support can", {
  event<-somatic_test_event(); peptide<-somatic_test_peptide()
  raw<-Neoantimon2:::somatic_raw_construct(list(event),peptide,
    list(somatic_test_source("dna_reads",list(somatic_test_row()))))
  policy<-Neoantimon2:::somatic_evidence_policy(list())$eligibility
  link<-raw$candidate_links[[1L]]; path<-link$origin_paths[[1L]]
  resolved<-Neoantimon2:::somatic_origin_eligibility(path,link$candidate_id,
    setNames(raw$event_evidence,"event:a"),policy)
  expect_identical(resolved$status,"eligible")
  path$phase_condition<-list(status="CONDITIONAL")
  conditional<-Neoantimon2:::somatic_origin_eligibility(path,link$candidate_id,
    setNames(raw$event_evidence,"event:a"),policy)
  expect_identical(conditional$status,"not_evaluated")
  expect_true(any(vapply(conditional$reasons,function(r)
    r$reason_code=="molecular_phase_not_resolved",TRUE)))
  link$origin_paths<-list(link$origin_paths[[1L]],path);raw$candidate_links<-list(link)
  expect_identical(Neoantimon2:::somatic_eligibility_construct(raw,policy)$candidate_eligibility[[1L]]$eligibility_status,"eligible")
  raw$candidate_links[[1L]]$origin_paths<-list(path)
  expect_identical(Neoantimon2:::somatic_eligibility_construct(raw,policy)$candidate_eligibility[[1L]]$eligibility_status,"not_evaluated")
  raw$event_evidence[[1L]]$tumor$alt_reads<-0L
  expect_identical(Neoantimon2:::somatic_eligibility_construct(raw,policy)$candidate_eligibility[[1L]]$eligibility_status,"ineligible")
})

test_that("vaccine conditional check respects actual per-Origin eligibility", {
  records<-vaccine_records();f<-vaccine_fixture();binding<-vaccine_binding()
  r<-records[[1L]];r$trajectory<-NULL
  r$evidence$origins[[1L]]$origin_warning_codes<-list("PHASE_UNCERTAIN")
  project<-function(rec)Neoantimon2:::vaccine_project_records(list(rec),f$contexts,f$relations,binding$normalized)$candidates[[1L]]$eligibility
  expect_identical(project(r)$status,"ineligible")
  expect_true(any(vapply(project(r)$checks,function(x)x$constraint=="molecular_phase_evaluable"&&x$status=="insufficient_evidence",TRUE)))
  o<-r$evidence$origins[[1L]];o$origin_id<-"resolved:other";o$origin_warning_codes<-list()
  r$evidence$origins[[2L]]<-o
  expect_identical(project(r)$status,"eligible")
  r$evidence$somatic_evidence$eligibility<-list(origin_results=list(list(origin_id=o$origin_id,status="ineligible"),list(origin_id="conditional",status="not_evaluated")))
  expect_identical(project(r)$status,"ineligible")
  r$evidence$somatic_evidence$eligibility$origin_results[[1L]]$status<-"eligible"
  expect_identical(project(r)$status,"eligible")
})

test_that("explicit cis UTR substitution is included before identity and trans/unknown remain guarded", {
  s<-reconstruction_test_row("s","chrF",38L,"T","C",side="LEFT")
  c<-reconstruction_test_row("c","chrF",41L,"C","T",side="LEFT")
  cis<-bounded_one(list(s,c),"TF.1")
  expect_identical(cis$translation$mutant_protein,"MKPGKPGKPQSG")
  expect_setequal(cis$somatic_event_ids,c("s","c"))
  expect_identical(cis$translation$mutant_stop_origin,37:39)
  expect_identical(cis$translation$mutant_protein,bounded_one(list(c,s),"TF.1")$translation$mutant_protein)
  u<-reconstruction_test_row("c","chrF",41L,"C","T")
  expect_identical(bounded_one(list(s,u),"TF.1")$translation$reason,"EXTENSION_PHASE_UNRESOLVED")
  t<-reconstruction_test_row("c","chrF",41L,"C","T",side="RIGHT")
  trans<-bounded_one(list(s,t),"TF.1")
  expect_identical(trans$translation$mutant_protein,"MKPGKPGKPQPG")
  expect_identical(trans$somatic_event_ids,"s")
  # Genome position50 -> transcript40, downstream of the recovered stop at37:39.
  far<-reconstruction_test_row("far","chrF",50L,"C","T",side="LEFT")
  alone<-bounded_one(list(s),"TF.1")
  after<-bounded_one(list(s,far),"TF.1")
  expect_identical(after$translation$mutant_protein,alone$translation$mutant_protein)
  expect_identical(after$somatic_event_ids,"s")
})
