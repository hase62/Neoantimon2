priority_fixture <- function(grades=c("A","B","C","D","E","INVALID","MISSING_DNA",
                                     "MISSING_PREDICTOR","II","SEC")) {
  spec <- list(A=list(rank=.5,rna=1L),B=list(rank=.5,rna=0L),
    C=list(rank=2,rna=1L),D=list(rank=2,rna=NULL),E=list(rank=2.01,rna=1L),
    INVALID=list(rank=.1,rna=1L,elig="ineligible"),
    MISSING_DNA=list(rank=.1,rna=1L,elig="not_evaluated"),
    MISSING_PREDICTOR=list(rank=NULL,rna=1L),
    II=list(rank=1,rna=1L,cl="II"),SEC=list(rank=NULL,rna=1L,sec=TRUE))
  spec <- spec[grades]
  candidates <- lapply(names(spec),function(id) {
    x <- spec[[id]]
    list(candidate_id=id,sample_id="sample-1",hla_class=if (is.null(x$cl)) "I" else x$cl,
      mhc_key=if (identical(x$cl,"II")) "HLA-DRA*01:01~HLA-DRB1*04:01" else "HLA-A*02:01",
      ligand_sequence=if (isTRUE(x$sec)) "ACDEFGHUK" else "ACDEFGHIK",
      peptide_sequence=if (isTRUE(x$sec)) "ACDEFGHUK" else "ACDEFGHIK",
      origin_ids=list(paste0("origin:",id)),self_exact_match=FALSE,self_status="EVALUATED",
      final_grade="NOT_GRADED",grading_status="NOT_GRADED",prediction_status="NOT_RUN",
      not_evaluated_codes=list("NOT_RUN"))
  })
  event <- lapply(names(spec),function(id) {
    x <- spec[[id]]
    list(event_id=paste0("event:",id),evidence_id=paste0("evidence:",id),
      tumor=list(status="observed_support",alt_reads=5L,vaf=.5),
      rna=list(status=if (is.null(x$rna)) "unavailable" else if (x$rna>0)
        "observed_support" else "observed_no_alt",alt_reads=x$rna,
        vaf=if (is.null(x$rna)) NULL else x$rna/10),
      ccf=list(status="ambiguous",bounded_value=NULL,solutions=list(list(
        multiplicity=1,compatible=TRUE),list(multiplicity=2,compatible=TRUE))),
      purity=list(value=.8),cn=list(status="mapped",total_cn=2,normal_cn=2),
      loh=list(status="no_evidence"))
  })
  links <- lapply(names(spec),function(id) list(candidate_id=id,sample_id="sample-1",
    origin_paths=list(list(origin_id=paste0("origin:",id),
      somatic_event_ids=list(paste0("event:",id)),
      event_evidence_ids=list(paste0("evidence:",id)),expression_evidence_ids=list()))))
  eligibility <- lapply(names(spec),function(id) {
    x <- spec[[id]]; status <- if (is.null(x$elig)) "eligible" else x$elig
    list(candidate_id=id,eligibility_status=status,reason_records=list(),
      origin_results=list(list(origin_id=paste0("origin:",id),status=status)))
  })
  evidence <- lapply(names(spec),function(id) {
    x <- spec[[id]]
    if (is.null(x$rank)) return(NULL)
    cl <- if (is.null(x$cl)) "I" else x$cl
    list(candidate_id=id,predictor=if (cl=="I") "netmhcpan" else "netmhciipan",
      task=if (cl=="I") "mhc1_presentation" else "mhc2_presentation",
      metric="EL_rank",status="evaluated",percentile_rank=x$rank,raw_score=x$rank,
      run_id=if (cl=="I") "run:I" else "run:II",raw_row_index=1L,
      predictor_version="fixture",model_version="fixture")
  })
  evidence <- Filter(Negate(is.null),evidence)
  runs <- list(list(provider="netmhcpan",run_id="run:I",model_version="fixture",
    evidence=Filter(function(x)x$predictor=="netmhcpan",evidence),statuses=list()),
    list(provider="netmhciipan",run_id="run:II",model_version="fixture",
      evidence=Filter(function(x)x$predictor=="netmhciipan",evidence),statuses=list()))
  list(candidates=candidates,somatic=list(raw=list(event_evidence=event,
    candidate_links=links,expression_evidence=list()),
    eligibility=list(candidate_eligibility=eligibility)),runs=runs)
}

test_that("inclusive class thresholds and all grade states preserve missingness", {
  p <- Neoantimon2:::somatic_priority_policy()
  band <- Neoantimon2:::somatic_pmhc_band
  expect_identical(band(.5,"I",p)$band,"STRONG")
  expect_identical(band(.500001,"I",p)$band,"BINDER")
  expect_identical(band(2,"I",p)$band,"BINDER")
  expect_identical(band(2.000001,"I",p)$band,"OUTSIDE")
  expect_identical(band(1,"II",p)$band,"STRONG")
  expect_identical(band(1.000001,"II",p)$band,"BINDER")
  expect_identical(band(5,"II",p)$band,"BINDER")
  expect_identical(band(5.000001,"II",p)$band,"OUTSIDE")
  f <- priority_fixture()
  rows <- Neoantimon2:::somatic_priority_construct(f$candidates,f$somatic,f$runs,p)
  by_id <- stats::setNames(rows,vapply(rows,`[[`,"","candidate_id"))
  expect_identical(vapply(by_id[c("A","B","C","D","E","INVALID","MISSING_DNA",
    "MISSING_PREDICTOR","II","SEC")],`[[`,"","final_grade"),
    c(A="A",B="B",C="C",D="D",E="E",INVALID="INVALID",
      MISSING_DNA="NOT_GRADED",MISSING_PREDICTOR="NOT_GRADED",II="A",SEC="NOT_GRADED"))
  expect_identical(by_id$B$rna_status,"observed_no_alt")
  expect_identical(by_id$D$rna_status,"unavailable")
  expect_null(by_id$D$rna_vaf)
  expect_identical(by_id$A$ccf_status,"ambiguous")
  expect_null(by_id$A$ccf_value)
  expect_length(by_id$A$ccf_solutions,2L)
  expect_null(by_id$INVALID$rank_within_sample_class)
  expect_identical(by_id$SEC$sec_status,"contains_U")
  expect_identical(by_id$II$rank_within_sample_class,1L)
})

test_that("selection rejects duplicate primary rows and leaves secondary predictions parallel", {
  p <- Neoantimon2:::somatic_priority_policy()
  f <- priority_fixture("A")
  second <- f$runs[[1L]]$evidence[[1L]]
  second$raw_row_index <- 2L
  f$runs[[1L]]$evidence <- c(f$runs[[1L]]$evidence,list(second))
  x <- Neoantimon2:::somatic_priority_construct(f$candidates,f$somatic,f$runs,p)[[1L]]
  expect_identical(x$final_grade,"NOT_GRADED")
  expect_identical(x$primary_prediction_status,"ambiguous")
  expect_length(x$primary_applicable_rows,2L)
  f$runs[[1L]]$evidence <- f$runs[[1L]]$evidence[1L]
  alt <- second; alt$predictor <- "mhcflurry"; alt$metric <- "presentation_percentile"
  alt$percentile_rank <- 50; alt$run_id <- "run:secondary"
  f$runs[[3L]] <- list(provider="mhcflurry",run_id="run:secondary",
    model_version="fixture",evidence=list(alt),statuses=list())
  x <- Neoantimon2:::somatic_priority_construct(f$candidates,f$somatic,f$runs,p)[[1L]]
  expect_identical(x$final_grade,"A")
  expect_true(x$predictor_disagreement)
  expect_identical(x$evaluable_provider_count,2L)
  expect_length(x$secondary_predictions,1L)
})

test_that("lexicographic rank is stable under row permutation and separates classes", {
  p <- Neoantimon2:::somatic_priority_policy()
  f <- priority_fixture(c("A","B","C","D","E","INVALID","II"))
  first <- Neoantimon2:::somatic_priority_construct(f$candidates,f$somatic,f$runs,p)
  second <- Neoantimon2:::somatic_priority_construct(rev(f$candidates),f$somatic,f$runs,p)
  mapping <- function(rows) stats::setNames(vapply(rows,function(x)
    if (is.null(x$rank_within_sample_class)) NA_integer_ else x$rank_within_sample_class,0L),
    vapply(rows,`[[`,"","candidate_id"))
  expect_identical(mapping(first),mapping(second))
  expect_identical(unname(mapping(first)[c("A","B","C","D","E")]),1:5)
  expect_identical(unname(mapping(first)["II"]),1L)
  expect_true(is.na(mapping(first)["INVALID"]))
})

test_that("each named ranking component breaks only its own tie", {
  p <- Neoantimon2:::somatic_priority_policy()
  base <- function(id) list(candidate_id=id,sample_id="s",hla_class="I",
    final_grade="B",primary_rank=.5,rna_status="unavailable",
    transcript_tpm=NULL,gene_tpm=NULL,ccf_value=NULL,tumor_alt_reads=NULL,
    exact_self_match=NULL,rank_within_sample_class=NULL)
  pair <- function(change_a,change_b) {
    a <- base("a"); b <- base("b")
    a[names(change_a)] <- change_a; b[names(change_b)] <- change_b
    rows <- Neoantimon2:::somatic_priority_rank(list(b,a),p)
    rows[[1L]]$candidate_id
  }
  expect_identical(pair(list(primary_rank=.4),list(primary_rank=.5)),"a")
  expect_identical(pair(list(rna_status="observed_support"),
    list(rna_status="observed_no_alt")),"a")
  expect_identical(pair(list(rna_status="partial_support"),
    list(rna_status="unavailable")),"a")
  expect_identical(pair(list(rna_status="unavailable"),
    list(rna_status="observed_no_alt")),"a")
  expect_identical(pair(list(transcript_tpm=0),list(transcript_tpm=NULL)),"a")
  expect_identical(pair(list(transcript_tpm=2),list(transcript_tpm=1)),"a")
  expect_identical(pair(list(gene_tpm=2),list(gene_tpm=1)),"a")
  expect_identical(pair(list(ccf_value=.8),list(ccf_value=.4)),"a")
  expect_identical(pair(list(tumor_alt_reads=10L),list(tumor_alt_reads=5L)),"a")
  expect_identical(pair(list(exact_self_match=FALSE),list(exact_self_match=TRUE)),"a")
  expect_identical(pair(list(),list()),"a")
  expect_null(base("x")$transcript_tpm)
})

priority_compound_fixture <- function(rank=.5,second_rna="observed_no_alt",
                                      second_tumor=5L) {
  f <- priority_fixture("A")
  f$runs[[1L]]$evidence[[1L]]$percentile_rank <- rank
  second <- f$somatic$raw$event_evidence[[1L]]
  second$event_id <- "event:second"
  second$evidence_id <- "evidence:second"
  second$rna$status <- second_rna
  second$rna$alt_reads <- if (second_rna=="unavailable") NULL else
    if (second_rna=="observed_no_alt") 0L else 1L
  second$tumor$alt_reads <- second_tumor
  f$somatic$raw$event_evidence[[1L]]$tumor$alt_reads <- 20L
  f$somatic$raw$event_evidence[[2L]] <- second
  path <- f$somatic$raw$candidate_links[[1L]]$origin_paths[[1L]]
  path$somatic_event_ids <- list("event:A","event:second")
  path$event_evidence_ids <- list("evidence:A","evidence:second")
  f$somatic$raw$candidate_links[[1L]]$origin_paths[[1L]] <- path
  f
}

test_that("compound support requires every event and preserves partial evidence", {
  p <- Neoantimon2:::somatic_priority_policy()
  assess <- function(f) Neoantimon2:::somatic_priority_construct(
    f$candidates,f$somatic,f$runs,p)[[1L]]
  partial <- assess(priority_compound_fixture())
  expect_identical(partial$final_grade,"B")
  expect_identical(partial$rna_status,"partial_support")
  expect_false(partial$rna_support)
  expect_identical(partial$tumor_alt_reads,5L)
  expect_identical(unname(vapply(partial$origin_priority_summaries[[1L]]$event_summaries,
    `[[`,0L,"tumor_alt_reads")),c(20L,5L))
  expect_identical(unname(vapply(partial$origin_priority_summaries[[1L]]$event_summaries,
    `[[`,"","rna_status")),c("observed_support","observed_no_alt"))
  expect_identical(partial$origin_priority_summaries[[1L]]$ccf_status,"multiple_events")
  expect_null(partial$ccf_value)
  binder <- assess(priority_compound_fixture(rank=2))
  expect_identical(binder$final_grade,"D")
  full <- assess(priority_compound_fixture(second_rna="observed_support"))
  expect_true(full$rna_support)
  expect_identical(full$rna_status,"observed_support")
  expect_identical(full$final_grade,"A")
  expect_identical(assess(priority_compound_fixture(rank=2,
    second_rna="observed_support"))$final_grade,"C")
  unavailable <- assess(priority_compound_fixture(second_rna="unavailable"))
  expect_identical(unavailable$rna_status,"partial_support")
  expect_false(unavailable$rna_support)
  missing_tumor <- assess(priority_compound_fixture(second_tumor=NULL))
  expect_null(missing_tumor$tumor_alt_reads)
  expect_null(missing_tumor$origin_priority_summaries[[1L]]$event_summaries[[2L]]$tumor_alt_reads)
})

priority_alternative_fixture <- function() {
  f <- priority_fixture(c("A","B"))
  first <- f$somatic$raw$event_evidence[[1L]]
  first$rna$status <- "observed_no_alt"
  first$rna$alt_reads <- 0L
  f$somatic$raw$event_evidence[[1L]] <- first
  second <- first
  second$event_id <- "event:z"
  second$evidence_id <- "evidence:z"
  second$rna$status <- "observed_support"
  second$rna$alt_reads <- 1L
  f$somatic$raw$event_evidence[[3L]] <- second
  a <- f$somatic$raw$candidate_links[[1L]]$origin_paths[[1L]]
  a$origin_id <- "origin:a"
  a$expression_evidence_ids <- list("expression:a")
  z <- a
  z$origin_id <- "origin:z"
  z$somatic_event_ids <- list("event:z")
  z$event_evidence_ids <- list("evidence:z")
  z$expression_evidence_ids <- list("expression:z")
  f$somatic$raw$candidate_links[[1L]]$origin_paths <- list(a,z)
  f$somatic$raw$expression_evidence <- list(
    list(evidence_id="expression:a",transcript_tpm=1,gene_tpm=2),
    list(evidence_id="expression:z",transcript_tpm=10,gene_tpm=20))
  f$somatic$eligibility$candidate_eligibility[[1L]]$origin_results <- list(
    list(origin_id="origin:a",status="eligible"),
    list(origin_id="origin:z",status="eligible"))
  f$candidates[[1L]]$origin_ids <- list("origin:a","origin:z")
  f
}

test_that("eligible alternative origins use evidence before lexical ID", {
  p <- Neoantimon2:::somatic_priority_policy()
  f <- priority_alternative_fixture()
  run <- function(f) Neoantimon2:::somatic_priority_construct(
    f$candidates,f$somatic,f$runs,p)
  first <- run(f)
  selected <- Filter(function(x)x$candidate_id=="A",first)[[1L]]
  expect_identical(selected$priority_origin_id,"origin:z")
  expect_identical(selected$priority_origin_selection_policy,
    "best-supported-eligible-origin-v1")
  expect_identical(selected$priority_source_event_ids,list("event:z"))
  expect_identical(selected$source_event_ids,list("event:A","event:z"))
  expect_identical(selected$transcript_tpm,10)
  expect_identical(selected$gene_tpm,20)
  expect_identical(selected$final_grade,"A")
  expect_length(selected$origin_priority_summaries,2L)
  expect_identical(vapply(selected$origin_priority_summaries,`[[`,"","origin_id"),
    c("origin:a","origin:z"))
  f$somatic$raw$candidate_links[[1L]]$origin_paths <-
    rev(f$somatic$raw$candidate_links[[1L]]$origin_paths)
  f$somatic$eligibility$candidate_eligibility[[1L]]$origin_results <-
    rev(f$somatic$eligibility$candidate_eligibility[[1L]]$origin_results)
  f$candidates <- rev(f$candidates)
  second <- Filter(function(x)x$candidate_id=="A",run(f))[[1L]]
  expect_identical(second,selected)
  expect_identical(vapply(first,`[[`,"","candidate_id"),
    vapply(run(f),`[[`,"","candidate_id"))
  separate <- f
  separate$somatic$raw$expression_evidence[[1L]]$transcript_tpm <- 100
  separate$somatic$raw$expression_evidence[[2L]]$transcript_tpm <- NULL
  from_z <- Filter(function(x)x$candidate_id=="A",run(separate))[[1L]]
  expect_identical(from_z$priority_origin_id,"origin:z")
  expect_null(from_z$transcript_tpm)
  expect_identical(from_z$gene_tpm,20)
  tied <- f
  tied$somatic$raw$event_evidence[[3L]]$rna$status <- "observed_no_alt"
  tied$somatic$raw$event_evidence[[3L]]$rna$alt_reads <- 0L
  tied$somatic$raw$event_evidence[[3L]]$tumor$alt_reads <- 5L
  tied$somatic$raw$expression_evidence[[2L]]$transcript_tpm <- 1
  tied$somatic$raw$expression_evidence[[2L]]$gene_tpm <- 2
  expect_identical(Filter(function(x)x$candidate_id=="A",run(tied))[[1L]]$
    priority_origin_id,"origin:a")
  ineligible <- tied
  ineligible$somatic$eligibility$candidate_eligibility[[1L]]$origin_results <-
    lapply(ineligible$somatic$eligibility$candidate_eligibility[[1L]]$origin_results,
      function(x) { x$status <- "ineligible"; x })
  ineligible$somatic$eligibility$candidate_eligibility[[1L]]$eligibility_status <-
    "ineligible"
  invalid <- Filter(function(x)x$candidate_id=="A",run(ineligible))[[1L]]
  expect_identical(invalid$final_grade,"INVALID")
  expect_null(invalid$priority_origin_id)
  expect_length(invalid$origin_priority_summaries,2L)
})

test_that("policy changes alter only priority identity", {
  p <- Neoantimon2:::somatic_priority_policy()
  q <- Neoantimon2:::somatic_priority_policy(list(thresholds=list(I=list(strong=.4))))
  expect_false(identical(Neoantimon2:::infra_hash(p),Neoantimon2:::infra_hash(q)))
  expect_identical(p$primary$I$metric,"EL_rank")
  expect_identical(p$primary$II$metric,"EL_rank")
  expect_identical(p$origin_selection_policy,"best-supported-eligible-origin-v1")
  expect_error(Neoantimon2:::somatic_priority_policy(list(
    origin_selection_policy="unknown")),"origin selection")
  expect_error(Neoantimon2:::somatic_priority_policy(list(
    primary=list(I=list(metric="EL_score")))),"percentile")
})
