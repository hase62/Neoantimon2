test_that("paper equations reproduce independent direct evaluation and discrepancies", {
  f <- Neoantimon2:::priority_paper_value
  I <- 120;W <- 800;M <- 2;A <- .2;E <- 3;CCF <- .95
  LI <- function(x) 1/(1+exp(.015*(x-500)))
  LR <- function(x) 1/(1+exp(5*(x-2)))
  C <- 1/(1-exp(-30*(CCF-.8)))
  expected <- LI(I)*tanh(tanh(5*A)*E)*(1-2^-M*LI(W))*C
  expect_equal(f(I,W,M,A,E,CCF,"IC50")$with_ccf,expected,tolerance=1e-12)
  rank_expected <- LR(.5)*tanh(tanh(5*A)*E)*(1-2^-M*LR(4))*C
  expect_equal(f(.5,4,M,A,E,CCF,"BA_rank")$with_ccf,rank_expected,tolerance=1e-12)
  old <- LI(I)*tanh(A*tanh(5*A)*E)*(1-2^-M*LI(W))*C
  expect_false(isTRUE(all.equal(old,expected)))
  expect_equal(f(I,W,M,NA_real_,E,CCF,"IC50")$with_ccf,0)
  expect_equal(f(I,W,M,0,E,CCF,"IC50")$with_ccf,0)
  expect_true(is.na(f(I,W,M,A,E,.8,"IC50")$with_ccf))
  expect_identical(f(I,W,M,A,E,.8,"IC50")$status,"singular_or_nonfinite")
  expect_lt(f(I,W,M,A,E,.5,"IC50")$with_ccf,0)
  expect_identical(f(I,W,M,A,E,.5,"IC50")$status,"finite_negative_do_not_prioritize")
  expect_identical(f(I,W,M,A,E,.8+1e-12,"IC50")$status,"finite_unstable_do_not_prioritize")
  expect_true(is.na(f(I,W,M,A,E,NA_real_,"IC50")$with_ccf))
  expect_equal(f(I,W,M,A,E,NA_real_,"IC50")$without_ccf,expected/C)
  expect_true(is.na(f(I,W,M,A,NA_real_,CCF,"IC50")$with_ccf))
  grid <- vapply(c(0,1,100,500,1000,1e10),function(x)f(x,W,M,A,E,CCF,"IC50")$with_ccf,0)
  expect_true(all(diff(grid)<=0));expect_true(all(is.finite(grid)))
  expect_error(f(I,W,M,2,E,CCF,"IC50"),"VAF")
  expect_error(f(I,W,-1,A,E,CCF,"IC50"),"Invalid")
})

test_that("one-run input/priority demonstrations retain identity and native outputs", {
  source(system.file("examples","source-evidence.R",package="Neoantimon2"),local=TRUE)
  demo <- run_source_evidence_demo(tempfile("source-evidence-test-"))
  expect_equal(unique(demo$views$A$dna_vaf),.2)
  expect_true(all(demo$views$A$prediction_status=="NOT_RUN"))
  if(!is.null(demo$B)) {
    b <- Neoantimon2:::input_evidence_read(demo$B)$event_evidence[[1]]$rna
    expect_equal(b$ref_reads,1);expect_equal(b$alt_reads,3);expect_equal(b$usable_depth,4)
    expect_equal(b$vaf,.75)
  }
  expect_true(all(demo$views$D_EL$category=="EL"));expect_true(all(demo$views$D_IM$category=="IM"))
  expect_true(all(demo$views$E$vaf_used==0));expect_true(all(is.na(demo$views$E$dna_vaf)))
  expect_true(all(demo$views$E$evidence_adjusted==0))
  expect_true(all(is.finite(demo$views$E$prediction_only)))
  results <- neo_results(demo$A)
  expect_identical(results,neo_results(neo_project(path=demo$A$path)))
  q <- neo_explain(demo$A,candidate_id=results$candidate_id[[1]])
  expect_equal(q$candidates[[1]]$input_evidence$event_evidence[[1]]$tumor$vaf,.2)
  expect_setequal(unique(demo$views$D_EL$origin_id),unique(unlist(neo_results(demo$A)$origin_ids)))
  for (id in unique(demo$views$D_EL$candidate_id)) {
    z <- demo$views$D_EL[demo$views$D_EL$candidate_id==id,]
    expect_length(unique(z$prediction_rank),1L)
  }
  expect_error(neo_priorities(demo$A,method="neoantimon-paper-2020-v1"),"expression_unit")
  expect_identical(neo_results(demo$A),results)
  if (!is.null(demo$B)) {
    source <- Neoantimon2:::somatic_source_files(demo$B)[[1]]
    event <- Neoantimon2:::read_runtime_snapshot(demo$B)$events[[1]]
    source$kind <- "dna_bam"
    expect_length(Neoantimon2:::input_bam_observations(source,list(event),demo$B),0L)
    source$options$mode <- "independent"
    counts <- Neoantimon2:::input_bam_observations(source,list(event),demo$B)
    expect_equal(counts[[event$event_id]]$tumor$vaf,.75)
    source$kind <- "rna_bam";source$options$sample_name <- "wrong"
    expect_error(Neoantimon2:::input_bam_observations(source,list(event),demo$B),"read-group")
    source$options$sample_name <- NULL
    event$ref <- "T"
    expect_error(Neoantimon2:::input_bam_observations(source,list(event),demo$B),"REF mismatch")
    event$ref <- "A";event$alt <- "C";event$pos <- 32L;event$end <- 32L
    no <- Neoantimon2:::input_bam_observations(source,list(event),demo$B)[[event$event_id]]$rna
    expect_equal(no$alt_reads,0);expect_null(no$vaf);expect_identical(no$status,"insufficient_depth")
    event$event_type <- "SV"
    no <- Neoantimon2:::input_bam_observations(source,list(event),demo$B)[[event$event_id]]$rna
    expect_null(no$alt_reads);expect_identical(no$status,"unsupported")
    # Two independently specified full-contig CIGAR fixtures, one REF and one ALT.
    for (type in c("INS","DEL")) {
      seq <- "ATGGAACCAAAAAAAAAAAATAAACCCGGGTAA"
      alt <- if(type=="INS") paste0(substr(seq,1,4),"T",substring(seq,5)) else
        paste0(substr(seq,1,4),substring(seq,6))
      cigar <- if(type=="INS") "4M1I29M" else "4M1D28M"
      root <- tempfile();dir.create(root);sam <- file.path(root,"indel.sam")
      row <- function(id,sequence,cigar)paste(id,0,"chr1",1,60,cigar,"*",0,0,sequence,
        paste(rep("I",nchar(sequence)),collapse=""),"RG:Z:rg1",sep="\t")
      writeLines(c("@HD\tVN:1.6\tSO:coordinate","@SQ\tSN:chr1\tLN:33","@RG\tID:rg1\tSM:sample-1",
        row("ref",seq,"33M"),row("alt",alt,cigar)),sam)
      source$resolved_path <- Rsamtools::asBam(sam,destination=file.path(root,"indel"),indexDestination=TRUE)
      source$index_path <- paste0(source$resolved_path,".bai")
      event$event_type <- type;event$pos <- 4L;event$end <- if(type=="INS") 4L else 5L
      event$ref <- if(type=="INS") "G" else "GA";event$alt <- if(type=="INS") "GT" else "G"
      counts <- Neoantimon2:::input_bam_observations(source,list(event),demo$B)[[event$event_id]]$rna
      expect_equal(counts$ref_reads,1);expect_equal(counts$alt_reads,1);expect_equal(counts$vaf,.5)
      writeLines("corrupt index",file.path(root,"bad.bai"));source$index_path <- file.path(root,"bad.bai")
      expect_error(suppressWarnings(Neoantimon2:::input_bam_observations(source,list(event),demo$B)))
    }
  }
})
