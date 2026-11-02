# Independent oracle: exhaustive substring scan of explicit proteins and changed
# positions/junctions. Does not call peptide_hits/map/enumerate to derive truth.
range_oracle <- function(sequence,lengths,positions=integer(),junctions=integer()) {
  out <- character()
  for (len in lengths) if (len<=nchar(sequence)) for (start in seq_len(nchar(sequence)-len+1L)) {
    end <- start+len-1L
    if (any(positions>=start & positions<=end) || any(junctions>start & junctions<=end))
      out <- c(out,paste(len,start,substring(sequence,start,end),sep=":"))
  }
  sort(out)
}
range_observed <- function(context,lengths) {
  x <- Neoantimon2:::peptide_enumerate(context,peptide_test_policy(lengths,lengths))$peptides
  sort(vapply(x,function(p)paste(p$peptide_length,p$peptide_start_aa,p$peptide_sequence,sep=":"),""))
}
test_that("each 8--12mer and their simultaneous union exhaust missense windows", {
  for (position in c(1L,12L,29L)) {
    context <- peptide_test_long(position)
    for (lengths in c(as.list(8:12),list(8:12))) {
      expected <- range_oracle(context$translation$mutant_protein,lengths,position)
      expect_identical(range_observed(context,lengths),expected)
    }
  }
  expect_identical(Neoantimon2:::peptide_policy(list(class_i_lengths=12:8))$class_i_lengths,as.list(8:12))
  for (bad in list(numeric(),NA_real_,8.5,0,-8,Inf,"8:12"))
    expect_error(Neoantimon2:::peptide_policy(list(class_i_lengths=bad)),"Invalid peptide lengths")
})
test_that("insertions, deletion junctions, novel tails and termination enumerate independently", {
  normal <- "MACDEFGHIKLMNPQRSTVWYACDEFGHIKL"
  cases <- list(
    list(first=12L,last=11L,insert="AAA",positions=12:14,junctions=integer(),frame="unchanged"),
    list(first=12L,last=14L,insert="",positions=integer(),junctions=12L,frame="unchanged"),
    list(first=12L,last=nchar(normal),insert="AWYAWYAWYAWYAWY",positions=12:26,junctions=integer(),frame="shifted_unrestored"),
    list(first=12L,last=nchar(normal),insert="",positions=integer(),junctions=12L,frame="unchanged"))
  for (x in cases) {
    context <- peptide_test_long(); t <- context$translation
    t$mutant_protein <- paste0(substring(normal,1,x$first-1L),x$insert,substring(normal,x$last+1L))
    interval <- list(normal_aa_start=x$first,normal_aa_end=x$last,
      mutant_aa_start=x$first,mutant_aa_end=x$first+nchar(x$insert)-1L,
      normal_sequence=if(x$last<x$first) "" else substring(normal,x$first,x$last),mutant_sequence=x$insert)
    t$normal_tumor_difference_intervals <- list(interval)
    t$reference_tumor_difference_intervals <- list(interval);t$frame_effect <- x$frame
    context$translation <- t
    for (lengths in c(as.list(8:12),list(8:12)))
      expect_identical(range_observed(context,lengths),range_oracle(t$mutant_protein,lengths,x$positions,x$junctions))
  }
  context <- peptide_test_long(); context$translation$mutant_protein <- "MAC"
  context$translation$patient_normal_protein <- "MAT";context$translation$reference_protein <- "MAT"
  interval <- list(normal_aa_start=3L,normal_aa_end=3L,mutant_aa_start=3L,mutant_aa_end=3L,normal_sequence="T",mutant_sequence="C")
  context$translation$normal_tumor_difference_intervals <- list(interval)
  context$translation$reference_tumor_difference_intervals <- list(interval)
  expect_identical(range_observed(context,8:12),character())
})

test_that("8--12mer deduplication keeps distinct transcript/Event Origins and unique HLA pairs", {
  base <- reconstruction_test_project();cfg <- base$config
  cfg$peptides <- peptide_test_policy(8:12,15L)
  project <- neo_project(samples=base$samples,config=cfg,
    hla=list("patient-1"=c("HLA-A*02:01","HLA-A*02:01","HLA-B*07:02")))
  event <- reconstruction_test_row("e1","chr1",4L,"G","A");event$event_key <- "key:e1"
  context <- peptide_test_contexts(list(event),"TA.1")[[1]]
  # These are protein-window unit fixtures, not assertions about the short
  # genomic reconstruction fixture. Translation truth is tested separately.
  context$translation <- peptide_test_long()$translation
  context$gene_symbol <- "GENE_A"
  other <- context;other$reconstruction_id <- "reconstruction:other"
  other$transcript_id <- "OTHER.1";other$gene_id <- "GOTHER.1";other$gene_symbol <- "GENE_B"
  other$somatic_event_ids <- "e2"
  event2 <- event;event2$event_id <- "e2";event2$event_key <- "key:e2"
  x <- Neoantimon2:::peptide_construct(project,list(context,other),list(event,event2))
  expected <- range_oracle(context$translation$mutant_protein,8:12,12L)
  expected_sequences <- sub("^[0-9]+:[0-9]+:","",expected)
  expect_setequal(vapply(x$candidates,`[[`,"","peptide_sequence"),expected_sequences)
  expect_length(Filter(function(o)o$peptide_length %in% 8:12,x$origins),2L*length(expected))
  expect_length(Filter(function(o)o$peptide_length==15L,x$origins),
    2L*length(range_oracle(context$translation$mutant_protein,15L,12L)))
  expect_length(x$candidates,2L*length(unique(expected_sequences)))
  for (c in x$candidates) {
    expect_length(c$origin_ids,2L)
    expect_setequal(unlist(c$gene_symbols),c("GENE_A","GENE_B"))
    expect_true(c$mhc_key %in% c("HLA-A*02:01","HLA-B*07:02"))
    own <- Filter(function(o)o$origin_id %in% unlist(c$origin_ids),x$origins)
    expect_setequal(unlist(lapply(own,`[[`,"source_event_ids")),c("e1","e2"))
    expect_true(all(vapply(own,function(o)o$peptide_sequence==c$peptide_sequence,TRUE)))
  }
  expect_false(anyDuplicated(vapply(x$candidates,`[[`,"","candidate_id"))>0)
  # Candidate creation is independent of installed provider support.
  expect_length(x$predictor_inputs,length(x$candidates))
  expect_true(all(vapply(x$predictor_inputs,function(q)q$supported_residue_check=="ADAPTER_PENDING",TRUE)))
})
