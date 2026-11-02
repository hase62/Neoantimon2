# Additional independent sequence truth, not full-sample biological validation.
args <- commandArgs(TRUE);root <- normalizePath(args[[1]]);out <- args[[2]]
.libPaths(c(args[[3]],.libPaths()));library(Neoantimon2);library(testthat)
setwd(file.path(root,"tests/testthat"))
helper <- readLines("helper-reconstruction.R")
helper <- gsub('testthat::test_path("..", "fixtures", "reconstruction")',
 deparse(file.path(root,"tests/fixtures/reconstruction")),helper,fixed=TRUE)
eval(parse(text=helper));source("helper-peptide.R")
# Literal proteins independently established from the fixture's spliced DNA
# (see test-translation-t009.R), not inferred from enumerator output.
cases <- list(
 minus_frameshift=list(rows=list(reconstruction_test_row("minus","chrR",46L,"GT","G")),tx="TR.1",protein="MNPGNPGNPNPGK",positions=2:13),
 phased_restore=list(rows=list(reconstruction_test_row("i","chrF",6L,"A","AC",side="LEFT"),
 reconstruction_test_row("d","chrF",26L,"CC","C",side="LEFT")),tx="TF.1",protein="MKPRETGKP",positions=4:6),
 phased_stop=list(rows=list(reconstruction_test_row("a","chrF",10L,"G","T",side="LEFT"),
 reconstruction_test_row("b","chrF",12L,"G","A",side="LEFT")),tx="TF.1",protein="MKP",positions=integer()))
counts <- list();expected <- list()
for(name in names(cases)) {
 x <- cases[[name]];ct <- peptide_test_contexts(x$rows,x$tx)[[1]]
 stopifnot(identical(ct$translation$mutant_protein,x$protein))
 for(len in 8:12) {
  truth <- character()
  if(nchar(x$protein)>=len)for(start in seq_len(nchar(x$protein)-len+1L))
   if(any(x$positions>=start & x$positions<=start+len-1L))truth <- c(truth,substring(x$protein,start,start+len-1L))
  actual <- Neoantimon2:::peptide_enumerate(ct,peptide_test_policy(len,len))$peptides
  stopifnot(identical(vapply(actual,`[[`,"","peptide_sequence"),truth))
  counts[[length(counts)+1L]] <- data.frame(case=name,length=len,expected=length(truth),actual=length(actual),equal=TRUE)
  if(length(truth))expected[[length(expected)+1L]] <- data.frame(case=name,length=len,sequence=truth)
 }
 combined <- Neoantimon2:::peptide_enumerate(ct,peptide_test_policy(8:12,8:12))$peptides
 stopifnot(length(combined)==sum(vapply(Filter(function(r)r$case==name,counts),function(r)r$expected,0L)))
}
write.csv(do.call(rbind,counts),file.path(out,"tables/strand-phase-lengths.csv"),row.names=FALSE)
write.csv(do.call(rbind,expected),file.path(out,"fixtures/strand-phase-expected-peptides.csv"),row.names=FALSE)
cat("15 per-length comparisons and 3 simultaneous-union comparisons passed\n")
