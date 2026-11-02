# SYNTHETIC A--E demonstration plus strictly read-only reopening checks.
# No real predictor, ASCAT execution, biological or clinical validation implied.
library(Neoantimon2)
run_reopen_evidence_demo <- function(output) {
  source(system.file("examples","source-evidence.R",package="Neoantimon2"),local=TRUE)
  demo <- run_source_evidence_demo(output)
  hashes <- function(root) {
    f <- sort(list.files(root,recursive=TRUE,full.names=TRUE,all.files=TRUE))
    f <- f[!dir.exists(f)]
    stats::setNames(vapply(f,function(x)digest::digest(file=x,algo="sha256",serialize=FALSE),""),f)
  }
  verification <- list()
  for (name in c("A","B","C","E")) {
    p <- demo[[name]];if(is.null(p)) next
    before <- neo_results(p);id <- before$candidate_id[[1L]]
    explain <- neo_explain(p,candidate_id=id);numeric <- neo_priorities(p)
    checksums <- hashes(p$path)
    q <- neo_project(path=p$path) # No neo_run() is required for these reads.
    stopifnot(identical(neo_results(q),before),
      identical(neo_explain(q,candidate_id=id),explain),
      identical(neo_priorities(q),numeric),identical(hashes(p$path),checksums))
    verification[[name]] <- list(rows=nrow(before),files_verified=length(checksums),
      read_only=TRUE,results_explanation_and_priorities_identical=TRUE)
  }
  jsonlite::write_json(verification,file.path(output,"reopen-verification.json"),
    auto_unbox=TRUE,pretty=TRUE)
  invisible(demo)
}
if (sys.nframe()==0L) {
  args <- commandArgs(TRUE)
  if(length(args)!=1L) stop("Usage: Rscript reopen-evidence.R /new/output")
  run_reopen_evidence_demo(args[[1L]])
}
