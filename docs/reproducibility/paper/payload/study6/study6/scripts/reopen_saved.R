args<-commandArgs(TRUE);B<-normalizePath(args[[1L]]);library(Neoantimon2);library(jsonlite)
write_full<-function(x,path)write_json(x,path,auto_unbox=TRUE,null='null',na='null',digits=NA,pretty=TRUE,dataframe='rows')
receipts<-list.files(file.path(B,'observed'),pattern='receipt.json',recursive=TRUE,full.names=TRUE)
summary<-list();fail<-FALSE
for(file in receipts) {
 x<-read_json(file);if(!identical(x$status,'COMPLETED'))next
 out<-dirname(file);started<-format(Sys.time(),'%Y-%m-%dT%H:%M:%OS6%z');rr<-list(case_id=x$case_id,started_at=started,R_process_pid=Sys.getpid(),Core_package_version=as.character(packageVersion('Neoantimon2')),Core_package_path=system.file(package='Neoantimon2'),library_paths=as.list(.libPaths()),session=capture.output(sessionInfo()),public_calls=list('neo_project(path)','neo_events','neo_origins','neo_results','neo_export_bundle','neo_external_results(id,project)','neo_external_explain(id,project)','neo_priorities(analysis=id)'))
 tryCatch({
  before<-readRDS(file.path(out,'first_process_views.rds'));p<-neo_project(path=x$project_path)
  # NO neo_run, predictor, reconstruction or translation is invoked in this separate R process.
  graph<-neo_export_bundle(p)$graph;native<-neo_results(p);origins<-neo_origins(p);events<-neo_events(p)
  external<-if(!is.null(x$analysis_id))neo_external_results(x$analysis_id,project=p) else NULL
  priorities<-if(!is.null(x$analysis_id))neo_priorities(p,analysis=x$analysis_id) else NULL
  rr$comparisons<-list(events=identical(before$events,events),origins=identical(before$origins,origins),native=identical(before$native,native),graph=identical(before$graph,graph),external=identical(before$external,external),priorities=identical(before$priorities,priorities))
  stopifnot(all(unlist(rr$comparisons)))
  if(nrow(native)) {
   explanation<-neo_external_explain(x$analysis_id,candidate_id=head(native$candidate_id,3L),project=p)
   write_full(explanation,file.path(out,'reopened_external_explain.json'))
  }
  rr$status<-'PASS';rr$scientific_recomputation<-FALSE
  write_full(list(events=events,origins=origins,native=native,external=external,priorities=priorities),file.path(out,'reopened_public_views.json'))
 },error=function(e){rr$status<<-'FAILED';rr$error<<-conditionMessage(e);fail<<-TRUE})
 rr$finished_at<-format(Sys.time(),'%Y-%m-%dT%H:%M:%OS6%z');write_full(rr,file.path(out,'separate_process_reopen.json'));summary[[length(summary)+1L]]<-rr
}
write_full(summary,file.path(B,'observed/SEPARATE_PROCESS_REOPEN_RESULTS.json'));quit(status=if(fail)1 else 0)
