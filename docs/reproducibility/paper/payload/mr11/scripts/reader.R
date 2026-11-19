a<-commandArgs(TRUE);C<-normalizePath(a[[1L]]);W<-normalizePath(a[[2L]]);library(Neoantimon2);library(jsonlite)
full<-function(x,p)write_json(x,p,auto_unbox=TRUE,null='null',na='null',digits=NA,pretty=TRUE,dataframe='rows')
rows<-list();failure<-FALSE
for(id in c('M01','M02','M03','M04')){D<-file.path(W,'observed',id);first<-read_json(file.path(D,'receipt.json'));rr<-list(case_id=id,R_process_pid=Sys.getpid(),first_R_process_pid=first$runtime$R_process_pid,session=capture.output(sessionInfo()),package_path=system.file(package='Neoantimon2'),scientific_reconstruction=FALSE)
 tryCatch({
  stopifnot(first$status=='COMPLETED',rr$R_process_pid!=rr$first_R_process_pid);v<-readRDS(file.path(D,'first_views.rds'));p<-neo_project(path=first$project_path)
  rr$comparisons<-list(events=identical(v$events,neo_events(p)),origins=identical(v$origins,neo_origins(p)),native=identical(v$native,neo_results(p)),graph=identical(v$graph,neo_export_bundle(p)$graph),external=identical(v$external,neo_external_results(first$analysis_id,project=p)),priorities=identical(v$priorities,neo_priorities(p,analysis=first$analysis_id)))
  # Compare the complete saved translation/peptide artifacts bytewise as well as public graph views.
  artifacts<-c('translation_resolution','peptide_candidates');rr$artifact_checks<-lapply(artifacts,function(n){f<-file.path(p$path,'intermediate',paste0(n,'.rds'));list(name=n,sha256=digest::digest(file=f,algo='sha256',serialize=FALSE))})
  stopifnot(all(unlist(rr$comparisons)),identical(rr$artifact_checks,first$stage_artifacts_sha));full(neo_explain(p,candidate_id=head(v$native$candidate_id,3L)) |> unclass(),file.path(D,'reopened_explain.json'));rr$status<-'PASS'
 },error=function(e){rr$status<<-'FAILED';rr$error<<-conditionMessage(e);failure<<-TRUE})
 full(rr,file.path(D,'separate_reader.json'));rows[[length(rows)+1L]]<-rr
}
full(rows,file.path(W,'SEPARATE_READER.json'));quit(status=if(failure)1 else 0)
