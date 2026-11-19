args <- commandArgs(TRUE); B <- normalizePath(args[[1L]]); case_id <- args[[2L]]
library(Neoantimon2); library(jsonlite)
stopifnot(length(getNamespaceExports('Neoantimon2'))==20L)
now <- function() format(Sys.time(),'%Y-%m-%dT%H:%M:%OS6%z')
write_json_full <- function(x,path) { if(inherits(x,'neo_explanation'))x<-unclass(x); write_json(x,path,auto_unbox=TRUE,null='null',na='null',digits=NA,pretty=TRUE,dataframe='rows') }
case <- Filter(function(x)x$case_id==case_id,read_json(file.path(B,'inputs/cases.json')))[[1L]]
resume <- if(length(args)>=3L)args[[3L]] else NULL
out <- file.path(B,'observed',paste0(case_id,if(!is.null(resume))'_resume02' else '')); stopifnot(!dir.exists(out)); dir.create(out)
receipt <- list(runtime=list(package_path=system.file(package='Neoantimon2'),package_version=as.character(packageVersion('Neoantimon2')),library_paths=as.list(.libPaths()),R_process_pid=Sys.getpid(),R_version=R.version.string,session=capture.output(sessionInfo()),loaded_namespace_paths=lapply(loadedNamespaces(),function(n)list(name=n,version=as.character(packageVersion(n)),path=system.file(package=n)))),case_id=case_id,started_at=now(),api_calls=list(),status='RUNNING',input_case=case)
write_json_full(receipt,file.path(out,'receipt.json'))
status <- 0L
tryCatch({
 ref <- file.path(B,case$reference)
 roles <- c(genome_fasta='genome.fa',genome_fai='genome.fa.fai',annotation_gtf='annotation.gtf',transcript_fasta='transcripts.fa',protein_coding_transcript_fasta='coding.fa',protein_fasta='proteins.fa',selenocysteine_metadata='sec.tsv',mane_summary='mane.tsv',mane_gtf='mane.gtf')
 resources <- lapply(names(roles),function(role){path<-file.path(ref,roles[[role]]);list(path=path,required=TRUE,type='file',resource_id=if(role=='annotation_gtf') 'gencode_primary_assembly_all_catalogs_mane_supplemental_v1' else role,version=if(startsWith(role,'genome')) 'GRCh38.p14' else if(startsWith(role,'mane')) '1.5' else '50',source='NEW_SYNTHETIC_CONTROL_NOT_REAL_REFERENCE',checksum=list(algorithm='sha256',scope='file',value=digest::digest(file=path,algo='sha256',serialize=FALSE)))})
 names(resources)<-paste0('reference.',names(roles))
 inputs<-lapply(case$vcfs,function(v){path<-file.path(B,v$path);list(input_id=paste(case_id,v$role,sep='-'),input_role=paste0(v$role,'_vcf'),sample_id='NV_SAMPLE',path=path,checksum=list(algorithm='sha256',scope='file',value=digest::digest(file=path,algo='sha256',serialize=FALSE)))})
 lens<-unlist(case$lengths)
 targets <- if(!is.null(case$targets))unlist(case$targets) else 'MKACDFGH'
 target_keys<-as.list(unlist(lapply(targets,function(seq)vapply(c('HLA-A*02:01','HLA-B*07:02'),function(h)Neoantimon2:::biological_id('candidate_key',list(species='NCBITaxon:9606',hla_class='I',ligand_sequence=seq,mhc_key=h),Neoantimon2:::peptide_contract()),'',USE.NAMES=FALSE))))
 extra <- if(isTRUE(case$native)) list(somatic_evidence=list(),vaccine_set=list(version='vaccine-set-policy-1.0.0',max_peptides=as.integer(case$max_selected),fill_to_capacity=TRUE,allowlist=target_keys)) else list()
 cfg<-do.call(neo_config,c(list(env=character(),paths=list(data_root=file.path(B,'inputs'),input_root=file.path(B,'inputs'),run_root=file.path(B,'runs',case_id)),scientific=list(genome_assembly='GRCh38.p14'),resources=resources,checksum='sha256',peptides=list(class_i_lengths=as.integer(lens[lens<15]),class_ii_lengths=15L)),extra))
 receipt$api_calls<-c(receipt$api_calls,list('neo_config','neo_project',if(is.null(resume))'neo_run' else 'neo_project(path=completed_saved_run)'))
 p<-neo_project(samples=data.frame(sample_id='NV_SAMPLE',patient_id='NV_SUBJECT'),config=cfg,hla=list(NV_SUBJECT=c('HLA-A*02:01','HLA-B*07:02','HLA-DRA*01:01~HLA-DRB1*01:01')),inputs=inputs)
 receipt$project_path<-p$path;write_json_full(receipt,file.path(out,'receipt.json'))
 if(isTRUE(case$native)) {
   # Only input identity preparation uses this internal importer, also used by package regression fixtures.
   # Molecular reconstruction, eligibility and vaccine selection below run through public neo_run/results/explain.
   imported <- Neoantimon2:::event_import_project(p)
   dna <- file.path(out,'synthetic_dna_counts.tsv')
   dna_rows <- lapply(imported,function(e){
    label <- Filter(function(v) identical(v$contig,e$contig) && as.integer(v$pos)==e$pos && identical(v$ref,e$ref) && identical(v$alt,e$alt),case$events)
    if(length(label)!=1L)stop('Controlled genomic identity not unique')
    alt<-as.integer(case$dna_counts[[label[[1L]]$event_id]])
    data.frame(event_id=e$event_id,tumor_ref_reads=20L-alt,tumor_alt_reads=alt,tumor_depth=20L,normal_ref_reads=20L,normal_alt_reads=0L,normal_depth=20L)
   })
   write.table(do.call(rbind,dna_rows),dna,sep='\t',quote=FALSE,row.names=FALSE)
   cfg$somatic_evidence <- list(sources=list(list(kind='dna_reads',sample_id='NV_SAMPLE',path=dna,source='DECLARED_SYNTHETIC_COUNTS_NOT_OBSERVED_PATIENT',source_version='1')))
   p<-neo_project(samples=data.frame(sample_id='NV_SAMPLE',patient_id='NV_SUBJECT'),config=cfg,hla=list(NV_SUBJECT=c('HLA-A*02:01','HLA-B*07:02','HLA-DRA*01:01~HLA-DRB1*01:01')),inputs=inputs)
   receipt$input_preparation_internal_helper<-'event_import_project: exact genomic identity to generated event_id only; no status injection'
 }
 p<-if(is.null(resume))neo_run(p) else neo_project(path=resume)
 receipt$project_path<-p$path;receipt$resumed_saved_pipeline_without_neo_run<-!is.null(resume)
 receipt$pipeline_completed_at<-now()
 e<-neo_events(p);o<-neo_origins(p);r<-neo_results(p)
 selection_paths<-list.files(p$path,pattern='selection.rds$',recursive=TRUE,full.names=TRUE)
 write_json_full(lapply(selection_paths,function(f)list(path=f,sha256=digest::digest(file=f,algo='sha256',serialize=FALSE),value=readRDS(f))),file.path(out,'actual_native_selection.json'))
 receipt$api_calls<-c(receipt$api_calls,list('neo_events','neo_origins','neo_results'))
 write_json_full(e,file.path(out,'events.json'));write_json_full(o,file.path(out,'origins.json'));write_json_full(r,file.path(out,'native_results.json'));write_json_full(attr(r,'vaccine_set'),file.path(out,'vaccine_summary.json'))
 saveRDS(if(nrow(r))neo_explain(p,candidate_id=r$candidate_id) else NULL,file.path(out,'all_native_explain.rds'))
 # Reading completed native stage artifacts is observation, never a separate reconstruction call.
 stages<-c(sequence_reconstruction='sequence_reconstruction.rds',translation_resolution='translation_resolution.rds',peptide_candidates='peptide_candidates.rds',somatic_priority='somatic_priority.rds',somatic_eligibility='somatic_eligibility.rds',somatic_evidence='somatic_evidence_raw.rds')
 for(n in names(stages)) {path<-file.path(p$path,'intermediate',stages[[n]]);if(file.exists(path))write_json_full(readRDS(path),file.path(out,paste0(n,'.json')))}
 file.copy(file.path(p$path,'intermediate/phase_resolution.json'),file.path(out,'phase_resolution.json'))
 if(nrow(r))write_json_full(neo_explain(p,candidate_id=head(r$candidate_id,3L)),file.path(out,'native_explain.json'))
 files<-list.files(p$path,recursive=TRUE,full.names=TRUE);files<-files[!dir.exists(files)];before<-vapply(files,function(f)digest::digest(file=f,algo='sha256',serialize=FALSE),'')
 bundle<-neo_export_bundle(p);receipt$api_calls<-c(receipt$api_calls,list('neo_export_bundle'))
 receipt$export_path<-bundle$path;write_json_full(bundle$graph,file.path(out,'public_graph.json'));write_json_full(bundle$manifest,file.path(out,'bundle_manifest.json'))
 receipt$counts<-list(Events=nrow(e),Origins=nrow(o),unique_sequences=length(unique(o$peptide_sequence)),candidates=nrow(r),links=length(bundle$graph$links))
 analysis<-NULL
 if(nrow(r)) {
  q<-read.delim(file.path(bundle$path,'pmhc_queries.tsv'),colClasses='character',check.names=FALSE)
  q<-q[order(q$peptide_sequence,q$mhc_key,method='radix'),,drop=FALSE]
  control<-lapply(seq_len(nrow(q)),function(i){
   if(i%%13L==0L)return(NULL)
   st<-if(i%%17L==0L)'failed' else if(i%%11L==0L)'unsupported' else 'evaluated'
   data.frame(schema_version='1.0.0',peptide_sequence=q$peptide_sequence[[i]],scope='peptide_hla',mhc_key=q$mhc_key[[i]],task='external:controlled_nv',provider='NV-CONTROL-NOT-A-PREDICTOR',predictor_version='NV1',model_version='FIXTURE',run_id=case_id,metric='controlled_score',metric_kind='raw_score',metric_semantics='raw',score=if(st=='evaluated')as.character(((i-1L)%%4L)/4) else '',score_direction='higher_is_better',unit='unitless',status=st,reason=if(st=='evaluated')'' else paste('CONTROLLED',st,'NOT_REAL_MODEL'),source_row=paste0('CONTROL:',i))})
  predictions<-do.call(rbind,Filter(Negate(is.null),control));predfile<-file.path(out,'controlled_predictions.tsv');write.table(predictions,predfile,sep='\t',quote=TRUE,qmethod='double',row.names=FALSE,na='')
  write_json_full(list(identity='CONTROLLED_IMPORT_NOT_REAL_MODEL_OUTPUT',rule='lexicographic sequence/HLA; index%%4 score including ties+zero; every11unsupported/17failed/13omitted; no predictor inference',query_count=nrow(q)),file.path(out,'controlled_prediction_policy.json'))
  analysis<-neo_analyze_external(p,predfile,policy=list(mode='external_rank',provider='NV-CONTROL-NOT-A-PREDICTOR',run_id=case_id,task='external:controlled_nv',metric='controlled_score',score_direction='higher_is_better',predictor_version='NV1',model_version='FIXTURE'))
  er<-neo_external_results(analysis);pr<-neo_priorities(p,analysis=analysis)
  write_json_full(er,file.path(out,'external_results.json'));write_json_full(pr,file.path(out,'priorities.json'));write_json_full(neo_external_explain(analysis,head(r$candidate_id,3L)),file.path(out,'external_explain.json'))
  receipt$analysis_id<-analysis$analysis_id;receipt$analysis_path<-analysis$path
  receipt$api_calls<-c(receipt$api_calls,list('neo_analyze_external(external_rank)','neo_external_results','neo_external_explain','neo_priorities(model-v1)'))
  stopifnot(identical(neo_results(p),r),identical(neo_origins(p),o),identical(neo_events(p),e))
 }
 after<-vapply(files,function(f)digest::digest(file=f,algo='sha256',serialize=FALSE),'');stopifnot(identical(before,after))
 write_json_full(list(paths=as.list(files),before_sha256=as.list(before),after_sha256=as.list(after),unchanged=identical(before,after)),file.path(out,'source_file_preservation.json'))
 saveRDS(list(events=e,origins=o,native=r,external=if(!is.null(analysis))neo_external_results(analysis) else NULL,priorities=if(!is.null(analysis))neo_priorities(p,analysis=analysis) else NULL,graph=bundle$graph),file.path(out,'first_process_views.rds'),version=3)
 receipt$runtime_finished<-list(package_path=system.file(package='Neoantimon2'),session=capture.output(sessionInfo()),loaded_namespace_paths=lapply(loadedNamespaces(),function(n)list(name=n,version=as.character(packageVersion(n)),path=system.file(package=n))))
 receipt$status<-'COMPLETED';receipt$native_unchanged_by_external_import<-TRUE;receipt$original_project_files_unchanged<-TRUE
},error=function(e){receipt$status<<-'FAILED';receipt$error<<-conditionMessage(e);receipt$error_class<<-class(e);receipt$error_code<<-attr(e,'neo_error_code');status<<-1L;message('NV CASE FAILURE: ',conditionMessage(e))})
receipt$finished_at<-now();receipt$exit_status<-status;write_json_full(receipt,file.path(out,'receipt.json'));print(receipt[c('case_id','status','counts','error')]);quit(status=status)
