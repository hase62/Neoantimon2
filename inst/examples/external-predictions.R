# SYNTHETIC / TEST-SCALE / NON-CLINICAL. No predictor, network or checkout.
library(Neoantimon2)
source(system.file('examples','quickstart.R',package='Neoantimon2'),local=TRUE)
native_results <- neo_results(p)
bundle <- neo_export_bundle(p)
queries <- read.delim(file.path(bundle$path,'pmhc_queries.tsv'),colClasses='character')
stopifnot(nrow(queries)>=3L)
selected_candidate <- bundle$graph$candidates[[1L]]
i <- which(queries$peptide_sequence==selected_candidate$peptide_sequence & queries$mhc_key==selected_candidate$mhc_key)[[1L]]
j <- setdiff(seq_len(nrow(queries)),i)[[1L]]
native_explanation <- neo_explain(p,candidate_id=selected_candidate$candidate_id)
original_files <- list.files(p$path,recursive=TRUE,full.names=TRUE)
original_hashes <- vapply(original_files,function(f)digest::digest(file=f,algo='sha256',serialize=FALSE),'')
row <- function(query,score='0',status='evaluated',scope='peptide_hla',mhc=queries$mhc_key[[query]],task='external:offline_example',source='X:1') {
  data.frame(schema_version='1.0.0',peptide_sequence=queries$peptide_sequence[[query]],
    scope=scope,mhc_key=mhc,task=task,provider='synthetic-model',predictor_version='1',model_version='1',
    run_id='offline-example',metric='score',metric_kind='raw_score',metric_semantics='raw',score=score,
    score_direction='higher_is_better',unit='unitless',status=status,
    reason=if(status=='evaluated') '' else 'synthetic unsupported query',source_row=source)
}
predictions <- rbind(row(i),row(j,score='',status='unsupported',source='X:2'),
  row(i,score='0.9',mhc='HLA-B*07:02',source='X:3'),
  row(i,score='0.5',scope='peptide',mhc='',task='external:intrinsic_example',source='X:4'))
prediction_file <- file.path(input_root,'external-standard.tsv')
write.table(predictions,prediction_file,sep='\t',quote=TRUE,qmethod='double',row.names=FALSE,na='')
external <- neo_analyze_external(p,prediction_file,policy=list(mode='external_rank',
  provider='synthetic-model',run_id='offline-example',task='external:offline_example',metric='score',score_direction='higher_is_better'))
external_results <- neo_external_results(external)
external_explanation <- neo_external_explain(external,selected_candidate$candidate_id)
selected <- external_results[external_results$candidate_id==selected_candidate$candidate_id,,drop=FALSE]
missing <- external_results$external_status=='NOT_PREDICTED'
unsupported <- external_results$external_status=='unsupported'
stopifnot(nrow(selected)==1L,selected$external_score==0,selected$external_rank==1L,
  any(missing),any(unsupported),all(is.na(external_results$external_rank[missing|unsupported])),
  all(is.na(external_results$external_score[missing|unsupported])),
  length(external_explanation$candidates[[1L]]$origins)==length(selected_candidate$origin_ids),
  identical(neo_results(p),native_results),
  identical(neo_explain(p,candidate_id=selected_candidate$candidate_id),native_explanation),
  identical(vapply(original_files,function(f)digest::digest(file=f,algo='sha256',serialize=FALSE),''),original_hashes),
  identical(neo_external_results(external$analysis_id,project=neo_project(path=p$path)),external_results))
unmatched <- read.delim(file.path(external$path,'unmatched.tsv'),colClasses='character')
stopifnot(nrow(unmatched)==1L,unmatched$mhc_key=='HLA-B*07:02')
external_example <- list(status='PASS',candidates=nrow(external_results),
  origin_count=length(selected_candidate$origin_ids),zero_score=selected$external_score,
  not_predicted=sum(missing),unsupported=sum(unsupported),unmatched=nrow(unmatched),
  native_unchanged=TRUE,project_path=p$path,external_analysis_id=external$analysis_id)
print(external_example)
