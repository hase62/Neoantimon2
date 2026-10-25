# Offline documentation acceptance against an existing installed v0.3.0 package.
args <- commandArgs(TRUE); root <- normalizePath(args[[1L]]); out <- args[[2L]]
stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE);out<-normalizePath(out)
library(Neoantimon2);stopifnot(as.character(packageVersion('Neoantimon2'))=='0.3.0')
helper <- new.env();sys.source(file.path(root,'docs/ai/run-analysis.R'),helper)
checks <- list()
pass <- function(name) {checks[[name]] <<- 'PASS';cat('PASS:',name,'\n')}
reject <- function(expr,pattern) {
  error<-tryCatch({force(expr);NULL},error=conditionMessage)
  stopifnot(!is.null(error),grepl(pattern,error));invisible(error)
}
# Validate entrypoint traversal, not the behavior of a separately launched agent.
for (entry in c('AGENTS.md','CLAUDE.md')) {
 text<-readLines(file.path(root,entry));stopifnot(any(grepl('docs/ai/AI_INSTRUCTIONS.md',text,fixed=TRUE)),
  any(grepl('docs/ai/ANALYSIS_PLAYBOOK.md',text,fixed=TRUE)))
 pass(paste(entry,'entrypoint reading/link traversal'))
}
e<-new.env();source(system.file('examples','quickstart.R',package='Neoantimon2'),local=e)
p<-e$p;r<-neo_results(p)
stopifnot(nrow(r)>0L,nrow(neo_events(p))>0L,nrow(neo_origins(p))>0L,
 all(r$prediction_status=='NOT_RUN'),all(r$grading_status=='NOT_GRADED'))
x<-neo_explain(p,candidate_id=r$candidate_id[[1L]])
stopifnot(identical(neo_results(neo_project(path=p$path)),r),
 identical(neo_explain(neo_project(path=p$path),candidate_id=r$candidate_id[[1L]]),x));pass('offline Quick Start, explanation and identical reopening')
pass('NOT_RUN and NOT_GRADED remain unavailable states, not biological negatives')
# Re-execute all four existing sample/patient/HLA examples from the installed release.
cohort<-new.env();source(system.file('examples','multi-sample.R',package='Neoantimon2'),local=cohort)
for(i in 1:4) {
 q<-get(paste0('p',i),cohort);z<-neo_results(q);ss<-get(paste0('samples_',i),cohort)
 stopifnot(setequal(z$sample_id,ss$sample_id),all(z$patient_id==ss$patient_id[match(z$sample_id,ss$sample_id)]))
 if(i>=3)stopifnot(!any(z$mhc_key[z$patient_id=='P1']=='HLA-B*07:02'),any(z$mhc_key[z$patient_id=='P2']=='HLA-B*07:02'))
 pass(paste('cohort',i,'sample/patient/HLA associations'))
}
# Translate a real local YAML workflow declaration; no developer paths in source.
# Keep inputs/references outside the new output session.
input<-file.path(out,'inputs');dir.create(input)
vcf<-file.path(input,'S1.vcf');file.copy(e$vcf,vcf);lines<-readLines(vcf);writeLines(sub('sample-1$','S1',lines),vcf)
config<-unclass(e$cfg);config$paths$input_root<-input;config$paths$data_root<-input;config$paths$run_root<-file.path(out,'unused')
configfile<-file.path(input,'config.yaml');yaml::write_yaml(config,configfile)
m<-list(workflow_declaration='neoantimon2-agent-v1',mode='somatic',reference_assembly='GRCh38.p14',
 resource_config=configfile,run_directory=file.path(out,'execution'),
 samples=list(list(sample_id='S1',patient_id='P1',vcf=vcf)),hla=list(P1=c('HLA-A*02:01','HLA-B*07:02')),
 prediction_providers=list(),optional_evidence=list(),external_predictions=list())
manifest<-file.path(input,'manifest.yaml');yaml::write_yaml(m,manifest)
pre<-helper$agent_execute(manifest);stopifnot(!dir.exists(m$run_directory));pass('manifest translation and read-only preflight')

# Exercise the checked-in three-sample/two-patient template after path substitution.
template<-yaml::read_yaml(file.path(root,'docs/ai/INPUT_MANIFEST.example.yaml'),eval.expr=FALSE)
template$resource_config<-configfile;template$run_directory<-file.path(out,'template-not-executed')
for(i in seq_along(template$samples)) {
 id<-template$samples[[i]]$sample_id;f<-file.path(input,paste0(id,'.vcf'))
 template$samples[[i]]$vcf<-f
}
# Each synthetic file has its authoritative sample header; no source header is rewritten.
for(i in seq_along(template$samples)) {
 id<-template$samples[[i]]$sample_id
 if(id!='S1')writeLines(sub('S1$',id,readLines(vcf)),template$samples[[i]]$vcf)
}
translated<-helper$agent_declare(template)
stopifnot(length(translated$project$samples)==3L,setequal(names(translated$project$hla),c('P1','P2')),
 !dir.exists(template$run_directory))
pass('checked-in two-patient three-sample manifest translates with explicit Class II pair')
bad<-m;bad$optional_evidence<-list(list(kind='rna_expression',path='unconfigured.tsv'))
reject(helper$agent_declare(bad),'optional_evidence')
bad<-m;bad$external_predictions<-list(list(path='unconfigured.tsv'))
reject(helper$agent_declare(bad),'external_predictions')
pass('unconfigured optional/external inputs rejected instead of ignored')
# Agent policy is intentionally stricter than core's allowed missing-HLA path.
for (value in list(NULL,list(),list(P1='A02'),list(P1='HLA-DRB1*15:01'),list(S1='HLA-A*02:01'))) {
 bad<-m;bad$hla<-value;reject(helper$agent_declare(bad),'HLA');
}
pass('agent rejects missing/ambiguous HLA and wrong patient associations')
# Multi-sample headers and identity mismatch fail before execution, not rewritten.
old<-readBin(vcf,'raw',n=file.info(vcf)$size)
writeLines(sub('\tS1$','\tS1\tS2',readLines(vcf)),vcf)
reject(helper$agent_declare(m),'one VCF sample');writeBin(old,vcf)
writeLines(sub('\tS1$','\tOTHER',readLines(vcf)),vcf)
reject(helper$agent_declare(m),'identity');writeBin(old,vcf)
pass('multi-sample VCF and inconsistent header rejected without inference')
missing_cfg<-neo_config(env=character(),predictors=list(list(provider='netmhcpan',task='mhc1_binding',
 executable=file.path(input,'absent-netmhcpan'))))
probe<-neo_check_predictors(config=missing_cfg,providers='netmhcpan',smoke_test=FALSE)
stopifnot(!isTRUE(probe$details[[1L]]$capability$available),probe$details[[1L]]$smoke_test$status=='NOT_RUN')
pass('missing predictor safely reported; no model run')
original<-vapply(c(manifest,configfile,vcf),function(f)digest::digest(file=f,algo='sha256',serialize=FALSE),'')
p2<-helper$agent_execute(manifest,execute=TRUE);summary<-jsonlite::read_json(file.path(m$run_directory,'report.json'))
stopifnot(summary$inputs_unchanged,summary$reopened_identical,
 file.exists(file.path(m$run_directory,'sessionInfo.txt')),file.exists(file.path(m$run_directory,'invocation.txt')),
 identical(original,vapply(c(manifest,configfile,vcf),function(f)digest::digest(file=f,algo='sha256',serialize=FALSE),'')),
 identical(neo_results(p2),neo_results(neo_project(path=p2$path))))
reject(helper$agent_execute(manifest,execute=TRUE),'already exists')
pass('authorized helper execution, saved report, original hashes and completed run preservation')
# Exact public export inventory and documented signature references.
expected<-c('neo_project','neo_config','neo_run','neo_results','neo_events','neo_origins','neo_multisample','neo_explain',
 'neo_export_bundle','neo_analyze_external','neo_external_results','neo_external_explain',
 'neo_download_resources','neo_check_resources','neo_check_predictors','neo_setup_predictor')
stopifnot(setequal(getNamespaceExports('Neoantimon2'),expected))
signatures<-lapply(expected,function(n) names(formals(getExportedValue('Neoantimon2',n))));names(signatures)<-expected
pass('16 released public API signatures inspected')
jsonlite::write_json(list(status='PASS',checks=checks,pass=length(checks),fail=0L,
 agent_sessions_launched=FALSE,agent_specific_limit='Entrypoint traversal and R execution tested; no independent Codex/Claude Code session',
 package_version='0.3.0',package_commit=packageDescription('Neoantimon2')$RemoteSha,
 signatures=signatures,manifest=manifest,helper_report=file.path(m$run_directory,'report.json')),
 file.path(out,'validation.json'),pretty=TRUE,auto_unbox=TRUE)
