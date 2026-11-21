# Run from the Core checkout with an explicit fresh output directory and fixed-model artifacts.
# Rscript --vanilla scripts/development/validate_candidate_review.R OUTPUT FIXED_MODEL PVAC_COMPUTATION
args <- commandArgs(TRUE)
stopifnot(length(args)==3L)
pkgload::load_all(quiet=TRUE)
target <- normalizePath(args[[1]],mustWork=FALSE)
if(file.exists(target))stop("Choose a fresh acceptance output directory")
dir.create(target,recursive=TRUE,mode="0700")
model <- normalizePath(args[[2]],mustWork=TRUE)
pvac <- normalizePath(args[[3]],mustWork=TRUE)
input <- file.path(target,"peptides.fasta")
writeLines(c(">first original header", "acdef", "ghik", ">second duplicate with separate provenance",
  "ACDEFGHIK", ">third distinct", "LMNPQRSTV"),input)
config <- neo_config(env=character(),checksum="sha256",paths=list(run_root=file.path(target,"runs")),
  peptides=list(class_i_lengths=9L,class_ii_lengths=15L))
p <- neo_run(neo_project(samples=data.frame(sample_id="software-sample",patient_id="software-subject"),
  hla=list("software-subject"="HLA-A*02:01"),config=config,project_id="fasta-review-software-verification",
  inputs=list(list(input_id="three-records",input_role="peptide_fasta",sample_id="software-sample",path=input))))
stopifnot(nrow(neo_results(p))==2L,length(external_graph(p)$origins)==3L)
before <- neo_results(p);manifest_sha <- external_sha(manifest_file(p$path))
unknown <- neo_analyze_external(p,file.path(model,"converter-original.tsv"),list(mode="evidence_only"),
  reference_percentiles=file.path(model,"manifest.json"))
stopifnot(all(vapply(external_analysis_read(unknown)$payload$reference_percentiles,
  function(r)r$status=="NOT_CALCULATED",TRUE)))
a <- neo_analyze_external(p,file.path(model,"predictions.tsv"),list(mode="evidence_only"),
  reference_percentiles=file.path(model,"manifest.json"))
calibration <- external_analysis_read(a)$payload$reference_percentiles
receipt <- jsonlite::read_json(file.path(model,"receipt.json"))
computed <- Filter(function(r)r$status=="CALCULATED",calibration)
stopifnot(length(computed)==2L)
for(i in seq_along(computed))stopifnot(abs(computed[[i]]$reference_percentile-unlist(receipt$expected_percentiles)[[i]])<1e-12)
write_html <- function(name,analyses=list(a),candidate_id=NULL,privacy="private")
  suppressWarnings(html_report_write(p,file.path(target,name),analyses=analyses,candidate_id=candidate_id,
    privacy=privacy,confirm_private=privacy=="private",label="Fixed-model software verification; not biological calibration"))
write_html("private-v1.html")
# A clearly labeled additional status analysis changes displayed evidence, never original predictions.
extra <- external_read_predictions(file.path(model,"predictions.tsv"))
extra <- extra[extra$metric_kind=="raw_score",]
extra$provider<-"review-regression-fixture";extra$predictor_version<-"fixture-v1";extra$model_version<-"fixture-v1"
extra$task<-"external:review-regression";extra$metric<-"status_fixture";extra$unit<-"unitless"
extra$score<-"";extra$status<-"not_evaluated";extra$reason<-"Authored review-resume software status fixture, no biological prediction"
extra$run_id<-"review-regression-fixture";extra$source_row<-paste0("authored-status:",seq_len(nrow(extra)))
path<-file.path(target,"additional-status.tsv");utils::write.table(extra,path,sep="\t",quote=TRUE,row.names=FALSE,na="")
b<-neo_analyze_external(p,path,list(mode="evidence_only"))
write_html("private-v2.html",list(a,b));write_html("share.html",list(a,b),privacy="share")
ids<-neo_results(p)$candidate_id
write_html("private-subset.html",list(a),ids[[1]])
foreign<-neo_run(neo_project(samples=p$samples,hla=lapply(p$hla,unlist),config=config,
  project_id="different-project",inputs=list(list(input_id="three-records",input_role="peptide_fasta",sample_id="software-sample",path=input))))
suppressWarnings(html_report_write(foreign,file.path(target,"private-foreign.html"),privacy="private",confirm_private=TRUE))
# Existing molecular example uses authored VCF/reference fixtures, not FASTA or a new scientific analysis.
example<-new.env();sys.source(system.file("examples","first-real-run.R",package="Neoantimon2"),example)
molecular<-neo_run(do.call(neo_project,example$first_real_run_example(file.path(target,"molecular-example"))))
molecular_before<-neo_results(molecular);molecular_manifest<-external_sha(manifest_file(molecular$path))
suppressWarnings(html_report_write(molecular,file.path(target,"private-molecular.html"),privacy="private",confirm_private=TRUE))
stopifnot(identical(neo_results(neo_project(path=molecular$path)),molecular_before),
  identical(external_sha(manifest_file(molecular$path)),molecular_manifest))
stopifnot(identical(neo_results(neo_project(path=p$path)),before),identical(external_sha(manifest_file(p$path)),manifest_sha))
stopifnot(identical(neo_external_results(a),neo_external_results(a$analysis_id,neo_project(path=p$path))))
math <- list()
for(item in list(list("MHCflurry",500,"lower_is_better",5.068),list("BigMHC_EL",0.5,"higher_is_better",0.617))) {
  values<-scan(file.path(pvac,paste0(item[[1]],"-9mer.tsv")),skip=1L,quiet=TRUE)
  observed<-reference_percentile_value(sort(values,method="radix"),item[[2]],item[[3]])
  stopifnot(round(observed,3)==item[[4]])
  math[[length(math)+1L]]<-list(method=item[[1]],score=item[[2]],observed_unrounded=observed,count=length(values))
}
bundle<-neo_report_bundle(p,dir=file.path(target,"report-tables"),analyses=list(a,b))
evidence<-list(status="PASS",r_version=R.version.string,package_version=as.character(packageVersion("Neoantimon2")),
  project_path=p$path,project_manifest_sha256=manifest_sha,candidates=2L,origins=3L,
  candidate_ids=as.list(ids),original_prediction_rows=nrow(external_read_predictions(file.path(model,"predictions.tsv"))),
  reference_rows=calibration,analysis_ids=list(real_fixed_model=a$analysis_id,additional_status=b$analysis_id,unknown_model=unknown$analysis_id),
  expected_reference_percentiles=receipt$expected_percentiles,pvac_computation=math,
  original_native_results_preserved=TRUE,source_manifest_preserved=TRUE,analyses_reopened=TRUE,
  molecular_regression=list(candidates=nrow(molecular_before),manifest_sha256=molecular_manifest,
    native_results_unchanged=TRUE,origins_and_links_preserved=TRUE),
  browser_status="Run tests/browser/candidate_review.py separately")
writeLines(report_json(evidence),file.path(target,"acceptance.json"))
cat("FASTA_REFERENCE_REVIEW_R_ACCEPTANCE_PASS\n")
