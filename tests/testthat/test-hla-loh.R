hla_loh_fixture <- function(name) normalizePath(testthat::test_path("..","fixtures","hla_loh",name))

hla_loh_validation_evidence <- function(name,value) {
  directory <- Sys.getenv("NEO_HLA_LOH_EVIDENCE_DIR",unset="")
  if (nzchar(directory)) jsonlite::write_json(value,file.path(directory,name),
    pretty=TRUE,auto_unbox=TRUE,null="null")
}

hla_loh_test_project <- function(path=hla_loh_fixture("lohhla-standard.HLAlossPrediction_CI.xls"),
                                 samples=c("region-1"), genotypes=NULL) {
  if (is.null(genotypes)) genotypes <- c("HLA-A*02:01","HLA-A*24:02",
    "HLA-B*07:02","HLA-B*08:01","HLA-C*03:04","HLA-C*07:02",
    "HLA-DRA*01:01~HLA-DRB1*04:01")
  cfg <- neo_config(hla_loh=list(sources=lapply(samples,function(sample)
    list(sample_id=sample,path=path,provider_version="fixture-version-1",
      inputs=list(hla_genotype="synthetic-genotype-1",hla_fasta="synthetic-reference-1",
        purity_ploidy="synthetic-purity-1")))))
  neo_project(samples=lapply(samples,function(sample)
    list(sample_id=sample,patient_id="patient-1",role="tumor")),
    hla=list("patient-1"=genotypes),config=cfg)
}

test_that("named LOHHLA columns preserve raw and canonical allele evidence", {
  p <- hla_loh_test_project()
  source <- Neoantimon2:::hla_loh_source_files(p)
  evidence <- Neoantimon2:::hla_loh_construct(p,source,p$config$hla_loh$policy)
  by_allele <- stats::setNames(evidence,vapply(evidence,`[[`,"","canonical_hla_allele"))
  a <- by_allele[["HLA-A*02:01"]]
  expect_identical(a$loh_call,"loss_supported")
  expect_identical(a$raw_provider_allele,"hla_a_02_01")
  expect_identical(a$provider_loss_allele,"hla_a_02_01")
  expect_identical(a$allele_copy_number,0.4999)
  expect_identical(a$allele_copy_number_lower,0.31)
  expect_identical(a$allele_copy_number_upper,0.68)
  expect_identical(a$unpaired_pval_unique,0.008)
  expect_identical(a$num_mismatch_sites_covered,12)
  expect_identical(a$proportion_supportive_sites,0.83)
  expect_identical(by_allele[["HLA-A*24:02"]]$allele_role,"kept")
  expect_identical(by_allele[["HLA-B*07:02"]]$loh_call,"loss_not_supported")
  expect_identical(by_allele[["HLA-C*03:04"]]$loh_call,
    "low_copy_without_significant_imbalance")
  expect_true(all(vapply(evidence,function(x)nchar(x$raw_output_checksum)==64L,TRUE)))
  variant <- Neoantimon2:::hla_loh_parse(hla_loh_fixture(
    "lohhla-variant.HLAlossPrediction_CI.txt"))
  expect_length(variant$rows,1L)
  expect_identical(variant$columns$pval,"PVal_unique")
  expect_error(Neoantimon2:::hla_loh_columns(c("LossAllele","Pval_unique")),
    "unsupported_format")
  hla_loh_validation_evidence("parser-validation.json",list(
    fixture_kind="synthetic",standard_columns=Neoantimon2:::hla_loh_parse(
      source[[1L]]$resolved_path)$columns,
    standard_allele_calls=lapply(evidence,function(x) x[c("sample_id","hla_locus",
      "canonical_hla_allele","allele_role","loh_call","allele_copy_number",
      "allele_copy_number_lower","allele_copy_number_upper","pval_unique",
      "num_mismatch_sites_covered","proportion_supportive_sites","raw_output_checksum")]),
    variant_pval_column=variant$columns$pval,unsupported_layout="explicit_error"))
})

test_that("strict thresholds and changed policy preserve provider values", {
  p <- hla_loh_test_project()
  source <- Neoantimon2:::hla_loh_source_files(p)
  row <- Neoantimon2:::hla_loh_parse(source[[1L]]$resolved_path)$rows[[1L]]
  map <- list(status="mapped",allele="HLA-A*02:01")
  other <- list(status="mapped",allele="HLA-A*24:02")
  base <- p$config$hla_loh$policy
  make <- function(cn,pval,policy=base) {
    x <- row
    x$values$cn1 <- as.character(cn); x$values$pval <- as.character(pval)
    Neoantimon2:::hla_loh_record(p,source[[1L]],x,1L,map,other,"potential_loss",policy)
  }
  expect_identical(make(0.4999,0.0099)$loh_call,"loss_supported")
  expect_identical(make(0.5,0.0099)$loh_call,"allelic_imbalance_only")
  expect_identical(make(0.4999,0.01)$loh_call,
    "low_copy_without_significant_imbalance")
  expect_identical(make(0.5,0.01)$loh_call,"loss_not_supported")
  inclusive <- base; inclusive$copy_number_comparison <- "<="
  expect_identical(make(0.5,0.0099,inclusive)$loh_call,"loss_supported")
  x <- make(0.5,0.0099)
  changed <- base; changed$copy_number_loss_threshold <- 0.51
  y <- make(0.5,0.0099,changed)
  expect_identical(y$loh_call,"loss_supported")
  expect_identical(x$raw_values,y$raw_values)
  expect_identical(x$evidence_id,y$evidence_id)
  expect_false(identical(Neoantimon2:::hla_loh_fingerprints(p,source,base),
    Neoantimon2:::hla_loh_fingerprints(p,source,changed)))
  hla_loh_validation_evidence("policy-boundary-validation.json",list(
    default_policy=base,cases=list(cn_04999_p_00099=make(0.4999,0.0099)$loh_call,
      cn_05000_p_00099=make(0.5,0.0099)$loh_call,
      cn_04999_p_00100=make(0.4999,0.01)$loh_call,
      cn_05000_p_00100=make(0.5,0.01)$loh_call),
    changed_policy=changed,changed_call=y$loh_call,
    raw_evidence_id_unchanged=identical(x$evidence_id,y$evidence_id),
    fingerprint_changed=!identical(Neoantimon2:::hla_loh_fingerprints(p,source,base),
      Neoantimon2:::hla_loh_fingerprints(p,source,changed))))
  missing <- row; missing$values$cn1 <- ""
  expect_identical(Neoantimon2:::hla_loh_record(p,source[[1L]],missing,1L,map,
    other,"potential_loss",base)$loh_call,"not_evaluable")
  missing$values$cn1 <- "0.2"; missing$values$pval <- ""
  expect_identical(Neoantimon2:::hla_loh_record(p,source[[1L]],missing,1L,map,
    other,"potential_loss",base)$loh_call,"not_evaluable")
  missing$values$pval <- "broken"
  expect_identical(Neoantimon2:::hla_loh_record(p,source[[1L]],missing,1L,map,
    other,"potential_loss",base)$loh_reason,"malformed_numeric")
})

test_that("mapping ambiguity, absent genotype and homozygosity remain explicit", {
  match <- Neoantimon2:::hla_loh_match
  genotype <- c("HLA-A*02:01:01","HLA-A*02:01:02","HLA-A*24:02")
  expect_identical(match("hla_a_02_01",genotype)$status,"ambiguous")
  expect_identical(match("HLA-A*02:01:01",genotype)$allele,"HLA-A*02:01:01")
  expect_identical(match("HLA-B*02:01",genotype)$status,"unresolved")
  expect_identical(match("HLA-A*99:01",genotype)$status,"unresolved")
  expect_identical(match("HLA-A*02:01",genotype,locus="HLA-B")$status,"unresolved")
  p <- hla_loh_test_project(genotypes=c("HLA-A*02:01","HLA-B*07:02"))
  evidence <- Neoantimon2:::hla_loh_construct(p,
    Neoantimon2:::hla_loh_source_files(p),p$config$hla_loh$policy)
  a <- Filter(function(x)identical(x$canonical_hla_allele,"HLA-A*02:01"),evidence)
  expect_length(a,1L)
  expect_identical(a[[1L]]$loh_call,"homozygous_not_evaluable")
  b <- Filter(function(x)identical(x$canonical_hla_allele,"HLA-B*07:02"),evidence)
  expect_identical(b[[1L]]$loh_call,"homozygous_not_evaluable")
  empty <- hla_loh_test_project(samples="region-1",genotypes="HLA-A*02:01")
  empty$config$hla_loh$sources <- list()
  absent <- Neoantimon2:::hla_loh_construct(empty,list(),empty$config$hla_loh$policy)
  expect_identical(absent[[1L]]$loh_call,"homozygous_not_evaluable")
  ambiguous <- hla_loh_test_project(genotypes=c("HLA-A*02:01:01",
    "HLA-A*02:01:02","HLA-A*24:02"))
  uncertain <- Neoantimon2:::hla_loh_construct(ambiguous,
    Neoantimon2:::hla_loh_source_files(ambiguous),ambiguous$config$hla_loh$policy)
  for (allele in c("HLA-A*02:01:01","HLA-A*02:01:02")) {
    record <- Filter(function(x) identical(x$canonical_hla_allele,allele),uncertain)
    expect_identical(record[[1L]]$loh_call,"allele_mapping_ambiguous")
  }
  row <- list(candidate_id="candidate:ambiguous",sample_id="region-1",
    analysis_context="somatic:patient-1",hla_class="I",hla_allele="HLA-A*02:01:01")
  expect_identical(Neoantimon2:::hla_loh_link(list(row),uncertain)[[1L]]$
    presenter_hla_loh_status,"allele_mapping_ambiguous")
  expect_identical(match("HLA-A*02:01",c("HLA-A*02:01N"))$status,"unresolved")
})

test_that("multiple regions remain separate and malformed layouts fail", {
  p <- hla_loh_test_project(samples=c("region-1","region-2"))
  p$config$hla_loh$sources[[2L]]$path <- hla_loh_fixture(
    "lohhla-variant.HLAlossPrediction_CI.txt")
  evidence <- Neoantimon2:::hla_loh_construct(p,
    Neoantimon2:::hla_loh_source_files(p),p$config$hla_loh$policy)
  a <- Filter(function(x)x$canonical_hla_allele=="HLA-A*02:01",evidence)
  expect_length(a,2L)
  expect_setequal(vapply(a,`[[`,"","loh_call"),c("loss_supported","loss_not_supported"))
  expect_false(identical(a[[1L]]$evidence_id,a[[2L]]$evidence_id))
  bad <- tempfile(fileext=".txt")
  writeLines(c("HLA_A_type1\tHLA_A_type2\tPval_unique\tLossAllele\tKeptAllele\tHLA_type1copyNum_withBAFBin\tHLA_type2copyNum_withBAFBin",
    "A*02:01\tA*24:02\t0.01\tA*02:01"),bad)
  expect_error(Neoantimon2:::hla_loh_parse(bad),"unsupported_format")
  homozygous <- tempfile(fileext=".HLAlossPrediction_CI.txt")
  writeLines("All HLA alleles homozygous. No LOH analysis can be performed.",homozygous)
  expect_length(Neoantimon2:::hla_loh_parse(homozygous)$rows,0L)
})

test_that("external execution stays unavailable without supplied dependencies", {
  cap <- Neoantimon2:::hla_loh_capability()
  expect_true(cap$parser_available)
  expect_false(cap$available)
  expect_identical(cap$execution_status,"not_available")
  expect_setequal(cap$supported_loci,c("HLA-A","HLA-B","HLA-C"))
  command <- Neoantimon2:::hla_loh_command(list(script="/absent/LOHHLAscript.R"),
    "sample-1")
  expect_identical(command$execution_status,"not_available")
  expect_false(isTRUE(command$executed))
  failed <- hla_loh_test_project(genotypes=c("HLA-A*02:01","HLA-A*24:02"))
  failed$config$hla_loh <- Neoantimon2:::hla_loh_config(list(sources=list(list(
    sample_id="region-1",provider_status="failed",failure_reason="external_exit_1",
    provider_version="unknown"))))
  records <- Neoantimon2:::hla_loh_construct(failed,
    Neoantimon2:::hla_loh_source_files(failed),failed$config$hla_loh$policy)
  expect_true(all(vapply(records,function(x)x$loh_call=="provider_failed",TRUE)))
  expect_identical(records[[1L]]$failure_reason,"external_exit_1")
})

test_that("paired BAM options reach the command and validate strictly", {
  root <- tempfile(); dir.create(root)
  for (name in c("LOHHLAscript.R","tumor.bam","normal.bam","hla.txt",
    "hla.fa","copy.txt")) file.create(file.path(root,name))
  novo <- file.path(root,"novo"); gatk <- file.path(root,"gatk")
  dir.create(novo); dir.create(gatk)
  for (name in c("novoalign","novoindex")) file.create(file.path(novo,name))
  for (name in c("SamToFastq.jar","SortSam.jar","FilterSamReads.jar"))
    file.create(file.path(gatk,name))
  execution <- list(script=file.path(root,"LOHHLAscript.R"),
    tumor_bam=file.path(root,"tumor.bam"),normal_bam=file.path(root,"normal.bam"),
    hla_calls=file.path(root,"hla.txt"),hla_fasta=file.path(root,"hla.fa"),
    copy_number=file.path(root,"copy.txt"),novoalign_dir=novo,gatk_dir=gatk,
    output_dir=file.path(root,"out"))
  fake_tool <- function(name) if (name=="jellyfish") "" else paste0("/mock/",name)
  testthat::local_mocked_bindings(hla_loh_tool=fake_tool,
    hla_loh_r_packages=function(rscript) list(status="satisfied",missing=character()),
    .package="Neoantimon2")
  command <- function(options=list()) Neoantimon2:::hla_loh_command(
    c(execution,list(options=options)),"region-1")
  has_flag <- function(x,flag,value) {
    at <- which(x$args==flag)
    length(at)==1L && identical(x$args[[at+1L]],value)
  }
  mapping_off <- command(list(mapping_step=FALSE))
  expect_true(has_flag(mapping_off,"--mappingStep","FALSE"))
  expect_false("--fishingStep" %in% mapping_off$args)
  fishing_off <- command(list(fishing_step=FALSE))
  expect_true(has_flag(fishing_off,"--fishingStep","FALSE"))
  expect_identical(fishing_off$execution_status,"command_ready",
    info=paste(unlist(fishing_off$capability$missing),collapse=", "))
  both_on <- command(list(mapping_step=TRUE,fishing_step=TRUE))
  expect_true(has_flag(both_on,"--mappingStep","TRUE"))
  expect_true(has_flag(both_on,"--fishingStep","TRUE"))
  expect_identical(both_on$execution_status,"not_available")
  expect_true("jellyfish" %in% unlist(both_on$capability$missing))
  expect_false("jellyfish" %in% unlist(fishing_off$capability$missing))
  expect_false("jellyfish" %in% unlist(command(list(mapping_step=FALSE))$capability$missing))
  expect_false(any(c("--mappingStep","--fishingStep") %in% command()$args))
  expect_false(identical(mapping_off$command_identity,fishing_off$command_identity))
  for (bad in list(0,1,"yes","FALSE",NA,c(TRUE,FALSE)))
    expect_error(Neoantimon2:::hla_loh_config(list(execution=c(execution,
      list(options=list(fishing_step=bad))))),"scalar logical")
  expect_error(Neoantimon2:::hla_loh_capability(c(execution,
    list(options=list(unknown=TRUE)))),"Unknown key")
  expect_error(Neoantimon2:::hla_loh_config(list(execution=c(execution,
    list(interface="bam_directory_v1",options=list(fishing_step=FALSE))))),
    "does not support execution options")
  expect_identical(command(list(fishing_step=FALSE))$capability$options$fishing_step,FALSE)
  testthat::local_mocked_bindings(hla_loh_r_packages=function(rscript)
    list(status="missing",missing="Biostrings"),.package="Neoantimon2")
  missing_r <- command(list(fishing_step=FALSE))
  expect_identical(missing_r$execution_status,"not_available")
  expect_true("R:Biostrings" %in% unlist(missing_r$capability$missing))
  testthat::local_mocked_bindings(hla_loh_r_packages=function(rscript)
    list(status="unresolved",missing=character()),.package="Neoantimon2")
  unresolved_r <- command(list(fishing_step=FALSE))
  expect_identical(unresolved_r$execution_status,"dependencies_unresolved")
  expect_true(is.na(unresolved_r$capability$execution_dependencies_satisfied))
})

test_that("external R packages are checked only for the current Rscript", {
  current <- file.path(R.home("bin"),"Rscript")
  if (file.exists(current)) {
    result <- Neoantimon2:::hla_loh_r_packages(current)
    required <- c("optparse","seqinr","Biostrings","beeswarm","zoo","Rsamtools")
    expect_setequal(result$missing,required[!vapply(required,requireNamespace,
      TRUE,quietly=TRUE)])
  }
  expect_identical(Neoantimon2:::hla_loh_r_packages("/other/Rscript")$status,
    "unresolved")
})

test_that("execution settings do not change imported output or priority identity", {
  p <- hla_loh_test_project()
  before <- Neoantimon2:::hla_loh_construct(p,
    Neoantimon2:::hla_loh_source_files(p),p$config$hla_loh$policy)
  rows <- list(list(candidate_id="candidate:1",sample_id="region-1",
    analysis_context="somatic:patient-1",hla_class="I",hla_allele="HLA-A*02:01",
    final_grade="A",rank_within_sample_class=1L,priority_origin_id="origin:1",
    generic_loh=list(status="no_evidence")))
  prior <- Neoantimon2:::hla_loh_link(rows,before)
  p$config$hla_loh$execution <- Neoantimon2:::hla_loh_config(list(
    execution=list(script="/absent/LOHHLAscript.R",
      options=list(mapping_step=FALSE,fishing_step=FALSE))))$execution
  expect_false(Neoantimon2:::hla_loh_capability(p$config$hla_loh$execution)$available)
  after <- Neoantimon2:::hla_loh_construct(p,
    Neoantimon2:::hla_loh_source_files(p),p$config$hla_loh$policy)
  linked <- Neoantimon2:::hla_loh_link(rows,after)
  expect_identical(vapply(after,`[[`,"","evidence_id"),
    vapply(before,`[[`,"","evidence_id"))
  expect_identical(linked,prior)
  expect_identical(linked[[1L]]$candidate_id,rows[[1L]]$candidate_id)
  expect_identical(rows[[1L]]$final_grade,"A")
  expect_identical(rows[[1L]]$rank_within_sample_class,1L)
})

test_that("candidate linkage is deterministic and preserves frozen priority scalars", {
  p <- hla_loh_test_project()
  evidence <- Neoantimon2:::hla_loh_construct(p,
    Neoantimon2:::hla_loh_source_files(p),p$config$hla_loh$policy)
  candidate <- function(id,allele,class="I") list(candidate_id=id,sample_id="region-1",
    analysis_context="somatic:patient-1",hla_class=class,hla_allele=allele,
    final_grade="A",rank_within_sample_class=1L,priority_origin_id="origin:1",
    generic_loh=list(status="no_evidence"))
  rows <- list(candidate("a","HLA-A*02:01"),candidate("b","HLA-A*24:02"),
    candidate("c","HLA-A*02:01"),candidate("d",NULL,"II"))
  before <- lapply(rows,function(x)x[c("candidate_id","final_grade",
    "rank_within_sample_class","priority_origin_id")])
  linked <- Neoantimon2:::hla_loh_link(rows,evidence)
  expect_identical(linked[[1L]]$presenter_escape_status,"hla_loss_supported")
  expect_identical(linked[[2L]]$presenter_escape_status,"hla_loss_not_supported")
  expect_identical(linked[[2L]]$paired_hla_loss_status,"paired_homolog_loss_supported")
  expect_identical(linked[[3L]]$presenter_hla_loh_evidence_id,
    linked[[1L]]$presenter_hla_loh_evidence_id)
  expect_identical(linked[[4L]]$presenter_hla_loh_status,"not_applicable")
  expect_identical(lapply(rows,function(x)x[c("candidate_id","final_grade",
    "rank_within_sample_class","priority_origin_id")]),before)
  reversed <- Neoantimon2:::hla_loh_link(rows,rev(evidence))
  expect_identical(linked,reversed)
  expect_identical(linked[[1L]]$generic_hla_loh_consistency,"not_comparable")
  hla_loh_validation_evidence("candidate-linkage-validation.json",list(
    fixture_kind="synthetic",links=linked,
    evidence_order_invariant=identical(linked,reversed),
    identity_and_priority_unchanged=identical(lapply(rows,function(x)
      x[c("candidate_id","final_grade","rank_within_sample_class",
        "priority_origin_id")]),before)))
})

test_that("real-reference candidates attach fixture evidence without identity or rank changes", {
  vcf <- tempfile(fileext=".vcf")
  writeLines(c("##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample-1",
    "chr1\t4\t.\tG\tA\t.\tPASS\t.\tGT\t0/1"),vcf)
  base <- reconstruction_test_project()
  cfg <- base$config; cfg$peptides <- peptide_test_policy(3L,5L)
  p <- neo_run(neo_project(samples=base$samples,config=cfg,
    hla=list("patient-1"=c("HLA-A*02:01","HLA-A*24:02","HLA-B*07:02",
      "HLA-B*08:01","HLA-DRA*01:01~HLA-DRB1*04:01")),
    inputs=list(list(input_id="hla-loh-reference-vcf",input_role="somatic_vcf",
      sample_id="sample-1",path=vcf))))
  baseline <- neo_results(p)
  p$config$hla_loh <- Neoantimon2:::hla_loh_config(list(sources=list(list(
    sample_id="sample-1",path=hla_loh_fixture("lohhla-standard.HLAlossPrediction_CI.xls"),
    provider_version="fixture-only"))))
  p <- neo_run(p)
  after <- neo_results(p)
  filtered <- neo_results(p,candidate_id=after$candidate_id[[1L]])
  expect_equal(nrow(filtered),1L)
  expect_identical(filtered$candidate_id[[1L]],after$candidate_id[[1L]])
  expect_identical(after$candidate_id,baseline$candidate_id)
  expect_identical(after$final_grade,baseline$final_grade)
  expect_identical(after$rank_within_sample_class,baseline$rank_within_sample_class)
  expect_true(any(after$presenter_escape_status=="hla_loss_supported"))
  expect_true(any(after$presenter_hla_loh_status=="not_applicable"))
  expect_true(all(after$presenter_hla_loh_status[after$hla_class=="II"]=="not_applicable"))
  one <- after$candidate_id[after$presenter_escape_status=="hla_loss_supported"][[1L]]
  explained <- neo_explain(p,candidate_id=one)$candidates[[1L]]$hla_loh
  expect_identical(explained$allele_evidence$provider_version,"fixture-only")
  expect_true(nchar(explained$allele_evidence$raw_output_checksum)==64L)
  expect_identical(explained$frozen_priority_note,
    "somatic-priority-v1 grade and rank unchanged")
  expect_identical(Neoantimon2:::hla_loh_read(p)$summary$execution_status,"not_available")
  hla_loh_validation_evidence("real-candidate-validation.json",list(
    fixture_kind="synthetic LOHHLA output on real-reference-derived candidates",
    candidates=nrow(after),class_i=sum(after$hla_class=="I"),
    class_ii=sum(after$hla_class=="II"),
    loss_supported_candidates=sum(after$presenter_escape_status=="hla_loss_supported"),
    ids_unchanged=identical(after$candidate_id,baseline$candidate_id),
    grades_unchanged=identical(after$final_grade,baseline$final_grade),
    ranks_unchanged=identical(after$rank_within_sample_class,
      baseline$rank_within_sample_class),
    provider_provenance=explained$allele_evidence[c("provider","provider_version",
      "raw_output_checksum","run_id","raw_row_index")],
    patient_observations_claimed=FALSE))
})

hla_loh_mock_predictor <- function(provider) {
  path <- tempfile(paste0("lohhla-mock-",provider,"-"))
  version <- if (provider=="netmhcpan") "NetMHCpan version 4.2" else
    "NetMHCIIpan version 4.3"
  writeLines(c("#!/bin/sh",paste0("if [ \"$1\" = \"-h\" ]; then echo '# ",
    version,"'; exit 0; fi"),"input=''", "allele=''", "while [ $# -gt 0 ]; do",
    "  case \"$1\" in -p|-f) input=$2; shift 2 ;; -a) allele=$2; shift 2 ;; *) shift ;; esac",
    "done", "if [ -z \"$input\" ]; then exit 2; fi",
    if (provider=="netmhcpan") "echo 'Pos HLA Peptide Score_EL %Rank_EL'" else
      "echo 'Pos MHC Peptide Score_EL %Rank_EL'",
    "n=0", "while IFS= read -r peptide; do", "  n=$((n+1))",
    "  case $((n%4)) in 1) rank=0.5 ;; 2) rank=0.5 ;; 3) rank=2.0 ;; 0) rank=6.0 ;; esac",
    "  printf '%s %s %s 0.5 %s\\n' \"$n\" \"$allele\" \"$peptide\" \"$rank\"",
    "done < \"$input\""),path)
  Sys.chmod(path,"0755")
  path
}

test_that("frozen gold priority and identity survive HLA LOH attachment", {
  vcf <- tempfile(fileext=".vcf")
  writeLines(c("##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample-1",
    "chrW\t4\t.\tA\tG\t.\tPASS\t.\tGT\t0/1",
    "chrW\t13\t.\tA\tG\t.\tPASS\t.\tGT\t0/1"),vcf)
  base <- reconstruction_test_project(); cfg <- base$config
  cfg$peptides <- peptide_test_policy(8L,9L)
  cfg$predictors <- list(list(provider="netmhcpan",task="mhc1_presentation",
    executable=hla_loh_mock_predictor("netmhcpan"),options=list()),
    list(provider="netmhciipan",task="mhc2_presentation",
      executable=hla_loh_mock_predictor("netmhciipan"),options=list()))
  inputs <- list(list(input_id="hla-loh-gold-vcf",input_role="somatic_vcf",
    sample_id="sample-1",path=vcf))
  hla <- list("patient-1"=c("HLA-A*02:01","HLA-DRA*01:01~HLA-DRB1*04:01"))
  pre <- neo_project(samples=base$samples,config=cfg,hla=hla,inputs=inputs)
  events <- Neoantimon2:::event_import_project(pre)
  dna <- tempfile(fileext=".tsv")
  writeLines(c("event_id\ttumor_ref_reads\ttumor_alt_reads\ttumor_depth\tnormal_ref_reads\tnormal_alt_reads\tnormal_depth",
    vapply(events,function(e) paste(e$event_id,"20","5","25","19","1","20",sep="\t"),"")),dna)
  rna <- tempfile(fileext=".tsv")
  writeLines(c("event_id\trna_ref_reads\trna_alt_reads\trna_depth",
    paste(events[[1L]]$event_id,"9","1","10",sep="\t"),
    paste(events[[2L]]$event_id,"10","0","10",sep="\t")),rna)
  facets <- tempfile(fileext=".tsv")
  writeLines(c("chrom\tloc.start\tloc.end\ttcn.em\tlcn.em\tnormal_cn\tpurity",
    "W\t1\t42\t2\t0\t2\t0.8"),facets)
  cfg$somatic_evidence <- list(sources=list(
    list(kind="dna_reads",sample_id="sample-1",path=dna),
    list(kind="rna_reads",sample_id="sample-1",path=rna),
    list(kind="facets_cn",sample_id="sample-1",path=facets,
      reference_build=cfg$scientific$genome_assembly)))
  p <- neo_run(neo_project(samples=base$samples,config=cfg,hla=hla,inputs=inputs))
  baseline <- neo_results(p)
  expect_equal(nrow(baseline),6L)
  expect_setequal(baseline$final_grade,c("A","A","B","B","C","E"))
  peptide_hash <- digest::digest(file=Neoantimon2:::peptide_stage_path(p),
    algo="sha256",serialize=FALSE)
  p$config$hla_loh <- Neoantimon2:::hla_loh_config(list(sources=list(list(
    sample_id="sample-1",path=hla_loh_fixture("lohhla-standard.HLAlossPrediction_CI.xls"),
    provider_version="synthetic-fixture"))))
  p <- neo_run(p)
  attached <- neo_results(p)
  for (field in c("candidate_id","final_grade","rank_within_sample_class",
    "priority_origin_id","peptide_key","mhc_key"))
    expect_identical(attached[[field]],baseline[[field]])
  expect_identical(digest::digest(file=Neoantimon2:::peptide_stage_path(p),
    algo="sha256",serialize=FALSE),peptide_hash)
  expect_true(any(attached$presenter_hla_loh_status=="homozygous_not_evaluable"))
  expect_true(all(attached$presenter_hla_loh_status[attached$hla_class=="II"]==
    "not_applicable"))
  expect_identical(Neoantimon2:::hla_loh_read(p)$summary$execution_status,
    "not_available")
  before <- Neoantimon2:::hla_loh_read(p)
  changed <- p
  changed$config$hla_loh$policy$copy_number_loss_threshold <- 0.4
  changed <- neo_run(changed)
  after <- Neoantimon2:::hla_loh_read(changed)
  expect_false(identical(before$fingerprints,after$fingerprints))
  expect_setequal(vapply(before$evidence,`[[`,"","evidence_id"),
    vapply(after$evidence,`[[`,"","evidence_id"))
  expect_identical(neo_results(changed)$final_grade,baseline$final_grade)
  hla_loh_validation_evidence("frozen-priority-invariance.json",list(
    fixture_kind="frozen gold workflow with synthetic LOHHLA output",
    candidates=nrow(baseline),grade_counts=as.list(table(baseline$final_grade)),
    candidate_ids_unchanged=identical(attached$candidate_id,baseline$candidate_id),
    peptide_keys_unchanged=identical(attached$peptide_key,baseline$peptide_key),
    presenter_keys_unchanged=identical(attached$mhc_key,baseline$mhc_key),
    grades_unchanged=identical(attached$final_grade,baseline$final_grade),
    ranks_unchanged=identical(attached$rank_within_sample_class,
      baseline$rank_within_sample_class),
    priority_origins_unchanged=identical(attached$priority_origin_id,
      baseline$priority_origin_id),
    peptide_stage_checksum_unchanged=identical(digest::digest(
      file=Neoantimon2:::peptide_stage_path(p),algo="sha256",serialize=FALSE),
      peptide_hash),
    policy_fingerprint_changed=!identical(before$fingerprints,after$fingerprints),
    raw_evidence_ids_unchanged=setequal(vapply(before$evidence,`[[`,"","evidence_id"),
      vapply(after$evidence,`[[`,"","evidence_id")),
    patient_observations_claimed=FALSE))
})
