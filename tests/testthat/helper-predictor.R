predictor_test_candidate <- function(id, sequence="ACDEFGHIK", presenter="HLA-A*02:01") {
  ii <- grepl("~",presenter,fixed=TRUE)
  chains <- strsplit(presenter,"~",fixed=TRUE)[[1L]]
  list(candidate_id=id,sample_peptide_id=paste0("sample-peptide:",id),
    ligand_sequence=sequence,mhc_key=presenter,hla_class=if (ii) "II" else "I",
    hla_allele=if (ii) NULL else presenter,
    hla_alpha=if (ii) chains[[1L]] else NULL,
    hla_beta=if (ii) chains[[2L]] else NULL,
    species="NCBITaxon:9606",analysis_context="somatic:patient-1",sample_id="sample-1")
}

predictor_test_capability <- function(provider) {
  x <- Neoantimon2:::predictor_capability(provider)
  x$available <- TRUE; x$availability_status <- "available"
  x$predictor_version <- switch(provider,netmhcpan="4.1",netmhciipan="4.3k",
    mhcflurry="2.0.1",fixture="fixture-1")
  x$executable_path <- if (provider=="fixture") NULL else Sys.which("true")
  x
}

predictor_test_fixture <- function(name) paste(readLines(testthat::test_path("..","fixtures",
  "predictor",name),warn=FALSE),collapse="\n")

predictor_test_project <- function() {
  vcf <- tempfile(fileext=".vcf")
  writeLines(c("##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample-1",
    "chr1\t4\t.\tG\tA\t.\tPASS\t.\tGT\t0/1"),vcf)
  base <- reconstruction_test_project()
  cfg <- base$config
  cfg$peptides <- Neoantimon2:::peptide_policy(list(class_i_lengths=3L,class_ii_lengths=5L))
  cfg$predictors <- list(list(provider="netmhcpan",task="mhc1_binding",
    executable=file.path(tempdir(),"absent-netmhcpan"),options=list()))
  neo_project(samples=base$samples,config=cfg,hla=list("patient-1"="HLA-A*02:01"),
    inputs=list(list(input_id="predictor-fixture",input_role="somatic_vcf",
      sample_id="sample-1",path=vcf)))
}
