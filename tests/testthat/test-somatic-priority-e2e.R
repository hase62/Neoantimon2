priority_mock_executable <- function(provider) {
  path <- tempfile(paste0("mock-",provider,"-"))
  version <- if (provider=="netmhcpan") "NetMHCpan version 4.2" else
    "NetMHCIIpan version 4.3"
  writeLines(c("#!/bin/sh",
    paste0("if [ \"$1\" = \"-h\" ]; then echo '# ",version,"'; exit 0; fi"),
    "input=''", "allele=''", "while [ $# -gt 0 ]; do",
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

test_that("gold project runs through final result and explanation artifacts", {
  vcf <- tempfile(fileext=".vcf")
  writeLines(c("##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample-1",
    "chrW\t4\t.\tA\tG\t.\tPASS\t.\tGT\t0/1",
    "chrW\t13\t.\tA\tG\t.\tPASS\t.\tGT\t0/1"),vcf)
  base <- reconstruction_test_project()
  cfg <- base$config
  cfg$peptides <- peptide_test_policy(8L,9L)
  cfg$predictors <- list(
    list(provider="netmhcpan",task="mhc1_presentation",
      executable=priority_mock_executable("netmhcpan"),options=list()),
    list(provider="netmhciipan",task="mhc2_presentation",
      executable=priority_mock_executable("netmhciipan"),options=list()))
  inputs <- list(list(input_id="priority-gold-vcf",input_role="somatic_vcf",
    sample_id="sample-1",path=vcf))
  hla <- list("patient-1"=c("HLA-A*02:01","HLA-DRA*01:01~HLA-DRB1*04:01"))
  pre <- neo_project(samples=base$samples,config=cfg,hla=hla,inputs=inputs)
  events <- Neoantimon2:::event_import_project(pre)
  expect_length(events,2L)
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
  p <- neo_project(samples=base$samples,config=cfg,hla=hla,inputs=inputs)
  p <- neo_run(p)
  final <- neo_results(p)
  expect_gt(nrow(final),0L)
  expect_true(all(c("A","B","C","E") %in% final$final_grade))
  expect_true(all(final$final_grade %in% c("A","B","C","D","E","INVALID","NOT_GRADED")))
  expect_true(any(final$primary_prediction_status=="evaluated"))
  expect_true(all(c("I","II") %in% final$hla_class))
  expect_true(all(final$policy_version=="somatic-priority-v1"))
  expect_true(all(final$priority_origin_selection_policy==
    "best-supported-eligible-origin-v1"))
  expect_true(all(!is.na(final$priority_origin_id)))
  expect_true(all(vapply(final$priority_source_event_ids,length,0L)>=1L))
  expect_true(any(final$ccf_status=="ambiguous"))
  expect_true(all(is.na(final$ccf_value[final$ccf_status=="ambiguous"])))
  expect_setequal(final$candidate_id,Neoantimon2:::read_runtime_snapshot(p)$neoantigen_results |>
    lapply(`[[`,"candidate_id") |> unlist())
  explanation <- neo_explain(p,candidate_id=final$candidate_id[[1L]])
  expect_identical(explanation$candidates[[1L]]$priority$final_grade,
    final$final_grade[[1L]])
  expect_gte(length(explanation$candidates[[1L]]$events),1L)
  expect_identical(explanation$candidates[[1L]]$priority$priority_origin_id,
    final$priority_origin_id[[1L]])
  expect_identical(explanation$candidates[[1L]]$priority$priority_source_event_ids,
    final$priority_source_event_ids[[1L]])
  expect_identical(explanation$candidates[[1L]]$priority$
    priority_origin_selection_policy,"best-supported-eligible-origin-v1")
  expect_length(explanation$candidates[[1L]]$priority$origin_priority_summaries,
    length(explanation$candidates[[1L]]$origins))
  expect_gte(length(explanation$candidates[[1L]]$somatic_evidence$event_evidence),1L)
  materialized <- jsonlite::read_json(file.path(p$path,"results","candidates.json"))
  expect_true(all(vapply(materialized,function(x)
    !is.null(x$priority_origin_id) && length(x$origin_priority_summaries)>0L,TRUE)))
  expect_match(neo_explain(p,candidate_id=final$candidate_id[[1L]],format="text"),
    "somatic-priority-v1")
  for (file in Neoantimon2:::somatic_priority_files())
    expect_true(file.exists(file.path(p$path,"results",file)))
  expect_true(any(vapply(Neoantimon2:::read_manifest(p$path)$checkpoints,
    function(x)x$step_id=="somatic_priority" && x$status=="COMPLETED",TRUE)))
  raw_before <- digest::digest(file=Neoantimon2:::somatic_raw_path(p),
    algo="sha256",serialize=FALSE)
  changed <- p
  changed$config$somatic_priority <- Neoantimon2:::somatic_priority_policy(list(
    thresholds=list(I=list(strong=.4))))
  changed <- neo_run(changed)
  changed_results <- neo_results(changed)
  expect_setequal(changed_results$candidate_id,final$candidate_id)
  expect_false(identical(changed_results$final_grade,final$final_grade))
  expect_identical(digest::digest(file=Neoantimon2:::somatic_raw_path(p),
    algo="sha256",serialize=FALSE),raw_before)
  variant_fp <- Neoantimon2:::somatic_priority_fingerprints(changed,
    changed$config$somatic_priority,Neoantimon2:::somatic_prediction_runs(changed))
  expect_true(file.exists(file.path(Neoantimon2:::somatic_priority_variant_dir(changed,
    variant_fp),"checkpoint.json")))

  missing_dna <- cfg
  missing_dna$somatic_evidence$sources <- Filter(function(x)x$kind!="dna_reads",
    cfg$somatic_evidence$sources)
  p_missing_dna <- neo_run(neo_project(samples=base$samples,config=missing_dna,
    hla=hla,inputs=inputs))
  missing_dna_results <- neo_results(p_missing_dna)
  expect_true(all(missing_dna_results$final_grade=="NOT_GRADED"))
  expect_true(all(missing_dna_results$technical_eligibility=="not_evaluated"))

  unavailable <- cfg
  unavailable$predictors <- list(list(provider="netmhcpan",task="mhc1_presentation",
    executable=file.path(tempdir(),"absent-netmhcpan"),options=list()))
  p_unavailable <- neo_run(neo_project(samples=base$samples,config=unavailable,
    hla=hla,inputs=inputs))
  unavailable_results <- neo_results(p_unavailable)
  expect_true(all(unavailable_results$final_grade=="NOT_GRADED"))
  expect_true(any(unavailable_results$technical_eligibility=="eligible"))
  expect_true(any(unavailable_results$primary_prediction_status=="unavailable"))

  bad_dna <- tempfile(fileext=".tsv")
  writeLines(c("event_id\ttumor_ref_reads\ttumor_alt_reads\ttumor_depth\tnormal_ref_reads\tnormal_alt_reads\tnormal_depth",
    vapply(events,function(e)paste(e$event_id,"20","4","24","18","2","20",sep="\t"),"")),bad_dna)
  invalid <- cfg
  for (i in seq_along(invalid$somatic_evidence$sources))
    if (invalid$somatic_evidence$sources[[i]]$kind=="dna_reads")
      invalid$somatic_evidence$sources[[i]]$path <- bad_dna
  p_invalid <- neo_run(neo_project(samples=base$samples,config=invalid,
    hla=hla,inputs=inputs))
  invalid_results <- neo_results(p_invalid)
  expect_true(all(invalid_results$final_grade=="INVALID"))
  expect_true(all(is.na(invalid_results$rank_within_sample_class)))

  missing_rna <- cfg
  missing_rna$somatic_evidence$sources <- Filter(function(x)x$kind!="rna_reads",
    cfg$somatic_evidence$sources)
  p_missing_rna <- neo_run(neo_project(samples=base$samples,config=missing_rna,
    hla=hla,inputs=inputs))
  expect_true(all(neo_results(p_missing_rna)$rna_status=="unavailable"))
  expect_true(any(final$rna_status=="observed_no_alt"))

  sec_vcf <- tempfile(fileext=".vcf")
  writeLines(c("##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample-1",
    "chrS2\t4\t.\tA\tG\t.\tPASS\t.\tGT\t0/1"),sec_vcf)
  sec_inputs <- list(list(input_id="priority-sec-vcf",input_role="somatic_vcf",
    sample_id="sample-1",path=sec_vcf))
  sec_pre <- neo_project(samples=base$samples,config=cfg,hla=hla,inputs=sec_inputs)
  sec_event <- Neoantimon2:::event_import_project(sec_pre)[[1L]]
  sec_dna <- tempfile(fileext=".tsv")
  writeLines(c("event_id\ttumor_ref_reads\ttumor_alt_reads\ttumor_depth\tnormal_ref_reads\tnormal_alt_reads\tnormal_depth",
    paste(sec_event$event_id,"20","5","25","19","1","20",sep="\t")),sec_dna)
  sec_cfg <- cfg
  sec_cfg$somatic_evidence$sources <- list(list(kind="dna_reads",
    sample_id="sample-1",path=sec_dna))
  p_sec <- neo_run(neo_project(samples=base$samples,config=sec_cfg,
    hla=hla,inputs=sec_inputs))
  sec_results <- neo_results(p_sec)
  sec_rows <- sec_results[grepl("U",sec_results$ligand_sequence,fixed=TRUE),,drop=FALSE]
  expect_gt(nrow(sec_rows),0L)
  expect_true(all(sec_rows$technical_eligibility=="eligible"))
  expect_true(all(sec_rows$final_grade=="NOT_GRADED"))
  expect_true(all(sec_rows$primary_prediction_status=="unsupported"))

  phased_vcf <- tempfile(fileext=".vcf")
  writeLines(c("##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>",
    "##FORMAT=<ID=PS,Number=1,Type=String,Description=Phase set>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample-1",
    "chrW\t4\t.\tA\tG\t.\tPASS\t.\tGT:PS\t0|1:gold",
    "chrW\t13\t.\tA\tG\t.\tPASS\t.\tGT:PS\t0|1:gold"),phased_vcf)
  phased_inputs <- list(list(input_id="priority-phased-vcf",input_role="somatic_vcf",
    sample_id="sample-1",path=phased_vcf))
  phased_pre <- neo_project(samples=base$samples,config=cfg,hla=hla,inputs=phased_inputs)
  phased_events <- Neoantimon2:::event_import_project(phased_pre)
  phased_dna <- tempfile(fileext=".tsv")
  writeLines(c("event_id\ttumor_ref_reads\ttumor_alt_reads\ttumor_depth\tnormal_ref_reads\tnormal_alt_reads\tnormal_depth",
    vapply(phased_events,function(e) paste(e$event_id,"20","5","25","19","1","20",sep="\t"),"")),phased_dna)
  phased_cfg <- cfg
  for (i in seq_along(phased_cfg$somatic_evidence$sources))
    if (phased_cfg$somatic_evidence$sources[[i]]$kind=="dna_reads")
      phased_cfg$somatic_evidence$sources[[i]]$path <- phased_dna
  p_phased <- neo_run(neo_project(samples=base$samples,config=phased_cfg,
    hla=hla,inputs=phased_inputs))
  compound <- Filter(function(x)length(x$source_event_ids)>1L,
    Neoantimon2:::read_runtime_snapshot(p_phased)$peptide_origins)
  expect_gt(length(compound),0L)
  expect_true(any(vapply(neo_results(p_phased)$source_event_ids,length,0L)>1L))
  compound_rows <- Filter(function(x)length(x$priority_source_event_ids)>1L,
    Neoantimon2:::somatic_priority_read(p_phased)$rows)
  expect_gt(length(compound_rows),0L)
  expect_true(any(vapply(compound_rows,function(x)
    identical(x$rna_status,"partial_support") && !x$rna_support,TRUE)))
  expect_true(all(vapply(compound_rows,function(x)
    x$final_grade %in% c("B","D","E","NOT_GRADED"),TRUE)))

  evidence_dir <- Sys.getenv("NEO_FINAL_E2E_EVIDENCE_DIR","")
  if (nzchar(evidence_dir)) {
    dir.create(evidence_dir,recursive=TRUE,showWarnings=FALSE)
    file.copy(file.path(p$path,"results",Neoantimon2:::somatic_priority_files()),
      evidence_dir,overwrite=TRUE)
    jsonlite::write_json(list(status="PASS",candidates=nrow(final),
      grades=as.list(table(final$final_grade)),
      classes=as.list(table(final$hla_class)),
      invalid_candidates=nrow(invalid_results),
      not_graded_missing_dna=nrow(missing_dna_results),
      not_graded_unavailable_predictor=nrow(unavailable_results),
      missing_rna_candidates=nrow(neo_results(p_missing_rna)),
      unsupported_sec_candidates=nrow(sec_rows),
      compound_origins=length(compound),
      candidate_ids=as.list(sort(final$candidate_id))),
      file.path(evidence_dir,"gold-e2e-summary.json"),auto_unbox=TRUE,pretty=TRUE)
  }
})
