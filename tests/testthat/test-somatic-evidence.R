test_that("DNA counts, VAF, missingness, and configurable technical thresholds are separate", {
  event <- somatic_test_event()
  peptide <- somatic_test_peptide()
  source <- somatic_test_source("dna_reads",list(somatic_test_row()))
  raw <- Neoantimon2:::somatic_raw_construct(list(event),peptide,list(source))
  Neoantimon2:::somatic_validate_raw(raw,list(event),peptide)
  dna <- raw$event_evidence[[1L]]
  expect_identical(dna$tumor$ref_reads,5L)
  expect_identical(dna$tumor$alt_reads,5L)
  expect_identical(dna$tumor$depth,10L)
  expect_equal(dna$tumor$vaf,0.5)
  expect_identical(dna$normal$alt_reads,1L)
  expect_equal(dna$normal$vaf,0.05)
  expect_identical(dna$rna$status,"unavailable")
  policy <- Neoantimon2:::somatic_evidence_policy(list())$eligibility
  eligible <- Neoantimon2:::somatic_eligibility_construct(raw,policy)
  expect_true(all(vapply(eligible$candidate_eligibility,function(x)
    x$eligibility_status=="eligible",TRUE)))
  comparisons <- eligible$candidate_eligibility[[1L]]$origin_results[[1L]]$dna_criteria
  expect_setequal(vapply(comparisons,`[[`,"","evidence_field"),
    c("tumor_alt_reads","normal_alt_reads"))
  expect_true(all(vapply(comparisons,`[[`,TRUE,"satisfied")))
  expect_setequal(vapply(comparisons,`[[`,0L,"policy_threshold"),c(5L,1L))
  expect_null(eligible$candidate_eligibility[[1L]]$origin_results[[1L]]$rna_criteria[[1L]]$satisfied)
  expect_identical(eligible$candidate_eligibility[[1L]]$policy_version,"somatic-technical-v1")
  changed <- source; changed$rows <- list(somatic_test_row(tumor_alt="4",normal_alt="2",
    tumor_ref="6",normal_ref="18"))
  failed <- Neoantimon2:::somatic_eligibility_construct(
    Neoantimon2:::somatic_raw_construct(list(event),peptide,list(changed)),policy)
  reasons <- failed$candidate_eligibility[[1L]]$reason_records
  expect_identical(failed$candidate_eligibility[[1L]]$eligibility_status,"ineligible")
  expect_setequal(vapply(reasons,`[[`,"","reason_code"),
    c("insufficient_tumor_alt_reads","excessive_normal_alt_reads"))
  expect_setequal(vapply(reasons,`[[`,0L,"observed_value"),c(4L,2L))
  permissive <- policy; permissive$tumor_min_alt_reads <- 4L; permissive$normal_max_alt_reads <- 2L
  expect_identical(Neoantimon2:::somatic_eligibility_construct(
    Neoantimon2:::somatic_raw_construct(list(event),peptide,list(changed)),permissive)$candidate_eligibility[[1L]]$eligibility_status,
    "eligible")
  zero <- source; zero$rows <- list(somatic_test_row(tumor_alt="0",tumor_ref="0",depth="0"))
  zero_assay <- Neoantimon2:::somatic_raw_construct(list(event),peptide,list(zero))$event_evidence[[1L]]$tumor
  expect_identical(zero_assay$depth,0L)
  expect_null(zero_assay$vaf)
  expect_identical(zero_assay$status,"insufficient_depth")
  missing <- source; missing$rows <- list(somatic_test_row(tumor_alt=NULL,tumor_ref=NULL,depth=NULL))
  missing_assay <- Neoantimon2:::somatic_raw_construct(list(event),peptide,list(missing))$event_evidence[[1L]]$tumor
  expect_null(missing_assay$depth)
  expect_null(missing_assay$alt_reads)
  expect_identical(missing_assay$status,"not_evaluated")
  alt_only <- source; alt_only$rows <- list(somatic_test_row(tumor_alt="5",
    tumor_ref=NULL,depth=NULL))
  alt_record <- Neoantimon2:::somatic_raw_construct(list(event),peptide,list(alt_only))
  expect_identical(alt_record$event_evidence[[1L]]$tumor$alt_reads,5L)
  expect_null(alt_record$event_evidence[[1L]]$tumor$vaf)
  expect_identical(Neoantimon2:::somatic_eligibility_construct(alt_record,policy)$candidate_eligibility[[1L]]$eligibility_status,
    "eligible")
  canonical <- event; canonical$tumor_ref_count <- 5L; canonical$tumor_alt_count <- 5L
  canonical$tumor_raw_depth <- 10L; canonical$normal_ref_count <- 19L
  canonical$normal_alt_count <- 1L; canonical$normal_raw_depth <- 20L
  inherited <- Neoantimon2:::somatic_raw_construct(list(canonical),peptide,list())$event_evidence[[1L]]
  expect_equal(inherited$tumor$vaf,0.5)
  expect_identical(inherited$tumor$provenance$source,"frozen_event")
  inconsistent <- source; inconsistent$rows <- list(somatic_test_row(tumor_alt="8",tumor_ref="5",depth="10"))
  bad <- Neoantimon2:::somatic_raw_construct(list(event),peptide,list(inconsistent))
  expect_identical(bad$event_evidence[[1L]]$tumor$status,"invalid_counts")
  expect_null(bad$event_evidence[[1L]]$tumor$vaf)
})

test_that("optional RNA support and expression retain distinct statuses", {
  event <- somatic_test_event(); peptide <- somatic_test_peptide()
  dna <- somatic_test_source("dna_reads",list(somatic_test_row()))
  rna <- somatic_test_source("rna_reads",list(list(event_id="event:a",
    rna_ref_reads="9",rna_alt_reads="1",rna_depth="10")))
  expression <- somatic_test_source("expression",list(list(sample_id="sample-1",
    transcript_id="TX1",gene_id="GENE1",gene_tpm="3.5",transcript_tpm="1.2")))
  raw <- Neoantimon2:::somatic_raw_construct(list(event),peptide,list(dna,rna,expression))
  expect_identical(raw$event_evidence[[1L]]$rna$status,"observed_support")
  expect_equal(raw$event_evidence[[1L]]$rna$vaf,0.1)
  expect_equal(raw$expression_evidence[[1L]]$gene_tpm,3.5)
  expect_equal(raw$expression_evidence[[1L]]$transcript_tpm,1.2)
  zero_expression <- expression; zero_expression$rows[[1L]]$gene_tpm <- "0"
  zero_expression$rows[[1L]]$transcript_tpm <- "0"
  observed_zero <- Neoantimon2:::somatic_raw_construct(list(event),peptide,
    list(dna,zero_expression))$expression_evidence[[1L]]
  expect_identical(observed_zero$status,"available")
  expect_equal(observed_zero$gene_tpm,0)
  expect_length(raw$candidate_links[[1L]]$origin_paths[[1L]]$expression_evidence_ids,1L)
  expect_length(raw$candidate_links[[2L]]$origin_paths[[1L]]$expression_evidence_ids,0L)
  expect_identical(raw$candidate_links[[1L]]$origin_paths[[1L]]$expression_status,"available")
  expect_identical(raw$candidate_links[[2L]]$origin_paths[[1L]]$expression_status,"not_evaluated")
  criterion <- Neoantimon2:::somatic_eligibility_construct(raw,
    Neoantimon2:::somatic_evidence_policy(list())$eligibility)$candidate_eligibility[[1L]]$origin_results[[1L]]$rna_criteria[[1L]]
  expect_true(criterion$satisfied)
  expect_identical(criterion$policy_threshold,1L)
  rna$rows[[1L]]$rna_alt_reads <- "0"; rna$rows[[1L]]$rna_ref_reads <- "10"
  zero <- Neoantimon2:::somatic_raw_construct(list(event),peptide,list(dna,rna))
  expect_identical(zero$candidate_links[[1L]]$origin_paths[[1L]]$expression_status,"unavailable")
  expect_identical(zero$event_evidence[[1L]]$rna$status,"observed_no_alt")
  expect_identical(Neoantimon2:::somatic_eligibility_construct(zero,
    Neoantimon2:::somatic_evidence_policy(list())$eligibility)$candidate_eligibility[[1L]]$eligibility_status,
    "eligible")
  expect_false(Neoantimon2:::somatic_eligibility_construct(zero,
    Neoantimon2:::somatic_evidence_policy(list())$eligibility)$candidate_eligibility[[1L]]$origin_results[[1L]]$rna_criteria[[1L]]$satisfied)
  rna$rows[[1L]]$rna_ref_reads <- "0"; rna$rows[[1L]]$rna_depth <- "0"
  expect_identical(Neoantimon2:::somatic_raw_construct(list(event),peptide,list(dna,rna))$event_evidence[[1L]]$rna$status,
    "insufficient_depth")
})

test_that("FACETS segment mapping preserves boundaries, ambiguity, and build safety", {
  path <- testthat::test_path("..","fixtures","somatic","facets-cncf.tsv")
  source <- somatic_test_source("facets_cn",build="GRCh38")
  source$resolved_path <- path
  source$rows <- Neoantimon2:::somatic_read_table(source)
  segments <- Neoantimon2:::somatic_facets_segments(source)
  expect_length(segments,3L)
  expect_equal(segments[[1L]]$total_cn,2)
  expect_equal(segments[[1L]]$minor_cn,0)
  expect_equal(segments[[1L]]$major_cn,2)
  expect_equal(segments[[1L]]$cellular_fraction,1)
  expect_equal(Neoantimon2:::somatic_purity(somatic_test_event(),NULL,source)$value,0.8)
  expect_equal(segments[[1L]]$ploidy,2.1)
  expect_identical(Neoantimon2:::somatic_cn_map(somatic_test_event(pos=150L),segments,source)$segment_id,
    segments[[1L]]$segment_id)
  expect_identical(Neoantimon2:::somatic_cn_map(somatic_test_event(pos=151L),segments,source)$segment_id,
    segments[[2L]]$segment_id)
  expect_identical(Neoantimon2:::somatic_cn_map(somatic_test_event(pos=250L),segments,source)$status,
    "no_overlap")
  expect_identical(Neoantimon2:::somatic_cn_map(somatic_test_event(contig="chr7"),segments,source)$status,
    "contig_mismatch")
  expect_identical(Neoantimon2:::somatic_cn_map(somatic_test_event(build="GRCh37"),segments,source)$status,
    "build_incompatible")
  expect_identical(Neoantimon2:::somatic_cn_map(somatic_test_event(build=NULL),segments,source)$status,
    "build_unknown")
  ambiguous <- source
  ambiguous$resolved_path <- testthat::test_path("..","fixtures","somatic","facets-ambiguous.tsv")
  ambiguous$rows <- Neoantimon2:::somatic_read_table(ambiguous)
  hits <- Neoantimon2:::somatic_cn_map(somatic_test_event(pos=155L),
    Neoantimon2:::somatic_facets_segments(ambiguous),ambiguous)
  expect_identical(hits$status,"ambiguous_overlap")
  expect_null(hits$segment_id)
  expect_identical(Neoantimon2:::somatic_loh(somatic_test_event(pos=151L),NULL,NULL,
    Neoantimon2:::somatic_cn_map(somatic_test_event(pos=151L),segments,source))$status,
    "no_evidence")
})

test_that("purity, conditional CCF, multiplicity ambiguity, and generic LOH stay distinct", {
  event <- somatic_test_event(); peptide <- somatic_test_peptide()
  dna <- somatic_test_source("dna_reads",list(somatic_test_row()))
  purity <- somatic_test_source("purity",list(list(sample_id="sample-1",purity="0.8",
    ploidy="2.1",method="FACETS")))
  cn <- somatic_test_source("facets_cn",build="GRCh38")
  cn$resolved_path <- testthat::test_path("..","fixtures","somatic","facets-cncf.tsv")
  cn$rows <- Neoantimon2:::somatic_read_table(cn)
  event$mutation_multiplicity <- 2
  raw <- Neoantimon2:::somatic_raw_construct(list(event),peptide,list(dna,purity,cn))
  x <- raw$event_evidence[[1L]]
  expect_equal(x$purity$value,0.8)
  expect_equal(x$purity$ploidy,2.1)
  expect_identical(x$cn$status,"mapped")
  expect_identical(x$ccf$status,"conditional")
  expect_equal(x$ccf$raw_estimate,0.625)
  expect_equal(x$ccf$bounded_value,0.625)
  expect_equal(x$ccf$inputs$normal_cn,2)
  expect_identical(x$loh$status,"copy_neutral_loh")
  low_purity <- purity; low_purity$rows[[1L]]$purity <- "0.1"
  low_raw <- Neoantimon2:::somatic_raw_construct(list(event),peptide,list(dna,low_purity,cn))
  expect_identical(Neoantimon2:::somatic_eligibility_construct(low_raw,
    Neoantimon2:::somatic_evidence_policy(list())$eligibility)$candidate_eligibility[[1L]]$eligibility_status,
    "eligible")
  event$mutation_multiplicity <- NULL
  dna$rows <- list(somatic_test_row(tumor_alt="2",tumor_ref="8"))
  ambiguous <- Neoantimon2:::somatic_raw_construct(list(event),peptide,list(dna,purity,cn))$event_evidence[[1L]]$ccf
  expect_identical(ambiguous$status,"ambiguous")
  expect_length(ambiguous$solutions,2L)
  expect_null(ambiguous$bounded_value)
  supplied <- somatic_test_source("ccf",list(list(event_id="event:a",ccf="0.7")))
  external <- Neoantimon2:::somatic_raw_construct(list(event),peptide,list(dna,supplied))$event_evidence[[1L]]
  expect_identical(external$ccf$status,"external")
  expect_equal(external$ccf$bounded_value,0.7)
  expect_identical(external$purity$status,"unavailable")
  expect_identical(external$loh$status,"unavailable")
  missing <- Neoantimon2:::somatic_raw_construct(list(event),peptide,list(dna))$event_evidence[[1L]]$ccf
  expect_identical(missing$reason,"purity_unavailable")
  germline <- event; germline$event_role <- "germline"
  expect_identical(Neoantimon2:::somatic_event_record(germline,list(dna,purity,cn),
    list("sample-1"=Neoantimon2:::somatic_facets_segments(cn)))$ccf$reason,
    "non_somatic_event")
  no_normal <- cn; no_normal$rows <- lapply(cn$rows,function(r) {r$normal_cn <- NULL; r})
  without_normal <- Neoantimon2:::somatic_raw_construct(list(event),peptide,
    list(dna,purity,no_normal))$event_evidence[[1L]]$ccf
  expect_identical(without_normal$reason,"normal_cn_unavailable")
  impossible <- somatic_test_source("ccf",list(list(event_id="event:a",ccf_raw="1.2")))
  outside <- Neoantimon2:::somatic_raw_construct(list(event),peptide,list(dna,impossible))$event_evidence[[1L]]$ccf
  expect_equal(outside$raw_estimate,1.2)
  expect_null(outside$bounded_value)
  expect_identical(outside$status,"out_of_range")
  loh <- somatic_test_source("loh",list(list(event_id="event:a",status="allele_specific_loss")))
  expect_identical(Neoantimon2:::somatic_raw_construct(list(event),peptide,list(dna,loh))$event_evidence[[1L]]$loh$status,
    "allele_specific_loss")
})

test_that("event evidence is canonical once across candidates and compound origins", {
  events <- list(somatic_test_event(),somatic_test_event("event:b",pos=151L))
  peptide <- somatic_test_peptide(compound=TRUE)
  dna <- somatic_test_source("dna_reads",list(somatic_test_row(),
    somatic_test_row("event:b",tumor_alt="4",tumor_ref="6",normal_alt="2",normal_ref="18")))
  raw <- Neoantimon2:::somatic_raw_construct(events,peptide,list(dna))
  expect_length(raw$event_evidence,2L)
  expect_length(raw$candidate_links,2L)
  first <- raw$candidate_links[[1L]]$origin_paths[[1L]]
  compound <- raw$candidate_links[[2L]]$origin_paths[[1L]]
  expect_length(first$event_evidence_ids,1L)
  expect_length(compound$event_evidence_ids,2L)
  expect_identical(first$event_evidence_ids[[1L]],compound$event_evidence_ids[[1L]])
  expect_setequal(unlist(compound$somatic_event_ids),c("event:a","event:b"))
  expect_length(compound$reconstruction_ids,1L)
  policy <- Neoantimon2:::somatic_evidence_policy(list())$eligibility
  result <- Neoantimon2:::somatic_eligibility_construct(raw,policy)
  expect_identical(result$candidate_eligibility[[1L]]$eligibility_status,"eligible")
  expect_identical(result$candidate_eligibility[[2L]]$eligibility_status,"ineligible")
  expect_length(result$candidate_eligibility[[2L]]$reason_records,2L)
  expect_false(any(c("final_grade","rank","composite_score") %in% names(result)))
  expect_identical(raw$event_evidence[[1L]]$evidence_id,
    Neoantimon2:::somatic_raw_construct(events,peptide,list(dna))$event_evidence[[1L]]$evidence_id)
  expect_identical(raw$source_provenance[[1L]]$sha256,dna$sha256)
  alternate <- peptide
  alternate$candidates[[1L]]$origin_ids <- list("origin:a","origin:b")
  alt_raw <- Neoantimon2:::somatic_raw_construct(events,alternate,list(dna))
  alt_result <- Neoantimon2:::somatic_eligibility_construct(alt_raw,policy)
  expect_identical(alt_result$candidate_eligibility[[1L]]$eligibility_status,"eligible")
  expect_length(alt_result$candidate_eligibility[[1L]]$origin_results,2L)
})

test_that("raw checkpoint identity excludes eligibility thresholds", {
  root <- tempfile("somatic-fingerprint-")
  dir.create(file.path(root,"intermediate"),recursive=TRUE)
  writeLines("phase",file.path(root,"intermediate","phase_resolution.json"))
  writeLines("peptide",file.path(root,"intermediate","peptide_candidates.rds"))
  writeLines("raw",file.path(root,"intermediate","somatic_evidence_raw.rds"))
  project <- list(path=root)
  source <- somatic_test_source("dna_reads",list())
  raw_before <- Neoantimon2:::somatic_raw_fingerprints(project,list(source))
  p1 <- Neoantimon2:::somatic_evidence_policy(list())$eligibility
  p2 <- p1; p2$tumor_min_alt_reads <- 6L
  raw_after <- Neoantimon2:::somatic_raw_fingerprints(project,list(source))
  expect_identical(raw_before,raw_after)
  expect_false(identical(Neoantimon2:::somatic_eligibility_fingerprints(project,p1),
    Neoantimon2:::somatic_eligibility_fingerprints(project,p2)))
  source$sha256 <- paste(rep("b",64),collapse="")
  expect_false(identical(raw_before,Neoantimon2:::somatic_raw_fingerprints(project,list(source))))
})

test_that("source checksum and policy validation are explicit", {
  path <- tempfile(fileext=".tsv")
  writeLines(c("event_id\ttumor_alt_reads","event:a\t5"),path)
  sha <- digest::digest(file=path,algo="sha256",serialize=FALSE)
  policy <- Neoantimon2:::somatic_evidence_policy(list(sources=list(list(
    kind="dna_reads",sample_id="sample-1",path=path,checksum=sha))))
  project <- list(config=list(somatic_evidence=policy,paths=list(input_root=tempdir())),
    samples=list(list(sample_id="sample-1")))
  bound <- Neoantimon2:::somatic_source_files(project)[[1L]]
  expect_identical(bound$sha256,sha)
  expect_identical(bound$source,"dna_reads")
  writeLines(c("event_id\ttumor_alt_reads","event:a\t4"),path)
  expect_error(Neoantimon2:::somatic_source_files(project),"checksum mismatch")
  expect_error(Neoantimon2:::somatic_evidence_policy(list(eligibility=list(
    tumor_min_alt_reads=-1))),"threshold")
})

test_that("configured somatic sidecar preserves frozen candidate IDs and predictor sidecar", {
  vcf <- tempfile(fileext=".vcf")
  writeLines(c("##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample-1",
    "chr1\t4\t.\tG\tA\t.\tPASS\t.\tGT\t0/1"),vcf)
  base <- reconstruction_test_project()
  cfg <- base$config
  cfg$peptides <- peptide_test_policy()
  hla <- list("patient-1"="HLA-A*02:01")
  inputs <- list(list(input_id="somatic-phase-vcf",input_role="somatic_vcf",
    sample_id="sample-1",path=vcf))
  original <- neo_project(samples=base$samples,config=cfg,hla=hla,inputs=inputs)
  event <- Neoantimon2:::event_import_project(original)[[1L]]
  dna <- tempfile(fileext=".tsv")
  writeLines(c("event_id\ttumor_ref_reads\ttumor_alt_reads\ttumor_depth\tnormal_ref_reads\tnormal_alt_reads\tnormal_depth",
    paste(event$event_id,"5","5","10","19","1","20",sep="\t")),dna)
  rna <- tempfile(fileext=".tsv")
  writeLines(c("event_id\trna_ref_reads\trna_alt_reads\trna_depth",
    paste(event$event_id,"9","1","10",sep="\t")),rna)
  facets <- tempfile(fileext=".tsv")
  writeLines(c("chrom\tloc.start\tloc.end\ttcn.em\tlcn.em\tnormal_cn",
    "1\t1\t10\t2\t0\t2"),facets)
  purity <- tempfile(fileext=".tsv")
  writeLines(c("sample_id\tpurity\tploidy\tmethod","sample-1\t0.8\t2.0\tFACETS"),purity)
  cfg$somatic_evidence <- list(sources=list(
    list(kind="dna_reads",sample_id="sample-1",path=dna),
    list(kind="rna_reads",sample_id="sample-1",path=rna),
    list(kind="facets_cn",sample_id="sample-1",path=facets,
      reference_build=original$config$scientific$genome_assembly),
    list(kind="purity",sample_id="sample-1",path=purity)))
  cfg$predictors <- list(list(provider="netmhcpan",task="mhc1_binding",
    executable=file.path(tempdir(),"absent-netmhcpan"),options=list()))
  with_evidence <- neo_project(samples=base$samples,config=cfg,hla=hla,inputs=inputs)
  expect_identical(original$analysis_id,with_evidence$analysis_id)
  original <- neo_run(original)
  with_evidence <- neo_run(with_evidence)
  before <- Neoantimon2:::read_runtime_snapshot(original)
  after <- Neoantimon2:::read_runtime_snapshot(with_evidence)
  expect_identical(vapply(before$neoantigen_results,`[[`,"","candidate_id"),
    vapply(after$neoantigen_results,`[[`,"","candidate_id"))
  sidecar <- Neoantimon2:::somatic_read_sidecar(with_evidence)
  expect_length(sidecar$raw$event_evidence,1L)
  expect_equal(sidecar$raw$event_evidence[[1L]]$tumor$vaf,0.5)
  expect_true(all(vapply(sidecar$eligibility$candidate_eligibility,function(x)
    x$eligibility_status=="eligible",TRUE)))
  expect_true(all(vapply(after$neoantigen_results,function(x) x$final_grade=="NOT_GRADED",TRUE)))
  cps <- Neoantimon2:::read_manifest(with_evidence$path)$checkpoints
  expect_true(all(c("somatic_evidence_raw","somatic_eligibility") %in%
    vapply(cps,`[[`,"","step_id")))
  expect_true(any(startsWith(vapply(cps,`[[`,"","step_id"),"prediction_")))
  candidate_id <- after$neoantigen_results[[1L]]$candidate_id
  explained <- neo_explain(with_evidence,candidate_id=candidate_id)$candidates[[1L]]
  expect_identical(explained$somatic_evidence$eligibility$candidate_id,candidate_id)
  expect_length(explained$somatic_evidence$event_evidence,1L)
})
