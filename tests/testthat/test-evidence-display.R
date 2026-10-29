test_that("read-only DNA VAF projection preserves origin and measurement semantics", {
  f <- Neoantimon2:::input_display_candidates
  graph <- somatic_test_peptide()
  row <- graph$candidates[[1]]; row$grade <- "NOT_GRADED"; row$rank <- 1L
  row$origin_ids <- list("origin:a","origin:b")
  evidence <- function(value) list(list(event_id="event:a",sample_id="sample-1",
    evidence_id="evidence:a",tumor=list(vaf=value,provenance=list(source="reported pooled TVAF"))))
  out <- f(list(row),graph$origins,evidence(.2))[[1]]
  expect_equal(out$tumor_vaf,.2)
  expect_identical(out$grade,row$grade);expect_identical(out$rank,row$rank)
  expect_length(out$dna_vaf_display$origin_observations,2L)
  expect_identical(out$dna_vaf_display$origin_observations[[1]]$events[[1]]$observation$provenance$source,
    "reported pooled TVAF")
  expect_equal(f(list(row),graph$origins,evidence(0))[[1]]$tumor_vaf,0)
  expect_null(f(list(row),graph$origins,evidence(NULL))[[1]]$tumor_vaf)
  graph$origins[[2]]$source_event_ids <- list("event:b")
  ev <- c(evidence(.2),list(list(event_id="event:b",sample_id="sample-1",tumor=list(vaf=.4))))
  conflict <- f(list(row),graph$origins,ev)[[1]]
  expect_null(conflict$tumor_vaf)
  expect_true("tumor_vaf" %in% names(conflict))
  expect_identical(conflict$dna_vaf_display$status,"conflicting_origins")
  row$priority_origin_id <- "origin:b"
  expect_equal(f(list(row),graph$origins,ev)[[1]]$tumor_vaf,.4)
  graph$origins[[2]]$source_event_ids <- list("event:a","event:b")
  expect_null(f(list(row),graph$origins,ev)[[1]]$tumor_vaf)
  ev[[1]]$sample_id <- "other"
  expect_error(f(list(row),graph$origins,ev),"sample mismatch")
})

test_that("count-derived and multiallelic fractions retain their selected ALT provenance", {
  f <- Neoantimon2:::input_display_candidates
  graph <- somatic_test_peptide()
  for (i in 1:2) {
    e <- somatic_test_event();e$event_type <- "SNV";e$ref <- "G";e$alt <- c("A","C")[[i]]
    e$event_details_json <- list(source_records=list(list(format="VCF",input_id="fixture",line=5L,
      original=list(alts=list("A","C")),selected_alt_index=i,
      raw_format=list(keys="GT:AD:DP",sample="1/2:50,20,30:100"))))
    raw <- Neoantimon2:::somatic_raw_construct(list(e),graph,list())
    z <- f(graph$candidates,graph$origins,raw$event_evidence)
    expect_true(all(vapply(z,function(x)x$tumor_vaf==c(.2,.3)[i],TRUE)))
    expect_equal(z[[1]]$dna_vaf_display$origin_observations[[1]]$events[[1]]$observation$provenance$selected_alt_index,i)
  }
})

test_that("ordinary results reopen with VAF without changing stored Grade/rank", {
  base <- reconstruction_test_project()
  vcf <- tempfile(fileext=".vcf")
  writeLines(c("##fileformat=VCFv4.2","##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allelic depths>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample-1",
    "chr1\t4\t.\tG\tA\t.\tPASS\t.\tAD\t80,20"),vcf)
  cfg <- base$config;cfg$peptides <- peptide_test_policy(3L,5L)
  p <- neo_run(neo_project(samples=base$samples,config=cfg,hla=list("patient-1"="HLA-A*02:01"),
    inputs=list(list(input_id="display",input_role="somatic_vcf",sample_id="sample-1",path=vcf))))
  before <- Neoantimon2:::somatic_priority_original_read(p)
  result <- neo_results(p)
  expect_gt(nrow(result),0L);expect_true(all(result$tumor_vaf==.2))
  expect_identical(neo_results(neo_project(path=p$path)),result)
  expect_identical(Neoantimon2:::somatic_priority_original_read(p),before)
  stored <- if(is.null(before)) Neoantimon2:::read_runtime_snapshot(p)$neoantigen_results else before$rows
  expect_identical(result$final_grade,vapply(stored,`[[`,"","final_grade"))
  expect_true(all(vapply(result$dna_vaf_provenance,function(x)x$status=="observed",TRUE)))
})

test_that("gene annotation is carried from GTF to Origins without changing biological identity", {
  base <- reconstruction_test_project(); bundle <- Neoantimon2:::reference_bundle(base)
  original <- bundle$resources$annotation_gtf
  path <- tempfile(fileext=".gtf")
  writeLines(sub('gene_id "G1.1";', 'gene_id "G1.1"; gene_name "SYNTHETIC_GENE";',
    readLines(original$resolved_path)),path)
  bundle$resources$annotation_gtf$resolved_path <- path
  bundle$resources$annotation_gtf$integrity <- Neoantimon2:::integrity_record(path,"sha256",NULL)
  models <- Neoantimon2:::reconstruction_models(bundle,c("TA.1","TB.1"))
  expect_true(all(vapply(models,function(x)identical(x$gene_symbol,"SYNTHETIC_GENE"),TRUE)))
  expect_true(all(vapply(models,function(x)identical(x$gene_id,"G1.1"),TRUE)))
  ctx <- peptide_test_contexts(list(reconstruction_test_row("e1","chr1",4L,"G","A")),"TA.1")[[1]]
  event <- reconstruction_test_row("e1","chr1",4L,"G","A");event$event_key <- "key:e1"
  peptide <- Neoantimon2:::peptide_enumerate(ctx,peptide_test_policy())$peptides[[1]]
  old <- Neoantimon2:::peptide_origin(peptide,ctx,base,list(event))
  ctx$gene_symbol <- "SYNTHETIC_GENE"
  new <- Neoantimon2:::peptide_origin(peptide,ctx,base,list(event))
  expect_identical(new$origin_id,old$origin_id)
  expect_identical(new$gene_symbol,"SYNTHETIC_GENE")
  expect_identical(new$genes_json[[1]]$gene_id,"G1.1")
  expect_null(old$gene_symbol)
})

test_that("coordinate evidence validates one bundle per import and retains strict matching", {
  project <- reconstruction_test_project();project$path <- tempfile("coordinate-test-")
  path <- tempfile(fileext=".tsv")
  writeLines(c("contig\tpos\tref\talt\ttumor_vaf",rep("chr1\t4\tG\tA\t0.2",10)),path)
  source <- somatic_test_source("dna_reads",build="GRCh38.p14");source$resolved_path <- path
  event <- reconstruction_test_row("e1","chr1",4L,"G","A")
  calls <- 0L; original <- Neoantimon2:::reference_bundle
  testthat::local_mocked_bindings(reference_bundle=function(...) { calls <<- calls+1L;original(...) },.package="Neoantimon2")
  out <- Neoantimon2:::input_source_rows(list(source),list(event),project)
  expect_equal(calls,1L)
  expect_length(out[[1]]$rows,10L)
  expect_true(all(vapply(out[[1]]$rows,function(x)x$event_id=="e1" && x$tumor_vaf=="0.2",TRUE)))
  event$sample_id <- "different"
  expect_error(Neoantimon2:::input_source_rows(list(source),list(event),project),"exactly one sample Event")
  source$reference_build <- "GRCh37"
  expect_error(Neoantimon2:::input_source_rows(list(source),list(event),project),"reference_build")
})
