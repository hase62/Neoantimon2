input_fixture_event <- function(values="0/1:80,20:100",keys="GT:AD:DP",alts=list("A"),index=1L) {
  event <- somatic_test_event();event$ref <- "G";event$alt <- alts[[index]];event$event_type <- "SNV"
  event$event_details_json <- list(source_records=list(list(format="VCF",input_id="fixture",
    line=5L,original=list(alts=alts),selected_alt_index=index,
    raw_format=list(keys=keys,sample=values),evidence_format_definitions=list(
      AF=list(Number="A",Description="Allele fraction in the sample")))))
  event
}

test_that("tumor FORMAT counts/estimated fractions use correct ALT and denominator", {
  f <- Neoantimon2:::input_vcf_assay
  e <- input_fixture_event();x <- f(e)
  expect_equal(x$vaf,.2);expect_identical(x$alt_reads,20L);expect_identical(x$vaf_denominator,100L)
  multi <- input_fixture_event("1/2:50,20,30:100",alts=list("A","C"),index=2L)
  expect_equal(f(multi)$vaf,.3);expect_identical(f(multi)$other_reads,20L)
  multi$event_details_json$source_records[[1]]$selected_alt_index <- 1L
  expect_equal(f(multi)$vaf,.2)
  expect_equal(f(input_fixture_event("20:80","AO:RO"))$vaf,.2)
  partial <- f(input_fixture_event("50,20,.","AD",alts=list("A","C")))
  expect_equal(partial$alt_reads,20L);expect_null(partial$vaf);expect_null(partial$vaf_denominator)
  expect_error(f(input_fixture_event(".2:0","AF:DP")),"zero")
  expect_null(f(input_fixture_event("100","DP")))
  pop <- input_fixture_event("100","DP");pop$event_details_json$source_records[[1]]$raw_info <- "AF=0.9"
  expect_null(f(pop))
  expect_equal(f(input_fixture_event("0,20","AD"))$vaf,1)
  expect_equal(f(input_fixture_event("20,0","AD"))$vaf,0)
  expect_null(f(input_fixture_event("0,0:0","AD:DP"))$vaf)
  expect_identical(f(input_fixture_event("0,0:0","AD:DP"))$status,"insufficient_depth")
  af <- f(input_fixture_event(".3","AF"));expect_equal(af$vaf,.3);expect_identical(af$status,"fraction_only")
  mismatch <- f(input_fixture_event("80,20:.3","AD:AF"))
  expect_equal(mismatch$vaf,.2);expect_equal(mismatch$reported_vaf,.3);expect_true(mismatch$provenance$discrepancy)
  expect_error(f(input_fixture_event("80,20:90","AD:DP")),"exceed")
  expect_error(f(input_fixture_event("-1,20","AD")),"integer")
  expect_error(f(input_fixture_event("Inf","AF")),"numeric")
  expect_error(f(input_fixture_event("1.1","AF")),"numeric")
  e$event_role <- "germline";expect_null(f(e))
})

test_that("explicit evidence retains precedence and all alternative provenance", {
  event <- input_fixture_event()
  source <- somatic_test_source("dna_reads")
  x <- Neoantimon2:::input_assay(event,list(tumor_ref_reads="50",tumor_alt_reads="50"),"tumor",source)
  expect_equal(x$vaf,.5);expect_equal(x$alternatives$vcf$vaf,.2);expect_true(x$source_disagreement)
  expect_equal(Neoantimon2:::input_assay(event,NULL,"tumor",NULL)$vaf,.2)
  partial <- input_fixture_event("50,20,.","AD",alts=list("A","C"))
  source$bam_observations <- list();source$bam_observations[[partial$event_id]] <- list(tumor=list(vaf=.5,status="observed_support"))
  expect_equal(Neoantimon2:::input_assay(partial,NULL,"tumor",source)$vaf,.5)
})

test_that("CIGAR haplotype classification handles indels, splicing and quality", {
  f <- Neoantimon2:::input_bam_haplotype
  expect_equal(f(1L,"5M","ACGTA","IIIII",1L,5L,20L),"ACGTA")
  expect_equal(f(1L,"2M1I3M","ACTGTA","IIIIII",1L,5L,20L),"ACTGTA")
  expect_equal(f(1L,"2M1D2M","ACTA","IIII",1L,5L,20L),"ACTA")
  expect_null(f(1L,"2M1N2M","ACTA","IIII",1L,5L,20L))
  expect_equal(f(1L,"2M5N3M","ACGTA","IIIII",8L,10L,20L),"GTA")
  expect_null(f(1L,"5M","ACGTA","II!II",1L,5L,20L))
  expect_null(f(2L,"4M","CGTA","IIII",1L,5L,20L))
})

test_that("native ASCAT saved object and segment tables preserve absolute CN", {
  root <- tempfile();dir.create(root)
  src <- somatic_test_source("ascat");src$sample_id <- "S1";src$reference_build <- "GRCh38"
  src$format <- "ascat_rds";src$options <- list(coordinate_system="1-based-inclusive",normal_cn=2)
  src$resolved_path <- file.path(root,"ascat.rds")
  segments <- data.frame(sample="S1",chr="1",startpos=1,endpos=200,nMajor=2,nMinor=1)
  saveRDS(list(segments=segments,purity=c(S1=.6),ploidy=c(S1=3)),src$resolved_path)
  project <- list(config=list(scientific=list(genome_assembly="GRCh38")))
  out <- Neoantimon2:::input_ascat(src,project)
  cn <- Neoantimon2:::somatic_facets_segments(out)[[1]]
  expect_equal(cn$total_cn,3);expect_equal(cn$major_cn,2);expect_equal(cn$minor_cn,1)
  expect_equal(cn$purity,.6);expect_equal(cn$ploidy,3);expect_equal(cn$normal_cn,2)
  expect_identical(Neoantimon2:::somatic_cn_map(somatic_test_event(),list(cn),out)$status,"mapped")
  src$format <- "ascat_segments";src$resolved_path <- file.path(root,"segments.txt")
  utils::write.table(segments,src$resolved_path,sep="\t",quote=FALSE,row.names=FALSE)
  expect_equal(Neoantimon2:::input_ascat(src,project)$rows[[1]]$tcn.em,"3")
  src$options$coordinate_system <- "0-based";expect_error(Neoantimon2:::input_ascat(src,project),"coordinate")
  src$options$coordinate_system <- "1-based-inclusive";src$reference_build <- "GRCh37"
  expect_error(Neoantimon2:::input_ascat(src,project),"reference_build")
})

test_that("Salmon RSEM and featureCounts keep sample IDs and units", {
  root <- tempfile();dir.create(root);source <- somatic_test_source("expression")
  source$format <- "salmon";source$resolved_path <- file.path(root,"quant.sf")
  writeLines(c("Name\tLength\tEffectiveLength\tTPM\tNumReads","TX1.1\t100\t80\t4.5\t9"),source$resolved_path)
  x <- Neoantimon2:::input_expression(source,NULL)[[1]]
  expect_equal(x$transcript_id,"TX1.1");expect_equal(x$transcript_tpm,"4.5");expect_null(x$gene_tpm)
  source$format <- "featurecounts"
  writeLines(c("# Program:featureCounts","Geneid\tChr\tStart\tEnd\tStrand\tLength\tsample-1","G1\tchr1\t1\t100\t+\t100\t8"),source$resolved_path)
  x <- Neoantimon2:::input_expression(source,NULL)[[1]]
  expect_equal(x$estimated_count,"8");expect_equal(x$expression_unit,"count");expect_null(x$gene_tpm)
  source$sample_id <- "wrong";expect_error(Neoantimon2:::input_expression(source,NULL),"column")
  source$format <- "rsem_genes";source$sample_id <- "sample-1"
  writeLines(c("gene_id\ttranscript_id(s)\tTPM\tFPKM\texpected_count","G1\tTX1\t2\t3\t4"),source$resolved_path)
  x <- Neoantimon2:::input_expression(source,NULL)[[1]]
  expect_equal(x$gene_tpm,"2");expect_equal(x$fpkm,"3")
})
