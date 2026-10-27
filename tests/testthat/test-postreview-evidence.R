test_that("depth-only assertions do not conceal VCF allele evidence", {
  e <- somatic_test_event();e$event_role <- "somatic"
  e$event_details_json <- list(source_records=list(list(format="VCF",input_id="truth",
    line=5L,original=list(alts=list("A")),selected_alt_index=1L,
    raw_format=list(keys="GT:AD:DP",sample="0/1:80,20:100"))))
  source <- somatic_test_source("dna_reads")
  f <- Neoantimon2:::input_assay
  x <- f(e,list(tumor_depth="100"),"tumor",source)
  expect_equal(x$vaf,.2);expect_equal(x$ref_reads,80);expect_equal(x$alt_reads,20)
  expect_equal(x$depth,100);expect_equal(x$vaf_denominator,100)
  expect_identical(x$measurement_sources$alleles$method,"FORMAT/AD")
  expect_equal(x$alternatives$explicit$depth,100);expect_false(x$source_disagreement)
  x <- f(e,list(tumor_depth="90"),"tumor",source)
  expect_equal(x$vaf,.2);expect_equal(x$depth,100);expect_true(x$source_disagreement)
  expect_equal(x$alternatives$explicit$depth,90)
  x <- f(e,list(tumor_alt_reads="40"),"tumor",source)
  expect_null(x$vaf);expect_identical(x$status,"ambiguous")
  expect_null(x$alt_reads);expect_equal(x$ambiguous_allele_observation$alt_reads,20)
  expect_equal(x$alternatives$explicit$alt_reads,40)
  expect_equal(f(e,list(tumor_alt_reads="20"),"tumor",source)$vaf,.2)
  e$event_details_json$source_records[[1]]$raw_format <- list(keys="GT:AD",sample="0/1:80,20")
  x <- f(e,list(tumor_depth="100"),"tumor",source)
  expect_equal(x$vaf,.2);expect_equal(x$depth,100)
  expect_identical(x$measurement_sources$depth$source_id,source$source_id)
  expect_identical(x$measurement_sources$alleles$method,"FORMAT/AD")
  x <- f(e,list(tumor_depth="0"),"tumor",source)
  expect_equal(x$vaf,.2);expect_true(x$depth_count_disagreement);expect_true(x$source_disagreement)
})

postreview_bam <- function(rows,extra=character(),length=33L) {
  root <- tempfile("postreview-bam-");dir.create(root)
  sam <- file.path(root,"reads.sam")
  writeLines(c("@HD\tVN:1.6\tSO:coordinate",paste0("@SQ\tSN:chr1\tLN:",length),
    extra,"@RG\tID:rg1\tSM:sample-1",rows),sam)
  bam <- Rsamtools::asBam(sam,destination=file.path(root,"reads"),indexDestination=TRUE)
  source <- somatic_test_source("rna_bam");source$reference_build <- "GRCh38.p14"
  source$resolved_path <- bam;source$index_path <- paste0(bam,".bai")
  source$sha256 <- digest::digest(file=bam,algo="sha256",serialize=FALSE)
  source$index_sha256 <- digest::digest(file=source$index_path,algo="sha256",serialize=FALSE)
  source
}
postreview_sam <- function(name,sequence,cigar="33M",flag=0L) paste(name,flag,"chr1",1,60,
  cigar,"*",0,0,sequence,paste(rep("I",nchar(sequence)),collapse=""),"RG:Z:rg1",sep="\t")
postreview_event <- function(type="SNV",pos=4L,ref="G",alt="A") {
  e <- somatic_test_event();e$event_type <- type;e$pos <- pos
  e$end <- pos+nchar(ref)-1L;e$ref <- ref;e$alt <- alt;e
}

test_that("only queried BAM contigs must match the selected reference", {
  skip_if_not_installed("Rsamtools")
  p <- reconstruction_test_project();e <- postreview_event()
  ref <- "ATGGAACCAAAAAAAAAAAATAAACCCGGGTAA"
  rows <- c(postreview_sam("ref",ref),postreview_sam("alt",sub("ATGG","ATGA",ref)))
  src <- postreview_bam(rows,extra=c("@SQ\tSN:chrUn_decoy\tLN:1000","@SQ\tSN:chr1_ALT\tLN:80"))
  count <- Neoantimon2:::input_bam_observations(src,list(e),p)[[e$event_id]]$rna
  expect_equal(count$vaf,.5);expect_equal(count$ref_reads,1);expect_equal(count$alt_reads,1)
  src <- postreview_bam(rows,extra="@SQ\tSN:chrUn_decoy\tLN:1000",length=34L)
  expect_error(Neoantimon2:::input_bam_observations(src,list(e),p),"Queried BAM contig length")
  src <- postreview_bam(rows);e$contig <- "chrUnknown"
  expect_error(Neoantimon2:::input_bam_observations(src,list(e),p),"Contig")
  e <- postreview_event();e$ref <- "T"
  expect_error(Neoantimon2:::input_bam_observations(src,list(e),p),"REF mismatch")
})

test_that("indel target edits tolerate nearby SNVs and preserve repeat and fragment truth", {
  skip_if_not_installed("Rsamtools")
  p <- reconstruction_test_project();ref <- "ATGGAACCAAAAAAAAAAAATAAACCCGGGTAA"
  get <- function(rows,e) Neoantimon2:::input_bam_observations(postreview_bam(rows),list(e),p)[[e$event_id]]$rna
  for(type in c("INS","DEL")) {
    e <- if(type=="INS") postreview_event(type,4L,"G","GT") else postreview_event(type,4L,"GA","G")
    alt <- if(type=="INS") paste0(substr(ref,1,4),"T",substring(ref,5)) else paste0(substr(ref,1,4),substring(ref,6))
    # A separate SNV at genomic position 2 is independent of the position-4 indel.
    substring(alt,2L,2L) <- "C"
    # Also mutate genomic position 6 immediately next to the indel. Query
    # coordinates differ by +1 for INS and -1 for DEL.
    at <- if(type=="INS") 7L else 5L;substring(alt,at,at) <- "C"
    cigar <- if(type=="INS") "4M1I29M" else "4M1D28M"
    rows <- c(postreview_sam("ref",ref),postreview_sam("alt-pair",alt,cigar,65L),postreview_sam("alt-pair",alt,cigar,129L))
    x <- get(rows,e)
    expect_equal(x$ref_reads,1);expect_equal(x$alt_reads,1);expect_equal(x$vaf,.5)
    expect_equal(x$usable_depth,2)
    x <- get(postreview_sam("only-ref",ref),e)
    expect_equal(x$vaf,0);expect_equal(x$alt_reads,0)
    # The pair makes contradictory target-allele assertions: withhold VAF.
    x <- get(c(postreview_sam("ref",ref),postreview_sam("pair",alt,cigar,65L),postreview_sam("pair",ref,"33M",129L)),e)
    expect_null(x$vaf);expect_equal(x$ambiguous_reads,1);expect_identical(x$status,"ambiguous")
    expect_null(x$alt_reads);expect_null(x$ref_reads);expect_equal(x$resolved_ref_fragments,1)
    # A splice beyond the required context does not negate the allele.
    splice_cigar <- if(type=="INS") "4M1I16M3N10M" else "4M1D15M3N10M"
    splice_seq <- if(type=="INS") paste0(substr(ref,1,4),"T",substr(ref,5,20),substring(ref,24)) else
      paste0(substr(ref,1,4),substr(ref,6,20),substring(ref,24))
    expect_equal(get(postreview_sam("splice",splice_seq,splice_cigar),e)$vaf,1)
  }
  # Independent truth: A insertion/deletion in the A9..A20 repeat can be placed
  # later by CIGAR; both normalize to the position-8 C anchor.
  for(type in c("INS","DEL")) {
    e <- if(type=="INS") postreview_event(type,8L,"C","CA") else postreview_event(type,8L,"CA","C")
    alt <- if(type=="INS") paste0(substr(ref,1,12),"A",substring(ref,13)) else
      paste0(substr(ref,1,12),substring(ref,14))
    cigar <- if(type=="INS") "12M1I21M" else "12M1D20M"
    substring(alt,2L,2L) <- "C"
    x <- get(c(postreview_sam("ref",ref),postreview_sam("repeat-alt",alt,cigar)),e)
    expect_equal(x$alt_reads,1);expect_equal(x$ref_reads,1);expect_equal(x$vaf,.5)
  }
  # RNA splice interrupts required indel context: missing, never observed zero.
  e <- postreview_event("INS",4L,"G","GT")
  x <- get(postreview_sam("gap","ATGGAAA","4M3N3M"),e)
  expect_null(x$vaf);expect_equal(x$usable_depth,0);expect_identical(x$status,"insufficient_depth")
  # A different insertion is not fabricated as target REF or ALT.
  alt <- paste0(substr(ref,1,4),"C",substring(ref,5))
  x <- get(c(postreview_sam("ref",ref),postreview_sam("different-indel",alt,"4M1I29M")),e)
  expect_equal(x$vaf,0);expect_equal(x$other_reads,1);expect_equal(x$alt_reads,0)
  # Multiple interacting edits cannot safely become negative ALT evidence.
  complex <- paste0(substr(ref,1,4),"T",substr(ref,5,5),"C",substring(ref,6))
  x <- get(c(postreview_sam("ref",ref),postreview_sam("complex",complex,"4M1I1M1I28M")),e)
  expect_null(x$vaf);expect_null(x$alt_reads);expect_equal(x$ambiguous_reads,1)
  # A deletion beginning before the queried window spans the target: it cannot
  # be misclassified as REF merely because its CIGAR anchor is outside it.
  e <- postreview_event("INS",24L,"A","AT")
  spanning <- paste0(substr(ref,1,3),substring(ref,29))
  x <- get(c(postreview_sam("ref",ref),postreview_sam("spanning",spanning,"3M25D5M")),e)
  expect_null(x$vaf);expect_null(x$ref_reads);expect_equal(x$ambiguous_reads,1)
})

test_that("expression disagreements remain unavailable with every assertion retained", {
  f <- Neoantimon2:::priority_expression_value
  rows <- list(list(gene_tpm=3,provenance=list(source="A")),list(gene_tpm=9,provenance=list(source="B")))
  expect_identical(f(rows,"gene_tpm")$status,"conflict");expect_null(f(rows,"gene_tpm")$value)
  expect_identical(f(rev(rows),"gene_tpm"),f(rows,"gene_tpm"))
  rows[[2]]$gene_tpm <- 3
  expect_equal(f(rows,"gene_tpm")$value,3)
  expect_identical(f(rows,"transcript_tpm")$status,"unavailable")
  rows[[2]]$gene_tpm <- Inf;expect_error(f(rows,"gene_tpm"),"Invalid expression")
})

test_that("ASCAT upstream sample selection is explicit and metadata stays bound", {
  src <- somatic_test_source("ascat");src$sample_id <- "S1";src$reference_build <- "GRCh38"
  src$format <- "ascat_rds";src$options <- list(coordinate_system="1-based-inclusive")
  src$resolved_path <- tempfile(fileext=".rds")
  obj <- list(segments=data.frame(sample=c("tumor_A","tumor_B"),chr="1",startpos=1,endpos=200,nMajor=c(2,4),nMinor=1),
    purity=c(tumor_A=.6,tumor_B=.4),ploidy=c(tumor_A=3,tumor_B=5))
  saveRDS(obj,src$resolved_path);p <- list(config=list(scientific=list(genome_assembly="GRCh38")))
  expect_error(Neoantimon2:::input_ascat(src,p),"sample identities")
  src$options$sample_name <- "tumor_A"
  mapped <- Neoantimon2:::input_ascat(src,p)
  expect_length(mapped$rows,1);expect_identical(mapped$rows[[1]]$sample_id,"S1")
  expect_equal(mapped$rows[[1]]$purity,"0.6");expect_equal(mapped$rows[[1]]$tcn.em,"3")
  expect_identical(Neoantimon2:::somatic_source_provenance(mapped)$sample_mapping$upstream_sample,"tumor_A")
  src$options$sample_name <- "unknown"
  expect_error(Neoantimon2:::input_ascat(src,p),"sample identities")
  src$options$sample_name <- "tumor_A";obj$purity <- .5;saveRDS(obj,src$resolved_path)
  expect_error(Neoantimon2:::input_ascat(src,p),"must be named")
})

test_that("skipped DNA fallback retains full cryptographic source verification", {
  skip_if_not_installed("Rsamtools")
  ref <- "ATGGAACCAAAAAAAAAAAATAAACCCGGGTAA"
  src <- postreview_bam(postreview_sam("ref",ref));p <- reconstruction_test_project()
  p$config$somatic_evidence <- Neoantimon2:::somatic_evidence_policy(list(sources=list(
    list(kind="dna_bam",sample_id="sample-1",path=src$resolved_path,reference_build="GRCh38.p14"))))
  bound <- Neoantimon2:::somatic_source_files(p)[[1L]]
  expect_identical(bound$sha256,digest::digest(file=src$resolved_path,algo="sha256",serialize=FALSE))
  expect_identical(bound$index_sha256,digest::digest(file=src$index_path,algo="sha256",serialize=FALSE))
  e <- postreview_event();e$event_details_json <- list(source_records=list(list(format="VCF",
    input_id="vcf",line=5L,original=list(alts=list("A")),selected_alt_index=1L,
    raw_format=list(keys="AD",sample="80,20"))))
  expect_length(Neoantimon2:::input_bam_observations(bound,list(e),p),0L)
  p$config$somatic_evidence$sources[[1L]]$checksum <- bound$sha256
  # Only the independently authored temporary test BAM is modified.
  con <- file(src$resolved_path,"ab");writeBin(as.raw(1),con);close(con)
  expect_error(Neoantimon2:::somatic_source_files(p),"checksum mismatch")
})
