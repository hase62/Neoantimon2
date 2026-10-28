# SYNTHETIC measurements and predictions, not biological or model validation.
# Run: Rscript --vanilla source-evidence.R /new/output/directory
library(Neoantimon2)
run_source_evidence_demo <- function(output) {
  if (file.exists(output)) stop("Use a new output directory")
  dir.create(output,recursive=TRUE)
  source_dir <- file.path(output,"inputs");dir.create(source_dir)
  refs <- system.file("extdata","quickstart",package="Neoantimon2")
  files <- c(genome_fasta="genome.fa",genome_fai="genome.fa.fai",annotation_gtf="annotation.gtf",
    transcript_fasta="transcripts.fa",protein_coding_transcript_fasta="coding.fa",
    protein_fasta="proteins.fa",selenocysteine_metadata="sec.tsv",mane_summary="mane.tsv",mane_gtf="mane.gtf")
  resources <- lapply(names(files),function(role) list(path=file.path(refs,files[[role]]),
    required=TRUE,resource_id=if(role=="annotation_gtf") "gencode_primary_assembly_all_catalogs_mane_supplemental_v1" else role,
    version=if(startsWith(role,"genome")) "GRCh38.p14" else if(startsWith(role,"mane")) "1.5" else "50",source="synthetic"))
  names(resources) <- paste0("reference.",names(files))
  write_vcf <- function(name,ad=TRUE) {
    file <- file.path(source_dir,paste0(name,".vcf"))
    writeLines(c("##fileformat=VCFv4.2",
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>",
      if(ad) "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allelic depths>" else character(),
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample-1",
      if(ad) "chr1\t4\t.\tG\tA\t.\tPASS\t.\tGT:AD\t0/1:80,20" else
        "chr1\t4\t.\tG\tA\t.\tPASS\t.\tGT\t0/1"),file)
    file
  }
  variants <- write_vcf("tumor");missing <- write_vcf("missing",FALSE)
  project <- function(name,vcf=variants,sources=list()) {
    cfg <- neo_config(env=character(),resources=resources,checksum="sha256",
      paths=list(input_root=source_dir,data_root=source_dir,run_root=file.path(output,"runs")),
      scientific=list(genome_assembly="GRCh38.p14"),
      peptides=list(class_i_lengths=3L,class_ii_lengths=5L),
      somatic_evidence=if(length(sources)) list(sources=sources) else NULL)
    neo_run(neo_project(samples=data.frame(sample_id="sample-1",patient_id="P1"),
      hla=list(P1=c("HLA-A*02:01","HLA-B*07:02")),config=cfg,project_id=name,
      inputs=list(list(input_id=name,input_role="somatic_vcf",sample_id="sample-1",path=vcf))))
  }
  # A: exactly one run, no BAM or internal IDs.
  a <- project("A-minimal")
  za <- neo_priorities(a)
  stopifnot(all(za$dna_vaf==.2),all(is.na(za$raw_score)))
  b <- NULL
  if (requireNamespace("Rsamtools",quietly=TRUE)) {
    sam <- file.path(output,"rna.sam")
    record <- function(id,flag,seq,cigar="8M",mapq=60,qual="IIIIIIII")
      paste(id,flag,"chr1",1,mapq,cigar,"*",0,0,seq,qual,"RG:Z:rg1",sep="\t")
    writeLines(c("@HD\tVN:1.6\tSO:coordinate","@SQ\tSN:chr1\tLN:33","@RG\tID:rg1\tSM:sample-1",
      record("ref",0,"ATGGAACC"),record("alt",0,"ATGAAACC"),
      record("pair",65,"ATGAAACC"),record("pair",129,"ATGAAACC"),
      record("splice",0,"ATGAAAAA","5M3N3M"),record("dup",1024,"ATGAAACC"),
      record("secondary",256,"ATGAAACC"),record("supp",2048,"ATGAAACC"),
      record("lowmapq",0,"ATGAAACC",mapq=0),record("lowbase",0,"ATGAAACC",qual="III!IIII")),sam)
    bam <- Rsamtools::asBam(sam,destination=file.path(output,"rna"),overwrite=FALSE,indexDestination=TRUE)
    b <- project("B-rna",sources=list(list(kind="rna_bam",sample_id="sample-1",path=bam,
      reference_build="GRCh38.p14",source="independently authored synthetic alignments")))
    zb <- neo_priorities(b)
    stopifnot(all(zb$rna_vaf==.75),all(zb$dna_vaf==.2))
  }
  # C: unchanged external expression/CN outputs, no evidence-table preparation.
  quant <- file.path(output,"quant.sf")
  writeLines(c("Name\tLength\tEffectiveLength\tTPM\tNumReads",
    "TA.1\t21\t15\t4.5\t9","TB.1\t21\t15\t3\t6"),quant)
  genes <- file.path(output,"genes.results")
  writeLines(c("gene_id\ttranscript_id(s)\tTPM\tFPKM\texpected_count",
    "G1.1\tTA.1,TB.1\t4.5\t3\t9"),genes)
  supplied_ccf <- file.path(output,"ccf.tsv")
  writeLines(c("contig\tpos\tref\talt\tccf\tmethod",
    "chr1\t4\tG\tA\t1\tsynthetic externally supplied CCF"),supplied_ccf)
  asc <- file.path(output,"ascat.rds")
  saveRDS(list(segments=data.frame(sample="sample-1",chr="1",startpos=1,endpos=33,nMajor=1,nMinor=1),
    purity=c("sample-1"=.4),ploidy=c("sample-1"=2)),asc)
  cc <- project("C-expression-CN",sources=list(
    list(kind="expression",format="rsem_genes",sample_id="sample-1",path=genes,source="synthetic RSEM genes schema"),
    list(kind="ccf",sample_id="sample-1",path=supplied_ccf,reference_build="GRCh38.p14",source="synthetic supplied estimate"),
    list(kind="ascat",format="ascat_rds",sample_id="sample-1",path=asc,reference_build="GRCh38.p14",
      source="synthetic ASCAT schema",options=list(coordinate_system="1-based-inclusive",normal_cn=2))))
  zc <- neo_priorities(cc)
  stopifnot(any(!is.na(zc$gene_tpm)),all(vapply(zc$biological_evidence,function(e)e[[1]]$cn$total_cn==2,TRUE)))
  # D/E: explicitly synthetic external predictions for this exact candidate set.
  predictions <- function(p,task,metric,provider) {
    candidates <- neo_results(p)
    data.frame(schema_version="1.0.0",peptide_sequence=candidates$peptide_sequence,
      scope="peptide_hla",mhc_key=candidates$mhc_key,task=task,provider=provider,
      predictor_version="synthetic-1",model_version="synthetic-1",run_id="synthetic-run",
      metric=metric,metric_kind="raw_score",metric_semantics="raw",
      score=seq_len(nrow(candidates))/10,score_direction="higher_is_better",unit="unitless",
      status="evaluated",reason="",source_row=paste0("synthetic:",seq_len(nrow(candidates))))
  }
  import <- function(p,category,task) {
    file <- file.path(output,paste0(category,"-predictions.tsv"));provider <- paste0("synthetic-",category)
    data <- predictions(p,task,if(category=="BA") "affinity_nM" else "score",provider)
    if(category=="BA") {data$unit <- "nM";data$score_direction <- "lower_is_better";data$score <- seq_len(nrow(data))*100}
    utils::write.table(data,file,sep="\t",quote=FALSE,row.names=FALSE)
    neo_analyze_external(p,file,policy=list(mode="external_rank",provider=provider,run_id="synthetic-run",
      task=task,metric=if(category=="BA") "affinity_nM" else "score",
      score_direction=if(category=="BA") "lower_is_better" else "higher_is_better"))
  }
  native <- neo_results(a)
  el <- import(a,"EL","mhc1_presentation");im <- import(a,"IM","external:immunogenicity")
  zel <- neo_priorities(a,analysis=el);zim <- neo_priorities(a,analysis=im)
  stopifnot(identical(neo_results(a),native),all(zel$category=="EL"),all(zim$category=="IM"),
    all(!is.na(zim$prediction_only)),all(zim$vaf_used==.2))
  e <- project("E-missing",missing)
  ei <- import(e,"IM-missing","external:immunogenicity");ze <- neo_priorities(e,analysis=ei)
  stopifnot(all(is.na(ze$dna_vaf)),all(ze$vaf_missing_substituted),all(ze$evidence_adjusted==0),
    all(!is.na(ze$prediction_only)),identical(neo_results(a),native),
    identical(neo_results(neo_project(path=a$path)),native))
  ba <- import(cc,"BA","mhc1_binding")
  zba <- neo_priorities(cc,analysis=ba)
  candidates <- neo_results(cc);origins <- neo_origins(cc)
  wt <- predictions(cc,"mhc1_binding","affinity_nM","synthetic-BA")[FALSE,]
  for(i in seq_len(nrow(candidates))) for(id in unlist(candidates$origin_ids[[i]])) {
    o <- origins[origins$origin_id==id,,drop=FALSE]
    if(nrow(o)==1L && !is.na(o$wt_peptide_sequence) && nzchar(o$wt_peptide_sequence)) {
      row <- predictions(cc,"mhc1_binding","affinity_nM","synthetic-BA")[i,,drop=FALSE]
      row$peptide_sequence <- o$wt_peptide_sequence;row$score <- 800;row$unit <- "nM"
      row$score_direction <- "lower_is_better";row$run_id <- "synthetic-run";row$source_row <- "synthetic-WT"
      wt <- rbind(wt,row)
    }
  }
  wt <- unique(wt);wt_file <- file.path(output,"wild-type.tsv")
  utils::write.table(wt,wt_file,sep="\t",quote=FALSE,row.names=FALSE)
  paper <- neo_priorities(cc,analysis=ba,method="neoantimon-paper-2020-v1",
    wild_type_predictions=wt_file,expression_unit="TPM")
  stopifnot(any(is.finite(paper$paper_with_ccf)),identical(neo_results(a),native))
  # Keep the full provenance in RDS; flat CSV intentionally excludes list columns.
  views <- list(A=za,B=if(is.null(b)) NULL else neo_priorities(b),C=zc,D_BA=zba,D_EL=zel,D_IM=zim,D_paper=paper,E=ze)
  saveRDS(views,file.path(output,"numerical-views.rds"))
  for (name in names(views)) if(!is.null(views[[name]])) {
    flat <- views[[name]][,!vapply(views[[name]],is.list,TRUE),drop=FALSE]
    utils::write.csv(flat,file.path(output,paste0(name,".csv")),row.names=FALSE)
  }
  result <- list(A=a,B=b,C=cc,E=e,external_BA=ba,external_EL=el,external_IM=im,views=views)
  saveRDS(result,file.path(output,"demo.rds"))
  cat("A VCF VAF=.2; B RNA VAF=",if(is.null(b)) "SKIP (Rsamtools unavailable)" else ".75",
    "; C expression/CN imported; D EL/IM independent; E missing VAF adjusted=0.\n",sep="")
  result
}
if (sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  if(length(args)!=1L) stop("Usage: Rscript source-evidence.R /new/output")
  run_source_evidence_demo(args[[1L]])
}
