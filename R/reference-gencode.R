reference_gencode <- function(paths, records, scientific, progress = function(stage) invisible(NULL)) {
  if (!identical(scientific$species, "NCBITaxon:9606")) stop("GENCODE human profile species conflict")
  version <- function(role) records[[role]]$version
  assembly <- scientific$genome_assembly
  if (is.null(assembly) || is.null(version("genome_fasta")) ||
      !identical(assembly, version("genome_fasta")) || !identical(assembly, version("genome_fai")))
    stop("Declared assembly conflict/missing declaration")
  if (!grepl("^GRCh38([.]p[0-9]+)?$", assembly))
    stop("MANE summary requires GRCh38 assembly context")
  release <- version("annotation_gtf"); mane_release <- version("mane_summary")
  if (is.null(release) || is.null(mane_release)) stop("Missing declared annotation/MANE release")
  for (r in c("transcript_fasta", "protein_coding_transcript_fasta", "protein_fasta", "selenocysteine_metadata"))
    if (!identical(version(r), release)) stop("Declared GENCODE release conflict")
  if (!identical(version("mane_gtf"), mane_release)) stop("Declared MANE release conflict")
  progress("Genome FASTA/FAI")
  genome <- reference_fai(paths[["genome_fai"]], paths[["genome_fasta"]])
  progress("Primary GTF")
  gtf <- reference_gtf(paths[["annotation_gtf"]], genome)
  headers <- attr(gtf, "source_headers")
  providers <- sub("^##provider:[[:space:]]*", "", grep("^##provider:", headers, value = TRUE))
  if (length(providers) && any(providers != "GENCODE")) stop("Observed annotation provider conflict")
  descriptions <- grep("^##description:", headers, value = TRUE)
  observed_release <- sub(".*version ([0-9]+).*", "\\1", grep("version [0-9]+", descriptions, value = TRUE))
  observed_assembly <- sub(".*\\((GRCh[0-9]+(?:[.]p[0-9]+)?)\\).*", "\\1",
    grep("\\(GRCh[0-9]+(?:[.]p[0-9]+)?\\)", descriptions, value = TRUE, perl = TRUE), perl = TRUE)
  if (length(observed_release) && any(observed_release != release)) stop("Observed GENCODE release conflict")
  if (length(observed_assembly) && any(observed_assembly != sub("[.]p[0-9]+$", "", assembly) &
      observed_assembly != assembly)) stop("Observed assembly conflict")
  progress("Transcript ALL catalog")
  tx <- reference_fasta(paths[["transcript_fasta"]], "transcript")
  progress("Coding transcript ALL catalog")
  coding <- reference_fasta(paths[["protein_coding_transcript_fasta"]], "transcript")
  progress("Protein ALL catalog")
  proteins <- reference_fasta(paths[["protein_fasta"]], "protein")
  annotated <- unique(gtf[!is.na(gtf$transcript_id), c("transcript_id", "gene_id")])
  if (any(!annotated$transcript_id %in% tx$id)) stop("GTF/transcript FASTA mismatch")
  if (any(!coding$id %in% tx$id)) stop("Coding/transcript FASTA mismatch")
  if (anyDuplicated(proteins$transcript_id) || any(!proteins$transcript_id %in% coding$id)) stop("Protein/transcript mismatch")
  match_gene <- function(a, b, key_a, key_b) {
    idx <- match(a[[key_a]], b[[key_b]]); present <- !is.na(idx)
    if (any(a$gene_id[present] != b$gene_id[idx[present]])) stop("Reference transcript/gene relationship mismatch")
  }
  match_gene(annotated, tx, "transcript_id", "id")
  match_gene(coding, tx, "id", "id")
  match_gene(proteins, tx, "transcript_id", "id")
  if (any(coding$length != tx$length[match(coding$id, tx$id)])) stop("Coding/transcript length mismatch")
  links <- unique(gtf[!is.na(gtf$protein_id), c("transcript_id", "protein_id", "gene_id")])
  if (anyDuplicated(links$transcript_id)) stop("Contradictory annotation protein relationship")
  idx <- match(links$transcript_id, proteins$transcript_id)
  if (anyNA(idx) || any(links$protein_id != proteins$protein_id[idx])) stop("Annotation/protein relationship mismatch")
  progress("MANE metadata")
  mane_gtf <- reference_gtf(paths[["mane_gtf"]])
  mane <- reference_table(paths[["mane_summary"]], 14L, "MANE summary")
  header <- as.character(mane[1, ]); header[1] <- sub("^#", "", header[1])
  expected <- c("NCBI_GeneID", "Ensembl_Gene", "HGNC_ID", "symbol", "name", "RefSeq_nuc",
    "RefSeq_prot", "Ensembl_nuc", "Ensembl_prot", "MANE_status", "GRCh38_chr", "chr_start", "chr_end", "chr_strand")
  if (!identical(header, expected) || nrow(mane) < 2L) stop("Malformed MANE header/record")
  mane <- mane[-1, ]; names(mane) <- header
  if (any(!nzchar(as.matrix(mane[setdiff(names(mane), c("HGNC_ID", "RefSeq_prot", "Ensembl_prot"))]))) || anyDuplicated(mane$Ensembl_nuc) ||
      any(!mane$MANE_status %in% c("MANE Select", "MANE Plus Clinical")) ||
      any(!mane$chr_strand %in% c("+", "-"))) stop("Malformed MANE record")
  mane$chr_start <- reference_integer(mane$chr_start, "MANE start", 1)
  mane$chr_end <- reference_integer(mane$chr_end, "MANE end", 1)
  if (any(mane$chr_start > mane$chr_end)) stop("Invalid MANE interval")
  mt <- unique(mane_gtf[!is.na(mane_gtf$transcript_id), c("transcript_id", "gene_id")])
  if (any(!mane$Ensembl_nuc %in% mt$transcript_id) || any(!mt$transcript_id %in% mane$Ensembl_nuc))
    stop("MANE identifier mismatch between summary and GTF")
  if (any(mane$Ensembl_Gene != mt$gene_id[match(mane$Ensembl_nuc, mt$transcript_id)]))
    stop("MANE gene identifier mismatch")
  mane <- reference_mane_scope(mane, mane_gtf, genome, gtf, tx, proteins)
  classify <- function(ids) ifelse(ids %in% annotated$transcript_id,
    "PRIMARY_RECONSTRUCTION_ELIGIBLE", "CATALOG_ONLY_OUT_OF_SCOPE")
  tx$reconstruction_scope <- classify(tx$id)
  coding$reconstruction_scope <- classify(coding$id)
  proteins$reconstruction_scope <- classify(proteins$transcript_id)
  proteins$self_catalog_member <- TRUE
  progress("Selenocysteine metadata")
  sec <- reference_table(paths[["selenocysteine_metadata"]], 2L, "Sec metadata")
  names(sec) <- c("transcript_id", "position")
  sec$position <- reference_integer(sec$position, "Sec position", 1)
  if (any(!sec$transcript_id %in% tx$id)) stop("Sec metadata references unknown entity")
  sec_idx <- match(sec$transcript_id, proteins$transcript_id)
  if (anyNA(sec_idx) || any(sec$position > proteins$length[sec_idx])) stop("Sec/protein relationship mismatch")
  sec$reconstruction_scope <- classify(sec$transcript_id)
  sec$catalog_traceability <- "GENCODE_ALL_TRANSCRIPT_AND_PROTEIN"
  list(profile = list(assembly_release = assembly, coordinate_scope = "primary_assembly",
    annotation_source = "GENCODE", annotation_release = release, annotation_scope = "primary_assembly",
    sequence_catalog_scope = "all_regions", mane_release = mane_release,
    scope_policy = reference_scope_policy()),
    genome = genome, annotation = gtf, transcripts = tx, coding_transcripts = coding,
    proteins = proteins, mane = reference_sort(mane), mane_annotation = mane_gtf, sec = reference_sort(sec),
    self_proteome_resource = "reference.protein_fasta", self_proteome_policy = "all supplied protein isoforms",
    coordinate_convention = "1-based inclusive", mitochondrial_contigs = intersect("chrM", genome$contig),
    observed = list(annotation_headers = headers, annotation_release = observed_release,
      assembly = observed_assembly, mane_headers = attr(mane_gtf, "source_headers")),
    mismatches = list(gtf_without_transcript = setdiff(annotated$transcript_id, tx$id),
      transcript_without_gtf = setdiff(tx$id, annotated$transcript_id),
      coding_without_transcript = setdiff(coding$id, tx$id),
      protein_without_gtf = setdiff(proteins$transcript_id, annotated$transcript_id),
      coding_without_protein = setdiff(coding$id, proteins$transcript_id),
      mane_without_exact_gencode = setdiff(mane$Ensembl_nuc, annotated$transcript_id),
      sec_without_gtf = setdiff(sec$transcript_id, annotated$transcript_id)),
    capabilities = list(genome_sequence_available = TRUE, transcript_sequence_available = TRUE,
      protein_reference_available = TRUE, mane_available = TRUE, sec_metadata_available = TRUE,
      assembly_verified = if (assembly %in% observed_assembly) TRUE else NA))
}
