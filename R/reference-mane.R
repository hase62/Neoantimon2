reference_mane_scope <- function(mane, annotation, genome, primary, transcripts, proteins) {
  idx <- match(mane$Ensembl_nuc, annotation$transcript_id)
  if (anyNA(idx)) stop("MANE transcript declaration missing")
  a <- annotation[idx, , drop = FALSE]
  has_protein <- nzchar(mane$Ensembl_prot)
  if (any(has_protein != nzchar(mane$RefSeq_prot)) ||
      any(!has_protein & (!is.na(a$protein_id) | is.na(a$transcript_type) | a$transcript_type == "protein_coding")))
    stop("MANE coding/noncoding annotation conflict")
  protein_conflict <- has_protein & (is.na(a$protein_id) | mane$Ensembl_prot != a$protein_id)
  if (any(mane$Ensembl_Gene != a$gene_id | mane$chr_start != a$start |
          mane$chr_end != a$end | mane$chr_strand != a$strand |
          protein_conflict))
    stop("MANE summary/GTF relationship mismatch")
  # Evidence comes from the bound MANE sources jointly: explicit typed GTF
  # sequence-region name, matching versioned transcript/gene/protein and exact
  # coordinates/strand, plus a RefSeq nonchromosomal genomic accession in summary.
  nonprimary_accession <- grepl("^(NW|NT)_[0-9]+[.][0-9]+$", mane$GRCh38_chr)
  patch <- grepl("^chr([1-9][0-9]*|X|Y)_[A-Z]{2}[0-9]+v[0-9]+_fix$", a$seqname)
  alternate <- grepl("^chr([1-9][0-9]*|X|Y)_[A-Z]{2}[0-9]+v[0-9]+_alt$", a$seqname)
  inside <- a$seqname %in% genome$contig
  unexplained <- !inside & !(nonprimary_accession & (patch | alternate))
  if (any(unexplained)) stop("UNEXPLAINED_REFERENCE_MISMATCH: ",
    paste(unique(a$seqname[unexplained]), collapse = ", "))
  # No unknown gene-only region can evade classification via transcript omission.
  if (length(setdiff(attr(annotation, "seqnames"), a$seqname)))
    stop("UNEXPLAINED_REFERENCE_MISMATCH: MANE region without summary transcript")
  if (any(inside & a$end > genome$length[match(a$seqname, genome$contig)], na.rm = TRUE))
    stop("MANE coordinate exceeds primary genome")
  # A region maps to exactly one supplied genomic accession, and vice versa.
  mapping <- unique(data.frame(contig = a$seqname, accession = mane$GRCh38_chr))
  if (anyDuplicated(mapping$contig) || anyDuplicated(mapping$accession)) stop("MANE sequence-region identity conflict")
  t <- match(mane$Ensembl_nuc, transcripts$id)
  if (anyNA(t) || any(transcripts$gene_id[t] != mane$Ensembl_Gene))
    stop("MANE transcript sequence relationship mismatch")
  p <- match(mane$Ensembl_prot[has_protein], proteins$id)
  if (anyNA(p) || any(proteins$transcript_id[p] != mane$Ensembl_nuc[has_protein] |
      proteins$gene_id[p] != mane$Ensembl_Gene[has_protein])) stop("MANE protein sequence relationship mismatch")
  pi <- match(mane$Ensembl_nuc, primary$transcript_id)
  if (any(inside & is.na(pi))) stop("In-scope MANE transcript absent from primary annotation")
  in_idx <- which(inside)
  if (any(primary$seqname[pi[in_idx]] != a$seqname[in_idx] |
      primary$gene_id[pi[in_idx]] != mane$Ensembl_Gene[in_idx] |
      primary$strand[pi[in_idx]] != a$strand[in_idx])) stop("MANE/primary context mismatch")
  mane$protein_reference_available <- has_protein
  mane$seqname <- a$seqname
  mane$scope <- ifelse(inside, "IN_SCOPE_PRIMARY", "OUT_OF_SCOPE_REFERENCE_REGION")
  mane$scope_reason <- ifelse(inside, "PRIMARY_SEQUENCE", ifelse(patch, "MANE_PATCH_SEQUENCE", "MANE_ALTERNATE_LOCUS"))
  mane$scope_evidence <- ifelse(inside, "EXACT_GENOME_AND_PRIMARY_GTF_SEQUENCE_RELATIONSHIP",
    "MANE_TYPED_REGION_AND_SUMMARY_ACCESSION_IDENTIFIERS_COORDINATES")
  mane$reconstruction_eligible <- inside
  mane$representative_eligible <- inside & mane$MANE_status == "MANE Select"
  reference_sort(mane)
}
