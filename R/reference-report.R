reference_summary <- function(bundle) {
  b <- bundle
  if (!is.null(b$validation_summary)) return(b$validation_summary)
  scope_counts <- function(x) {
    groups <- split(seq_len(nrow(x)), x$scope)
    lapply(groups, function(i) list(records = length(i), transcripts = length(unique(x$Ensembl_nuc[i])),
      genes = length(unique(x$Ensembl_Gene[i])), contigs = length(unique(x$seqname[i]))))
  }
  m <- b$mane
  region <- attr(b$mane_annotation, "region_feature_counts")
  outside <- m$scope == "OUT_OF_SCOPE_REFERENCE_REGION"
  reasons <- lapply(c("MANE_PATCH_SEQUENCE", "MANE_ALTERNATE_LOCUS", "OTHER_NON_PRIMARY"), function(reason) {
    rows <- which(m$scope_reason == reason); contigs <- unique(m$seqname[rows])
    list(reason = reason, records = length(rows), genes = length(unique(m$Ensembl_Gene[rows])),
      contigs = length(contigs), feature_rows = sum(region[contigs]))
  })
  list(status = "PASS", family = b$family, analysis_id = b$analysis_id,
    profile = b$profile, declared = b$declared, observed = b$observed,
    resource_paths = lapply(b$resources, function(r) r[c("logical_name", "declared_path", "resolved_path", "version", "source")]),
    checksum_state = as.list(b$checksum_state),
    checksum_records = lapply(b$resources, `[[`, "integrity"),
    counts = list(genome_contigs = nrow(b$genome), primary_gtf_seqnames = length(attr(b$annotation, "seqnames")),
      primary_gtf_transcripts = nrow(b$annotation), transcripts = nrow(b$transcripts),
      coding_transcripts = nrow(b$coding_transcripts), proteins = nrow(b$proteins),
      reconstruction_eligible_transcripts = sum(b$transcripts$reconstruction_scope == "PRIMARY_RECONSTRUCTION_ELIGIBLE"),
      catalog_only_transcripts = sum(b$transcripts$reconstruction_scope == "CATALOG_ONLY_OUT_OF_SCOPE"),
      sec_records = nrow(b$sec), sec_catalog_only = sum(b$sec$reconstruction_scope == "CATALOG_ONLY_OUT_OF_SCOPE"),
      mane_records = nrow(m),
      mane_select = sum(m$MANE_status == "MANE Select"),
      mane_plus_clinical = sum(m$MANE_status == "MANE Plus Clinical"),
      mane_in_scope_select = sum(!outside & m$MANE_status == "MANE Select"),
      mane_in_scope_plus_clinical = sum(!outside & m$MANE_status == "MANE Plus Clinical"),
      mane_out_of_scope_select = sum(outside & m$MANE_status == "MANE Select"),
      mane_out_of_scope_plus_clinical = sum(outside & m$MANE_status == "MANE Plus Clinical"),
      mane_seqnames = length(attr(b$mane_annotation, "seqnames")),
      mane_feature_rows = sum(region), mane_out_of_scope_feature_rows = sum(region[unique(m$seqname[outside])])),
    mane_scope = scope_counts(m), mane_reasons = reasons, unexplained_mane_mismatches = 0L,
    mismatches = lapply(b$mismatches, length), self_catalog_scope = b$self_proteome_policy,
    mitochondrial_contigs = as.list(b$mitochondrial_contigs),
    limitations = c("No sequence reconstruction or translation", "No self matching",
      "Assembly patch release is declared unless explicitly observed", "Source integrity state is reported per resource"))
}

# Storage projection only: validation has already completed. Scientific scope,
# all bound resources and the complete MANE/Sec metadata remain unchanged.
reference_bundle_metadata <- function(bundle) {
  summary <- reference_summary(bundle)
  b <- bundle[setdiff(names(bundle), c("annotation", "transcripts", "coding_transcripts",
    "proteins", "mane_annotation"))]
  b$validation_summary <- summary
  b$indexes_materialized <- FALSE
  b$catalog_locators <- list(annotation = "reference.annotation_gtf",
    transcripts = "reference.transcript_fasta", coding_transcripts = "reference.protein_coding_transcript_fasta",
    proteins = "reference.protein_fasta", mane_annotation = "reference.mane_gtf")
  class(b) <- class(bundle)
  b
}
