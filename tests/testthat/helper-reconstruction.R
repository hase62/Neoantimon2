reconstruction_test_project <- function(inputs = list()) {
  root <- normalizePath(testthat::test_path("..", "fixtures", "reconstruction"))
  roles <- c(genome_fasta = "genome.fa", genome_fai = "genome.fa.fai",
    annotation_gtf = "annotation.gtf", transcript_fasta = "transcripts.fa",
    protein_coding_transcript_fasta = "coding.fa", protein_fasta = "proteins.fa",
    selenocysteine_metadata = "sec.tsv", mane_summary = "mane.tsv", mane_gtf = "mane.gtf")
  resources <- lapply(names(roles), function(role) list(path = file.path(root, roles[[role]]),
    required = TRUE, resource_id = if (role == "annotation_gtf")
      Neoantimon2:::reference_scope_policy() else NULL,
    version = if (startsWith(role, "genome")) "GRCh38.p14" else
      if (startsWith(role, "mane")) "1.5" else "50", source = "synthetic"))
  names(resources) <- paste0("reference.", names(roles))
  cfg <- neo_config(env = character(), resources = resources, checksum = "sha256",
    paths = list(run_root = tempfile("reconstruction-runs-")),
    scientific = list(genome_assembly = "GRCh38.p14"))
  neo_project(samples = data.frame(sample_id = "sample-1", patient_id = "patient-1"),
    inputs = inputs, config = cfg)
}

reconstruction_test_row <- function(id, contig, pos, ref, alt, role = "somatic",
                                    side = NULL, block = "block:fixture", both = FALSE) {
  list(event_id = id, event_role = role, event_type = if (nchar(ref) == nchar(alt))
    if (nchar(ref) == 1L) "SNV" else "MNV" else if (nchar(ref) == 1L) "INS" else
      if (nchar(alt) == 1L) "DEL" else "INDEL",
    patient_id = "patient-1", sample_id = "sample-1", contig = contig,
    pos = as.integer(pos), ref = ref, alt = alt, phase_status = "EVALUATED",
    phase_links_json = list(), event_details_json = list(phase_assertions =
      if (is.null(side) && !both) list() else list(list(status = "INTERPRETABLE",
        genotype_state = if (both) "BOTH" else side, local_side = side,
        block_id = block))))
}
