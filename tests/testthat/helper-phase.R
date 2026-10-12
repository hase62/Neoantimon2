phase_event_fixture <- function(name) normalizePath(testthat::test_path("..", "fixtures", "events", name))

phase_event_test_project <- function(input_names, roles, ids = basename(input_names)) {
  root <- tempfile("t007-"); dir.create(root)
  files <- list.files(testthat::test_path("..", "fixtures", "reference"), full.names = TRUE)
  file.copy(files, root)
  ref_roles <- c(genome_fasta = "genome.fa", genome_fai = "genome.fa.fai",
    annotation_gtf = "annotation.gtf", transcript_fasta = "transcripts.fa",
    protein_coding_transcript_fasta = "coding.fa", protein_fasta = "proteins.fa",
    selenocysteine_metadata = "sec.tsv", mane_summary = "mane.tsv", mane_gtf = "mane.gtf")
  resources <- lapply(names(ref_roles), function(role) list(path = file.path(root, ref_roles[[role]]),
    required = TRUE, resource_id = if (role == "annotation_gtf") Neoantimon2:::reference_scope_policy() else NULL,
    version = if (startsWith(role, "genome")) "GRCh38.p14" else
      if (startsWith(role, "mane")) "1.5" else "50", source = "synthetic"))
  names(resources) <- paste0("reference.", names(ref_roles))
  inputs <- unname(Map(function(name, role, id) list(input_id = id, input_role = role,
    sample_id = "sample-1", path = if (file.exists(name)) normalizePath(name) else phase_event_fixture(name)),
    input_names, roles, ids))
  config <- neo_config(env = character(), resources = resources, checksum = "sha256",
    paths = list(run_root = file.path(root, "runs")),
    scientific = list(genome_assembly = "GRCh38.p14"))
  neo_project(samples = data.frame(sample_id = "sample-1", patient_id = "patient-1"),
              inputs = inputs, config = config)
}
