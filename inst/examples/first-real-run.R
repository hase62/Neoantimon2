# Sourceable installed example helpers, NOT new package exports.
# Sourcing defines helpers only: no files, project, network or provider process.
first_real_run_check <- function(config, samples = NULL, hla = NULL, inputs = list(),
                                 output_dir = config$paths$run_root, include_private = FALSE) {
  utils::getFromNamespace("first_run_preflight", "Neoantimon2")(
    config, samples, hla, inputs, output_dir, include_private)
}

first_real_run_example <- function(run_root) {
  # Authored synthetic DNA/annotation, but the ordinary somatic execution route.
  # These compatibility labels do NOT turn the tiny fixture into human references.
  root <- system.file("extdata", "first-real-run", package = "Neoantimon2")
  if (!nzchar(root)) stop("Install a revision containing the First Real-Data Workflow example")
  files <- c(genome_fasta = "genome.fa", genome_fai = "genome.fa.fai",
    annotation_gtf = "annotation.gtf", transcript_fasta = "transcripts.fa",
    protein_coding_transcript_fasta = "coding.fa", protein_fasta = "proteins.fa",
    selenocysteine_metadata = "sec.tsv", mane_summary = "mane.tsv", mane_gtf = "mane.gtf")
  resources <- lapply(names(files), function(role) {
    path <- file.path(root, files[[role]])
    list(path = path, required = TRUE, type = "file",
      resource_id = if (role == "annotation_gtf")
        "gencode_primary_assembly_all_catalogs_mane_supplemental_v1" else role,
      version = if (startsWith(role, "genome")) "GRCh38.p14" else
        if (startsWith(role, "mane")) "1.5" else "50",
      source = "independently authored synthetic first-run fixture; not human resources",
      checksum = list(algorithm = "sha256", scope = "file",
        value = digest::digest(file = path, algo = "sha256", serialize = FALSE)))
  })
  names(resources) <- paste0("reference.", names(files))
  config <- Neoantimon2::neo_config(env = character(),
    paths = list(data_root = root, input_root = root, run_root = run_root),
    scientific = list(genome_assembly = "GRCh38.p14"), resources = resources,
    checksum = "sha256", peptides = list(class_i_lengths = 8:12, class_ii_lengths = 15L),
    predictors = NULL)
  list(config = config, samples = data.frame(sample_id = "demo-sample", patient_id = "demo-subject"),
    hla = list("demo-subject" = c("HLA-A*02:01", "HLA-DRA*01:01~HLA-DRB1*01:01")),
    inputs = list(list(input_id = "synthetic-somatic", input_role = "somatic_vcf",
      sample_id = "demo-sample", path = file.path(root, "somatic.vcf"),
      checksum = list(algorithm = "sha256", scope = "file",
        value = digest::digest(file = file.path(root, "somatic.vcf"), algo = "sha256", serialize = FALSE)))))
}
