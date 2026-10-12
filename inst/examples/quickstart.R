# SYNTHETIC / TEST-SCALE / NON-CLINICAL. No network or external predictor.
# Reference bytes: committed tests/fixtures/reconstruction/ at
# 635823fc4a41308d010396fdababb2ca5d88878b; see extdata/quickstart/README.md.
library(Neoantimon2)

example_root <- tempfile("neoantimon2-example-")
dir.create(example_root)
input_root <- file.path(example_root, "inputs")
dir.create(input_root)
reference_dir <- system.file("extdata", "quickstart", package = "Neoantimon2")
stopifnot(nzchar(reference_dir))
reference_files <- c(
  genome_fasta = "genome.fa", genome_fai = "genome.fa.fai",
  annotation_gtf = "annotation.gtf", transcript_fasta = "transcripts.fa",
  protein_coding_transcript_fasta = "coding.fa", protein_fasta = "proteins.fa",
  selenocysteine_metadata = "sec.tsv", mane_summary = "mane.tsv",
  mane_gtf = "mane.gtf"
)
resources <- lapply(names(reference_files), function(role) {
  path <- file.path(reference_dir, reference_files[[role]])
  list(path = path, required = TRUE, type = "file",
    resource_id = if (role == "annotation_gtf")
      "gencode_primary_assembly_all_catalogs_mane_supplemental_v1" else role,
    version = if (startsWith(role, "genome")) "GRCh38.p14" else
      if (startsWith(role, "mane")) "1.5" else "50",
    source = "synthetic reconstruction fixture",
    checksum = list(algorithm = "sha256", scope = "file",
      value = digest::digest(file = path, algo = "sha256", serialize = FALSE)))
})
names(resources) <- paste0("reference.", names(reference_files))

vcf <- file.path(input_root, "synthetic.vcf")
writeLines(c(
  "##fileformat=VCFv4.2",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample-1",
  "chr1\t4\t.\tG\tA\t.\tPASS\t.\tGT\t0/1"
), vcf)
cfg <- neo_config(env = character(),
  paths = list(data_root = input_root, input_root = input_root,
    run_root = file.path(example_root, "runs")),
  scientific = list(genome_assembly = "GRCh38.p14"),
  checksum = "sha256", resources = resources,
  peptides = list(class_i_lengths = 3L, class_ii_lengths = 5L))
p <- neo_project(
  samples = data.frame(sample_id = "sample-1", patient_id = "synthetic-subject"),
  config = cfg, hla = list("synthetic-subject" = "HLA-A*02:01"),
  inputs = list(list(input_id = "synthetic-variant", input_role = "somatic_vcf",
    sample_id = "sample-1", path = vcf,
    checksum = list(algorithm = "sha256", scope = "file",
      value = digest::digest(file = vcf, algo = "sha256", serialize = FALSE)))))
p <- neo_run(p)
events <- neo_events(p)
origins <- neo_origins(p)
results <- neo_results(p)
multisample <- neo_multisample(p)
explanation <- neo_explain(p, candidate_id = results$candidate_id[[1L]])
reopened <- neo_project(path = p$path)
stopifnot(nrow(events) == 1L, nrow(origins) == 7L, nrow(results) == 4L,
  identical(unique(results$prediction_status), "NOT_RUN"),
  identical(neo_events(reopened), events),
  identical(neo_origins(reopened), origins),
  identical(neo_results(reopened), results),
  identical(neo_results(neo_run(reopened)), results))
cat("Synthetic quick start: 1 event, 7 origins, 4 candidates; prediction NOT_RUN; reopen identical.\n")
print(results[, c("sample_id", "ligand_sequence", "mhc_key",
  "prediction_status", "grading_status", "final_grade")])
cat("Example files and run directory:", example_root, "\n")
