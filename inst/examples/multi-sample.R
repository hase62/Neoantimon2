# SYNTHETIC / TEST-SCALE / NON-CLINICAL. No network or predictor required.
# Run from a source checkout with Neoantimon2 >= 0.2.0 installed, or source
# this file from an installed package that includes it.
# The same code is printed in docs/user/multi-sample-examples.md.

# --- Shared setup ---
library(Neoantimon2)

example_root <- tempfile("neoantimon2-cohort-")
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

for (sample_id in c("S1", "S2", "S3")) {
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>",
    paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", sample_id),
    "chr1\t4\t.\tG\tA\t.\tPASS\t.\tGT\t0/1"
  ), file.path(input_root, paste0(sample_id, ".vcf")))
}
cfg <- neo_config(env = character(),
  paths = list(data_root = input_root, input_root = input_root,
    run_root = file.path(example_root, "runs")),
  scientific = list(genome_assembly = "GRCh38.p14"),
  checksum = "sha256", resources = resources,
  peptides = list(class_i_lengths = 3L, class_ii_lengths = 5L))
cat("Synthetic input and run root:", example_root, "\n")

# --- Scenario 1: one patient, one sample ---
samples_1 <- data.frame(sample_id = "S1", patient_id = "P1")
hla_1 <- list(P1 = "HLA-A*02:01")
inputs_1 <- list(
  list(input_id = "S1-somatic", input_role = "somatic_vcf",
    sample_id = "S1", path = "S1.vcf"))
p1 <- neo_project(samples = samples_1, hla = hla_1, inputs = inputs_1,
  config = cfg, project_id = "one-patient-one-sample")
p1 <- neo_run(p1)
events_1 <- neo_events(p1, sample_id = "S1")
origins_1 <- neo_origins(p1, sample_id = "S1")
results_1 <- neo_results(p1, sample_id = "S1")
print(results_1[, c("patient_id", "sample_id", "ligand_sequence", "mhc_key",
  "prediction_status", "grading_status")])

# --- Scenario 2: one patient, two samples ---
samples_2 <- data.frame(sample_id = c("S1", "S2"), patient_id = c("P1", "P1"))
hla_2 <- list(P1 = "HLA-A*02:01")
inputs_2 <- list(
  list(input_id = "S1-somatic", input_role = "somatic_vcf",
    sample_id = "S1", path = "S1.vcf"),
  list(input_id = "S2-somatic", input_role = "somatic_vcf",
    sample_id = "S2", path = "S2.vcf"))
p2 <- neo_project(samples = samples_2, hla = hla_2, inputs = inputs_2,
  config = cfg, project_id = "one-patient-two-samples")
p2 <- neo_run(p2)
events_2 <- neo_events(p2)
origins_2 <- neo_origins(p2)
results_2_S1 <- neo_results(p2, sample_id = "S1")
results_2_S2 <- neo_results(p2, sample_id = "S2")
print(neo_results(p2)[, c("patient_id", "sample_id", "mhc_key")])
stopifnot(nrow(neo_multisample(p2)) == 0L)

# --- Scenario 3: two patients, one sample each ---
samples_3 <- data.frame(sample_id = c("S1", "S2"), patient_id = c("P1", "P2"))
hla_3 <- list(P1 = "HLA-A*02:01", P2 = c("HLA-A*02:01", "HLA-B*07:02"))
inputs_3 <- list(
  list(input_id = "S1-somatic", input_role = "somatic_vcf",
    sample_id = "S1", path = "S1.vcf"),
  list(input_id = "S2-somatic", input_role = "somatic_vcf",
    sample_id = "S2", path = "S2.vcf"))
p3 <- neo_project(samples = samples_3, hla = hla_3, inputs = inputs_3,
  config = cfg, project_id = "two-patients-two-samples")
p3 <- neo_run(p3)
events_3 <- neo_events(p3)
origins_3 <- neo_origins(p3)
results_3_S1 <- neo_results(p3, sample_id = "S1")
results_3_S2 <- neo_results(p3, sample_id = "S2")
print(neo_results(p3)[, c("patient_id", "sample_id", "mhc_key")])

# --- Scenario 4: two patients, three samples ---
samples_4 <- data.frame(sample_id = c("S1", "S2", "S3"),
  patient_id = c("P1", "P1", "P2"))
hla_4 <- list(P1 = "HLA-A*02:01", P2 = c("HLA-A*02:01", "HLA-B*07:02"))
inputs_4 <- list(
  list(input_id = "S1-somatic", input_role = "somatic_vcf",
    sample_id = "S1", path = "S1.vcf"),
  list(input_id = "S2-somatic", input_role = "somatic_vcf",
    sample_id = "S2", path = "S2.vcf"),
  list(input_id = "S3-somatic", input_role = "somatic_vcf",
    sample_id = "S3", path = "S3.vcf"))
p4 <- neo_project(samples = samples_4, hla = hla_4, inputs = inputs_4,
  config = cfg, project_id = "two-patients-three-samples")
p4 <- neo_run(p4)
events_4_S2 <- neo_events(p4, sample_id = "S2")
origins_4_S2 <- neo_origins(p4, sample_id = "S2")
results_4_S1 <- neo_results(p4, sample_id = "S1")
results_4_S2 <- neo_results(p4, sample_id = "S2")
results_4_S3 <- neo_results(p4, sample_id = "S3")
results_4_P1 <- neo_results(p4, sample_id = samples_4$sample_id[samples_4$patient_id == "P1"])
print(neo_results(p4)[, c("patient_id", "sample_id", "mhc_key")])
