# Multiple patients and multiple samples

[Getting started](getting-started.md) · [Inputs/configuration](inputs.md) ·
[Longitudinal workflows](workflows.md#longitudinal-and-multi-sample-analysis)

One somatic project can declare several patients and samples. Each sample has
one explicit patient, each patient's HLA is declared once, and each input file is
bound to a sample. `neo_run()` processes all these declared inputs in one run.
This guide uses **separate single-sample VCF files**. A VCF containing multiple
sample columns is rejected; this guide does not add support for it.

## Run the examples

Install the [public v0.3.0 package](getting-started.md#install-and-run-the-offline-example).
Then execute the shared setup below followed by any scenario. Each scenario is
complete after that setup and can run independently of the other scenarios.
All R code blocks on this page, in order, are also available in
[`inst/examples/multi-sample.R`](../../inst/examples/multi-sample.R). From the root
of a source checkout containing that file, run
`source("inst/examples/multi-sample.R", local = TRUE)` to execute all four.
The script uses only the existing v0.2.0 API and installed synthetic references;
it requires no internet or licensed predictor after package installation.

These are synthetic technical examples, not patient data or biological
validation. All three files deliberately contain the **same** chr1:4 G>A variant.
This makes it possible to check that identical sequences and variant keys do not
merge different samples. HLA declarations are minimal test inputs, not complete
clinical genotypes. P2 shares HLA-A*02:01 with P1 and additionally declares
HLA-B*07:02, allowing patient-specific HLA resolution to be checked.

## Shared offline setup

The setup creates a fresh temporary directory with the layout below. Scenario
runs are generated under `runs/`; references are read from the installed package.
The directory is printed for inspection and is not durable storage.

```text
<example_root>/
  inputs/
    S1.vcf  # header sample S1
    S2.vcf  # header sample S2
    S3.vcf  # header sample S3
  runs/     # generated run directories, one per executed scenario
<installed Neoantimon2>/extdata/quickstart/  # nine tiny reference files
```

The synthetic references have release/build labels for contract testing; they
are not the human reference genome. Their
[provenance](../../inst/extdata/quickstart/README.md) is unchanged. The 3-aa/5-aa
windows below are test settings. For your own data use the
[complete real-reference profile and production lengths](inputs.md).

```r
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
```

## Scenario 1: one patient, one sample

```text
Patient P1
  Sample S1 -> inputs/S1.vcf
```

```r
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
```

Expect one P1/S1 event, seven P1/S1 origins and four P1/S1 candidates for HLA-A*02:01. `events_1`, `origins_1` and `results_1` are sample-filtered data frames.

## Scenario 2: one patient, two samples

```text
Patient P1
  Sample S1 -> inputs/S1.vcf
  Sample S2 -> inputs/S2.vcf
```

```r
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
```

Expect one event, seven origins and four HLA-A*02:01 candidates **per sample**. P1's HLA is declared once and used by both samples. `results_2_S1` and `results_2_S2` stay separate despite identical sequences. `neo_multisample(p2)` is empty because no longitudinal settings were requested.

## Scenario 3: two patients, one sample each

```text
Patient P1
  Sample S1 -> inputs/S1.vcf
Patient P2
  Sample S2 -> inputs/S2.vcf
```

```r
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
```

Expect one event and seven origins per sample. S1/P1 has four HLA-A*02:01 candidates; S2/P2 has eight candidates across its two declared alleles. HLA-B*07:02 must not appear in S1. Shared HLA-A*02:01 and peptide sequences do not merge the two patients' candidate records.

## Scenario 4: two patients, three samples

```text
Patient P1
  Sample S1 -> inputs/S1.vcf
  Sample S2 -> inputs/S2.vcf
Patient P2
  Sample S3 -> inputs/S3.vcf
```

```r
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
```

Expect three events and 21 origins, preserving the declared patient/sample relationships. S1/P1 and S2/P1 each have four candidates; S3/P2 has eight. The total is 16 candidates. `results_4_P1` selects both of P1's samples through the supported `sample_id` filter; `neo_results()` does not take a `patient_id` argument.

For every scenario, prediction is `NOT_RUN`, grading and final grade are
`NOT_GRADED`: no predictor or quantitative somatic evidence was requested.
Completion is an execution result, not evidence of immunogenicity. Missing
predictions are not zero or evaluated negatives. Counts describe only these
fixed synthetic fixtures; real inputs need not yield candidates.

## Identity, validation and result selection

- `samples` is a data frame or unnamed list of records. Each record needs a
  nonempty string `sample_id` and `patient_id`; sample IDs are unique across the
  project, including across patients. Repeated patient IDs group samples.
- `hla` is a named list keyed by patients declared in `samples`. A sample ID or
  unknown patient is not a valid HLA key. A patient's genotype is reused across
  its samples. HLA is not inferred; omitted/empty or unresolved declarations can
  leave events and origins intact while producing no candidates for that patient.
  Do not interpret that as a negative prediction. Class II requires an explicit
  supported alpha~beta pair; see [HLA formats](inputs.md#genomic-and-hla-input).
- `inputs` is an unnamed list of records. Input IDs must be unique and every
  `sample_id` must exist in `samples`. A relative path resolves under the
  configured `input_root`. Missing required input files are rejected at project
  construction. Unknown declaration keys are rejected.
- `neo_project()` validates declarations and binds file integrity; `neo_run()`
  parses VCF contents. A VCF sample column must match its authoritative input
  binding. A mismatched header or multiple sample columns causes execution to
  fail. Filenames do not assign sample or patient identity.
- Ordinary somatic contexts are created per patient. Event, origin and candidate
  records retain patient/sample/context associations. Shared event, peptide or
  candidate *keys* describe shared biological content; they do not authorize
  merging sample-specific IDs or dropping origins.
- Use `neo_events(p, sample_id = "S1")`, `neo_origins(p, sample_id = "S1")` and
  `neo_results(p, sample_id = "S1")` for sample-specific output. To select a
  patient, first obtain its sample IDs from your `samples` declaration, as in
  Scenario 4. There is no patient filter on these three accessors.

A project declaration is not a guarantee that every event can be reconstructed,
or that every optional workflow accepts every cohort arrangement. Unsupported
sequence contexts retain their documented status; see
[scientific limitations](methods-and-limitations.md).

## Ordinary samples versus longitudinal analysis

The four scenarios require one ordinary project/run each. They do not require
one run per sample and do not establish a temporal order or a trajectory.
Adding a sample `timepoint` string alone does not enable longitudinal analysis.
`neo_multisample()` is **not** a general inventory of project samples: the normal
runtime leaves its summary empty without a configured longitudinal overlay.
Use the sample-specific accessors above for ordinary projects.

The [longitudinal workflow](workflows.md#longitudinal-and-multi-sample-analysis)
requires a `longitudinal_settings` input, explicit subjects, sample contexts,
within-subject relations and evaluability declarations as needed. Its overlay
can operate on the current project's sample graphs, or an aggregate project can
bind a completed run for **every** declared sample using `sample_runs`. The
existing workflow example demonstrates the latter route; inspect the source runs
for their individual events/origins. Comparable explicit time axes establish
chronology; filenames, row order and informal time labels do not.

Subject-scoped trajectories require the relevant multiple-tumor context. A
single sample or matched-normal-only comparison does not create that overlay.
Failure to observe an event/candidate is not an evaluated absence without an
explicit adequate observation opportunity. External prediction import also has
its own [integration boundaries](external-predictions.md); it does not
automatically populate every longitudinal or other optional module.
