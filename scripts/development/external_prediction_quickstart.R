# Run from the repository. Original synthetic inputs survive outside runs.
args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default) {
  i <- match(name, args)
  if (is.na(i)) default else {
    if (i == length(args)) stop(paste("Missing value for", name))
    args[[i + 1L]]
  }
}
repo <- normalizePath(".", mustWork = TRUE)
if (!file.exists(file.path(repo, "DESCRIPTION"))) stop("Run from the Neoantimon2 repository")
if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(repo, quiet = TRUE) else library(Neoantimon2)
data_root <- path.expand(arg("--data-root", "~/codex/data/Neoantimon2"))
run_root <- path.expand(arg("--run-root", "~/codex/runs/Neoantimon2"))
input_root <- file.path(data_root, "external-interop-synthetic-v1")
dir.create(input_root, recursive = TRUE, showWarnings = FALSE)
preserve <- function(name, bytes) {
  path <- file.path(input_root, name)
  if (file.exists(path)) {
    old <- readBin(path, "raw", n = file.info(path)$size)
    if (!identical(old, bytes)) stop(paste("Existing original input differs; refusing overwrite:", path))
  } else writeBin(bytes, path)
  path
}
model <- preserve("shared-model.gtf", charToRaw(paste0(
  'chr1\tx\texon\t1\t7\t.\t+\t.\ttranscript_id "shared";\n',
  'chr1\tx\texon\t21\t33\t.\t+\t.\ttranscript_id "shared";\n')))
source_fixture <- file.path(repo, "tests", "fixtures", "source", "splice-structure.tsv")
splice <- preserve("splice-structure.tsv", readBin(source_fixture, "raw", n = file.info(source_fixture)$size))
reference <- file.path(repo, "tests", "fixtures", "long_read", "reference")
roles <- c(genome_fasta = "genome.fa", genome_fai = "genome.fa.fai", annotation_gtf = "annotation.gtf",
  transcript_fasta = "transcripts.fa", protein_coding_transcript_fasta = "coding.fa", protein_fasta = "proteins.fa",
  selenocysteine_metadata = "sec.tsv", mane_summary = "mane.tsv", mane_gtf = "mane.gtf")
resources <- lapply(names(roles), function(role) list(path = file.path(reference, roles[[role]]),
  required = TRUE, resource_id = if (role == "annotation_gtf") Neoantimon2:::reference_scope_policy() else NULL,
  version = if (startsWith(role, "genome")) "GRCh38.p14" else if (startsWith(role, "mane")) "1.5" else "50",
  source = "synthetic"))
names(resources) <- paste0("reference.", names(roles))
cfg <- neo_config(env = character(), paths = list(repository = repo, data_root = data_root, run_root = run_root),
  resources = resources, checksum = "sha256", scientific = list(genome_assembly = "GRCh38.p14"),
  peptides = list(class_i_lengths = list(6L), class_ii_lengths = list(13L)))
p <- neo_run(neo_project(samples = data.frame(sample_id = "sample-1", patient_id = "patient-1"),
  hla = list("patient-1" = "HLA-A*02:01"), config = cfg,
  inputs = list(list(input_id = "long-read", input_role = "long_read_gtf", sample_id = "sample-1", path = model, version = "synthetic-v1"),
    list(input_id = "splice", input_role = "splice_structure", sample_id = "sample-1", path = splice, version = "synthetic-v1"))))
native <- neo_results(p)
bundle <- neo_export_bundle(p)
multi <- Filter(function(c) length(c$origin_ids) > 1L, bundle$graph$candidates)
stopifnot(length(multi) > 0L)
priority_path <- Neoantimon2:::somatic_priority_path(p)
originals <- unique(c(Neoantimon2:::manifest_file(p$path), Neoantimon2:::peptide_stage_path(p),
  Neoantimon2:::runtime_snapshot_path(p), if (file.exists(priority_path)) priority_path,
  list.files(file.path(p$path, "results"), recursive = TRUE, full.names = TRUE)))
hashes <- vapply(originals, Neoantimon2:::external_sha, "")
native_explanation <- neo_explain(p, candidate_id = multi[[1L]]$candidate_id)
queries <- read.delim(file.path(bundle$path, "pmhc_queries.tsv"), colClasses = "character")
# Simulate X reordered and deduplicated by an external sequence/HLA model.
raw <- data.frame(peptide = queries$peptide_sequence, hla = queries$mhc_key,
  value = seq_len(nrow(queries)) / max(1L, nrow(queries)))
raw <- raw[rev(seq_len(nrow(raw))), , drop = FALSE]
staged <- tempfile(tmpdir = p$path, fileext = ".tsv")
write.table(raw, staged, sep = "\t", quote = TRUE, row.names = FALSE)
raw_path <- preserve("output-X.tsv", readBin(staged, "raw", n = file.info(staged)$size)); unlink(staged)
converter <- file.path(repo, "inst", "extdata", "external-example-converter.R")
converter_path <- preserve("converter.R", readBin(converter, "raw", n = file.info(converter)$size))
policy <- list(mode = "external_rank", provider = "synthetic-model", run_id = "demo-run-1",
  task = "external:synthetic_score", metric = "score", score_direction = "higher_is_better")
a <- neo_analyze_external(p, raw_path, policy, converter = converter_path)
r <- neo_external_results(a)
explanation <- neo_external_explain(a, multi[[1L]]$candidate_id)
stopifnot(length(explanation$candidates[[1L]]$origins) >= 2L,
  all(r$external_rank_status == "USER_SELECTED_EXTERNAL_RANK"),
  all(r$final_grade == "NOT_GRADED"), identical(neo_results(p), native),
  identical(neo_explain(p, candidate_id = multi[[1L]]$candidate_id), native_explanation),
  identical(vapply(originals, Neoantimon2:::external_sha, ""), hashes))
# The converted table is directly importable without executing any converter.
direct <- neo_analyze_external(p, file.path(a$path, "normalized.tsv"), policy)
stopifnot(identical(neo_external_results(direct), r),
  identical(neo_external_results(a$analysis_id, neo_project(path = p$path)), r))
summary <- list(status = "PASS", run_path = p$path, analysis_id = a$analysis_id,
  direct_analysis_id = direct$analysis_id, candidates = nrow(r),
  multi_origin_candidate = multi[[1L]]$candidate_id,
  origin_count = length(explanation$candidates[[1L]]$origins),
  native_hashes_unchanged = TRUE, original_priority_present = file.exists(priority_path),
  native_artifact_hashes = as.list(stats::setNames(hashes, substring(originals, nchar(p$path) + 2L))),
  import_report = explanation$import_report)
Neoantimon2:::atomic_json(summary, file.path(p$path, "quickstart-summary.json"))
cat(jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE), "\n")
