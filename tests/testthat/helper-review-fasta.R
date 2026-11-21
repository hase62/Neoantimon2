review_fasta_project <- function(sequences = c("ACDEFGHIK", "ACDEFGHIK", "LMNPQRSTV"),
    sample_ids = "s", config = NULL) {
  f <- tempfile(fileext = ".fasta")
  writeLines(unlist(lapply(seq_along(sequences), function(i)
    c(paste0(">record", i, " supplied description ", i), sequences[[i]]))), f)
  if (is.null(config)) config <- neo_config(env = character(), paths = list(run_root = tempfile()),
    peptides = list(class_i_lengths = 9L, class_ii_lengths = 15L))
  p <- neo_project(samples = data.frame(sample_id = sample_ids, patient_id = "p"),
    hla = list(p = "HLA-A*02:01"), config = config,
    inputs = unname(lapply(sample_ids, function(s) list(input_id = paste0("f-",s),
      input_role = "peptide_fasta", sample_id = s, path = f))))
  neo_run(p)
}

review_reference <- function(scores = c(0,1,1,2), direction = "lower_is_better", ...) {
  root <- tempfile(); dir.create(root); file <- file.path(root,"scores.tsv")
  writeLines(c("score", sprintf("%.17g", scores)), file)
  d <- list(reference_id = "software-fixture-v1", provider = "new-model", predictor_version = "1",
    model_version = "1", task = "external:example", metric = "score", unit = "unitless",
    score_direction = direction, mhc_key = "HLA-A*02:01", peptide_length = 9L,
    species = "NCBITaxon:9606", hla_class = "I", reference_set = "Authored software-test scores, no biological calibration",
    source = "test fixture", source_version = "1", original_sha256 = external_sha(file),
    scores_file = "scores.tsv", scores_sha256 = external_sha(file), count = length(scores), license = "repository license")
  changes <- list(...); d[names(changes)] <- changes
  path <- file.path(root,"manifest.json"); atomic_json(list(schema_version="reference-percentiles-1.0.0", distributions=list(d)), path)
  path
}
