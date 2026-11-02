# Independently authored synthetic evidence and score fixtures, NOT biology.
yield_test_scores <- function(scores = c(.1, .4, .7, 1)) {
  n <- length(scores); aa <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]];i <- seq_len(n) - 1L
  peptides <- paste0("ACDEF",aa[(i %/% 400) %% 20 + 1L],aa[(i %/% 20) %% 20 + 1L],aa[i %% 20 + 1L],"K")
  x <- data.frame(candidate_id = paste0("c",seq_len(n)), sample_id = "s1",
    peptide_sequence = peptides, mhc_key = "HLA-A*02:01", raw_score = scores,
    prediction_status = "evaluated", source = "External", provider = "test-model",
    predictor_version = "1", model_version = "1", run_id = "test-run",
    task = "external:immunogenicity", metric = "score", metric_kind = "raw_score",
    metric_semantics = "raw", unit = "unitless", score_direction = "higher_is_better",
    native_grade = rep(c("A","B","C","D","E"), length.out = n), native_rank = seq_len(n))
  x$origin_id <- paste0("o",seq_len(n));x$origin_ids <- I(lapply(x$origin_id,as.list));x
}
yield_test_profile <- function(rates = c(.01,.02,.03,.04), direction = "higher_is_better")
  list(schema_version = "experimental-yield-profile-0.1.0", id = "synthetic-assumptions", version = "1",
    mode = "scenario", endpoint = "synthetic-T-cell", assay = "synthetic-assay",
    provenance = list(source = "independently authored synthetic fixture"),
    population = "synthetic", limitations = "Not experimental observations or calibration",
    binning = list(basis = "fixed_thresholds", thresholds = c(.3,.6,.9), score_direction = direction), rates = rates)
yield_test_record <- function(id, pep = "ACDEFGHIK", hla = "HLA-A*02:01", result = "positive",
    patient = "p1", study = "synthetic-study", assay = "synthetic-assay")
  list(record_id = id, peptide_sequence = pep, mhc_key = hla, endpoint = "synthetic-T-cell",
    assay = assay, result = result, patient_id = patient, study_id = study,
    citation = "independently authored synthetic fixture; not a biological result")
yield_test_reference <- function(records = list()) list(schema_version = "experimental-reference-0.1.0",
  id = "synthetic-reference", version = "1", provenance = list(source = "synthetic"),
  patient_identity = "study_local", records = records)
yield_test_evidence <- function(x, ref, ...) neo_yield_evidence(x, ref,
  endpoint = "synthetic-T-cell", assay = "synthetic-assay", ...)
