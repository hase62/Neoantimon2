external_test_row <- function(sequence = "ACDEFGHIK", mhc = "HLA-A*02:01", ...) {
  row <- list(schema_version = "1.0.0", peptide_sequence = sequence,
    scope = "peptide_hla", mhc_key = mhc, task = "external:example",
    provider = "new-model", predictor_version = "1", model_version = "1",
    run_id = "external-run-1", metric = "score", metric_kind = "raw_score",
    metric_semantics = "raw", score = "0.9", score_direction = "higher_is_better",
    unit = "unitless", status = "evaluated", reason = "", source_row = "X:1")
  changes <- list(...); row[names(changes)] <- changes
  as.data.frame(row, stringsAsFactors = FALSE)
}

external_test_graph <- function() {
  candidate <- function(id, sample, mhc, origins) list(candidate_id = id,
    sample_id = sample, peptide_sequence = "ACDEFGHIK", mhc_key = mhc,
    hla_class = "I", origin_ids = as.list(origins), final_grade = "NOT_GRADED")
  list(candidates = list(candidate("c1", "s1", "HLA-A*02:01", c("o1", "o2")),
    candidate("c2", "s2", "HLA-A*02:01", "o3"),
    candidate("c3", "s1", "HLA-B*07:02", "o4")))
}

external_test_write <- function(x) {
  path <- tempfile(fileext = ".tsv")
  utils::write.table(x, path, sep = "\t", quote = TRUE, qmethod = "double", row.names = FALSE, na = "")
  path
}
