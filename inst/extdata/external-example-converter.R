# Synthetic example only: X contains peptide, hla, value; no biological IDs.
args <- commandArgs(trailingOnly = TRUE)
arg <- function(name) {
  i <- match(name, args)
  if (is.na(i) || i == length(args)) stop(paste("Missing", name))
  args[[i + 1L]]
}
x <- read.delim(arg("--input"), colClasses = "character", check.names = FALSE)
if (!all(c("peptide", "hla", "value") %in% names(x))) stop("Expected peptide, hla, value")
n <- nrow(x)
repn <- function(s) rep(s, n)
out <- data.frame(schema_version = repn("1.0.0"), peptide_sequence = x$peptide,
  scope = repn("peptide_hla"), mhc_key = x$hla, task = repn("external:synthetic_score"),
  provider = repn("synthetic-model"), predictor_version = repn("1"), model_version = repn("demo-1"),
  run_id = repn("demo-run-1"), metric = repn("score"), metric_kind = repn("raw_score"),
  metric_semantics = repn("raw"), score = x$value, score_direction = repn("higher_is_better"),
  unit = repn("unitless"), status = repn("evaluated"), reason = repn(""),
  source_row = paste0("X:", seq_len(n)), stringsAsFactors = FALSE)
write.table(out, arg("--output"), sep = "\t", quote = TRUE, row.names = FALSE)
