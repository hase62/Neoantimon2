peptide_presenters <- function(inputs) {
  if (!length(inputs)) return(list())
  out <- lapply(sort(unique(unlist(inputs)), method = "radix"), function(raw) {
    canonical <- toupper(trimws(raw))
    chains <- strsplit(canonical, "~", fixed = TRUE)[[1L]]
    chains <- vapply(chains, function(x) {
      x <- trimws(x)
      if (!startsWith(x, "HLA-")) paste0("HLA-", x) else x
    }, "")
    valid <- grepl("^HLA-[A-Z0-9]+\\*[0-9]{2,3}(:[0-9]{2,3}){1,3}[NLSCAQ]?$", chains)
    locus <- sub("\\*.*$", "", sub("^HLA-", "", chains))
    result <- list(raw_input = raw, canonical_input = paste(chains, collapse = "~"),
      status = "UNRESOLVED", reason = "UNSUPPORTED_OR_INCOMPLETE_PRESENTER",
      hla_class = NULL, mhc_key = NULL, hla_allele = NULL,
      hla_alpha = NULL, hla_beta = NULL, pairing_status = NULL)
    if (length(chains) == 1L && all(valid) && locus %in% c("A", "B", "C", "E", "F", "G")) {
      result$status <- "RESOLVED"; result$reason <- NULL
      result$hla_class <- "I"; result$mhc_key <- chains[[1L]]
      result$hla_allele <- chains[[1L]]; result$pairing_status <- "SINGLE_CHAIN"
    } else if (length(chains) == 2L && all(valid)) {
      alpha <- which(locus %in% c("DRA", "DQA1", "DPA1"))
      beta <- which(locus %in% c("DRB1", "DRB3", "DRB4", "DRB5", "DQB1", "DPB1"))
      if (length(alpha) == 1L && length(beta) == 1L &&
          substr(locus[alpha], 1L, 2L) == substr(locus[beta], 1L, 2L)) {
        result$status <- "RESOLVED"; result$reason <- NULL
        result$hla_class <- "II"; result$hla_alpha <- chains[[alpha]]
        result$hla_beta <- chains[[beta]]
        result$mhc_key <- paste(result$hla_alpha, result$hla_beta, sep = "~")
        result$pairing_status <- "EXPLICIT_HETERODIMER"
      }
    }
    result
  })
  out
}

# This is an input boundary, not prediction evidence. Adapters must declare
# supported alphabets and return unsupported_residue without modifying input.
peptide_prediction_tasks <- function() c("mhc1_binding", "mhc1_presentation",
  "mhc2_binding", "mhc2_presentation", "processing.class1", "transport.tap",
  "stability.pmhc1", "recognition.peptide_intrinsic", "recognition.tcr_specific")

peptide_prediction_evidence <- function() data.frame(peptide_id = character(),
  sample_peptide_id = character(), candidate_id = character(), presenter = character(), task = character(),
  predictor = character(), predictor_version = character(), metric = character(),
  raw_score = numeric(), percentile_rank = numeric(), score_direction = character(),
  model_version = character(), run_id = character(), stringsAsFactors = FALSE)
