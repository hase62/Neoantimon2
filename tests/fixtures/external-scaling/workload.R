# Deterministic, authored AA20+U strings and provenance; no real patient data.
external_scaling_workload <- function(n) {
  stopifnot(length(n) == 1L, n >= 0L, n == as.integer(n))
  alphabet <- strsplit("ACDEFGHIKLMNPQRSTVWYU", "", fixed = TRUE)[[1L]]
  sequence <- function(i) {
    digits <- integer(7L)
    for (j in seq_along(digits)) { digits[[j]] <- i %% length(alphabet) + 1L; i <- i %/% length(alphabet) }
    paste0("AC", paste0(alphabet[digits], collapse = ""))
  }
  keys <- c("HLA-A*02:01", "HLA-B*07:02", "HLA-DQA1*05:01~HLA-DQB1*02:01")
  candidates <- lapply(seq_len(n), function(i) list(candidate_id = sprintf("c%08d", i),
    sample_id = paste0("s", i %% 3L), peptide_sequence = sequence(i),
    mhc_key = keys[[i %% 3L + 1L]], hla_class = if (i %% 3L == 2L) "II" else "I",
    origin_ids = list(sprintf("o%08da", i), sprintf("o%08db", i)), final_grade = "NOT_GRADED"))
  origins <- unlist(lapply(candidates, function(c) lapply(c$origin_ids, function(id)
    list(origin_id = id, sample_id = c$sample_id, peptide_sequence = c$peptide_sequence))), recursive = FALSE)
  records <- c(candidates, lapply(candidates[seq_len(n) %% 10L == 0L], function(c) {
    c$origin_ids <- c$origin_ids[2L]; c
  }))
  columns <- list(schema_version = "1.0.0", peptide_sequence = "", scope = "peptide_hla", mhc_key = "",
    task = "external:scaling", provider = "synthetic-scaling", predictor_version = "1", model_version = "1",
    run_id = "scaling-v1", metric = "score", metric_kind = "raw_score", metric_semantics = "raw", score = "0.9",
    score_direction = "higher_is_better", unit = "unitless", status = "evaluated", reason = "", source_row = "")
  x <- as.data.frame(lapply(columns, rep, times = n), stringsAsFactors = FALSE)
  x$peptide_sequence <- vapply(candidates, `[[`, "", "peptide_sequence")
  x$mhc_key <- vapply(candidates, `[[`, "", "mhc_key")
  x$source_row <- if (n) paste0("X:", seq_len(n)) else character()
  x$score <- as.character(seq_len(n) %% 101L / 100)
  x$status[seq_len(n) %% 17L == 0L] <- "failed"
  x$status[seq_len(n) %% 19L == 0L] <- "unsupported"
  x$status[seq_len(n) %% 29L == 0L] <- "not_evaluated"
  x$score[x$status != "evaluated"] <- ""; x$reason[x$status != "evaluated"] <- "synthetic_reason"
  # Unknown keys stay diagnostics; some candidates intentionally have no score.
  x$peptide_sequence[seq_len(n) %% 23L == 0L] <- paste0("YYY", x$peptide_sequence[seq_len(n) %% 23L == 0L])
  intrinsic <- x[seq_len(n) %% 10L == 0L, , drop = FALSE]
  intrinsic$scope <- rep("peptide", nrow(intrinsic)); intrinsic$mhc_key <- rep("", nrow(intrinsic))
  intrinsic$metric <- rep("intrinsic", nrow(intrinsic))
  x <- Neoantimon2:::external_validate_predictions(rbind(x, intrinsic))
  policy <- list(mode = "external_rank", provider = "synthetic-scaling", run_id = "scaling-v1",
    task = "external:scaling", metric = "score", score_direction = "higher_is_better")
  list(records = records, origins = origins, events = list(), annotations = list(), predictions = x,
    selection = if (n) Neoantimon2:::external_policy(policy, x) else
      Neoantimon2:::external_policy(list(mode = "evidence_only"), x))
}
