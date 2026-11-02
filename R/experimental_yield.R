yield_channel_fields <- function() c("source", "provider", "predictor_version", "model_version",
  "run_id", "task", "metric", "metric_kind", "metric_semantics", "unit", "score_direction")
yield_model_fields <- function() setdiff(yield_channel_fields(), c("source", "run_id"))

yield_profile <- function(profile, identity) {
  if (identical(profile, "improve-exploratory-v0.1"))
    profile <- system.file("extdata", "yield", "improve-scenario.json", package = "Neoantimon2")
  obj <- yield_object(profile, "experimental-yield-profile-0.1.0"); p <- obj$data
  if (is.null(p$mode) || !p$mode %in% c("scenario", "empirical")) stop("Yield mode must be scenario or empirical")
  for (f in c("endpoint", "assay", "population", "limitations")) yield_text(p[[f]], paste("Profile", f))
  if (is.null(p$binning$basis) || !p$binning$basis %in% c("fixed_thresholds", "cohort_relative")) stop("Explicit binning basis is required")
  same <- !is.null(p$model) && all(vapply(yield_model_fields(), function(f)
    is.character(p$model[[f]]) && length(p$model[[f]]) == 1L &&
      !is.na(p$model[[f]]) && nzchar(p$model[[f]]) && identical(p$model[[f]], identity[[f]]), TRUE))
  if (p$mode == "empirical") {
    if (!same) stop("Empirical calibration predictor/model/metric/unit/direction mismatch")
    if (p$binning$basis != "fixed_thresholds") stop("Empirical profiles require frozen reference thresholds")
    if (length(p$bins) != 4L) stop("Empirical profile requires four tested bins")
    rates <- vapply(seq_len(4L), function(i) {
      b <- p$bins[[i]]; yield_number(b$bin, "bin", 1, 4, TRUE)
      if (b$bin != i) stop("Empirical bins must be ordered 1 to 4")
      for (f in c("tested", "positive", "negative")) yield_number(b[[f]], f, 0, integer = TRUE)
      if (b$tested < 1 || b$positive + b$negative != b$tested) stop("Inconsistent empirical tested/positive/negative counts")
      yield_number(b$rate, "empirical rate", 0, 1)
      if (abs(b$rate - b$positive / b$tested) > 1e-12) stop("Empirical rate differs from observed counts")
      b$rate
    }, 0)
  } else {
    rates <- unlist(p$rates, use.names = FALSE)
    if (!is.numeric(rates) || length(rates) != 4L || anyNA(rates) ||
        any(!is.finite(rates) | rates < 0 | rates > 1)) stop("Scenario requires four valid probabilities in [0,1]")
  }
  if (p$binning$basis == "fixed_thresholds") {
    t <- unlist(p$binning$thresholds, use.names = FALSE)
    if (!is.numeric(t) || length(t) != 3L || any(!is.finite(t))) stop("Require three finite frozen raw-score thresholds")
    if (!same && !is.null(p$model)) stop("Cross-model scenario must use exploratory cohort-relative ranks, not incompatible frozen score thresholds")
    if (is.null(p$model) && !identical(p$binning$score_direction, identity$score_direction))
      stop("Scenario threshold score_direction does not match the actual predictor")
    q <- if (identity$score_direction == "higher_is_better") t else -t
    if (any(diff(q) < 0)) stop("Frozen thresholds must follow increasing prediction quality (30th,60th,90th)")
  }
  list(data = p, artifact = obj$artifact, rates = rates, model_match = same,
    interpretation = if (p$mode == "empirical") "empirical_model_endpoint_matched; transportability_not_established" else
      if (!is.null(p$model) && !same) "EXPLORATORY_CROSS_MODEL_ASSUMPTION_NOT_CALIBRATED" else "USER_SCENARIO_NOT_CALIBRATED")
}

yield_score_input <- function(x, analysis, channel) {
  if (inherits(x, "neo_project")) x <- neo_priorities(x, analysis = analysis)
  else if (!is.null(analysis)) stop("analysis is supported only with a neo_project")
  required <- c("candidate_id", "sample_id", "peptide_sequence", "mhc_key", "raw_score",
    "prediction_status", yield_channel_fields())
  if (!is.data.frame(x) || !all(required %in% names(x)))
    stop("Yield requires a neo_project or neo_priorities-compatible data frame; missing: ",
      paste(setdiff(required, names(x)), collapse = ", "))
  if (!nrow(x)) stop("No candidate rows supplied")
  for (f in c("candidate_id", "sample_id", "prediction_status")) {
    if (!is.character(x[[f]]) || anyNA(x[[f]]) || any(!nzchar(x[[f]]))) stop("Invalid ", f)
  }
  if (length(unique(x$sample_id)) != 1L) stop("Select one sample per yield plan; never pool samples silently")
  if ("patient_id" %in% names(x)) {
    patients <- unique(x$patient_id[!is.na(x$patient_id) & nzchar(x$patient_id)])
    if (length(patients) > 1L) stop("Conflicting patient identities within sample")
  }
  id_keys <- paste(x$sample_id, yield_peptide(x$peptide_sequence), yield_hla(x$mhc_key), sep = "\037")
  for (id in unique(x$candidate_id)) if (length(unique(id_keys[x$candidate_id == id])) != 1L)
    stop("candidate_id maps to conflicting sample/peptide/HLA identities")
  if (!is.numeric(x$raw_score) || any(is.nan(x$raw_score)) || any(is.infinite(x$raw_score)))
    stop("raw_score must be finite numeric or genuinely missing NA")
  fields <- yield_channel_fields()
  for (f in fields) if (!is.character(x[[f]])) stop("Channel identity must be character: ", f)
  complete <- Reduce(`&`, lapply(fields, function(f) !is.na(x[[f]]) &
    (f %in% c("predictor_version", "model_version") | nzchar(x[[f]]))))
  available <- x[complete, fields, drop = FALSE]
  if (!is.null(channel)) {
    if (!is.list(channel) || !length(channel) || is.null(names(channel)) ||
        any(!names(channel) %in% fields)) stop("channel must select explicit predictor identity fields")
    for (f in names(channel)) {
      yield_text(channel[[f]], paste("channel", f))
      available <- available[available[[f]] == channel[[f]], , drop = FALSE]
    }
  }
  identities <- unique(available)
  if (nrow(identities) != 1L) stop("Select exactly one compatible predictor/model/run/metric channel")
  identity <- as.list(identities[1, , drop = FALSE])
  if (!identity$score_direction %in% c("higher_is_better", "lower_is_better")) stop("Undefined score direction")
  match <- Reduce(`&`, lapply(fields, function(f) !is.na(x[[f]]) & x[[f]] == identity[[f]]))
  # Unpredicted candidates with no channel cannot be silently dropped from pool size.
  no_prediction <- !complete & (is.na(x$provider) | !nzchar(x$provider)) &
    x$prediction_status %in% c("NOT_RUN", "NOT_PREDICTED", "not_run")
  keep <- match | no_prediction; selected <- x[keep, , drop = FALSE]
  selected <- yield_pair_input(selected)
  pairs <- paste(selected$sample_id, selected$peptide_sequence, selected$mhc_key, sep = "\037")
  observations <- lapply(unique(pairs), function(key) {
    at <- which(pairs == key); z <- selected[at, , drop = FALSE]
    evaluated <- z$prediction_status == "evaluated" & !is.na(z$raw_score)
    if (length(unique(z$prediction_status)) > 1L || length(unique(z$raw_score)) > 1L)
      stop("Conflicting duplicate peptide-HLA measurements; select one source channel")
    if (any(evaluated) && identity$unit == "nM" && any(z$raw_score[evaluated] < 0)) stop("Negative affinity is invalid")
    ids <- unique(as.character(unlist(c(z[intersect(c("origin_id", "origin_ids"), names(z))]), use.names = FALSE)))
    ids <- ids[!is.na(ids) & nzchar(ids)]
    list(row = z[1, , drop = FALSE], candidate_ids = sort(unique(z$candidate_id), method = "radix"),
      origin_ids = sort(ids, method = "radix"), source_rows = z,
      evaluated = all(evaluated), input_position = min(at))
  })
  list(observations = observations, identity = identity, input = x,
    excluded_channel_rows = sum(!keep), selected_input_rows = nrow(selected))
}

yield_bins <- function(scores, good, identity, profile) {
  b <- rep(NA_integer_, length(scores)); q <- if (identity$score_direction == "lower_is_better") -scores else scores
  if (profile$data$binning$basis == "fixed_thresholds") {
    cuts <- unlist(profile$data$binning$thresholds)
    if (identity$score_direction == "lower_is_better") cuts <- -cuts
    b[good] <- vapply(q[good], function(v) 1L + sum(v > cuts), 0L)
  } else if (any(good)) {
    # Mid-rank percentile over distinct scored peptide-HLA pairs, not Origin rows.
    u <- (rank(q[good], ties.method = "average") - .5) / sum(good)
    b[good] <- vapply(u, function(v) 1L + sum(v > c(.3, .6, .9)), 0L)
  }
  b
}

#' Optional assumption-based experimental yield, separate from native Grades
#' @param x A completed neo_project or a data frame returned by neo_priorities.
#' @param profile Versioned experimental-yield-profile-0.1.0 JSON path/list, or
#'   "improve-exploratory-v0.1" for the explicitly uncalibrated illustrative scenario.
#' @param reference Versioned experimental-reference-0.1.0 path/list, or NULL.
#' @param analysis Optional existing external analysis handle or ID, for project input.
#' @param channel Named list selecting one provider/model/run/task/metric channel.
#' @param unit Unique synthesized peptide proxy or exact peptide-HLA assays.
#' @param ordering Direction-aware score order, or first appearance in the input.
#' @param primary_hla Deterministic best-score proxy (ties lexical HLA) or lexical HLA.
#' @param min_records,min_patients Provisional C evidence thresholds.
#' @return A neo_experimental_yield list containing units, pair_evidence, curve,
#'   targets, summary, composition, explicit profile/reference and original provenance.
#' @export
neo_experimental_yield <- function(x, profile, reference = NULL, analysis = NULL,
    channel = NULL, unit = c("peptide", "peptide_hla"), ordering = c("score", "input"),
    primary_hla = c("best_score", "lexicographic"), min_records = 100L, min_patients = 3L) {
  unit <- match.arg(unit); ordering <- match.arg(ordering); primary_hla <- match.arg(primary_hla)
  input <- yield_score_input(x, analysis, channel); obs <- input$observations
  p <- yield_profile(profile, input$identity)
  pairs <- do.call(rbind, lapply(obs, `[[`, "row")); rownames(pairs) <- NULL
  pairs$yield_evaluated <- vapply(obs, `[[`, TRUE, "evaluated")
  pairs$score_bin <- yield_bins(pairs$raw_score, pairs$yield_evaluated, input$identity, p)
  pairs$probability <- p$rates[pairs$score_bin]
  pairs$probability_status <- ifelse(is.na(pairs$probability), "UNSCORED_NO_ASSUMED_NEGATIVE", p$interpretation)
  # Original measurements and links are never replaced by these calculated values.
  pairs$all_candidate_ids <- I(lapply(obs, `[[`, "candidate_ids"))
  pairs$all_origin_ids <- I(lapply(obs, `[[`, "origin_ids"))
  pairs$source_rows <- I(lapply(obs, `[[`, "source_rows"))
  pairs$input_position <- vapply(obs, `[[`, 0L, "input_position")
  pairs <- neo_yield_evidence(pairs, reference, p$data$endpoint, p$data$assay,
    min_records = min_records, min_patients = min_patients)
  q <- if (input$identity$score_direction == "lower_is_better") -pairs$raw_score else pairs$raw_score
  q[!pairs$yield_evaluated] <- NA_real_
  groups <- if (unit == "peptide") split(seq_len(nrow(pairs)), pairs$peptide_sequence) else
    as.list(seq_len(nrow(pairs)))
  chosen <- vapply(groups, function(at) {
    ord <- if (primary_hla == "best_score") order(-q[at], pairs$mhc_key[at],
      pairs$candidate_id[at], na.last = TRUE, method = "radix") else
      order(pairs$mhc_key[at], pairs$candidate_id[at], method = "radix")
    at[ord[1L]]
  }, 0L)
  units <- pairs[chosen, , drop = FALSE]; units$alternative_pairs <- I(lapply(groups,
    function(at) pairs[at, , drop = FALSE]))
  units$all_candidate_ids <- I(lapply(groups, function(at) sort(unique(unlist(pairs$all_candidate_ids[at])), method = "radix")))
  units$all_origin_ids <- I(lapply(groups, function(at) sort(unique(unlist(pairs$all_origin_ids[at])), method = "radix")))
  units$synthesized_peptide <- units$peptide_sequence
  units$unit_definition <- unit; units$primary_HLA_proxy <- unit == "peptide"
  units$primary_HLA_rule <- if (unit == "peptide") primary_hla else "not_applicable"
  order_at <- if (ordering == "score") order(-q[chosen], units$peptide_sequence,
    units$mhc_key, method = "radix", na.last = TRUE) else order(units$input_position, method = "radix")
  units <- units[order_at, , drop = FALSE]; rownames(units) <- NULL
  n <- nrow(units); rates <- units$probability
  duplicate_peptide <- duplicated(units$peptide_sequence)
  missing <- is.na(rates); unknown_prefix <- cumsum(missing) > 0L
  repeated_prefix <- cumsum(duplicate_peptide) > 0L
  expected <- cumsum(rates); prob <- -expm1(cumsum(log1p(-rates)))
  prob[unknown_prefix | repeated_prefix] <- NA_real_
  curve <- data.frame(N_tested = 0:n,
    distinct_synthesized_peptides = c(0L, cumsum(!duplicate_peptide)),
    expected_positive_units = c(0, expected), probability_at_least_one = c(0, prob),
    probability_status = c("EMPTY_POOL", ifelse(unknown_prefix, "UNSCORED_PREFIX_UNAVAILABLE",
      ifelse(repeated_prefix, "SHARED_PEPTIDE_HLA_JOINT_PROBABILITY_UNAVAILABLE", "INDEPENDENCE_ASSUMPTION"))))
  targets <- do.call(rbind, lapply(c(.5, .8, .9, .95), function(target) {
    at <- which(!is.na(curve$probability_at_least_one) & curve$probability_at_least_one >= target)
    hit <- if (length(at)) at[1L] else NA_integer_
    data.frame(target = target, minimum_N = if (is.na(hit)) NA_integer_ else curve$N_tested[hit],
      distinct_synthesized_peptides = if (is.na(hit)) NA_integer_ else curve$distinct_synthesized_peptides[hit],
      reachable = if (length(at)) TRUE else if (any(missing | duplicate_peptide)) NA else FALSE,
      reason = if (length(at)) "REACHED_WITH_AVAILABLE_ORDERED_UNITS" else if (any(missing))
        "UNSCORED_UNITS_PREVENT_COMPLETE_ESTIMATE" else if (any(duplicate_peptide))
        "SHARED_PEPTIDE_HLA_JOINT_PROBABILITY_UNAVAILABLE" else "TARGET_UNREACHABLE_AVAILABLE_POOL")
  }))
  good <- !missing
  subset_prob <- if (any(duplicated(units$peptide_sequence[good]))) NA_real_ else -expm1(sum(log1p(-rates[good])))
  summary <- data.frame(sample_id = units$sample_id[1L], unit = unit, N_available = n,
    N_scored = sum(good), N_unscored = sum(missing), distinct_synthesized_peptides = length(unique(units$peptide_sequence)),
    distinct_scored_peptides = length(unique(units$peptide_sequence[good])),
    expected_positive_units = if (any(missing)) NA_real_ else sum(rates),
    probability_at_least_one = utils::tail(curve$probability_at_least_one, 1L),
    expected_positive_scored_subset = sum(rates[good]), probability_scored_subset = subset_prob,
    profile_id = p$data$id, profile_version = p$data$version, mode = p$data$mode,
    probability_interpretation = p$interpretation, endpoint = p$data$endpoint, assay = p$data$assay)
  warnings <- c("Candidate outcomes are assumed independent; shared patient/tumor biology can invalidate this assumption.",
    "Evidence A/B/C/D and B1/B2 never multiply or modify probabilities. A prior result is not a guarantee in a new patient.",
    "Endpoint/assay rates do not establish ELISpot, tumor-recognition or clinical validity for other endpoints.")
  if (p$data$mode == "scenario") warnings <- c(warnings, p$interpretation)
  if (p$data$binning$basis == "cohort_relative") warnings <- c(warnings,
    "EXPLORATORY candidate-set-dependent mid-rank bins; not frozen empirical calibration. Ties can create empty bins.")
  if (p$data$mode == "empirical") warnings <- c(warnings,
    "User-supplied matching empirical counts are not independent validation; population transportability remains unestablished.")
  if (unit == "peptide") warnings <- c(warnings,
    paste("Single-HLA", primary_hla, "proxy for each synthesized peptide; alternatives retained, no cross-HLA probability combination."))
  if (any(missing)) warnings <- c(warnings, "Missing predictions remain unavailable, never assumed negative or zero probability.")
  if (any(duplicate_peptide)) warnings <- c(warnings,
    "Multiple HLA assays sharing a peptide: expectation counts assays, not positive peptides; joint success/targets unavailable after repeated peptide.")
  structure(list(summary = summary, units = units, pair_evidence = pairs, curve = curve, targets = targets,
    composition = list(score_bins = table(factor(units$score_bin, levels = 1:4), useNA = "ifany"),
      evidence_classes = table(factor(units$evidence_category, levels = c("A", "B", "C", "D")), useNA = "ifany"),
      B_subtypes = table(factor(unlist(units$b_subtypes), levels = c("B1", "B2")))),
    profile = p, reference = attr(pairs, "reference"), evidence_policy = attr(pairs, "evidence_policy"),
    predictor_identity = input$identity, input_provenance = input$input,
    excluded_channel_rows = input$excluded_channel_rows,
    assumptions = list(version = "experimental-yield-0.1.0", independence = TRUE,
      evidence_probability_multiplier = FALSE, ordering = ordering,
      bin_reference_unit = "distinct scored sample/peptide/exact-HLA pairs",
      fixed_threshold_ties = "lower bin; intervals (-Inf,t30],(t30,t60],(t60,t90],(t90,Inf) in quality direction"),
    warnings = warnings), class = "neo_experimental_yield")
}

#' @export
plot.neo_experimental_yield <- function(x, y = NULL, ...) {
  graphics::plot(x$curve$N_tested, x$curve$probability_at_least_one, type = "l",
    ylim = c(0, 1), xlab = paste("Number of", x$summary$unit, "units tested"),
    ylab = "P(at least one positive), independence assumed", ...)
  graphics::abline(h = x$targets$target, lty = 3, col = "grey70")
  invisible(x$curve)
}
