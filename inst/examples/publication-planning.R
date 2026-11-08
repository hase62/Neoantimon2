# Optional, sourceable composition example. Not an exported package API.
# No predictors, downloads, project writes or reconstruction are performed.
# Existing neo_experimental_yield supplies ALL proxy selection, bins and maths.
publication_model_selection <- function(x, channels, N, profile,
    pool = c("common_evaluated", "model_available"), allow_different_tasks = FALSE) {
  pool <- match.arg(pool)
  if (!is.data.frame(x) || !all(c("hla_class", "prediction_scope", "candidate_id") %in% names(x)))
    stop("Supply priorities-compatible records with explicit hla_class and prediction_scope")
  if (length(unique(x$sample_id)) != 1L || length(unique(x$hla_class)) != 1L ||
      anyNA(x$hla_class)) stop("Select exactly one sample and MHC class")
  if ("prediction_scope" %in% names(x) && anyNA(x$prediction_scope)) stop("Unknown prediction scope")
  if ("prediction_scope" %in% names(x) && any(x$prediction_scope != "peptide_hla"))
    stop("This comparison requires exact peptide-HLA predictions; peptide-only scope is separate")
  if (!is.numeric(N) || !length(N) || anyNA(N) || any(!is.finite(N) | N < 0 | N != floor(N)))
    stop("N must contain nonnegative integer budgets")
  if (!is.list(channels) || length(channels) < 2L || is.null(names(channels)) ||
      any(!nzchar(names(channels))) || anyDuplicated(names(channels))) stop("Name at least two channels")
  # Validate identities across models, not just within each call to Yield.
  ids <- unique(x[c("candidate_id", "sample_id", "peptide_sequence", "mhc_key")])
  if (anyDuplicated(ids$candidate_id)) stop("candidate_id maps to conflicting identity")
  base <- lapply(channels, function(ch) {
    # The public API rejects ambiguous model/version/run/metric selectors.
    y <- Neoantimon2::neo_experimental_yield(x, profile, channel = ch)
    z <- y$input_provenance
    for (f in names(y$predictor_identity)) z <- z[!is.na(z[[f]]) & z[[f]] == y$predictor_identity[[f]], , drop=FALSE]
    for (f in intersect(c("analysis_id", "prediction_scope"), names(z)))
      if (length(unique(z[[f]])) != 1L) stop("Select one analysis and prediction scope per channel")
    y
  })
  signatures <- vapply(base, function(y) jsonlite::toJSON(y$predictor_identity, auto_unbox=TRUE), "")
  if (anyDuplicated(signatures)) stop("Repeated identical channels are not distinct comparisons")
  tasks <- vapply(base, function(y) y$predictor_identity$task, "")
  if (length(unique(tasks)) > 1L && !isTRUE(allow_different_tasks))
    stop("Different tasks require explicit allow_different_tasks; this is model-and-objective selection, not accuracy")
  evaluated <- lapply(base, function(y) {
    p <- y$pair_evidence
    sort(unique(as.character(unlist(p$all_candidate_ids[p$yield_evaluated], use.names=FALSE))), method="radix")
  })
  common <- Reduce(intersect, evaluated)
  plans <- lapply(seq_along(base), function(i) {
    if (pool == "model_available") return(base[[i]])
    if (!length(common)) return(NULL)
    # Restrict the exact common pMHC pool BEFORE proxy selection, never just top-N.
    Neoantimon2::neo_experimental_yield(x[x$candidate_id %in% common, , drop=FALSE],
      profile, channel=channels[[i]])
  })
  names(plans) <- names(channels)
  coverage <- do.call(rbind, lapply(seq_along(base), function(i) {
    y <- base[[i]]; p <- y$pair_evidence
    available <- unique(unlist(p$all_candidate_ids, use.names=FALSE))
    data.frame(model=names(channels)[i],candidate_id=ids$candidate_id,
      status=ifelse(ids$candidate_id %in% evaluated[[i]], "evaluated",
        ifelse(ids$candidate_id %in% available, "recorded_unavailable", "no_record")),
      common_evaluated=ids$candidate_id %in% common, stringsAsFactors=FALSE)
  }))
  selection <- list(); summary <- list(); sets <- list()
  for (i in seq_along(plans)) for (n in sort(unique(N))) {
    y <- plans[[i]]
    u <- if (is.null(y)) data.frame() else y$units[!is.na(y$units$raw_score) &
      y$units$prediction_status == "evaluated", , drop=FALSE]
    k <- min(n, nrow(u)); chosen <- if (k) u[seq_len(k), , drop=FALSE] else u
    key <- paste(i,n,sep=":"); sets[[key]] <- if (k) chosen$peptide_sequence else character()
    if (k) {
      chosen$model <- names(channels)[i]; chosen$requested_N <- n
      chosen$selection_rank <- seq_len(k); chosen$comparison_pool <- pool
      selection[[length(selection)+1L]] <- chosen
    }
    summary[[length(summary)+1L]] <- data.frame(model=names(channels)[i], requested_N=n,
      selected_peptides=k, available_peptides=nrow(u), common_evaluated_candidates=length(common),
      model_evaluated_candidates=length(evaluated[[i]]), declared_candidates=nrow(ids),
      selection_coverage=length(evaluated[[i]])/nrow(ids), comparison_pool=pool,
      status=if(k<n) "FEWER_THAN_REQUESTED" else "REQUEST_FILLED")
  }
  overlap <- list()
  for (i in seq_len(length(plans)-1L)) for (j in seq.int(i+1L,length(plans))) for (n in sort(unique(N))) {
    a <- sets[[paste(i,n,sep=":")]]; b <- sets[[paste(j,n,sep=":")]]
    both <- intersect(a,b); union <- union(a,b)
    row <- data.frame(model_1=names(channels)[i],model_2=names(channels)[j],requested_N=n,
      intersection=length(both),union=length(union),jaccard=if(length(union))length(both)/length(union) else NA_real_,
      comparison_pool=pool,comparison_interpretation=if(length(unique(tasks))>1L) "MODEL_AND_TASK_CHOICE_NOT_ACCURACY" else "MODEL_SELECTION_NOT_ACCURACY")
    row$common_peptides <- I(list(both));row$model_1_only <- I(list(setdiff(a,b)));row$model_2_only <- I(list(setdiff(b,a)))
    overlap[[length(overlap)+1L]] <- row
  }
  list(summary=do.call(rbind,summary),selected=if(length(selection)) do.call(rbind,selection) else data.frame(),
    pairwise=do.call(rbind,overlap),coverage=coverage,common_candidate_ids=common,
    predictor_identities=lapply(base,`[[`,"predictor_identity"),full_pool_plans=base,plans=plans,
    input_records=x,assumptions=c("No cross-model score combination; no accuracy inference",
      "Existing best-HLA proxy and unique-peptide order; all alternatives/source_rows retained",
      "Profile rates do not affect selection; probabilities in selection plans are not used to compare models",
      "Common pool is exact evaluated candidate intersection; model-available pools may differ"))
}

publication_scenarios <- function(x, profiles, budgets, channel=NULL, reference=NULL) {
  if (!is.list(profiles) || !length(profiles) || is.null(names(profiles)) ||
      any(!nzchar(names(profiles))) || anyDuplicated(names(profiles))) stop("Name explicit profiles")
  if (!is.numeric(budgets) || !length(budgets) || anyNA(budgets) ||
      any(!is.finite(budgets) | budgets < 0 | budgets != floor(budgets))) stop("Invalid budgets")
  plans <- lapply(profiles, function(p) Neoantimon2::neo_experimental_yield(x,p,
    channel=channel,reference=reference,unit="peptide",primary_hla="best_score",ordering="score"))
  # Call on the FULL identical pool for every scenario. Never re-bin a top-N subset.
  first <- plans[[1]]$units
  for (y in plans) if (!identical(y$profile$data$binning, plans[[1]]$profile$data$binning) || !identical(y$units[c("peptide_sequence","mhc_key","score_bin")],
      first[c("peptide_sequence","mhc_key","score_bin")]))
    stop("Scenario comparison requires identical peptide order, primary HLA and binning")
  annotate <- function(z,i) {
    z$scenario <- names(plans)[i];z$profile_id <- plans[[i]]$profile$data$id
    z$profile_version <- plans[[i]]$profile$data$version
    z$probability_interpretation <- plans[[i]]$summary$probability_interpretation;z
  }
  curves <- targets <- summaries <- list()
  for(i in seq_along(plans)) {
    y <- plans[[i]];curves[[i]]<-annotate(y$curve,i);targets[[i]]<-annotate(y$targets,i)
    for(n in sort(unique(budgets))) {
      k <- min(n,y$summary$N_available)
      z <- y$curve[y$curve$N_tested==k,,drop=FALSE]
      z$requested_N<-n;z$available_pool<-y$summary$N_available
      z$budget_status<-if(n>k) "FEWER_THAN_REQUESTED" else "AVAILABLE"
      z$assumed_rates_json<-jsonlite::toJSON(y$profile$rates,auto_unbox=TRUE)
      summaries[[length(summaries)+1L]]<-annotate(z,i)
    }
  }
  list(budgets=do.call(rbind,summaries),targets=do.call(rbind,targets),
    curves=do.call(rbind,curves),plans=plans,
    assumptions=c("User-supplied sensitivity assumptions, not confidence intervals",
      "Expected positives differ from independence-based P(at least one)",
      "No evidence-category probability multiplier; no extrapolated candidates"))
}
