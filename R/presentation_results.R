presentation_candidate_view <- function(x, candidate_id) {
  empty <- is.null(x)
  links <- if (empty) list() else Filter(function(l) l$candidate_id == candidate_id, x$processing_links)
  applicable_links <- Filter(function(l) l$applicability == "applicable", links)
  query_ids <- vapply(applicable_links, `[[`, "", "processing_query_id")
  e <- if (empty) list() else Filter(function(e) e$processing_query_id %in% query_ids, x$processing_evidence)
  queries <- if (empty) list() else Filter(function(q) q$processing_query_id %in%
    vapply(links, `[[`, "", "processing_query_id"), x$queries)
  processing <- Filter(function(e) e$task != "pmhc_stability", e)
  stability <- Filter(function(e) e$task == "pmhc_stability", e)
  ms <- if (empty) list() else Filter(function(l) l$candidate_id == candidate_id, x$ms_links)
  direct <- Filter(function(l) l$status == "exact_unmodified_ligand_match", ms)
  observations <- if (empty) list() else Filter(function(o)
    o$immunopeptidomics_observation_id %in% vapply(ms, `[[`, "", "observation_id"), x$observations)
  psms <- if (empty) list() else Filter(function(p) p$psm_id %in% unlist(lapply(ms, `[[`, "psm_ids")), x$psms)
  apm_links <- if (empty) list() else Filter(function(l) l$candidate_id == candidate_id, x$apm_links)
  apm <- if (empty) list() else Filter(function(e)
    e$apm_observation_id %in% vapply(apm_links, `[[`, "", "apm_observation_id"), x$apm_observations)
  applicable <- Filter(function(l) l$applicability == "applicable", apm_links)
  expression <- Filter(function(l) l$measurement_type %in% c("gene_expression", "allele_expression") &&
    grepl("^HLA-", l$gene), applicable)
  states <- unique(vapply(applicable, `[[`, "", "state"))
  summary <- list(processing_evidence_count = length(processing),
    processing_origin_count = length(unique(unlist(lapply(Filter(function(l) l$processing_query_id %in%
      vapply(processing, `[[`, "", "processing_query_id"), applicable_links), function(l)
        vapply(l$origin_contexts, `[[`, "", "origin_id"))))),
    pmhc_stability_evidence_count = length(stability),
    immunopeptidomics_status = if (length(direct)) "exact_unmodified_ligand_observed" else
      if (length(ms)) "other_evidence_present" else "no_evidence_supplied",
    exact_observed_ligand_count = length(direct), ms_psm_count = length(unique(unlist(lapply(direct, `[[`, "psm_ids")))),
    ms_run_count = length(unique(vapply(Filter(function(o) o$immunopeptidomics_observation_id %in%
      vapply(direct, `[[`, "", "observation_id"), observations), `[[`, "", "ms_run"))),
    direct_hla_assignment_status = if (any(vapply(direct, function(l) l$direct_hla_assignment == "direct_allele_resolved", TRUE)))
      "direct_allele_resolved" else "hla_unresolved",
    presenter_expression_status = if (length(expression)) "evidence_present" else "not_evaluated",
    apm_evidence_status = if (!length(applicable)) "not_evaluated" else if ("direct_loss_supported" %in% states)
      "direct_loss_supported" else if ("provider_loss_call" %in% states) "provider_loss_call" else
        if ("abnormality_observed" %in% states) "abnormality_observed" else "evidence_present_unresolved")
  list(summary = summary, processing_links = links, processing_queries = queries, antigen_processing = processing, pmhc_stability = stability,
    immunopeptidomics = list(links = ms, observations = observations, psms = psms),
    hla_expression = Filter(function(e) e$measurement_type %in% c("gene_expression", "allele_expression") &&
      grepl("^HLA-", e$gene), apm), antigen_presentation_machinery = list(observations = apm, links = apm_links),
    processing_statuses = lapply(links, function(l) {
      status <- Filter(function(s) s$processing_query_id == l$processing_query_id, x$statuses)
      list(candidate_id = candidate_id, processing_query_id = l$processing_query_id,
        status = if (l$applicability != "applicable") l$status else status[[1L]]$status)
    }),
    interpretation = "Independent evidence only; no composite score, pathway-intact claim, or grade/rank change.")
}
presentation_results_append <- function(out, project, snapshot) {
  x <- presentation_read(project, snapshot)
  if (is.null(x)) return(out)
  views <- lapply(out$candidate_id, function(id) presentation_candidate_view(x, id)$summary)
  template <- presentation_candidate_view(x, "")$summary
  for (field in names(template)) out[[field]] <- vapply(views, `[[`, template[[field]], field)
  out
}
presentation_explain_append <- function(candidates, project, snapshot) {
  x <- presentation_read(project, snapshot)
  if (is.null(x)) return(candidates)
  lapply(candidates, function(c) { c$presentation_evidence <- presentation_candidate_view(x, c$result$candidate_id); c })
}
presentation_explain_text <- function(item) {
  x <- item$presentation_evidence
  if (is.null(x)) return("")
  render <- function(e) if (!length(e)) "no evidence supplied" else paste(vapply(e, function(r)
    paste(r$provider, r$metric, r$raw_value, r$unit, "query", r$processing_query_id), ""), collapse = "; ")
  paste0("\nANTIGEN PROCESSING: ", render(x$antigen_processing),
    "\nPROCESSING ORIGIN LINKS: ", paste(unique(unlist(lapply(x$processing_links, function(l)
      vapply(l$origin_contexts, `[[`, "", "origin_id")))), collapse = "; "),
    "\nPMHC BINDING / PRESENTATION: ", item$result$prediction_status,
    "\nPMHC STABILITY: ", render(x$pmhc_stability),
    "\nIMMUNOPEPTIDOMICS: ", x$summary$exact_observed_ligand_count, " exact unmodified ligand observations; ",
    x$summary$ms_psm_count, " PSMs; HLA assignment=", x$summary$direct_hla_assignment_status,
    "; other/modified/ambiguous matches=", length(x$immunopeptidomics$links) - x$summary$exact_observed_ligand_count,
    "\nHLA AVAILABILITY / EXPRESSION: ", x$summary$presenter_expression_status,
    "; gene and allele/heterodimer-chain measurements retained separately",
    "\nANTIGEN-PRESENTATION MACHINERY: ", x$summary$apm_evidence_status,
    "\nHLA-SPECIFIC LOH: ", if (is.null(item$hla_loh)) "not_evaluated" else item$hla_loh$candidate_link$presenter_hla_loh_status)
}
