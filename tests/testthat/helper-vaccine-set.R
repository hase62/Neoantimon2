vaccine_fixture <- function(file = "evidence.json") jsonlite::read_json(testthat::test_path("..", "fixtures", "vaccine_set", file))
vaccine_binding <- function(...) Neoantimon2:::vaccine_policy_bind(c(list(version = "vaccine-set-policy-1.0.0", max_peptides = 4L), list(...)))
vaccine_records <- function(fixture = vaccine_fixture()) {
  records <- list()
  for (spec in fixture$candidates) {
    presenter <- Neoantimon2:::peptide_presenters(list(spec$presenter))[[1L]]
    key <- Neoantimon2:::biological_id("candidate_key", list(species = "NCBITaxon:9606", hla_class = presenter$hla_class,
      ligand_sequence = spec$ligand, mhc_key = presenter$mhc_key), Neoantimon2:::peptide_contract())
    cells <- lapply(seq_along(fixture$contexts), function(i) {
      ct <- fixture$contexts[[i]]; state <- spec$states[[i]]
      assay <- if (state == "not_observed_evaluable") longitudinal_assay(ct$sample_id, key) else NULL
      clonal <- if (isTRUE(spec$clonal) && state == "observed") list(list(clonal_evidence_id = paste("clonal", ct$sample_id, spec$event),
        mapping_status = "mapped", sample_id = ct$sample_id, event_ids = list(spec$event), provider_cluster_id = "1",
        provider_mutation_id = spec$event, cluster_scope = list(provider = "external-PyClone", run_id = "fit:1", subject_id = fixture$subject_id),
        cellular_prevalence = .72, cellular_prevalence_std = .08, cluster_assignment_prob = .65)) else list()
      list(sample_id = ct$sample_id, state = state, evaluability = assay, observations = list(), clonal_evidence = clonal,
        domain_evaluability = list(list(domain = "ms", cell = list(state = if (i == 1L) "observed" else "not_evaluated"))))
    })
    t <- list(subject_id = fixture$subject_id, feature_key = key, ordering = Neoantimon2:::longitudinal_order(fixture$contexts)$ordering,
      cells = cells, summary = Neoantimon2:::longitudinal_summary(cells, fixture$contexts, fixture$relations),
      trajectory_evidence_id = paste0("authoritative-fixture:", spec$label))
    for (i in seq_along(cells)) if (cells[[i]]$state == "observed") {
      ct <- fixture$contexts[[i]]; sample <- ct$sample_id
      row <- list(candidate_id = Neoantimon2:::biological_id("candidate_id", list(context_id = "somatic:patient-1", sample_id = sample, candidate_key = key),
        Neoantimon2:::peptide_contract()), candidate_key = key, sample_id = sample, ligand_sequence = spec$ligand,
        peptide_sequence = spec$ligand, mhc_key = presenter$mhc_key, hla_allele = presenter$hla_allele, hla_class = presenter$hla_class,
        final_grade = if (is.null(spec$grade)) "B" else spec$grade, grading_status = if (identical(spec$grade, "NOT_GRADED")) "NOT_GRADED" else "GRADED",
        technical_eligibility = if (identical(spec$grade, "INVALID")) "invalid" else "eligible", policy_version = if (is.null(spec$source)) "somatic-priority-v1" else "source-not-graded-v1",
        rank_within_sample_class = if (is.null(spec$source)) 2L else NULL, priority_origin_id = paste0("origin:", spec$label, ":", sample),
        epitope_key = spec$epitope, self_status = if (is.null(spec$exact_self)) "NOT_EVALUABLE" else "EVALUATED",
        self_exact_match = spec$exact_self, ccf_value = .31)
      origin <- list(origin_id = row$priority_origin_id, source_event_ids = list(spec$event), transcript_id = "TX:fixture",
        peptide_start_aa = spec$start, peptide_end_aa = spec$start + nchar(spec$ligand) - 1L,
        origin_validity = "VALID", source_type = spec$source)
      ms <- list(summary = list(immunopeptidomics_status = if (i == 1L) "exact_unmodified_ligand_observed" else "no_evidence_supplied"),
        immunopeptidomics = list(observations = if (i == 1L) list(list(context = list(sample_id = sample), ms_run = "ms:S1")) else list()))
      apm <- list()
      if (sample %in% unlist(spec$apm_samples)) apm <- list(list(apm_observation_id = paste0("apm:", sample), gene = "B2M", role = "class_i_structural",
        hla_class = "I", state = "direct_loss_supported", context = list(sample_id = sample), interpretation_reference = "fixture:validated-loss"))
      ms$antigen_presentation_machinery <- list(observations = apm, links = lapply(apm, function(a) list(apm_observation_id = a$apm_observation_id,
        applicability = if (presenter$hla_class == "I") "applicable" else "not_applicable", state = a$state)))
      loh <- if (sample %in% unlist(spec$loh_samples)) list(presenter_allele = presenter$mhc_key, candidate_link = list(
        presenter_hla_loh_call = "loss_supported", presenter_hla_loh_evidence_id = paste0("loh:", sample, ":", presenter$mhc_key))) else NULL
      e <- list(result = row, origins = list(origin), hla_loh = loh, presentation_evidence = ms,
        somatic_evidence = list(event_evidence = list(list(event_id = spec$event, ccf_value = .31))),
        receptor_evidence = list(status = "no_observation"))
      records <- c(records, list(list(result = row, candidate_key = key, subject_id = fixture$subject_id, trajectory = t, evidence = e)))
    }
  }
  records
}
vaccine_fixture_projection <- function(binding = vaccine_binding(), fixture = vaccine_fixture())
  Neoantimon2:::vaccine_project_records(vaccine_records(fixture), fixture$contexts, fixture$relations, binding$normalized)
vaccine_fixture_select <- function(binding = vaccine_binding(), fixture = vaccine_fixture())
  Neoantimon2:::vaccine_construct(vaccine_fixture_projection(binding, fixture), binding)
vaccine_keys <- function(x) sort(vapply(x$selections[[1L]]$selected_candidates, `[[`, "", "candidate_key"))
vaccine_annotate <- function(p, binding) {
  p$config$vaccine_set <- binding
  neo_run(p)
}
vaccine_real_e2e <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    original <- longitudinal_e2e()
    runs <- original$runs
    # Frozen integrity tests deliberately tamper with their own input handles.
    # This fixture owns a separate artifact copied from immutable imported data.
    clonal_path <- tempfile(fileext = ".tsv")
    clonal_fields <- c("mutation_id", "sample_id", "cluster_id", "cellular_prevalence", "cellular_prevalence_std", "cluster_assignment_prob")
    writeLines(c(paste(clonal_fields, collapse = "\t"), vapply(original$overlay$clonal_evidence, function(e)
      paste(vapply(clonal_fields, function(f) if (is.null(e$raw_fields[[f]])) "." else e$raw_fields[[f]], ""), collapse = "\t"), "")), clonal_path)
    extra <- lapply(original$extra, function(i) { if (i$input_role == "pyclone_vi") i$path <- clonal_path; i })
    # Genuine source imports use the same synthetic primary reference as the small variants.
    for (i in seq_along(runs)) runs[[i]] <- neo_run(recognition_clone_with(runs[[i]], list(list(
      input_id = "vaccine-fusion", input_role = "star_fusion", sample_id = paste0("S", i),
      path = normalizePath(testthat::test_path("..", "fixtures", "source", "star-fusion.tsv")), version = "synthetic-v1"))))
    settings <- original$settings
    settings$sample_runs <- lapply(seq_along(runs), function(i) list(sample_id = paste0("S", i), path = runs[[i]]$path))
    before <- neo_run(longitudinal_project(settings, extra))
    pool <- neo_results(before)
    # This declared allowlist reserves A*24:02 coverage for source candidates,
    # demonstrating NOT_GRADED participation without any grade reinterpretation.
    allowed <- pool$candidate_key[pool$final_grade == "NOT_GRADED" | pool$hla_class == "II" | pool$mhc_key == "HLA-A*02:01"]
    binding <- vaccine_binding(min_class_i = 2L, min_class_ii = 1L,
      allowlist = as.list(unique(allowed)),
      required_presenters = list("HLA-A*02:01", "HLA-A*24:02", "HLA-DRA*01:01~HLA-DRB1*04:01"))
    after <- vaccine_annotate(before, binding)
    cache <<- list(before = before, after = after, selection = Neoantimon2:::vaccine_read(after), runs = runs)
    cache
  }
})

vaccine_review_duplicate_policy <- function() {
  raw <- vaccine_fixture("review-support.json")$duplicate_policy
  key <- vaccine_records()[[1L]]$candidate_key
  raw$allowlist <- list(key, "unused:b", key)
  raw
}
