ia_simple_allo <- function(x) {
  allo_keys(x, c("donor", "recipient"), label = "donor/recipient typing")
  loci <- unlist(indirect_loci()); typing <- list()
  for (role in c("donor", "recipient")) {
    own <- x[[role]]; required <- c(loci, if (role == "recipient") "DRA")
    if (!is.list(own) || is.null(names(own)) || anyDuplicated(names(own)) ||
        !setequal(names(own), required)) stop("Supply explicit donor/recipient loci A/B/C/DRB1/DQB1 and recipient DRA")
    for (locus in required) {
      alleles <- own[[locus]]
      if (!is.character(alleles) || length(alleles) != 2L || anyNA(alleles))
        stop("Supply two explicit allele copies per locus (repeat for declared homozygosity); no genotype imputation")
      alleles <- vapply(alleles, allo_allele, "", locus = locus)
      for (a in sort(unique(alleles), method = "radix")) typing[[length(typing) + 1L]] <- list(
        subject_id = role, typing_scenario_id = "typed:1", locus = locus, alleles = list(a),
        complete = TRUE, copy_number = as.integer(sum(alleles == a)), phase = NULL)
    }
  }
  list(version = allo_version(), subjects = list(list(subject_id = "donor"), list(subject_id = "recipient")),
    contexts = list(list(transplant_donor_subject_id = "donor", transplant_recipient_subject_id = "recipient",
      source_subject_id = "donor", self_subject_id = "recipient", presenter_subject_id = "recipient", responding_subject_id = "recipient",
      source_subject_role = "donor", self_subject_role = "recipient", presenter_subject_role = "recipient", responding_subject_role = "recipient",
      transplant_type = "SOT", scenario_id = "sot_recipient_response", direction = "donor_to_recipient",
      tissue = "unspecified_solid_organ", apc_compartment = "recipient_APC", occasion_id = "user_unspecified",
      claim_level = "hypothetical_molecular_scenario")), typing = typing,
    self_search = list(list(scenario_id = "sot_recipient_response", self_subject_id = "recipient", typing_scenario_id = "typed:1",
      scope = "all_self_hla", required_loci = as.list(loci), required_domains = list("full_precursor"),
      catalog_complete = TRUE, isoform_policy = "explicit_mapped_isoforms", sequence_policy = "precursor")))
}
ia_plan <- function(foundation, presenters) {
  if (length(foundation$contexts) != 1L) stop("Autorun v0.1 requires one explicit donor-to-recipient SOT context")
  context <- foundation$contexts[[1L]]
  if (context$semantics$scenario_id != "sot_recipient_response") stop("Only SOT donor-to-recipient is supported")
  generation <- indirect_generate(foundation, context)
  catalog <- indirect_self_catalog(foundation, context, NULL)
  if (is.character(presenters)) presenters <- lapply(presenters, function(p) list(
    scenario_id = "sot_recipient_response", subject_id = context$semantics$presenter_subject_id,
    presenter = p, resolution_policy = "explicit_typed_alpha_beta"))
  presenters <- allo_records(lapply(presenters, function(p) {
    if (p$scenario_id != "sot_recipient_response" || p$subject_id != context$semantics$presenter_subject_id ||
        p$resolution_policy != "explicit_typed_alpha_beta") stop("Explicit recipient presenter binding required")
    p$presenter <- peptide_presenters(list(p$presenter))[[1L]]$mhc_key
    ia_vendor_allele(p$presenter); p
  }))
  if (!length(presenters)) stop("Explicit recipient DRA~DRB1 presenters are required")
  # No substitution through reference protein equivalence is accepted by autorun.
  if (any(vapply(foundation$typing, function(t) t$evaluability$status == "evaluated" && is.null(t$resolved_allele), TRUE)))
    stop("Exact pinned allele mapping required; no allele completion or sequence-equivalent substitution")
  scope <- indirect_presenters(list(presenters = presenters), foundation, context)
  if (any(vapply(scope, function(p) p$status != "evaluated", TRUE)))
    stop("Recipient presenter is unsupported or not exactly represented in supplied typing/reference")
  evidence <- list(); add <- function(peptide, presenter, role, extra) {
    id <- indirect_id("indirect-autorun-query-v1", list(peptide_sequence = peptide, presenter = presenter))
    evidence[[length(evidence) + 1L]] <<- c(list(query_id = id, peptide_sequence = peptide, presenter = presenter, role = role), extra)
  }
  for (p in presenters) {
    for (o in generation$occurrences) if (o$status == "evaluated") add(o$peptide_sequence, p$presenter, "donor", list(
      occurrence_id = o$indirect_peptide_occurrence_id, source_occurrence_id = o$semantics$source_occurrence_id,
      allele = o$source_allele_hypothesis, sequence_id = o$hla_source_sequence_id,
      start = o$semantics$start, end = o$semantics$end))
    for (w in catalog$windows) add(w$peptide_sequence, p$presenter, "self", list(
      allele = w$allele, sequence_id = w$sequence_id, start = w$start, end = w$end,
      self_universe_id = catalog$universe$self_universe_id))
  }
  evidence <- allo_records(evidence)
  queries <- allo_records(lapply(evidence, function(x) x[c("query_id", "peptide_sequence", "presenter")]))
  list(version = "indirect-autorun-plan-0.2.0", profile = indirect_profile(), presenters = presenters,
    queries = queries, query_links = evidence, generation = generation,
    self_catalog_complete = catalog$complete, self_window_count = length(catalog$windows),
    self_reason_codes = catalog$reason_codes, self_catalog = catalog)
}

#' Run indirect-HLA Class II analysis with a pinned predictor profile
#'
#' Reuses the frozen SOT donor-to-recipient molecular baseline; not official
#' PIRCHE reproduction, Frost, or a clinical risk score. No vendor acquisition.
#' @param allo Existing allo-foundation config, or donor/recipient named typing
#'   lists with two explicit copies at A/B/C/DRB1/DQB1 and recipient DRA.
#' @param presenters Explicit recipient DRA*01:01~DRB1 pairs (two-field DRB1).
#' @param reference List with path and SHA-256 of a versioned literal precursor
#'   hla-reference-1.0.0 manifest. No nearest-allele resolution is performed.
#' @param executable Authorized local NetMHCIIpan 4.3k (standard) or 3.2 command; never downloaded.
#' @param model List with full installation path containing data/, identity and
#'   optional expected tree SHA-256. 4.3k additionally permits runtime_subdir and pinned internal data links; 3.2 keeps physical-tree validation.
#' @param output_dir New output directory. Existing files are never overwritten.
#' @param dry_run Plan donor and self queries without executing any command.
#' @param resume Reopen a completed run/plan only after exact request/tool and
#'   artifact hash verification. Failed/interrupted runs are retained, not resumed.
#' @param authorization Explicit list(status, terms). status must be
#'   authorized_local_installation, or mock_only for controlled software tests.
#' @param timeout Per-command timeout in seconds.
#' @param batch_size Maximum unique peptides per presenter invocation.
#' @param predictor_version Exact 4.3k (default), explicit legacy 3.2, or opt-in auto detection.
#' @param inversion Vendor default (DP inversion only), or all loci for 4.3k.
#' @return A list with status, manifest, plan, output_dir and, when complete,
#'   existing public result and explanation projections. Missing evidence never
#'   becomes a definitive zero. Raw output, imports and checksums are persisted.
#' @export
neo_indirect_autorun <- function(allo, presenters, reference, executable = NULL,
                                model = NULL, output_dir, dry_run = FALSE, resume = FALSE,
                                authorization = NULL, timeout = 600L, batch_size = 256L,
                                predictor_version = c("4.3k", "3.2", "auto"), inversion = c("vendor_default", "all")) {
  predictor_version <- match.arg(predictor_version); inversion <- match.arg(inversion)
  if (predictor_version == "3.2" && inversion != "vendor_default") stop("3.2 cannot use 4.3k inversion settings")
  selected_profile <- if (predictor_version == "3.2") indirect_profile() else if (predictor_version == "4.3k") ia43_profile() else "auto-resolve-pinned-NetMHCIIpan-profile"
  setup_flag(dry_run, "dry_run"); setup_flag(resume, "resume")
  for (n in list(timeout = timeout, batch_size = batch_size))
    if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1 || n != floor(n)) stop("timeout/batch_size must be positive integers")
  allo_keys(reference, c("path", "sha256"), label = "pinned reference")
  if (!is.character(reference$sha256) || !grepl("^[0-9a-f]{64}$", reference$sha256)) stop("Explicit reference SHA-256 required")
  if (is.null(allo$version)) allo <- ia_simple_allo(allo)
  if (!is.null(allo$reference$sha256) && allo$reference$sha256 != reference$sha256) stop("Conflicting reference bindings")
  allo$reference <- reference; bound <- allo_config(allo, getwd()); foundation <- allo_foundation(bound)
  plan <- ia_plan(foundation, presenters)
  plan$profile <- selected_profile; plan$requested_predictor_version <- predictor_version
  tool <- NULL; pinned <- NULL; mock <- FALSE
  if (!dry_run) {
    allo_keys(authorization, c("status", "terms"), label = "explicit local authorization")
    if (!authorization$status %in% c("authorized_local_installation", "mock_only")) stop("Explicit local authorization status required")
    scalar_text(authorization$terms, "local license/authorization terms")
    mock <- authorization$status == "mock_only"
    exe <- setup_executable(executable, "netMHCIIpan")
    if (is.null(exe)) stop("NetMHCIIpan unavailable: install a licensed local executable and model data; use dry_run=TRUE to plan")
    tool <- list(path = exe, sha256 = ia_hash(exe), authorization = authorization)
    pinned <- if (predictor_version == "3.2") ia_model(model) else ia43_model(model)
  }
  request <- list(version = "indirect-autorun-0.2.0", code_sha256 = ia_code_hash(), profile = selected_profile,
    predictor_version = predictor_version, inversion = inversion,
    foundation_id = foundation$allo_analysis_id, reference_sha256 = reference$sha256,
    presenters = plan$presenters, query_ids = lapply(plan$queries, `[[`, "query_id"), dry_run = dry_run,
    executable = tool, model = if (is.null(pinned)) NULL else pinned[intersect(c("path", "rdir", "identity", "tree_sha256"), names(pinned))],
    batch_size = as.integer(batch_size), timeout = as.integer(timeout))
  request_id <- infra_hash(request); root <- setup_dest(output_dir)
  if (dir.exists(root)) {
    if (!resume) stop("Output directory exists; use resume=TRUE only for an identical completed run/plan")
    f <- file.path(root, "run.json")
    if (!file.exists(f)) stop("Existing directory has no completed autorun manifest; preserved without changes")
    m <- jsonlite::read_json(f)
    if (!identical(m$request_id, request_id)) stop("Resume request/tool/reference identity mismatch")
    if (!m$status %in% c("COMPLETED", "DRY_RUN")) stop("Failed/interrupted run is retained; resume supports completed runs/plans only; choose a new output_dir")
    ia_verify_files(root, m); return(ia_return(root))
  }
  if (resume) stop("No existing run to resume; choose resume=FALSE for a fresh directory")
  if (!dir.create(root, recursive = TRUE, showWarnings = FALSE)) stop("Cannot create output directory")
  lock <- file.path(root, ".autorun-lock")
  if (!dir.create(lock, showWarnings = FALSE)) stop("Autorun output is locked")
  writeLines(as.character(Sys.getpid()), file.path(lock, "owner"))
  on.exit(unlink(lock, recursive = TRUE), add = TRUE) # Our own transient lock only.
  manifest <- list(version = "indirect-autorun-run-0.2.0", request_id = request_id, request = request,
    status = "RUNNING", real_predictor_executed = FALSE, execution_kind = if (dry_run) "NOT_EXECUTED" else if (mock) "MOCK_EXECUTION" else "AUTHORIZED_LOCAL_PROVIDER",
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), claim = "Literature-derived molecular baseline; not official PIRCHE or clinical risk")
  ia_json(manifest, file.path(root, "run.json"))
  tryCatch({
    file.copy(bound$reference$path, file.path(root, "reference.json"))
    if (ia_hash(file.path(root, "reference.json")) != reference$sha256) stop("Reference changed during setup")
    allo$reference <- list(path = "reference.json", sha256 = reference$sha256)
    ia_json(allo, file.path(root, "allo-config.json")); ia_json(plan, file.path(root, "plan.json"))
    ia_json(foundation, file.path(root, "foundation.json"))
    ia_table(plan$queries, file.path(root, "queries.tsv")); ia_table(plan$query_links, file.path(root, "query-origin-links.tsv"))
    ia_table(plan$generation$occurrences, file.path(root, "donor-occurrences.tsv"))
    ia_table(plan$generation$ledger, file.path(root, "generation-coverage.tsv"))
    ia_table(plan$self_catalog$windows, file.path(root, "self-windows.tsv"))
    ia_table(lapply(plan$queries, function(q) c(q, list(status = "not_evaluated"))), file.path(root, "query-coverage.tsv"))
    if (!dry_run) {
      ia_json(tool, file.path(root, "execution/executable.json")); ia_json(pinned, file.path(root, "execution/model-tree.json"))
      pp <- vapply(plan$presenters, `[[`, "", "presenter")
      probe <- if (predictor_version == "3.2") ia_probe(tool$path, pinned, pp, root, timeout, mock) else
        ia43_probe(tool$path, pinned, pp, root, timeout, mock, predictor_version, inversion)
      modern <- identical(probe$version, "4.3k")
      plan$profile <- if (modern) ia43_profile() else indirect_profile()
      ia_json(plan, file.path(root, "plan.json"))
      probe$presenter_map <- stats::setNames(vapply(pp, ia_vendor_allele, ""), pp)
      rows <- list(); commands <- probe$commands; raw_index <- list(); batch <- 0L
      for (p in pp) {
        peptides <- sort(vapply(Filter(function(q) q$presenter == p, plan$queries), `[[`, "", "peptide_sequence"), method = "radix")
        groups <- if (!length(peptides)) list() else split(peptides, ceiling(seq_along(peptides) / batch_size))
        for (g in groups) {
          batch <- batch + 1L; stem <- sprintf("raw/batch-%05d", batch)
          input <- file.path(root, paste0(stem, ".pep")); dir.create(dirname(input), recursive = TRUE, showWarnings = FALSE); writeLines(g, input)
          args <- c("-rdir", if (is.null(pinned$rdir)) pinned$path else pinned$rdir, "-inptype", "1", "-f", input, "-a", probe$presenter_map[[p]])
          if (modern) args <- c(args, "-BA", if (inversion == "all") "-inv_all" else character())
          # 4.3k keeps its official TMPDIR template: long output-derived -tdir paths
          # overflow the native nnmultipred command buffer on Darwin arm64.
          if (!modern && probe$temp_flag) { dir.create(file.path(root, "temporary"), showWarnings = FALSE); args <- c(args, "-tdir", file.path(root, "temporary")) }
          ran <- ia_exec(tool$path, args, root, stem, timeout); commands <- c(commands, list(ran$record))
          parsed <- if (modern) ia43_parse(ran$output, g, p) else ia_parse(ran$output, g, p)
          parsed <- lapply(parsed, function(r) { r$raw_output$artifact <- paste0(stem, ".stdout"); r })
          rows <- c(rows, parsed); raw_index[[length(raw_index) + 1L]] <- list(path = paste0(stem, ".stdout"), sha256 = ia_hash(file.path(root, paste0(stem, ".stdout"))))
        }
      }
      if (!identical(ia_hash(tool$path), tool$sha256) || !identical((if (predictor_version == "3.2") ia_model(model) else ia43_model(model))$tree_sha256, pinned$tree_sha256)) stop("Executable/model changed during execution")
      ia_json(raw_index, file.path(root, "execution/raw-output-index.json")); ia_json(commands, file.path(root, "execution/commands.json"))
      imported <- (if (modern) ia43_import else ia_import)(rows, file.path(root, "queries.tsv"), file.path(root, "execution/raw-output-index.json"), tool, pinned, probe, commands, mock)
      ia_json(imported, file.path(root, "binding.json"))
      config <- list(version = if (modern) "indirect-hla-netmhciipan43k-ba-1.0.0" else "indirect-hla-baseline-1.0.0", profile = plan$profile, presenters = plan$presenters,
        provider = list(path = "binding.json", sha256 = ia_hash(file.path(root, "binding.json"))))
      ia_json(config, file.path(root, "indirect-config.json"))
      existing <- indirect_baseline(foundation, indirect_config(config, root))
      # Only exported physical bindings become relative; semantic IDs/measurements are untouched.
      foundation$reference$provenance$artifact$path <- "reference.json"
      existing$provenance$artifact$path <- "binding.json"
      ia_json(indirect_public_result(existing, foundation), file.path(root, "result.json"))
      explanation <- indirect_public_explanation(existing, foundation)
      ia_json(explanation, file.path(root, "explanation.json"))
      writeLines(indirect_public_text(explanation), file.path(root, "result.txt"))
      ia_table(existing$provider$rows, file.path(root, "query-coverage.tsv"))
      ia_table(existing$summaries, file.path(root, "summary.tsv"))
      manifest$real_predictor_executed <- !mock && length(plan$queries) > 0L
      manifest$provider_version <- probe$version
      manifest$resolved_profile <- existing$profile; manifest$threshold_policy <- existing$policy$threshold
      manifest$counting_policy <- existing$policy$counting
      if (modern) {
        manifest$execution_settings <- list(BA = TRUE, context_encoding = FALSE, inversion = inversion,
          temporary_policy = "vendor TMPDIR default; resolved template retained in raw stdout",
          decision_metric = "Affinity(nM)", retained_nondecision_metrics = list("Score_EL", "%Rank_EL", "Core_Rel"))
        ia_table(ia43_metrics_table(rows), file.path(root, "model-metrics.tsv"))
      }
    }
    writeLines(c("Neoantimon2 indirect autorun v0.2; not official PIRCHE or clinical risk.",
      "Inputs: allo-config.json + pinned reference.json; queries.tsv includes donor AND recipient self queries.",
      "binding.json uses the explicit version-specific BA import contract; EL is retained separately, never used as nM.",
      "See execution/*.command.json for exact argv/exit/stdout and execution/model-tree.json for data hashes.",
      "Reopen with neo_indirect_autorun using identical inputs/tool and resume=TRUE; changed inputs/artifacts fail closed.",
      "Failed/interrupted runs are preserved and cannot resume in place; use a new output directory.",
      if (mock) "MOCK EXECUTION ONLY; REAL_PREDICTOR_E2E_NOT_TESTED" else if (dry_run) "DRY RUN: no executable was invoked" else "Authorized local provider output; no clinical validation claimed."), file.path(root, "README.txt"))
    manifest$status <- if (dry_run) "DRY_RUN" else "COMPLETED"
    ia_finish(root, manifest); ia_return(root)
  }, error = function(e) {
    manifest$status <- "FAILED"; manifest$error <- conditionMessage(e)
    ia_finish(root, manifest); stop(e)
  })
}
