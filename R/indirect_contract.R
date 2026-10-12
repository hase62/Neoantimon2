# T046-T048: local imported evidence, literature-derived molecular baseline.
indirect_profile <- function() "neo-literature-derived-indirect-hla-ii-dr-strict-v1"
indirect_loci <- function() list("A", "B", "C", "DQB1", "DRB1")
indirect_id <- function(namespace, payload) {
  paste0(namespace, ":", infra_hash(list(namespace = namespace, payload = payload)))
}
indirect_occurrence_payload <- function(source_occurrence_id, start, end, peptide_key, generation_policy_id) {
  list(source_occurrence_id = source_occurrence_id, start = as.integer(start),
    end = as.integer(end), peptide_key = peptide_key, generation_policy_id = generation_policy_id)
}
indirect_pmhc_payload <- function(context_id, candidate_key, compartment, scenario) {
  list(allo_context_id = context_id, candidate_key = candidate_key,
    presenting_compartment = compartment, presenting_scenario = scenario)
}
indirect_comparison_payload <- function(query_identity, query_type, universe_id, algorithm, self_policy) {
  list(query_identity = query_identity, query_type = query_type, self_universe_id = universe_id,
    algorithm = algorithm, self_presentation_policy = self_policy)
}
indirect_generation_policy <- function() list(profile = indirect_profile(), version = "1.0.0",
  loci = indirect_loci(), sequence_coordinates = "precursor", coverage = "full_protein",
  length = 15L, stride = 1L, padding = FALSE, completion = "none", imputation = "none")
indirect_threshold <- function() list(version = "strict-affinity-1000nm-1.0.0",
  profile = indirect_profile(), dimension = "binding_affinity", comparator = "<",
  value = 1000, unit = "nM", precision = "raw_no_rounding",
  citation = "doi:10.1007/s00251-019-01140-x")
indirect_binding <- function(x, base) {
  allo_keys(x, "path", c("sha256", "manifest"), "indirect evidence binding")
  path <- resolve_path(x$path, base)
  if (!file.exists(path) || dir.exists(path)) stop("Pinned indirect provider artifact unavailable")
  hash <- digest::digest(file = path, algo = "sha256", serialize = FALSE)
  if ("sha256" %in% names(x) && !identical(x$sha256, hash)) stop("Indirect provider integrity mismatch")
  m <- jsonlite::read_json(path); reject_duplicate_keys(m)
  if ("manifest" %in% names(x) && !identical(canonical_bytes(x$manifest), canonical_bytes(m)))
    stop("Pinned indirect provider content mismatch")
  indirect_provider(m)
  list(path = path, sha256 = hash, manifest = m)
}
indirect_config <- function(x, base = NULL) {
  allo_keys(x, c("version", "profile", "presenters", "provider"), label = "indirect baseline config")
  if (x$version != "indirect-hla-baseline-1.0.0" || x$profile != indirect_profile())
    stop("Unsupported literature-derived indirect-HLA profile")
  allo_records(x$presenters)
  for (p in x$presenters) {
    allo_keys(p, c("scenario_id", "subject_id", "presenter", "resolution_policy"), label = "DR presenter declaration")
    for (v in p) scalar_text(v, "presenter declaration")
    if (p$resolution_policy != "explicit_typed_alpha_beta") stop("DRA resolution must be explicit; no inferred alpha chain")
  }
  x$provider <- indirect_binding(x$provider, base)
  x
}
indirect_enabled <- function(project) !is.null(project$config$indirect_hla)
indirect_read <- function(project) {
  if (!indirect_enabled(project)) return(NULL)
  x <- indirect_config(project$config$indirect_hla)
  expected <- indirect_baseline(allo_read(project), x)
  if (!identical(canonical_bytes(expected), canonical_bytes(project$indirect_hla)))
    stop("Indirect baseline derived contract mismatch")
  expected
}
indirect_summary <- function(x) list(schema_version = x$schema_version, profile = x$profile,
  claim = x$claim, evidence_kind = x$provider$semantics$evidence_kind,
  allo_analysis_id = x$allo_analysis_id,
  summaries = lapply(x$summaries, function(s) s[setdiff(names(s), c("member_ids", "coverage"))]))
