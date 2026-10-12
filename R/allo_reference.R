# Reuse the frozen presenter syntax normalizer; reference resolution is a
# separate operation and never treats a canonical allele string as a sequence.
allo_allele <- function(raw, locus = NULL) {
  scalar_text(raw, "HLA typing")
  canonical <- peptide_presenters(list(raw))[[1L]]$canonical_input
  if (!grepl("^HLA-[A-Z0-9]+\\*[0-9]{2,3}(:[0-9]{2,3}){0,3}[NLSCAQGP]?$", canonical))
    stop("Invalid single-chain HLA typing syntax")
  found <- sub("\\*.*$", "", sub("^HLA-", "", canonical))
  if (!found %in% c("A", "B", "C", "E", "F", "G", "DRA", "DRB1", "DRB3", "DRB4", "DRB5", "DQA1", "DQB1", "DPA1", "DPB1"))
    stop("Unsupported HLA locus")
  if (!is.null(locus) && !identical(found, locus)) stop("Typing locus conflicts with allele")
  canonical
}
allo_intervals <- function(x, length_aa, label) {
  out <- lapply(allo_records(x), function(d) {
    allo_keys(d, c("name", "start", "end"), label = label)
    d$name <- allo_text(d$name, "domain/region name")
    for (k in c("start", "end")) {
      v <- d[[k]]
      if (!is.numeric(v) || length(v) != 1L || is.na(v) || v != floor(v) || v < 1L || v > length_aa)
        stop("Invalid inclusive sequence coordinates")
      d[[k]] <- as.integer(v)
    }
    if (d$start > d$end) stop("Invalid sequence interval order")
    d
  })
  allo_records(out)
}
# Molecular content identity excludes sequence-status and reference evidence.
allo_sequence_identity_payload <- function(species, r) list(species = species,
  amino_acids = r$protein_sequence, unknown_positions = r$unknown_positions,
  protein_definition = r$protein_definition, isoform = r$isoform,
  coordinate_map = r$coordinate_map, domains = r$domains,
  covered_regions = r$covered_regions)

allo_reference <- function(binding) {
  m <- binding$manifest
  required <- c("version", "status", "reference_source", "reference_release", "nomenclature_authority",
    "nomenclature_release", "species", "provider", "parser", "coordinate_convention",
    "isoform_policy", "completion_policy", "sequence_policy", "records")
  allo_keys(m, required, c("provenance", "reason_codes"), "HLA reference manifest")
  if ("provenance" %in% names(m)) allo_provenance(m$provenance)
  if ("reason_codes" %in% names(m)) allo_strings(m$reason_codes, "reference reason codes")
  if (m$version != "hla-reference-1.0.0") stop("Unsupported HLA reference schema")
  for (k in setdiff(required, c("provider", "parser", "records"))) scalar_text(m[[k]], k)
  if (!m$status %in% c("evaluated", "failed") || m$species != "NCBITaxon:9606" ||
      m$completion_policy != "none" || m$isoform_policy != "explicit_mapped_isoforms" ||
      !m$sequence_policy %in% c("precursor", "mature") ||
      m$coordinate_convention != "one_based_inclusive") stop("Unsupported reference policy; completion/imputation forbidden")
  for (k in c("provider", "parser")) {
    allo_keys(m[[k]], c("name", "version"), label = k)
    for (s in m[[k]]) scalar_text(s, k)
  }
  if (m$parser$name != "literal_amino_acids" || m$parser$version != "1.0.0")
    stop("Unsupported reference parser; alignment notation must not be stripped")
  records <- lapply(allo_records(m$records), function(r) {
    allo_keys(r, c("allele", "sequence_identifier", "protein_sequence", "sequence_status", "supported",
      "protein_definition", "isoform", "coordinate_map", "domains", "covered_regions", "resolution_provenance"),
      label = "HLA reference sequence", nullable = "protein_sequence")
    r$allele <- allo_allele(r$allele)
    if (!grepl(":", r$allele) || grepl("[GP]$", r$allele)) stop("Reference mapping requires an explicit allele, not a group")
    for (k in c("sequence_identifier", "protein_definition", "isoform", "resolution_provenance"))
      r[[k]] <- allo_text(r[[k]], k)
    r$supported <- allo_boolean(r$supported, "supported")
    if (!r$sequence_status %in% c("complete", "partial", "unavailable")) stop("Invalid sequence status")
    allo_keys(r$coordinate_map, c("sequence_coordinates", "precursor_start", "mature_start"), label = "coordinate map", nullable = c("precursor_start", "mature_start"))
    if (!identical(r$coordinate_map$sequence_coordinates, m$sequence_policy)) stop("Sequence coordinate policy mismatch")
    for (k in c("precursor_start", "mature_start")) {
      v <- r$coordinate_map[[k]]
      if (!is.null(v) && (!is.numeric(v) || length(v) != 1L || is.na(v) || v != floor(v) || v < 1L)) stop("Invalid coordinate map")
      if (!is.null(v)) r$coordinate_map[[k]] <- as.integer(v)
    }
    if (is.null(r$coordinate_map[[paste0(m$sequence_policy, "_start")]]) ||
      r$coordinate_map[[paste0(m$sequence_policy, "_start")]] != 1L) stop("Declared sequence coordinate origin required")
    if (r$sequence_status == "unavailable") {
      if (!is.null(r$protein_sequence) || length(r$domains) || length(r$covered_regions)) stop("Unavailable sequence must not fabricate residues/coverage")
      r["hla_source_sequence_id"] <- list(NULL)
      r$unknown_positions <- list()
    } else {
      scalar_text(r$protein_sequence, "protein sequence")
      if (!grepl("^[ACDEFGHIKLMNPQRSTVWYX]+$", r$protein_sequence)) stop("Invalid literal protein sequence; no alignment notation")
      n <- nchar(r$protein_sequence)
      r$domains <- allo_intervals(r$domains, n, "selected domain")
      r$covered_regions <- allo_intervals(r$covered_regions, n, "covered region")
      if (!length(r$domains)) stop("Selected domain intervals required")
      unknown <- which(strsplit(r$protein_sequence, "", fixed = TRUE)[[1L]] == "X")
      covered <- unique(unlist(lapply(r$covered_regions, function(d) seq.int(d$start, d$end))))
      if (any(unknown %in% covered)) stop("Unknown residues cannot be covered")
      if (r$sequence_status == "complete" && (length(unknown) || !setequal(covered, seq_len(n)))) stop("Complete sequence requires literal full coverage")
      r$unknown_positions <- as.list(as.integer(unknown))
      r$hla_source_sequence_id <- allo_id("sequence", allo_sequence_identity_payload(m$species, r))
    }
    r["expression_annotation"] <- list(if (grepl("[NLSCAQ]$", r$allele)) substr(r$allele, nchar(r$allele), nchar(r$allele)) else NULL)
    r
  })
  records <- allo_records(records)
  mapping_keys <- vapply(records, function(r) canonical_bytes(r[c("allele", "sequence_identifier")]), "")
  if (anyDuplicated(mapping_keys)) stop("Conflicting duplicate reference allele/sequence mapping")
  pin <- m[setdiff(required, c("records", "status"))]
  pin$records <- allo_records(lapply(records, function(r) r[setdiff(names(r), "resolution_provenance")]))
  # Sequence aliases are mapping configuration; unlike the content ID above,
  # this policy identity pins nomenclature/database release and actual mapping.
  policy <- list(profile = "neo-allo-phase-a-reference-v1", version = "1.0.0", reference = pin,
    typing_imputation = "none", resolution = "all_compatible_pinned_sequences_v1")
  list(schema_version = "hla-reference-1.0.0", reference_policy_id = allo_id("policy", policy),
    semantics = pin[setdiff(names(pin), "records")], records = records,
    status = m$status, reason_codes = if (m$status == "failed") allo_set(c(list("reference_provider_failure"), m$reason_codes)) else list(),
    provenance = list(artifact = list(path = binding$path, sha256 = binding$sha256),
      metadata = if (is.null(m$provenance)) list() else allo_provenance(m$provenance)))
}
allo_resolve_typing <- function(config, reference) {
  subjects <- vapply(config$subjects, function(s) allo_text(s$subject_id, "subject_id"), "")
  out <- lapply(config$typing, function(t) {
    allo_keys(t, c("subject_id", "typing_scenario_id", "locus", "alleles", "complete", "copy_number", "phase"), label = "HLA typing", nullable = c("copy_number", "phase"))
    for (k in c("subject_id", "typing_scenario_id", "locus")) t[[k]] <- allo_text(t[[k]], k)
    if (!t$subject_id %in% subjects) stop("Unknown typing subject")
    if (!t$locus %in% c("A", "B", "C", "E", "F", "G", "DRA", "DRB1", "DRB3", "DRB4", "DRB5", "DQA1", "DQB1", "DPA1", "DPB1")) stop("Unsupported typing locus")
    t$complete <- allo_boolean(t$complete, "typing complete")
    if (!is.null(t$copy_number) && (!is.numeric(t$copy_number) || length(t$copy_number) != 1L || is.na(t$copy_number) || t$copy_number != floor(t$copy_number) || t$copy_number < 1L)) stop("Invalid typing copy number")
    if (!is.null(t$copy_number)) t$copy_number <- as.integer(t$copy_number)
    if (!is.null(t$phase)) t$phase <- allo_text(t$phase, "typing phase")
    raw <- allo_strings(t$alleles, "typing alleles")
    canonical <- allo_set(lapply(raw, allo_allele, locus = t$locus))
    matches <- allo_records(unlist(lapply(canonical, function(a) {
      # A group label does not define membership by allele-prefix matching.
      # This provider has no approved group expansion/imputation operation.
      if (grepl("[GP]$", a)) return(list())
      stem <- sub("[NLSCAQGP]$", "", a)
      suffix <- if (grepl("[NLSCAQ]$", a)) substr(a, nchar(a), nchar(a)) else NULL
      Filter(function(r) {
        rs <- sub("[NLSCAQ]$", "", r$allele)
        (rs == stem || startsWith(rs, paste0(stem, ":"))) &&
          (is.null(suffix) || identical(r$expression_annotation, suffix))
      }, reference$records)
    }), recursive = FALSE))
    # Every explicit alternative, including its suffix, must have a mapping.
    missing_alternative <- any(vapply(canonical, function(a) {
      stem <- sub("[NLSCAQGP]$", "", a)
      suffix <- if (grepl("[NLSCAQ]$", a)) substr(a, nchar(a), nchar(a)) else NULL
      !any(vapply(matches, function(r) {
        rs <- sub("[NLSCAQ]$", "", r$allele)
        (rs == stem || startsWith(rs, paste0(stem, ":"))) &&
          (is.null(suffix) || identical(r$expression_annotation, suffix))
      }, TRUE))
    }, TRUE))
    status <- "evaluated"; reasons <- list()
    if (reference$status == "failed") { status <- "failed"; reasons <- reference$reason_codes }
    else if (!length(canonical) || any(!grepl(":", unlist(canonical))) || any(grepl("[GP]$", unlist(canonical)))) {
      status <- "insufficient"; reasons <- list("hla_resolution_insufficient")
    } else if (!length(matches) || missing_alternative || any(!vapply(matches, `[[`, TRUE, "supported"))) {
      status <- "unsupported"; reasons <- list("reference_allele_unsupported")
    } else if (any(vapply(matches, function(r) is.null(r$hla_source_sequence_id), TRUE))) {
      status <- "unavailable"; reasons <- list("source_sequence_unavailable")
    } else if (length(unique(vapply(matches, `[[`, "", "hla_source_sequence_id"))) > 1L ||
      length(unique(vapply(matches, function(r) if (is.null(r$expression_annotation)) "none" else r$expression_annotation, ""))) > 1L) {
      status <- "ambiguous"; reasons <- list("hla_typing_ambiguous")
    } else if (!t$complete || any(vapply(matches, function(r) r$sequence_status != "complete", TRUE))) {
      status <- "insufficient"; reasons <- if (!t$complete) list("hla_resolution_insufficient") else list("source_sequence_partial")
    }
    exact <- length(canonical) == 1L && length(matches) == 1L && identical(canonical[[1L]], matches[[1L]]$allele)
    list(subject_id = t$subject_id, typing_scenario_id = t$typing_scenario_id, locus = t$locus,
      hla_class = if (t$locus %in% c("A", "B", "C", "E", "F", "G")) "I" else "II",
      raw_alleles = t$alleles, normalized_alleles = canonical, complete = t$complete,
      copy_number = t$copy_number, phase = t$phase, compatible_alleles = allo_set(lapply(matches, `[[`, "allele")),
      sequence_ids = allo_set(lapply(matches, `[[`, "hla_source_sequence_id")),
      expression_annotations = allo_set(lapply(matches, `[[`, "expression_annotation")),
      resolved_allele = if (status == "evaluated" && exact) canonical[[1L]] else NULL,
      resolution_operation = if (status != "evaluated") "unresolved" else if (exact) "exact_pinned_mapping" else "pinned_sequence_equivalence_without_genotype_imputation",
      evaluability = allo_state(status, reasons, as.list(canonical), lapply(matches, `[[`, "allele"), "pinned_hla_protein_resolution"))
  })
  # Raw spelling/order remains provenance but never controls identity or lookup.
  allo_records(out)
}
