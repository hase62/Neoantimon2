# PSM provenance, peptide observations, and experimental restriction stay separate.
ms_experiment <- function(input, settings) {
  found <- Filter(function(x) identical(x$input_id, input$input_id), settings$experiments)
  x <- if (length(found)) found[[1L]] else list()
  strict_keys(x, c("input_id", "ms_run", "region", "timepoint", "tissue", "cell_context", "ip_type",
    "hla_class", "presenter", "assignment_basis", "assignment_reference", "antibody", "reagent",
    "analysis_software", "software_version", "search_database", "database_location", "database_checksum"), "MS experiment")
  if (is.null(x$ms_run)) x$ms_run <- input$input_id
  if (is.null(x$ip_type)) x$ip_type <- "unknown"
  if (!x$ip_type %in% c("class_i_pan", "class_ii_pan", "pan_hla", "allele_specific", "monoallelic", "other", "unknown"))
    stop("Unsupported IP context")
  if (is.null(x$hla_class) && x$ip_type %in% c("class_i_pan", "class_ii_pan"))
    x$hla_class <- if (x$ip_type == "class_i_pan") "I" else "II"
  class <- recognition_mhc_class(x$hla_class)
  if (!is.null(x$hla_class) && is.null(class)) stop("Invalid MS HLA class")
  x$hla_class <- class
  for (field in c("ms_run", "assignment_basis", "assignment_reference"))
    if (!is.null(x[[field]])) scalar_text(x[[field]], paste("MS experiment", field))
  x$context <- presentation_context(x, input$sample_id)
  x
}
ms_restriction <- function(experiment) {
  p <- recognition_presenter(experiment$presenter)
  basis <- experiment$assignment_basis
  explicit <- !is.null(basis) && basis %in% c("monoallelic_system", "allele_specific_isolation", "explicit_source_assignment") &&
    !is.null(experiment$assignment_reference)
  canonical <- if (p$status == "known") peptide_presenters(list(p$mhc_key))[[1L]] else NULL
  conflict <- !is.null(experiment$hla_class) && !is.null(canonical) && experiment$hla_class != canonical$hla_class
  list(status = if (explicit && p$status == "known" && !conflict) "direct_allele_resolved" else "hla_unresolved",
    presenter = if (explicit && p$status == "known" && !conflict) p$mhc_key else NULL,
    hla_class = experiment$hla_class, raw_presenter = experiment$presenter,
    assignment_basis = basis, assignment_reference = experiment$assignment_reference,
    conflict = conflict, ip_type = experiment$ip_type)
}
ms_modifications <- function(x) {
  if (is.null(x)) return(list())
  if (is.character(x)) x <- jsonlite::fromJSON(x, simplifyVector = FALSE)
  if (!is.list(x) || !is.null(names(x))) stop("Modifications must be an array of records")
  lapply(x, function(m) {
    if (!is.list(m)) stop("Invalid modification")
    if (!is.null(m$location)) {
      pos <- presentation_number(m$location, TRUE)
      if (pos != floor(pos)) stop("Invalid modification location")
      m$location <- pos
    }
    if (!is.null(m$monoisotopicMassDelta)) m$monoisotopicMassDelta <- presentation_number(m$monoisotopicMassDelta)
    m
  })
}
ms_psm <- function(row, input, experiment) {
  if (!is.null(row[["sample_id"]]) && row[["sample_id"]] != input$sample_id) stop("PSM sample binding conflict")
  for (field in c("region", "timepoint", "tissue", "cell_context", "ms_run")) if (!is.null(row[[field]])) {
    if (!is.null(experiment[[field]]) && experiment[[field]] != row[[field]]) stop("PSM experiment context conflict")
    experiment[[field]] <- row[[field]]
  }
  experiment$context <- presentation_context(experiment, input$sample_id)
  peptide <- presentation_sequence(row[["peptide_sequence"]])
  modifications <- ms_modifications(row[["modifications"]])
  modified <- recognition_nullable(row[["modified_sequence"]])
  p <- list(psm_source_id = row[["psm_id"]], peptide_source_id = row[["peptide_id"]],
    spectrum_id = row[["spectrum_id"]], spectrum_file = row[["spectrum_file"]],
    peptide_sequence = peptide, modifications = modifications, modified_sequence = modified,
    modification_status = if (length(modifications) || (!is.null(modified) && modified != peptide)) "modified" else "unmodified",
    peptide_evidence = if (is.null(row[["peptide_evidence"]])) list() else row[["peptide_evidence"]],
    protein_accessions = if (is.null(row[["protein_accessions"]])) list() else row[["protein_accessions"]],
    decoy = presentation_bool(row[["decoy"]]), rank = presentation_number(row[["rank"]], TRUE),
    passThreshold = presentation_bool(row[["passThreshold"]]), q_value = presentation_number(row[["q_value"]], TRUE),
    pep = presentation_number(row[["pep"]], TRUE), fdr = presentation_number(row[["fdr"]], TRUE),
    fdr_local = presentation_number(row[["fdr_local"]], TRUE),
    identification_list_metadata = row[["identification_list_metadata"]],
    search_engine_score = row[["search_engine_score"]], cv_params = if (is.null(row[["cv_params"]])) list() else row[["cv_params"]],
    charge_state = presentation_number(row[["charge_state"]], TRUE),
    experimental_mz = presentation_number(row[["experimental_mz"]], TRUE), calculated_mz = presentation_number(row[["calculated_mz"]], TRUE),
    analysis_software = experiment$analysis_software, software_version = experiment$software_version,
    search_database = experiment$search_database, database_location = experiment$database_location,
    database_checksum = experiment$database_checksum, experiment = experiment,
    context = experiment$context, ms_run = experiment$ms_run, restriction = ms_restriction(experiment),
    raw_artifact = recognition_artifact(input), raw_fields = row)
  if (is.null(p$psm_source_id) || is.null(p$spectrum_id)) stop("PSM and spectrum IDs required")
  for (m in modifications) if (!is.null(m$location) && m$location > nchar(peptide) + 1L) stop("Modification beyond peptide")
  for (k in c("q_value", "pep", "fdr")) if (!is.null(p[[k]]) && p[[k]] > 1) stop("PSM quality proportion outside [0,1]")
  presentation_seal(p, "psm", "psm_id")
}
ms_generic_import <- function(input, experiment) lapply(presentation_table(input$resolved_path,
  c("psm_id", "spectrum_id", "peptide_sequence")), ms_psm, input = input, experiment = experiment)
ms_xml_attrs <- function(node) as.list(xml2::xml_attrs(node))
ms_xml_find <- function(node, tag) xml2::xml_find_all(node, paste0("./*[local-name()='", tag, "']"))
ms_xml_params <- function(node) lapply(c(ms_xml_find(node, "cvParam"), ms_xml_find(node, "userParam")), ms_xml_attrs)
ms_mzidentml_import <- function(input, experiment) {
  # NONET; DTDs/entities are outside the identification subset and rejected.
  content <- readLines(input$resolved_path, warn = FALSE)
  if (any(grepl("<!DOCTYPE|<!ENTITY", content))) stop("DTD/entity declarations not supported")
  doc <- xml2::read_xml(input$resolved_path, options = "NONET")
  if (xml2::xml_name(doc) != "MzIdentML" || !grepl("^1\\.3(\\.|$)", xml2::xml_attr(doc, "version")))
    stop("Expected mzIdentML 1.3 identification subset")
  xml_all <- function(tag) xml2::xml_find_all(doc, paste0("//*[local-name()='", tag, "']"))
  index <- function(tag) {
    nodes <- xml_all(tag); ids <- xml2::xml_attr(nodes, "id")
    if (anyNA(ids) || anyDuplicated(ids)) stop("Missing/duplicate XML ID: ", tag)
    stats::setNames(as.list(nodes), ids)
  }
  peptides <- index("Peptide"); evidence <- index("PeptideEvidence"); db <- index("DBSequence")
  databases <- index("SearchDatabase"); spectra <- index("SpectraData")
  software <- lapply(xml_all("AnalysisSoftware"), function(n) list(attributes = ms_xml_attrs(n),
    software_name = lapply(ms_xml_find(n, "SoftwareName"), ms_xml_params)))
  database_records <- lapply(databases, function(n) list(attributes = ms_xml_attrs(n), cv_params = ms_xml_params(n),
    database_name = lapply(ms_xml_find(n, "DatabaseName"), ms_xml_params)))
  out <- list()
  for (result in xml_all("SpectrumIdentificationResult")) {
    ra <- ms_xml_attrs(result); sd <- spectra[[ra$spectraData_ref]]
    if (is.null(sd)) stop("Unresolved spectrum file reference")
    for (item in ms_xml_find(result, "SpectrumIdentificationItem")) {
      a <- ms_xml_attrs(item); peptide <- peptides[[a$peptide_ref]]
      if (is.null(peptide)) stop("Unresolved peptide reference")
      refs <- lapply(ms_xml_find(item, "PeptideEvidenceRef"), ms_xml_attrs)
      pe <- lapply(refs, function(ref) {
        n <- evidence[[ref$peptideEvidence_ref]]
        if (is.null(n)) stop("Unresolved peptide evidence reference")
        v <- ms_xml_attrs(n)
        if (!identical(v$peptide_ref, a$peptide_ref)) stop("Conflicting peptide evidence")
        d <- db[[v$dBSequence_ref]]
        if (is.null(d)) stop("Unresolved DBSequence reference")
        da <- ms_xml_attrs(d)
        if (is.null(databases[[da$searchDatabase_ref]])) stop("Unresolved search database reference")
        list(attributes = v, db_sequence = da)
      })
      flags <- lapply(pe, function(p) presentation_bool(p$attributes$isDecoy))
      decoy <- if (length(flags) && all(vapply(flags, identical, TRUE, TRUE))) TRUE else
        if (length(flags) && all(vapply(flags, identical, TRUE, FALSE))) FALSE else NULL
      params <- ms_xml_params(item)
      quality <- function(accessions) {
        found <- Filter(function(p) p$accession %in% accessions, params)
        if (length(found) > 1L) stop("Multiple noncomparable quality values")
        if (length(found)) found[[1L]]$value else NULL
      }
      mods <- lapply(ms_xml_find(peptide, "Modification"), function(n) c(ms_xml_attrs(n), list(cv_params = ms_xml_params(n))))
      # SubstitutionModification is retained as modification, never erased.
      mods <- c(mods, lapply(ms_xml_find(peptide, "SubstitutionModification"), function(n)
        c(ms_xml_attrs(n), list(kind = "substitution"))))
      exp <- experiment; exp$analysis_software <- software; exp$search_database <- database_records
      row <- list(psm_id = a$id, peptide_id = a$peptide_ref, spectrum_id = ra$spectrumID,
        spectrum_file = c(list(spectra_data_ref = ra$spectraData_ref), ms_xml_attrs(sd)),
        peptide_sequence = xml2::xml_text(ms_xml_find(peptide, "PeptideSequence")), modifications = mods,
        peptide_evidence = pe, protein_accessions = lapply(pe, function(p) p$db_sequence$accession),
        decoy = decoy, rank = a$rank, passThreshold = a$passThreshold,
        charge_state = a$chargeState, experimental_mz = a$experimentalMassToCharge,
        calculated_mz = a$calculatedMassToCharge, cv_params = params,
        q_value = quality(c("MS:1002354", "MS:1001491")), pep = quality("MS:1001493"),
        fdr = NULL, fdr_local = quality("MS:1002351"),
        identification_list_metadata = list(attributes = ms_xml_attrs(xml2::xml_parent(result)),
          cv_params = ms_xml_params(xml2::xml_parent(result))),
        item_attributes = a, result_attributes = ra)
      out[[length(out) + 1L]] <- ms_psm(row, input, exp)
    }
  }
  out
}
ms_observations <- function(psms) {
  if (anyDuplicated(vapply(psms, `[[`, "", "psm_id"))) stop("Duplicate PSM identity")
  # Unknown/mixed target-decoy assignments stay unresolved, not direct observations.
  spectrum_key <- function(p) infra_hash(list(p$context, p$ms_run, p$spectrum_file, p$spectrum_id))
  groups <- split(psms, vapply(psms, spectrum_key, ""))
  for (g in groups) {
    plausible <- Filter(function(p) !isTRUE(p$decoy) && !identical(p$passThreshold, FALSE), g)
    forms <- unique(vapply(plausible, function(p) infra_hash(list(p$peptide_sequence, p$modifications, p$modified_sequence)), ""))
    for (p in g) {
      state <- if (isTRUE(p$decoy)) "decoy" else if (identical(p$passThreshold, FALSE)) "threshold_failed" else
        if (is.null(p$decoy)) "unresolved" else if (length(forms) > 1L) "ambiguous_assignment" else "unique_assignment"
      id <- match(p$psm_id, vapply(psms, `[[`, "", "psm_id"))
      psms[[id]]$assignment_status <- state
    }
  }
  key <- function(p) infra_hash(list(context = p$context, ms_run = p$ms_run, peptide = p$peptide_sequence,
    modifications = p$modifications, modified_sequence = p$modified_sequence,
    restriction = p$restriction, assignment = p$assignment_status, artifact = p$raw_artifact))
  observations <- lapply(split(psms, vapply(psms, key, "")), function(g) {
    p <- g[[1L]]
    observation <- p[c("peptide_sequence", "modifications", "modified_sequence", "modification_status",
      "context", "ms_run", "restriction", "assignment_status", "experiment", "raw_artifact")]
    observation$psm_ids <- as.list(sort(vapply(g, `[[`, "", "psm_id")))
    observation$biological_peptide_id <- biological_id("peptide_key", list(ligand_sequence = p$peptide_sequence), peptide_contract())
    observation$origin_assignment <- "unresolved_shared_peptide"
    presentation_seal(observation, "observation", "immunopeptidomics_observation_id")
  })
  list(psms = recognition_sort(psms, "psm_id"), observations = recognition_sort(observations, "immunopeptidomics_observation_id"))
}
ms_candidate_links <- function(candidates, observations) {
  links <- list()
  for (c in candidates) for (o in observations) {
    if (!identical(c$ligand_sequence, o$peptide_sequence) || !presentation_context_matches(o$context, c)) next
    restriction <- o$restriction
    compatible_class <- is.null(restriction$hla_class) || restriction$hla_class == c$hla_class
    compatible_hla <- is.null(restriction$presenter) || identical(restriction$presenter, c$mhc_key)
    status <- if (!compatible_class || !compatible_hla) "restriction_mismatch" else
      if (o$assignment_status != "unique_assignment") o$assignment_status else
        if (o$modification_status == "modified") "backbone_match_with_modification" else "exact_unmodified_ligand_match"
    links[[length(links) + 1L]] <- list(candidate_id = c$candidate_id,
      observation_id = o$immunopeptidomics_observation_id, psm_ids = o$psm_ids,
      status = status, direct_hla_assignment = if (status == "exact_unmodified_ligand_match" &&
        restriction$status == "direct_allele_resolved") "direct_allele_resolved" else "hla_unresolved",
      origin_ids = c$origin_ids, origin_assignment = "unresolved_shared_peptide")
  }
  sorted_records(links)
}
