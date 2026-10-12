apm_component_version <- function() "apm-component-map-1.0.0"
apm_component_map <- function() {
  groups <- list(class_i_hla = c("HLA-A", "HLA-B", "HLA-C"), class_i_structural = "B2M",
    class_i_transport = c("TAP1", "TAP2"), class_i_loading = c("TAPBP", "CALR", "PDIA3"),
    class_i_trimming = c("ERAP1", "ERAP2"), class_i_immunoproteasome = c("PSMB8", "PSMB9", "PSMB10"),
    class_i_regulation = "NLRC5", class_ii_hla = c("HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5",
      "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1"), class_ii_loading = c("CD74", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB"),
    class_ii_regulation = "CIITA")
  unlist(lapply(names(groups), function(role) lapply(groups[[role]], function(g) list(gene = g,
    role = role, hla_class = if (startsWith(role, "class_i_")) "I" else "II"))), recursive = FALSE)
}
apm_import <- function(input) {
  rows <- if (input$input_role == "arcas_hla_quant") arcas_quant_rows(input) else
    presentation_table(input$resolved_path, if (input$input_role == "hla_expression") "allele" else c("gene", "measurement_type"))
  lapply(rows, function(r) {
    allele <- recognition_nullable(r[["allele"]])
    canonical <- if (is.null(allele)) NULL else peptide_presenters(list(allele))[[1L]]
    if (!is.null(canonical) && (is.null(canonical$canonical_input) || grepl("~", canonical$canonical_input, fixed = TRUE)))
      stop("Allele expression requires a single canonical HLA chain")
    gene <- r[["gene"]]
    if (!is.null(canonical)) {
      allele_gene <- sub("\\*.*$", "", canonical$canonical_input)
      if (!is.null(gene) && gene != allele_gene) stop("HLA gene/allele conflict")
      gene <- allele_gene
    }
    scalar_text(gene, "APM gene")
    found <- Filter(function(c) c$gene == gene, apm_component_map())
    role <- if (length(found)) found[[1L]]$role else "outside_core_apm_map"
    class <- if (length(found)) found[[1L]]$hla_class else NULL
    type <- if (input$input_role == "hla_expression") "allele_expression" else r[["measurement_type"]]
    if (!type %in% c("gene_expression", "allele_expression", "copy_number", "somatic_alteration", "protein_ihc", "functional_annotation"))
      stop("Unsupported APM measurement type")
    if (type == "allele_expression" && is.null(allele)) stop("Allele expression identity missing")
    raw_count <- presentation_number(r[["raw_count"]], TRUE)
    abundance <- presentation_number(r[["abundance"]], TRUE); fraction <- presentation_number(r[["allele_fraction"]], TRUE)
    if (!is.null(fraction) && fraction > 1) stop("Allele fraction outside [0,1]")
    if (!is.null(abundance) && is.null(r[["unit"]])) stop("Abundance unit required")
    strong <- r[["loss_support"]] %in% c("documented_biallelic_loss", "documented_homozygous_deletion",
      "validated_protein_loss", "provider_functional_null") && !is.null(r[["interpretation_reference"]])
    if (is.null(r[["loss_support"]])) strong <- FALSE
    state <- if (strong) "direct_loss_supported" else if (!is.null(r[["provider_loss_call"]])) "provider_loss_call" else
      if (!is.null(r[["abnormality_annotation"]])) "abnormality_observed" else "evidence_present_unresolved"
    e <- list(component_map_version = apm_component_version(), gene = gene, gene_id = r[["gene_id"]],
      allele = if (is.null(canonical)) NULL else canonical$canonical_input, role = role, hla_class = class,
      measurement_type = type, raw_count = raw_count, abundance = abundance, unit = r[["unit"]],
      unit_definition = r[["unit_definition"]], allele_fraction = fraction, copy_number = presentation_number(r[["copy_number"]], TRUE),
      event_id = r[["event_id"]], variant_type = r[["variant_type"]], consequence = r[["consequence"]],
      zygosity = r[["zygosity"]], loh = presentation_bool(r[["loh"]]), functional_annotation = r[["functional_annotation"]],
      pathogenicity = r[["pathogenicity"]], loss_support = r[["loss_support"]],
      provider_loss_call = r[["provider_loss_call"]], abnormality_annotation = r[["abnormality_annotation"]],
      interpretation_reference = r[["interpretation_reference"]], state = state,
      provider = r[["provider"]], provider_version = r[["provider_version"]], reference_identity = r[["reference_identity"]],
      mapping_qc = r[["mapping_qc"]], context = presentation_context(r, input$sample_id),
      raw_artifact = recognition_artifact(input), raw_fields = r)
    presentation_seal(e, "apm", "apm_observation_id")
  })
}
apm_candidate_links <- function(candidates, observations) {
  links <- list()
  for (c in candidates) for (e in observations) {
    if (!presentation_context_matches(e$context, c)) next
    chains <- strsplit(c$mhc_key, "~", fixed = TRUE)[[1L]]
    genes <- sub("\\*.*$", "", chains)
    applicable <- !is.null(e$hla_class) && e$hla_class == c$hla_class &&
      (!grepl("_hla$", e$role) || e$gene %in% genes) && (is.null(e$allele) || e$allele %in% chains)
    links[[length(links) + 1L]] <- list(candidate_id = c$candidate_id, apm_observation_id = e$apm_observation_id,
      gene = e$gene, allele = e$allele, measurement_type = e$measurement_type,
      applicability = if (applicable) "applicable" else "not_applicable",
      state = if (applicable) e$state else "not_applicable")
  }
  sorted_records(links)
}
# Documented quant.alleles.tsv / quant.genes.tsv variants; no positional guesses.
arcas_quant_rows <- function(input) {
  rows <- presentation_table(input$resolved_path, "gene")
  unlist(lapply(rows, function(r) {
    gene <- if (startsWith(r[["gene"]], "HLA-")) r[["gene"]] else paste0("HLA-", r[["gene"]])
    common <- list(gene = gene, provider = "arcasHLA", provider_version = input$version,
      reference_identity = r[["reference_identity"]], mapping_qc = r[["mapping_qc"]], unit = "TPM",
      unit_definition = "provider-reported transcripts per million", provider_fields = r)
    if (all(c("allele1", "allele1_count", "allele1_tpm") %in% names(r))) {
      if (any(c("allele2", "allele2_count", "allele2_tpm") %in% names(r)) &&
          !all(c("allele2", "allele2_count", "allele2_tpm") %in% names(r))) stop("Incomplete arcasHLA allele2 header")
      slots <- if ("allele2" %in% names(r)) 1:2 else 1L
      Filter(Negate(is.null), lapply(slots, function(i) {
        allele <- r[[paste0("allele", i)]]
        if (is.null(allele)) {
          if (!is.null(r[[paste0("allele", i, "_count")]]) || !is.null(r[[paste0("allele", i, "_tpm")]]))
            stop("arcasHLA quantity without allele identity")
          return(NULL)
        }
        c(common, list(allele = allele, measurement_type = "allele_expression",
          raw_count = r[[paste0("allele", i, "_count")]], abundance = r[[paste0("allele", i, "_tpm")]],
          provider_allele_slot = i))
      }))
    } else if (all(c("count", "tpm", "abundance") %in% names(r)))
      list(c(common, list(measurement_type = "gene_expression", raw_count = r[["count"]], abundance = r[["tpm"]]))) else
        stop("Unsupported arcasHLA quant named-field variant")
  }), recursive = FALSE)
}
