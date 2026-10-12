# Cross-record constraints supplement (never duplicate) field definitions.
validate_snapshot <- function(x, contract = schema_contract()) {
  fail <- function(msg) stop(msg, call. = FALSE)
  eq <- function(a, b, msg) if (!identical(canonical_bytes(a), canonical_bytes(b))) fail(msg)
  present <- function(v) !is_missing(v)
  lookup <- function(rows, key, value) {
    hits <- Filter(function(r) identical(r[[key]], value), rows)
    if (length(hits) != 1L) fail(paste("Unresolved/ambiguous foreign key", key, value))
    hits[[1L]]
  }
  expect_id <- function(value, kind, vals) eq(value, biological_id(kind, vals, contract), paste("Invalid", kind))
  for (t in c("events", "peptide_origins", "neoantigen_results", "run_manifest"))
    if (!t %in% names(x)) fail(paste("Missing output", t))
  m <- x$run_manifest
  validate_record(m, contract$tables$run_manifest, contract, "run_manifest")
  for (t in setdiff(names(contract$tables), "run_manifest")) {
    if (t %in% names(x)) validate_table(x[[t]], t, contract)
    expected_version <- if (is.null(contract$output_schema_versions))
      contract$schema_version else contract$output_schema_versions[[t]]
    if (!identical(m$output_schema_versions[[t]], expected_version)) fail("Output schema version mismatch")
  }
  contexts <- m$analysis_contexts
  ids <- vapply(contexts, `[[`, "", "context_id")
  if (anyDuplicated(ids)) fail("Duplicate context")
  for (ct in contexts) {
    for (role in c("source", "self", "presenter")) {
      if (!ct[[paste0(role, "_subject_role")]] %in% unlist(contract$role_bindings[[ct$mode]][[role]])) fail("Invalid subject role for mode")
    }
    if (ct$mode == "somatic") {
      eq(ct$source_subject_id, ct$self_subject_id, "Somatic subjects must identify the same patient")
      eq(ct$self_subject_id, ct$presenter_subject_id, "Somatic presenter must identify patient")
      if (ct$direction != "within_patient") fail("Invalid somatic direction")
    } else {
      if (!ct$direction %in% c("donor_to_recipient", "recipient_to_donor")) fail("Explicit allo direction required")
      forward <- ct$direction == "donor_to_recipient"
      eq(ct$source_subject_role, if (forward) "donor" else "recipient", "Source/direction mismatch")
      eq(ct$self_subject_role, if (forward) "recipient" else "donor", "Self/direction mismatch")
      eq(ct$presenter_subject_role, ct$self_subject_role, "Presenter/direction mismatch")
      if (ct$source_subject_id == ct$self_subject_id) fail("Allo source and self must differ")
      eq(ct$self_subject_id, ct$presenter_subject_id, "Allo presenter/self mismatch")
      if (is.null(m$allo_policy)) fail("Allo policy metadata required")
    }
    for (n in c("organ", "tissue", "timepoint")) if (is_missing(ct[[n]]) && is.null(ct$missing_reasons[[n]])) fail("Missing context reason")
  }
  if (!is.null(m$allo_policy)) {
    p <- m$allo_policy
    if (p$policy_status == "NOT_RUN" && is_missing(p$reason)) fail("Missing policy reason")
    if (p$policy_status == "DEFINED" && any(vapply(p[setdiff(names(p), "reason")], is_missing, TRUE))) fail("Incomplete defined allo policy")
  }
  for (t in c("events", "peptide_origins", "neoantigen_results")) for (r in x[[t]]) {
    ct <- lookup(contexts, "context_id", r$analysis_context)
    eq(r$analysis_id, m$analysis_id, "Analysis mismatch")
    if (!r$sample_id %in% unlist(ct$sample_bindings)) fail("Sample not bound to context")
    eq(r$patient_id, if (ct$mode == "somatic") ct$self_subject_id else
      if (ct$self_subject_role == "recipient") ct$self_subject_id else ct$source_subject_id, "Patient compatibility field mismatch")
    for (n in grep("_subject_(id|role)$", names(r), value = TRUE)) if (present(r[[n]])) eq(r[[n]], ct[[n]], "Role projection mismatch")
  }
  for (e in x$events) {
    expect_id(e$event_key, "event_key", e[c("genome_assembly", "canonical_event")])
    expect_id(e$event_id, "event_id", list(context_id=e$analysis_context, patient_id=e$patient_id, sample_id=e$sample_id, event_key=e$event_key))
    for (assay in c("tumor", "normal", "rna")) {
      names_count <- paste0(assay, "_", c("raw_depth", "usable_depth", "ref_count", "alt_count", "other_count", "ambiguous_count", "vaf_denominator"))
      measured <- any(vapply(e[names_count], present, TRUE)) || present(e[[paste0(assay,"_vaf")]])
      if (measured && (is_missing(e[[paste0(assay,"_count_unit")]]) || !identical(e[[paste0(assay,"_assessment_status")]], "EVALUATED"))) fail("Count requires unit and EVALUATED status")
      if (identical(e[[paste0(assay,"_vaf_denominator")]], 0L) || identical(e[[paste0(assay,"_vaf_denominator")]], 0))
        if (present(e[[paste0(assay,"_vaf")]])) fail("Zero denominator requires missing VAF")
    }
    for (n in c("dna_split", "dna_spanning", "rna_junction", "rna_spanning"))
      if (present(e[[paste0(n,"_support")]]) && is_missing(e[[paste0(n,"_count_unit")]])) fail("Support count unit missing")
  }
  for (o in x$peptide_origins) {
    eq(o$ligand_sequence, o$peptide_sequence, "Ligand alias mismatch")
    if (o$peptide_length != nchar(o$ligand_sequence) || o$peptide_end_aa-o$peptide_start_aa+1 != o$peptide_length) fail("Peptide coordinates/length mismatch")
    expect_id(o$peptide_key, "peptide_key", o["ligand_sequence"])
    expect_id(o$sample_peptide_id, "sample_peptide_id", list(context_id=o$analysis_context,patient_id=o$patient_id,sample_id=o$sample_id,peptide_key=o$peptide_key))
    expect_id(o$origin_id, "origin_id", o[c("sample_peptide_id", "origin_identity")])
    ii <- o$origin_identity
    for (n in c("transcript_id","transcript_version","peptide_start_aa","peptide_end_aa","phase_namespace","phase_set","haplotype_id","source_antigen_key")) eq(o[[n]], ii[[n]], "Origin identity projection mismatch")
    ct <- lookup(contexts,"context_id",o$analysis_context)
    eq(ii$source_subject_id,ct$source_subject_id,"Origin source mismatch")
    eq(ii$scenario_id,ct$scenario_id,"Origin scenario mismatch")
    keys <- character()
    for (n in c("source_event_ids", "context_event_ids", "germline_context_ids")) for (eid in o[[n]]) {
      e <- lookup(x$events, "event_id", eid)
      eq(e$sample_id, o$sample_id, "Origin event sample mismatch")
      eq(e$analysis_context, o$analysis_context, "Origin event context mismatch")
      if (n == "source_event_ids") keys <- c(keys, e$event_key)
    }
    eq(as.list(sort(unique(keys),method="radix")),ii$source_event_keys,"Origin event identity set mismatch")
    if (present(o$rna_origin_support_count) && is_missing(o$rna_origin_count_unit)) fail("Origin count unit missing")
    if (length(o$source_metadata_json)) {
      validate_record(o$source_metadata_json,contract$records$source_sequence_detail,contract,"source_metadata_json")
      eq(o$source_metadata_json$source_antigen_key,o$source_antigen_key,"Source metadata antigen mismatch")
    }
    if (ct$mode == "allo_hla") {
      if (is_missing(o$source_antigen_key) || is_missing(o$source_sequence_ref)) fail("Allo source identity missing")
      if (any(vapply(o[c("source_protein_id","source_allele","source_locus")],is_missing,TRUE)) && !length(o$source_metadata_json$uncertainty)) fail("Unresolved allo source needs uncertainty")
    }
  }
  for (r in x$neoantigen_results) {
    ct <- lookup(contexts,"context_id",r$analysis_context)
    eq(r$ligand_sequence,r$peptide_sequence,"Ligand alias mismatch")
    if (r$peptide_length != nchar(r$ligand_sequence)) fail("Ligand length mismatch")
    expect_id(r$peptide_key,"peptide_key",r["ligand_sequence"])
    expect_id(r$candidate_key,"candidate_key",r[c("species","hla_class","ligand_sequence","mhc_key")])
    expect_id(r$candidate_id,"candidate_id",list(context_id=r$analysis_context,sample_id=r$sample_id,candidate_key=r$candidate_key))
    if (ct$mode != "somatic" && r$final_grade != "NOT_GRADED") fail("Allo grade not implemented")
    os <- Filter(function(o) identical(o$sample_peptide_id,r$sample_peptide_id),x$peptide_origins)
    if (!length(os) || !setequal(unlist(r$origin_ids),vapply(os,`[[`,"","origin_id")) || r$n_origins != length(os)) fail("All origins must be preserved")
    for (o in os) {
      eq(o$analysis_context,r$analysis_context,"Origin/result context mismatch")
      eq(o$sample_id,r$sample_id,"Origin/result sample mismatch")
      eq(o$peptide_key,r$peptide_key,"Origin/result peptide mismatch")
    }
    for (n in c("grading_origin_id","representative_origin_id")) if (present(r[[n]]) && !r[[n]] %in% unlist(r$origin_ids)) fail("Selected origin not in candidate")
    # Display projections must never assemble independent origins' best evidence.
    projected <- c("tumor_alt_count","normal_alt_count","rna_alt_count","tumor_vaf","normal_vaf","rna_vaf","ccf","gene_tpm","transcript_tpm")
    if (any(vapply(r[projected],present,TRUE))) {
      o <- lookup(os,"origin_id",r$grading_origin_id)
      for (n in projected) if (present(r[[n]])) {
        if (n %in% c("gene_tpm","transcript_tpm")) eq(r[[n]],o[[n]],"Cross-origin expression projection")
        else {
          if (length(o$source_event_ids) != 1L) fail("Multi-event projection requires future explicit evaluation; scalar must be null")
          e <- lookup(x$events,"event_id",o$source_event_ids[[1L]])
          eq(r[[n]],e[[n]],"Cross-origin evidence projection")
        }
      }
    }
    if (r$final_grade == "INVALID" && !all(vapply(os,function(o) identical(o$origin_validity,"INVALID"),TRUE))) fail("Valid candidate is not INVALID")
    expected_status <- if (r$final_grade %in% c("INVALID","NOT_GRADED")) r$final_grade else "GRADED"
    eq(r$grading_status,expected_status,"Grade/status mismatch")
    if (r$final_grade == "NOT_GRADED" && !length(r$not_evaluated_codes)) fail("NOT_GRADED requires reason")
    if (r$prediction_status %in% c("NOT_RUN","FAILED","UNSUPPORTED","NOT_APPLICABLE")) {
      if (present(r$mhc_rank_pct) || present(r$mhc_affinity_nm)) fail("Unperformed prediction is not a negative value")
      if (!r$final_grade %in% c("NOT_GRADED","INVALID")) fail("Unperformed prediction requires NOT_GRADED")
    }
    if (any(vapply(r[c("self_exact_match","self_match_canonical","self_match_protein_count")],present,TRUE)) && !identical(r$self_status,"EVALUATED")) fail("Self match value requires evaluated search")
    assignments <- r$epitope_assignments_json
    if (anyDuplicated(vapply(assignments,`[[`,"","assignment_id"))) fail("Duplicate assignment")
    primary <- Filter(function(z) isTRUE(z$is_primary),assignments)
    if (length(primary)>1L) fail("Multiple primary assignments")
    for (z in assignments) {
      eq(z$candidate_id,r$candidate_id,"Assignment candidate mismatch")
      pr <- lookup(m$predictions,"prediction_ref",z$prediction_ref)
      for (n in c("request_id","predictor_version","model_checksum")) eq(z[[n]],pr[[n]],"Prediction provenance mismatch")
      if (z$status == "EVALUATED") {
        if (r$mhc_key == "UNSPECIFIED" || r$hla_class == "UNSPECIFIED") fail("Epitope requires resolved presenter")
        if (any(vapply(z[c("binding_core","binding_core_start","binding_core_length","binding_orientation","epitope_key")],is_missing,TRUE))) fail("Evaluated assignment lacks core")
        if (z$binding_orientation != "forward") fail("Noncanonical mapping requires NOT_EVALUABLE")
        start <- z$binding_core_start; len <- z$binding_core_length
        if (start+len-1 > nchar(r$ligand_sequence) || substr(r$ligand_sequence,start,start+len-1) != z$binding_core || nchar(z$binding_core)!=len) fail("Invalid core coordinates")
        if (length(z$position_mapping) && !identical(as.integer(unlist(z$position_mapping)),seq.int(start,length.out=len))) fail("Invalid position mapping")
        vals <- c(r[c("species","hla_class")],z["binding_core"],r["mhc_key"])
        expect_id(z$epitope_key,"epitope_key",vals)
      } else if (present(z$epitope_key) || is_missing(z$reason)) fail("Unresolved assignment needs null key and reason")
    }
    for (n in c("binding_core","binding_core_start","binding_core_length","binding_orientation","epitope_key")) {
      if (length(primary)) eq(r[[n]],primary[[1L]][[n]],"Primary core projection mismatch")
      else if (present(r[[n]])) fail("Scalar core requires primary assignment")
    }
    if (length(primary)) eq(r$primary_prediction_ref,primary[[1L]]$prediction_ref,"Primary prediction mismatch")
    if (length(assignments)>1L && !length(primary) && !"AMBIGUOUS_CORE" %in% unlist(r$not_evaluated_codes)) fail("Ambiguous core reason required")
    if (anyDuplicated(vapply(r$evidence_records_json,`[[`,"","evidence_id"))) fail("Duplicate evidence")
    for (ev in r$evidence_records_json) {
      eq(ev$analysis_context,r$analysis_context,"Evidence context mismatch")
      targets <- switch(ev$target_type,candidate=r$candidate_id,epitope_assignment=vapply(assignments,`[[`,"","assignment_id"),origin=unlist(r$origin_ids),source_antigen=unlist(lapply(os,`[[`,"source_antigen_key")),presenter=r$mhc_key)
      if (!ev$target_id %in% targets) fail("Evidence target foreign key mismatch")
      if (ev$evidence_type %in% c("source_antigen_expression","processing_plausibility","observed_proteomics","structural_b_cell","structural_accessibility") && ev$target_type == "candidate" && length(os)>1L) fail("Origin-dependent evidence needs explicit origin/source antigen target")
      if (ev$assessment_status != "EVALUATED") {
        if (present(ev$value) || present(ev$finding) || is_missing(ev$reason)) fail("Missing evidence cannot become negative")
      } else {
        if (is_missing(ev$value) || is_missing(ev$method) || is_missing(ev$source_ref) || is_missing(ev$metric_id) || is_missing(ev$unit)) fail("Evaluated evidence requires value/unit/provenance")
        metric <- contract$metrics[[ev$evidence_type]]
        if (!is.null(metric)) {
          validate_value(ev$value,c(metric,list(nullable=FALSE)),contract,"evidence.value")
          eq(ev$unit,metric$unit,"Evidence metric unit mismatch")
        }
        if (identical(ev$finding,"NEGATIVE") && (is.null(ev$details$scope) || !nzchar(ev$details$scope))) fail("Negative evidence requires evaluated scope")
        if (ev$evidence_type %in% c("observed_proteomics","observed_immunopeptidomics","tcr_evidence") && identical(ev$finding,"NEGATIVE") && is_missing(ev$assay)) fail("Observed negative needs assay")
        if (grepl("self_match$",ev$evidence_type)) {
          if ((identical(ev$value,FALSE) && !identical(ev$finding,"NEGATIVE")) || (identical(ev$value,TRUE) && !identical(ev$finding,"POSITIVE"))) fail("Self value/finding mismatch")
          detail <- ev$details$self_match
          validate_record(detail,contract$records$self_match_detail,contract,"self_match")
          eq(detail$self_subject_id,ct$self_subject_id,"Self comparator mismatch")
        }
      }
    }
    for (link in r$linked_recognition_json) {
      eq(link$analysis_context,r$analysis_context,"Link context mismatch")
      if (present(link$source_antigen_complex_key)) fail("Complex membership requires future verified registry")
      for (oid in link$t_cell_origin_ids) {
        o <- lookup(os,"origin_id",oid)
        eq(o$source_antigen_key,link$source_antigen_key,"Link crosses source antigens")
      }
      if (!all(unlist(link$candidate_ids) %in% r$candidate_id) || !all(unlist(link$assignment_ids) %in% vapply(assignments,`[[`,"","assignment_id"))) fail("Link foreign key mismatch")
      for (eid in link$b_cell_evidence_ids) {
        ev <- lookup(r$evidence_records_json,"evidence_id",eid)
        if (ev$target_type == "source_antigen") eq(ev$target_id,link$source_antigen_key,"B evidence antigen mismatch")
        else if (ev$target_type == "origin") eq(lookup(os,"origin_id",ev$target_id)$source_antigen_key,link$source_antigen_key,"B evidence origin mismatch")
        else fail("B evidence must target source antigen/origin")
      }
    }
  }
  for (s in x$multisample_summary) {
    ct <- lookup(contexts,"context_id",s$analysis_context)
    if (ct$mode != "somatic") fail("Allo multisample is unsupported")
    eq(s$patient_id,ct$self_subject_id,"Summary patient/context mismatch")
    eq(s$analysis_id,m$analysis_id,"Summary analysis mismatch")
    if (s$n_samples_peptide_assessed+s$n_samples_not_assessed != s$n_samples_total || s$n_samples_detected>s$n_samples_peptide_assessed) fail("Summary denominators inconsistent")
  }
  invisible(TRUE)
}
