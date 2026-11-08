# Optional presentation only. All science and graph construction stay in neo_report_bundle.
html_report_id_family <- function(key) {
  if(grepl("candidate_(id|ids|key)$",key)) return("Candidate")
  if(grepl("origin_(id|ids)$",key)) return("Origin")
  if(grepl("event_refs?$",key)) return("EventRef")
  if(grepl("event_ids?$|context_ids$",key)) return("Event")
  if(grepl("sample_id$|patient_id$|subject_id$",key)) return("Subject")
  if(grepl("analysis_id$|run_id$|input_id$|evidence_id$|prediction_record_id$|scope_id$|source_set_id$|peptide_unit_id$|peptide_key$|sample_peptide_id$|sequence_context_id$|haplotype_id$|phase_namespace$|sequence_ref$",key)) return("Record")
  NULL
}
html_report_fields <- function() c(
  "candidate_id","candidate_key","sample_id","patient_id","peptide_sequence","peptide_length","ligand_sequence","mhc_key","hla_allele","hla_class","hla_alpha","hla_beta","origin_ids","n_origins","peptide_unit_id",
  "html_native_prediction_status","html_native_grade_status","html_native_reasons","html_native_status_source","native_grade","native_rank","native_rank_status","native_policy_version","final_grade","grading_status","grade_before_warnings","rule_version","prediction_status","primary_predictor","primary_metric","mhc_affinity_nm","mhc_rank_pct","mhc_raw_score","not_evaluated_codes","exclusion_reasons","downgrade_reasons","warning_codes","default_include",
  "origin_id","origin_type","origin_validity","origin_warning_codes","origin_evidence_status","source_event_ids","context_event_ids","germline_context_ids","source_event_mapping_status","report_source_event_refs","report_context_event_refs","report_candidate_ids","compound_source_context","source_domain","source_locus",
  "transcript_id","transcript_version","protein_id","protein_version","gene_symbol","gene_symbols","consequence","translation_status","canonical_status","mane_status","isoform_source","phase_status","phase_namespace","haplotype_id","sequence_certainty","sequence_context_id","peptide_start_aa","peptide_end_aa","wt_peptide_sequence","wt_comparison_status","normal_comparison_status","wt_definition","mutation_positions_json","junction_positions_json","upstream_flank","downstream_flank","source_allele",
  "event_ref","event_domain","event_id","event_role","event_type","contig","pos","end","ref","alt","record_status","link_role","is_source","source_field","source_set_id","source_set_event_refs","independent_causation_claimed",
  "source","analysis_id","run_id","prediction_record_id","provider","predictor_version","model_version","task","metric","metric_kind","metric_semantics","unit","score_direction","prediction_scope","model_input_mhc_key","record_kind","raw_score","original_score_text","reason","evidence_id","status",
  "evidence_category","evidence_status","B_subtypes","A_experimental_status","B1_experimental_status","B2_experimental_status","reference_endpoint","reference_assay","experimental_reference_id","experimental_reference_version",
  "tumor_vaf","normal_vaf","rna_vaf","tumor_alt_count","normal_alt_count","rna_alt_count","tumor_evidence_status","normal_evidence_status","rna_evidence_status","gene_tpm","transcript_tpm","ccf",
  "native_priority_origin_id","best_candidate_id","best_origin_id","scope_id","n_distinct_peptides","n_candidate_peptides","n_pmhc_candidates","n_origins","n_candidate_origin_links","n_context_candidates","candidate_ids","context_origin_ids","source_event_refs","context_event_refs","alternative_HLAs","synthesized_peptide_unit","planning_status","primary_HLA_proxy","primary_candidate_id","score_bin","assumed_or_empirical_rate","probability_interpretation","planning_profile_id","planning_profile_version")
html_report_share <- function(tables, subjects, metadata=list()) {
  maps <- new.env(parent=emptyenv())
  alias <- function(value,family) {
    if(is.na(value))return(NA_character_)
    if(!exists(family,maps,inherits=FALSE))assign(family,character(),maps)
    m<-get(family,maps);at<-match(value,names(m))
    if(is.na(at)) { m<-c(m,stats::setNames(paste0(family,"-",sprintf("%04d",length(m)+1L)),value));assign(family,m,maps);at<-length(m) }
    unname(m[[at]])
  }
  scrub <- function(x,key="") {
    if(is.list(x)) {
      if(is.null(names(x)))return(lapply(x,scrub,key=key))
      return(stats::setNames(lapply(names(x),function(k)scrub(x[[k]],k)),names(x)))
    }
    if(!is.character(x))return(x)
    family<-html_report_id_family(key)
    if(!is.null(family))return(vapply(x,alias,"",family=family,USE.NAMES=FALSE))
    # Sequence strings are scientific values, even if a subject used that string as an ID.
    sequence<-grepl("peptide_sequence$|ligand_sequence$|flank$|^ref$|^alt$",key)
    if(sequence)return(x)
    vapply(x,function(s) {
      if(is.na(s))return(NA_character_)
      s<-gsub("file://[^[:space:]\"']+","[path withheld]",s,ignore.case=TRUE,perl=TRUE)
      s<-gsub("(^|[[:space:]\"'])(/[[:alnum:]_.~][^[:space:]\"']*|[A-Za-z]:[\\\\/][^[:space:]\"']*|\\\\\\\\[^[:space:]\"']+)","\\1[path withheld]",s,perl=TRUE)
      # Typed scientific values are never masked because a subject happened to
      # use the same word/number as an ID. Paths remain withheld even here.
      if(key%in%c(setdiff(html_report_fields(),"reason"),"code")||grepl("sequence|^mhc_key$|^hla_|provider|predictor|model_|metric|^unit$|score|^task$|status|grade|rank|version|gene|transcript|protein|codes",key))return(s)
      if(s %in% subjects)return("[identifier withheld]")
      for(id in subjects[nchar(subjects)>=3L])s<-gsub(id,"[identifier withheld]",s,fixed=TRUE)
      s
    },"",USE.NAMES=FALSE)
  }
  out<-list();withheld<-list()
  for(n in names(tables)) {
    d<-tables[[n]]
    if(!is.data.frame(d))next
    fields<-intersect(names(d),html_report_fields())
    rows<-report_rows(d[fields]);out[[n]]<-lapply(rows,scrub)
    withheld[[n]]<-length(setdiff(names(d),fields))
  }
  list(tables=out,withheld_field_counts=withheld,metadata=scrub(metadata),private_id_map=unlist(lapply(ls(maps),function(f) lapply(seq_along(get(f,maps)),function(i)list(family=f,original=names(get(f,maps))[[i]],alias=unname(get(f,maps)[[i]])))),recursive=FALSE))
}
html_report_validate <- function(t) {
  ids<-t$pmhc_candidates$candidate_id;oi<-t$origins$origin_id;er<-t$events$event_ref
  if(anyDuplicated(ids)||anyDuplicated(oi)||anyDuplicated(er))stop("HTML requires unique saved record identities")
  links<-t$candidate_origin_links;el<-t$origin_event_links
  if(any(!links$candidate_id%in%ids)||any(!links$origin_id%in%oi)||any(!el$origin_id%in%oi)||any(!el$event_ref%in%er))stop("HTML linkage target missing; no record was omitted")
  for(c in report_rows(t$pmhc_candidates))if(!setequal(report_chars(c$origin_ids),links$origin_id[links$candidate_id==c$candidate_id]))stop("HTML candidate-Origin associations differ from reporting table")
  invisible(TRUE)
}
html_report_payload <- function(bundle, project_metadata=list(), privacy="share", selected=FALSE) {
  t<-bundle$tables;html_report_validate(t)
  cc<-t$pmhc_candidates
  nn<-if("native_public_result"%in%names(cc))cc$native_public_result else rep(list(NULL),nrow(cc))
  scalar<-function(x)if(is.null(x)||!length(x))NA_character_ else as.character(x[[1L]])
  cc$html_native_prediction_status<-vapply(seq_len(nrow(cc)),function(i)scalar(if(is.null(nn[[i]]$prediction_status))cc$prediction_status[[i]] else nn[[i]]$prediction_status),"")
  cc$html_native_grade_status<-vapply(seq_len(nrow(cc)),function(i)scalar(if(!is.null(nn[[i]]$grade_status))nn[[i]]$grade_status else if(!is.null(nn[[i]]$grading_status))nn[[i]]$grading_status else cc$grading_status[[i]]),"")
  cc$html_native_reasons<-I(lapply(seq_len(nrow(cc)),function(i)list(not_evaluated_codes=nn[[i]]$not_evaluated_codes,grade_reason_codes=nn[[i]]$grade_reason_codes)))
  cc$html_native_status_source<-vapply(nn,function(n)if(is.null(n$prediction_status))"candidate snapshot (public status not supplied)" else "native_public_result.prediction_status","")
  t$pmhc_candidates<-cc
  # A bounded candidate view uses the exact report links, never new molecular inference.
  if(selected) {
    refs<-unique(t$origin_event_links$event_ref)
    for(n in c("events","event_summary","event_best_by_channel"))t[[n]]<-t[[n]][t[[n]]$event_ref%in%refs,,drop=FALSE]
  }
  counts<-lapply(t,nrow)
  counts$unique_peptide_sequences<-length(unique(t$peptide_summary$peptide_sequence))
  subjects<-unique(c(report_chars(lapply(report_rows(t$pmhc_candidates),function(x)x[c("sample_id","patient_id")])),project_metadata$subjects))
  meta<-list(contract="offline-report-view-1.0.0",reporting_contract=bundle$manifest$contract,
    package_version=bundle$manifest$package_version,source_manifest_sha256=bundle$manifest$source_manifest$sha256,
    original_project_graph_counts=bundle$manifest$original_project_graph_counts,
    output_scope=if(selected)"EXPLICIT_CANDIDATE_SELECTION; all linked Origins and associated Events retained" else "ALL_REPORT_RECORDS",
    scoped_table_rows=counts,reference_versions=project_metadata$references,input_versions=project_metadata$inputs,
    privacy=privacy,label=project_metadata$label,source_notices=bundle$manifest$warnings,probability_note="No probability or new score is calculated by this viewer.")
  if(privacy=="share") {
    share<-html_report_share(t,subjects,meta);data<-share$tables;meta<-share$metadata
    meta$withheld_field_counts<-share$withheld_field_counts
    # Safe version/checksum metadata only. Private run paths/options remain in full mode.
    meta$privacy_note<-"Sharing view: identifiers pseudonymized; paths and non-allowlisted/free-form metadata withheld. All selected Origin rows and links remain. Biological sequences, mutations and annotations can still identify a person: review before sharing. Original identifiers are not embedded."
  } else {data<-lapply(t,report_rows);meta$original_manifest<-bundle$manifest;meta$privacy_note<-"PRIVATE FULL DATA: identifiers, absolute paths and original nested records may be present. Do not share without review."}
  # Same alias dictionary must cover project version identifiers too; project metadata
  # intentionally carries no paths or personal identifiers, only role/version/hash.
  out<-list(metadata=meta,tables=data)
  if(privacy=="share")attr(out,"private_id_map")<-share$private_id_map
  out
}
html_report_json <- function(x) {
  s<-report_json(x)
  for(k in c("<",">","&","\u2028","\u2029"))s<-gsub(k,switch(k,"<"="\\u003c",">"="\\u003e","&"="\\u0026","\u2028"="\\u2028","\u2029"="\\u2029"),s,fixed=TRUE)
  s
}
html_report_write <- function(project, file, analyses=list(), candidate_id=NULL,
    privacy="share", confirm_private=FALSE, max_candidates=500L,
    max_bytes=10*1024^2, report_options=list(), label="Saved scientific results") {
  if(!is.character(privacy)||length(privacy)!=1L||is.na(privacy)||!privacy%in%c("share","private"))stop("privacy must be share or private")
  setup_flag(confirm_private,"confirm_private")
  if(privacy=="private"&&!confirm_private)stop("Private HTML contains identifiers and paths; set confirm_private=TRUE only after reviewing sharing risks")
  if(privacy=="private")warning("PRIVATE HTML contains embedded scientific data, identifiers and possibly paths; review before sharing",call.=FALSE)
  if(!is.numeric(max_candidates)||length(max_candidates)!=1L||!is.finite(max_candidates)||max_candidates<1||max_candidates>5000||max_candidates!=floor(max_candidates))stop("max_candidates must be an integer from 1 to 5000")
  if(!is.numeric(max_bytes)||length(max_bytes)!=1L||!is.finite(max_bytes)||max_bytes<1024||max_bytes>50*1024^2)stop("max_bytes must be between 1024 and 50 MiB")
  yield_text(label,"HTML label")
  strict_keys(report_options,c("yield_plan","external_channels","reference","endpoint","assay"),"HTML report_options")
  native<-neo_results(project,candidate_id=candidate_id)
  if(nrow(native)>max_candidates)stop("Candidate count exceeds explicit HTML limit; select exact candidate_id values or raise max_candidates deliberately. No candidate was silently truncated.")
  b<-do.call(neo_report_bundle,c(list(project=project,analyses=analyses,candidate_id=candidate_id),report_options))
  if(nrow(b$tables$pmhc_candidates)>max_candidates)stop("Reporting candidate count exceeds HTML limit; no candidate was silently truncated")
  version_record<-function(r)list(role=if(is.null(r$input_role))r$logical_name else r$input_role,resource_id=r$resource_id,version=r$version,integrity=r$integrity[c("mode","observed","declared","verified")])
  metadata<-list(label=label,subjects=unique(c(report_chars(lapply(project$samples,`[[`,"sample_id")),report_chars(lapply(project$samples,`[[`,"patient_id")))),
    references=unname(lapply(project$resources,version_record)),inputs=unname(lapply(project$inputs,version_record)))
  payload<-html_report_payload(b,metadata,privacy,!is.null(candidate_id))
  json<-html_report_json(payload)
  root<-system.file("report-html",package="Neoantimon2")
  css<-paste(readLines(file.path(root,"viewer.css"),warn=FALSE),collapse="\n")
  js<-paste(readLines(file.path(root,"viewer.js"),warn=FALSE),collapse="\n")
  nonce<-digest::digest(paste(json,css,js),algo="sha256",serialize=FALSE)
  template<-paste(readLines(file.path(root,"viewer.html"),warn=FALSE),collapse="\n")
  replace<-function(text,token,value) {
    # literal replacement without gsub replacement-string interpretation
    at<-regexpr(token,text,fixed=TRUE);if(at[[1L]]<0L)stop("Missing HTML template token")
    paste0(substr(text,1L,at[[1L]]-1L),value,substring(text,at[[1L]]+nchar(token)))
  }
  for(pair in list(c("{{NONCE}}",nonce),c("{{CSS}}",css),c("{{SCRIPT}}",js))) {
    while(grepl(pair[[1L]],template,fixed=TRUE))template<-replace(template,pair[[1L]],pair[[2L]])
  }
  template<-replace(template,"{{DATA}}",json)
  if(nchar(template,type="bytes")>max_bytes)stop("HTML byte limit exceeded; select fewer candidates or raise max_bytes deliberately. All Origins are retained; no HTML was written.")
  yield_text(file,"HTML file");target<-resolve_path(file)
  if(file.exists(target))stop("HTML target already exists; use a fresh file")
  forbidden<-c(project$path,project$config$paths$repository,vapply(project$inputs,`[[`,"","resolved_path"),vapply(Filter(function(r)!is.null(r$resolved_path),project$resources),`[[`,"","resolved_path"))
  forbidden<-forbidden[!is.na(forbidden)&nzchar(forbidden)]
  if(any(vapply(forbidden,function(p)within_path(target,resolve_path(p))||within_path(resolve_path(p),target),TRUE)))stop("HTML output must not overlap source project/repository/inputs/references")
  if(!dir.exists(dirname(target)))stop("Create the HTML parent directory explicitly before writing")
  tmp<-tempfile(".html-",tmpdir=dirname(target));on.exit(unlink(tmp),add=TRUE)
  writeLines(template,tmp,useBytes=TRUE);Sys.chmod(tmp,"0600")
  if(!file.link(tmp,target))stop("Cannot publish fresh HTML exclusively; choose a local filesystem supporting file links. Existing targets are never overwritten.")
  invisible(list(path=target,sha256=external_sha(target),privacy=privacy,scope=payload$metadata$output_scope,counts=payload$metadata$scoped_table_rows,private_id_map=attr(payload,"private_id_map")))
}
