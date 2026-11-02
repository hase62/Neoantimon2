# Optional reporting projection only. Biological/scoring identity is never changed.
report_version <- function() "result-report-0.1.0"
report_chars <- function(x) unique(as.character(unlist(x, use.names = FALSE)))
# Explicit 17-digit finite-double encoding supports round-trip TSV/list JSON.
# The original upstream decimal score string is retained separately, unchanged.
report_json <- function(x) {
  if(is.null(x))return("null")
  if(is.data.frame(x))return(paste0("[",paste(vapply(report_rows(x),report_json,""),collapse=","),"]"))
  if(is.list(x)) {
    if(!length(x))return(if(is.null(names(x)))"[]" else "{}")
    values<-vapply(x,report_json,"")
    if(is.null(names(x)))return(paste0("[",paste(values,collapse=","),"]"))
    if(anyNA(names(x))||any(!nzchar(names(x)))||anyDuplicated(names(x)))stop("JSON object keys must be unique and nonempty")
    keys<-vapply(names(x),function(k)as.character(jsonlite::toJSON(k,auto_unbox=TRUE)),"")
    return(paste0("{",paste(paste0(keys,":",values),collapse=","),"}"))
  }
  if(is.object(x)&&!inherits(x,"AsIs"))return(as.character(jsonlite::toJSON(x,auto_unbox=TRUE,null="null",na="null",digits=NA,force=TRUE)))
  values<-if(is.numeric(x))vapply(x,function(v)if(is.na(v)||!is.finite(v))"null" else sprintf("%.17g",v),"") else
    if(is.logical(x))ifelse(is.na(x),"null",ifelse(x,"true","false")) else
      vapply(x,function(v)if(is.na(v))"null" else as.character(jsonlite::toJSON(v,auto_unbox=TRUE)),"")
  if(length(values)==1L)return(unname(values))
  paste0("[",paste(values,collapse=","),"]")
}
report_hash <- function(x) digest::digest(report_json(x),algo="sha256",serialize=FALSE)
report_rows <- function(x) {
  if (is.null(x)) return(list())
  if (!is.data.frame(x)) return(x)
  lapply(seq_len(nrow(x)), function(i) lapply(x, function(col)
    if (is.list(col)) col[[i]] else col[i]))
}
report_frame <- function(rows, required = character()) {
  fields <- sort(unique(c(required, unlist(lapply(rows, names), use.names = FALSE))), method = "radix")
  cols <- lapply(fields, function(f) {
    values <- lapply(rows, `[[`, f); present <- Filter(Negate(is.null), values)
    if (any(vapply(present, is.list, TRUE)) || any(lengths(present) != 1L)) return(I(values))
    numeric <- length(present) && all(vapply(present, is.numeric, TRUE))
    logical <- length(present) && all(vapply(present, is.logical, TRUE))
    vapply(values, function(v) if (is.null(v)) {
      if (numeric) NA_real_ else if (logical) NA else NA_character_
    } else if (numeric) as.numeric(v) else if (logical) v else as.character(v),
      if (numeric) 0 else if (logical) TRUE else "")
  }); names(cols) <- fields
  as.data.frame(cols, check.names = FALSE, stringsAsFactors = FALSE)
}
report_with_raw <- function(r, raw = r) {
  r$source_record_json <- report_json(raw); r
}
report_index <- function(rows, field) {
  if(!length(rows))return(list())
  split(rows,vapply(rows,`[[`,"",field))
}
report_gather <- function(index, keys) unname(unlist(index[keys],recursive=FALSE))
report_key <- function(...) paste(..., sep = "\037")
report_event_ref <- function(domain, id) paste0("report-event:", report_hash(list(domain = domain, event_id = id)))
report_peptide_ref <- function(sample, sequence) paste0("report-peptide:",
  report_hash(list(sample_id = sample, peptide_sequence = sequence)))
report_scalar <- function(x, default = NA_character_) if (is.null(x) || !length(x)) default else x

report_catalog <- function(graph) {
  catalog <- list()
  add <- function(domain, id, raw, sample = NULL, available = TRUE) {
    if (!is.character(id) || length(id) != 1L || !nzchar(id)) stop("Reporting Event identity missing")
    ref <- report_event_ref(domain, id)
    node <- list(event_ref = ref, event_domain = domain, event_id = id,
      sample_id = report_scalar(if (is.null(raw$sample_id)) sample else raw$sample_id),
      event_role = report_scalar(raw$event_role, if (domain == "core") "unknown" else "structural_source"),
      record_status = if (available) "PRESERVED_SOURCE_RECORD" else "SOURCE_EVENT_RECORD_UNAVAILABLE",
      event_type = report_scalar(raw$event_type, report_scalar(raw$source_type,
        if (domain == "long_read") "long_read_isoform" else NA_character_)),
      contig = report_scalar(raw$contig), pos = report_scalar(raw$pos, NA_real_),
      end = report_scalar(raw$end, NA_real_), ref = report_scalar(raw$ref), alt = report_scalar(raw$alt),
      original_location = if (!is.null(raw$structure)) raw$structure else raw[intersect(c("contig", "pos", "end", "ref", "alt", "mate_contig", "mate_pos"),names(raw))],
      raw = raw)
    if (!is.null(catalog[[ref]]) && !identical(catalog[[ref]], node)) stop("Conflicting Event records in the same namespace")
    catalog[[ref]] <<- node; ref
  }
  for (e in graph$events) add("core", e$event_id, e)
  source <- external_annotation_records(graph, "source", "canonical_events")
  lr <- external_annotation_records(graph, "long_read", "events")
  for (e in source) add("source", e$source_event_id, e, e$structure$sample_id)
  for (e in lr) add("long_read", e$long_read_source_event_id, e, e$context$sample_id)
  mappings <- list()
  for (family in c("source", "long_read")) {
    field <- if (family == "source") "source_event_id" else "long_read_source_event_id"
    for (p in external_annotation_records(graph, family, "peptides"))
      if (!is.null(p$origin_id) && !is.null(p[[field]])) mappings[[length(mappings) + 1L]] <-
        list(origin_id = p$origin_id, domain = family, id = p[[field]], field = paste(family, "peptides", field, sep = "."))
  }
  links <- list(); unlinked <- character()
  for (o in graph$origins) {
    own <- list()
    put <- function(domain, id, field, context = FALSE) {
      ref <- report_event_ref(domain, id)
      if (is.null(catalog[[ref]])) add(domain, id, list(), o$sample_id, available = FALSE)
      e <- catalog[[ref]]
      is_source <- !context && (domain != "core" || identical(e$event_role, "somatic") ||
        identical(e$record_status, "SOURCE_EVENT_RECORD_UNAVAILABLE"))
      if (is_source && e$record_status == "PRESERVED_SOURCE_RECORD" &&
          !is.na(e$sample_id) && !identical(e$sample_id, o$sample_id))
        stop("Causative Event/Origin sample mismatch")
      own[[length(own) + 1L]] <<- list(origin_id = o$origin_id, sample_id = o$sample_id,
        event_ref = ref, event_domain = domain, event_id = id, source_field = field,
        link_role = if (context) "germline_context" else if (is_source) "generation_source_set_member" else "non_somatic_context",
        is_source = is_source, source_record_status = e$record_status,
        source_event_sample_id = e$sample_id)
    }
    for (id in report_chars(o$source_event_ids)) put("core", id, "source_event_ids")
    for (id in report_chars(o$germline_context_ids)) put("core", id, "germline_context_ids", TRUE)
    ext <- Filter(function(m) identical(m$origin_id, o$origin_id), mappings)
    for (m in ext) put(m$domain, m$id, m$field)
    if (!length(ext) && !is.null(o$source_antigen_key)) {
      # Exact registered identity only. Never infer a namespace from an ID prefix.
      at <- Filter(function(e) e$event_domain != "core" && identical(e$event_id, o$source_antigen_key), catalog)
      if (length(at) == 1L) put(at[[1L]]$event_domain, at[[1L]]$event_id, "source_antigen_key_exact_registered")
      else if (length(at) > 1L) stop("Ambiguous structural Event namespace; explicit origin mapping is required")
      else put("unresolved_structural", o$source_antigen_key, "source_antigen_key_unresolved")
    }
    if (!length(own)) unlinked <- c(unlinked, o$origin_id)
    if (length(own)) {
      keys <- vapply(own, function(l) report_key(l$event_ref, l$link_role), "")
      own <- own[!duplicated(keys)]
      sources <- sort(unique(vapply(Filter(function(l) l$is_source, own), `[[`, "", "event_ref")), method = "radix")
      set_id <- if (length(sources)) paste0("report-source-set:", report_hash(as.list(sources))) else NA_character_
      own <- lapply(own, function(l) {
        l$source_set_id <- set_id; l$source_set_event_refs <- as.list(sources)
        l$compound_source_context <- length(sources) > 1L
        l$independent_causation_claimed <- FALSE; l
      })
      links <- c(links, own)
    }
  }
  list(events = unname(catalog), links = links, origins_without_declared_events = unlinked)
}

report_prediction_fields <- function() c("source", "analysis_id", "prediction_scope", "provider", "predictor_version", "model_version",
  "run_id", "task", "metric", "metric_kind", "metric_semantics", "unit", "score_direction")
report_prediction_rows <- function(graph, runs, external) {
  ids <- vapply(graph$candidates, `[[`, "", "candidate_id"); rows <- list()
  native <- function(e, run, kind) {
    list(source = "Native", analysis_id = "Native", prediction_scope = "peptide_hla", candidate_id = e$candidate_id,
      provider = report_scalar(e$predictor, report_scalar(run$provider)),
      predictor_version = report_scalar(e$predictor_version, report_scalar(run$predictor_version, "")),
      model_version = report_scalar(e$model_version, report_scalar(run$model_version, "")),
      run_id = report_scalar(e$run_id, report_scalar(run$run_id)), task = report_scalar(e$task, report_scalar(run$task)),
      metric = report_scalar(e$metric), metric_kind = if (identical(e$unit, "percent")) "percentile_rank" else "raw_score",
      metric_semantics = if (identical(e$unit, "percent")) if (endsWith(report_scalar(e$task, ""), "_binding")) "mhc_binding_percentile" else
        if (endsWith(report_scalar(e$task, ""), "_presentation")) "mhc_presentation_percentile" else "raw" else "raw",
      raw_score = if (is.null(e$raw_score)) NA_real_ else e$raw_score,
      original_score_text = NA_character_, unit = report_scalar(e$unit), score_direction = report_scalar(e$score_direction),
      prediction_status = report_scalar(e$status, "unavailable"), reason = report_scalar(e$reason, ""),
      record_kind = kind, evidence_id = report_scalar(e$evidence_id), source_record = e)
  }
  for (run in runs) {
    for (e in run$evidence) if (e$candidate_id %in% ids)
      rows[[length(rows) + 1L]] <- native(e, run, "measurement")
    for (e in run$statuses) if (!is.null(e$candidate_id) && e$candidate_id %in% ids)
      rows[[length(rows) + 1L]] <- native(e, run, "provider_status")
  }
  for (a in external) for (e in a$payload$evidence) if (e$candidate_id %in% ids) {
    if (identical(e$status, "evaluated") && (is.null(e$numeric_score) || !is.finite(e$numeric_score)))
      stop("Evaluated external observation lacks a finite score")
    if (any(!report_chars(e$origin_ids) %in% report_chars(graph$candidates[[match(e$candidate_id, ids)]]$origin_ids)))
      stop("External observation has an Origin outside its candidate")
    rows[[length(rows) + 1L]] <- list(source = "External", analysis_id = a$manifest$analysis_id,
      candidate_id = e$candidate_id, prediction_scope = e$scope, model_input_mhc_key = e$mhc_key, provider = e$provider, predictor_version = e$predictor_version,
      model_version = e$model_version, run_id = e$run_id, task = e$task, metric = e$metric,
      metric_kind = e$metric_kind, metric_semantics = e$metric_semantics,
      raw_score = if (is.null(e$numeric_score)) NA_real_ else e$numeric_score,
      original_score_text = e$score, unit = e$unit, score_direction = e$score_direction,
      prediction_status = e$status, reason = e$reason, record_kind = "measurement",
      evidence_id = e$evidence_id, source_record = e)
  }
  # These are explicit coverage records, not invented measurements or negatives.
  for (id in ids) {
    if (!any(vapply(rows, function(r) identical(r$source, "Native") && identical(r$candidate_id, id), TRUE)))
      rows[[length(rows) + 1L]] <- list(source = "Native", analysis_id = "Native", candidate_id = id,
        prediction_status = "NOT_RECORDED", record_kind = "coverage_status", raw_score = NA_real_,
        reason = "No persisted Native prediction record for this candidate")
    for (a in external) if (!any(vapply(rows, function(r) identical(r$analysis_id, a$manifest$analysis_id) && identical(r$candidate_id, id), TRUE)))
      rows[[length(rows) + 1L]] <- list(source = "External", analysis_id = a$manifest$analysis_id,
        candidate_id = id, prediction_status = "NOT_PREDICTED", record_kind = "coverage_status", raw_score = NA_real_,
        reason = "No matched observation in explicitly selected analysis")
  }
  rows
}

report_best <- function(events, event_links, graph, native, predictions, selectors) {
  candidates <- stats::setNames(graph$candidates, vapply(graph$candidates, `[[`, "", "candidate_id"))
  origins <- stats::setNames(graph$origins, vapply(graph$origins, `[[`, "", "origin_id"))
  nr <- stats::setNames(native, vapply(native, `[[`, "", "candidate_id"))
  # Native policy ranks are consumed, never reconstructed from Grade labels.
  native_scopes <- unique(vapply(graph$candidates, function(c) report_key(c$sample_id, c$hla_class,
    report_scalar(nr[[c$candidate_id]]$policy_version), report_scalar(c$analysis_context)), ""))
  external_rows <- Filter(function(p) p$source == "External" && p$record_kind == "measurement", predictions)
  eligible_selector <- function(p) {
    if (is.null(selectors)) return(TRUE)
    any(vapply(selectors, function(s) all(vapply(names(s), function(f)
      identical(p[[f]], s[[f]]), TRUE)), TRUE))
  }
  if (!is.null(selectors)) {
    if (!is.list(selectors)) stop("external_channels must be a list of named selectors")
    for (s in selectors) {
      if (!is.list(s) || !length(s) || is.null(names(s)) || any(!names(s) %in% report_prediction_fields()))
        stop("Invalid external channel selector")
      for (f in names(s)) yield_text(s[[f]], paste("External selector", f))
      if (!any(vapply(external_rows, function(p) all(vapply(names(s), function(f) identical(p[[f]], s[[f]]), TRUE)), TRUE)))
        stop("Requested external representative channel is not in the selected analyses")
    }
  }
  external_rows <- Filter(eligible_selector, external_rows)
  scopes <- list()
  for (key in native_scopes) {
    cs <- Filter(function(c) identical(report_key(c$sample_id,c$hla_class,report_scalar(nr[[c$candidate_id]]$policy_version),report_scalar(c$analysis_context)),key),graph$candidates)
    first <- cs[[1L]];scopes[[length(scopes)+1L]] <- list(source="Native",sample_id=first$sample_id,hla_class=first$hla_class,
      analysis_context=report_scalar(first$analysis_context),policy_version=report_scalar(nr[[first$candidate_id]]$policy_version),
      candidate_ids=vapply(cs,`[[`,"","candidate_id"),scope_id=paste0("report-rank-scope:",report_hash(list(source="Native",key=key))))
  }
  fields <- report_prediction_fields()
  keys <- vapply(external_rows,function(p) report_json(p[fields]),"")
  for (key in unique(keys)) {
    ps <- external_rows[keys==key]
    for (sample in unique(vapply(ps,function(p)candidates[[p$candidate_id]]$sample_id,"")))
      for (cl in unique(vapply(ps,function(p)candidates[[p$candidate_id]]$hla_class,""))) {
        own <- Filter(function(p) identical(candidates[[p$candidate_id]]$sample_id,sample) && identical(candidates[[p$candidate_id]]$hla_class,cl),ps)
        if(!length(own))next
        sc <- ps[[1L]][fields];sc$sample_id<-sample;sc$hla_class<-cl
        sc$candidate_ids<-unique(vapply(own,`[[`,"","candidate_id"));sc$predictions<-own
        sc$scope_id<-paste0("report-rank-scope:",report_hash(sc[c(fields,"sample_id","hla_class")]))
        scopes[[length(scopes)+1L]]<-sc
      }
  }
  best<-list()
  event_index<-report_index(event_links,"event_ref");origin_link_index<-report_index(graph$links,"origin_id")
  source_origin_index<-report_index(Filter(function(l)isTRUE(l$is_source),event_links),"origin_id")
  scopes<-lapply(scopes,function(sc) {
    if(!is.null(sc$predictions))sc$prediction_index<-report_index(sc$predictions,"candidate_id")
    sc
  })
  scope_of<-list()
  for(sc in scopes)for(id in sc$candidate_ids)scope_of[[id]]<-c(scope_of[[id]],sc$scope_id)
  for(e in events) {
    el<-Filter(function(l)isTRUE(l$is_source),event_index[[e$event_ref]])
    oids<-unique(vapply(el,`[[`,"","origin_id"))
    linked<-report_gather(origin_link_index,oids)
    cids<-unique(vapply(linked,`[[`,"","candidate_id"))
    own_scopes<-Filter(function(sc)identical(sc$sample_id,e$sample_id),scopes)
    if(!length(own_scopes))own_scopes<-list(list(source="Native",sample_id=e$sample_id,hla_class=NA_character_,scope_id="no_candidate_scope",candidate_ids=character()))
    for(sc in own_scopes) {
      ids<-cids[vapply(cids,function(id)sc$scope_id %in% scope_of[[id]],TRUE)];available<-character();pick<-NULL
      if(sc$source=="Native") {
        available<-Filter(function(id) {
          r<-nr[[id]]; !is.null(r$rank_within_sample_class)&&is.finite(r$rank_within_sample_class)&&
            identical(r$rank_status,"RANKED")&&identical(r$primary_prediction_status,"evaluated")
        },as.list(ids));available<-unlist(available,use.names=FALSE)
        if(length(available))pick<-available[order(vapply(available,function(id)nr[[id]]$rank_within_sample_class,0),available,method="radix")][1L]
      } else {
        ps<-Filter(function(p)identical(p$prediction_status,"evaluated")&&is.finite(p$raw_score),report_gather(sc$prediction_index,ids))
        if(length(ps)) {
          direction<-sc$score_direction;if(!direction %in% c("higher_is_better","lower_is_better"))stop("Unknown external score direction")
          byid<-split(ps,vapply(ps,`[[`,"","candidate_id"))
          if(any(vapply(byid,function(z)length(unique(vapply(z,`[[`,0,"raw_score")))>1L,TRUE)))stop("Conflicting scores within exact analysis/channel")
          available<-names(byid);scores<-vapply(byid,function(z)z[[1L]]$raw_score,0)
          pick<-available[order(if(direction=="higher_is_better")-scores else scores,available,method="radix")][1L]
        }
      }
      if(e$record_status!="PRESERVED_SOURCE_RECORD")pick<-NULL
      r<-c(e[c("sample_id","event_ref","event_domain","event_id","event_type","contig","pos","end","ref","alt")],
        sc[setdiff(names(sc),c("candidate_ids","predictions","prediction_index","sample_id"))])
      r$associated_candidate_ids<-as.list(sort(ids,method="radix"));r$associated_origin_ids<-as.list(sort(oids,method="radix"))
      r$selection_policy<-if(sc$source=="Native")"existing_native_rank_within_sample_class_policy_context_v1" else "exact_analysis_model_run_task_metric_unit_direction_then_candidate_id_v1"
      r$status<-if(!is.null(pick))"SELECTED" else "UNAVAILABLE"
      r$reason<-if(e$record_status!="PRESERVED_SOURCE_RECORD")"SOURCE_EVENT_RECORD_UNAVAILABLE" else if(!length(cids))"NO_ASSOCIATED_SOURCE_CANDIDATE_IN_EXPORT_SCOPE" else if(!length(ids))"NO_ASSOCIATED_CANDIDATE_IN_THIS_SCOPE" else "NO_EVALUABLE_RANKED_CANDIDATE_IN_THIS_SCOPE"
      r$best_candidate_id<-if(is.null(pick))NA_character_ else pick
      if(!is.null(pick)) {
        c<-candidates[[pick]];n<-nr[[pick]];own_o<-sort(intersect(oids,report_chars(c$origin_ids)),method="radix")
        if(!length(own_o))stop("Representative has no Event-associated Origin")
        preferred<-n$priority_origin_id
        oid<-if(!is.null(preferred)&&!is.na(preferred)&&preferred %in% own_o)preferred else own_o[1L]
        o<-origins[[oid]]
        r$best_peptide_sequence<-c$peptide_sequence;r$best_peptide_length<-nchar(c$peptide_sequence);r$best_HLA<-c$mhc_key;r$best_origin_id<-oid
        r$gene_symbol<-o$gene_symbol;r$genes<-o$genes_json;r$transcript_id<-o$transcript_id;r$transcript_version<-o$transcript_version;r$protein_id<-o$protein_id
        r$native_grade<-n$final_grade;r$native_rank<-n$rank_within_sample_class;r$native_rank_status<-n$rank_status
        r$native_policy_version<-n$policy_version;r$native_priority_origin_id<-n$priority_origin_id
        r$best_origin_is_native_priority_origin<-!is.null(preferred)&&!is.na(preferred)&&identical(oid,preferred)
        r$native_grade_scope<-"persisted_candidate_global_grade; not recomputed for the Event-specific Origin"
        r$origin_selection_policy<-"native_priority_origin_if_Event_associated_else_lexical_associated_origin_v1"
        ol<-source_origin_index[[oid]]
        r$joint_source_event_refs<-as.list(sort(unique(vapply(ol,`[[`,"","event_ref")),method="radix"))
        r$compound_source_context<-length(r$joint_source_event_refs)>1L;r$independent_causation_claimed<-FALSE
        if(sc$source=="External") {
          selected<-Filter(function(p)identical(p$prediction_status,"evaluated")&&is.finite(p$raw_score),sc$prediction_index[[pick]])[[1L]]
          r$external_score<-selected$raw_score;r$external_original_score_text<-selected$original_score_text;r$prediction_status<-selected$prediction_status
        } else r$prediction_status<-n$primary_prediction_status
        r$reason<-"Representative only; every alternative remains exported"
      }
      best[[length(best)+1L]]<-r
    }
  }
  best
}

report_construct <- function(graph, native, runs = list(), external = list(), raw_evidence = NULL,
    yield_plan = NULL, candidate_id = NULL, external_channels = NULL, reference = NULL,
    endpoint = NULL, assay = NULL) {
  original_counts <- list(events = length(graph$events), origins = length(graph$origins),
    candidates = length(graph$candidates), candidate_origin_links = length(graph$links))
  cids <- vapply(graph$candidates, `[[`, "", "candidate_id")
  if (anyDuplicated(cids)) stop("Duplicate candidate identity")
  if (!is.null(candidate_id)) {
    validate_filter(candidate_id, "candidate_id")
    if (any(!candidate_id %in% cids)) stop("Requested candidate_id is absent from saved graph")
    graph$candidates <- Filter(function(c) c$candidate_id %in% candidate_id, graph$candidates)
    graph$links <- Filter(function(l) l$candidate_id %in% candidate_id, graph$links)
    oids <- unique(vapply(graph$links, `[[`, "", "origin_id"))
    graph$origins <- Filter(function(o) o$origin_id %in% oids, graph$origins)
  }
  cids <- vapply(graph$candidates, `[[`, "", "candidate_id")
  origins <- stats::setNames(graph$origins, vapply(graph$origins, `[[`, "", "origin_id"))
  candidates <- stats::setNames(graph$candidates, cids)
  for (l in graph$links) if (!l$candidate_id %in% cids || !l$origin_id %in% names(origins)) stop("Unresolved candidate/Origin linkage")
  if (anyDuplicated(names(origins))) stop("Duplicate Origin identity")
  candidate_link_index<-report_index(graph$links,"candidate_id")
  origin_link_index<-report_index(graph$links,"origin_id")
  for (c in candidates) {
    own <- candidate_link_index[[c$candidate_id]]
    actual <- vapply(own, `[[`, "", "origin_id")
    if (anyDuplicated(actual) || !setequal(actual, report_chars(c$origin_ids))) stop("Candidate/Origin linkage mismatch")
    for (id in actual) {
      o <- origins[[id]]
      if (is.null(o) || !identical(o$sample_id,c$sample_id) || !identical(o$peptide_sequence,c$peptide_sequence))
        stop("Candidate/Origin sample or peptide identity mismatch")
    }
  }
  catalog <- report_catalog(graph); event_links <- catalog$links
  event_link_index<-report_index(event_links,"event_ref");origin_event_index<-report_index(event_links,"origin_id")
  native <- report_rows(native); native <- Filter(function(n) n$candidate_id %in% cids, native)
  if (anyDuplicated(vapply(native, `[[`, "", "candidate_id"))) stop("Native result identity ambiguous")
  ni <- stats::setNames(native, vapply(native, `[[`, "", "candidate_id"))
  predictions <- report_prediction_rows(graph, runs, external)
  warnings <- character()
  if (!is.null(yield_plan)) {
    if (!inherits(yield_plan,"neo_experimental_yield")) stop("yield_plan must be an existing neo_experimental_yield result")
    yp <- report_rows(yield_plan$pair_evidence)
    for (p in yp) {
      ids <- report_chars(p$all_candidate_ids)
      if (!length(ids) || any(!ids %in% cids)) stop("Yield candidate lies outside export scope")
      for (id in ids) {
        c <- candidates[[id]]
        if (!identical(c$sample_id,p$sample_id) || !identical(c$peptide_sequence,p$peptide_sequence) || !identical(c$mhc_key,p$mhc_key))
          stop("Yield peptide/HLA/sample identity mismatch")
        if (!setequal(report_chars(c$origin_ids),report_chars(p$all_origin_ids))) stop("Yield Origin provenance is incomplete or mismatched")
      }
    }
    if (is.null(reference) && isTRUE(yield_plan$reference$available)) {
      a <- yield_plan$reference$artifact
      if (!is.null(a$path) && file.exists(a$path)) {
        if (!identical(external_sha(a$path),a$sha256)) stop("Yield experimental reference changed")
        reference <- a$path
      } else {
        reference <- yield_plan$reference$data
        reference$records <- lapply(reference$records,function(r) {
          if(!is.null(r$source_mhc_key))r$mhc_key<-r$source_mhc_key
          if(!is.null(r$source_peptide_sequence))r$peptide_sequence<-r$source_peptide_sequence
          r
        })
      }
    }
    if (is.null(endpoint)) endpoint <- yield_plan$evidence_policy$endpoint
    if (is.null(assay)) assay <- yield_plan$evidence_policy$assay
  }
  pairs <- report_frame(lapply(candidates,function(c)c[c("candidate_id","sample_id","peptide_sequence","mhc_key")]),
    c("candidate_id","sample_id","peptide_sequence","mhc_key"))
  if(nrow(pairs)) evidence <- neo_yield_evidence(pairs,reference,endpoint,assay,
    min_records=if(is.null(yield_plan))100L else yield_plan$evidence_policy$min_records,
    min_patients=if(is.null(yield_plan))3L else yield_plan$evidence_policy$min_patients) else {
    evidence <- pairs; evidence$evidence_category <- character();evidence$evidence_status<-character()
    evidence$b_subtypes<-I(list());evidence$matches<-I(list())
  }
  if (!is.null(yield_plan) && isTRUE(yield_plan$reference$available) && isTRUE(attr(evidence,"reference")$available)) {
    er<-attr(evidence,"reference")$data;yr<-yield_plan$reference$data
    if(!identical(er$id,yr$id)||!identical(er$version,yr$version)||!identical(report_hash(er$records),report_hash(yr$records))||!identical(endpoint,yield_plan$evidence_policy$endpoint)||!identical(assay,yield_plan$evidence_policy$assay))
      stop("Report reference/endpoint/assay differs from the supplied yield plan")
    a<-attr(evidence,"reference")$artifact;b<-yield_plan$reference$artifact
    if(!is.null(a$path)&&!is.null(b$path)&&!identical(a$sha256,b$sha256))stop("Reference content differs from supplied yield plan")
  }
  evrows <- report_rows(evidence);ei<-stats::setNames(evrows,vapply(evrows,`[[`,"","candidate_id"))
  candidate_rows <- lapply(candidates,function(c) {
    r<-c;n<-ni[[c$candidate_id]];e<-ei[[c$candidate_id]]
    r$peptide_unit_id<-report_peptide_ref(c$sample_id,c$peptide_sequence)
    r$native_grade<-n$final_grade;r$native_rank<-n$rank_within_sample_class;r$native_rank_status<-n$rank_status
    r$native_policy_version<-n$policy_version;r$native_priority_origin_id<-n$priority_origin_id
    r$native_public_result<-n
    r$evidence_category<-e$evidence_category;r$evidence_status<-e$evidence_status;r$B_subtypes<-e$b_subtypes
    r$A_experimental_status<-e$A_status;r$B1_experimental_status<-e$B1_status;r$B2_experimental_status<-e$B2_status
    r$experimental_matches<-e$matches;r$hla_tested_units<-e$hla_tested_units;r$hla_verified_patients<-e$hla_verified_patients
    r$reference_endpoint<-endpoint;r$reference_assay<-assay
    rd<-attr(evidence,"reference")$data
    r$experimental_reference_id<-rd$id;r$experimental_reference_version<-rd$version;r$experimental_reference_provenance<-rd$provenance
    report_with_raw(r,c)
  })
  ci<-stats::setNames(candidate_rows,cids)
  predictions <- lapply(seq_along(predictions),function(i) {
    p<-predictions[[i]];c<-ci[[p$candidate_id]]
    p$prediction_record_id<-paste0("report-prediction:",report_hash(list(index=i,record=p)))
    p$sample_id<-c$sample_id;p$peptide_sequence<-c$peptide_sequence;p$peptide_length<-nchar(c$peptide_sequence)
    p$mhc_key<-c$mhc_key;p$hla_class<-c$hla_class;p$peptide_unit_id<-c$peptide_unit_id;p$origin_ids<-c$origin_ids
    for(f in c("native_grade","native_rank","native_rank_status","native_policy_version","evidence_category","evidence_status",
      "B_subtypes","A_experimental_status","B1_experimental_status","B2_experimental_status","experimental_reference_id","experimental_reference_version","reference_endpoint","reference_assay"))p[f]<-c[f]
    report_with_raw(p,p$source_record)
  })
  best<-report_best(catalog$events,event_links,graph,native,predictions,external_channels)
  event_rows<-lapply(catalog$events,function(e) report_with_raw(c(e[setdiff(names(e),"raw")],list(original_event_record=e$raw)),e$raw))
  best_event_index<-report_index(best,"event_ref")
  event_summary<-lapply(catalog$events,function(e) {
    ol<-event_link_index[[e$event_ref]]
    so<-unique(vapply(Filter(function(l)isTRUE(l$is_source),ol),`[[`,"","origin_id"))
    co<-unique(vapply(Filter(function(l)!isTRUE(l$is_source),ol),`[[`,"","origin_id"))
    ls<-report_gather(origin_link_index,so);lc<-report_gather(origin_link_index,co)
    cs<-sort(unique(vapply(ls,`[[`,"","candidate_id")),method="radix")
    ps<-unique(vapply(origins[so],`[[`,"","peptide_sequence"))
    reps<-best_event_index[[e$event_ref]]
    c(e[setdiff(names(e),"raw")],list(n_distinct_peptides=length(ps),n_candidate_peptides=length(unique(vapply(candidates[cs],`[[`,"","peptide_sequence"))),
      n_pmhc_candidates=length(cs),n_origins=length(so),n_candidate_origin_links=length(ls),n_context_candidates=length(unique(vapply(lc,`[[`,"","candidate_id"))),
      candidate_ids=as.list(cs),origin_ids=as.list(sort(so,method="radix")),context_origin_ids=as.list(co),
      compound_origin_count=length(unique(vapply(Filter(function(l)l$origin_id %in% so&&isTRUE(l$compound_source_context),ol),`[[`,"","origin_id"))),
      n_selected_representative_scopes=sum(vapply(reps,function(r)identical(r$status,"SELECTED"),TRUE)),
      representative_scope_ids=as.list(vapply(reps,`[[`,"","scope_id")),
      status=if(e$record_status!="PRESERVED_SOURCE_RECORD")"SOURCE_RECORD_UNAVAILABLE" else if(length(cs))"ASSOCIATED_CANDIDATES_RETAINED" else if(length(lc))"CONTEXT_ONLY_NOT_CAUSATIVE" else "NO_CANDIDATE_IN_EXPORT_SCOPE",
      source_association_policy="declared joint generation source set; not independent per-Event causation"))
  })
  origin_rows<-lapply(origins,function(o) {
    r<-o;ol<-origin_event_index[[o$origin_id]]
    r$peptide_unit_id<-report_peptide_ref(o$sample_id,o$peptide_sequence)
    r$report_source_event_refs<-as.list(vapply(Filter(function(l)isTRUE(l$is_source),ol),`[[`,"","event_ref"))
    r$report_context_event_refs<-as.list(vapply(Filter(function(l)!isTRUE(l$is_source),ol),`[[`,"","event_ref"))
    r$source_event_mapping_status<-if(length(ol))"DECLARED_LINKS_RETAINED" else "NO_DECLARED_SOURCE_EVENT"
    r$compound_source_context<-length(r$report_source_event_refs)>1L
    r$report_candidate_ids<-as.list(vapply(origin_link_index[[o$origin_id]],`[[`,"","candidate_id"))
    r$molecular_evidence_join<-"source_event_ids -> molecular_event_evidence.event_id; measurements are not collapsed across a compound Origin"
    report_with_raw(r,o)
  })
  groups<-split(origin_rows,vapply(origin_rows,function(o)report_key(o$sample_id,o$peptide_sequence),""))
  for(c in candidates) {
    key<-report_key(c$sample_id,c$peptide_sequence)
    if(!key %in% names(groups))groups[[key]]<-list()
  }
  peptide_rows<-lapply(names(groups),function(key) {
    oo<-groups[[key]]
    first<-if(length(oo))oo[[1L]] else {
      c<-Filter(function(c)identical(report_key(c$sample_id,c$peptide_sequence),key),candidates)[[1L]]
      list(sample_id=c$sample_id,peptide_sequence=c$peptide_sequence,peptide_unit_id=report_peptide_ref(c$sample_id,c$peptide_sequence))
    };oids<-vapply(oo,`[[`,"","origin_id")
    cs<-Filter(function(c)identical(c$sample_id,first$sample_id)&&identical(c$peptide_sequence,first$peptide_sequence),candidates)
    r<-list(peptide_unit_id=first$peptide_unit_id,sample_id=first$sample_id,peptide_sequence=first$peptide_sequence,
      peptide_length=nchar(first$peptide_sequence),synthesized_peptide_unit="unique_sequence_per_sample",
      candidate_ids=as.list(sort(vapply(cs,`[[`,"","candidate_id"),method="radix")),origin_ids=as.list(sort(oids,method="radix")),
      alternative_HLAs=as.list(sort(unique(vapply(cs,`[[`,"","mhc_key")),method="radix")),
      source_event_refs=as.list(sort(unique(report_chars(lapply(oo,`[[`,"report_source_event_refs"))),method="radix")),
      context_event_refs=as.list(sort(unique(report_chars(lapply(oo,`[[`,"report_context_event_refs"))),method="radix")),
      gene_symbols=as.list(sort(unique(report_chars(lapply(oo,`[[`,"gene_symbol"))),method="radix")),
      transcript_associations=lapply(oo,function(o)o[intersect(c("origin_id","transcript_id","transcript_version","protein_id","genes_json"),names(o))]),
      planning_status="NOT_REQUESTED_FOR_THIS_PEPTIDE",primary_HLA_proxy=NA_character_,assumed_or_empirical_rate=NA_real_)
    if(!is.null(yield_plan)) {
      y<-yield_plan$units[yield_plan$units$sample_id==first$sample_id&yield_plan$units$peptide_sequence==first$peptide_sequence,,drop=FALSE]
      if(nrow(y)==1L&&identical(yield_plan$summary$unit,"peptide")) {
        r$planning_status<-"SINGLE_HLA_PROXY_FROM_UNCHANGED_YIELD_PLAN";r$primary_HLA_proxy<-y$mhc_key[[1L]]
        r$primary_candidate_id<-y$candidate_id[[1L]];r$score_bin<-y$score_bin[[1L]]
        r$assumed_or_empirical_rate<-y$probability[[1L]];r$probability_interpretation<-y$probability_status[[1L]]
        r$evidence_category<-y$evidence_category[[1L]];r$B_subtypes<-y$b_subtypes[[1L]]
        r$planning_profile_id<-yield_plan$profile$data$id;r$planning_profile_version<-yield_plan$profile$data$version
        r$profile_provenance<-yield_plan$profile$data$provenance;r$profile_artifact<-yield_plan$profile$artifact
        r$primary_HLA_rule<-y$primary_HLA_rule[[1L]]
      } else if(nrow(y))r$planning_status<-"ASSAY_UNIT_PLAN_SEE_UNCHANGED_YIELD_TABLES"
    }
    r
  })
  prediction_index<-report_index(predictions,"candidate_id")
  prediction_links<-list()
  if(!is.null(yield_plan)) for(p in report_rows(yield_plan$pair_evidence)) {
    ids<-report_chars(p$all_candidate_ids)
    for(id in ids) {
      same<-Filter(function(e) identical(e$candidate_id,id)&&identical(e$record_kind,"measurement")&&
        all(vapply(setdiff(yield_channel_fields(),"source"),function(f) {
          target<-if(f=="raw_score")p$raw_score else p[[f]]
          identical(e[[f]],target)
        },TRUE))&&identical(e$source,p$source),prediction_index[[id]])
      evaluated<-isTRUE(p$yield_evaluated)
      if(length(same)&&evaluated&&!any(vapply(same,function(e)identical(as.numeric(e$raw_score),as.numeric(p$raw_score)),TRUE)))
        stop("Yield score differs from persisted selected model/run/channel")
      if(!length(same)&&evaluated)warnings<-c(warnings,"YIELD_PREDICTION_NOT_BOUND_TO_SELECTED_PERSISTED_ANALYSES")
      for(e in if(length(same))same else list(NULL))prediction_links[[length(prediction_links)+1L]]<-list(candidate_id=id,
        prediction_record_id=if(is.null(e))NA_character_ else e$prediction_record_id,
        analysis_id=if(is.null(e))NA_character_ else e$analysis_id,
        status=if(is.null(e))"SOURCE_CHANNEL_NOT_IN_SELECTED_ANALYSES" else "PERSISTED_CHANNEL_AND_SCORE_VERIFIED")
    }
  }
  tables<-list(events=report_frame(event_rows,c("event_ref","event_domain","event_id","sample_id")),
    event_summary=report_frame(event_summary,c("event_ref","event_id","n_pmhc_candidates")),
    event_best_by_channel=report_frame(best,c("event_ref","event_id","scope_id","source","status","best_candidate_id")),
    peptide_summary=report_frame(peptide_rows,c("peptide_unit_id","sample_id","peptide_sequence")),
    pmhc_candidates=report_frame(candidate_rows,c("candidate_id","sample_id","peptide_sequence","mhc_key")),
    pmhc_predictions=report_frame(predictions,c("prediction_record_id","candidate_id","source","analysis_id","raw_score","prediction_status")),
    origins=report_frame(origin_rows,c("origin_id","sample_id","peptide_sequence")),
    candidate_origin_links=report_frame(graph$links,c("candidate_id","origin_id","sample_id")),
    origin_event_links=report_frame(event_links,c("origin_id","event_ref","event_domain","event_id","link_role","is_source")),
    candidate_peptide_links=report_frame(lapply(candidates,function(c)list(candidate_id=c$candidate_id,sample_id=c$sample_id,
      peptide_unit_id=report_peptide_ref(c$sample_id,c$peptide_sequence))),c("candidate_id","sample_id","peptide_unit_id")),
    molecular_event_evidence=report_frame(lapply(raw_evidence$event_evidence,report_with_raw),"event_id"),
    expression_evidence=report_frame(lapply(raw_evidence$expression_evidence,report_with_raw),c("sample_id","transcript_id")),
    source_observations=report_frame(lapply(external_annotation_records(graph,"source","observations"),report_with_raw),"source_event_id"),
    long_read_observations=report_frame(lapply(external_annotation_records(graph,"long_read","observations"),report_with_raw),"long_read_source_event_id"))
  tables$native_run_statuses<-report_frame(unlist(lapply(runs,function(run) lapply(run$statuses,function(s) {
    r<-s;r$run_id<-run$run_id;r$provider<-run$provider;report_with_raw(r,s)
  })),recursive=FALSE),c("run_id","provider","candidate_id","status","reason"))
  if(!is.null(yield_plan)) {
    for(n in c("units","pair_evidence","curve","targets","summary"))tables[[paste0("experimental_yield_",n)]]<-yield_plan[[n]]
    tables$experimental_yield_prediction_links<-report_frame(prediction_links,c("candidate_id","prediction_record_id","analysis_id","status"))
  }
  list(tables=tables,metadata=list(original_project_graph_counts=original_counts,
    scope=list(candidate_ids=if(is.null(candidate_id))"ALL" else as.list(sort(unique(candidate_id),method="radix")),
      Event_policy="all original Event nodes retained, with export-scope counts",Origin_policy=if(is.null(candidate_id))"all Origins including orphans" else "all Origins linked to selected candidates"),
    warnings=as.list(unique(warnings)),unmapped_origin_ids=as.list(catalog$origins_without_declared_events),
    reference=if(is.null(attr(evidence,"reference")))list(status="UNAVAILABLE") else {
      rr<-attr(evidence,"reference");if(!is.null(rr$data))rr$data$records<-NULL;rr
    },
    yield=if(is.null(yield_plan))list(status="NOT_REQUESTED") else list(status="SUPPLIED_UNCHANGED",
      profile=yield_plan$profile,evidence_policy=yield_plan$evidence_policy,assumptions=yield_plan$assumptions,warnings=yield_plan$warnings,
      source_digest=digest::digest(yield_plan,algo="sha256")),
    policies=list(native="existing evaluated Native rank; separate sample/class/policy/context scopes; no regrading",
      external="separate analysis/provider/version/model/run/task/metric/kind/semantics/unit/direction/sample/class scopes; raw score direction then candidate ID",
      origin="native priority Origin if Event-associated, else lexical associated Origin; global candidate Grade is never assigned as an Event-specific Origin Grade",
      compound="joint declared generation source sets; never claim independent per-Event causation",list_cells="JSON with null/empty arrays preserved; source_record_json is authoritative original record")))
}

report_write_table <- function(x, path) {
  enc <- lapply(x,function(col) {
    if(is.list(col))return(vapply(col,report_json,""))
    if(is.numeric(col))return(vapply(col,function(v)if(is.na(v))NA_character_ else
      sprintf("%.17g",v),""))
    if(is.logical(col))return(ifelse(is.na(col),NA_character_,ifelse(col,"true","false")))
    as.character(col)
  })
  flat<-as.data.frame(enc,check.names=FALSE,stringsAsFactors=FALSE)
  utils::write.table(flat,path,sep="\t",quote=TRUE,qmethod="double",row.names=FALSE,na="")
  list(file=basename(path),rows=nrow(x),sha256=external_sha(path),bytes=unname(file.info(path)$size),
    columns=lapply(names(x),function(f)list(name=f,encoding=if(is.list(x[[f]]))"JSON" else
      if(is.numeric(x[[f]]))"JSON_number" else if(is.logical(x[[f]]))"boolean" else "string")))
}

#' Complete read-only Event/peptide/pMHC/Origin reporting
#' @param project A completed local neo_project. Ordinary multi-sample projects
#'   and core/source/long-read namespaces are supported. For longitudinal aggregate
#'   handles, report the underlying sample projects individually.
#' @param dir New directory outside all source project/repository/input paths;
#'   NULL returns the reporting tables in memory only.
#' @param analyses Explicit list of existing external analysis handles or IDs.
#'   No analyses are discovered, imported or executed automatically.
#' @param candidate_id Optional exact candidate IDs for a bounded report. All
#'   original Events remain visible; Origins/links are bounded to these candidates.
#' @param yield_plan Optional existing neo_experimental_yield object, retained
#'   unchanged, including all alternative HLA pairs.
#' @param external_channels Optional list of named exact channel selectors for
#'   representatives only. NULL reports each distinct selected-analysis channel
#'   separately. Long predictions always retain every selected-analysis record.
#' @param reference,endpoint,assay Optional explicit experimental reference and
#'   exact evidence scope, or reuse the reference already supplied in yield_plan.
#' @return A neo_report_bundle list with tables, versioned manifest and optional path.
#' @export
neo_report_bundle <- function(project, dir = NULL, analyses = list(), candidate_id = NULL,
    yield_plan = NULL, external_channels = NULL, reference = NULL, endpoint = NULL, assay = NULL) {
  if(longitudinal_aggregate(project))stop("Report underlying longitudinal sample projects individually; aggregate export is not supported in result-report-0.1.0")
  graph<-external_graph(project)
  if(!is.list(analyses))stop("analyses must be an explicit list of existing analysis handles or IDs")
  external<-lapply(analyses,function(a)external_analysis_read(a,project))
  ids<-vapply(external,function(a)a$manifest$analysis_id,"")
  if(anyDuplicated(ids))stop("Select each external analysis once")
  if(length(ids))external<-external[order(ids,method="radix")]
  native<-neo_results(project)
  runs<-somatic_prediction_runs(project)
  raw<-input_evidence_read(project)
  if(is.null(raw))raw<-somatic_read_sidecar(project)$raw
  out<-report_construct(graph,native,runs,external,raw,yield_plan,candidate_id,
    external_channels,reference,endpoint,assay)
  manifest<-c(list(contract=report_version(),schema_version="0.1.0",
    package_version=as.character(utils::packageVersion("Neoantimon2")),
    project_run_id=project$run_id,project_analysis_id=project$analysis_id,
    source_manifest=list(path=manifest_file(project$path),sha256=external_sha(manifest_file(project$path))),
    external_analyses=lapply(external,function(a)list(analysis_id=a$manifest$analysis_id,
      manifest=a$manifest,path=a$path,policy=a$payload$policy)),
    native_run_metadata=lapply(runs,function(r)r[intersect(c("run_id","provider","task","predictor_version","model_version","model_identity","capability"),names(r))]),
    explicit_external_selectors=external_channels),out$metadata)
  for(n in names(out$tables)) {
    d<-out$tables[[n]]
    if(nrow(d)&&ncol(d)) {
      fields<-intersect(c("sample_id","event_domain","event_id","scope_id","peptide_sequence","candidate_id","origin_id","analysis_id","prediction_record_id"),names(d))
      if(length(fields)&&all(!vapply(d[fields],is.list,TRUE))) {
        ord<-do.call(order,c(lapply(d[fields],as.character),list(method="radix",na.last=TRUE)))
        # Yield curve/target/unit order is part of the plan and never sorted away.
        if(!startsWith(n,"experimental_yield_"))out$tables[[n]]<-d[ord,,drop=FALSE]
      }
    }
  }
  if(is.null(dir))return(structure(list(tables=out$tables,manifest=manifest,path=NULL),class="neo_report_bundle"))
  yield_text(dir,"Reporting target directory")
  target<-resolve_path(dir)
  forbidden<-c(project$path,project$config$paths$repository,
    vapply(project$inputs,`[[`,"","resolved_path"),
    vapply(Filter(function(r)!is.null(r$resolved_path),project$resources),`[[`,"","resolved_path"))
  forbidden<-forbidden[!is.na(forbidden)&nzchar(forbidden)]
  if(any(vapply(forbidden,function(p)within_path(target,resolve_path(p))||within_path(resolve_path(p),target),TRUE)))
    stop("Reporting output must not overlap source project/repository/inputs/references")
  if(file.exists(target))stop("Reporting target already exists; use a fresh directory")
  dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
  tmp<-tempfile(".report-",tmpdir=dirname(target));if(!dir.create(tmp))stop("Cannot create report staging directory")
  on.exit(unlink(tmp,recursive=TRUE),add=TRUE)
  manifest$tables<-lapply(names(out$tables),function(n)report_write_table(out$tables[[n]],file.path(tmp,paste0(n,".tsv"))))
  manifest$report_id<-paste0("report-bundle:",report_hash(manifest))
  writeLines(report_json(manifest),file.path(tmp,"manifest.json"),useBytes=TRUE)
  jsonlite::read_json(file.path(tmp,"manifest.json"))
  for(t in manifest$tables)if(!identical(external_sha(file.path(tmp,t$file)),t$sha256))stop("Reporting integrity failure")
  if(!file.rename(tmp,target))stop("Cannot store report atomically")
  structure(list(path=target,tables=out$tables,manifest=manifest),class="neo_report_bundle")
}
