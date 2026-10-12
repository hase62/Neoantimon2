somatic_eligibility_reason <- function(candidate_id,origin_id,event_id,code,field,
                                       observed,threshold,operator,hard) {
  list(candidate_id=candidate_id,origin_id=origin_id,event_id=event_id,
    reason_code=code,evidence_field=field,observed_value=observed,
    policy_threshold=threshold,policy_operator=operator,
    policy_version=somatic_eligibility_version(),hard_failure=hard)
}

somatic_origin_eligibility <- function(path,candidate_id,event_index,policy) {
  reasons <- list(); criteria <- list(); dna_criteria <- list(); failures <- FALSE; missing <- FALSE
  ids <- unlist(path$somatic_event_ids)
  if (!length(ids)) {
    missing <- TRUE
    reasons[[1L]] <- somatic_eligibility_reason(candidate_id,path$origin_id,NULL,
      "no_somatic_source_event","source_event_ids",NULL,NULL,NULL,FALSE)
  }
  for (id in ids) {
    event <- event_index[[id]]
    if (is.null(event)) stop("Somatic origin event evidence missing")
    for (assay in c("tumor","normal")) {
      x <- event[[assay]]
      field <- paste0(assay,"_alt_reads")
      threshold <- if (assay=="tumor") policy$tumor_min_alt_reads else
        policy$normal_max_alt_reads
      operator <- if (assay=="tumor") ">=" else "<="
      count <- x$alt_reads
      satisfies <- if (is.null(count) || x$status %in% c("invalid_counts",
        "unavailable","insufficient_depth")) NULL else
        if (assay=="tumor") count>=threshold else count<=threshold
      dna_criteria[[length(dna_criteria)+1L]] <- list(event_id=id,
        evidence_field=field,observed_value=count,policy_threshold=threshold,
        policy_operator=operator,satisfied=satisfies,assessment_status=x$status,
        policy_version=somatic_eligibility_version())
      if (x$status=="invalid_counts") {
        failures <- TRUE
        reasons[[length(reasons)+1L]] <- somatic_eligibility_reason(candidate_id,
          path$origin_id,id,"invalid_dna_counts",field,count,threshold,operator,TRUE)
      } else if (is.null(count) || x$status %in% c("unavailable","insufficient_depth")) {
        missing <- TRUE
        reasons[[length(reasons)+1L]] <- somatic_eligibility_reason(candidate_id,
          path$origin_id,id,paste0(assay,"_alt_not_evaluated"),field,count,
          threshold,operator,FALSE)
      } else if ((assay=="tumor" && count<threshold) ||
                 (assay=="normal" && count>threshold)) {
        failures <- TRUE
        reasons[[length(reasons)+1L]] <- somatic_eligibility_reason(candidate_id,
          path$origin_id,id,if (assay=="tumor") "insufficient_tumor_alt_reads" else
            "excessive_normal_alt_reads",field,count,threshold,operator,TRUE)
      }
    }
    rna <- event$rna
    criteria[[length(criteria)+1L]] <- list(event_id=id,evidence_field="rna_alt_reads",
      observed_value=rna$alt_reads,policy_threshold=policy$rna_min_alt_reads,
      policy_operator=">=",satisfied=if (is.null(rna$alt_reads) ||
      rna$status %in% c("unavailable","insufficient_depth",
          "invalid_counts")) NULL else rna$alt_reads>=policy$rna_min_alt_reads,
      rna_status=rna$status,policy_version=somatic_eligibility_version())
  }
  list(origin_id=path$origin_id,status=if (failures) "ineligible" else
    if (missing) "not_evaluated" else "eligible",reasons=reasons,
    dna_criteria=dna_criteria,rna_criteria=criteria)
}

somatic_eligibility_construct <- function(raw,policy) {
  event_index <- stats::setNames(raw$event_evidence,
    vapply(raw$event_evidence,`[[`,"","event_id"))
  records <- lapply(raw$candidate_links,function(link) {
    paths <- lapply(link$origin_paths,somatic_origin_eligibility,
      candidate_id=link$candidate_id,event_index=event_index,policy=policy)
    statuses <- vapply(paths,`[[`,"","status")
    status <- if ("eligible" %in% statuses) "eligible" else
      if (length(statuses) && all(statuses=="ineligible")) "ineligible" else "not_evaluated"
    list(candidate_id=link$candidate_id,eligibility_status=status,
      reason_records=unlist(lapply(paths,`[[`,"reasons"),recursive=FALSE),
      origin_results=paths,policy_version=somatic_eligibility_version())
  })
  list(schema_version=somatic_evidence_version(),policy_version=somatic_eligibility_version(),
    policy=policy,candidate_eligibility=records)
}

somatic_validate_raw <- function(raw,events,peptide) {
  if (!identical(raw$schema_version,somatic_evidence_version()))
    stop("Unsupported somatic evidence schema version")
  ids <- vapply(events,`[[`,"","event_id")
  if (!setequal(vapply(raw$event_evidence,`[[`,"","event_id"),ids))
    stop("Somatic event evidence identity mismatch")
  candidates <- vapply(peptide$candidates,`[[`,"","candidate_id")
  if (!setequal(vapply(raw$candidate_links,`[[`,"","candidate_id"),candidates))
    stop("Somatic candidate evidence identity mismatch")
  event_evidence_ids <- vapply(raw$event_evidence,`[[`,"","evidence_id")
  for (link in raw$candidate_links) for (path in link$origin_paths)
    if (any(!unlist(path$event_evidence_ids) %in% event_evidence_ids))
      stop("Candidate source evidence link missing")
  invisible(TRUE)
}

somatic_validate_eligibility <- function(result,raw) {
  if (!identical(result$schema_version,somatic_evidence_version()) ||
      !identical(result$policy_version,somatic_eligibility_version()))
    stop("Unsupported somatic eligibility version")
  if (!setequal(vapply(result$candidate_eligibility,`[[`,"","candidate_id"),
      vapply(raw$candidate_links,`[[`,"","candidate_id")))
    stop("Somatic eligibility candidate identity mismatch")
  for (record in result$candidate_eligibility) {
    if (!record$eligibility_status %in% c("eligible","ineligible","not_evaluated"))
      stop("Invalid somatic eligibility status")
    for (reason in record$reason_records) {
      if (!identical(reason$candidate_id,record$candidate_id) ||
          !identical(reason$policy_version,result$policy_version))
        stop("Somatic eligibility reason identity mismatch")
    }
  }
  invisible(TRUE)
}
