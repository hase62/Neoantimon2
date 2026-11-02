# Additive source-file evidence. No alteration of somatic grade/rank policies.
input_evidence_version <- function() "source-evidence-v2-postreview"

input_vcf_assay <- function(event,allow_conflict=FALSE) {
  if (!identical(event$event_role,"somatic")) return(NULL)
  observations <- list()
  for (record in event$event_details_json$source_records) {
    if (!identical(record$format,"VCF") || is.null(record$raw_format)) next
    fields <- event_vcf_format(record$raw_format)
    n <- length(record$original$alts); i <- record$selected_alt_index
    # Only the bound single tumor sample's FORMAT fields; never INFO/AF.
    parse <- function(key,number) {
      x <- fields[[key]]
      if (is.null(x) || x==".") return(NULL)
      v <- event_vcf_allele_values(x,number,n,paste("tumor FORMAT",key))
      if (is.null(v)) return(NULL)
      lapply(v,function(x) if(x==".") NA_integer_ else somatic_integer(x,paste("tumor FORMAT",key)))
    }
    ad <- parse("AD","R"); ao <- parse("AO","A")
    ro <- if (is.null(fields$RO) || fields$RO==".") NULL else
      somatic_integer(fields$RO,"tumor FORMAT RO")
    dp <- if (is.null(fields$DP) || fields$DP==".") NULL else
      somatic_integer(fields$DP,"tumor FORMAT DP")
    counted <- NULL; method <- NULL
    if (!is.null(ad)) { counted <- unlist(ad); method <- "FORMAT/AD" }
    if ((is.null(counted) || anyNA(counted)) && !is.null(ao) && !anyNA(unlist(ao)) && !is.null(ro)) {
      counted <- c(ro,unlist(ao)); method <- "FORMAT/AO+RO"
    }
    reported <- NULL
    # A declared AF must actually describe sample allele fraction (not a
    # population frequency). Undeclared/ambiguous AF is retained as raw data.
    definition <- record$evidence_format_definitions$AF
    if (!is.null(fields$AF) && fields$AF!="." && !is.null(definition) &&
        identical(definition$Number,"A") &&
        grepl("allele.*(fraction|frequency)",definition$Description,ignore.case=TRUE) &&
        !grepl("population|gnomad|1000.genome",definition$Description,ignore.case=TRUE)) {
      values <- event_vcf_allele_values(fields$AF,"A",n,"tumor FORMAT AF")
      if (!is.null(values) && values[[i]]!=".")
        reported <- somatic_number(values[[i]],"tumor FORMAT AF",0,1)
    }
    if (is.null(counted) && is.null(reported) && is.null(dp)) next
    if (!is.null(counted) && sum(counted,na.rm=TRUE)>.Machine$integer.max)
      stop("Tumor VCF allele depths exceed supported integer range: ",event$event_id)
    # DP and AD can use different filters (e.g. GATK). Preserve a mismatch,
    # but never replace the AD denominator by DP or reject valid caller data.
    depth_discrepancy <- !is.null(counted) && !is.null(dp) && sum(counted,na.rm=TRUE)>dp
    denominator <- if (is.null(counted) || anyNA(counted)) NULL else as.integer(sum(counted))
    count_at <- function(j) if(is.null(counted) || is.na(counted[[j]])) NULL else as.integer(counted[[j]])
    value <- if (!is.null(denominator) && denominator>0) counted[[i+1L]]/denominator else NULL
    if (!is.null(reported) && ((!is.null(denominator) && denominator==0L) || (!is.null(dp) && dp==0L)))
      stop("Tumor FORMAT AF conflicts with zero allele depth")
    conflict <- !is.null(value) && !is.null(reported) && abs(value-reported)>1e-8
    obs <- list(ref_reads=count_at(1L),
      alt_reads=count_at(i+1L),
      depth=dp,usable_depth=denominator,other_reads=if (is.null(counted) || anyNA(counted)) NULL else
        as.integer(sum(counted)-counted[[1L]]-counted[[i+1L]]),
      ambiguous_reads=NULL,vaf_denominator=denominator,
      vaf=if (!is.null(value)) value else reported,reported_vaf=reported,
      vaf_status=if (!is.null(value)) "derived" else if (!is.null(reported)) "supplied" else "unavailable",
      count_unit=if (is.null(counted) && is.null(dp)) NULL else "read",
      status=if ((!is.null(denominator) && denominator==0L) || (is.null(counted) && !is.null(dp) && dp==0L)) "insufficient_depth" else
        if (is.null(count_at(i+1L))) {if(!is.null(reported)) "fraction_only" else "not_evaluated"} else
        if (counted[[i+1L]]>0) "observed_support" else "observed_no_alt",
      provenance=list(source="tumor_vcf",input_id=record$input_id,line=record$line,
        selected_alt_index=i,method=if (is.null(method)) if(is.null(reported)) "FORMAT/DP_only" else "FORMAT/AF" else method,
        denominator_policy="sum_all_allele_depths",raw_format=record$raw_format,
        counted_vaf=value,caller_estimated_vaf=reported,discrepancy=conflict,depth_discrepancy=depth_discrepancy))
    observations[[length(observations)+1L]] <- obs
  }
  if (!length(observations)) return(NULL)
  conflict <- FALSE
  for (field in c("ref_reads","alt_reads","vaf_denominator","vaf")) {
    known <- unlist(lapply(observations,`[[`,field))
    if (length(known)>1L && any(abs(known-known[[1L]])>1e-8)) conflict <- TRUE
  }
  # Missing fields do not contradict a more complete source observation.
  strength <- vapply(observations,function(x) if(!is.null(x$vaf))
    if(!is.null(x$alt_reads)) 3L else 2L else if(!is.null(x$alt_reads)) 1L else 0L,0L)
  result <- observations[[which.max(strength)]]
  result$provenance$source_observation_conflict <- conflict
  depths <- unlist(lapply(observations,`[[`,"depth"))
  result$provenance$depth_source_disagreement <- length(unique(depths))>1L
  result$provenance$source_records <- lapply(observations,`[[`,"provenance")
  if (conflict && !allow_conflict) {
    # Two valid caller observations must not stop otherwise supported peptide
    # reconstruction, nor be pooled or silently chosen as a measurement.
    for(field in c("ref_reads","alt_reads","depth","usable_depth","other_reads",
      "ambiguous_reads","vaf_denominator","vaf","reported_vaf")) result[field] <- list(NULL)
    result$status <- "ambiguous";result$vaf_status <- "ambiguous"
  }
  result
}

input_assay <- function(event,row,assay,source) {
  explicit <- somatic_assay(event,row,assay,source)
  # Allele measurements are selected as a coherent bundle. Depth alone is an
  # independent observation, not evidence that ALT fraction was measured.
  complete <- !is.null(explicit$vaf) ||
    (!is.null(explicit$ref_reads) && !is.null(explicit$alt_reads)) ||
    identical(explicit$status,"invalid_counts")
  vcf <- if (assay=="tumor") input_vcf_assay(event,allow_conflict=complete) else NULL
  bam <- if (is.null(source)) NULL else source$bam_observations[[event$event_id]][[assay]]
  selected <- if (complete) explicit else if (!is.null(vcf$vaf)) vcf else
    if (!is.null(bam)) bam else if (!is.null(vcf)) vcf else explicit
  if (!is.null(vcf) || !is.null(bam)) {
    known <- function(x) any(vapply(x[c("ref_reads","alt_reads","depth","vaf")],
      function(v)!is.null(v),TRUE))
    alternatives <- list(vcf=vcf,bam=bam)
    if (known(explicit)) alternatives$explicit <- explicit
    differs <- function(x,y,fields) any(vapply(fields,function(f)
      !is.null(x[[f]]) && !is.null(y[[f]]) && abs(x[[f]]-y[[f]])>1e-8,TRUE))
    comparable <- function(x,y) c("vaf",if(!is.null(x$count_unit) &&
      identical(x$count_unit,y$count_unit)) c("ref_reads","alt_reads","depth"))
    conflict <- !complete && differs(explicit,selected,intersect(
      comparable(explicit,selected),c("ref_reads","alt_reads","vaf")))
    selected$measurement_sources <- list(alleles=selected$provenance,depth=selected$provenance)
    # Attach an independently sourced depth only when the selected bundle lacks
    # it and its unit agrees. Never use this depth to derive or replace VAF.
    if (!complete && is.null(selected$depth) && !is.null(explicit$depth) &&
        identical(explicit$count_unit,selected$count_unit)) {
      selected$depth <- explicit$depth
      selected$measurement_sources$depth <- explicit$provenance
    }
    selected$depth_count_disagreement <- !is.null(selected$depth) &&
      !is.null(selected$vaf_denominator) && selected$vaf_denominator>selected$depth
    selected$alternatives <- alternatives
    selected$source_precedence <- "complete_explicit_alleles > tumor_VCF > BAM_fallback; independent_depth"
    selected$source_disagreement <- isTRUE(vcf$provenance$source_observation_conflict) ||
      selected$depth_count_disagreement ||
      any(vapply(alternatives,function(x)differs(x,selected,comparable(x,selected)),TRUE))
    if (conflict) {
      selected$ambiguous_allele_observation <- selected[c("ref_reads","alt_reads","vaf_denominator","vaf")]
      selected$ref_reads <- NULL;selected$alt_reads <- NULL
      selected$vaf <- NULL;selected$vaf_status <- "ambiguous";selected$status <- "ambiguous"
    }
  }
  selected
}

input_evidence_path <- function(project) file.path(project$path,"intermediate","input_evidence.rds")
input_evidence_stage_run <- function(project,events,peptide) {
  sources <- somatic_source_files(project)
  fp <- checkpoint_fingerprints(inputs=list(events=phase_stage_output(project),
    peptides=peptide_stage_output(project),sources=lapply(sources,function(s)
      s[intersect(names(s),c("kind","sample_id","source_id","sha256","format","options","index_sha256"))])),
    config=list(version=input_evidence_version()),dependencies=list())
  step <- "source_input_evidence"
  if (checkpoint_reusable(project,step,fp) && file.exists(input_evidence_path(project)))
    return(readRDS(input_evidence_path(project)))
  checkpoint_update(project,step,fingerprints=fp)
  checkpoint_update(project,step,"RUNNING",fp)
  tryCatch({
    parsed <- input_source_rows(sources,events,project)
    raw <- somatic_raw_construct(events,peptide,parsed,project$run_id)
    # Retain bound file identities without rehashing large BAMs at the next stage.
    raw$bound_source_files <- sources
    raw$vcf_sources <- lapply(Filter(function(x)x$input_role=="somatic_vcf",project$inputs),
      function(x) list(input_id=x$input_id,sample_id=x$sample_id,path=x$resolved_path,
        sha256=digest::digest(file=x$resolved_path,algo="sha256",serialize=FALSE)))
    somatic_validate_raw(raw,events,peptide)
    temp <- tempfile("input-evidence-",tmpdir=dirname(input_evidence_path(project)))
    on.exit(if(file.exists(temp)) unlink(temp),add=TRUE)
    saveRDS(raw,temp,version=3L)
    if (!file.rename(temp,input_evidence_path(project))) stop("Could not store input evidence")
    checkpoint_update(project,step,"COMPLETED",fp,output_records=list(
      somatic_stage_output(input_evidence_path(project),"intermediate/input_evidence.rds")))
    raw
  },error=function(e) {
    checkpoint_update(project,step,"FAILED",fp,error=list(code="INPUT_EVIDENCE_FAILED"));stop(e)
  })
}
input_evidence_read <- function(project) {
  cp <- Filter(function(x) x$step_id=="source_input_evidence",read_manifest(project$path)$checkpoints)
  if (!length(cp)) return(NULL) # Old completed runs remain readable unchanged.
  if (length(cp)!=1L || cp[[1L]]$status!="COMPLETED") stop("Input evidence incomplete")
  verify_checkpoint_outputs(cp[[1L]],project$path)
  raw <- readRDS(input_evidence_path(project))
  if (!identical(raw$run_id,project$run_id)) stop("Input evidence run identity mismatch")
  raw
}
