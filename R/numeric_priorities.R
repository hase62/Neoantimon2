# Additional, versioned numerical views. Never writes native or external results.
priority_paper_value <- function(mutant,wild,mismatch,vaf,expression,ccf,metric) {
  scalar <- function(x) length(x)==1L && !is.na(x) && is.finite(x)
  inputs <- list(mutant,wild,mismatch,expression)
  if (!all(vapply(inputs,scalar,TRUE))) return(list(without_ccf=NA_real_,
    with_ccf=NA_real_,ccf_factor=NA_real_,status="missing_required_prediction_expression_or_mismatch"))
  if (any(c(mutant,wild,mismatch,expression)<0) || mismatch!=floor(mismatch)) stop("Invalid paper scoring input")
  a <- if (!scalar(vaf)) 0 else vaf
  if (a<0 || a>1) stop("Invalid scoring VAF")
  L <- switch(metric,IC50=function(x)stats::plogis(-.015*(x-500)),
    BA_rank=function(x)stats::plogis(-5*(x-2)),stop("Paper equations require IC50 nM or binding percent rank"))
  base <- L(mutant)*tanh(tanh(5*a)*expression)*(1-2^(-mismatch)*L(wild))
  if (!scalar(ccf)) return(list(without_ccf=base,with_ccf=NA_real_,
    ccf_factor=NA_real_,status="CCF_omitted_not_full_paper_score"))
  if (ccf<0 || ccf>1) stop("Invalid CCF")
  den <- -expm1(-30*(ccf-.8))
  factor <- if (den==0) NA_real_ else 1/den
  score <- base*factor
  status <- if (den==0 || !is.finite(score)) "singular_or_nonfinite" else
    if (abs(den)<1e-6) "finite_unstable_do_not_prioritize" else
      if (factor<0) "finite_negative_do_not_prioritize" else "finite"
  list(without_ccf=base,with_ccf=if(is.finite(score)) score else NA_real_,
    ccf_factor=factor,status=status)
}

priority_metric_category <- function(task) {
  if (is.null(task)) return(NA_character_)
  if (grepl("_binding$",task)) "BA" else if (grepl("_presentation$",task)) "EL" else
    if (identical(task,"external:immunogenicity")) "IM" else "unclassified"
}

priority_paper_metric <- function(category,unit,kind,semantics) {
  if (!identical(category,"BA")) return(NULL)
  if (identical(unit,"nM")) return("IC50")
  if (identical(unit,"percent") && identical(kind,"percentile_rank") &&
      identical(semantics,"mhc_binding_percentile")) return("BA_rank")
  NULL
}

#' Additional model-specific priorities and source-file evidence
#' @param project A completed neo_project.
#' @param analysis Optional completed external analysis handle or ID.
#' @param method model-v1 or neoantimon-paper-2020-v1.
#' @param vaf DNA or RNA VAF channel for the designated adjusted score.
#' @param wild_type_predictions Optional standard prediction TSV for patient-normal
#'   peptides, matched by exact peptide, HLA, provider, model, metric and unit.
#' @param expression_unit Explicit TPM acknowledgement for historical scoring.
#' @return A data frame with one candidate-origin-prediction row, biological
#'   observations, model-specific numerical priorities and full provenance lists.
#' @export
neo_priorities <- function(project, analysis = NULL,
    method = c("model-v1","neoantimon-paper-2020-v1"), vaf = c("dna","rna"),
    wild_type_predictions = NULL, expression_unit = NULL) {
  method <- match.arg(method);vaf <- match.arg(vaf)
  snapshot <- read_runtime_snapshot(project)
  raw <- input_evidence_read(project)
  if (is.null(raw)) raw <- somatic_read_sidecar(project)$raw
  graph <- external_graph(project)
  if (!is.null(analysis)) {
    external <- external_analysis_read(analysis,project)
    if (!identical(external$manifest$fingerprint$source_manifest_sha256,
        external_sha(manifest_file(project$path)))) stop("External analysis belongs to another source project")
  }
  predictions <- if (is.null(analysis)) unlist(lapply(somatic_prediction_runs(project),
    `[[`,"evidence"),recursive=FALSE) else external$payload$evidence
  native <- is.null(analysis)
  wt <- if (is.null(wild_type_predictions)) NULL else external_read_predictions(wild_type_predictions)
  if (method=="neoantimon-paper-2020-v1" && !identical(expression_unit,"TPM"))
    stop("Paper scoring requires explicit expression_unit='TPM'; counts/FPKM cannot be silently substituted")
  event_index <- stats::setNames(raw$event_evidence,vapply(raw$event_evidence,`[[`,"","event_id"))
  rows <- list()
  for (candidate in graph$candidates) {
    own <- Filter(function(p)identical(p$candidate_id,candidate$candidate_id),predictions)
    if (!length(own)) own <- list(NULL)
    for (origin_id in unlist(candidate$origin_ids)) {
      origins <- Filter(function(o)o$origin_id==origin_id,graph$origins)
      if (length(origins)!=1L) stop("Priority origin identity is ambiguous")
      origin <- origins[[1L]]
      ids <- as.character(unlist(origin$source_event_ids))
      ev <- event_index[intersect(ids,names(event_index))]
      single <- if (length(ids)==1L && length(ev)==1L) ev[[1L]] else NULL
      dna <- single$tumor$vaf;rna <- single$rna$vaf
      a <- if(vaf=="dna") dna else rna
      a_missing <- is.null(a)
      expr <- Filter(function(x)x$sample_id==candidate$sample_id &&
        !is.null(x$transcript_id) && sub("\\.[0-9]+$","",x$transcript_id)==origin$transcript_id,raw$expression_evidence)
      best <- function(field) {x <- unlist(lapply(expr,`[[`,field));if(length(x)) max(x) else NULL}
      gene <- best("gene_tpm");tx <- best("transcript_tpm")
      expression <- gene # Paper E is corresponding gene expression, never silently transcript abundance.
      ccf <- if(!is.null(single) && single$ccf$status %in% c("external","conditional")) single$ccf$bounded_value else NULL
      for (prediction in own) {
        score <- if (native) prediction$raw_score else prediction$numeric_score
        evaluated <- !is.null(prediction) && identical(prediction$status,"evaluated") &&
          length(score)==1L && is.finite(score)
        direction <- prediction$score_direction
        if (evaluated && !direction %in% c("higher_is_better","lower_is_better"))
          stop("Numerical priority requires an explicit score direction")
        scalar <- function(x,default=NA_real_) if(is.null(x)) default else x
        text <- function(x) scalar(x,NA_character_)
        row <- list(candidate_id=candidate$candidate_id,sample_id=candidate$sample_id,
          patient_id=origin$patient_id,origin_id=origin_id,origin_ids=candidate$origin_ids,
          event_ids=as.list(ids),peptide_sequence=candidate$peptide_sequence,mhc_key=candidate$mhc_key,
          hla_class=candidate$hla_class,source=if(native) "Native" else "External",
          native_grade=text(candidate$final_grade),
          native_rank_within_sample_class=scalar(candidate$rank_within_sample_class),
          native_technical_eligibility=text(candidate$technical_eligibility),
          native_priority_origin_id=text(candidate$priority_origin_id),
          provider=text(if(native) prediction$predictor else prediction$provider),
          predictor_version=scalar(prediction$predictor_version,""),model_version=scalar(prediction$model_version,""),
          run_id=text(prediction$run_id),task=text(prediction$task),category=priority_metric_category(prediction$task),
          metric=text(prediction$metric),metric_kind=if(native)
            if(identical(prediction$unit,"percent")) "percentile_rank" else "raw_score" else text(prediction$metric_kind),
          metric_semantics=if(native) {
            if(identical(prediction$unit,"percent") && grepl("_binding$",prediction$task)) "mhc_binding_percentile" else
              if(identical(prediction$unit,"percent") && grepl("_presentation$",prediction$task)) "mhc_presentation_percentile" else "raw"
          } else text(prediction$metric_semantics),raw_score=if(evaluated) score else NA_real_,
          unit=text(prediction$unit),score_direction=text(direction),
          prediction_status=if(is.null(prediction)) "NOT_RUN" else prediction$status,
          prediction_only=if(!evaluated) NA_real_ else if(direction=="lower_is_better") -score else score,
          prediction_transform="direction_only; original units retained; within-model only",
          prediction_rank=NA_real_,within_channel_quantile=NA_real_,evidence_adjusted=NA_real_,
          algorithm="model-v1",algorithm_version="1.0.0",requested_method=method,
          prediction_algorithm="direction-only-v1",adjusted_algorithm="within-channel-order-vaf-v1",
          paper_algorithm=if(method=="neoantimon-paper-2020-v1") method else NA_character_,
          paper_equation=NA_character_,ccf_input_type=if(is.null(single)) "unavailable" else single$ccf$status,
          dna_vaf=scalar(dna),rna_vaf=scalar(rna),
          dna_ref_reads=scalar(single$tumor$ref_reads),dna_alt_reads=scalar(single$tumor$alt_reads),
          dna_vaf_denominator=scalar(single$tumor$vaf_denominator),dna_count_unit=text(single$tumor$count_unit),
          rna_ref_reads=scalar(single$rna$ref_reads),rna_alt_reads=scalar(single$rna$alt_reads),
          rna_usable_depth=scalar(single$rna$usable_depth),rna_count_unit=text(single$rna$count_unit),
          purity=scalar(single$purity$value),ploidy=scalar(single$purity$ploidy),
          total_cn=scalar(single$cn$total_cn),major_cn=scalar(single$cn$major_cn),minor_cn=scalar(single$cn$minor_cn),
          vaf_used=if(a_missing) 0 else a,
          vaf_channel=vaf,vaf_missing_substituted=a_missing,
          vaf_note=if(length(ids)>1L) "compound_origin_no_single_VAF" else if(a_missing) "unavailable_observation_remains_NA" else "observed",
          gene_tpm=scalar(gene),transcript_tpm=scalar(tx),ccf=scalar(ccf),
          evidence_adjustment="within-channel empirical order quantile * VAF; no probability interpretation",
          biological_evidence=ev,expression_evidence=expr,prediction_provenance=prediction,
          paper_without_ccf=NA_real_,paper_with_ccf=NA_real_,paper_ccf_factor=NA_real_,
          paper_status="not_requested",paper_wt_provenance=NULL)
        if (method=="neoantimon-paper-2020-v1") {
          row$paper_status <- "missing_or_incompatible_WT_BA_prediction"
          # Modern EL/IM outputs cannot be fed into an IC50/rank equation.
          paper_metric <- priority_paper_metric(row$category,row$unit,row$metric_kind,row$metric_semantics)
          if (!is.null(paper_metric)) row$paper_equation <- if(paper_metric=="IC50") "P_I" else "P_R"
          if (!is.null(paper_metric) && evaluated && !is.null(wt) &&
              !is.null(origin$wt_peptide_sequence) &&
              nchar(origin$wt_peptide_sequence)==nchar(candidate$peptide_sequence)) {
            target <- wt[wt$status=="evaluated" & wt$scope=="peptide_hla" &
              wt$peptide_sequence==origin$wt_peptide_sequence & wt$mhc_key==candidate$mhc_key &
              wt$provider==row$provider & wt$metric==row$metric & wt$task==row$task &
              wt$unit==row$unit & wt$score_direction==row$score_direction &
              wt$metric_kind==row$metric_kind & wt$metric_semantics==row$metric_semantics &
              wt$predictor_version==row$predictor_version & wt$model_version==row$model_version,,drop=FALSE]
            if (nrow(target)>1L) stop("Ambiguous wild-type channel: retain one matching model/run")
            if (nrow(target)==1L) {
              mm <- sum(strsplit(candidate$peptide_sequence,"",fixed=TRUE)[[1L]]!=
                strsplit(origin$wt_peptide_sequence,"",fixed=TRUE)[[1L]])
              p <- priority_paper_value(score,as.numeric(target$score),mm,
                scalar(a),scalar(expression),scalar(ccf),paper_metric)
              row$paper_without_ccf <- p$without_ccf;row$paper_with_ccf <- p$with_ccf
              row$paper_ccf_factor <- p$ccf_factor;row$paper_status <- p$status
              row$paper_wt_provenance <- list(row=as.list(target[1,,drop=FALSE]),
                artifact_sha256=external_sha(wild_type_predictions),mismatch=mm,
                expression_unit="TPM",expression_level=if(is.null(gene)) "transcript" else "gene")
            }
          }
        }
        rows[[length(rows)+1L]] <- row
      }
    }
  }
  if (!length(rows)) return(data.frame())
  # Never mix providers, versions, runs, tasks, metrics, units or samples.
  group_fields <- c("sample_id","hla_class","source","provider","predictor_version",
    "model_version","run_id","task","metric","unit","score_direction")
  keys <- vapply(rows,function(x)paste(unlist(x[group_fields]),collapse="\037"),"")
  for (key in unique(keys)) {
    ix <- which(keys==key & vapply(rows,function(x)is.finite(x$prediction_only),TRUE))
    if (!length(ix)) next
    candidates <- vapply(rows[ix],`[[`,"","candidate_id")
    for (id in unique(candidates)) if (length(unique(vapply(rows[ix[candidates==id]],`[[`,0,"prediction_only")))!=1L)
      stop("Conflicting numerical predictions for candidate/model channel")
    distinct <- ix[!duplicated(candidates)]
    scores <- vapply(rows[distinct],`[[`,0,"prediction_only")
    ranks <- rank(-scores,ties.method="min")
    # Ties remain ties. Quantile is conditional on this exact submitted set.
    q <- if(length(scores)==1L) 1 else (length(scores)-ranks)/(length(scores)-1)
    for(i in ix) {
      at <- match(rows[[i]]$candidate_id,vapply(rows[distinct],`[[`,"","candidate_id"))
      rows[[i]]$prediction_rank <- ranks[[at]];rows[[i]]$within_channel_quantile <- q[[at]]
      rows[[i]]$evidence_adjusted <- q[[at]]*rows[[i]]$vaf_used
    }
  }
  # List columns preserve complete event/origin/source provenance.
  fields <- names(rows[[1L]])
  columns <- lapply(fields,function(f) {
    x <- lapply(rows,`[[`,f)
    if (f %in% c("origin_ids","event_ids","biological_evidence","expression_evidence",
      "prediction_provenance","paper_wt_provenance")) return(I(x))
    if (f=="vaf_missing_substituted") return(vapply(x,isTRUE,TRUE))
    numeric <- f %in% c("raw_score","prediction_only","prediction_rank","within_channel_quantile",
      "evidence_adjusted","dna_vaf","rna_vaf","vaf_used","gene_tpm","transcript_tpm",
      "ccf","paper_without_ccf","paper_with_ccf","paper_ccf_factor",
      "native_rank_within_sample_class","dna_ref_reads","dna_alt_reads","dna_vaf_denominator",
      "rna_ref_reads","rna_alt_reads","rna_usable_depth","purity","ploidy","total_cn","major_cn","minor_cn")
    vapply(x,function(v)if(is.null(v)) if(numeric) NA_real_ else NA_character_ else v,
      if(numeric) 0 else "")
  });names(columns) <- fields
  as.data.frame(columns,check.names=FALSE,stringsAsFactors=FALSE)
}
