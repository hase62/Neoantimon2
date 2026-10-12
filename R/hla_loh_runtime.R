hla_loh_source_files <- function(project) {
  cfg <- project$config$hla_loh
  if (is.null(cfg)) return(list())
  installed <- hla_loh_capability(cfg$execution)
  sample_ids <- vapply(project$samples,`[[`,"","sample_id")
  lapply(cfg$sources,function(s) {
    if (!s$sample_id %in% sample_ids) stop("LOHHLA source sample not in project")
    path <- if (is.null(s$path)) NULL else resolve_path(s$path,project$config$paths$input_root)
    if (s$provider_status=="completed" && (is.null(path) ||
        !file.exists(path) || dir.exists(path))) stop("LOHHLA output file missing")
    sha <- if (!is.null(path) && file.exists(path) && !dir.exists(path))
      digest::digest(file=path,algo="sha256",serialize=FALSE) else NULL
    if (!is.null(s$checksum) && !identical(tolower(s$checksum),sha))
      stop("LOHHLA output checksum mismatch")
    if (is.null(s$provider_version)) s$provider_version <- installed$provider_version
    s$resolved_path <- path; s$sha256 <- sha
    s$run_id <- if (is.null(s$run_id)) paste0("lohhla-run:",infra_hash(list(
      sample=s$sample_id,context=s$context_type,sha256=sha,inputs=s$inputs,
      provider_status=s$provider_status,failure_reason=s$failure_reason))) else s$run_id
    s
  })
}

hla_loh_fingerprints <- function(project,sources,policy) checkpoint_fingerprints(
  inputs=list(samples=project$samples,contexts=project$analysis_context,
    genotype=project$hla,raw_outputs=lapply(sources,function(s) list(
      sample_id=s$sample_id,context_type=s$context_type,sha256=s$sha256,
      run_id=s$run_id,provider_version=s$provider_version,inputs=s$inputs,
      provider_status=s$provider_status,failure_reason=s$failure_reason))),
  config=list(schema=hla_loh_schema_version(),parser=hla_loh_parser_version(),policy=policy),
  dependencies=list(hla_reference=lapply(sources,function(s)s$inputs$hla_fasta),
    purity_ploidy=lapply(sources,function(s)s$inputs$purity_ploidy)))

hla_loh_context <- function(project,sample_id) {
  found <- Filter(function(x) sample_id %in% unlist(x$sample_bindings),
    project$analysis_context)
  if (length(found)!=1L) stop("Ambiguous HLA LOH sample context")
  found[[1L]]
}

hla_loh_genotype <- function(project,sample_id) {
  sample <- Filter(function(x)x$sample_id==sample_id,project$samples)
  if (length(sample)!=1L) stop("Unknown HLA LOH sample")
  alleles <- project$hla[[sample[[1L]]$patient_id]]
  if (is.null(alleles)) return(character())
  found <- unique(vapply(alleles,function(x) {
    p <- peptide_presenters(list(x))[[1L]]
    if (p$status=="RESOLVED" && p$hla_class=="I" &&
        grepl("^HLA-[ABC]\\*",p$hla_allele)) p$hla_allele else ""
  },""))
  found[nzchar(found)]
}

hla_loh_raw_number <- function(value,lower=0,upper=Inf) {
  n <- hla_loh_number(value,lower,upper)
  if (inherits(n,"hla_loh_invalid")) return(n)
  n
}

hla_loh_record <- function(project,source,row,index,match,other_match,role,policy) {
  v <- row$values
  own <- if (index==1L) "1" else "2"
  other <- if (index==1L) "2" else "1"
  get <- function(stem,which) hla_loh_raw_number(v[[paste0(stem,which)]])
  context <- hla_loh_context(project,source$sample_id)
  token <- v[[paste0("allele",index)]]
  canonical <- match$allele
  locus <- hla_loh_token(token)
  locus <- if (is.null(locus)) NULL else sub("\\*.*$","",locus)
  record <- list(schema_version=hla_loh_schema_version(),sample_id=source$sample_id,
    context_type=source$context_type,analysis_context=context$context_id,
    hla_locus=locus,canonical_hla_allele=canonical,
    other_hla_allele=other_match$allele,allele_role=role,provider="lohhla",
    provider_version=source$provider_version,run_id=source$run_id,
    raw_output_artifact=source$resolved_path,raw_output_checksum=source$sha256,
    raw_row_index=row$row_index,raw_values=row$raw,
    raw_provider_allele=token,raw_other_allele=v[[paste0("allele",3L-index)]],
    allele_copy_number=get("cn",own),
    allele_copy_number_lower=get("cn",paste0(own,"lower")),
    allele_copy_number_upper=get("cn",paste0(own,"upper")),
    other_allele_copy_number=get("cn",other),
    other_allele_copy_number_lower=get("cn",paste0(other,"lower")),
    other_allele_copy_number_upper=get("cn",paste0(other,"upper")),
    pval_unique=hla_loh_raw_number(v$pval,0,1),
    unpaired_pval_unique=hla_loh_raw_number(v$unpaired,0,1),
    num_mismatch_sites_covered=hla_loh_raw_number(v$mismatch),
    proportion_supportive_sites=hla_loh_raw_number(v$supportive,0,1),
    provider_loss_allele=v$loss,provider_kept_allele=v$kept,
    mapping_status=match$status,input_provenance=source$inputs,
    parser_version=hla_loh_parser_version())
  record$evidence_id <- paste0("hla-loh-evidence:",infra_hash(record))
  hla_loh_interpret(record,policy)
}

hla_loh_absent_record <- function(project,sample_id,allele,role,policy,
                                  mapping_status="mapped",source=NULL) {
  context <- hla_loh_context(project,sample_id)
  record <- list(schema_version=hla_loh_schema_version(),sample_id=sample_id,
    context_type=context$mode,analysis_context=context$context_id,
    hla_locus=sub("\\*.*$","",allele),canonical_hla_allele=allele,
    other_hla_allele=NULL,allele_role=role,provider="lohhla",
    provider_version=if (is.null(source)) NULL else source$provider_version,
    run_id=if (is.null(source)) NULL else source$run_id,
    raw_output_artifact=if (is.null(source)) NULL else source$resolved_path,
    raw_output_checksum=if (is.null(source)) NULL else source$sha256,
    allele_copy_number=NULL,allele_copy_number_lower=NULL,
    allele_copy_number_upper=NULL,other_allele_copy_number=NULL,
    other_allele_copy_number_lower=NULL,other_allele_copy_number_upper=NULL,
    pval_unique=NULL,unpaired_pval_unique=NULL,num_mismatch_sites_covered=NULL,
    proportion_supportive_sites=NULL,provider_loss_allele=NULL,
    provider_kept_allele=NULL,mapping_status=mapping_status,
    provider_status=if (is.null(source)) "not_available" else source$provider_status,
    failure_reason=if (is.null(source)) NULL else source$failure_reason,
    input_provenance=if (is.null(source)) NULL else source$inputs,
    parser_version=hla_loh_parser_version())
  record$evidence_id <- paste0("hla-loh-evidence:",infra_hash(record))
  hla_loh_interpret(record,policy)
}

hla_loh_construct <- function(project,sources,policy) {
  out <- list()
  for (sample in project$samples) {
    genotype <- hla_loh_genotype(project,sample$sample_id)
    source <- Filter(function(s)s$sample_id==sample$sample_id,sources)
    source <- if (length(source)) source[[1L]] else NULL
    parsed <- if (is.null(source) || source$provider_status=="failed") NULL else
      hla_loh_parse(source$resolved_path)
    if (!is.null(parsed) && length(parsed$rows)) {
      regions <- unique(vapply(parsed$rows,function(r) {
        value <- r$values$region
        if (is.null(value)) "" else value
      },""))
      if (length(regions)>1L) stop("unsupported_format: multiple regions in one LOHHLA source")
    }
    if (!is.null(source) && !identical(source$context_type,
        hla_loh_context(project,sample$sample_id)$mode))
      stop("LOHHLA source context_type does not match sample context")
    observed <- character()
    if (!is.null(parsed)) for (row in parsed$rows) {
      v <- row$values
      tokens <- list(v$allele1,v$allele2)
      first <- hla_loh_token(tokens[[1L]]); second <- hla_loh_token(tokens[[2L]])
      locus1 <- if (is.null(first)) NULL else sub("\\*.*$","",first)
      locus2 <- if (is.null(second)) NULL else sub("\\*.*$","",second)
      if (is.null(locus1) || !identical(locus1,locus2))
        stop("unsupported_format: LOHHLA row has invalid or mixed HLA loci")
      if (!locus1 %in% hla_loh_registry()$lohhla$supported_loci)
        stop("unsupported_format: LOHHLA locus outside Class I scope")
      mapping <- lapply(tokens,hla_loh_match,genotype=genotype,locus=locus1)
      if (!is.null(mapping[[1L]]$allele) && identical(mapping[[1L]]$allele,mapping[[2L]]$allele))
        stop("unsupported_format: duplicate provider alleles in locus")
      loss <- hla_loh_token(v$loss); kept <- hla_loh_token(v$kept)
      roles_valid <- !is.null(loss) && !is.null(kept) && !identical(loss,kept) &&
        setequal(c(loss,kept),vapply(tokens,hla_loh_token,""))
      germline_alleles <- genotype[startsWith(genotype,paste0(locus1,"*"))]
      for (i in 1:2) {
        role <- if (length(germline_alleles)==1L) "homozygous" else
          if (!roles_valid) "unresolved" else
          if (identical(hla_loh_token(tokens[[i]]),loss))
          "potential_loss" else if (!is.null(kept) &&
            identical(hla_loh_token(tokens[[i]]),kept)) "kept" else "unresolved"
        record <- hla_loh_record(project,source,row,i,mapping[[i]],mapping[[3L-i]],role,policy)
        out[[length(out)+1L]] <- record
        if (!is.null(record$canonical_hla_allele)) observed <- c(observed,record$canonical_hla_allele)
      }
    }
    for (locus in hla_loh_registry()$lohhla$supported_loci) {
      at_locus <- genotype[startsWith(genotype,paste0(locus,"*"))]
      unresolved <- Filter(function(x) x$sample_id==sample$sample_id &&
        identical(x$hla_locus,locus) && x$mapping_status!="mapped",out)
      unresolved_status <- if (any(vapply(unresolved,function(x)
        x$mapping_status=="ambiguous",TRUE))) "ambiguous" else
        if (length(unresolved)) "unresolved" else "mapped"
      for (allele in setdiff(at_locus,observed)) out[[length(out)+1L]] <-
        hla_loh_absent_record(project,sample$sample_id,allele,
          if (length(at_locus)==1L) "homozygous" else
            if (!is.null(source) && source$provider_status=="failed")
              "provider_failed" else "unreported",policy,
          mapping_status=if (length(at_locus)==1L) "mapped" else unresolved_status,
          source=source)
    }
  }
  keys <- vapply(out,function(r) paste(r$sample_id,r$analysis_context,
    if (is.null(r$canonical_hla_allele)) paste0("raw:",r$raw_row_index,":",r$raw_provider_allele)
    else r$canonical_hla_allele,sep="\037"),"")
  if (anyDuplicated(keys)) stop("Duplicate HLA LOH evidence for sample/context/allele")
  if (length(keys)) out <- out[order(keys,method="radix")]
  out
}

hla_loh_link <- function(rows,evidence) {
  mapped <- Filter(function(e)!is.null(e$canonical_hla_allele),evidence)
  by_key <- stats::setNames(mapped,vapply(mapped,function(e)
    paste(e$sample_id,e$analysis_context,e$canonical_hla_allele,sep="\037"),""))
  lapply(rows,function(row) {
    if (row$hla_class!="I" || is.null(row$hla_allele) ||
        !grepl("^HLA-[ABC]\\*",row$hla_allele))
      return(list(candidate_id=row$candidate_id,sample_id=row$sample_id,
        hla_class=row$hla_class,
        presenter_hla_loh_status=if (row$hla_class=="II") "not_applicable" else
          "unsupported_locus",presenter_hla_loh_call=NULL,
        presenter_hla_copy_number=NULL,presenter_hla_loh_pval_unique=NULL,
        presenter_hla_loh_evidence_id=NULL,presenter_escape_status="not_applicable",
        paired_hla_loss_status=NULL,generic_hla_loh_consistency="not_comparable"))
    key <- paste(row$sample_id,row$analysis_context,row$hla_allele,sep="\037")
    e <- if (key %in% names(by_key)) by_key[[key]] else NULL
    if (is.null(e)) return(list(candidate_id=row$candidate_id,sample_id=row$sample_id,
      hla_class=row$hla_class,
      presenter_hla_loh_status="not_evaluable",presenter_hla_loh_call="not_evaluable",
      presenter_hla_copy_number=NULL,presenter_hla_loh_pval_unique=NULL,
      presenter_hla_loh_evidence_id=NULL,presenter_escape_status="unresolved",
      paired_hla_loss_status=NULL,generic_hla_loh_consistency="insufficient_evidence"))
    paired <- Filter(function(x)x$sample_id==e$sample_id &&
      x$analysis_context==e$analysis_context && x$hla_locus==e$hla_locus &&
      !is.null(x$canonical_hla_allele) && x$canonical_hla_allele!=e$canonical_hla_allele,
      evidence)
    paired_loss <- any(vapply(paired,function(x)x$loh_call=="loss_supported",TRUE))
    # Frozen generic_loh is tied to contributing variant events and may be
    # elsewhere in the genome. No HLA-region equivalence is established here.
    consistency <- "not_comparable"
    list(candidate_id=row$candidate_id,sample_id=row$sample_id,
      hla_class=row$hla_class,
      presenter_hla_loh_status=e$evaluation_status,presenter_hla_loh_call=e$loh_call,
      presenter_hla_copy_number=if (is.numeric(e$allele_copy_number))
        e$allele_copy_number else NULL,
      presenter_hla_loh_pval_unique=if (is.numeric(e$pval_unique))
        e$pval_unique else NULL,
      presenter_hla_loh_evidence_id=e$evidence_id,
      presenter_escape_status=if (e$loh_call=="loss_supported")
        "hla_loss_supported" else if (e$evaluation_status=="evaluated")
        "hla_loss_not_supported" else "unresolved",
      paired_hla_loss_status=if (paired_loss) "paired_homolog_loss_supported" else NULL,
      generic_hla_loh_consistency=consistency)
  })
}

hla_loh_summary <- function(project,evidence,links,capability) {
  statuses <- vapply(evidence,`[[`,"","loh_call")
  link_status <- vapply(links,`[[`,"","presenter_escape_status")
  class_ii <- sum(vapply(links,function(x)x$hla_class=="II",TRUE))
  evaluated <- Filter(function(x)x$evaluation_status=="evaluated",evidence)
  list(provider="lohhla",execution_status=capability$execution_status,
    samples_with_evaluable_class_i=as.list(sort(unique(vapply(evaluated,`[[`,"","sample_id")),
      method="radix")),alleles_evaluated=length(evaluated),
    locus_evaluated_counts=as.list(table(vapply(evaluated,`[[`,"","hla_locus"))),
    allele_status_counts=as.list(table(statuses)),candidate_escape_counts=as.list(table(link_status)),
    class_ii_candidates_not_applicable=class_ii,
    unresolved_or_not_evaluable_alleles=sum(!vapply(evidence,function(x)
      x$evaluation_status=="evaluated",TRUE)),
    provider_failure_count=sum(vapply(evidence,function(x)
      x$evaluation_status=="provider_failed",TRUE)),
    generic_hla_consistency_counts=as.list(table(vapply(links,`[[`,"",
      "generic_hla_loh_consistency"))))
}

hla_loh_sidecar_dir <- function(project,fp) file.path(project$path,"results",
  "hla_loh",infra_hash(fp))

hla_loh_write <- function(project,fp,evidence,links,capability) {
  root <- hla_loh_sidecar_dir(project,fp)
  if (!dir.exists(root) && !dir.create(root,recursive=TRUE,showWarnings=FALSE))
    stop("Cannot create HLA LOH result directory")
  summary <- hla_loh_summary(project,evidence,links,capability)
  payload <- list(schema_version=hla_loh_schema_version(),run_id=project$run_id,
    fingerprints=fp,evidence=evidence,links=links,summary=summary,capability=capability)
  saveRDS(payload,file.path(root,"hla_loh.rds"),version=3L)
  atomic_json(evidence,file.path(root,"allele_evidence.json"))
  atomic_json(links,file.path(root,"candidate_evidence.json"))
  atomic_json(summary,file.path(root,"summary.json"))
  somatic_write_tsv(links,c("candidate_id","sample_id","presenter_hla_loh_status",
    "presenter_hla_loh_call","presenter_hla_copy_number",
    "presenter_hla_loh_pval_unique","presenter_hla_loh_evidence_id",
    "presenter_escape_status","paired_hla_loss_status","generic_hla_loh_consistency"),
    file.path(root,"candidate_evidence.tsv"))
  writeLines(c("# HLA-specific LOH evidence", "",
    "LOHHLA fixture or user-supplied output is interpreted separately from generic LOH.",
    "The frozen somatic-priority-v1 grade, rank and candidate identity are unchanged.","",
    paste0("Evaluated Class I alleles: ",summary$alleles_evaluated),
    paste0("Class II candidates not applicable: ",summary$class_ii_candidates_not_applicable),
    paste0("Execution status: ",summary$execution_status)),
    file.path(root,"report.md"))
  files <- c("hla_loh.rds","allele_evidence.json","candidate_evidence.json",
    "candidate_evidence.tsv","summary.json","report.md")
  hashes <- lapply(files,function(f) list(file=f,sha256=digest::digest(
    file=file.path(root,f),algo="sha256",serialize=FALSE)))
  atomic_json(list(schema_version=hla_loh_schema_version(),fingerprints=fp,
    outputs=hashes),file.path(root,"checkpoint.json"))
  payload
}

hla_loh_read <- function(project,priority_rows=NULL) {
  if (is.null(project$config$hla_loh)) return(NULL)
  sources <- hla_loh_source_files(project)
  fp <- hla_loh_fingerprints(project,sources,project$config$hla_loh$policy)
  root <- hla_loh_sidecar_dir(project,fp)
  checkpoint <- file.path(root,"checkpoint.json")
  if (!file.exists(checkpoint)) stop("HLA LOH sidecar missing; run neo_run")
  cp <- jsonlite::read_json(checkpoint,simplifyVector=FALSE)
  if (!identical(cp$fingerprints,fp) || !identical(cp$schema_version,hla_loh_schema_version()))
    stop("HLA LOH fingerprint mismatch")
  for (record in cp$outputs) if (!identical(record$sha256,digest::digest(
    file=file.path(root,record$file),algo="sha256",serialize=FALSE)))
      stop("HLA LOH sidecar checksum mismatch")
  payload <- readRDS(file.path(root,"hla_loh.rds"))
  if (!identical(payload$run_id,project$run_id) || !identical(payload$fingerprints,fp))
    stop("HLA LOH sidecar identity mismatch")
  if (!is.null(priority_rows) && !all(vapply(priority_rows,`[[`,"","candidate_id") %in%
      vapply(payload$links,`[[`,"","candidate_id")))
      stop("HLA LOH candidate identity mismatch")
  payload
}

hla_loh_postrun <- function(project,snapshot) {
  if (is.null(project$config$hla_loh)) return(invisible(NULL))
  sources <- hla_loh_source_files(project)
  fp <- hla_loh_fingerprints(project,sources,project$config$hla_loh$policy)
  root <- hla_loh_sidecar_dir(project,fp)
  if (file.exists(file.path(root,"checkpoint.json")))
    return(invisible(hla_loh_read(project)))
  if (dir.exists(root)) stop("Incomplete HLA LOH sidecar directory")
  priority <- somatic_priority_read(project,snapshot)
  rows <- if (is.null(priority)) snapshot$neoantigen_results else priority$rows
  evidence <- hla_loh_construct(project,sources,project$config$hla_loh$policy)
  links <- hla_loh_link(rows,evidence)
  cap <- hla_loh_capability(project$config$hla_loh$execution)
  hla_loh_write(project,fp,evidence,links,cap)
  verified <- hla_loh_read(project,rows)
  if (!setequal(vapply(verified$links,`[[`,"","candidate_id"),
      vapply(rows,`[[`,"","candidate_id")))
    stop("HLA LOH sidecar candidate set mismatch")
  invisible(verified)
}
