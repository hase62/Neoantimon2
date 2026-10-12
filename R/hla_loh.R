# HLA-specific evidence is an optional sidecar. It does not participate in
# somatic-priority-v1 or the frozen peptide/presenter identity contracts.
hla_loh_schema_version <- function() "hla-allele-specific-loh-v1"
hla_loh_parser_version <- function() "lohhla-named-columns-v1"

hla_loh_config <- function(x) {
  strict_keys(x,c("sources","policy","execution"),"HLA LOH configuration")
  sources <- if (is.null(x$sources)) list() else x$sources
  if (!is.list(sources) || !is.null(names(sources))) stop("HLA LOH sources must be a record array")
  for (i in seq_along(sources)) {
    s <- sources[[i]]
    strict_keys(s,c("sample_id","path","checksum","provider_version","run_id",
      "context_type","inputs","provider_status","failure_reason"),"LOHHLA source")
    scalar_text(s$sample_id,"LOHHLA sample_id")
    if (is.null(s$provider_status)) s$provider_status <- "completed"
    if (!s$provider_status %in% c("completed","failed"))
      stop("Invalid LOHHLA provider status")
    if (s$provider_status=="completed") scalar_text(s$path,"LOHHLA output path")
    if (!is.null(s$path)) scalar_text(s$path,"LOHHLA output path")
    for (key in c("checksum","provider_version","run_id","context_type"))
      if (!is.null(s[[key]])) scalar_text(s[[key]],key)
    if (!is.null(s$checksum) && !grepl("^[a-fA-F0-9]{64}$",s$checksum))
      stop("LOHHLA checksum must be SHA256 hex")
    if (is.null(s$context_type)) s$context_type <- "somatic"
    if (!is.null(s$failure_reason)) scalar_text(s$failure_reason,"LOHHLA failure_reason")
    if (is.null(s$inputs)) s$inputs <- list()
    strict_keys(s$inputs,c("tumor_bam","normal_bam","hla_genotype","hla_fasta",
      "purity_ploidy","copy_number","command","dependencies"),"LOHHLA input provenance")
    sources[[i]] <- s
  }
  keys <- vapply(sources,function(s) paste(s$sample_id,s$context_type,sep="\037"),"")
  if (anyDuplicated(keys)) stop("Duplicate LOHHLA source sample/context")
  if (length(keys)) sources <- sources[order(keys,method="radix")]
  p <- if (is.null(x$policy)) list() else x$policy
  strict_keys(p,c("policy_version","copy_number_loss_threshold","pval_unique_threshold",
    "copy_number_comparison","pval_unique_comparison"),"HLA LOH policy")
  defaults <- list(policy_version="hla-loh-policy-v1",copy_number_loss_threshold=0.5,
    pval_unique_threshold=0.01,copy_number_comparison="<",pval_unique_comparison="<")
  for (key in names(p)) defaults[[key]] <- p[[key]]
  scalar_text(defaults$policy_version,"HLA LOH policy_version")
  for (key in c("copy_number_loss_threshold","pval_unique_threshold")) {
    v <- defaults[[key]]
    if (!is.numeric(v) || length(v)!=1L || !is.finite(v) || v<0 ||
        (key=="pval_unique_threshold" && v>1)) stop("Invalid HLA LOH threshold")
  }
  for (key in c("copy_number_comparison","pval_unique_comparison"))
    if (!identical(defaults[[key]],"<") && !identical(defaults[[key]],"<="))
      stop("HLA LOH comparison must be < or <=")
  execution <- if (is.null(x$execution)) NULL else x$execution
  if (!is.null(execution)) execution <- hla_loh_execution(execution)
  list(sources=sources,policy=defaults,execution=execution)
}

hla_loh_execution <- function(execution) {
  strict_keys(execution,c("script","tumor_bam","normal_bam","hla_calls",
    "hla_fasta","copy_number","novoalign_dir","gatk_dir","output_dir",
    "purity_ploidy","options","interface"),"LOHHLA execution")
  if (is.null(execution$interface)) execution$interface <- "paired_bam_v1"
  if (!identical(execution$interface,"paired_bam_v1") &&
      !identical(execution$interface,"bam_directory_v1"))
    stop("Unsupported LOHHLA execution interface")
  if (is.null(execution$options)) execution$options <- list()
  strict_keys(execution$options,c("mapping_step","fishing_step"),"LOHHLA options")
  if (identical(execution$interface,"bam_directory_v1") && length(execution$options))
    stop("LOHHLA bam_directory_v1 does not support execution options")
  for (key in names(execution$options)) {
    value <- execution$options[[key]]
    if (!is.logical(value) || length(value)!=1L || is.na(value))
      stop("LOHHLA ",key," must be a scalar logical TRUE or FALSE")
  }
  execution
}

hla_loh_tool <- function(name) unname(Sys.which(name))

hla_loh_r_packages <- function(rscript) {
  expected <- file.path(R.home("bin"),"Rscript")
  if (!nzchar(rscript) || !file.exists(expected) ||
      !identical(normalizePath(rscript,mustWork=FALSE),normalizePath(expected)))
    return(list(status="unresolved",missing=character()))
  required <- c("optparse","seqinr","Biostrings","beeswarm","zoo","Rsamtools")
  missing <- required[!vapply(required,requireNamespace,TRUE,quietly=TRUE)]
  list(status=if (length(missing)) "missing" else "satisfied",missing=missing)
}

hla_loh_registry <- function() list(lohhla=list(
  provider_type=c("hla_allele_specific_copy_number","hla_loh"),
  supported_loci=c("HLA-A","HLA-B","HLA-C"),
  input_concepts=c("tumor_bam","normal_bam","hla_genotype","hla_fasta",
    "copy_number","purity_ploidy","novoalign","gatk_picard"),
  parser_capability="named_loss_prediction_columns",
  execution_capability="explicit_local_command_construction",
  license="external_user_supplied"))

hla_loh_capability <- function(execution=NULL) {
  spec <- hla_loh_registry()$lohhla
  if (is.null(execution)) return(c(spec,list(parser_available=TRUE,
    command_constructible=FALSE,execution_dependencies_satisfied=FALSE,
    available=FALSE,execution_status="not_available",missing=list(),
    provider_version=NULL,version_status="unknown")))
  execution <- hla_loh_execution(execution)
  required <- c("script","tumor_bam","normal_bam","hla_calls","hla_fasta",
    "copy_number","novoalign_dir","gatk_dir","output_dir")
  missing <- required[vapply(required,function(k) is.null(execution[[k]]) ||
    !is.character(execution[[k]]) || length(execution[[k]])!=1L ||
    is.na(execution[[k]]) || !nzchar(execution[[k]]),TRUE)]
  command_constructible <- !length(missing)
  dependency_status <- "unresolved"
  dependency_paths <- list()
  r_packages <- list(status="unresolved",missing=character())
  if (command_constructible) {
    files <- setdiff(required,c("novoalign_dir","gatk_dir","output_dir"))
    missing <- c(missing,files[!vapply(execution[files],file.exists,TRUE)])
    dirs <- c("novoalign_dir","gatk_dir")
    missing <- c(missing,dirs[!vapply(execution[dirs],dir.exists,TRUE)])
    binaries <- c("novoalign","novoindex")
    missing <- c(missing,binaries[!vapply(binaries,function(x)
      file.exists(file.path(execution$novoalign_dir,x)),TRUE)])
    jars <- c("SamToFastq.jar","SortSam.jar","FilterSamReads.jar")
    missing <- c(missing,jars[!vapply(jars,function(x)
      file.exists(file.path(execution$gatk_dir,x)),TRUE)])
    dependencies <- c("samtools","bedtools","java")
    if (!identical(execution$options$mapping_step,FALSE) &&
        !identical(execution$options$fishing_step,FALSE))
      dependencies <- c(dependencies,"jellyfish")
    dependency_paths <- as.list(stats::setNames(vapply(dependencies,hla_loh_tool,""),
      dependencies))
    missing <- c(missing,dependencies[!nzchar(unlist(dependency_paths))])
    if (identical(execution$interface,"bam_directory_v1") &&
        !identical(dirname(normalizePath(execution$tumor_bam,mustWork=TRUE)),
          dirname(normalizePath(execution$normal_bam,mustWork=TRUE))))
      missing <- c(missing,"bam_directory_mismatch")
    rscript <- hla_loh_tool("Rscript")
    dependency_paths$Rscript <- rscript
    if (!nzchar(rscript)) missing <- c(missing,"Rscript")
    else {
      r_packages <- hla_loh_r_packages(rscript)
      if (length(r_packages$missing))
        missing <- c(missing,paste0("R:",r_packages$missing))
    }
    dependency_status <- if (length(missing)) "missing" else r_packages$status
  }
  script <- execution$script
  version <- NULL
  if (is.character(script) && length(script)==1L && file.exists(script)) {
    # A user-supplied checkout may carry an exact git commit. No command runs
    # the LOHHLA script merely to determine its version.
    root <- dirname(normalizePath(script,mustWork=TRUE))
    if (dir.exists(file.path(root,".git")) && nzchar(Sys.which("git"))) {
      found <- suppressWarnings(system2("git",c("-C",shQuote(root),"rev-parse","HEAD"),
        stdout=TRUE,stderr=FALSE))
      if (length(found)==1L && grepl("^[a-f0-9]{40,64}$",found)) version <- found
    }
  }
  ready <- command_constructible && identical(dependency_status,"satisfied")
  c(spec,list(parser_available=TRUE,interface=execution$interface,
    options=execution$options,script=execution$script,
    command_constructible=command_constructible,
    execution_dependencies_satisfied=if (identical(dependency_status,"unresolved"))
      NA else identical(dependency_status,"satisfied"),
    dependency_status=dependency_status,dependency_paths=dependency_paths,
    r_package_status=r_packages$status,r_packages_missing=as.list(r_packages$missing),
    available=ready,execution_status=if (ready) "command_ready" else if (
      command_constructible && identical(dependency_status,"unresolved"))
      "dependencies_unresolved" else "not_available",missing=as.list(unique(missing)),
    provider_version=version,version_status=if (is.null(version)) "unknown" else "commit"))
}

hla_loh_command <- function(execution, sample_id) {
  execution <- hla_loh_execution(execution)
  cap <- hla_loh_capability(execution)
  scalar_text(sample_id,"LOHHLA sample_id")
  if (!cap$command_constructible) return(list(execution_status=cap$execution_status,
    capability=cap,executed=FALSE))
  bam_args <- if (execution$interface=="paired_bam_v1")
    c("--normalBAMfile",execution$normal_bam,"--tumorBAMfile",execution$tumor_bam) else
    c("--normalBAMfile",execution$normal_bam,"--BAMDir",dirname(execution$tumor_bam))
  args <- c(execution$script,"--patientId",sample_id,"--outputDir",execution$output_dir,
    bam_args,
    "--hlaPath",execution$hla_calls,"--HLAfastaLoc",execution$hla_fasta,
    "--CopyNumLoc",execution$copy_number,"--novoDir",execution$novoalign_dir,
    "--gatkDir",execution$gatk_dir)
  flags <- c(mapping_step="--mappingStep",fishing_step="--fishingStep")
  for (key in names(execution$options))
    args <- c(args,flags[[key]],if (execution$options[[key]]) "TRUE" else "FALSE")
  list(execution_status=cap$execution_status,interface=execution$interface,
    options=execution$options,executable=cap$dependency_paths$Rscript,
    script=execution$script,args=args,command_identity=infra_hash(list(
      interface=execution$interface,options=execution$options,
      executable=cap$dependency_paths$Rscript,script=execution$script,
      args=as.list(args),dependencies=cap$dependency_paths)),
    capability=cap,executed=FALSE)
}

hla_loh_header_aliases <- function() list(
  allele1=c("HLA_A_type1","HLA_allele1","HLA_allele_1"),
  allele2=c("HLA_A_type2","HLA_allele2","HLA_allele_2"),
  pval=c("Pval_unique","PVal_unique"),
  unpaired=c("UnPairedPval_unique","UnpairedPval_unique"),
  loss=c("LossAllele"),kept=c("KeptAllele"),
  cn1=c("HLA_type1copyNum_withBAFBin"),
  cn1lower=c("HLA_type1copyNum_withBAFBin_lower"),
  cn1upper=c("HLA_type1copyNum_withBAFBin_upper"),
  cn2=c("HLA_type2copyNum_withBAFBin"),
  cn2lower=c("HLA_type2copyNum_withBAFBin_lower"),
  cn2upper=c("HLA_type2copyNum_withBAFBin_upper"),
  mismatch=c("numMisMatchSitesCov"),supportive=c("propSupportiveSites"),
  region=c("region","Region","sample_id","Sample"))

hla_loh_columns <- function(header) {
  if (anyDuplicated(header) || any(!nzchar(header))) stop("unsupported_format: duplicate or empty columns")
  aliases <- hla_loh_header_aliases()
  selected <- lapply(aliases,function(xs) {
    hit <- intersect(xs,header)
    if (length(hit)>1L) stop("unsupported_format: duplicate semantic columns")
    if (length(hit)) hit[[1L]] else NULL
  })
  required <- c("allele1","allele2","pval","loss","kept","cn1","cn2")
  if (any(vapply(selected[required],is.null,TRUE)))
    stop("unsupported_format: required LOHHLA columns missing")
  selected
}

hla_loh_parse <- function(path) {
  if (!file.exists(path) || dir.exists(path)) stop("LOHHLA output file missing")
  lines <- readLines(path,warn=FALSE)
  if (!length(lines)) stop("unsupported_format: empty LOHHLA output")
  if (length(lines)==1L && grepl("^All HLA alleles homozygous\\. No LOH analysis can be performed\\.$",
      lines[[1L]])) return(list(parser_version=hla_loh_parser_version(),
        columns=list(),rows=list(),provider_status="homozygous_not_evaluable"))
  if (length(lines)==1L && grepl("^No suitable HLA alleles! No LOH analysis can be performed\\.$",
      lines[[1L]])) return(list(parser_version=hla_loh_parser_version(),
        columns=list(),rows=list(),provider_status="not_evaluable"))
  header <- strsplit(sub("\r$","",lines[[1L]]),"\t",fixed=TRUE)[[1L]]
  columns <- hla_loh_columns(header)
  parsed <- tryCatch(utils::read.delim(path,header=TRUE,check.names=FALSE,
    colClasses="character",na.strings=character(),quote="",comment.char="",
    fill=FALSE,stringsAsFactors=FALSE),error=function(e)
      stop("unsupported_format: malformed LOHHLA rows"))
  if (!identical(names(parsed),header)) stop("unsupported_format: header mismatch")
  rows <- lapply(seq_len(nrow(parsed)),function(i) {
    row <- as.list(parsed[i,,drop=FALSE])
    names(row) <- names(parsed)
    list(row_index=i,raw=row,values=lapply(columns,function(name)
      if (is.null(name)) NULL else row[[name]]))
  })
  list(parser_version=hla_loh_parser_version(),columns=columns,rows=rows)
}

hla_loh_token <- function(raw) {
  if (is.null(raw) || length(raw)!=1L || is.na(raw)) return(NULL)
  x <- toupper(trimws(raw))
  x <- sub("^HLA[-_]", "", x)
  x <- gsub("_", ":", x, fixed=TRUE)
  x <- sub("^([ABC]):", "\\1*", x)
  if (!grepl("^[ABC]\\*[0-9]{2,3}(:[0-9]{2,3}){1,3}[NLSCAQ]?$",x)) return(NULL)
  paste0("HLA-",x)
}

hla_loh_match <- function(token, genotype, locus=NULL) {
  normalized <- hla_loh_token(token)
  if (is.null(normalized)) return(list(status="unresolved",allele=NULL))
  actual_locus <- sub("\\*.*$","",normalized)
  if (!is.null(locus) && !identical(actual_locus,locus))
    return(list(status="unresolved",allele=NULL))
  canonical <- vapply(genotype,function(x) {
    p <- peptide_presenters(list(x))[[1L]]
    if (p$status=="RESOLVED" && p$hla_class=="I") p$hla_allele else ""
  },"")
  canonical <- unique(canonical[nzchar(canonical) & startsWith(canonical,paste0(actual_locus,"*"))])
  compatible <- canonical[vapply(canonical,function(x) {
    a <- sub("^[^*]+\\*","",x); b <- sub("^[^*]+\\*","",normalized)
    af <- strsplit(sub("[NLSCAQ]$","",a),":",fixed=TRUE)[[1L]]
    bf <- strsplit(sub("[NLSCAQ]$","",b),":",fixed=TRUE)[[1L]]
    suffix_a <- if (grepl("[NLSCAQ]$",a)) sub(".*([NLSCAQ])$","\\1",a) else ""
    suffix_b <- if (grepl("[NLSCAQ]$",b)) sub(".*([NLSCAQ])$","\\1",b) else ""
    n <- min(length(af),length(bf))
    identical(af[seq_len(n)],bf[seq_len(n)]) &&
      (!nzchar(suffix_b) || identical(suffix_a,suffix_b)) &&
      (!nzchar(suffix_a) || nzchar(suffix_b))
  },TRUE)]
  if (length(compatible)==1L) return(list(status="mapped",allele=compatible[[1L]]))
  list(status=if (length(compatible)) "ambiguous" else "unresolved",allele=NULL)
}

hla_loh_number <- function(x,lower=0,upper=Inf) {
  if (is.null(x) || length(x)!=1L || is.na(x) || !nzchar(trimws(x)) ||
      trimws(x) %in% c("NA",".")) return(NULL)
  n <- suppressWarnings(as.numeric(x))
  if (!is.finite(n) || n<lower || n>upper) return(structure(list(raw=x),class="hla_loh_invalid"))
  n
}

hla_loh_compare <- function(value,threshold,operator) if (operator=="<")
  value<threshold else value<=threshold

hla_loh_interpret <- function(record,policy) {
  record$policy_version <- policy$policy_version
  record$policy_comparisons <- list(copy_number_loss_threshold=policy$copy_number_loss_threshold,
    copy_number_comparison=policy$copy_number_comparison,
    pval_unique_threshold=policy$pval_unique_threshold,
    pval_unique_comparison=policy$pval_unique_comparison)
  if (record$mapping_status!="mapped") {
    record$evaluation_status <- if (record$mapping_status=="ambiguous")
      "allele_mapping_ambiguous" else "not_evaluable"
    record$loh_call <- record$evaluation_status; record$loh_reason <- record$mapping_status
    return(record)
  }
  if (identical(record$allele_role,"homozygous")) {
    record$evaluation_status <- "homozygous_not_evaluable"
    record$loh_call <- "homozygous_not_evaluable"
    record$loh_reason <- "germline_homozygous_locus"
    return(record)
  }
  if (identical(record$allele_role,"provider_failed")) {
    record$evaluation_status <- "provider_failed"; record$loh_call <- "provider_failed"
    record$loh_reason <- "provider_run_failed"; return(record)
  }
  if (identical(record$allele_role,"unreported")) {
    record$evaluation_status <- "not_evaluable"; record$loh_call <- "not_evaluable"
    record$loh_reason <- "no_provider_row"; return(record)
  }
  if (any(vapply(record[c("allele_copy_number","pval_unique")],
    function(v) inherits(v,"hla_loh_invalid"),TRUE))) {
    record$evaluation_status <- "not_evaluable"; record$loh_call <- "not_evaluable"
    record$loh_reason <- "malformed_numeric"; return(record)
  }
  if (is.null(record$allele_copy_number) || is.null(record$pval_unique)) {
    record$evaluation_status <- "not_evaluable"; record$loh_call <- "not_evaluable"
    record$loh_reason <- "missing_required_metric"; return(record)
  }
  if (!record$allele_role %in% c("potential_loss","kept")) {
    record$evaluation_status <- "unresolved"; record$loh_call <- "unresolved"
    record$loh_reason <- "provider_allele_role_unresolved"; return(record)
  }
  record$evaluation_status <- "evaluated"
  low <- hla_loh_compare(record$allele_copy_number,
    policy$copy_number_loss_threshold,policy$copy_number_comparison)
  significant <- hla_loh_compare(record$pval_unique,
    policy$pval_unique_threshold,policy$pval_unique_comparison)
  record$policy_comparisons$copy_number_below_threshold <- low
  record$policy_comparisons$pval_below_threshold <- significant
  if (record$allele_role=="kept") {
    record$loh_call <- "loss_not_supported"; record$loh_reason <- "provider_kept_allele"
  } else if (low && significant) {
    record$loh_call <- "loss_supported"; record$loh_reason <- "cn_and_pval_support"
  } else if (!low && significant) {
    record$loh_call <- "allelic_imbalance_only"; record$loh_reason <- "cn_not_below_threshold"
  } else if (low && !significant) {
    record$loh_call <- "low_copy_without_significant_imbalance"
    record$loh_reason <- "pval_not_below_threshold"
  } else {
    record$loh_call <- "loss_not_supported"; record$loh_reason <- "neither_threshold_met"
  }
  record
}
