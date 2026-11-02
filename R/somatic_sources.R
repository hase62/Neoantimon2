somatic_source_files <- function(project) {
  policy <- project$config$somatic_evidence
  if (is.null(policy)) return(list())
  samples <- vapply(project$samples,`[[`,"","sample_id")
  lapply(policy$sources,function(s) {
    if (!s$sample_id %in% samples) stop("Evidence source sample is not in project")
    path <- resolve_path(s$path,project$config$paths$input_root)
    if (!file.exists(path) || dir.exists(path)) stop("Evidence source file missing")
    sha <- digest::digest(file=path,algo="sha256",serialize=FALSE)
    if (!is.null(s$checksum) && tolower(s$checksum)!=sha)
      stop("Evidence source checksum mismatch")
    s$resolved_path <- path
    s$sha256 <- sha
    if (s$kind %in% c("dna_bam","normal_dna_bam","rna_bam")) {
      indexes <- c(paste0(path,".bai"),sub("[.]bam$",".bai",path),paste0(path,".csi"))
      indexes <- unique(indexes[file.exists(indexes)])
      if (length(indexes)!=1L) stop("BAM requires exactly one existing unambiguous BAI/CSI index; index it externally")
      s$index_path <- indexes[[1L]]
      s$index_sha256 <- digest::digest(file=s$index_path,algo="sha256",serialize=FALSE)
    }
    if (s$kind=="ascat" && !is.null(s$options$metadata_path)) {
      if (!identical(s$format,"ascat_segments")) stop("ASCAT companion metadata is only for ascat_segments")
      if (length(s$options$metadata_format)!=1L ||
          !s$options$metadata_format %in% c("ascat_rds","ascat_rdata")) stop("Declare ASCAT metadata_format explicitly")
      scalar_text(s$options$metadata_path,"ASCAT metadata_path")
      s$metadata_path <- resolve_path(s$options$metadata_path,project$config$paths$input_root)
      if (!file.exists(s$metadata_path) || dir.exists(s$metadata_path)) stop("ASCAT companion file missing")
      s$metadata_sha256 <- digest::digest(file=s$metadata_path,algo="sha256",serialize=FALSE)
      if (!is.null(s$options$metadata_checksum) &&
          !identical(tolower(s$options$metadata_checksum),s$metadata_sha256)) stop("ASCAT companion checksum mismatch")
    }
    s$source_id <- paste0("somatic-source:",infra_hash(list(kind=s$kind,
      sample_id=s$sample_id,source=s$source,version=s$source_version,
      reference_build=s$reference_build,sha256=sha)))
    if (!is.null(s$metadata_sha256)) s$source_id <- paste0("somatic-source:",
      infra_hash(list(primary_source_id=s$source_id,metadata_sha256=s$metadata_sha256)))
    s
  })
}

somatic_read_table <- function(source) {
  reader <- if (grepl("\\.csv$",source$resolved_path,ignore.case=TRUE)) utils::read.csv else
    utils::read.delim
  table <- tryCatch(reader(source$resolved_path,check.names=FALSE,
    colClasses="character",na.strings=c("","NA","."),stringsAsFactors=FALSE),
    error=function(e) stop("Malformed somatic evidence table: ",conditionMessage(e)))
  if (anyDuplicated(names(table)) || any(!nzchar(names(table))))
    stop("Duplicate or empty somatic evidence column")
  lapply(seq_len(nrow(table)),function(i) {
    row <- as.list(table[i,,drop=FALSE])
    lapply(row,function(x) if (length(x)!=1L || is.na(x)) NULL else x)
  })
}

somatic_source_rows <- function(sources) {
  lapply(sources,function(s) {
    s$rows <- somatic_read_table(s)
    s
  })
}

somatic_integer <- function(value,label) {
  if (is.null(value)) return(NULL)
  if (is.numeric(value) && length(value)==1L && is.finite(value) &&
      value>=0 && value==floor(value) && value<=.Machine$integer.max) return(as.integer(value))
  if (is.character(value) && length(value)==1L && grepl("^[0-9]+$",value)) {
    number <- suppressWarnings(as.integer(value))
    if (!is.na(number)) return(number)
  }
  stop("Invalid nonnegative integer in ",label)
}

somatic_number <- function(value,label,minimum=-Inf,maximum=Inf) {
  if (is.null(value)) return(NULL)
  number <- suppressWarnings(as.numeric(value))
  if (length(number)!=1L || !is.finite(number) || number<minimum || number>maximum)
    stop("Invalid numeric value in ",label)
  number
}

somatic_row_value <- function(row,name) {
  value <- row[[name]]
  if (is.null(value) || length(value)!=1L || is.na(value)) NULL else value
}

somatic_source_provenance <- function(source) {
  if (is.null(source)) return(NULL)
  result <- list(source_id=source$source_id,kind=if(is.null(source$original_kind)) source$kind else source$original_kind,source=source$source,
    source_version=source$source_version,artifact=source$resolved_path,
    sha256=source$sha256,reference_build=source$reference_build,
    sample_id=source$sample_id)
  if (!is.null(source$metadata_sha256)) result$companion <- list(artifact=source$metadata_path,
    sha256=source$metadata_sha256,format=source$options$metadata_format)
  if (!is.null(source$sample_mapping)) result$sample_mapping <- source$sample_mapping
  result
}

somatic_source_index <- function(sources,kind,sample_id) {
  found <- Filter(function(x) x$kind==kind && x$sample_id==sample_id,sources)
  if (!length(found)) NULL else found[[1L]]
}

somatic_match_event_row <- function(source,event) {
  if (is.null(source)) return(NULL)
  rows <- Filter(function(x) identical(x$event_id,event$event_id),source$rows)
  if (length(rows)>1L) stop("Duplicate evidence event_id")
  if (!length(rows)) return(NULL)
  row <- rows[[1L]]
  if (!is.null(row$sample_id) && row$sample_id!=event$sample_id)
    stop("Evidence row sample/event mismatch")
  row
}
