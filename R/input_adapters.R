# Explicit upstream schemas; no implicit count/TPM/FPKM conversion.
input_table <- function(path) {
  tab <- utils::read.delim(path,check.names=FALSE,colClasses="character",quote="\"",
    comment.char="#",na.strings=c("","NA","."),stringsAsFactors=FALSE)
  if (!nrow(tab) || anyDuplicated(names(tab))) stop("Empty/duplicate-column upstream table")
  tab
}
input_expression <- function(source,project) {
  format <- source$format
  if (is.null(format) || format=="generic") return(somatic_read_table(source))
  tab <- input_table(source$resolved_path)
  required <- switch(format,salmon=c("Name","Length","EffectiveLength","TPM","NumReads"),
    rsem_genes=c("gene_id","transcript_id(s)","TPM","FPKM","expected_count"),
    rsem_transcripts=c("transcript_id","gene_id","TPM","FPKM","expected_count"),
    featurecounts=c("Geneid","Chr","Start","End","Strand","Length"),
    stop("Unsupported expression format; choose salmon, rsem_genes, rsem_transcripts, featurecounts or generic"))
  if (!all(required %in% names(tab))) stop("Expression output does not match declared ",format," schema")
  tx <- switch(format,salmon=tab$Name,rsem_transcripts=tab$transcript_id,NULL)
  gene <- switch(format,rsem_genes=tab$gene_id,rsem_transcripts=tab$gene_id,featurecounts=tab$Geneid,NULL)
  count <- if (format=="featurecounts") {
    column <- source$options$column
    if (is.null(column)) column <- source$sample_id
    if (!column %in% setdiff(names(tab),required)) stop("featureCounts requires the exact declared sample count column")
    tab[[column]]
  } else if (format=="salmon") tab$NumReads else tab$expected_count
  lapply(seq_len(nrow(tab)),function(i) {
    value <- function(x) if(is.null(x) || is.na(x[[i]])) NULL else x[[i]]
    list(sample_id=source$sample_id,transcript_id=value(tx),gene_id=value(gene),
      gene_tpm=if (format=="rsem_genes") value(tab$TPM) else NULL,
      transcript_tpm=if (format %in% c("salmon","rsem_transcripts")) value(tab$TPM) else NULL,
      fpkm=value(tab$FPKM),estimated_count=value(count),quantification_method=format,
      expression_unit=if (format=="featurecounts") "count" else "TPM")
  })
}

# Native saveRDS(ascat.output) or the documented ASCAT_objects.Rdata workspace.
# Load workspaces into a private environment; never guess which object to use.
input_ascat_object <- function(path,format) {
  if (identical(format,"ascat_rds")) object <- readRDS(path) else if (identical(format,"ascat_rdata")) {
    env <- new.env(parent=emptyenv())
    names <- load(path,envir=env)
    if (!"ascat.output" %in% names) stop("ASCAT workspace requires exact object name 'ascat.output'")
    object <- get("ascat.output",envir=env,inherits=FALSE)
  } else stop("ASCAT metadata_format must be ascat_rds or ascat_rdata")
  if (!is.list(object) || !is.data.frame(object$segments))
    stop("Expected saved ASCAT output list with segments data.frame")
  object
}

input_ascat <- function(source,project) {
  if (!identical(source$reference_build,project$config$scientific$genome_assembly))
    stop("ASCAT reference_build mismatch; genome lift-over is not implicit")
  opts <- source$options
  strict_keys(if(is.null(opts)) list() else opts,c("coordinate_system","normal_cn","sample_name","metadata_path","metadata_format","metadata_checksum"),"ASCAT options")
  if (is.null(opts$coordinate_system) || opts$coordinate_system!="1-based-inclusive")
    stop("ASCAT requires explicit coordinate_system='1-based-inclusive' for SNP-position segments")
  format <- source$format
  if (is.null(format)) stop("ASCAT format must be ascat_segments, ascat_rds or ascat_rdata")
  object <- NULL
  if (format %in% c("ascat_rds","ascat_rdata")) {
    object <- input_ascat_object(source$resolved_path,format)
    tab <- object$segments
  } else if (format=="ascat_segments") tab <- input_table(source$resolved_path)
  else stop("Unsupported ASCAT schema; raw LogR/BAF is not absolute copy number")
  if (!is.null(opts$metadata_path)) {
    if (format!="ascat_segments") stop("ASCAT companion metadata is only for ascat_segments")
    if (is.null(source$metadata_path) || is.null(source$metadata_sha256))
      stop("ASCAT companion must be bound and checksummed before import")
    object <- input_ascat_object(source$metadata_path,opts$metadata_format)
  } else if (!is.null(opts$metadata_format) || !is.null(opts$metadata_checksum))
    stop("ASCAT metadata_format/checksum require metadata_path")
  required <- c("sample","chr","startpos","endpos","nMajor","nMinor")
  if (!all(required %in% names(tab))) stop("ASCAT segments need sample, chr, startpos, endpos, nMajor, nMinor")
  if (anyNA(tab[required])) stop("ASCAT segments contain missing values or conflicting sample identities")
  upstream_samples <- unique(as.character(tab$sample))
  expected <- if(is.null(opts$sample_name)) source$sample_id else opts$sample_name
  scalar_text(expected,"ASCAT sample_name")
  if (!expected %in% upstream_samples ||
      (is.null(opts$sample_name) && any(upstream_samples!=expected)))
    stop("ASCAT sample identities mismatch; declare a verified options$sample_name mapping explicitly")
  tab <- tab[as.character(tab$sample)==expected,,drop=FALSE]
  source$sample_mapping <- list(upstream_sample=expected,neo_sample=source$sample_id,
    policy=if(is.null(opts$sample_name)) "exact_identity" else "explicit_verified_declaration")
  if (!is.null(opts$metadata_path)) {
    # A companion from another fit must not supply purity/ploidy to these segments.
    other <- object$segments
    if (!all(required %in% names(other))) stop("ASCAT companion segments schema mismatch")
    other <- other[as.character(other$sample)==expected,,drop=FALSE]
    canonical <- function(x) {
      x <- x[required]
      x[] <- lapply(x,as.character)
      for (key in c("startpos","endpos","nMajor","nMinor"))
        x[[key]] <- suppressWarnings(as.numeric(x[[key]]))
      x <- x[do.call(order,x),,drop=FALSE];rownames(x) <- NULL;x
    }
    if (!identical(canonical(tab),canonical(other)))
      stop("ASCAT companion segments do not match the declared sample/fit")
  }
  metadata <- function(name,old=NULL) {
    value <- object[[name]]
    if (is.null(value) && !is.null(old)) value <- object[[old]]
    if (is.null(value)) return(NULL)
    metadata_samples <- unique(as.character(object$segments$sample))
    if (is.null(names(value)) && length(metadata_samples)>1L)
      stop("ASCAT multi-sample metadata must be named by upstream sample")
    if (!is.null(names(value))) {
      if (anyDuplicated(names(value)) || !expected %in% names(value)) stop("ASCAT metadata sample mismatch")
      value <- value[[expected]]
    }
    if (length(value)==1L && is.na(value)) return(NULL)
    somatic_number(value,paste("ASCAT",name),minimum=0,maximum=if(name=="purity") 1 else Inf)
  }
  purity <- metadata("purity","aberrantcellfraction");ploidy <- metadata("ploidy")
  normal <- somatic_number(opts$normal_cn,"explicit normal_cn",minimum=0)
  rows <- lapply(seq_len(nrow(tab)),function(i) {
    major <- somatic_number(tab$nMajor[[i]],"ASCAT nMajor",minimum=0)
    minor <- somatic_number(tab$nMinor[[i]],"ASCAT nMinor",minimum=0)
    if (major<minor || major!=floor(major) || minor!=floor(minor))
      stop("Rounded ASCAT nMajor/nMinor must be nonnegative integer copies with major >= minor")
    list(chrom=as.character(tab$chr[[i]]),loc.start=as.character(tab$startpos[[i]]),
      loc.end=as.character(tab$endpos[[i]]),tcn.em=as.character(major+minor),
      lcn.em=as.character(minor),normal_cn=normal,
      purity=if(is.null(purity)) NULL else as.character(purity),
      ploidy=if(is.null(ploidy)) NULL else as.character(ploidy),
      sample_id=source$sample_id,reference_build=source$reference_build)
  })
  source$original_kind <- "ascat";source$kind <- "facets_cn";source$adapter <- "ASCAT-segments-v1";source$rows <- rows
  segments <- somatic_facets_segments(source)
  for (i in seq_along(segments)) for (j in seq_len(i-1L))
    if (segments[[i]]$contig==segments[[j]]$contig &&
        segments[[i]]$start<=segments[[j]]$end && segments[[j]]$start<=segments[[i]]$end)
      stop("ASCAT segments overlap; no arbitrary segment selection")
  source
}

input_source_rows <- function(sources,events,project) {
  parsed <- list();bams <- list(); coordinate_ref <- NULL; coordinate_bundle <- NULL
  on.exit(if (!is.null(coordinate_ref)) event_reference_close(coordinate_ref), add=TRUE)
  coordinate_key <- function(sample,contig,pos,ref,alt)
    paste(sample,contig,pos,ref,alt,sep="\037")
  keys <- vapply(events,function(e) if(is.null(e$contig) || is.null(e$ref) || is.null(e$alt))
    paste0("non-coordinate:",e$event_id) else
    coordinate_key(e$sample_id,e$contig,e$pos,e$ref,e$alt),"")
  event_lookup <- external_group_index(keys)
  for (s in sources) {
    if (s$kind %in% c("dna_bam","normal_dna_bam","rna_bam")) { bams[[length(bams)+1L]] <- s;next }
    if (s$kind=="ascat") s <- input_ascat(s,project) else
      s$rows <- if(s$kind=="expression") input_expression(s,project) else somatic_read_table(s)
    if (s$kind %in% c("dna_reads","rna_reads","ccf","loh")) {
      for (j in seq_along(s$rows)) if (is.null(s$rows[[j]]$event_id)) {
        row <- s$rows[[j]]
        if (!all(c("contig","pos","ref","alt") %in% names(row)))
          stop("Evidence rows require event_id or contig,pos,ref,alt")
        if (!identical(s$reference_build,project$config$scientific$genome_assembly))
          stop("Coordinate evidence requires matching reference_build")
        if (is.null(coordinate_ref)) {
          coordinate_bundle <- reference_bundle(project)
          coordinate_ref <- event_reference_open(coordinate_bundle)
        }
        canonical <- event_normalize(coordinate_ref,row$contig,
          somatic_integer(row$pos,"evidence pos"),row$ref,row$alt)
        matches <- external_index_get(event_lookup,coordinate_key(s$sample_id,canonical$contig,
          canonical$pos,canonical$ref,canonical$alt))
        if(length(matches)!=1L) stop("Coordinate evidence must match exactly one sample Event")
        s$rows[[j]]$event_id <- events[[matches]]$event_id
      }
    }
    parsed[[length(parsed)+1L]] <- s
  }
  kinds <- vapply(parsed,function(s)paste(s$kind,s$sample_id),"")
  if (anyDuplicated(kinds)) stop("Choose one copy-number source per sample; ASCAT and FACETS cannot silently override")
  for (s in bams) {
    query_events <- events
    if (s$kind=="dna_bam" && input_bam_options(s$options)$mode=="fallback") {
      supplied <- somatic_source_index(parsed,"dna_reads",s$sample_id)
      query_events <- Filter(function(e) is.null(somatic_assay(e,
        somatic_match_event_row(supplied,e),"tumor",supplied)$vaf),events)
    }
    observations <- if(length(query_events)) input_bam_observations(s,query_events,project) else list()
    kind <- if(s$kind=="rna_bam") "rna_reads" else "dna_reads"
    ix <- which(vapply(parsed,function(x)x$kind==kind && x$sample_id==s$sample_id,TRUE))
    if (!length(ix)) {
      carrier <- s;carrier$kind <- kind;carrier$rows <- list();carrier$bam_observations <- list()
      parsed[[length(parsed)+1L]] <- carrier;ix <- length(parsed)
    }
    at <- ix[[1L]]
    for (id in names(observations)) parsed[[at]]$bam_observations[[id]] <-
      c(parsed[[at]]$bam_observations[[id]],observations[[id]])
  }
  # Resolve gene-only estimates through the reference-derived transcript models,
  # not guessed gene symbols, and retain both original IDs in the evidence.
  contexts <- if(file.exists(reconstruction_stage_path(project))) readRDS(reconstruction_stage_path(project)) else list()
  for (i in seq_along(parsed)) if(parsed[[i]]$kind=="expression") {
    expanded <- list()
    for (row in parsed[[i]]$rows) {
      if (is.null(row$transcript_id) && !is.null(row$gene_id)) {
        hits <- Filter(function(x)x$sample_id==parsed[[i]]$sample_id &&
          sub("\\.[0-9]+$","",x$gene_id)==sub("\\.[0-9]+$","",row$gene_id),contexts)
        txs <- unique(vapply(hits,`[[`,"","transcript_id"))
        if (length(txs)) for (tx in txs) {mapped <- row;mapped$transcript_id <- tx;expanded[[length(expanded)+1L]] <- mapped}
        else expanded[[length(expanded)+1L]] <- row
      } else expanded[[length(expanded)+1L]] <- row
    }
    parsed[[i]]$rows <- expanded
  }
  # Revalidate the bound reference at the operation boundary; reuse is local
  # to this call, never a persistent cache or a weaker checksum policy.
  if (!is.null(coordinate_bundle)) for (r in coordinate_bundle$resources) {
    current <- integrity_record(r$resolved_path,r$integrity$mode,r$integrity$declared)
    if (!identical(current,r$integrity)) stop("Reference changed during coordinate evidence import")
    if (!is.null(r$integrity$declared) && r$integrity$mode!="sha256")
      integrity_record(r$resolved_path,"sha256",r$integrity$declared)
  }
  parsed
}
