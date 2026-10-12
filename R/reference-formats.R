reference_fasta <- function(path, kind = c("transcript", "protein", "legacy")) {
  kind <- match.arg(kind)
  con <- reference_connection(path); on.exit(close(con))
  header_chunks <- list(); size_chunks <- list(); k <- 0L
  pending_header <- NULL; pending_size <- 0
  alphabet <- if (kind == "protein") "^[ABCDEFGHIKLMNPQRSTVWXYZUO*abcdefghiklmnpqrstvwxyzuo]+$" else
    "^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+$"
  repeat {
    lines <- readLines(con, n = 20000L, warn = FALSE)
    if (!length(lines)) break
    h <- startsWith(lines, ">")
    if (any(!grepl(alphabet, lines[!h]))) stop("Invalid FASTA sequence/record")
    at <- which(h)
    if (!length(at)) {
      if (is.null(pending_header)) stop("Sequence before FASTA header")
      pending_size <- pending_size + sum(nchar(lines, type = "bytes")); next
    }
    if (at[1] > 1L && is.null(pending_header)) stop("Sequence before FASTA header")
    if (!is.null(pending_header)) {
      size <- pending_size + sum(nchar(lines[seq_len(at[1] - 1L)], type = "bytes"))
      if (!size) stop("Empty FASTA record")
      k <- k + 1L; header_chunks[[k]] <- pending_header; size_chunks[[k]] <- size
    }
    cumulative <- c(0, cumsum(ifelse(h, 0, nchar(lines, type = "bytes"))))
    ends <- c(at[-1] - 1L, length(lines))
    sizes <- cumulative[ends + 1L] - cumulative[at + 1L]
    if (length(at) > 1L) {
      if (any(utils::head(sizes, -1L) == 0)) stop("Empty FASTA record")
      k <- k + 1L; header_chunks[[k]] <- substring(lines[utils::head(at, -1L)], 2)
      size_chunks[[k]] <- utils::head(sizes, -1L)
    }
    pending_header <- substring(lines[utils::tail(at, 1L)], 2); pending_size <- utils::tail(sizes, 1L)
  }
  if (is.null(pending_header) || !pending_size) stop("Empty FASTA")
  headers <- c(unlist(header_chunks, use.names = FALSE), pending_header)
  sizes <- c(unlist(size_chunks, use.names = FALSE), pending_size)
  token <- sub("[[:space:]].*$", "", headers)
  fields <- strsplit(token, "|", fixed = TRUE)
  ids <- if (kind == "legacy") token else vapply(fields, `[`, "", 1L)
  if (any(!nzchar(ids)) || anyDuplicated(ids)) stop("Duplicate/empty FASTA identifier")
  out <- data.frame(id = ids, raw_header = headers, length = sizes, stringsAsFactors = FALSE)
  if (kind != "legacy") {
    required <- if (kind == "protein") 8L else 7L
    if (any(lengths(fields) < required)) stop("Malformed GENCODE FASTA header")
    get <- function(i) vapply(fields, `[`, "", i)
    out$transcript_id <- if (kind == "protein") get(2L) else ids
    out$gene_id <- get(if (kind == "protein") 3L else 2L)
    if (kind == "protein") out$protein_id <- ids
    declared_length <- reference_integer(get(required), "FASTA declared length", 1)
    if (any(declared_length != sizes)) stop("FASTA declared length mismatch")
  }
  reference_sort(out)
}
reference_gtf <- function(path, genome = NULL, collect = c("transcripts", "all"), visitor = NULL) {
  collect <- match.arg(collect)
  con <- reference_connection(path); on.exit(close(con))
  chunks <- list(); evidence <- character(); k <- 0L
  feature_counts <- integer(); region_counts <- integer(); seqnames <- character(); associations <- list(); bounds <- list(); j <- 0L
  repeat {
    lines <- readLines(con, n = 10000L, warn = FALSE)
    if (!length(lines)) break
    comments <- startsWith(lines, "#")
    evidence <- c(evidence, lines[comments]); lines <- lines[!comments]
    if (!length(lines)) next
    fields <- stringi::stri_split_fixed(lines, "\t", omit_empty = FALSE)
    if (any(lengths(fields) != 9L)) stop("Malformed GTF field count")
    x <- as.data.frame(do.call(rbind, fields), stringsAsFactors = FALSE)
    names(x) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "phase", "raw_attributes")
    x$start <- reference_integer(x$start, "GTF start", 1)
    x$end <- reference_integer(x$end, "GTF end", 1)
    if (any(x$start > x$end)) stop("Invalid GTF coordinate interval")
    if (any(!x$strand %in% c("+", "-"))) stop("Invalid GTF strand")
    if (any(!x$phase %in% c(".", "0", "1", "2")) || any(x$feature == "CDS" & x$phase == "."))
      stop("Invalid GTF phase")
    if (!is.null(genome)) {
      idx <- match(x$seqname, genome$contig)
      if (anyNA(idx)) stop("GTF contig absent from genome (aliases are not repaired): ",
        paste(unique(x$seqname[is.na(idx)]), collapse = ", "))
      if (any(x$end > genome$length[idx])) stop("GTF coordinate exceeds genome")
    }
    # GENCODE quoted attributes and MANE's unquoted integer exon_number.
    pattern <- '[A-Za-z_][A-Za-z0-9_]*[[:space:]]+("[^"\\r\\n]*"|[0-9]+)[[:space:]]*;'
    residue <- trimws(gsub(pattern, "", x$raw_attributes, perl = TRUE))
    if (any(nzchar(residue))) stop("Malformed GTF attributes")
    keys <- c("gene_id", "transcript_id", "protein_id", "gene_type", "transcript_type",
      "gene_biotype", "transcript_biotype", "exon_id", "exon_number")
    for (key in keys) {
      pat <- paste0('(?:^|(?<=;))[[:space:]]*', key, '[[:space:]]+(?:"([^"]+)"|([0-9]+));')
      if (any(stringi::stri_count_regex(x$raw_attributes, pat) > 1L))
        stop("Duplicate GTF attribute: ", key)
      matches <- stringi::stri_match_first_regex(x$raw_attributes, pat)
      value <- ifelse(is.na(matches[, 2]), matches[, 3], matches[, 2])
      x[[key]] <- value
    }
    if (anyNA(x$gene_id) || any(is.na(x$transcript_id) & x$feature != "gene")) stop("Missing GTF identity")
    if (!is.null(visitor)) visitor(x)
    counts <- table(x$feature)
    for (f in names(counts)) feature_counts[f] <- as.integer(counts[f]) +
      if (f %in% names(feature_counts)) feature_counts[[f]] else 0L
    counts <- table(x$seqname)
    for (f in names(counts)) region_counts[f] <- as.integer(counts[f]) +
      if (f %in% names(region_counts)) region_counts[[f]] else 0L
    seqnames <- union(seqnames, x$seqname)
    j <- j + 1L
    z <- x[!is.na(x$transcript_id), c("transcript_id", "start", "end")]
    low <- z[order(z$transcript_id, z$start, method = "radix"), ]
    low <- low[!duplicated(low$transcript_id), ]
    high <- z[order(z$transcript_id, -z$end, method = "radix"), ]
    high <- high[!duplicated(high$transcript_id), ]
    bounds[[j]] <- data.frame(transcript_id = low$transcript_id, start = low$start, end = high$end)
    associations[[j]] <- unique(x[!is.na(x$transcript_id),
      c("transcript_id", "gene_id", "seqname", "strand", "protein_id")])
    if (collect == "transcripts") {
      # Raw attributes remain in the bound source and are exposed to the
      # visitor before this compaction. Do not retain hundreds of MB of
      # repeated attribute syntax in the transcript lookup index.
      x <- x[x$feature == "transcript", setdiff(names(x), "raw_attributes"), drop = FALSE]
    }
    if (nrow(x)) { k <- k + 1L; chunks[[k]] <- x }
  }
  if (!length(chunks)) stop("Empty GTF")
  x <- do.call(rbind, chunks)
  associations <- unique(do.call(rbind, associations))
  tx <- unique(associations[c("transcript_id", "gene_id", "seqname", "strand")])
  if (anyDuplicated(tx$transcript_id)) stop("Contradictory GTF transcript identity")
  links <- unique(associations[!is.na(associations$protein_id), c("transcript_id", "protein_id")])
  if (anyDuplicated(links$transcript_id)) stop("Contradictory GTF protein relationship")
  declared <- x[x$feature == "transcript", , drop = FALSE]
  if (anyDuplicated(declared$transcript_id) || any(!tx$transcript_id %in% declared$transcript_id))
    stop("Missing/duplicate GTF transcript declaration")
  bounds <- do.call(rbind, bounds)
  idx <- match(bounds$transcript_id, declared$transcript_id)
  if (any(bounds$start < declared$start[idx] | bounds$end > declared$end[idx]))
    stop("GTF feature outside transcript bounds")
  if (collect == "transcripts") {
    linked <- links$protein_id[match(x$transcript_id, links$transcript_id)]
    x$protein_id[is.na(x$protein_id)] <- linked[is.na(x$protein_id)]
  }
  structure(reference_sort(x), source_headers = sort(unique(evidence), method = "radix"),
    feature_counts = feature_counts[order(names(feature_counts), method = "radix")],
    seqnames = sort(seqnames, method = "radix"), source_path = path,
    region_feature_counts = region_counts[order(names(region_counts), method = "radix")])
}

# Explicit on-demand feature access preserves raw source attributes, phase,
# exon/CDS/start/stop/Sec rows without embedding the full GTF in the project.
reference_transcript_features <- function(bundle, transcript_ids, source = c("annotation", "mane")) {
  source <- match.arg(source)
  if (!inherits(bundle, "neo_reference_bundle") || !is.character(transcript_ids) || anyNA(transcript_ids))
    stop("Invalid reference feature query")
  role <- if (source == "annotation") "annotation_gtf" else "mane_gtf"
  if (!role %in% names(bundle$resources)) stop("GTF source unavailable")
  parts <- list(); k <- 0L
  r <- bundle$resources[[role]]
  if (!identical(integrity_record(r$resolved_path, r$integrity$mode, r$integrity$declared), r$integrity))
    stop("Reference changed since bundle validation")
  reference_gtf(r$resolved_path, visitor = function(x) {
    x <- x[!is.na(x$transcript_id) & x$transcript_id %in% transcript_ids, , drop = FALSE]
    if (nrow(x)) { k <<- k + 1L; parts[[k]] <<- x }
  })
  if (!length(parts)) return(data.frame())
  reference_sort(do.call(rbind, parts))
}
