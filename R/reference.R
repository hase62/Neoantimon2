# T005 is an internal, read-only reference validator. No stage or public API.
reference_scope_policy <- function() "gencode_primary_assembly_all_catalogs_mane_supplemental_v1"
reference_roles <- function(family) {
  if (family == "legacy") return(c("refFlat", "refMrna"))
  c("genome_fasta", "genome_fai", "annotation_gtf", "transcript_fasta",
    "protein_coding_transcript_fasta", "protein_fasta", "selenocysteine_metadata",
    "mane_summary", "mane_gtf")
}
reference_sort <- function(x) {
  if (!nrow(x)) return(x)
  keys <- x[!vapply(x, is.list, TRUE)]
  x <- x[do.call(order, c(unname(keys), list(method = "radix", na.last = TRUE))), , drop = FALSE]
  rownames(x) <- NULL
  x
}
reference_connection <- function(path) {
  con <- file(path, "rb"); magic <- readBin(con, "raw", 2L); close(con)
  compressed <- identical(magic, as.raw(c(31, 139)))
  if (grepl("[.]gz$", path) != compressed) stop("Compression/extension mismatch: ", basename(path))
  if (compressed) gzfile(path, "rt") else file(path, "rt")
}
reference_lines <- function(path) {
  con <- reference_connection(path); on.exit(close(con))
  readLines(con, warn = FALSE)
}
reference_integer <- function(x, label, minimum = 0) {
  if (anyNA(x) || any(!grepl("^[0-9]+$", x))) stop("Invalid integer: ", label)
  y <- as.numeric(x)
  if (any(!is.finite(y) | y < minimum | y > 2^53 - 1)) stop("Invalid integer: ", label)
  y
}
reference_table <- function(path, n, label) {
  lines <- reference_lines(path)
  if (!length(lines) || any(!nzchar(lines))) stop("Empty/malformed ", label)
  rows <- stringi::stri_split_fixed(lines, "\t", omit_empty = FALSE)
  if (any(lengths(rows) != n)) stop("Wrong field count in ", label)
  as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
}
reference_fai <- function(path, fasta) {
  x <- reference_table(path, 5L, "FAI")
  names(x) <- c("contig", "length", "offset", "line_bases", "line_width")
  if (any(!nzchar(x$contig)) || anyDuplicated(x$contig) || any(grepl("[[:space:]]", x$contig)))
    stop("Invalid/duplicate FAI contig")
  for (n in names(x)[-1]) x[[n]] <- reference_integer(x[[n]], paste("FAI", n), 1)
  if (any(x$line_bases > x$length | x$line_width <= x$line_bases |
          x$line_width > x$line_bases + 2)) stop("Invalid FAI line widths")
  # Seek through the uncompressed FASTA using FAI geometry; never load the genome.
  if (grepl("[.]gz$", fasta)) stop("Genome FASTA must be uncompressed for FAI seek validation")
  con <- file(fasta, "rb"); on.exit(close(con))
  size <- file.info(fasta)$size
  for (i in seq_len(nrow(x))) {
    r <- x[i, ]; boundary <- r$offset + floor((r$length - 1) / r$line_bases) *
      r$line_width + (r$length - 1) %% r$line_bases
    if (boundary >= size) stop("FAI extends beyond genome FASTA")
    seek(con, max(0, r$offset - 4096), origin = "start")
    prefix <- rawToChar(readBin(con, "raw", min(4096, r$offset)))
    lines <- strsplit(prefix, "\n", fixed = TRUE)[[1]]
    header <- utils::tail(lines[nzchar(lines)], 1L)
    if (!length(header) || !startsWith(header, ">") ||
        sub("[[:space:]].*$", "", substring(header, 2)) != r$contig)
      stop("FAI/genome header mismatch")
    seek(con, r$offset, origin = "start")
    first <- rawToChar(readBin(con, "raw", min(r$line_bases, 4096)))
    if (!grepl("^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+$", first)) stop("Invalid genome FASTA sequence")
    seek(con, r$offset + r$line_bases, origin = "start")
    eol <- rawToChar(readBin(con, "raw", r$line_width - r$line_bases))
    if (!eol %in% c("\n", "\r\n")) stop("FAI/genome line-width mismatch")
    seek(con, boundary, origin = "start")
    last <- rawToChar(readBin(con, "raw", 1L))
    if (!grepl("^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]$", last)) stop("FAI/genome boundary mismatch")
  }
  reference_sort(x)
}
reference_compare_dictionary <- function(supplied, genome) {
  for (x in list(supplied, genome)) {
    if (!is.data.frame(x) || !all(c("contig", "length") %in% names(x)) ||
        anyNA(x[c("contig", "length")]) || anyDuplicated(x$contig) ||
        any(!nzchar(x$contig)) || !is.numeric(x$length) ||
        any(!is.finite(x$length) | x$length <= 0 | x$length != floor(x$length)))
      stop("Invalid sequence dictionary")
  }
  idx <- match(supplied$contig, genome$contig)
  missing <- supplied$contig[is.na(idx)]
  lengths <- supplied$contig[which(!is.na(idx) & supplied$length != genome$length[idx])]
  alias <- function(x) sub("^MT$", "M", sub("^chr", "", x))
  aliases <- unlist(lapply(missing, function(n) {
    hits <- genome$contig[alias(genome$contig) == alias(n)]
    if (length(hits)) paste(n, hits, sep = " -> ") else character()
  }), use.names = FALSE)
  missing_from_supplied <- setdiff(genome$contig, supplied$contig)
  list(compatible = !length(missing) && !length(missing_from_supplied) && !length(lengths),
       missing_contigs = missing, missing_from_supplied = missing_from_supplied,
       length_mismatches = lengths, alias_candidates = aliases)
}
reference_bundle <- function(project, progress = function(stage) invisible(NULL), retain_indexes = FALSE) {
  if (!is.logical(retain_indexes) || length(retain_indexes) != 1L || is.na(retain_indexes))
    stop("Invalid retain_indexes")
  if (!inherits(project, "neo_project")) stop("Expected neo_project")
  records <- project$resources
  logical <- vapply(records, `[[`, "", "logical_name")
  modern <- paste0("reference.", reference_roles("gencode"))
  legacy <- paste0("reference.", reference_roles("legacy"))
  if (any(logical %in% modern) && any(logical %in% legacy)) stop("Mixed modern/legacy reference is forbidden")
  family <- if (any(logical %in% legacy)) "legacy" else "gencode"
  roles <- reference_roles(family); wanted <- paste0("reference.", roles)
  if (anyDuplicated(logical) || any(!wanted %in% logical)) stop("Missing/duplicate reference role")
  records <- records[match(wanted, logical)]; names(records) <- roles
  states <- character(length(roles)); names(states) <- roles
  for (role in roles) {
    r <- records[[role]]; p <- r$resolved_path
    if (is.null(p) || !file.exists(p) || dir.exists(p) || r$type != "file" || file.access(p, 4) != 0)
      stop("Missing/unreadable reference file: ", role)
    if (r$integrity$mode == "none" && is.null(r$integrity$declared))
      stop("Reference needs metadata identity or a bound content checksum")
    # Reobserve using exactly T003's policy and reject a stale project snapshot.
    current <- integrity_record(p, r$integrity$mode, r$integrity$declared)
    if (!identical(current, r$integrity)) stop("Reference changed since project binding: ", role)
    if (!is.null(r$integrity$declared)) {
      if (r$integrity$mode != "sha256") integrity_record(p, "sha256", r$integrity$declared)
      states[[role]] <- "verified_sha256"
    } else states[[role]] <- if (r$integrity$mode == "sha256") "observed_sha256" else
      if (r$integrity$mode == "metadata") "metadata_only_not_content_verified" else "unverified"
  }
  if (family == "gencode" && !identical(records$annotation_gtf$resource_id, reference_scope_policy()))
    stop("Missing/unsupported reference scope policy in annotation resource_id")
  paths <- vapply(records, `[[`, "", "resolved_path")
  b <- if (family == "legacy") reference_legacy(paths) else reference_gencode(paths, records, project$config$scientific, progress)
  b$family <- family; b$resources <- records; b$analysis_id <- project$analysis_id
  b$resource_identity <- identity_records(unname(records)); b$checksum_state <- states
  b$declared <- list(species = project$config$scientific$species,
    assembly = project$config$scientific$genome_assembly,
    versions = lapply(records, `[[`, "version"), sources = lapply(records, `[[`, "source"))
  class(b) <- c("neo_reference_bundle", "list")
  b$indexes_materialized <- TRUE
  if (family == "gencode" && !retain_indexes) b <- reference_bundle_metadata(b)
  b
}
