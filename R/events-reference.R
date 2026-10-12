# All coordinates in this module are 1-based, inclusive genomic coordinates.
event_reference_open <- function(bundle) {
  if (!inherits(bundle, "neo_reference_bundle") || bundle$family != "gencode" ||
      bundle$profile$coordinate_scope != "primary_assembly")
    stop("Event import requires the selected T005 primary genome reference")
  fasta <- bundle$resources$genome_fasta$resolved_path
  fai <- bundle$genome
  con <- file(fasta, "rb")
  list(con = con, fai = fai, assembly = bundle$profile$assembly_release)
}

event_reference_close <- function(reference) close(reference$con)

event_reference_fetch <- function(reference, contig, start, end) {
  if (!is.character(contig) || length(contig) != 1L || is.na(contig)) stop("Invalid contig")
  i <- match(contig, reference$fai$contig)
  if (is.na(i)) stop("Contig absent from selected primary reference: ", contig)
  r <- reference$fai[i, ]
  if (length(start) != 1L || length(end) != 1L || !is.numeric(start) || !is.numeric(end) ||
      anyNA(c(start, end)) || any(!is.finite(c(start, end))) || start != floor(start) ||
      end != floor(end) || start < 1 || end < start || end > r$length)
    stop("Event interval outside selected reference: ", contig, ":", start, "-", end)
  # FAI offsets are zero-based byte offsets; the caller API remains 1-based.
  byte <- function(p) r$offset + floor((p - 1) / r$line_bases) * r$line_width +
    (p - 1) %% r$line_bases
  seek(reference$con, byte(start), origin = "start")
  raw <- readBin(reference$con, "raw", n = as.integer(byte(end) - byte(start) + 1))
  seq <- gsub("[\r\n]", "", rawToChar(raw))
  if (nchar(seq, type = "bytes") != end - start + 1)
    stop("FAI seek returned incorrect sequence length")
  toupper(seq)
}

event_validate_ref <- function(reference, contig, pos, ref) {
  observed <- event_reference_fetch(reference, contig, pos, pos + nchar(ref) - 1L)
  if (!identical(observed, ref))
    stop("REF mismatch at ", contig, ":", pos, ": source ", ref,
         ", selected reference ", observed, call. = FALSE)
  invisible(TRUE)
}

event_normalize <- function(reference, contig, pos, ref, alt) {
  if (!is.character(ref) || !is.character(alt) || length(ref) != 1L || length(alt) != 1L ||
      anyNA(c(ref, alt)) || !grepl("^[ACGT]+$", ref) || !grepl("^[ACGT]+$", alt))
    stop("Malformed sequence allele; expected nonempty forward A/C/G/T REF and ALT")
  if (identical(ref, alt)) stop("No-op REF=ALT allele")
  event_validate_ref(reference, contig, pos, ref)
  original <- list(pos = pos, ref = ref, alt = alt)
  # VCF-style minimalization retains one anchor base for unequal lengths.
  trim <- function() {
    while (nchar(ref) > 1L && nchar(alt) > 1L &&
           substr(ref, nchar(ref), nchar(ref)) == substr(alt, nchar(alt), nchar(alt))) {
      ref <<- substr(ref, 1L, nchar(ref) - 1L)
      alt <<- substr(alt, 1L, nchar(alt) - 1L)
    }
    while (nchar(ref) > 1L && nchar(alt) > 1L && substr(ref, 1L, 1L) == substr(alt, 1L, 1L)) {
      ref <<- substring(ref, 2L); alt <<- substring(alt, 2L); pos <<- pos + 1L
    }
  }
  trim()
  if (nchar(ref) != nchar(alt)) {
    repeat {
      if (pos <= 1L || substr(ref, nchar(ref), nchar(ref)) !=
          substr(alt, nchar(alt), nchar(alt))) break
      prev <- event_reference_fetch(reference, contig, pos - 1L, pos - 1L)
      ref <- paste0(prev, substr(ref, 1L, nchar(ref) - 1L))
      alt <- paste0(prev, substr(alt, 1L, nchar(alt) - 1L))
      pos <- pos - 1L
      trim()
    }
  }
  event_validate_ref(reference, contig, pos, ref)
  type <- if (nchar(ref) == nchar(alt)) {
    if (nchar(ref) == 1L) "SNV" else "MNV"
  } else if (nchar(ref) == 1L && startsWith(alt, ref)) "INS"
  else if (nchar(alt) == 1L && startsWith(ref, alt)) "DEL"
  else "INDEL"
  list(contig = contig, pos = as.integer(pos), ref = ref, alt = alt,
       end = as.integer(pos + nchar(ref) - 1L), event_type = type,
       changed = !identical(original, list(pos = pos, ref = ref, alt = alt)))
}
