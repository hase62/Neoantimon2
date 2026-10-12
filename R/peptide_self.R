peptide_self_search <- function(sequences, path = NULL, reference_identity = NULL,
                                sec_encoding = "UNDECLARED") {
  sequences <- sort(unique(sequences), method = "radix")
  for (s in sequences) peptide_sequence_validate(s)
  available <- !is.null(path) && file.exists(path)
  checksum <- if (available && length(sequences)) digest::digest(file = path,
    algo = "sha256", serialize = FALSE) else NULL
  out <- stats::setNames(lapply(sequences, function(s) list(
    status = if (!available) "UNAVAILABLE" else if (grepl("U", s, fixed = TRUE) &&
      sec_encoding != "LITERAL_U") "UNRESOLVED_RESIDUE_REPRESENTATION" else "EVALUATED",
    exact_match = NULL, matches = list(), reference_identity = reference_identity,
    comparator = "BOUND_REFERENCE_PROTEOME_ALL_ISOFORMS", reference_checksum = checksum,
    sec_encoding = sec_encoding)), sequences)
  keys <- sequences[vapply(out, function(x) x$status == "EVALUATED", TRUE)]
  if (!length(keys)) return(out)
  con <- reference_connection(path); on.exit(close(con))
  ids <- character(); proteins <- character(); headers <- character()
  current <- NULL; pieces <- character()
  scan_batch <- function() {
    if (!length(ids)) return()
    for (key in keys) {
      hits <- which(grepl(key, proteins, fixed = TRUE))
      for (i in hits) {
        # Lookahead reports overlapping exact occurrences (e.g. homopolymers).
        positions <- gregexpr(paste0("(?=", key, ")"), proteins[[i]], perl = TRUE)[[1L]]
        for (p in positions[positions > 0L]) out[[key]]$matches[[length(out[[key]]$matches) + 1L]] <<-
          list(protein_id = ids[[i]], header = headers[[i]], start_aa = as.integer(p),
            end_aa = as.integer(p + nchar(key) - 1L))
      }
    }
    ids <<- character(); proteins <<- character(); headers <<- character()
  }
  flush <- function() {
    if (is.null(current)) return()
    ids <<- c(ids, strsplit(current, "|", fixed = TRUE)[[1L]][[1L]])
    headers <<- c(headers, current)
    proteins <<- c(proteins, paste0(pieces, collapse = ""))
    if (length(ids) >= 2000L) scan_batch()
  }
  repeat {
    lines <- readLines(con, n = 20000L, warn = FALSE)
    if (!length(lines)) break
    for (line in lines) if (startsWith(line, ">")) {
      flush(); current <- substring(line, 2L); pieces <- character()
    } else if (!is.null(current)) pieces <- c(pieces, line)
  }
  flush(); scan_batch()
  for (key in keys) {
    out[[key]]$exact_match <- length(out[[key]]$matches) > 0L
    if (length(out[[key]]$matches)) out[[key]]$matches <- sorted_records(out[[key]]$matches)
  }
  out
}
