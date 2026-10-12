# T008 edits are expressed against the unchanged, forward-genomic reference.
reconstruction_edit_core <- function(row) {
  ref <- row$ref; alt <- row$alt
  if (is.null(ref) || is.null(alt) || !grepl("^[ACGT]+$", ref) ||
      !grepl("^[ACGT]+$", alt)) return(NULL)
  a <- strsplit(ref, "", fixed = TRUE)[[1L]]
  b <- strsplit(alt, "", fixed = TRUE)[[1L]]
  left <- 0L
  while (left < min(length(a), length(b)) && identical(a[left + 1L], b[left + 1L]))
    left <- left + 1L
  right <- 0L
  while (right < min(length(a), length(b)) - left &&
         identical(a[length(a) - right], b[length(b) - right])) right <- right + 1L
  old <- if (left + right == length(a)) "" else
    paste(a[seq.int(left + 1L, length(a) - right)], collapse = "")
  new <- if (left + right == length(b)) "" else
    paste(b[seq.int(left + 1L, length(b) - right)], collapse = "")
  list(event_id = row$event_id, start = as.integer(row$pos + left),
       end = as.integer(row$pos + length(a) - right - 1L),
       old = old, new = new, insertion = !nzchar(old),
       delta = nchar(new) - nchar(old))
}

reconstruction_reverse_complement <- function(x) {
  if (!nzchar(x)) return("")
  paste(rev(chartr("ACGT", "TGCA", strsplit(x, "", fixed = TRUE)[[1L]])), collapse = "")
}

# Segments are in increasing genomic order. An edit must fit wholly in one
# segment; an insertion needs both flanking coding bases in that segment.
reconstruction_project_edit <- function(edit, segments) {
  if (is.null(edit)) return(list(status = "UNSUPPORTED_ALLELE"))
  hit <- if (edit$insertion)
    which(segments$start <= edit$start - 1L & segments$end >= edit$start)
  else which(segments$start <= edit$start & segments$end >= edit$end)
  if (length(hit) != 1L) return(list(status = "CDS_BOUNDARY_OR_SPLICE_INTERSECTION"))
  edit$segment <- hit
  edit$offset <- as.integer(edit$start - segments$start[hit])
  edit$status <- "PROJECTED"
  edit
}

# Every edit is validated against the original segment. Compatible duplicates
# collapse; intersecting assertions otherwise defer rather than taking row order.
reconstruction_apply_edits <- function(sequence, edits) {
  if (!length(edits)) return(list(status = "COMPLETE", sequence = sequence))
  atoms <- list()
  for (e in edits) {
    if (nchar(e$old) == nchar(e$new) && nchar(e$old) > 1L) {
      old <- strsplit(e$old, "", fixed = TRUE)[[1L]]
      new <- strsplit(e$new, "", fixed = TRUE)[[1L]]
      for (j in which(old != new)) {
        atom <- e; atom$offset <- e$offset + j - 1L
        atom$old <- old[[j]]; atom$new <- new[[j]]
        atoms[[length(atoms) + 1L]] <- atom
      }
    } else atoms[[length(atoms) + 1L]] <- e
  }
  edits <- atoms
  if (!length(edits)) return(list(status = "COMPLETE", sequence = sequence))
  edits <- edits[order(vapply(edits, `[[`, 0L, "offset"),
                       vapply(edits, function(e) as.integer(!e$insertion), 0L),
                       vapply(edits, `[[`, 0L, "end"),
                       vapply(edits, `[[`, "", "event_id"), method = "radix")]
  keep <- list()
  for (e in edits) {
    observed <- if (e$insertion) "" else substr(sequence, e$offset + 1L,
                                                   e$offset + nchar(e$old))
    if (!identical(observed, e$old)) return(list(status = "REFERENCE_ALLELE_MISMATCH"))
    if (length(keep)) {
      p <- keep[[length(keep)]]
      duplicate <- identical(c(e$offset, e$old, e$new), c(p$offset, p$old, p$new))
      if (duplicate) {
        keep[[length(keep)]]$event_ids <- sort(unique(c(p$event_ids, e$event_id)), method = "radix")
        next
      }
      # A gap at the edge of another edit is unambiguous; a gap within an
      # altered interval, two different insertions at one gap, or shared bases is not.
      p_end <- p$offset + nchar(p$old)
      if (e$offset < p_end || (e$offset == p$offset && e$insertion && p$insertion))
        return(list(status = "OVERLAPPING_EDIT_CONFLICT"))
    }
    e$event_ids <- e$event_id
    keep[[length(keep) + 1L]] <- e
  }
  out <- character(); at <- 1L
  for (e in keep) {
    if (e$offset >= at) out <- c(out, substr(sequence, at, e$offset))
    out <- c(out, e$new)
    at <- e$offset + nchar(e$old) + 1L
  }
  if (at <= nchar(sequence)) out <- c(out, substr(sequence, at, nchar(sequence)))
  list(status = "COMPLETE", sequence = paste0(out, collapse = ""),
       event_ids = sort(unique(unlist(lapply(keep, `[[`, "event_ids"))), method = "radix"))
}
