# Stream plain list containers while delegating every leaf's JSON semantics to
# jsonlite. This changes whitespace for large files, never scientific values.
json_write_chunks <- function(value, con, leaf_bytes = 8 * 1024^2) {
  emit <- function(text) writeBin(charToRaw(enc2utf8(text)), con)
  plain <- is.list(value) && !is.object(value) && is.null(dim(value))
  if (!plain || !length(value) || as.numeric(utils::object.size(value)) <= leaf_bytes) {
    emit(as.character(jsonlite::toJSON(value, auto_unbox = TRUE, null = "null",
      digits = NA, pretty = TRUE)))
    return(invisible(NULL))
  }
  labels <- names(value)
  object <- !is.null(labels)
  emit(if (object) "{\n" else "[\n")
  for (i in seq_along(value)) {
    if (i > 1L) emit(",\n")
    if (object) {
      label <- labels[[i]]
      # jsonlite uses one-based positions for unnamed members of named lists.
      if (is.na(label) || !nzchar(label)) label <- as.character(i)
      emit(as.character(jsonlite::toJSON(label, auto_unbox = TRUE)))
      emit(": ")
    }
    json_write_chunks(value[[i]], con, leaf_bytes)
  }
  emit(if (object) "\n}" else "\n]")
  invisible(NULL)
}

# Equivalent to the existing jsonlite round trip, without joining a whole graph
# into one character string. Parent list shape and order are retained.
json_plain_chunks <- function(value, leaf_bytes = 8 * 1024^2) {
  plain <- is.list(value) && !is.object(value) && is.null(dim(value))
  if (!plain || !length(value) || as.numeric(utils::object.size(value)) <= leaf_bytes)
    return(jsonlite::fromJSON(jsonlite::toJSON(value, auto_unbox = TRUE,
      null = "null", digits = NA), simplifyVector = FALSE))
  out <- lapply(value, json_plain_chunks, leaf_bytes = leaf_bytes)
  labels <- names(value)
  if (!is.null(labels)) {
    missing <- is.na(labels) | !nzchar(labels)
    labels[missing] <- as.character(which(missing))
    names(out) <- labels
  }
  out
}
