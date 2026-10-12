schema_contract <- function(path = system.file("schema", "contract.json", package = "Neoantimon2")) {
  if (!nzchar(path) || !file.exists(path)) stop("Contract not found")
  jsonlite::read_json(path, simplifyVector = FALSE)
}

is_missing <- function(x) is.null(x) || (is.atomic(x) && length(x) == 1L && is.na(x) && !is.nan(x))

validate_value <- function(x, f, contract, path) {
  fail <- function(msg) stop(paste0(path, ": ", msg), call. = FALSE)
  if (!is.null(f$vocabulary)) {
    if (!is.null(f$enum)) fail("enum and vocabulary cannot coexist")
    if (is.null(contract$vocabularies[[f$vocabulary]])) fail("Unknown vocabulary reference")
    f$enum <- contract$vocabularies[[f$vocabulary]]
  }
  if (!is.null(f$items_vocabulary)) {
    if (!is.null(f$items_enum)) fail("items_enum and items_vocabulary cannot coexist")
    if (is.null(contract$vocabularies[[f$items_vocabulary]])) fail("Unknown items vocabulary reference")
  }
  if (is_missing(x)) {
    if (!isTRUE(f$nullable)) fail("null is not allowed")
    return(invisible(TRUE))
  }
  scalar <- is.atomic(x) && length(x) == 1L
  valid <- switch(f$type,
    string = scalar && is.character(x),
    boolean = scalar && is.logical(x),
    integer = scalar && is.numeric(x) && is.finite(x) && x == floor(x) && abs(x) <= 2147483647,
    number = scalar && is.numeric(x) && is.finite(x),
    scalar = scalar && (is.character(x) || is.logical(x) || (is.numeric(x) && is.finite(x))),
    array = is.list(x) && is.null(names(x)),
    object = is.list(x) && !is.null(names(x)),
    json = is.list(x) || scalar,
    FALSE)
  if (!isTRUE(valid)) fail(paste("expected", f$type))
  if (!is.null(f$enum) && !x %in% unlist(f$enum)) fail("value outside enum")
  if (!is.null(f$minimum) && x < f$minimum) fail("below minimum")
  if (!is.null(f$maximum) && x > f$maximum) fail("above maximum")
  if (!is.null(f$pattern) && !grepl(f$pattern, x, perl = TRUE)) fail("invalid sequence/pattern")
  if (!is.null(f$record)) validate_record(x, contract$records[[f$record]], contract, path)
  if (f$type == "array") {
    if (isTRUE(f$unique_items) && anyDuplicated(vapply(x, canonical_bytes, ""))) fail("duplicate array item")
    if (isTRUE(f$sorted_set) && !identical(unlist(x), sort(unlist(x), method = "radix"))) fail("set must be sorted")
    for (i in seq_along(x)) {
      p <- paste0(path, "[", i, "]")
      if (!is.null(f$items_record)) validate_record(x[[i]], contract$records[[f$items_record]], contract, p)
      else validate_value(x[[i]], list(type = f$items_type, nullable = FALSE,
        enum = f$items_enum, vocabulary = f$items_vocabulary), contract, p)
    }
  }
  invisible(TRUE)
}

validate_record <- function(row, def, contract, path = "record") {
  if (!is.list(row) || (length(row) && is.null(names(row))) || anyDuplicated(names(row)))
    stop(path, ": expected unique named record")
  allowed <- vapply(def$fields, `[[`, "", "name")
  if (length(setdiff(names(row), allowed))) stop(path, ": unknown fields: ", paste(setdiff(names(row), allowed), collapse = ", "))
  for (f in def$fields) {
    if (!f$name %in% names(row)) {
      if (isTRUE(f$required)) stop(path, ": required field ", f$name)
    } else validate_value(row[[f$name]], f, contract, paste0(path, ".", f$name))
  }
  invisible(TRUE)
}

validate_table <- function(rows, name, contract = schema_contract()) {
  def <- contract$tables[[name]]
  if (is.null(def)) stop("Unknown table")
  if (is.data.frame(rows)) {
    required <- vapply(Filter(function(f) isTRUE(f$required), def$fields), `[[`, "", "name")
    if (length(setdiff(required, names(rows)))) stop("Missing required columns")
    if (length(setdiff(names(rows), vapply(def$fields, `[[`, "", "name")))) stop("Unknown columns")
    rows <- lapply(seq_len(nrow(rows)), function(i) lapply(rows, function(col) col[[i]]))
  }
  if (!is.list(rows) || !is.null(names(rows))) stop("Table must be array of records")
  keys <- character()
  for (i in seq_along(rows)) {
    row <- rows[[i]]
    validate_record(row, def, contract, paste0(name, "[", i, "]"))
    key <- canonical_bytes(unname(row[unlist(def$primary_key)]))
    if (key %in% keys) stop(name, ": duplicate primary key")
    keys <- c(keys, key)
  }
  invisible(TRUE)
}

# Field documentation is rendered from the authority, never maintained separately.
schema_field_markdown <- function(contract = schema_contract()) {
  out <- c("<!-- Generated from inst/schema/contract.json; do not edit. -->", "")
  for (name in names(c(contract$tables, contract$records))) {
    def <- c(contract$tables, contract$records)[[name]]
    out <- c(out, paste0("## ", name), "", "| Field | Type | Required | Nullable | Unit | Vocabulary | Definition |",
             "| --- | --- | --- | --- | --- | --- | --- |")
    for (f in def$fields) out <- c(out, paste0("| ", paste(c(f$name, f$type, f$required,
      f$nullable, if (is.null(f$unit)) "\u2014" else f$unit,
      paste(c(f$vocabulary, f$items_vocabulary, unlist(f$enum), unlist(f$items_enum)), collapse = ", "),
      gsub("|", "\\|", f$definition, fixed = TRUE)), collapse = " | "), " |"))
    out <- c(out, "")
  }
  out
}
