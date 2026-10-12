# Transient indices only: keys and values never participate in biological IDs.
external_value_index <- function(keys, values) {
  if (!length(keys)) return(new.env(hash = TRUE, parent = emptyenv()))
  list2env(stats::setNames(values, keys), hash = TRUE, parent = emptyenv())
}

external_group_index <- function(keys) {
  first <- unique(keys)
  # Integer grouping preserves original within-key order and avoids repeated
  # appends/copies for a large fan-out bucket. Environment lookup is hashed.
  groups <- split(seq_along(keys), match(keys, first))
  external_value_index(first, unname(groups))
}

external_index_get <- function(index, key) get0(key, envir = index,
  inherits = FALSE, ifnotfound = integer())
