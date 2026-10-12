# T007 interprets source assertions. All side labels are local to a block.
phase_hash <- function(x) digest::digest(charToRaw(canonical_bytes(x)),
  algo = "sha256", serialize = FALSE)

phase_format <- function(source) {
  raw <- source$raw_format
  if (is.null(raw)) return(list())
  keys <- strsplit(raw$keys, ":", fixed = TRUE)[[1L]]
  values <- strsplit(raw$sample, ":", fixed = TRUE)[[1L]]
  if (length(keys) != length(values) || anyDuplicated(keys)) return(list())
  stats::setNames(as.list(values), keys)
}

phase_genotype <- function(token, alt_count) {
  if (is.null(token) || identical(token, ".") || !nzchar(token))
    return(list(state = "MISSING", alleles = NULL))
  if (grepl("[/|]", token) && grepl("[/]", token) && grepl("[|]", token))
    return(list(state = "MALFORMED", alleles = NULL))
  sep <- if (grepl("|", token, fixed = TRUE)) "|" else if (grepl("/", token, fixed = TRUE)) "/" else NULL
  parts <- if (is.null(sep)) token else strsplit(token, sep, fixed = TRUE)[[1L]]
  if (length(parts) != 2L || any(!grepl("^[0-9]+$", parts)))
    return(list(state = "NON_DIPLOID_OR_INCOMPLETE", alleles = NULL))
  alleles <- suppressWarnings(as.numeric(parts))
  if (any(!is.finite(alleles)) || any(alleles > alt_count) || any(alleles > 2147483647))
    return(list(state = "INVALID_INDEX", alleles = NULL))
  list(state = if (identical(sep, "|")) "PHASED" else "UNPHASED",
       alleles = as.integer(alleles))
}

phase_side <- function(alleles, selected) {
  hits <- which(alleles == selected)
  if (!length(hits)) "ABSENT" else if (length(hits) == 2L) "BOTH"
  else if (hits == 1L) "LEFT" else "RIGHT"
}

phase_source_assertions <- function(row, source) {
  if (!identical(source$format, "VCF") || is.null(row$contig)) return(list())
  fmt <- phase_format(source)
  alt_count <- length(source$original$alts)
  selected <- source$selected_alt_index
  if (is.null(selected) || length(selected) != 1L || selected < 1L || selected > alt_count)
    return(list(list(status = "INVALID_INDEX", input_id = source$input_id)))
  gt <- phase_genotype(fmt$GT, alt_count)
  pgt <- phase_genotype(fmt$PGT, alt_count)
  has_pid <- !is.null(fmt$PID) && nzchar(fmt$PID) && fmt$PID != "."
  has_ps <- !is.null(fmt$PS) && nzchar(fmt$PS) && fmt$PS != "."
  record_key <- paste0("record:", phase_hash(list(input_id = source$input_id,
    original = source$original, raw_format = source$raw_format)))
  base <- list(input_id = source$input_id, source_record = record_key,
    caller = source$caller,
    selected_alt_index = as.integer(selected), raw_gt = fmt$GT,
    raw_ps = fmt$PS, raw_pid = fmt$PID, raw_pgt = fmt$PGT)
  invalid_pgt <- has_pid && !is.null(fmt$PGT) && pgt$state != "PHASED"
  incompatible <- has_pid && pgt$state == "PHASED" && !is.null(gt$alleles) &&
    !identical(sort(gt$alleles), sort(pgt$alleles))
  opposite <- has_pid && pgt$state == "PHASED" && gt$state == "PHASED" &&
    !identical(gt$alleles, pgt$alleles)
  if (invalid_pgt || incompatible || opposite)
    return(list(c(base, list(status = "CONFLICT", genotype_state = "UNKNOWN",
      mechanism = NULL, block_id = NULL, local_side = NULL))))
  make <- function(genotype, mechanism, tag) {
    side <- phase_side(genotype$alleles, selected)
    scope <- if (is.null(tag)) record_key else tag
    namespace <- list(input_id = source$input_id, patient_id = row$patient_id,
      sample_id = row$sample_id, vcf_sample = source$sample_id,
      contig = row$contig, mechanism = mechanism)
    block_id <- paste0("block:", phase_hash(c(namespace, list(scope = scope))))
    c(base, list(status = "INTERPRETABLE", genotype_state = side,
      mechanism = mechanism, block_id = block_id, local_side =
        if (side %in% c("LEFT", "RIGHT")) side else NULL,
      namespace = paste0("namespace:", phase_hash(namespace)), block_label = tag))
  }
  out <- list()
  if (gt$state == "PHASED") out[[length(out) + 1L]] <- make(gt, "GT_PS", if (has_ps) fmt$PS else NULL)
  if (has_pid && pgt$state == "PHASED") out[[length(out) + 1L]] <- make(pgt, "PGT_PID", fmt$PID)
  if (length(out)) return(out)
  state <- if (gt$state == "UNPHASED") "UNPHASED" else gt$state
  list(c(base, list(status = state, genotype_state = if (is.null(gt$alleles))
    "UNKNOWN" else phase_side(gt$alleles, selected), mechanism = NULL,
    block_id = NULL, local_side = NULL)))
}

phase_marker_key <- function(row) phase_hash(list(patient_id = row$patient_id,
  sample_id = row$sample_id, assembly = row$genome_assembly,
  contig = row$contig, pos = row$pos, ref = row$ref, alt = row$alt))

phase_assertion_key <- function(x) canonical_bytes(x[c("input_id", "source_record",
  "selected_alt_index", "mechanism", "block_id", "local_side", "status")])

phase_sort <- function(xs, key) {
  if (!length(xs)) return(list())
  xs[order(vapply(xs, key, ""), method = "radix")]
}

phase_conflicted_assertions <- function(assertions) {
  if (any(vapply(assertions, function(a) identical(a$status, "CONFLICT"), TRUE))) return(TRUE)
  usable <- Filter(function(a) identical(a$status, "INTERPRETABLE"), assertions)
  if (length(usable) < 2L) return(FALSE)
  groups <- split(usable, vapply(usable, `[[`, "", "input_id"))
  any(vapply(groups, function(group) {
    states <- unique(vapply(group, `[[`, "", "genotype_state"))
    if (any(states %in% c("BOTH", "ABSENT")) &&
        any(states %in% c("LEFT", "RIGHT"))) return(TRUE)
    by_block <- split(group, vapply(group, `[[`, "", "block_id"))
    any(vapply(by_block, function(block)
      length(unique(vapply(block, `[[`, "", "genotype_state"))) > 1L, TRUE))
  }, TRUE))
}

phase_orientations <- function(rows) {
  anchors <- list(); states <- list(); n <- 0L
  for (row in rows) {
    assertions <- row$event_details_json$phase_assertions
    for (a in assertions) {
      key <- paste(phase_marker_key(row), a$input_id, sep = "/")
      states[[key]] <- c(states[[key]], if (identical(a$status, "CONFLICT"))
        "CONFLICT" else if (identical(a$status, "INTERPRETABLE"))
        a$genotype_state else character())
      if (!identical(a$status, "INTERPRETABLE") || is.null(a$local_side)) next
      n <- n + 1L
      anchors[[n]] <- list(marker = phase_marker_key(row), event_id = row$event_id,
        block_id = a$block_id, input_id = a$input_id, side = a$local_side,
        source_record = a$source_record, mechanism = a$mechanism)
    }
  }
  anchors <- Filter(function(a) {
    observed <- unique(states[[paste(a$marker, a$input_id, sep = "/")]])
    length(observed) == 1L && observed %in% c("LEFT", "RIGHT")
  }, anchors)
  n <- length(anchors)
  if (n < 2L) return(list())
  by_marker <- split(anchors, vapply(anchors, `[[`, "", "marker"))
  evidence <- list(); k <- 0L
  for (group in by_marker) {
    group <- phase_sort(group, function(x) canonical_bytes(x))
    if (length(group) < 2L) next
    for (i in seq_len(length(group) - 1L)) for (j in seq.int(i + 1L, length(group))) {
      a <- group[[i]]; b <- group[[j]]
      # A record carrying both GT/PS and PGT/PID maps those two mechanisms.
      # Otherwise, only independent inputs with an exact shared marker map.
      same_record <- identical(a$input_id, b$input_id) &&
        identical(a$source_record, b$source_record) &&
        !identical(a$mechanism, b$mechanism)
      cross_source <- !identical(a$input_id, b$input_id)
      if ((!same_record && !cross_source) || identical(a$block_id, b$block_id)) next
      pair <- sort(c(a$block_id, b$block_id), method = "radix")
      k <- k + 1L
      evidence[[k]] <- list(source_block = pair[[1L]], target_block = pair[[2L]],
        orientation = if (identical(a$side, b$side)) "SAME" else "FLIPPED",
        evidence = if (same_record) "SOURCE_ASSERTION" else "SHARED_MARKER",
        marker = a$marker, anchor_event_ids = as.list(sort(unique(c(a$event_id,
          b$event_id)), method = "radix")), source_records = as.list(sort(c(
          a$source_record, b$source_record), method = "radix")))
    }
  }
  if (!length(evidence)) return(list())
  by_pair <- split(evidence, vapply(evidence, function(x)
    paste(x$source_block, x$target_block, sep = "/"), ""))
  edges <- lapply(by_pair, function(group) {
    group <- phase_sort(group, function(x) canonical_bytes(x))
    orientations <- unique(vapply(group, `[[`, "", "orientation"))
    list(kind = "block_orientation", source_block = group[[1L]]$source_block,
      target_block = group[[1L]]$target_block,
      orientation = if (length(orientations) == 1L) orientations[[1L]] else "CONFLICT",
      evidence = group[[1L]]$evidence, anchors = group)
  })
  edges <- phase_sort(unname(edges), function(x) paste(x$source_block, x$target_block))
  # A parity contradiction disables cross-source reconciliation throughout
  # its component. Source-local assertions retain their own orientation.
  nodes <- sort(unique(unlist(lapply(edges, function(e) c(e$source_block,
    e$target_block)))), method = "radix")
  parent <- stats::setNames(nodes, nodes)
  parity <- stats::setNames(integer(length(nodes)), nodes)
  find <- function(x) {
    p <- parent[[x]]
    if (identical(p, x)) return(list(root = x, side = 0L))
    f <- find(p)
    list(root = f$root, side = bitwXor(parity[[x]], f$side))
  }
  bad <- character()
  for (e in edges) {
    a <- find(e$source_block); b <- find(e$target_block)
    if (e$orientation == "CONFLICT") { bad <- c(bad, a$root, b$root); next }
    want <- as.integer(e$orientation == "FLIPPED")
    if (identical(a$root, b$root)) {
      if (bitwXor(a$side, b$side) != want) bad <- c(bad, a$root)
    } else {
      parent[[b$root]] <- a$root
      parity[[b$root]] <- bitwXor(bitwXor(a$side, b$side), want)
    }
  }
  # Include directly contradictory pairs when identifying poisoned nodes.
  bad_nodes <- unique(unlist(lapply(edges, function(e) if (e$orientation == "CONFLICT")
    c(e$source_block, e$target_block) else character())))
  roots <- unique(c(vapply(bad, function(x) find(x)$root, ""),
    vapply(bad_nodes, function(x) find(x)$root, "")))
  if (length(roots)) for (i in seq_along(edges))
    if (identical(edges[[i]]$evidence, "SHARED_MARKER") &&
        find(edges[[i]]$source_block)$root %in% roots)
      edges[[i]]$orientation <- "CONFLICT"
  edges
}

phase_resolve_events <- function(rows) {
  if (!length(rows)) return(rows)
  before <- lapply(rows, function(x) x[c("event_id", "event_key", "canonical_event",
    "event_role", "contig", "pos", "ref", "alt")])
  for (i in seq_along(rows)) {
    row <- rows[[i]]
    assertions <- unlist(lapply(row$event_details_json$source_records,
      function(s) phase_source_assertions(row, s)), recursive = FALSE)
    assertions <- phase_sort(assertions, phase_assertion_key)
    row$event_details_json$phase_assertions <- assertions
    usable <- Filter(function(a) identical(a$status, "INTERPRETABLE") &&
      !identical(a$genotype_state, "ABSENT"), assertions)
    conflict <- phase_conflicted_assertions(assertions)
    row["phase_status"] <- list(if (row$event_type == "SEQUENCE") "NOT_APPLICABLE"
      else if (conflict) "NOT_EVALUABLE"
      else if (length(usable)) "EVALUATED" else "NOT_PROVIDED")
    blocks <- sort(unique(vapply(usable, `[[`, "", "block_id")), method = "radix")
    namespaces <- sort(unique(vapply(usable, `[[`, "", "namespace")), method = "radix")
    row["phase_set"] <- list(if (length(blocks) == 1L) blocks else NULL)
    row["phase_namespace"] <- list(if (length(namespaces) == 1L) namespaces else NULL)
    if (conflict) row$event_warning_codes <- as.list(sort(unique(c(
      unlist(row$event_warning_codes), "PHASE_UNCERTAIN")), method = "radix"))
    rows[[i]] <- row
  }
  edges <- phase_orientations(rows)
  # Each graph edge is stored once, on its least anchor event ID. The block
  # IDs in the edge target the per-event memberships; no event-pair expansion.
  for (edge in edges) {
    ids <- sort(unique(unlist(lapply(edge$anchors, `[[`, "anchor_event_ids"))), method = "radix")
    owner <- match(ids[[1L]], vapply(rows, `[[`, "", "event_id"))
    if (is.na(owner)) stop("Phase edge anchor missing")
    rows[[owner]]$phase_links_json[[length(rows[[owner]]$phase_links_json) + 1L]] <- edge
  }
  for (i in seq_along(rows)) {
    rows[[i]]$phase_links_json <- phase_sort(rows[[i]]$phase_links_json,
      function(x) paste(x$source_block, x$target_block))
    if (any(vapply(rows[[i]]$phase_links_json,
      function(x) identical(x$orientation, "CONFLICT"), TRUE))) {
      rows[[i]]$event_warning_codes <- as.list(sort(unique(c(
        unlist(rows[[i]]$event_warning_codes), "PHASE_UNCERTAIN")), method = "radix"))
    }
  }
  after <- lapply(rows, function(x) x[c("event_id", "event_key", "canonical_event",
    "event_role", "contig", "pos", "ref", "alt")])
  if (!identical(before, after)) stop("Phase resolution changed event identity")
  validate_table(rows, "events")
  rows
}

# Internal T008 handoff: derive a relation without materializing all pairs.
phase_relative_events <- function(rows, first_id, second_id) {
  by_id <- stats::setNames(rows, vapply(rows, `[[`, "", "event_id"))
  a <- by_id[[first_id]]; b <- by_id[[second_id]]
  if (is.null(a) || is.null(b) || a$patient_id != b$patient_id ||
      a$sample_id != b$sample_id || is.null(a$contig) ||
      !identical(a$contig, b$contig) ||
      identical(a$phase_status, "NOT_EVALUABLE") ||
      identical(b$phase_status, "NOT_EVALUABLE"))
    return(list(relation = "UNKNOWN", evidence = "NONE"))
  candidates <- list()
  links <- unlist(lapply(rows, `[[`, "phase_links_json"), recursive = FALSE)
  graph <- Filter(function(e) identical(e$kind, "block_orientation") &&
    e$orientation != "CONFLICT", links)
  walk <- function(from, to) {
    if (from == to) return(list(side = 0L, evidence = "SOURCE_ASSERTION"))
    queue <- list(list(node = from, side = 0L, depth = 0L, evidence = "SOURCE_ASSERTION")); seen <- character()
    while (length(queue)) {
      item <- queue[[1L]]; queue <- queue[-1L]
      if (item$node %in% seen) next
      seen <- c(seen, item$node)
      if (item$node == to) return(list(side = item$side,
        evidence = if (item$depth > 1L) "TRANSITIVE" else item$evidence))
      for (e in graph) if (item$node %in% c(e$source_block, e$target_block))
        queue[[length(queue) + 1L]] <- list(node = if (item$node == e$source_block)
          e$target_block else e$source_block,
          side = bitwXor(item$side, as.integer(e$orientation == "FLIPPED")),
          depth = item$depth + 1L, evidence = e$evidence)
    }
    NULL
  }
  for (x in a$event_details_json$phase_assertions) for (y in b$event_details_json$phase_assertions) {
    if (is.null(x$local_side) || is.null(y$local_side)) next
    path <- walk(x$block_id, y$block_id)
    if (is.null(path)) next
    relation <- if (bitwXor(bitwXor(as.integer(x$local_side == "RIGHT"),
      as.integer(y$local_side == "RIGHT")), path$side) == 0L)
      "SAME_HAPLOTYPE" else "OPPOSITE_HAPLOTYPE"
    candidates[[length(candidates) + 1L]] <- list(relation = relation,
      evidence = path$evidence)
  }
  if (!length(candidates) || length(unique(vapply(candidates, `[[`, "", "relation"))) != 1L)
    return(list(relation = "UNKNOWN", evidence = "NONE"))
  candidates[[order(vapply(candidates, function(x) match(x$evidence,
    c("SOURCE_ASSERTION", "SHARED_MARKER", "TRANSITIVE")), 0L))[[1L]]]]
}
