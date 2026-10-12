# T028 overlay: structure is biological; provider models and counts are evidence.
long_read_version <- function() "long-read-source-1.0.0"
long_read_identity_version <- function() "long-read-isoform-id-v1"
long_read_id <- function(kind, payload) paste0(kind, ":", infra_hash(list(
  version = long_read_identity_version(), kind = kind, payload = payload)))
long_read_roles <- function() c("long_read_gtf", "sqanti3_gtf", "flair_gtf",
  "isoquant_gtf", "long_read_metadata", "sqanti3_metadata", "flair_metadata",
  "isoquant_metadata", "long_read_fasta", "long_read_settings")
long_read_inputs <- function(project) Filter(function(x)
  x$input_role %in% long_read_roles(), project$inputs)
long_read_artifacts <- function(project) lapply(long_read_inputs(project), function(x) {
  if (is.null(x$version)) stop("Long-read provider/parser version declaration required")
  observed <- integrity_record(x$resolved_path, project$config$checksum, x$integrity$declared)
  if (!identical(observed, x$integrity)) stop("Long-read input changed since project binding")
  list(input_id = x$input_id, input_role = x$input_role, sample_id = x$sample_id,
    version = x$version, sha256 = digest::digest(file = x$resolved_path,
      algo = "sha256", serialize = FALSE))
})

long_read_attributes <- function(text) {
  fields <- strsplit(text, ";", fixed = TRUE)[[1L]]
  out <- list()
  for (field in trimws(fields)) {
    if (!nzchar(field)) next
    if (grepl("=", field, fixed = TRUE)) {
      k <- regexpr("=", field, fixed = TRUE)[[1L]]
      key <- substring(field, 1L, k - 1L); value <- utils::URLdecode(substring(field, k + 1L))
    } else {
      match <- stringi::stri_match_first_regex(field, '^([^ ]+) +"([^"]*)"$')
      if (is.na(match[1L, 2L])) stop("Malformed long-read GTF attribute")
      key <- match[1L, 2L]; value <- match[1L, 3L]
    }
    if (!is.null(out[[key]])) stop("Duplicate long-read annotation attribute")
    out[[key]] <- value
  }
  out
}

long_read_structure <- function(build, contig, strand, exons) {
  if (!strand %in% c("+", "-") || !nzchar(contig) || !nrow(exons) ||
      anyNA(exons) || any(exons$start < 1L | exons$end < exons$start))
    stop("Invalid long-read exon structure")
  exons <- unique(exons[, c("start", "end"), drop = FALSE])
  exons <- exons[order(exons$start, exons$end), , drop = FALSE]
  if (any(exons$start[-1L] <= utils::head(exons$end, -1L)))
    stop("Overlapping long-read exon blocks")
  if (strand == "-") exons <- exons[rev(seq_len(nrow(exons))), , drop = FALSE]
  blocks <- lapply(seq_len(nrow(exons)), function(i) list(
    start = as.integer(exons$start[[i]]), end = as.integer(exons$end[[i]])))
  junctions <- lapply(seq_len(max(0L, length(blocks) - 1L)), function(i) list(
    donor = if (strand == "+") blocks[[i]]$end else blocks[[i]]$start,
    acceptor = if (strand == "+") blocks[[i + 1L]]$start else blocks[[i + 1L]]$end))
  list(reference_build = build, contig = contig, strand = strand,
    exons = blocks, junctions = junctions,
    transcript_start = if (strand == "+") blocks[[1L]]$start else blocks[[1L]]$end,
    transcript_end = if (strand == "+") utils::tail(blocks, 1L)[[1L]]$end else
      utils::tail(blocks, 1L)[[1L]]$start)
}
long_read_structure_id <- function(s) {
  payload <- s[c("reference_build", "contig", "strand", "exons", "junctions",
    "transcript_start", "transcript_end")]
  payload$exons <- lapply(s$exons, function(e) e[c("start", "end")])
  payload$junctions <- lapply(s$junctions, function(j) j[c("donor", "acceptor")])
  long_read_id("isoform_structure_id", payload)
}
long_read_integer_evidence <- function(x) {
  if (length(x) != 1L || is.na(x) || !grepl("^[0-9]+$", x)) return(NA_integer_)
  suppressWarnings(as.integer(x))
}

# One parser for GTF and GFF3; explicit Parent relationships, never row positions.
long_read_models <- function(path, build) {
  lines <- reference_lines(path)
  features <- list(); transcripts <- list()
  for (line in lines[!startsWith(lines, "#") & nzchar(lines)]) {
    f <- stringi::stri_split_fixed(line, "\t", omit_empty = FALSE)[[1L]]
    if (length(f) != 9L) stop("Malformed long-read GTF/GFF record")
    if (!f[[3L]] %in% c("transcript", "mRNA", "exon", "CDS", "start_codon", "stop_codon")) next
    if (!all(grepl("^[1-9][0-9]*$", f[4:5])) || !f[[7L]] %in% c("+", "-"))
      stop("Invalid long-read GTF/GFF coordinates or strand")
    a <- long_read_attributes(f[[9L]])
    ids <- if (!is.null(a$transcript_id)) a$transcript_id else
      if (f[[3L]] %in% c("transcript", "mRNA")) a$ID else a$Parent
    if (is.null(ids) || !nzchar(ids)) stop("Missing transcript/Parent structural input")
    ids <- strsplit(ids, ",", fixed = TRUE)[[1L]]
    for (id in ids) {
      record <- list(contig = f[[1L]], feature = f[[3L]], start = as.integer(f[[4L]]),
        end = as.integer(f[[5L]]), strand = f[[7L]], phase = f[[8L]], attributes = a)
      if (anyNA(c(record$start, record$end)) || record$end < record$start)
        stop("Invalid long-read feature interval")
      features[[id]] <- c(features[[id]], list(record))
      if (f[[3L]] %in% c("transcript", "mRNA")) {
        if (!is.null(transcripts[[id]])) stop("Duplicate long-read transcript declaration")
        transcripts[[id]] <- record
      }
    }
  }
  if (!length(features)) stop("No long-read transcript structures")
  lapply(sort(names(features)), function(id) {
    rows <- features[[id]]
    contigs <- unique(vapply(rows, `[[`, "", "contig"))
    strands <- unique(vapply(rows, `[[`, "", "strand"))
    exon <- Filter(function(x) x$feature == "exon", rows)
    if (!length(exon)) stop("Missing required long-read exon structure: ", id)
    genes <- unique(unlist(lapply(rows, function(x) x$attributes$gene_id)))
    tx <- transcripts[[id]]
    if (!is.null(tx) && !length(genes)) genes <- tx$attributes$Parent
    if (length(contigs) != 1L || length(strands) != 1L || length(genes) > 1L)
      return(list(provider_transcript_id = id, provider_gene_ids = as.list(genes),
        status = "CHIMERIC_STRUCTURE_DEFERRED", raw_features = rows, structure = NULL))
    e <- data.frame(start = vapply(exon, `[[`, 0L, "start"), end = vapply(exon, `[[`, 0L, "end"))
    s <- long_read_structure(build, contigs[[1L]], strands[[1L]], e)
    if (!is.null(tx) && (tx$start != min(e$start) || tx$end != max(e$end)))
      stop("Transcript bounds do not equal outer exon boundaries")
    list(provider_transcript_id = id, provider_gene_ids = as.list(genes),
      status = "STRUCTURE_PARSED", structure = s, raw_features = rows,
      attributes = if (is.null(tx)) exon[[1L]]$attributes else tx$attributes,
      cds_features = Filter(function(x) x$feature == "CDS", rows))
  })
}

long_read_metadata <- function(input) {
  sq <- input$input_role == "sqanti3_metadata"
  rows <- source_table(input$resolved_path, if (sq) c("isoform", "chrom", "strand",
    "length", "exons", "structural_category") else character(), "long-read metadata")
  lapply(rows, function(r) {
    keys <- intersect(c("isoform", "transcript_id", "feature_id", "ids"), names(r))
    if (length(keys) != 1L || !nzchar(r[[keys[[1L]]]]))
      stop("Metadata requires exactly one named transcript identifier column")
    list(provider_transcript_id = r[[keys[[1L]]]], raw_fields = r,
      input_id = input$input_id, provider_version = input$version,
      sha256 = digest::digest(file = input$resolved_path, algo = "sha256", serialize = FALSE))
  })
}

# Internal settings are an explicit versioned input artifact, not a public API.
long_read_bindings <- function(project) {
  inputs <- long_read_inputs(project)
  models <- Filter(function(x) endsWith(x$input_role, "_gtf"), inputs)
  if (length(inputs) && !length(models)) stop("Missing required long-read structural GTF/GFF input")
  settings <- unlist(lapply(Filter(function(x) x$input_role == "long_read_settings", inputs),
    function(x) {
      value <- jsonlite::read_json(x$resolved_path)
      reject_duplicate_keys(value)
      if (!identical(value$version, "long-read-run-settings-v1")) stop("Unsupported long-read settings")
      value$bindings
    }), recursive = FALSE)
  ids <- vapply(settings, `[[`, "", "model_input_id")
  if (anyDuplicated(ids) || any(!ids %in% vapply(models, `[[`, "", "input_id")))
    stop("Invalid/duplicate long-read model binding")
  get_input <- function(id, sample, roles) {
    found <- Filter(function(x) x$input_id == id && x$sample_id == sample &&
      x$input_role %in% roles, inputs)
    if (length(found) != 1L) stop("Unresolved long-read auxiliary input binding")
    found[[1L]]
  }
  lapply(models, function(model) {
    provider <- sub("_gtf$", "", model$input_role)
    hit <- settings[ids == model$input_id]
    s <- if (length(hit)) hit[[1L]] else list()
    strict_keys(s, c("model_input_id", "metadata_input_ids", "fasta_input_id",
      "de_novo_aug", "sequencing_technology", "observation_role", "matched_normal_status",
      "region", "timepoint", "cell_context", "reference_build"), "long-read binding")
    if (!is.null(s$de_novo_aug) && (!is.logical(s$de_novo_aug) ||
        length(s$de_novo_aug) != 1L || is.na(s$de_novo_aug))) stop("Invalid de novo AUG setting")
    if (!is.null(s$reference_build) && s$reference_build != project$config$scientific$genome_assembly)
      stop("Long-read reference build conflict")
    sample <- Filter(function(x) x$sample_id == model$sample_id, project$samples)[[1L]]
    role <- if (is.null(s$observation_role)) if (identical(sample$role, "normal")) "normal" else "tumor" else s$observation_role
    if (!role %in% c("tumor", "normal")) stop("Invalid long-read observation role")
    state <- if (is.null(s$matched_normal_status)) "normal_not_profiled" else s$matched_normal_status
    if (!state %in% c("tumor_observed", "normal_observed", "tumor_only_observed",
        "normal_not_profiled", "normal_not_detected", "unresolved")) stop("Invalid matched-normal status")
    metadata <- if (!is.null(s$metadata_input_ids)) lapply(s$metadata_input_ids, get_input,
      sample = model$sample_id, roles = paste0(provider, "_metadata")) else
      Filter(function(x) x$sample_id == model$sample_id &&
        x$input_role == paste0(provider, "_metadata"), inputs)
    fasta <- if (!is.null(s$fasta_input_id)) list(get_input(s$fasta_input_id,
      model$sample_id, "long_read_fasta")) else Filter(function(x)
        x$sample_id == model$sample_id && x$input_role == "long_read_fasta", inputs)
    if (length(fasta) > 1L) stop("Ambiguous FASTA binding; declare fasta_input_id")
    list(input = model, provider = provider, metadata = metadata, fasta = fasta,
      de_novo_aug = isTRUE(s$de_novo_aug), sequencing_technology = s$sequencing_technology,
      matched_normal_status = state, context = list(sample_id = model$sample_id,
        observation_role = role, region = s$region, timepoint = if (is.null(s$timepoint)) sample$timepoint else s$timepoint,
        cell_context = s$cell_context))
  })
}

long_read_fasta <- function(path) {
  out <- list(); id <- NULL; sequence <- character()
  flush <- function() {
    if (is.null(id)) return()
    if (!is.null(out[[id]])) stop("Duplicate provider FASTA identifier")
    seq <- toupper(paste0(sequence, collapse = ""))
    if (!grepl("^[ACGTN]+$", seq)) stop("Invalid provider transcript FASTA")
    out[[id]] <<- seq
  }
  for (line in reference_lines(path)) {
    if (startsWith(line, ">")) {
      flush(); id <- strsplit(substring(line, 2L), "[[:space:]]")[[1L]][[1L]]
      sequence <- character()
    } else {
      if (is.null(id)) stop("Provider FASTA sequence before header")
      sequence <- c(sequence, line)
    }
  }
  flush(); out
}
