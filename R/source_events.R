# T025 source-event overlay. Provider records are evidence; biological identity
# contains only the explicitly mapped structure and bound reference identity.
source_extension_version <- function() "source-event-2.0.0"
source_identity_schema_version <- function() "source-biological-id-v1"
source_id <- function(kind, payload) paste0(kind, ":", infra_hash(list(
  version = source_identity_schema_version(), kind = kind, payload = payload)))

source_table <- function(path, required, label) {
  lines <- reference_lines(path)
  if (!length(lines)) stop(label, " table is empty", call. = FALSE)
  header <- stringi::stri_split_fixed(sub("^#", "", lines[[1L]]), "\t",
    omit_empty = FALSE)[[1L]]
  if (anyDuplicated(header) || any(!nzchar(header)) ||
      !all(required %in% header)) stop("Unsupported ", label, " header", call. = FALSE)
  rows <- lapply(seq_along(lines)[-1L], function(i) {
    fields <- stringi::stri_split_fixed(lines[[i]], "\t", omit_empty = FALSE)[[1L]]
    if (length(fields) != length(header)) stop("Malformed ", label, " record ", i)
    stats::setNames(as.list(fields), header)
  })
  rows
}

source_breakpoint <- function(x, strand = NULL) {
  m <- regexec("^([^:]+):([1-9][0-9]*)(?::([+-]))?$", x, perl = TRUE)
  parts <- regmatches(x, m)[[1L]]
  if (!length(parts)) stop("Unsupported fusion breakpoint: ", x)
  found <- if (length(parts) >= 4L && nzchar(parts[[4L]])) parts[[4L]] else strand
  if (is.null(found) || !found %in% c("+", "-"))
    stop("Fusion breakpoint strand is unresolved")
  position <- suppressWarnings(as.integer(parts[[3L]]))
  if (is.na(position) || position < 1L) stop("Invalid fusion breakpoint coordinate")
  list(contig = parts[[2L]], position = position, strand = found)
}

source_gene <- function(x) {
  parts <- strsplit(x, "^", fixed = TRUE)[[1L]]
  list(name = parts[[1L]], id = if (length(parts) > 1L) parts[[2L]] else NULL)
}

source_tumor_specificity <- function(value) {
  if (is.null(value) || !nzchar(value)) return("normal_not_evaluated")
  allowed <- c("tumor_specific_supported", "observed_in_tumor_only",
    "normal_not_evaluated", "normal_observed", "unresolved")
  if (!value %in% allowed) stop("Unsupported source tumor-specificity status")
  value
}

source_fusion_rows <- function(input, assembly) {
  star <- identical(input$input_role, "star_fusion")
  required <- if (star) c("FusionName", "JunctionReadCount", "SpanningFragCount",
    "LeftGene", "LeftBreakpoint", "RightGene", "RightBreakpoint") else
    c("gene1", "gene2", "strand1(gene/fusion)", "strand2(gene/fusion)",
      "breakpoint1", "breakpoint2", "split_reads1", "split_reads2",
      "discordant_mates", "confidence")
  rows <- source_table(input$resolved_path, required,
    if (star) "STAR-Fusion" else "Arriba")
  checksum <- digest::digest(file = input$resolved_path, algo = "sha256", serialize = FALSE)
  lapply(seq_along(rows), function(i) {
    r <- rows[[i]]
    if (star) {
      left <- source_breakpoint(r$LeftBreakpoint)
      right <- source_breakpoint(r$RightBreakpoint)
      genes <- list(source_gene(r$LeftGene), source_gene(r$RightGene))
      support <- list(junction_reads = as.integer(r$JunctionReadCount),
        spanning_fragments = as.integer(r$SpanningFragCount), ffpm = r$FFPM)
    } else {
      strand <- lapply(c(r[["strand1(gene/fusion)"]], r[["strand2(gene/fusion)"]]),
        function(s) {
          p <- strsplit(s, "/", fixed = TRUE)[[1L]]
          if (length(p) != 2L || !all(p %in% c("+", "-")))
            stop("Unsupported Arriba strand/orientation")
          p[[2L]]
        })
      left <- source_breakpoint(r$breakpoint1, strand[[1L]])
      right <- source_breakpoint(r$breakpoint2, strand[[2L]])
      genes <- list(list(name = r$gene1, id = r$gene_id1),
        list(name = r$gene2, id = r$gene_id2))
      support <- list(junction_reads = as.integer(r$split_reads1) +
        as.integer(r$split_reads2), split_reads1 = as.integer(r$split_reads1),
        split_reads2 = as.integer(r$split_reads2),
        spanning_fragments = as.integer(r$discordant_mates), confidence = r$confidence)
    }
    if (anyNA(unlist(support[c("junction_reads", "spanning_fragments")])) ||
        any(unlist(support[c("junction_reads", "spanning_fragments")]) < 0L))
      stop("Invalid fusion read support")
    structure <- list(sample_id = input$sample_id, reference_build = assembly,
      source_type = "fusion", left = left, right = right, genes = genes)
    gene_key <- function(g) if (!is.null(g$id) && nzchar(g$id) && g$id != ".")
      g$id else g$name
    event_id <- source_id("fusion_event_id", list(sample_id = input$sample_id,
      reference_build = assembly, left = left, right = right,
      gene_5prime = gene_key(genes[[1L]]),
      gene_3prime = gene_key(genes[[2L]])))
    list(source_type = "fusion", source_event_id = event_id, structure = structure,
      observation_id = source_id("source_observation_id", list(input = input$input_id,
        checksum = checksum, line = i + 1L)), provider = if (star) "STAR-Fusion" else "Arriba",
      provider_version = input$version, source_record_id = if (star) r$FusionName else
        paste0(input$input_id, ":", i + 1L), support = support,
      provider_frame = if (star) r$PROT_FUSION_TYPE else r$reading_frame,
      provider_gene_fusion_strands = if (star) NULL else list(
        r[["strand1(gene/fusion)"]], r[["strand2(gene/fusion)"]]),
      provider_classification = if (star) r$SpliceType else r$type,
      provider_filters = if (star) r$annots else r$filters,
      raw_fields = r, artifact = list(input_id = input$input_id, sha256 = checksum),
      mapping_status = "PENDING", transcript_hypothesis_status = "PENDING",
      tumor_specificity = source_tumor_specificity(r$tumor_specificity))
  })
}

source_splice_rows <- function(input, assembly) {
  required <- c("source_record_id", "source_type", "gene_id", "contig",
    "strand", "donor", "acceptor", "transcript_id", "observed_status")
  rows <- source_table(input$resolved_path, required, "splice structure")
  checksum <- digest::digest(file = input$resolved_path, algo = "sha256", serialize = FALSE)
  supported <- c("exon_skip", "alternative_donor", "alternative_acceptor",
    "retained_intron", "novel_junction")
  lapply(seq_along(rows), function(i) {
    r <- rows[[i]]
    if (!r$source_type %in% supported || !r$strand %in% c("+", "-") ||
        !r$observed_status %in% c("observed", "specified"))
      stop("Unsupported splice structure record")
    donor <- suppressWarnings(as.integer(r$donor))
    acceptor <- suppressWarnings(as.integer(r$acceptor))
    if (anyNA(c(donor, acceptor)) || any(c(donor, acceptor) < 1L))
      stop("Invalid splice junction coordinates")
    structure <- list(sample_id = input$sample_id, reference_build = assembly,
      source_type = r$source_type, gene_id = r$gene_id, contig = r$contig,
      strand = r$strand, donor = donor, acceptor = acceptor,
      transcript_id = if (nzchar(r$transcript_id)) r$transcript_id else NULL,
      retained_start = if (!is.null(r$retained_start) && nzchar(r$retained_start))
        as.integer(r$retained_start) else NULL,
      retained_end = if (!is.null(r$retained_end) && nzchar(r$retained_end))
        as.integer(r$retained_end) else NULL)
    if (anyNA(unlist(structure[c("retained_start", "retained_end")])) ||
        any(unlist(structure[c("retained_start", "retained_end")]) < 1L))
      stop("Invalid retained interval")
    list(source_type = r$source_type,
      source_event_id = source_id("source_event_id", structure), structure = structure,
      observation_id = source_id("source_observation_id", list(input = input$input_id,
        checksum = checksum, line = i + 1L)), provider = "explicit_splice_structure",
      provider_version = input$version, source_record_id = r$source_record_id,
      support = list(junction_reads = r$junction_reads, unique_junction_reads = r$unique_junction_reads,
        junction_fraction = r$junction_fraction,
        normal_junction_reads = r$normal_junction_reads,
        transcript_tpm = r$transcript_tpm), raw_fields = r,
      artifact = list(input_id = input$input_id, sha256 = checksum),
      observed_status = r$observed_status,
      tumor_specificity = source_tumor_specificity(r$tumor_specificity))
  })
}

source_spliceai_rows <- function(input) {
  lines <- reference_lines(input$resolved_path)
  checksum <- digest::digest(file = input$resolved_path, algo = "sha256", serialize = FALSE)
  out <- list()
  for (i in seq_along(lines)) {
    line <- lines[[i]]
    if (startsWith(line, "#")) next
    fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    if (length(fields) < 8L) stop("Malformed SpliceAI VCF record")
    info <- strsplit(fields[[8L]], ";", fixed = TRUE)[[1L]]
    annotation <- sub("^SpliceAI=", "", info[startsWith(info, "SpliceAI=")])
    annotation <- unlist(strsplit(annotation, ",", fixed = TRUE), use.names = FALSE)
    for (a in annotation) {
      values <- strsplit(a, "|", fixed = TRUE)[[1L]]
      if (length(values) != 10L) stop("Malformed SpliceAI annotation")
      names(values) <- c("ALLELE", "SYMBOL", "DS_AG", "DS_AL", "DS_DG",
        "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL")
      out[[length(out) + 1L]] <- list(sample_id = input$sample_id,
        prediction_id = source_id("splice_prediction_id", list(input = input$input_id,
          line = i, values = as.list(values))), values = as.list(values),
        locus = fields[1:5], artifact = list(input_id = input$input_id, sha256 = checksum),
        sequence_generating = FALSE)
    }
  }
  out
}

source_star_junction_rows <- function(input) {
  lines <- reference_lines(input$resolved_path)
  checksum <- digest::digest(file = input$resolved_path, algo = "sha256", serialize = FALSE)
  lapply(seq_along(lines), function(i) {
    fields <- stringi::stri_split_fixed(lines[[i]], "\t", omit_empty = FALSE)[[1L]]
    if (length(fields) != 9L || !all(grepl("^[0-9]+$", fields[2:9])))
      stop("Unsupported STAR SJ.out.tab record")
    values <- as.integer(fields[2:9])
    if (anyNA(values) || values[[1L]] > values[[2L]] ||
        !values[[3L]] %in% 0:2 || !values[[4L]] %in% 0:6 ||
        !values[[5L]] %in% 0:1)
      stop("Invalid STAR SJ.out.tab record")
    list(junction_evidence_id = source_id("junction_evidence_id", list(
      input = input$input_id, sha256 = checksum, line = i)),
      sample_id = input$sample_id, contig = fields[[1L]],
      intron_first_base = values[[1L]], intron_last_base = values[[2L]],
      strand = c("unknown", "+", "-")[[values[[3L]] + 1L]],
      motif_code = values[[4L]], annotation_code = values[[5L]],
      unique_junction_reads = values[[6L]], multi_mapping_reads = values[[7L]],
      maximum_overhang = values[[8L]], provider = "STAR",
      provider_version = input$version,
      artifact = list(input_id = input$input_id, sha256 = checksum),
      sequence_generating = FALSE)
  })
}
