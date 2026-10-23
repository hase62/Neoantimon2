event_input_roles <- function() c("somatic_vcf", "germline_vcf", "somatic_maf", "germline_maf", "canonical_events",
  "somatic_annovar", "germline_annovar", "somatic_vep", "germline_vep",
  "sequence_fragment", "sequence_entire_region", "germline_snp_table", "somatic_bnd_table")

event_source_record <- function(input, format, line, original, selected_alt_index = 1L,
                                caller = NULL, raw_info = NULL, raw_format = NULL,
                                projected_info = list(), projected_format = list()) {
  list(input_id = input$input_id, input_role = input$input_role, sample_id = input$sample_id,
       format = format, line = as.integer(line), caller = caller,
       original = original, selected_alt_index = as.integer(selected_alt_index),
       raw_info = raw_info, raw_format = raw_format, projected_info = projected_info,
       projected_format = projected_format)
}

event_record_row <- function(project, input, normalized, provenance, role, caller = NULL,
                             deferred = FALSE) {
  sample <- Filter(function(x) identical(x$sample_id, input$sample_id), project$samples)
  if (length(sample) != 1L) stop("Input sample binding is ambiguous")
  contexts <- Filter(function(x) input$sample_id %in% unlist(x$sample_bindings) &&
                       identical(x$source_subject_id, sample[[1L]]$patient_id), project$analysis_context)
  if (length(contexts) != 1L) stop("Input sample/context binding is ambiguous")
  ct <- contexts[[1L]]
  contract <- schema_contract()
  row <- stats::setNames(rep(list(NULL), length(contract$tables$events$fields)),
                  vapply(contract$tables$events$fields, `[[`, "", "name"))
  identity <- list(event_role = role, contig = normalized$contig,
                   pos = if (is.null(normalized$pos)) NULL else as.integer(normalized$pos),
                   ref = normalized$ref, alt = normalized$alt)
  if (deferred) identity$source_assertion <- normalized$source_assertion
  row$canonical_event <- list(namespace = if (deferred) "t006-deferred-source-v1" else
    "t006-primary-small-variant-v1", identity = identity)
  row$schema_version <- contract$schema_version
  row$analysis_id <- project$analysis_id
  row$analysis_context <- ct$context_id
  row$patient_id <- sample[[1L]]$patient_id
  row$sample_id <- input$sample_id
  row$event_key <- biological_id("event_key", list(genome_assembly = project$config$scientific$genome_assembly,
                                                     canonical_event = row$canonical_event), contract)
  row$event_id <- biological_id("event_id", list(context_id = ct$context_id,
    patient_id = row$patient_id, sample_id = row$sample_id, event_key = row$event_key), contract)
  row$event_type <- normalized$event_type
  row$event_role <- role
  row$genome_assembly <- project$config$scientific$genome_assembly
  row["contig"] <- list(normalized$contig)
  row["pos"] <- list(if (is.null(normalized$pos)) NULL else as.integer(normalized$pos))
  row["end"] <- list(if (is.null(normalized$end)) NULL else as.integer(normalized$end))
  row["ref"] <- list(normalized$ref)
  row["alt"] <- list(normalized$alt)
  if (deferred) {
    row["sv_type"] <- list(normalized$sv_type)
    row["mate_contig"] <- list(normalized$mate_contig)
    row["mate_pos"] <- list(normalized$mate_pos)
  }
  row["caller"] <- list(caller)
  row["caller_filter"] <- list(provenance$original$filter)
  row$original_record_ref <- paste0(input$input_id, ":", provenance$line)
  row$caller_evidence_json <- list()
  row$expression_json <- list()
  row$cn_by_breakpoint_json <- list()
  row$ccf_candidates_json <- list()
  row$phase_links_json <- list()
  row$event_details_json <- list(source_records = list(provenance))
  row$event_warning_codes <- list()
  row$event_detection_status <- "SOURCE_PROVIDED"
  row$reconstruction_status <- if (deferred) "UNSUPPORTED" else "NOT_RUN"
  row$candidate_generation_status <- "NOT_RUN"
  row
}

event_merge_rows <- function(rows) {
  if (!length(rows)) return(rows)
  by_id <- split(rows, vapply(rows, `[[`, "", "event_id"))
  out <- lapply(by_id, function(group) {
    row <- group[[1L]]
    records <- unlist(lapply(group, function(x) x$event_details_json$source_records), recursive = FALSE)
    keys <- vapply(records, canonical_bytes, "")
    records <- records[order(keys, method = "radix")]
    records <- records[!duplicated(vapply(records, canonical_bytes, ""))]
    row$event_details_json$source_records <- records
    callers <- sort(unique(vapply(records, function(x) if (is.null(x$caller)) "unknown" else x$caller, "")),
                    method = "radix")
    row["caller"] <- list(if (length(callers) == 1L && callers != "unknown") callers else NULL)
    row["caller_filter"] <- list(if (length(records) == 1L) records[[1L]]$original$filter else NULL)
    row["original_record_ref"] <- list(if (length(records) == 1L)
      paste0(records[[1L]]$input_id, ":", records[[1L]]$line) else NULL)
    row$caller_evidence_json <- lapply(records, function(x) list(input_id = x$input_id,
      caller = if (is.null(x$caller)) "unknown" else x$caller,
      source_record = paste0(x$input_id, ":", x$line)))
    row
  })
  unname(out[order(vapply(out, `[[`, "", "event_id"), method = "radix")])
}

event_vcf_meta <- function(line) {
  body <- sub("^##[^=]+=<", "", sub(">$", "", line))
  match <- regmatches(body, gregexpr('(?:^|,)([A-Za-z][A-Za-z0-9_]*)=("[^"]*"|[^,]*)', body,
                                      perl = TRUE))[[1L]]
  out <- list()
  for (m in match) {
    m <- sub("^,", "", m)
    i <- regexpr("=", m, fixed = TRUE)[1L]
    out[[substr(m, 1L, i - 1L)]] <- gsub('^"|"$', "", substring(m, i + 1L))
  }
  out
}

event_check_build_declaration <- function(value, assembly, label) {
  if (is.null(value) || !nzchar(value)) return(invisible(TRUE))
  if (grepl("GRCh37|hg19|b37", value, ignore.case = TRUE))
    stop(label, " contradicts selected reference: ", value)
  patch <- regmatches(value, regexpr("GRCh38[.]p[0-9]+", value, ignore.case = TRUE))
  if (length(patch) && nzchar(patch) && !identical(toupper(patch), toupper(assembly)))
    stop(label, " patch assembly contradicts selected reference: ", value)
  invisible(TRUE)
}

event_vcf_fields <- function(raw, defs, alt_index, alt_count) {
  if (is.null(raw) || raw == "." || !nzchar(raw)) return(list())
  fields <- event_vcf_info_map(raw)
  out <- list()
  for (key in names(fields)) {
    value <- fields[[key]]
    number <- defs[[key]]$Number
    if (is.null(number) || !number %in% c("A", "R")) next
    items <- event_vcf_allele_values(value, number, alt_count, paste("INFO", key))
    if (is.null(items)) next
    if (number == "A") out[[key]] <- items[[alt_index]]
    if (number == "R")
      out[[key]] <- list(ref = items[[1L]], alt = items[[alt_index + 1L]])
  }
  out
}

event_vcf_info_map <- function(raw) {
  if (is.null(raw) || raw == ".") return(list())
  out <- list()
  for (item in stringi::stri_split_fixed(raw, ";", omit_empty = FALSE)[[1L]]) {
    if (!nzchar(item)) stop("Malformed/duplicate VCF INFO key")
    p <- strsplit(item, "=", fixed = TRUE)[[1L]]
    if (!nzchar(p[[1L]]) || p[[1L]] %in% names(out)) stop("Malformed/duplicate VCF INFO key")
    out[[p[[1L]]]] <- if (length(p) > 1L) paste(p[-1L], collapse = "=") else TRUE
  }
  out
}

event_vcf_allele_values <- function(value, number, alt_count, label) {
  if (!is.character(value) || length(value) != 1L || !nzchar(value))
    stop("Malformed VCF ", label, " allele field")
  expected <- alt_count + if (number == "R") 1L else 0L
  # Preserve the released literal-dot projection when the expected length is
  # one; a whole missing array otherwise has no per-ALT projection.
  if (value == "." && expected > 1L) return(NULL)
  items <- stringi::stri_split_fixed(value, ",", omit_empty = FALSE)[[1L]]
  if (length(items) != expected || any(!nzchar(items)))
    stop("VCF ", label, " Number=", number, " cardinality mismatch")
  items
}

event_vcf_format <- function(raw) {
  if (is.null(raw)) return(list())
  keys <- stringi::stri_split_fixed(raw$keys, ":", omit_empty = FALSE)[[1L]]
  values <- stringi::stri_split_fixed(raw$sample, ":", omit_empty = FALSE)[[1L]]
  if (any(!nzchar(keys)) || anyDuplicated(keys) || any(!nzchar(values)) ||
      length(values) > length(keys)) stop("Malformed VCF FORMAT")
  # VCF permits omitted trailing sample fields. Keep raw_format unchanged.
  if (length(values) < length(keys)) values <- c(values, rep(".", length(keys) - length(values)))
  stats::setNames(as.list(values), keys)
}

event_parse_vcf <- function(project, input, reference) {
  con <- reference_connection(input$resolved_path); on.exit(close(con))
  info_defs <- list(); format_defs <- list(); contigs <- list(); sources <- character()
  declared_reference <- NULL; header <- NULL; line_number <- 0L; rows <- list(); n <- 0L
  repeat {
    chunk <- readLines(con, n = 10000L, warn = FALSE)
    if (!length(chunk)) break
    for (line in chunk) {
      line_number <- line_number + 1L
      if (startsWith(line, "##")) {
        if (startsWith(line, "##contig=<")) {
          d <- event_vcf_meta(line)
          if (is.null(d$ID) || d$ID %in% names(contigs)) stop("Malformed/duplicate VCF contig header")
          contigs[[d$ID]] <- d
        } else if (startsWith(line, "##INFO=<")) {
          d <- event_vcf_meta(line); if (!is.null(d$ID)) info_defs[[d$ID]] <- d
        } else if (startsWith(line, "##FORMAT=<")) {
          d <- event_vcf_meta(line); if (!is.null(d$ID)) format_defs[[d$ID]] <- d
        } else if (startsWith(line, "##reference=")) declared_reference <- substring(line, 13L)
        else if (startsWith(line, "##source=")) sources <- c(sources, substring(line, 10L))
        else if (startsWith(line, "##SomaticSeq=")) sources <- c(sources, paste0("SomaticSeq ", substring(line, 14L)))
        next
      }
      if (startsWith(line, "#CHROM")) {
        if (!is.null(header)) stop("Duplicate VCF column header")
        header <- strsplit(line, "\t", fixed = TRUE)[[1L]]
        if (length(header) < 8L || !identical(header[1:8],
          c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")))
          stop("Malformed VCF column header")
        if (length(header) > 9L && !identical(header[[10L]], input$sample_id))
          stop("VCF sample does not match authoritative input sample binding")
        if (length(header) > 10L) stop("Multi-sample VCF requires a single explicit sample input")
        if (!length(header) %in% c(8L, 10L) ||
            (length(header) == 10L && header[[9L]] != "FORMAT"))
          stop("Malformed VCF column header")
        event_check_build_declaration(declared_reference, reference$assembly, "VCF reference declaration")
        for (contig in names(contigs)) {
          d <- contigs[[contig]]; i <- match(contig, reference$fai$contig)
          event_check_build_declaration(d$assembly, reference$assembly, "VCF contig assembly")
          if (!is.na(i) && !is.null(d$length) &&
              (!grepl("^[0-9]+$", d$length) || as.numeric(d$length) != reference$fai$length[i]))
            stop("VCF contig length mismatch: ", contig)
        }
        next
      }
      if (startsWith(line, "#")) next
      if (is.null(header)) stop("VCF data before column header")
      fields <- stringi::stri_split_fixed(line, "\t", omit_empty = FALSE)[[1L]]
      if (length(fields) != length(header)) stop("Malformed VCF record at line ", line_number)
      if (any(!nzchar(fields))) stop("Malformed VCF record: empty field at line ", line_number)
      names(fields) <- header
      pos <- suppressWarnings(as.numeric(fields[["POS"]]))
      if (!grepl("^[0-9]+$", fields[["POS"]]) || !is.finite(pos) ||
          pos != floor(pos) || pos < 1 || pos > 2147483647)
        stop("Invalid VCF POS at line ", line_number)
      chrom <- fields[["#CHROM"]]
      if (!chrom %in% reference$fai$contig) stop("Contig absent from selected primary reference: ", chrom)
      ref <- fields[["REF"]]
      if (!grepl("^[ACGT]+$", ref)) stop("Malformed VCF REF at line ", line_number)
      alts <- stringi::stri_split_fixed(fields[["ALT"]], ",", omit_empty = FALSE)[[1L]]
      if (!length(alts) || any(!nzchar(alts))) stop("Malformed VCF ALT")
      caller <- if (length(unique(sources)) == 1L) unique(sources) else NULL
      raw_format <- if (length(fields) >= 10L) list(keys = fields[["FORMAT"]], sample = fields[[10L]]) else NULL
      format_values <- event_vcf_format(raw_format)
      event_vcf_info_map(fields[["INFO"]])
      called_alts <- seq_along(alts)
      if (identical(input$input_role, "germline_vcf") && !is.null(raw_format)) {
        gt <- format_values$GT
        if (!is.null(gt)) {
          if (!grepl("^[0-9]+([/|][0-9]+)*$", gt)) stop("Malformed germline GT")
          tokens <- strsplit(gt, "[/|]", perl = TRUE)[[1L]]
          if (any(!grepl("^[0-9]+$", tokens))) stop("Germline GT contains unresolved allele")
          indexes <- suppressWarnings(as.numeric(tokens))
          if (any(!is.finite(indexes)) || any(indexes > length(alts)))
            stop("Germline GT allele index exceeds ALT count")
          called_alts <- sort(unique(indexes[indexes > 0L]))
        }
      }
      projected_format <- function(index) {
        if (is.null(raw_format)) return(list())
        out <- list()
        for (key in names(format_values)) {
          number <- format_defs[[key]]$Number
          if (is.null(number) || !number %in% c("A", "R")) next
          items <- event_vcf_allele_values(format_values[[key]], number, length(alts), paste("FORMAT", key))
          if (is.null(items)) next
          if (number == "A") out[[key]] <- items[[index]]
          if (number == "R") out[[key]] <- list(ref = items[[1L]], alt = items[[index + 1L]])
        }
        out
      }
      original <- list(chrom = chrom, pos = as.integer(pos), id = fields[["ID"]], ref = ref,
                       alts = as.list(alts), filter = fields[["FILTER"]], qual = fields[["QUAL"]])
      for (i in seq_along(alts)) {
        if (!i %in% called_alts) next
        alt <- alts[[i]]
        if (alt == "*") stop("VCF spanning-deletion '*' allele is unsupported at line ", line_number)
        is_bnd <- grepl("[", alt, fixed = TRUE) || grepl("]", alt, fixed = TRUE)
        deferred <- grepl("^<[^>]+>$", alt) || is_bnd
        if (!deferred && !grepl("^[ACGT]+$", alt)) stop("Unsupported VCF ALT syntax at line ", line_number)
        if (deferred) {
          event_validate_ref(reference, chrom, pos, ref)
          info_map <- event_vcf_info_map(fields[["INFO"]])
          sv_end <- if (is.null(info_map$END)) pos + nchar(ref) - 1L else
            suppressWarnings(as.numeric(info_map$END))
          if (!is.finite(sv_end) || sv_end != floor(sv_end) || sv_end < pos ||
              sv_end > reference$fai$length[match(chrom, reference$fai$contig)])
            stop("Invalid symbolic event END at line ", line_number)
          mate <- if (is_bnd) {
            pieces <- strsplit(gsub("]", "[", alt, fixed = TRUE), "[", fixed = TRUE)[[1L]]
            parts <- strsplit(pieces[grepl(":", pieces, fixed = TRUE)], ":", fixed = TRUE)
            if (length(parts) != 1L || length(parts[[1L]]) != 2L ||
                !grepl("^[0-9]+$", parts[[1L]][[2L]]))
              stop("Malformed BND mate syntax at line ", line_number)
            parts[[1L]]
          } else character()
          if (is_bnd) {
            mate_pos <- suppressWarnings(as.numeric(mate[[2L]]))
            mate_index <- match(mate[[1L]], reference$fai$contig)
            if (is.na(mate_index)) stop("BND mate contig absent from selected primary reference at line ", line_number)
            if (!is.finite(mate_pos) || mate_pos != floor(mate_pos) || mate_pos < 1 ||
                mate_pos > reference$fai$length[[mate_index]] || mate_pos > 2147483647)
              stop("BND mate position outside selected primary reference at line ", line_number)
          }
          normalized <- list(contig = chrom, pos = as.integer(pos), end = as.integer(pos + nchar(ref) - 1L),
                             ref = ref, alt = alt, event_type = if (is_bnd) "BND" else "SV",
                             changed = FALSE)
          normalized$end <- as.integer(sv_end)
          normalized$sv_type <- if (is_bnd) "BND" else gsub("[<>]", "", alt)
          normalized$mate_contig <- if (is_bnd) mate[[1L]] else NULL
          normalized$mate_pos <- if (is_bnd) as.integer(mate[[2L]]) else NULL
          # Deferred identity denotes a source assertion. No SV equivalence is
          # claimed across records until T022 establishes structural semantics.
          normalized$source_assertion <- list(input_id = input$input_id,
            source_record_id = fields[["ID"]], info = info_map)
        } else normalized <- event_normalize(reference, chrom, pos, ref, alt)
        prov <- event_source_record(input, "VCF", line_number, original, i, caller,
          fields[["INFO"]], raw_format,
          event_vcf_fields(fields[["INFO"]], info_defs, i, length(alts)), projected_format(i))
        prov$declared_reference <- declared_reference
        if (identical(input$input_role, "germline_vcf"))
          prov$genotype_status <- if (is.null(raw_format) ||
            !"GT" %in% strsplit(raw_format$keys, ":", fixed = TRUE)[[1L]])
            "UNKNOWN_NO_GT" else "RAW_GT_PRESENT"
        prov$header_sources <- as.list(sources)
        prov$normalization <- list(source = list(pos = as.integer(pos), ref = ref, alt = alt),
                                   canonical = normalized[c("pos", "ref", "alt")], changed = normalized$changed)
        role <- if (startsWith(input$input_role, "somatic")) "somatic" else "germline"
        n <- n + 1L; rows[[n]] <- event_record_row(project, input, normalized, prov, role, caller, deferred)
      }
    }
  }
  if (is.null(header)) stop("Missing VCF column header")
  rows
}

event_parse_maf <- function(project, input, reference) {
  con <- reference_connection(input$resolved_path); on.exit(close(con))
  header <- NULL; rows <- list(); n <- 0L; line_number <- 0L
  required <- c("Chromosome", "Start_Position", "End_Position", "Reference_Allele",
                "Tumor_Seq_Allele1", "Tumor_Seq_Allele2")
  repeat {
    chunk <- readLines(con, n = 10000L, warn = FALSE)
    if (!length(chunk)) break
    for (line in chunk) {
      line_number <- line_number + 1L
      if (startsWith(line, "#") || !nzchar(line)) next
      fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
      if (is.null(header)) {
        header <- fields
        if (anyDuplicated(header) || length(setdiff(required, header))) stop("Unsupported MAF profile: required columns missing")
        next
      }
      if (length(fields) != length(header)) stop("Malformed MAF row at line ", line_number)
      names(fields) <- header
      chrom <- fields[["Chromosome"]]
      if (!chrom %in% reference$fai$contig) stop("Contig absent from selected primary reference: ", chrom)
      start <- suppressWarnings(as.numeric(fields[["Start_Position"]]))
      end <- suppressWarnings(as.numeric(fields[["End_Position"]]))
      if (any(!is.finite(c(start, end))) || any(c(start, end) != floor(c(start, end))) ||
          start < 1 || end < 0 || start > 2147483647 || end > 2147483647)
        stop("Invalid MAF coordinates at line ", line_number)
      if ("NCBI_Build" %in% header && fields[["NCBI_Build"]] != ".")
        event_check_build_declaration(fields[["NCBI_Build"]], reference$assembly, "MAF build")
      if ("Tumor_Sample_Barcode" %in% header && nzchar(fields[["Tumor_Sample_Barcode"]]) &&
          fields[["Tumor_Sample_Barcode"]] != input$sample_id)
        stop("MAF sample does not match authoritative input sample binding")
      source_ref <- fields[["Reference_Allele"]]
      alleles <- fields[c("Tumor_Seq_Allele1", "Tumor_Seq_Allele2")]
      if (!source_ref %in% c("-") && !grepl("^[ACGT]+$", source_ref)) stop("Malformed MAF REF")
      if (any(!vapply(alleles, function(x) x == "-" || grepl("^[ACGT]+$", x), TRUE)))
        stop("Malformed MAF tumor allele")
      alts <- unique(unname(alleles[alleles != source_ref]))
      if (!length(alts)) stop("MAF no-op tumor alleles")
      if (start > end) stop("MAF Start_Position exceeds End_Position")
      if (source_ref == "-") {
        if (end != start + 1L) stop("MAF insertion requires End_Position = Start_Position + 1")
      } else if (end != start + nchar(source_ref) - 1L) stop("MAF REF length/interval conflict")
      original <- as.list(fields)
      original$filter <- if ("FILTER" %in% header) fields[["FILTER"]] else NULL
      for (i in seq_along(alts)) {
        source_alt <- alts[[i]]
        variant_type <- if (source_ref == "-" || nchar(source_ref) < nchar(source_alt)) "INS" else
          if (source_alt == "-" || nchar(source_ref) > nchar(source_alt)) "DEL" else
            if (nchar(source_ref) > 3L) "ONP" else c("SNP", "DNP", "TNP")[[nchar(source_ref)]]
        declared_type <- if ("Variant_Type" %in% header) fields[["Variant_Type"]] else ""
        if (nzchar(declared_type) && declared_type != "." && declared_type != variant_type)
          stop("MAF Variant_Type contradicts genomic allele at line ", line_number)
        if (source_ref == "-") {
          anchor <- event_reference_fetch(reference, chrom, start, start)
          pos <- start; ref <- anchor; alt <- paste0(anchor, source_alt)
        } else if (source_alt == "-") {
          if (start <= 1L) stop("MAF deletion at contig start has no left anchor")
          anchor <- event_reference_fetch(reference, chrom, start - 1L, start - 1L)
          pos <- start - 1L; ref <- paste0(anchor, source_ref); alt <- anchor
        } else {
          pos <- start; ref <- source_ref; alt <- source_alt
        }
        normalized <- event_normalize(reference, chrom, pos, ref, alt)
        prov <- event_source_record(input, "MAF", line_number, original, i,
          if ("Caller" %in% header && nzchar(fields[["Caller"]])) fields[["Caller"]] else NULL)
        prov$normalization <- list(source = list(start = as.integer(start), end = as.integer(end),
          ref = source_ref, alt = source_alt), canonical = normalized[c("pos", "ref", "alt")],
          changed = normalized$changed || source_ref == "-" || source_alt == "-")
        role <- if (startsWith(input$input_role, "somatic")) "somatic" else "germline"
        n <- n + 1L; rows[[n]] <- event_record_row(project, input, normalized, prov, role, prov$caller)
      }
    }
  }
  if (is.null(header)) stop("Missing MAF header")
  rows
}

event_parse_annovar <- function(project, input, reference) {
  con <- reference_connection(input$resolved_path); on.exit(close(con))
  header <- strsplit(readLines(con, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1L]]
  if (anyDuplicated(header) || length(setdiff(c("Chr", "Start", "End", "Ref", "Alt"), header)))
    stop("Unsupported ANNOVAR table: required columns missing")
  rows <- list(); n <- 0L; line_number <- 1L
  repeat {
    chunk <- readLines(con, n = 10000L, warn = FALSE)
    if (!length(chunk)) break
    for (line in chunk) {
      line_number <- line_number + 1L
      fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
      if (length(fields) != length(header)) stop("Malformed ANNOVAR row at line ", line_number)
      names(fields) <- header
      chrom <- fields[["Chr"]]
      if (!chrom %in% reference$fai$contig) stop("Contig absent from selected primary reference: ", chrom)
      start <- suppressWarnings(as.numeric(fields[["Start"]]))
      end <- suppressWarnings(as.numeric(fields[["End"]]))
      if (any(!is.finite(c(start, end))) || any(c(start, end) != floor(c(start, end))) ||
          start < 1 || end < 0 || start > 2147483647 || end > 2147483647)
        stop("Invalid ANNOVAR coordinates at line ", line_number)
      source_ref <- fields[["Ref"]]; source_alt <- fields[["Alt"]]
      if (source_ref == "0")
        stop("ANNOVAR Ref=0 lacks an independently REF-validatable small-event allele at line ", line_number)
      if (!(source_ref == "-" || grepl("^[ACGT]+$", source_ref)) ||
          !(source_alt == "-" || grepl("^[ACGT]+$", source_alt)) || source_ref == source_alt)
        stop("Malformed ANNOVAR alleles at line ", line_number)
      if (source_ref == "-") {
        if (start != end) stop("ANNOVAR insertion requires Start = End")
        anchor <- event_reference_fetch(reference, chrom, start, start)
        pos <- start; ref <- anchor; alt <- paste0(anchor, source_alt)
      } else if (source_alt == "-") {
        if (end != start + nchar(source_ref) - 1L || start <= 1L)
          stop("ANNOVAR deletion interval conflict")
        anchor <- event_reference_fetch(reference, chrom, start - 1L, start - 1L)
        pos <- start - 1L; ref <- paste0(anchor, source_ref); alt <- anchor
      } else {
        if (end != start + nchar(source_ref) - 1L) stop("ANNOVAR interval/REF conflict")
        pos <- start; ref <- source_ref; alt <- source_alt
      }
      normalized <- event_normalize(reference, chrom, pos, ref, alt)
      original <- as.list(fields); original$filter <- NULL
      prov <- event_source_record(input, "ANNOVAR", line_number, original)
      prov$normalization <- list(source = list(start = as.integer(start), end = as.integer(end),
        ref = source_ref, alt = source_alt), canonical = normalized[c("pos", "ref", "alt")],
        changed = normalized$changed || source_ref == "-" || source_alt == "-")
      role <- if (startsWith(input$input_role, "somatic")) "somatic" else "germline"
      n <- n + 1L; rows[[n]] <- event_record_row(project, input, normalized, prov, role)
    }
  }
  rows
}

event_parse_snp_table <- function(project, input, reference) {
  con <- reference_connection(input$resolved_path); on.exit(close(con))
  header <- strsplit(readLines(con, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1L]]
  if (anyDuplicated(header) || length(setdiff(c("CHROM", "POS", "REF", "ALT"), header)))
    stop("Unsupported legacy SNP table")
  rows <- list(); n <- 0L; line_number <- 1L
  repeat {
    chunk <- readLines(con, n = 10000L, warn = FALSE)
    if (!length(chunk)) break
    for (line in chunk) {
      line_number <- line_number + 1L
      fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
      if (length(fields) != length(header)) stop("Malformed SNP table row")
      names(fields) <- header
      pos <- suppressWarnings(as.numeric(fields[["POS"]]))
      if (!is.finite(pos) || pos != floor(pos) || pos < 1 || pos > 2147483647)
        stop("Invalid SNP table POS")
      alts <- strsplit(fields[["ALT"]], ",", fixed = TRUE)[[1L]]
      if (!length(alts) || any(!grepl("^[ACGT]+$", alts))) stop("Malformed SNP table ALT")
      original <- as.list(fields)
      original$filter <- if ("FILTER" %in% header) fields[["FILTER"]] else NULL
      for (i in seq_along(alts)) {
        normalized <- event_normalize(reference, fields[["CHROM"]], as.integer(pos),
          fields[["REF"]], alts[[i]])
        prov <- event_source_record(input, "legacy_SNP_table", line_number, original, i)
        prov$genotype_status <- "UNKNOWN_NO_SAMPLE_GT"
        prov$normalization <- list(source = list(pos = as.integer(pos), ref = fields[["REF"]],
          alt = alts[[i]]), canonical = normalized[c("pos", "ref", "alt")], changed = normalized$changed)
        n <- n + 1L; rows[[n]] <- event_record_row(project, input, normalized, prov, "germline")
      }
    }
  }
  rows
}

event_parse_bnd_table <- function(project, input, reference) {
  con <- reference_connection(input$resolved_path); on.exit(close(con))
  header <- strsplit(readLines(con, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1L]]
  required <- c("Chr", "Start", "End", "Ref", "Alt", "mateID")
  if (anyDuplicated(header) || length(setdiff(required, header)))
    stop("Unsupported legacy BND table")
  rows <- list(); n <- 0L; line_number <- 1L
  repeat {
    chunk <- readLines(con, n = 10000L, warn = FALSE)
    if (!length(chunk)) break
    for (line in chunk) {
      line_number <- line_number + 1L
      fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
      if (length(fields) != length(header)) stop("Malformed BND table row")
      names(fields) <- header
      pos <- suppressWarnings(as.numeric(fields[["Start"]]))
      end <- suppressWarnings(as.numeric(fields[["End"]]))
      if (!is.finite(pos) || pos != floor(pos) || pos < 1 || pos > 2147483647 ||
          !identical(pos, end) || nchar(fields[["Ref"]]) != 1L)
        stop("Invalid BND table local coordinates")
      chrom <- fields[["Chr"]]; ref <- fields[["Ref"]]; alt <- fields[["Alt"]]
      event_validate_ref(reference, chrom, as.integer(pos), ref)
      pieces <- strsplit(gsub("]", "[", alt, fixed = TRUE), "[", fixed = TRUE)[[1L]]
      parts <- strsplit(pieces[grepl(":", pieces, fixed = TRUE)], ":", fixed = TRUE)
      if (length(parts) != 1L || length(parts[[1L]]) != 2L ||
          !grepl("^[0-9]+$", parts[[1L]][[2L]]) ||
          !(grepl("[", alt, fixed = TRUE) || grepl("]", alt, fixed = TRUE)))
        stop("Malformed BND table mate syntax")
      mate_index <- match(parts[[1L]][[1L]], reference$fai$contig)
      mate_pos <- suppressWarnings(as.numeric(parts[[1L]][[2L]]))
      if (is.na(mate_index)) stop("BND mate contig absent from selected primary reference")
      if (!is.finite(mate_pos) || mate_pos < 1 || mate_pos > 2147483647 ||
          mate_pos > reference$fai$length[[mate_index]]) stop("BND mate position outside selected primary reference")
      original <- as.list(fields); original$filter <- NULL
      normalized <- list(contig = chrom, pos = as.integer(pos), end = as.integer(end),
        ref = ref, alt = alt, event_type = "BND", changed = FALSE, sv_type = "BND",
        mate_contig = parts[[1L]][[1L]], mate_pos = as.integer(mate_pos),
        source_assertion = list(input_id = input$input_id, mateID = fields[["mateID"]],
          source_record = original))
      prov <- event_source_record(input, "legacy_BND_table", line_number, original)
      n <- n + 1L; rows[[n]] <- event_record_row(project, input, normalized, prov, "somatic", deferred = TRUE)
    }
  }
  rows
}

event_parse_sequence <- function(project, input) {
  con <- reference_connection(input$resolved_path); on.exit(close(con))
  header <- strsplit(readLines(con, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1L]]
  fragment <- identical(input$input_role, "sequence_fragment")
  required <- if (fragment) c("label", "sequence") else c("label", "input_nm_id", "reading_frame")
  allowed <- c(required, "group_id", "reference_nm_id", "reference_gene_symbol", "source_annotation", "sequence")
  if (anyDuplicated(header) || length(setdiff(required, header)) || length(setdiff(header, allowed)))
    stop("Unsupported legacy sequence source table")
  rows <- list(); n <- 0L
  repeat {
    chunk <- readLines(con, n = 10000L, warn = FALSE)
    if (!length(chunk)) break
    for (line in chunk) {
      fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
      if (length(fields) != length(header)) stop("Malformed legacy sequence row")
      names(fields) <- header
      if (!nzchar(fields[["label"]])) stop("Missing sequence label")
      if (fragment && !grepl("^[A-Za-z*]+$", fields[["sequence"]])) stop("Malformed sequence fragment")
      if (!fragment && (!nzchar(fields[["input_nm_id"]]) ||
          !fields[["reading_frame"]] %in% c("1", "2", "3")))
        stop("Invalid entire-region source")
      original <- as.list(fields); original$source_mode <- if (fragment) "fragment" else "entire-region"
      original$filter <- NULL
      prov <- event_source_record(input, "legacy_sequence_TSV", n + 2L, original)
      normalized <- list(contig = NULL, pos = NULL, end = NULL, ref = NULL, alt = NULL,
                         event_type = "SEQUENCE")
      row <- event_record_row(project, input, normalized, prov, "source", deferred = TRUE)
      row$canonical_event <- list(namespace = "t006-legacy-sequence-source-v1", identity = list(
        source_mode = original$source_mode, input_id = input$input_id, label = fields[["label"]],
        source_record = original))
      row["genome_assembly"] <- list(NULL)
      row["pos"] <- list(NULL); row["end"] <- list(NULL)
      row$event_key <- biological_id("event_key", list(genome_assembly = NULL,
        canonical_event = row$canonical_event), schema_contract())
      row$event_id <- biological_id("event_id", list(context_id = row$analysis_context,
        patient_id = row$patient_id, sample_id = row$sample_id, event_key = row$event_key), schema_contract())
      n <- n + 1L; rows[[n]] <- row
    }
  }
  rows
}

event_parse_vep <- function(project, input, reference) {
  con <- reference_connection(input$resolved_path); on.exit(close(con))
  header <- strsplit(readLines(con, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1L]]
  required <- c("Uploaded_variation", "Location", "Allele", "Consequence")
  if (anyDuplicated(header) || length(setdiff(required, header)))
    stop("Unsupported VEP table: required columns missing")
  rows <- list(); n <- 0L; line_number <- 1L
  repeat {
    chunk <- readLines(con, n = 10000L, warn = FALSE)
    if (!length(chunk)) break
    for (line in chunk) {
      line_number <- line_number + 1L
      fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
      if (length(fields) != length(header)) stop("Malformed VEP row at line ", line_number)
      names(fields) <- header
      location <- strsplit(fields[["Location"]], ":", fixed = TRUE)[[1L]]
      if (length(location) != 2L || !grepl("^[0-9]+(-[0-9]+)?$", location[[2L]]) ||
          !grepl("^(-|[ACGT]+)$", fields[["Allele"]]))
        stop("VEP source cannot be mapped to a genomic sequence assertion at line ", line_number)
      chrom <- location[[1L]]
      bounds <- suppressWarnings(as.numeric(strsplit(location[[2L]], "-", fixed = TRUE)[[1L]]))
      if (length(bounds) == 1L) bounds <- rep(bounds, 2L)
      start <- bounds[[1L]]; end <- bounds[[2L]]
      if (any(!is.finite(bounds)) || any(bounds > 2147483647) || start < 1 || end < 0 ||
          !(start <= end || start == end + 1L))
        stop("Invalid VEP Location at line ", line_number)
      # VEP's legacy table has no genomic REF. Resolve it from the selected
      # reference, and retain that fact separately from the source columns.
      source_alt <- fields[["Allele"]]
      if (start == end + 1L) {
        if (end < 1L || source_alt == "-") stop("Invalid VEP insertion assertion")
        anchor <- event_reference_fetch(reference, chrom, as.integer(end), as.integer(end))
        pos <- end; ref <- anchor; alt <- paste0(anchor, source_alt)
      } else if (source_alt == "-") {
        if (start <= 1L) stop("VEP deletion at contig start has no left anchor")
        anchor <- event_reference_fetch(reference, chrom, as.integer(start - 1L), as.integer(start - 1L))
        sequence <- event_reference_fetch(reference, chrom, as.integer(start), as.integer(end))
        pos <- start - 1L; ref <- paste0(anchor, sequence); alt <- anchor
      } else {
        pos <- start
        ref <- event_reference_fetch(reference, chrom, as.integer(start), as.integer(end))
        alt <- source_alt
      }
      normalized <- event_normalize(reference, chrom, as.integer(pos), ref, alt)
      original <- as.list(fields); original$filter <- NULL
      prov <- event_source_record(input, "VEP", line_number, original)
      prov$reference_allele_origin <- "selected_T005_reference_not_source_column"
      prov$normalization <- list(source = list(location = fields[["Location"]], alt = source_alt),
        canonical = normalized[c("pos", "ref", "alt")],
        changed = normalized$changed || start == end + 1L || source_alt == "-")
      role <- if (startsWith(input$input_role, "somatic")) "somatic" else "germline"
      n <- n + 1L; rows[[n]] <- event_record_row(project, input, normalized, prov, role)
    }
  }
  rows
}

event_parse_canonical <- function(project, input, reference) {
  # Input-only TSV profile: every T002 field is a column; structured cells are JSON.
  # This does not define final public result serialization.
  con <- reference_connection(input$resolved_path); on.exit(close(con))
  header <- strsplit(readLines(con, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1L]]
  contract <- schema_contract(); fields <- contract$tables$events$fields
  allowed <- vapply(fields, `[[`, "", "name")
  if (anyDuplicated(header) || length(setdiff(allowed, header)) || length(setdiff(header, allowed)))
    stop("Canonical event TSV must contain exactly T002 event columns")
  defs <- stats::setNames(fields, allowed); rows <- list(); n <- 0L
  repeat {
    chunk <- readLines(con, n = 10000L, warn = FALSE)
    if (!length(chunk)) break
    for (line in chunk) {
      values <- strsplit(line, "\t", fixed = TRUE)[[1L]]
      if (length(values) != length(header)) stop("Malformed canonical TSV row")
      names(values) <- header
      row <- stats::setNames(rep(list(NULL), length(allowed)), allowed)
      for (key in allowed) {
        value <- values[[key]]; type <- defs[[key]]$type
        if (value == "\\N") next
        row[key] <- list(switch(type,
          object = jsonlite::fromJSON(value, simplifyVector = FALSE),
          array = jsonlite::fromJSON(value, simplifyVector = FALSE),
          json = jsonlite::fromJSON(value, simplifyVector = FALSE),
          integer = as.integer(value), number = as.numeric(value),
          boolean = if (value == "true") TRUE else if (value == "false") FALSE else stop("Invalid boolean"),
          value))
      }
      validate_table(list(row), "events", contract)
      if (!identical(row$event_key, biological_id("event_key",
            row[c("genome_assembly", "canonical_event")], contract)) ||
          !identical(row$event_id, biological_id("event_id", list(context_id = row$analysis_context,
            patient_id = row$patient_id, sample_id = row$sample_id,
            event_key = row$event_key), contract))) stop("Canonical TSV source ID mismatch")
      sample <- Filter(function(x) identical(x$sample_id, input$sample_id), project$samples)
      if (length(sample) != 1L || !identical(row$sample_id, input$sample_id) ||
          !identical(row$patient_id, sample[[1L]]$patient_id) ||
          !identical(row$genome_assembly, reference$assembly)) stop("Canonical event input binding mismatch")
      if (!row$event_type %in% c("SNV", "MNV", "INS", "DEL", "INDEL"))
        stop("Canonical TSV profile supports sequence-resolved small events only")
      normalized <- event_normalize(reference, row$contig, row$pos, row$ref, row$alt)
      if (normalized$changed || !identical(normalized$event_type, row$event_type) ||
          !identical(normalized$end, row$end)) stop("Canonical TSV event is not normalized")
      if (!row$event_role %in% c("somatic", "germline"))
        stop("Canonical TSV profile requires somatic or germline event role")
      expected <- list(namespace = "t006-primary-small-variant-v1", identity = list(
        event_role = row$event_role, contig = normalized$contig, pos = normalized$pos,
        ref = normalized$ref, alt = normalized$alt))
      if (!identical(canonical_bytes(row$canonical_event), canonical_bytes(expected)))
        stop("Canonical TSV canonical_event disagrees with normalized allele")
      original <- list(source_analysis_id = row$analysis_id,
        source_analysis_context = row$analysis_context, source_patient_id = row$patient_id,
        source_sample_id = row$sample_id, source_canonical_event = row$canonical_event,
        source_event_id = row$event_id,
        source_event_key = row$event_key, source_details = row$event_details_json,
        chrom = row$contig, pos = row$pos, ref = row$ref,
        alts = list(row$alt), filter = row$caller_filter)
      prov <- event_source_record(input, "canonical_T002_TSV", n + 2L,
                                  original, caller = row$caller)
      prov$normalization <- list(source = list(pos = row$pos, ref = row$ref, alt = row$alt),
        canonical = normalized[c("pos", "ref", "alt")], changed = FALSE)
      n <- n + 1L; rows[[n]] <- event_record_row(project, input, normalized, prov,
                                                 row$event_role, row$caller)
    }
  }
  validate_table(rows, "events", contract)
  rows
}

event_import_project <- function(project) {
  if (!length(project$inputs)) stop("Real event run requires event input; no-input production is not implemented in T004")
  roles <- vapply(project$inputs, `[[`, "", "input_role")
  if (any(!roles %in% c(event_input_roles(), source_input_roles(), long_read_roles(), recognition_roles(), presentation_roles(), longitudinal_roles())))
    stop("Production input stage for this role is not implemented in T004/T006")
  inputs <- Filter(function(x) x$input_role %in% event_input_roles(), project$inputs)
  if (!length(inputs)) return(list())
  bundle <- reference_bundle(project)
  reference <- event_reference_open(bundle); on.exit(event_reference_close(reference))
  rows <- list()
  for (input in inputs) {
    current <- integrity_record(input$resolved_path, project$config$checksum, input$integrity$declared)
    if (!identical(current, input$integrity)) stop("Event source changed since project binding: ", input$input_id)
    parsed <- switch(input$input_role,
      somatic_vcf = event_parse_vcf(project, input, reference),
      germline_vcf = event_parse_vcf(project, input, reference),
      somatic_maf = event_parse_maf(project, input, reference),
      germline_maf = event_parse_maf(project, input, reference),
      canonical_events = event_parse_canonical(project, input, reference),
      somatic_annovar = event_parse_annovar(project, input, reference),
      germline_annovar = event_parse_annovar(project, input, reference),
      somatic_vep = event_parse_vep(project, input, reference),
      germline_vep = event_parse_vep(project, input, reference),
      sequence_fragment = event_parse_sequence(project, input),
      sequence_entire_region = event_parse_sequence(project, input),
      germline_snp_table = event_parse_snp_table(project, input, reference),
      somatic_bnd_table = event_parse_bnd_table(project, input, reference))
    rows <- c(rows, parsed)
  }
  event_merge_rows(rows)
}
