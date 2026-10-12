reconstruction_policy <- function() "modern-primary-cds-v1"

reconstruction_codons <- local({
  bases <- strsplit("TCAG", "", fixed = TRUE)[[1L]]
  codons <- unlist(lapply(bases, function(first) unlist(lapply(bases,
    function(second) paste0(first, second, bases)), use.names = FALSE)), use.names = FALSE)
  # NCBI translation table 1, in T/C/A/G order at each codon position.
  aa <- strsplit(paste0("FFLLSSSSYY**CC*W", "LLLLPPPPHHQQRRRR",
                        "IIIMTTTTNNKKSSRR", "VVVVAAAADDEEGGGG"), "", fixed = TRUE)[[1L]]
  stats::setNames(aa, codons)
})

reconstruction_translate <- function(cds) {
  if (!grepl("^[ACGT]+$", cds) || nchar(cds) %% 3L) return(NULL)
  starts <- seq.int(1L, nchar(cds), 3L)
  codons <- substring(cds, starts, starts + 2L)
  paste0(unname(reconstruction_codons[codons]), collapse = "")
}

# Read only requested catalog records. The full ALL-region FASTA is never kept.
reconstruction_fasta_select <- function(path, ids) {
  ids <- unique(ids); out <- list()
  if (!length(ids)) return(out)
  con <- reference_connection(path); on.exit(close(con))
  current <- NULL; pieces <- character()
  flush <- function() {
    if (!is.null(current) && current %in% ids) out[[current]] <<- paste0(pieces, collapse = "")
  }
  repeat {
    lines <- readLines(con, n = 20000L, warn = FALSE)
    if (!length(lines)) break
    for (line in lines) if (startsWith(line, ">")) {
      flush(); current <- strsplit(substring(line, 2L), "|", fixed = TRUE)[[1L]][1L]
      pieces <- character()
    } else if (!is.null(current) && current %in% ids) pieces <- c(pieces, toupper(line))
  }
  flush(); out
}

reconstruction_candidate_ids <- function(bundle, rows) {
  small <- Filter(function(x) x$event_type %in% c("SNV", "MNV", "INS", "DEL", "INDEL") &&
                    !is.null(x$contig), rows)
  if (!length(small)) return(character())
  intervals <- lapply(small, function(row) {
    e <- reconstruction_edit_core(row)
    if (is.null(e)) return(NULL)
    list(contig = row$contig, start = if (e$insertion) e$start - 1L else e$start,
         end = if (e$insertion) e$start else e$end)
  })
  intervals <- Filter(Negate(is.null), intervals)
  if (!length(intervals)) return(character())
  by_contig <- lapply(split(intervals, vapply(intervals, `[[`, "", "contig")), function(x) {
    start <- vapply(x, `[[`, 0L, "start"); end <- vapply(x, `[[`, 0L, "end")
    order <- order(start, method = "radix")
    list(start = start[order], max_end = cummax(end[order]))
  })
  con <- reference_connection(bundle$resources$annotation_gtf$resolved_path); on.exit(close(con))
  ids <- character()
  repeat {
    lines <- readLines(con, n = 10000L, warn = FALSE)
    if (!length(lines)) break
    lines <- lines[!startsWith(lines, "#")]
    if (!length(lines)) next
    parts <- stringi::stri_split_fixed(lines, "\t", omit_empty = FALSE)
    for (p in parts) {
      if (!p[[3L]] %in% c("CDS", "stop_codon") ||
          is.null(by_contig[[p[[1L]]]])) next
      start <- as.integer(p[[4L]]); end <- as.integer(p[[5L]])
      index <- by_contig[[p[[1L]]]]
      k <- findInterval(end, index$start)
      relevant <- k > 0L && index$max_end[[k]] >= start
      if (!relevant) next
      id <- stringi::stri_match_first_regex(p[[9L]], 'transcript_id "([^"]+)";')[, 2L]
      if (!is.na(id)) ids <- c(ids, id)
    }
  }
  sort(unique(ids), method = "radix")
}

reconstruction_models <- function(bundle, ids) {
  if (!length(ids)) return(list())
  features <- reference_transcript_features(bundle, ids)
  groups <- split(features, features$transcript_id)
  lapply(groups, function(x) {
    cds <- x[x$feature == "CDS", , drop = FALSE]
    if (!nrow(cds)) return(NULL)
    tx <- x[x$feature == "transcript", , drop = FALSE]
    if (nrow(tx) != 1L) stop("Transcript declaration missing/duplicated")
    if (!identical(tx$transcript_type[[1L]], "protein_coding") &&
        !identical(tx$transcript_biotype[[1L]], "protein_coding")) return(NULL)
    cds <- cds[order(cds$start, method = "radix"), , drop = FALSE]
    exons <- x[x$feature == "exon", , drop = FALSE]
    exons <- exons[order(exons$start, method = "radix"), , drop = FALSE]
    model_reason <- if (!nrow(exons) ||
      any(exons$start[-1L] <= utils::head(exons$end, -1L)))
      "EXON_STRUCTURE_UNASSESSABLE" else if
      (any(cds$start[-1L] <= utils::head(cds$end, -1L)))
      "CDS_STRUCTURE_UNASSESSABLE" else NULL
    protein <- unique(stats::na.omit(c(cds$protein_id, tx$protein_id)))
    if (length(protein) > 1L) stop("Contradictory transcript/protein relationship")
    list(transcript_id = tx$transcript_id[[1L]],
         protein_id = if (length(protein)) protein[[1L]] else NA_character_,
         gene_id = tx$gene_id[[1L]], contig = tx$seqname[[1L]],
         strand = tx$strand[[1L]], segments = cds[, c("start", "end", "phase")],
         exons = exons[, c("start", "end"), drop = FALSE],
         model_reason = model_reason,
         start_codon = any(x$feature == "start_codon"),
         stop_codon = any(x$feature == "stop_codon"),
         stop_segments = x[x$feature == "stop_codon", c("start", "end"), drop = FALSE],
         incomplete = any(grepl("cds_(start|end)_NF", x$raw_attributes, perl = TRUE)),
         mane = unique(bundle$mane$MANE_status[bundle$mane$Ensembl_nuc == tx$transcript_id[[1L]]]))
  })
}

reconstruction_reference <- function(model, reference, protein_catalog, sec_ids,
                                     transcript_catalog = NULL) {
  if (!is.null(model$model_reason)) return(list(cds = NULL, protein = NULL,
    transcript = NULL, terminal_codon = NULL, reason = model$model_reason,
    forward = NULL, exon_forward = NULL))
  seg <- model$segments
  forward <- vapply(seq_len(nrow(seg)), function(i)
    event_reference_fetch(reference, model$contig, seg$start[i], seg$end[i]), "")
  ordered <- if (model$strand == "-") rev(seq_len(nrow(seg))) else seq_len(nrow(seg))
  coding <- if (model$strand == "-")
    vapply(forward[ordered], reconstruction_reverse_complement, "") else forward[ordered]
  cds <- paste0(coding, collapse = "")
  exon_forward <- vapply(seq_len(nrow(model$exons)), function(i)
    event_reference_fetch(reference, model$contig, model$exons$start[i],
                          model$exons$end[i]), "")
  transcript <- if (model$strand == "-") paste0(vapply(rev(exon_forward),
    reconstruction_reverse_complement, ""), collapse = "") else
    paste0(exon_forward, collapse = "")
  terminal <- if (nrow(model$stop_segments)) {
    stop_forward <- vapply(seq_len(nrow(model$stop_segments)), function(i)
      event_reference_fetch(reference, model$contig, model$stop_segments$start[i],
                            model$stop_segments$end[i]), "")
    if (model$strand == "-") paste0(vapply(rev(stop_forward),
      reconstruction_reverse_complement, ""), collapse = "") else
      paste0(stop_forward, collapse = "")
  } else NULL
  phases <- as.integer(seg$phase[ordered]); sizes <- nchar(coding)
  expected <- (3L - c(0L, utils::head(cumsum(sizes), -1L)) %% 3L) %% 3L
  reason <- NULL
  if (anyNA(phases) || phases[[1L]] != 0L || !model$start_codon ||
      !model$stop_codon || model$incomplete)
    reason <- "INCOMPLETE_CDS"
  else if (!identical(phases, as.integer(expected)))
    reason <- "GTF_PHASE_INCONSISTENT"
  else if (model$contig == "chrM") reason <- "MITOCHONDRIAL_CODE"
  else if (model$transcript_id %in% sec_ids) reason <- "SELENOCYSTEINE"
  else if (nchar(cds) %% 3L) reason <- "INCOMPLETE_CDS"
  else if (!terminal %in% c("TAA", "TAG", "TGA")) reason <- "REFERENCE_TERMINAL_STOP_EXCEPTION"
  else if (is.na(model$protein_id)) reason <- "PROTEIN_RELATIONSHIP_MISSING"
  else if (!is.null(transcript_catalog) &&
           !identical(transcript_catalog[[model$transcript_id]], transcript))
    reason <- "REFERENCE_TRANSCRIPT_MISMATCH"
  protein <- if (is.null(reason)) reconstruction_translate(cds) else NULL
  if (is.null(reason) && (is.null(protein) || !startsWith(cds, "ATG") ||
                          grepl("*", protein, fixed = TRUE))) reason <- "REFERENCE_ORF_EXCEPTION"
  catalog <- if (is.na(model$protein_id)) NULL else protein_catalog[[model$protein_id]]
  if (is.null(reason) && (is.null(catalog) || !identical(sub("\\*$", "", catalog), protein)))
    reason <- "REFERENCE_PROTEIN_MISMATCH"
  list(cds = cds, protein = if (is.null(reason)) protein else NULL,
       transcript = transcript, terminal_codon = terminal, reason = reason,
       forward = forward, exon_forward = exon_forward)
}

reconstruction_protein_links <- function(path, transcript_ids) {
  con <- reference_connection(path); on.exit(close(con))
  links <- list()
  repeat {
    lines <- readLines(con, n = 20000L, warn = FALSE)
    if (!length(lines)) break
    for (line in lines[startsWith(lines, ">")]) {
      fields <- strsplit(substring(line, 2L), "|", fixed = TRUE)[[1L]]
      if (length(fields) >= 2L && fields[[2L]] %in% transcript_ids)
        links[[fields[[2L]]]] <- fields[[1L]]
    }
  }
  links
}

reconstruction_relevant <- function(model, rows) {
  Filter(function(row) {
    if (!identical(row$contig, model$contig) ||
        !row$event_type %in% c("SNV", "MNV", "INS", "DEL", "INDEL")) return(FALSE)
    e <- reconstruction_edit_core(row)
    if (is.null(e)) return(FALSE)
    seg <- rbind(model$segments[, c("start", "end"), drop = FALSE],
                 model$stop_segments)
    any(if (e$insertion) e$start >= seg$start & e$start <= seg$end + 1L else
      e$start <= seg$end & e$end >= seg$start)
  }, rows)
}

reconstruction_homozygous <- function(row) {
  a <- row$event_details_json$phase_assertions
  length(a) > 0L && !identical(row$phase_status, "NOT_EVALUABLE") &&
    all(vapply(a, function(x) identical(x$genotype_state, "BOTH"), TRUE))
}

reconstruction_layer <- function(model, reference, rows) {
  seg <- model$segments
  projected <- lapply(rows, function(row) reconstruction_project_edit(reconstruction_edit_core(row), seg))
  bad <- vapply(projected, function(e) e$status != "PROJECTED", TRUE)
  if (any(bad)) return(list(status = projected[[which(bad)[1L]]]$status))
  sequences <- reference$forward
  for (i in seq_along(sequences)) {
    edits <- Filter(function(e) e$segment == i, projected)
    result <- reconstruction_apply_edits(sequences[[i]], edits)
    if (result$status != "COMPLETE") return(list(status = result$status))
    sequences[[i]] <- result$sequence
  }
  if (model$strand == "-") sequences <- vapply(rev(sequences), reconstruction_reverse_complement, "")
  list(status = "COMPLETE", cds = paste0(sequences, collapse = ""), edits = projected)
}

reconstruction_context_id <- function(analysis_id, sample_id, transcript_id, germline_ids,
                                      somatic_ids, phase_state, phase_context) {
  paste0("reconstruction:", digest::digest(charToRaw(canonical_bytes(list(
    policy = reconstruction_policy(), analysis_id = analysis_id, sample_id = sample_id,
    transcript_id = transcript_id, germline_ids = as.list(sort(germline_ids, method = "radix")),
    somatic_ids = as.list(sort(somatic_ids, method = "radix")), phase_state = phase_state,
    phase_context = phase_context))),
    algo = "sha256", serialize = FALSE))
}

reconstruction_protein_change <- function(normal, tumor) {
  a <- strsplit(normal, "", fixed = TRUE)[[1L]]
  b <- strsplit(tumor, "", fixed = TRUE)[[1L]]
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
  list(first_position = if (identical(normal, tumor)) NULL else left + 1L,
       last_normal_position = if (identical(normal, tumor)) NULL else length(a) - right,
       last_tumor_position = if (identical(normal, tumor)) NULL else length(b) - right,
       removed = old, inserted = new)
}

reconstruction_orf_status <- function(cds, protein) {
  if (nchar(cds) %% 3L) return("FRAMESHIFT")
  if (!startsWith(cds, "ATG")) return("START_LOSS")
  if (is.null(protein)) return("TRANSLATION_EXCEPTION")
  if (grepl("*", protein, fixed = TRUE)) return("STOP_GAIN")
  "COMPLETE"
}

reconstruction_tumor_rows <- function(germline, somatic) {
  out <- germline
  if (!length(somatic)) return(list(status = "COMPLETE", rows = out))
  for (s in somatic) {
    se <- reconstruction_edit_core(s)
    for (i in rev(seq_along(out))) {
      if (!identical(out[[i]]$event_role, "germline")) next
      ge <- reconstruction_edit_core(out[[i]])
      intersect <- if (se$insertion && ge$insertion) se$start == ge$start else
        if (se$insertion) se$start > ge$start && se$start <= ge$end else
          if (ge$insertion) ge$start > se$start && ge$start <= se$end else
            se$start <= ge$end && ge$start <= se$end
      if (!intersect) next
      simple_same_base <- !se$insertion && !ge$insertion &&
        se$start == se$end && ge$start == ge$end && se$start == ge$start &&
        nchar(se$new) == 1L && nchar(ge$new) == 1L
      if (!simple_same_base && !identical(c(se$start, se$end, se$new),
                                          c(ge$start, ge$end, ge$new)))
        return(list(status = "LAYERED_OVERLAP_UNRESOLVED"))
      # Both alleles were asserted against the reference. The somatic allele
      # uniquely determines the tumor base at a simple shared substitution.
      out[[i]] <- NULL
    }
    out[[length(out) + 1L]] <- s
  }
  list(status = "COMPLETE", rows = out)
}

reconstruction_normal_backgrounds <- function(germline, anchors, rows,
                                              opposite = FALSE, sole_both = FALSE) {
  homo <- Filter(reconstruction_homozygous, germline)
  hetero <- Filter(function(x) !reconstruction_homozygous(x), germline)
  background <- function(known, uncertain = FALSE) list(
    rows = c(homo, known),
    phase_state = if (uncertain) "BACKGROUND_UNKNOWN" else "RESOLVED_OR_UNAMBIGUOUS")
  if (sole_both) {
    if (!length(hetero)) return(list(background(list())))
    if (length(hetero) == 1L)
      return(list(background(list()), background(hetero)))
    first <- hetero[[1L]]
    relations <- vapply(hetero[-1L], function(g)
      phase_relative_events(rows, first$event_id, g$event_id)$relation, "")
    known <- all(relations %in% c("SAME_HAPLOTYPE", "OPPOSITE_HAPLOTYPE"))
    if (known) for (i in seq_along(hetero)) for (j in seq_along(hetero)) {
      if (i >= j) next
      expected <- if (i == 1L) relations[[j - 1L]] else
        if (identical(relations[[i - 1L]], relations[[j - 1L]]))
          "SAME_HAPLOTYPE" else "OPPOSITE_HAPLOTYPE"
      if (!identical(phase_relative_events(rows, hetero[[i]]$event_id,
        hetero[[j]]$event_id)$relation, expected)) known <- FALSE
    }
    if (!known) return(list(background(list(), TRUE)))
    same <- c(list(first), hetero[-1L][relations == "SAME_HAPLOTYPE"])
    other <- hetero[-1L][relations == "OPPOSITE_HAPLOTYPE"]
    return(list(background(same), background(other)))
  }
  known <- list(); uncertain <- FALSE
  for (g in hetero) {
    relations <- vapply(anchors, function(s)
      phase_relative_events(rows, g$event_id, s$event_id)$relation, "")
    target <- if (opposite) "OPPOSITE_HAPLOTYPE" else "SAME_HAPLOTYPE"
    excluded <- if (opposite) "SAME_HAPLOTYPE" else "OPPOSITE_HAPLOTYPE"
    if (length(relations) && all(relations == target))
      known[[length(known) + 1L]] <- g
    else if (!length(relations) || !all(relations == excluded)) uncertain <- TRUE
  }
  list(background(known, uncertain))
}

reconstruction_contexts <- function(project, bundle, rows) {
  if (!length(rows)) return(list())
  ids <- reconstruction_candidate_ids(bundle, rows)
  if (!length(ids)) return(list())
  models <- Filter(Negate(is.null), reconstruction_models(bundle, ids))
  if (!length(models)) return(list())
  missing_protein <- vapply(models, function(x) is.na(x$protein_id), TRUE)
  if (any(missing_protein)) {
    links <- reconstruction_protein_links(bundle$resources$protein_fasta$resolved_path,
      vapply(models[missing_protein], `[[`, "", "transcript_id"))
    for (i in which(missing_protein)) {
      linked <- links[[models[[i]]$transcript_id]]
      if (!is.null(linked)) models[[i]]$protein_id <- linked
    }
  }
  catalog <- reconstruction_fasta_select(bundle$resources$protein_fasta$resolved_path,
    unique(stats::na.omit(vapply(models, `[[`, "", "protein_id"))))
  transcripts <- reconstruction_fasta_select(bundle$resources$transcript_fasta$resolved_path,
    unique(vapply(models, `[[`, "", "transcript_id")))
  reference <- event_reference_open(bundle); on.exit(event_reference_close(reference))
  output <- list(); oracle <- character()
  for (model in models) {
    relevant_all <- reconstruction_relevant(model, rows)
    somatic_all <- Filter(function(x) identical(x$event_role, "somatic"), relevant_all)
    if (!length(somatic_all)) next
    ref <- reconstruction_reference(model, reference, catalog,
      bundle$sec$transcript_id, transcripts)
    oracle <- c(oracle, if (is.null(ref$reason)) "exact" else ref$reason)
    sample_ids <- sort(unique(vapply(somatic_all, `[[`, "", "sample_id")), method = "radix")
    for (sample_id in sample_ids) {
    relevant <- Filter(function(x) identical(x$sample_id, sample_id), relevant_all)
    somatic <- Filter(function(x) identical(x$event_role, "somatic"), relevant)
    homo_somatic <- Filter(reconstruction_homozygous, somatic)
    seeds <- Filter(function(x) !reconstruction_homozygous(x), somatic)
    groups <- list()
    for (seed in seeds) {
      peers <- Filter(function(x) identical(x$event_id, seed$event_id) ||
        identical(phase_relative_events(rows, seed$event_id, x$event_id)$relation,
                  "SAME_HAPLOTYPE"), seeds)
      # A contradictory pair cannot be installed on one side.
      if (length(peers) > 1L) for (a in peers) for (b in peers)
        if (a$event_id != b$event_id &&
            !identical(phase_relative_events(rows, a$event_id, b$event_id)$relation,
                       "SAME_HAPLOTYPE")) peers <- list(seed)
      peers <- c(peers, Filter(function(x) !x$event_id %in%
        vapply(peers, `[[`, "", "event_id"), homo_somatic))
      key <- paste(sort(vapply(peers, `[[`, "", "event_id"), method = "radix"), collapse = "/")
      groups[[key]] <- list(events = peers, anchors = Filter(function(x)
        !reconstruction_homozygous(x), peers), opposite = FALSE, sole_both = FALSE)
    }
    if (length(homo_somatic)) {
      all_cis <- length(seeds) == 1L
      if (length(seeds) > 1L) {
        all_cis <- TRUE
        for (i in seq_len(length(seeds) - 1L))
          for (j in seq.int(i + 1L, length(seeds)))
            if (!identical(phase_relative_events(rows, seeds[[i]]$event_id,
              seeds[[j]]$event_id)$relation, "SAME_HAPLOTYPE")) all_cis <- FALSE
      }
      if (!length(seeds) || all_cis) {
        groups[[paste(sort(vapply(homo_somatic, `[[`, "", "event_id"),
          method = "radix"), collapse = "/")]] <- list(events = homo_somatic,
          anchors = seeds, opposite = length(seeds) > 0L,
          sole_both = !length(seeds))
      }
    }
    for (spec in groups) {
      group <- spec$events
      germline <- Filter(function(x) identical(x$event_role, "germline"), relevant)
      backgrounds <- reconstruction_normal_backgrounds(germline, spec$anchors, rows,
        opposite = spec$opposite, sole_both = spec$sole_both)
      for (background in backgrounds) {
      normal_rows <- background$rows
      phase_state <- background$phase_state
      germ_ids <- sort(vapply(normal_rows, `[[`, "", "event_id"), method = "radix")
      som_ids <- sort(vapply(group, `[[`, "", "event_id"), method = "radix")
      phase_rows <- relevant[order(vapply(relevant, `[[`, "", "event_id"), method = "radix")]
      phase_context <- lapply(phase_rows, function(x) list(event_id = x$event_id,
        relation = if (reconstruction_homozygous(x)) "BOTH" else
          if (x$event_id %in% som_ids) "SAME_HAPLOTYPE" else
            phase_relative_events(rows, group[[1L]]$event_id, x$event_id)$relation))
      result <- list(reconstruction_id = reconstruction_context_id(project$analysis_id,
        group[[1L]]$sample_id, model$transcript_id, germ_ids, som_ids, phase_state,
        phase_context),
        analysis_id = project$analysis_id, sample_id = group[[1L]]$sample_id,
        transcript_id = model$transcript_id, protein_id = model$protein_id,
        gene_id = model$gene_id, contig = model$contig, strand = model$strand,
        mane_status = model$mane, reference_cds = ref$cds,
        reference_transcript = ref$transcript,
        reference_protein = ref$protein, reference_terminal_codon = ref$terminal_codon,
        germline_event_ids = germ_ids,
        somatic_event_ids = som_ids, phase_state = phase_state,
        phase_context = phase_context, reconstruction_policy = reconstruction_policy(),
        reference_resources = bundle$resource_identity, status = "DEFERRED",
        reason = ref$reason)
      stop_event <- any(vapply(c(group, normal_rows), function(s) {
        e <- reconstruction_edit_core(s)
        nrow(model$stop_segments) && any(e$start <= model$stop_segments$end &
          e$end >= model$stop_segments$start)
      }, TRUE))
      if (stop_event && !is.null(ref$forward)) {
        result$reason <- if (is.null(ref$reason)) "TERMINAL_STOP_EDIT" else ref$reason
        tx_model <- model; tx_model$segments <- model$exons
        tx_ref <- list(forward = ref$exon_forward)
        normal_tx <- reconstruction_layer(tx_model, tx_ref, normal_rows)
        tumor_rows <- reconstruction_tumor_rows(normal_rows, group)
        tumor_tx <- if (tumor_rows$status == "COMPLETE")
          reconstruction_layer(tx_model, tx_ref, tumor_rows$rows) else
          list(status = tumor_rows$status)
        result$patient_normal_transcript <- normal_tx$cds
        result$tumor_transcript <- tumor_tx$cds
        stop_model <- model
        stop_model$segments <- model$stop_segments
        stop_ref <- list(forward = vapply(seq_len(nrow(model$stop_segments)), function(i)
          event_reference_fetch(reference, model$contig, model$stop_segments$start[i],
                                model$stop_segments$end[i]), ""))
        touches_stop <- function(x) {
          e <- reconstruction_edit_core(x)
          any(e$start <= model$stop_segments$end & e$end >= model$stop_segments$start)
        }
        touches_cds <- function(x) {
          e <- reconstruction_edit_core(x)
          any(if (e$insertion) e$start >= model$segments$start &
            e$start <= model$segments$end + 1L else
            e$start <= model$segments$end & e$end >= model$segments$start)
        }
        normal_stop <- reconstruction_layer(stop_model, stop_ref, Filter(touches_stop, normal_rows))
        tumor_stop <- if (tumor_rows$status == "COMPLETE")
          reconstruction_layer(stop_model, stop_ref,
            Filter(touches_stop, tumor_rows$rows)) else list(status = tumor_rows$status)
        if (normal_stop$status == "COMPLETE") result$patient_normal_terminal_codon <- normal_stop$cds
        if (tumor_stop$status == "COMPLETE") result$tumor_terminal_codon <- tumor_stop$cds
        if (normal_stop$status != "COMPLETE" || tumor_stop$status != "COMPLETE") {
          result$terminal_reason <- "TERMINAL_STOP_PROJECTION_UNRESOLVED"
          if (is.null(ref$reason)) result$reason <- result$terminal_reason
        }
        else if (normal_stop$cds %in% c("TAA", "TAG", "TGA") &&
                 !tumor_stop$cds %in% c("TAA", "TAG", "TGA")) {
          result$terminal_reason <- "STOP_LOSS"
          if (is.null(ref$reason)) result$reason <- result$terminal_reason
        }
        if (any(vapply(c(group, normal_rows), function(x)
          touches_stop(x) && touches_cds(x), TRUE))) {
          result$reason <- "CDS_TERMINAL_BOUNDARY_INTERSECTION"
        } else {
          normal_cds <- reconstruction_layer(model, ref,
            Filter(function(x) !touches_stop(x), normal_rows))
          tumor_cds <- if (tumor_rows$status == "COMPLETE")
            reconstruction_layer(model, ref,
              Filter(function(x) !touches_stop(x), tumor_rows$rows)) else
              list(status = tumor_rows$status)
          result$patient_normal_cds <- normal_cds$cds
          result$tumor_cds <- tumor_cds$cds
          if (normal_cds$status != "COMPLETE" || tumor_cds$status != "COMPLETE")
            result$reason <- "CDS_PROJECTION_UNRESOLVED"
        }
      }
      if (!is.null(ref$forward) && !stop_event) {
        normal <- reconstruction_layer(model, ref, normal_rows)
        tumor_rows <- reconstruction_tumor_rows(normal_rows, group)
        tumor <- if (tumor_rows$status == "COMPLETE")
          reconstruction_layer(model, ref, tumor_rows$rows) else list(status = tumor_rows$status)
        result$patient_normal_cds <- normal$cds
        result$tumor_cds <- tumor$cds
        tx_model <- model; tx_model$segments <- model$exons
        tx_ref <- list(forward = ref$exon_forward)
        normal_tx <- reconstruction_layer(tx_model, tx_ref, normal_rows)
        tumor_tx <- if (tumor_rows$status == "COMPLETE")
          reconstruction_layer(tx_model, tx_ref, tumor_rows$rows) else
          list(status = tumor_rows$status)
        result$patient_normal_transcript <- normal_tx$cds
        result$tumor_transcript <- tumor_tx$cds
        if (normal$status != "COMPLETE") result$reason <- normal$status
        else if (tumor$status != "COMPLETE") result$reason <- tumor$status
        else if (normal_tx$status != "COMPLETE") result$reason <- normal_tx$status
        else if (tumor_tx$status != "COMPLETE") result$reason <- tumor_tx$status
        else if (!is.null(ref$reason)) result$reason <- ref$reason
        else {
          # Opposing frame-changing edits can restore total CDS length while
          # leaving a shifted reading frame between their original positions.
          frame_edit <- any(vapply(c(normal$edits, tumor$edits), function(e)
            e$delta %% 3L != 0L, TRUE))
          if (frame_edit) result$reason <- "FRAMESHIFT" else {
            np <- reconstruction_translate(normal$cds)
            tp <- reconstruction_translate(tumor$cds)
            ns <- reconstruction_orf_status(normal$cds, np)
            ts <- reconstruction_orf_status(tumor$cds, tp)
            result$reason <- if (ns != "COMPLETE") ns else
              if (ts != "COMPLETE") ts else NULL
            if (is.null(result$reason)) {
              result$patient_normal_protein <- np
              result$tumor_protein <- tp
              result$protein_change <- reconstruction_protein_change(np, tp)
              result$status <- "COMPLETE"
            }
          }
        }
      }
      output[[length(output) + 1L]] <- result
      }
    }
    }
  }
  standard_oracle <- oracle[!oracle %in% c("INCOMPLETE_CDS", "MITOCHONDRIAL_CODE",
    "SELENOCYSTEINE", "REFERENCE_ORF_EXCEPTION", "REFERENCE_TERMINAL_STOP_EXCEPTION")]
  if (length(standard_oracle) >= 3L && all(standard_oracle %in% c(
      "REFERENCE_PROTEIN_MISMATCH", "REFERENCE_TRANSCRIPT_MISMATCH")))
    stop("Bundle-wide reference protein incompatibility")
  output[order(vapply(output, `[[`, "", "reconstruction_id"), method = "radix")]
}
