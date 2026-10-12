# Coordinates are one-based genomic bases in transcript orientation. The
# extension never uses caller-provided fusion sequence as sequence authority.
source_model_catalog <- function(bundle, events) {
  annotation <- bundle$annotation
  genes <- unique(unlist(lapply(events, function(e) {
    if (e$source_type == "fusion") vapply(e$structure$genes, function(g)
      if (is.null(g$id)) "" else g$id, "") else e$structure$gene_id
  })))
  ids <- unique(annotation$transcript_id[annotation$gene_id %in% genes])
  ids <- ids[!is.na(ids)]
  named <- unique(unlist(lapply(events, function(e) if (e$source_type == "fusion")
    vapply(e$structure$genes, `[[`, "", "name") else character())))
  if (length(named)) {
    gene_name <- stringi::stri_match_first_regex(annotation$raw_attributes,
      'gene_name "([^"]+)"')[, 2L]
    ids <- union(ids, annotation$transcript_id[gene_name %in% named])
  }
  ids <- sort(unique(ids[!is.na(ids)]), method = "radix")
  models <- reconstruction_models(bundle, ids)
  Filter(Negate(is.null), models)
}

source_model_sequence <- function(model, bundle, reference,
                                  transcript_catalog = NULL) {
  if (is.null(model$model_reason) && model$contig != "chrM") {
    coordinates <- translation_coordinates(model$exons, model$strand)
    if (!length(coordinates)) return(NULL)
    if (is.null(transcript_catalog)) transcript_catalog <-
      reconstruction_fasta_select(bundle$resources$transcript_fasta$resolved_path,
        model$transcript_id)
    sequence <- transcript_catalog[[model$transcript_id]]
    if (is.null(sequence) || nchar(sequence) != length(coordinates)) return(NULL)
    exonic <- vapply(seq_len(nrow(model$exons)), function(i)
      event_reference_fetch(reference, model$contig, model$exons$start[[i]],
        model$exons$end[[i]]), "")
    order <- if (model$strand == "-") rev(seq_along(exonic)) else seq_along(exonic)
    genome_sequence <- paste0(if (model$strand == "-")
      vapply(exonic[order], reconstruction_reverse_complement, "") else
      exonic[order], collapse = "")
    if (!identical(sequence, genome_sequence))
      stop("Source transcript FASTA/genome mismatch: ", model$transcript_id)
    cds <- translation_coordinates(model$segments, model$strand)
    start <- match(cds[[1L]], coordinates)
    if (is.na(start) || !isTRUE(model$start_codon) ||
        substring(sequence, start, start + 2L) != "ATG") return(NULL)
    sec_aa <- as.integer(bundle$sec$position[bundle$sec$transcript_id ==
      model$transcript_id])
    sec_sites <- lapply(sec_aa, function(aa) {
      p <- start + (aa - 1L) * 3L
      if (p + 2L > length(coordinates)) return(NULL)
      list(coordinates = as.integer(coordinates[p + 0:2]),
        id = paste0("GENCODE_SEC:", model$transcript_id, ":", aa),
        contig = model$contig)
    })
    sec_sites <- Filter(Negate(is.null), sec_sites)
    list(model = model, sequence = sequence, coordinates = coordinates,
      coding_start = as.integer(start), coding_start_coordinate = cds[[1L]],
      sec_sites = sec_sites)
  } else NULL
}

source_gene_matches <- function(model, gene, bundle) {
  if (!is.null(gene$id) && nzchar(gene$id) && gene$id != ".")
    return(identical(model$gene_id, gene$id))
  x <- bundle$annotation$raw_attributes[bundle$annotation$transcript_id ==
    model$transcript_id]
  any(stringi::stri_match_first_regex(x, 'gene_name "([^"]+)"')[, 2L] ==
    gene$name, na.rm = TRUE)
}

source_fusion_structure <- function(event, left, right, bundle) {
  s <- event$structure
  if (!source_gene_matches(left$model, s$genes[[1L]], bundle) ||
      !source_gene_matches(right$model, s$genes[[2L]], bundle) ||
      left$model$contig != s$left$contig || right$model$contig != s$right$contig ||
      left$model$strand != s$left$strand || right$model$strand != s$right$strand)
    return(NULL)
  a <- match(s$left$position, left$coordinates)
  b <- match(s$right$position, right$coordinates)
  if (is.na(a) || is.na(b) || b > length(right$coordinates))
    return(NULL)
  # Left breakpoint is the last retained exonic base; right is the first.
  coords <- c(left$coordinates[seq_len(a)], right$coordinates[seq.int(b,
    length(right$coordinates))])
  sequence <- paste0(substring(left$sequence, 1L, a),
    substring(right$sequence, b), collapse = "")
  if (identical(sequence, left$sequence) || identical(sequence, right$sequence))
    return(NULL)
  start <- left$coding_start
  retained <- a >= start + 2L &&
    identical(substring(sequence, start, start + 2L), "ATG")
  list(sequence = sequence, left_length = as.integer(a), coordinates = as.integer(coords),
    contigs = c(rep(left$model$contig, a), rep(right$model$contig,
      length(right$coordinates) - b + 1L)), coding_start = start,
    retained_start = retained,
    status = if (retained) "RECONSTRUCTED" else "ORF_UNAVAILABLE_START_LOST",
    junction_nt = as.integer(a),
    novel_nt_intervals = list(),
    blocks = list(list(partner = "5prime", transcript_id = left$model$transcript_id,
      reference_coordinates = as.list(left$coordinates[seq_len(a)])),
      list(partner = "3prime", transcript_id = right$model$transcript_id,
        reference_coordinates = as.list(right$coordinates[seq.int(b,
          length(right$coordinates))]))),
    source_models = list(left = left, right = right), source_start = b)
}

source_splice_structure <- function(event, model, reference) {
  s <- event$structure
  if (model$model$gene_id != s$gene_id || model$model$contig != s$contig ||
      model$model$strand != s$strand ||
      (!is.null(s$transcript_id) && model$model$transcript_id != s$transcript_id))
    return(NULL)
  a <- match(s$donor, model$coordinates)
  b <- match(s$acceptor, model$coordinates)
  extra_left <- integer(); extra_right <- integer()
  step <- if (s$strand == "+") 1L else -1L
  exons <- model$model$exons
  ordered <- if (s$strand == "-") rev(seq_len(nrow(exons))) else seq_len(nrow(exons))
  boundaries <- vapply(ordered, function(i) if (s$strand == "+")
    exons$end[[i]] else exons$start[[i]], 0)
  starts <- vapply(ordered, function(i) if (s$strand == "+")
    exons$start[[i]] else exons$end[[i]], 0)
  if (s$source_type == "alternative_donor" && is.na(a)) {
    for (j in seq_len(length(ordered) - 1L)) {
      begin <- boundaries[[j]] + step
      limit <- starts[[j + 1L]] - step
      if ((limit - begin) * step < 0L) next
      if (s$donor %in% seq.int(begin, limit, by = step)) {
        a <- match(boundaries[[j]], model$coordinates)
        extra_left <- seq.int(begin, s$donor, by = step)
        break
      }
    }
  }
  if (s$source_type == "alternative_acceptor" && is.na(b)) {
    for (j in seq_len(length(ordered) - 1L)) {
      begin <- boundaries[[j]] + step
      limit <- starts[[j + 1L]] - step
      if ((limit - begin) * step < 0L) next
      if (s$acceptor %in% seq.int(begin, limit, by = step)) {
        b <- match(starts[[j + 1L]], model$coordinates)
        extra_right <- seq.int(s$acceptor, limit, by = step)
        break
      }
    }
  }
  if (is.na(a) || is.na(b) || a >= b) return(NULL)
  if (s$source_type == "exon_skip") {
    x <- match(s$donor, boundaries); y <- match(s$acceptor, starts)
    if (is.na(x) || is.na(y) || y <= x + 1L) return(NULL)
  }
  if (s$source_type == "alternative_donor" &&
      !s$acceptor %in% starts) return(NULL)
  if (s$source_type == "alternative_acceptor" &&
      !s$donor %in% boundaries) return(NULL)
  if (s$source_type == "retained_intron") {
    if (b != a + 1L || is.null(s$retained_start) || is.null(s$retained_end) ||
        s$retained_start > s$retained_end) return(NULL)
    expected <- if (s$strand == "+") seq.int(s$donor + 1L, s$acceptor - 1L) else
      seq.int(s$donor - 1L, s$acceptor + 1L)
    if (!identical(sort(expected), seq.int(s$retained_start, s$retained_end))) return(NULL)
    intron <- event_reference_fetch(reference, s$contig, s$retained_start, s$retained_end)
    if (s$strand == "-") intron <- reconstruction_reverse_complement(intron)
    sequence <- paste0(substring(model$sequence, 1L, a), intron,
      substring(model$sequence, b))
    coords <- c(model$coordinates[seq_len(a)], expected,
      model$coordinates[seq.int(b, length(model$coordinates))])
    junctions <- c(a, a + nchar(intron))
  } else {
    if (b == a + 1L && !length(extra_left) && !length(extra_right))
      return(NULL) # unchanged annotated adjacency
    fetch <- function(coords) {
      if (!length(coords)) return("")
      seq <- event_reference_fetch(reference, s$contig, min(coords), max(coords))
      if (s$strand == "-") reconstruction_reverse_complement(seq) else seq
    }
    sequence <- paste0(substring(model$sequence, 1L, a), fetch(extra_left),
      fetch(extra_right), substring(model$sequence, b))
    coords <- c(model$coordinates[seq_len(a)], extra_left, extra_right,
      model$coordinates[seq.int(b, length(model$coordinates))])
    junctions <- a + length(extra_left)
  }
  start <- match(model$coding_start_coordinate, coords)
  retained <- !is.na(start) && start + 2L <= nchar(sequence) &&
    substring(sequence, start, start + 2L) == "ATG"
  upstream <- model$coordinates[seq_len(a)]
  downstream <- model$coordinates[seq.int(b, length(model$coordinates))]
  skipped <- if (b > a + 1L) model$coordinates[seq.int(a + 1L, b - 1L)] else integer()
  inserted <- if (s$source_type == "retained_intron") expected else
    c(extra_left, extra_right)
  # Newly present sequence content, in altered-transcript coordinates.
  # Deletion-only structures have no inserted interval.
  novel_nt_intervals <- if (length(inserted)) list(list(
    start_nt = as.integer(a + 1L), end_nt = as.integer(a + length(inserted)),
    contig = s$contig, strand = s$strand,
    reference_coordinates = as.list(as.integer(inserted)))) else list()
  blocks <- list(list(role = "upstream", transcript_id = model$model$transcript_id,
    contig = s$contig, strand = s$strand,
    reference_coordinates = as.list(as.integer(upstream))))
  if (length(inserted)) blocks[[length(blocks) + 1L]] <- list(
    role = if (s$source_type == "retained_intron") "retained_intron" else
      "alternative_boundary", contig = s$contig, strand = s$strand,
    reference_coordinates = as.list(as.integer(inserted)))
  blocks[[length(blocks) + 1L]] <- list(role = "downstream",
    transcript_id = model$model$transcript_id, contig = s$contig,
    strand = s$strand, reference_coordinates = as.list(as.integer(downstream)))
  list(sequence = sequence, coordinates = as.integer(coords),
    contigs = rep(s$contig, length(coords)),
    coding_start = if (retained) as.integer(start) else NULL,
    retained_start = retained,
    status = if (retained) "RECONSTRUCTED" else "ORF_UNAVAILABLE_START_LOST",
    junction_nt = as.integer(junctions),
    novel_nt_intervals = novel_nt_intervals,
    source_resume_nt = if (s$source_type == "retained_intron")
      as.integer(a + length(expected) + 1L) else
      as.integer(a + length(extra_left) + length(extra_right) + 1L),
    blocks = blocks, skipped_reference_coordinates = as.list(as.integer(skipped)),
    source_models = list(reference = model), source_start = b)
}

source_orf <- function(event, structure, hypothesis_id) {
  if (identical(structure$status, "ORF_UNAVAILABLE_START_LOST"))
    return(list(status = "UNAVAILABLE", reason = structure$status, peptides = list()))
  start <- structure$coding_start
  if (is.null(start) || !structure$retained_start)
    return(list(status = "UNAVAILABLE", reason = "ANNOTATED_START_NOT_RETAINED", peptides = list()))
  sites <- unlist(lapply(structure$source_models, `[[`, "sec_sites"),
    recursive = FALSE)
  sec_positions <- integer(); sec_ids <- character()
  for (site in sites) {
    possible <- which(structure$coordinates == site$coordinates[[1L]] &
      structure$contigs == site$contig)
    for (p in possible) if (p + 2L <= length(structure$coordinates) &&
      identical(as.integer(structure$coordinates[p + 0:2]),
        as.integer(site$coordinates)) &&
      all(structure$contigs[p + 0:2] == site$contig) &&
      (p - start) %% 3L == 0L) {
      sec_positions <- c(sec_positions, as.integer(p))
      sec_ids <- c(sec_ids, site$id)
    }
  }
  spec <- list(sec_tx_positions = sec_positions,
    translation_exception_ids = sec_ids,
    translation_system = "nuclear_standard", start_source = "GTF_START_CODON",
    terminal_completion = "NONE")
  read <- translation_read(structure$sequence, seq_len(nchar(structure$sequence)),
    start, spec)
  if (is.null(read$protein) || !nzchar(read$protein))
    return(list(status = "UNAVAILABLE", reason = read$boundary, peptides = list()))
  junction_nt <- structure$junction_nt
  junction_aa <- as.integer((junction_nt - start + 1L) %/% 3L)
  junction_aa <- junction_aa[junction_aa >= 1L & junction_aa < nchar(read$protein)]
  frame <- if (event$source_type == "fusion") {
    three <- structure$source_models$right
    original <- three$coding_start
    new_three <- structure$left_length + 1L
    if ((new_three - start - (structure$source_start - original)) %% 3L == 0L)
      "in_frame" else "out_of_frame"
  } else {
    reference <- structure$source_models$reference
    original <- reference$coding_start
    new_three <- structure$source_resume_nt
    if ((new_three - start - (structure$source_start - original)) %% 3L == 0L)
      "in_frame" else "out_of_frame"
  }
  # Only codons that actually translated to residues can contribute novelty.
  # This includes both boundary codons and codons wholly inside an insertion.
  codon_starts <- start + (seq_len(nchar(read$protein)) - 1L) * 3L
  inserted_aa <- integer()
  for (interval in structure$novel_nt_intervals) {
    inserted_aa <- c(inserted_aa, which(codon_starts <= interval$end_nt &
      codon_starts + 2L >= interval$start_nt))
  }
  inserted_aa <- sort(unique(as.integer(inserted_aa)))
  novel <- inserted_aa
  for (j in junction_nt) {
    cut <- as.integer((j - start + 1L) %/% 3L + 1L)
    if ((j - start + 1L) %% 3L != 0L) novel <- c(novel, cut)
    if (frame == "out_of_frame" && cut <= nchar(read$protein)) novel <- c(novel,
      seq.int(max(1L, cut), nchar(read$protein)))
  }
  novel <- sort(unique(novel[novel >= 1L & novel <= nchar(read$protein)]))
  list(status = "TRANSLATED", protein_origin_id = source_id("protein_origin_id",
    list(hypothesis = hypothesis_id, protein = read$protein)),
    orf_hypothesis_id = source_id("orf_hypothesis_id", list(hypothesis = hypothesis_id,
      start = start, policy = "fusion-splice-orf-v1")), protein = read$protein,
    frame_state = frame, termination = read$boundary,
    stop_nt = read$stop_nt, residual_nt = read$residual_nt,
    translation_exception_ids = as.list(read$sec_sites_used),
    junction_aa = as.list(junction_aa), novel_aa = as.list(novel),
    inserted_sequence_aa = as.list(inserted_aa),
    reference_comparison = if (event$source_type == "fusion")
      "parental_source_mapping" else "source_transcript_mapping")
}

source_peptide_windows <- function(orf, policy) {
  if (orf$status != "TRANSLATED") return(list())
  positions <- unlist(orf$novel_aa)
  inserted_positions <- unlist(orf$inserted_sequence_aa)
  cuts <- unlist(orf$junction_aa)
  n <- nchar(orf$protein)
  out <- list()
  for (length in sort(unique(unlist(policy)), method = "radix")) {
    if (length > n) next
    for (start in seq_len(n - length + 1L)) {
      end <- start + length - 1L
      novel <- any(positions >= start & positions <= end)
      crosses <- any(cuts >= start & cuts < end)
      inserted <- any(inserted_positions >= start & inserted_positions <= end)
      if (!novel && !crosses) next
      peptide <- substring(orf$protein, start, end)
      peptide_sequence_validate(peptide)
      out[[length(out) + 1L]] <- list(sequence = peptide,
        start_aa = as.integer(start), end_aa = as.integer(end),
        novel_residue = novel, junction_crossing = crosses,
        reference_comparison = if (inserted) "inserted_sequence_no_single_reference" else
          if (crosses) "junction_no_single_reference" else
          if (orf$frame_state == "out_of_frame") "novel_frameshift_no_reference" else
          "unresolved")
    }
  }
  out
}
